library(readr)
library(here)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(patchwork)
library(glue)
library(data.table)
library(ModelMetrics)
library(parallel)
library(qs)
library(logger)
library(plotROC)
library(MASS)
library(caret)
#library(ggh4x)

# Set fields from clean metadata to use for modelling
# Select desired fields
meta_select_vars = c(
  "multi_q", # Sample identifier
  "ethnicity", # Ethnicity
  "bmi", # Patient weight
  "map", # mean arterial pressure
  "insulin", # Factor for insulin. Not predictive.
  "metformin", # Factor for metformin. Not predictive.
  "diet", # Factor for diet. Not predictive.
  "aspirin" # Factor for aspirin. Not predictive.
)

# Refactor blood pressure field
blood_pressure_refactor <- function(x){
  # ::: Multiple blood pressure measurements were taken from the right and left arms.
  # ::: but blood pressure if ogten taken from just one arm. Diastoolic and systolic are
  # ::: calculated together and used to define thresholds. 
  # ::: First, get an average sys/dias for the patient and classify high BP using NHS guidelines.
  # ::: Then, calculate the difference between arms and classify using guidelines as arm
  # ::: difference in blood pressure is a suggested biomarker
  get_rowmeans <- function(x) x %>% dplyr::mutate_all(as.numeric) %>% as.matrix %>% rowMeans(na.rm=T)
  sbp = x[, c("rsbp1", "rsbp2", "lsbp1", "lsbp2")] %>% get_rowmeans()
  dbp = x[, c("rdbp1", "rdbp2", "ldbp1", "ldbp2")] %>% get_rowmeans()
  rsbp = x[, c("rsbp1", "rsbp2")] %>% get_rowmeans()
  lsbp = x[, c("lsbp1", "lsbp2")] %>% get_rowmeans()
  rdbp = x[, c("rdbp1", "rdbp2")] %>% get_rowmeans()
  ldbp = x[, c("ldbp1", "ldbp2")] %>% get_rowmeans()
  
  # Calculate mean arterial pressure
  # https://onlinelibrary.wiley.com/doi/full/10.1002/ccd.20217
  pp = (sbp-dbp)
  map = dbp + (1/3*pp)
  
  # Calculate MAP class
  # https://link.springer.com/article/10.1007/s00134-018-5446-8
  map_class = dplyr::case_when(
    map <65 ~ "low",
    between(map, 65, 75) ~ "low-normal",
    between(map, 80, 100) ~ "high-normal",
    map > 100 ~ "high"
    
  )
  
  # Calculate blood pressure class
  # https://www.nhs.uk/common-health-questions/lifestyle/what-is-blood-pressure/
  bp_class = dplyr::case_when(
    sbp < 90 & dbp < 60 ~ "low",
    between(sbp, 90, 120) & between(dbp, 60, 90) ~ "normal",
    sbp > 120 & dbp > 90 ~ "high",
    TRUE ~ "other"
  )
  
  
  # Calculate differences
  sbp_diff = rsbp - lsbp
  dbp_diff = rdbp - ldbp
  
  # Add these features to x
  x = dplyr::mutate(
    x,
    sbp=sbp, dbp=dbp, map=map, bp_class=bp_class, map_class=map_class,
    sbp_diff=sbp_diff, dbp_diff=dbp_diff
  )
  
  return(x)
}


medication_reformat <- function(meta){
  meta$other_medication = str_to_lower(meta$other_medication)
  meta$gdm_treatment = str_to_lower(meta$gdm_treatment)
  
  insulin = dplyr::transmute(
    meta,
    insulin = (
      !is.na(other_medication) & str_detect(other_medication, "insulin") |
        !is.na(gdm_treatment) & str_detect(gdm_treatment, "insulin")
    )
  ) %>% pull(insulin) %>% as.numeric %>% as.factor()
  metformin =dplyr::transmute(
    meta,
    metformin = (
      !is.na(other_medication) & str_detect(other_medication, "metformin") |
        !is.na(gdm_treatment) & str_detect(gdm_treatment, "metformin")
    )
  ) %>% pull(metformin) %>% as.numeric %>% as.factor()
  diet =dplyr::transmute(
    meta,
    diet = (
      !is.na(other_medication) & str_detect(other_medication, "diet") |
        !is.na(gdm_treatment) & str_detect(gdm_treatment, "diet")
    )
  ) %>% pull(diet) %>% as.numeric %>% as.factor()
  
  aspirin = as.numeric(meta$low_dose_aspirin) %>% ifelse(is.na(.), 0, .) %>% as.factor()
  
  meta = dplyr::mutate(
    meta,
    "insulin"=insulin, 
    "metformin"=metformin, 
    "diet"=diet,
    "aspirin"=aspirin
  )
  return(meta)
}


make_avpat <- function(x){
  # Make seuart object of patietn averages
  DefaultAssay(x) = "RNA"
  Idents(x) = "multi_q"
  avp = AverageExpression(x, assays="RNA", return.seurat=T, group.by="multi_q")
  avp$multi_q = rownames(avp@meta.data)
  return(avp)
}

getclusterprop <- function(x){
  x@meta.data %>% 
    dplyr::group_by(multi_q, seurat_clusters, cluster) %>% dplyr::count() %>%
    dplyr::ungroup() %>% dplyr::group_by(multi_q) %>%
    dplyr::mutate(n_cells = sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(prop=n/n_cells)
}

get_pat_average <- function(av_pat_markers, treg, cd4){
  # Get average experssion per patient
  avpat_cd4 <- make_avpat(cd4)
  avpat_treg <- make_avpat(treg)
  # Get average cluster proportion per patient
  propc_cd4 <- getclusterprop(cd4) %>% dplyr::mutate(dataset="cd4")
  propc_treg <- getclusterprop(treg) %>% dplyr::mutate(dataset="treg")
  
  # Create empty dataframe for patient meta results
  res = data.frame(
    multi_q = unique(avpat_cd4$multi_q)
  )
  # Use av_pat_marker names to add fields to dataframe
  for(i in av_pat_markers){
    spl = stringr::str_split(i, "_", n=3, simplify=T) %>% as.character()
    pdata = dplyr::case_when(
      spl[[1]]=="Treg" & spl[[3]]=="p" ~ "propc_treg",
      spl[[1]]=="Treg" & spl[[3]]=="" ~ "avpat_treg",
      spl[[1]]=="CD4" & spl[[3]]=="p" ~ "propc_cd4",
      spl[[1]]=="CD4" & spl[[3]]=="" ~ "avpat_cd4"
    )
    # Convert to symbol to dataset
    x = eval(sym(pdata))
    # Get field as required
    if(spl[[3]] == "p"){
      cluster_gene = spl[[2]]
      prop_d = dplyr::filter(x, str_detect(x$cluster, cluster_gene)) %>%
        dplyr::select(multi_q, prop)
      names(prop_d)[[2]] = str_to_lower(i)
      res = merge(res, prop_d, all=T)
    } else {
      # Extract gene 
      rna_gene = spl[[2]]
      v = x[rna_gene]@assays$RNA@counts
      rna_d = data.frame(
        multi_q = as.character(colnames(v)),
        rna = as.numeric(v)
      )
      names(rna_d)[[2]] = str_to_lower(i)
      res = merge(res, rna_d, all=T)
    }
  }
  
  return(res)
}


parse_coefs <- function(mod_list, label){
  v = lapply(seq_len(length(mod_list)),
             function(x){
               res = summary(mod_list[[x]])
               dt = as.data.frame(res$coefficients)
               names(dt) = c("estimate", "std_error", "z_value", "prob")
               dt$aic = res$aic
               dt$record = x
               dt$feature = rownames(dt)
               return(dt)
             })
  v = Reduce(rbind, v)
  v$model = label
  
  return(v)
}


extract_feature_response <- function(y, label, fits){
  ix = which(fits$model == label)
  x = fits[ix,]$fits[[1]]
  # Return predictions for the feature and response for a series of model fits
  # given a specific feature e.g. bmi as y
  v = lapply(
    seq_len(nrow(x)),
    function(i){
      o = x[i,]
      data.frame(
        control_b = o$pred,
        feature = sapply(o$data, function(x) x[[1]][[y]])
      )
    }
  )
  res = Reduce(rbind, v)
  res$model = label
  res$fname = y
  return(res)
}

plot_meta_lr <- function(meta, feat, extra, linecol="red"){
  names(extra)[[2]] = feat # Fix feature col
  ggplot(meta, aes(x=!!sym(feat), y=control_b)) + 
    geom_smooth(method="glm",  method.args = list(family="binomial"), data=extra, se=F, color=linecol) + 
    geom_jitter(height = 0.001, width=0.001, shape=21, fill="white", color="black", size=3) +
    geom_point(data=extra, color="grey") + theme_light() +
    labs(title=extra$model[[1]])
}

add_case_control_seurat <- function(meta, cd4, treg){
  # Use seurat object annotations to be consistent with patients filtered
  case_control <- unique(
    rbind(
      treg@meta.data[, c("multi_q", "cond_cln")],
      cd4@meta.data[, c("multi_q", "cond_cln")]
    )
  )
  meta <- merge(meta, case_control, all.x=T)
  meta <- dplyr::filter(meta, !is.na(cond_cln))
  meta <- dplyr::arrange(meta, cond_cln, ethnicity)
  return(meta)
}

recurrent_deg_filter <- function(x) x %>% dplyr::filter(is_sig & abs(avg_log2FC) > 0.1) %>% 
  dplyr::select(seurat_clusters, cluster, gene) %>% dplyr::group_by(gene) %>% 
  dplyr::count() %>% dplyr::filter(n>=3) %>% pull(gene)

cluster_deg_filter <- function(x, y) x %>% dplyr::filter(is_sig & abs(avg_log2FC) > 0.1) %>%
  dplyr::filter(stringr::str_detect(cluster, y)) %>% pull(gene) %>% unique()

get_coexp_genes <- function(coexp){
  # TODO: Hardcoded
  ## Validation: You can see coexpr network as top N selected
  ## Fix coexp TOM as it has now moved
  ## coexp@misc$`treg_Naive-Activation-2`$wgcna_net$TOMFiles = "data/TOM/treg_Naive-Activation-2_TOM.rda"
  ## You cansee hub network is selected by topN in the module :: #hdWGCNA::HubGeneNetworkPlot(coexp)
  coexp_topcells = colnames(coexp)[coexp$`Naive-Activation-21` > 5]
  coexp_means = rowMeans(coexp@assays$RNA@counts[,coexp_topcells], na.rm=T)
  head(coexp_means[order(-coexp_means)], 50)
  ma = hdWGCNA::GetModules(coexp)
  ma_genes = dplyr::filter(ma, module=="Naive-Activation-21") %>% arrange(-abs(`kME_Naive-Activation-21`)) %>% head(10) %>% pull(gene_name)
  return(ma_genes)
}

load_stirm_bulk <- function(mrna_tsv, series_matrix){
  # Load mrna tsv
  bulke <- suppressWarnings(data.table::fread(mrna_tsv))
  names(bulke)[1] = "gene"
  
  # Load series matrix
  # Read the data from the file
  data <- readLines(series_matrix)
  
  # Find the rows containing "individual id" and "health status"
  individual_rows <- grep("individual id:", data)
  health_status_rows <- grep("health status:", data)

  # Extract the individual ids and health statuses
  individual_ids <- stringr::str_extract_all(data[individual_rows], "GD\\d+")[[1]]
  health_statuses <- str_extract_all(data[health_status_rows], "GDM|NGT")[[1]]
  
  # Create a data frame
  result_df <- data.table(id = individual_ids, status = health_statuses)
  # Assignments are duplicated. Drop.
  result_df = unique(result_df)
  
  # Convert to long-format with annotations
  exprl <- tidyr::pivot_longer(bulke, -gene, names_to="id", values_to="expr")
  exprl <- merge(exprl, result_df, all.x=T, by=c("id"))
  return(data.table(exprl))
}


# Get LOO splits
loo_split <- function(meta){
  fold_ix <- seq_len(nrow(meta))
  get_fold <- function(i){
    list(
      "ix_eval"=i,
      "eval"=meta[i,],
      "train"=meta[-i,]
    )
  }
  folds <- lapply(fold_ix, get_fold)
  names(folds) = fold_ix
  return(folds)
}

# Source.: https://bookdown.org/ndphillips/YaRrr/logistic-regression-with-glmfamily-binomial.html
logit_to_prob <- function(x)  1 / (1 + exp((-1 * x)))

model_fold_with_metrics <- function(fi){
  # fi is a single fold object
  fold_no = fi$ix_eval
  fold_id = fi$id
  mod = glm(status ~ ., data = fi$train , family = "binomial")
  mod_aic = mod$aic
  x_truth = fi$train$status
  x_prob = mod$fitted.values
  x_pred = ifelse(x_prob > 0.5, 1, 0)
  x_acc = mean(x_pred == x_truth)
  x_auc = ModelMetrics::auc(x_truth, x_prob)
  x_recall = ModelMetrics::sensitivity(x_truth, x_pred)
  x_spec = ModelMetrics::specificity(x_truth, x_pred)
  x_prec = ModelMetrics::precision(x_truth, x_pred)
  x_f1 = ModelMetrics::f1Score(x_truth, x_pred)
  y_truth <- fi$eval$status
  y_prob <- as.numeric(predict(mod, fi$eval))
  y_prob <- logit_to_prob(y_prob)
  y_pred <- ifelse(y_prob > 0.5, 1, 0)
  
  fold_metrics = tibble(
    fold_no, fold_id, model=list(mod), model_aic = mod_aic, x_truth=list(x_truth), x_prob=list(x_prob),
    x_pred=list(x_pred), x_acc, x_auc, x_recall, x_spec, x_prec, x_f1, y_truth, y_prob, y_pred)
  
  return(fold_metrics)
  
}

add_test_metrics <- function(fold_table){
  test_prob = fold_table$y_prob
  test_pred = fold_table$y_pred
  test_truth = fold_table$y_truth
  
  fold_table$all_y_acc = mean(test_truth == test_pred)
  fold_table$all_y_auc = ModelMetrics::auc(test_truth, test_prob)
  fold_table$all_y_aic = mean(fold_table$model_aic)
  fold_table$all_y_recall = ModelMetrics::recall(test_truth, test_pred)
  fold_table$all_y_spec = ModelMetrics::specificity(test_truth, test_pred)
  fold_table$all_y_prec = ModelMetrics::precision(test_truth, test_pred)
  fold_table$all_y_f1 = ModelMetrics::f1Score(test_truth, test_pred)
  
  return(fold_table)
}

fit_and_collect_model_metrics <- function(dataset, dataset_name, var_selectors, merge=F){
  temp_dir <- fs::dir_create(tempfile(tmpdir = "tempdir"), recurse=T)

  all_folds_files <- mclapply(
    seq_along(var_selectors),
    function(i){
      
      ## Log
      var_field = var_selectors[[i]]
      model_name = paste0(var_field, collapse=";")
      log_info("Processing {model_name}")
      
      ## Return empty dataframe if the fields are not present
      no_fields_present = length(intersect(names(dataset), var_field)) < length(var_field)
      if(no_fields_present){
        ## Cache result to preserve memory
        temp_file <- file.path(temp_dir, paste0("cached_folds_", i, ".qs"))
        qsave(data.frame(), temp_file)
        return(temp_file)
      }
      
      ## Subset data for model
      data <- dplyr::select(dataset, all_of(var_field))
      
      ## Generate LOOCV folds
      folds <- loo_split(data)
      
      ## Add sample being evaluated as ID
      names(folds) <- dataset$id
      for (v in names(folds)){ folds[[v]]$id = v }
      
      ## Fit model and predict on each fold
      fold_table <- Reduce(rbind, lapply(folds, model_fold_with_metrics))
      fold_table <- add_test_metrics(fold_table)
      fold_table$data_name = dataset_name
      fold_table$model_name = model_name # Set of variables in this model
      fold_table$all_y_feats = stringr::str_count(unique(model_name), ';')
      
      ## Cache result to preserve memory
      temp_file <- file.path(temp_dir, paste0("cached_folds_", i, ".qs"))
      qsave(fold_table, temp_file)
      return(temp_file)
    }, mc.cores=10
  )
  
  if(merge){
    all_folds_data <- data.table::rbindlist(lapply(all_folds_files, qread))
    unlink(temp_dir, recursive = TRUE)
    return(all_folds_data)
  }
  # Else
  return(all_folds_files)
}

#' Select the top 2000 models based on AUC score and merge datasets.
#' a is a list of *.qs objects with model results
merge_metrics_mem_aware <- function(a, outfile, top=1000){
  # Take whichever is smaller of the given values (top or length of a)
  top = min(top, length(a))
  
  # Define a function to return auc and filename given a file
  # Expects auc score to be the same in the dataset
  get_auc <- function(x){
    xdata = qs::qread(x)[1, "all_y_auc"]
    xdata$file = x
    return(xdata)
  }
  
  # Apply to all files to create selector table
  log_info("Building file selector")
  sel = data.table::rbindlist(
    mclapply(
      a, get_auc, mc.cores=5
    )
  )
  # Rank selection by top 1000 model AUC scores
  ranks <- rank(-sel$all_y_auc, ties.method = "min")
  # Return files matching the top 1000 model AUC scores
  top_files = sel[which(ranks<=top)]$file
  
  # Load these files in a memory-aware manner
  top_fold = qs::qread(top_files[1])
  log_info("Loading top_files")
  for(i in seq(2, length(top_files))){
    if (i %% 100 == 0){
      log_info("Merge metrics mem_aware: {i}")
    }
    top_fold = rbind(top_fold, qs::qread(top_files[i]))
  }
  
  # Save data object to outfile
  log_info("Saving top folds")
  qs::qsave(top_fold, outfile)
  
  # Return data object
  return(top_fold)
}

build_variable_selectors <- function(x, sizes=NULL){
  if(is.null(sizes)){
    sizes = 1:length(x)
  }
  # Generate all combinations
  all_combinations <- unlist(lapply(sizes, function(i) combn(x, i, simplify = FALSE)), recursive = FALSE)
  # Add target variable 'status' to all combinations
  var_selectors <- lapply(all_combinations, function(comb) c("status", comb))
  return(var_selectors)
}


var_sel_to_ids = function(x, y){
  tibble(
    name_1= sapply(x, function(i) paste0(i, collapse=';')),
    name_2= sapply(y, function(i) paste0(i, collapse=';')),
    model_num = paste0("model_", seq_along(x))
  )
}


tally_model_feats <- function(mv, anno, tids, nameid){
  # Create a table tallying features in each model name
  mv = dplyr::select(mv, model_name)
  mv_genes = unique(mv$model_name) %>% stringr::str_split(";")
  uq_genes = setdiff(unique(unlist(mv_genes)), "status")
  pres_vec = lapply(uq_genes, function(x){
    sapply(mv_genes, function(y) any(stringr::str_detect(y, x)))
  })
  names(pres_vec) = (uq_genes)
  for(i in names(pres_vec)){
    mv[[i]] = as.numeric(pres_vec[[i]])
  }
  mv$model_data = anno
  mv = mv %>% tidyr::pivot_longer(cols=!contains("model"), names_to="gene", values_to='present')
  # Add nums
  tids = dplyr::select(tids, model_num, model_name=!!sym(nameid))
  mv = merge(mv, tids, by="model_name", all.x=T)
  return(mv)
}

build_roc_data <- function(strims_top_mods, cd4_mod, cd4s_mod, strims_mod){
  
  cd4_roc = cd4_mod %>% left_join(
    strims_top_mods %>% select(-all_y_auc, -model_name) %>%
      rename(model_name=model_name_cd4),
    by="model_name"
  ) %>% dplyr::filter(!is.na(model_num)) %>%
    dplyr::mutate(dataset="CD4+") %>%
    dplyr::select(dataset, model_num, model_name_cd4=model_name, fold_no, fold_id, all_y_auc, y_prob, y_truth) %>% left_join(
      strims_top_mods %>% select(model_name_strims=model_name, model_name_cd4)
    )
  
  cd4s_roc = cd4s_mod %>% left_join(
    strims_top_mods %>% select(-all_y_auc),
    by="model_name"
  ) %>% dplyr::filter(!is.na(model_num)) %>%
    dplyr::mutate(dataset="CD4+_Stirm") %>%
    dplyr::select(dataset, model_num, model_name_strims=model_name, fold_no, fold_id, all_y_auc, y_prob, y_truth) %>% left_join(
      strims_top_mods %>% select(model_name_strims=model_name, model_name_cd4)
    )
  
  strims_roc = strims_mod %>% left_join(
    strims_top_mods %>% select(-all_y_auc),
    by="model_name"
  ) %>% dplyr::filter(!is.na(model_num)) %>%
    dplyr::mutate(dataset="Strim") %>%
    dplyr::select(dataset, model_num, model_name_strims=model_name, fold_no, fold_id, all_y_auc, y_prob, y_truth) %>% left_join(
      strims_top_mods %>% select(model_name_strims=model_name, model_name_cd4)
    )
  
  roc_data = Reduce(rbind, list(cd4_roc, cd4s_roc, strims_roc))
  return(roc_data)
}

run_robust_reg_gene_status <- function(bp){
  bp_spl = split(bp, paste0(bp$gene, "_", bp$status_c))
  bp_pred = rbindlist(lapply(bp_spl, function(x){
    mapd = dplyr::filter(bp, gene==x$gene[[1]] & status_c == x$status_c[[1]])
    mod = rlm(expr ~ map, mapd)
    pred = predict(mod, mapd)
    mapd$pred = pred
    mapd$mod = list(mod)
    return(mapd)
  }))
  return(bp_pred)
}

make_scrna_model_dataset <- function(x, markers){
  avx <- make_avpat(x)
  avx <- subset(avx, multi_q %in% meta$id)
  avxm <- avx@assays$RNA@data[markers,]
  avxd <- t(avxm)
  avxdt <- data.frame(avxd)
  avxdt$id <- rownames(avxdt)
  rownames(avxdt) <- NULL
  return(avxdt)
}

load_wang_bulk <- function(x){
  wang = readr::read_csv(fs::path(x)) %>%
    dplyr::select(gene=Gene_Name, contains("GDM"), contains("NEG"))
  wangl = wang %>% tidyr::pivot_longer(cols=-gene, names_to="id", values_to="expression")
  wangl$status = ifelse(stringr::str_detect(wangl$id, "GDM"), 1, 0)
  wangl$status_c = ifelse(wangl$status==1, "GDM", "CONTROL")
  return(wangl)
}

metadata_odds_ratio <- function(meta){
  model <- glm(status ~ bmi + ethnicity + map, data = meta, family = "binomial")
  # Obtain odds ratios and CIs
  odds_ratios <- exp(coef(model))
  ci <- exp(confint(model))
  
  # Create a data frame for plotting
  odds_data = exp(cbind(coef(model), confint(model)))
  odds_data = data.frame(odds_data)
  names(odds_data) = c("odds_ratio", "lower_ci", "upper_ci")
  odds_data$variable = rownames(odds_data)
  rownames(odds_data) = NULL
  
  odds_data <- odds_data %>%
    mutate(
      clean = case_when(
        variable == "bmi" ~ "body\nmass\nindex",
        variable == "ethnicitywhite" ~ "ethnicity\n(white)",
        variable == "map" ~ "mean\narterial\npressure",
        TRUE ~ variable  # Keep other variable names as is
      )
    )
 
  return(odds_data) 
}

get_mod_tops <- function(x){
  list(
    'tops'= hget_mod_tops(x),
    'tops_s5'= hget_mod_tops(x, sub_5=T),
    'tops_aic' = hget_mod_tops(x, sub_5=F,'all_y_aic'),
    'tops_aic_s5' = hget_mod_tops(x, sub_5=T, 'all_y_aic')
  )
}

hget_mod_tops <- function(x, NN=50, sub_5=F, rankmet='all_y_auc'){
  
  if(sub_5){
    # Models with fewer than 5 genes
    x = x %>% dplyr::filter(all_y_feats < 5) 
  }
  
  mod_stats = unique(x %>% dplyr::select(contains("all_y"), data_name, model_name))
  if(rankmet=='all_y_auc'){
    mod_stats = mod_stats %>% dplyr::top_n(NN, get(rankmet))
  }else{ #all_y_aic
    mod_stats = mod_stats %>% dplyr::top_n(NN, -get(rankmet))
  }
  mod_stats$rankmet = rankmet
  
  return(mod_stats)
}

roc_mini <- function(x, data_name){ 
  x %>% mutate(dataset=data_name) %>% dplyr::select(dataset, model_name, fold_no, fold_id, all_y_auc, all_y_f1, y_prob, y_truth)
}

smod_to_roc_ss <- function(wang_smod, stirms_smod, cd4_smod, treg_smod){
  
  roc_data = Reduce(
    rbind,
    list(
      roc_mini(wang_smod, "WANG"),
      roc_mini(stirms_smod, "STIRM"),
      roc_mini(cd4_smod,"CD4+"),
      roc_mini(treg_smod, "Treg")
    )
  )
  
  top_ss_mod = roc_data %>% dplyr::filter(dataset!="cd4_wide") %>% 
    dplyr::select(dataset, model_name, all_y_auc) %>% unique() %>% 
    dplyr::group_by(model_name) %>% dplyr::summarise(mauc=median(all_y_auc)) %>% 
    dplyr::filter(mauc>=0.6)
  
  roc_data_ss = dplyr::filter(roc_data, model_name %in% top_ss_mod$model_name) %>%
    dplyr::mutate(model_name=stringr::str_remove(model_name, "status;")) %>%
    dplyr::mutate(dataset=stringr::str_to_upper(stringr::str_remove(dataset, "_.*"))) %>%
    dplyr::mutate(
      dataname = dplyr::case_when(
        dataset == "TREG" ~ "Treg scRNA",
        dataset == "CD4" ~ "CD4+ T scRNA",
        dataset == "STIRMS" ~ "Stirm et al. 2018",
        dataset == "WANG" ~ "Wang et al. 2021"
      )
    )
  
  return(roc_data_ss)
}

smod_simple_format <- function(wang_smod, stirms_smod, cd4_smod, treg_smod){
  
  roc_data = Reduce(
    rbind,
    list(
      roc_mini(wang_smod, "WANG"),
      roc_mini(stirms_smod, "STIRM"),
      roc_mini(cd4_smod,"CD4+"),
      roc_mini(treg_smod, "Treg")
    )
  )
  
  roc_data_ss = dplyr::filter(roc_data) %>%
    dplyr::mutate(model_name=stringr::str_remove(model_name, "status;")) %>%
    dplyr::mutate(dataset=stringr::str_to_upper(stringr::str_remove(dataset, "_.*"))) %>%
    dplyr::mutate(
      dataname = dplyr::case_when(
        dataset == "TREG" ~ "Treg scRNA",
        dataset == "CD4" ~ "CD4+ T scRNA",
        dataset == "STIRMS" ~ "Stirm et al. 2018",
        dataset == "WANG" ~ "Wang et al. 2021"
      )
    )
  
  return(roc_data_ss)
}


get_top_rgenes_bulk <- function(x){
  x = dplyr::filter(x, dataset %in% c("WANG", "STIRMS"))
  x = unique(dplyr::select(x, dataset, model_name, all_y_auc, all_y_f1))
  x$score = x$all_y_auc * x$all_y_f1
  xw = tidyr::pivot_wider(x, id_cols=c('model_name'), names_from='dataset', values_from=all_y_auc)
  xw$asc = rowMeans(xw[, c("WANG", "STIRMS")], na.rm=T)
  xw = dplyr::arrange(xw, -asc)
  return(xw)
}
