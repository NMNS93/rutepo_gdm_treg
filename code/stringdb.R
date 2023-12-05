library(STRINGdb)
library(png)
library(grid)
library(igraph)

make_sdb_object <- function(df, thresh=200){
  # Df must have 'gene' column with relevant genes
  # Package info says default threshold is 400
  sdb <- STRINGdb$new(version="11.5", score_threshold=thresh, species=9606)
  mapped <- sdb$map( df, "gene", removeUnmappedRows = TRUE )
  return(list(sdb, mapped))
}

plot_sdb_net <- function(sdb, mapped, dm){
  # Return ggplot object with stringdb plot
  
  # plot tempfile
  tmpf = tempfile()
  png(filename=tmpf, res=300, width=10, height=10, units="in")
  sdb$plot_network(mapped$STRING_id, )
  dev.off()
  
  # load iamge
  img <- readPNG(tmpf)
  g <- rasterGrob(img, interpolate=TRUE)
  
  # plot as ggplot
  ggplot(data.frame(seq(0, 0.1, by=0.1),seq(0, 0.5, by=0.1)), geom="blank") +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    theme_bw()
}

get_string_hubs <- function(sdb, mapped){
  # sdb is string db object
  fg = sdb$get_graph()
  deg = degree(fg)
  deg_df = data.frame(STRING_id = names(deg), degree=deg)
  mapped2 = merge(mapped, deg_df, all.x=T)
  topg = dplyr::top_n(mapped2, 10, degree) %>% pull(gene)
  return(topg)
}

