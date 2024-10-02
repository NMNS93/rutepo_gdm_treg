from cellphonedb.utils import db_utils
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import pickle
import os
# -- Version of the databse
cpdb_version = 'v5.0.0'
# -- Path where the input files to generate the database are located
cpdb_target_dir = os.path.join('/genomics/projects/B061-PanicosShangaris/nana/data/cellphonedb', cpdb_version)
db_utils.download_database(cpdb_target_dir, cpdb_version)

# Paths
cpdb_file_path = '/genomics/projects/B061-PanicosShangaris/nana/data/cellphonedb/v5.0.0/cellphonedb.zip'
deg_file_path = '/genomics/projects/B061-PanicosShangaris/nana/data/cellphonedb/cd4_degs.txt'

# Define output paths for pickled results
pickle_out_path_case = '/genomics/projects/B061-PanicosShangaris/nana/gdm/results_v3/cellphonedb/case/cpdb_results_case.pkl'
pickle_out_path_control = '/genomics/projects/B061-PanicosShangaris/nana/gdm/results_v3/cellphonedb/control/cpdb_results_control.pkl'
stat_out_path_case = '/genomics/projects/B061-PanicosShangaris/nana/gdm/results_v3/cellphonedb/case/cpdb_stat_case.pkl'
stat_out_path_control = '/genomics/projects/B061-PanicosShangaris/nana/gdm/results_v3/cellphonedb/control/cpdb_stat_control.pkl'

# For CASE
meta_file_path_case = '/genomics/projects/B061-PanicosShangaris/nana/data/cellphonedb/cd4_case_metadata.txt'
counts_file_path_case = '/genomics/projects/B061-PanicosShangaris/nana/data/cellphonedb/cd4_case_expression_matrix.txt'
out_path_case = '/genomics/projects/B061-PanicosShangaris/nana/gdm/results_v3/cellphonedb/case'

from cellphonedb.src.core.methods import cpdb_degs_analysis_method

cpdb_results_case = cpdb_degs_analysis_method.call(
    cpdb_file_path=cpdb_file_path,                             # CellphoneDB database zip file
    meta_file_path=meta_file_path_case,                        # Metadata file for CASE
    counts_file_path=counts_file_path_case,                    # Expression matrix for CASE
    degs_file_path=deg_file_path,                              # DEG file (same for both CASE and CONTROL)
    counts_data='hgnc_symbol',                                 # Defines the gene annotation in counts matrix
    score_interactions=True,                                   # Whether to score interactions or not
    threshold=0.1,                                             # Minimum % of cells expressing a gene for analysis
    result_precision=3,                                        # Rounding for mean values in results
    separator='|',                                             # Separator for cells in result dataframes
    debug=False,                                               # Save intermediate tables in pkl format
    output_path=out_path_case,                                 # Output path for CASE
    output_suffix=None,                                        # Optional output suffix
    threads=25                                                 # Number of threads
)

# Save cpdb_results_case to a pickle file
with open(pickle_out_path_case, 'wb') as f:
    pickle.dump(cpdb_results_case, f)

# For CONTROL
meta_file_path_control = '/genomics/projects/B061-PanicosShangaris/nana/data/cellphonedb/cd4_control_metadata.txt'
counts_file_path_control = '/genomics/projects/B061-PanicosShangaris/nana/data/cellphonedb/cd4_control_expression_matrix.txt'
out_path_control = '/genomics/projects/B061-PanicosShangaris/nana/gdm/results_v3/cellphonedb/control'

cpdb_results_control = cpdb_degs_analysis_method.call(
    cpdb_file_path=cpdb_file_path,                             # CellphoneDB database zip file
    meta_file_path=meta_file_path_control,                     # Metadata file for CONTROL
    counts_file_path=counts_file_path_control,                 # Expression matrix for CONTROL
    degs_file_path=deg_file_path,                              # DEG file (same for both CASE and CONTROL)
    counts_data='hgnc_symbol',                                 # Defines the gene annotation in counts matrix
    score_interactions=True,                                   # Whether to score interactions or not
    threshold=0.1,                                             # Minimum % of cells expressing a gene for analysis
    result_precision=3,                                        # Rounding for mean values in results
    separator='|',                                             # Separator for cells in result dataframes
    debug=False,                                               # Save intermediate tables in pkl format
    output_path=out_path_control,                              # Output path for CONTROL
    output_suffix=None,                                        # Optional output suffix
    threads=25                                                 # Number of threads
)

# Save cpdb_results_control to a pickle file
with open(pickle_out_path_control, 'wb') as f:
    pickle.dump(cpdb_results_control, f)


# ---- No significant interactions according to the plotting function, but this is false and likely not related to v5.

# ---- Stat method tried below but no longer needed as DEGs are more informative.

cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                 # mandatory: CellphoneDB database zip file.
    meta_file_path = meta_file_path_case,                 # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path_case,             # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
    active_tfs_file_path = None,           # optional: defines cell types and their active TFs.
    microenvs_file_path = None,       # optional (default: None): defines cells per microenvironment.
    score_interactions = True,                       # optional: whether to score interactions or not. 
    iterations = 1000,                               # denotes the number of shufflings performed in the analysis.
    threshold = 0.1,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
    threads = 5,                                     # number of threads to use in the analysis.
    debug_seed = 42,                                 # debug randome seed. To disable >=0.
    result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
    pvalue = 0.05,                                   # P-value threshold to employ for significance.
    subsampling = False,                             # To enable subsampling the data (geometri sketching).
    subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
    subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
    subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
    separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
    output_path = '/genomics/projects/B061-PanicosShangaris/nana/gdm/results_v3/cellphonedb_stat',                          # Path to save results.
    output_suffix = None                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
    )
    
  # Save cpdb_results_control to a pickle file
with open(stat_out_path_case, 'wb') as f:
    pickle.dump(cpdb_results, f)

cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path = cpdb_file_path,                 # mandatory: CellphoneDB database zip file.
    meta_file_path = meta_file_path_control,                 # mandatory: tsv file defining barcodes to cell label.
    counts_file_path = counts_file_path_control,             # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
    counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
    active_tfs_file_path = None,           # optional: defines cell types and their active TFs.
    microenvs_file_path = None,       # optional (default: None): defines cells per microenvironment.
    score_interactions = True,                       # optional: whether to score interactions or not. 
    iterations = 1000,                               # denotes the number of shufflings performed in the analysis.
    threshold = 0.1,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
    threads = 5,                                     # number of threads to use in the analysis.
    debug_seed = 42,                                 # debug randome seed. To disable >=0.
    result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
    pvalue = 0.05,                                   # P-value threshold to employ for significance.
    subsampling = False,                             # To enable subsampling the data (geometri sketching).
    subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
    subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
    subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
    separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
    debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
    output_path = '/genomics/projects/B061-PanicosShangaris/nana/gdm/results_v3/cellphonedb_stat/control',                          # Path to save results.
    output_suffix = None                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
    )
    
  # Save cpdb_results_control to a pickle file
with open(stat_out_path_control, 'wb') as f:
    pickle.dump(cpdb_results, f)

