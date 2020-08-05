#!/usr/bin/env Rscript

# Example of a call:
# Rscript ra_chris_wallace/Data_prep/matrix_transposer.R -a TRUE
# -d ./MDI/Data/Original\ data/ -w ./VSN_NA_min_data/ -n TRUE
# --na_people_threshold 0.1 --na_probe_threshold 0.1 -t TRUE -v TRUE


# Rscript to convert data from CEDAR cohorts (found here: http://cedar-web.giga.ulg.ac.be/)
# Call from the command line with input arguments
# Writes to current directory

# === Libraries ================================================================

# MDI specific codes (not sure if still used)
# source("/home/MINTS/sdc56/Desktop/MDI/mdipp-1.0.1/scripts/analysis.R") # install.packages("mcclust", dep = T)

# For posterioir similarity matrices
# Rcpp::sourceCpp("/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/posterior_sim_mat.cpp") # install.packages("Rcpp", dep = T)



# For tibbles
library(tibble, quietly = T) # for dataframe of lists

# For data wrangling
library(dplyr, quietly = T)

# Heatmapping
library(pheatmap, quietly = T) # install.packages("pheatmap", dep = T)

# Colour palettes
library(RColorBrewer, quietly = T)  # for pheatmap values
library(viridis, quietly = T)       # for annotation labels

# For Bayesian clustering output (possibly should use mcclust) and Adjusted Rand Index
library(mclust, quietly = T) # install.packages("mclust", dep = T)

# For plotting
library(ggplot2, quietly = T)

# For command line arguments
library(optparse, quietly = T) # install.packages("optparse")

# Load data.table to access fread and fwrite functions
library(data.table, quietly = T) # install.packages("data.table", dep = T)

# Load magrittr for the pipe %>%
library(magrittr, quietly = T)

# For defensive programming in plot_fused_samples
library(attempt, quietly = T)

# For having several plots on the same grid
library(cowplot, quietly = T)

# For praise
library(praise, quietly = T)

# My own helper function
library(mdiHelpR)

# function_dir <- "~/Documents/PhD/Year_1/Consensus_clustering/ra_chris_wallace/Analysis/Analysis_script_functions/"
#
# function_scripts <- c(
#   "plot_phi_series.R",
#   "plot_phi_densities.R",
#   "plot_phi_histograms.R",
#   "plot_similarity_matrices.R",
#   "plot_rand_index.R",
#   "plot_fused_samples.R",
#   "plot_comparison_expression_clustering.R",
#   "plot_mass_parameters.R",
#   "plot_clusters_series.R",
#   "create_psm.R",
#   "arandi_matrices.R"
# )
#
# for (f in paste0(function_dir, function_scripts)) {
#   source(f)
# }

# === Functions ================================================================

# print("Functions")

# Function to find the mode of a vector
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Function to plot and save tree for a given similarity matrix
plot_tree <- function(sim_mat, plot_name, plot_type = ".png") {
  if (!(plot_type %in% c(".pdf", ".png"))) {
    stop("Wrong plot type attempted: must be one of '.pdf' or '.png'")
  }
  
  # Open a file
  if (plot_type == ".png") {
    png(plot_name)
  } else {
    pdf(plot_name)
  }
  # 2. Create a plot
  sim_mat %>%
    dist() %>%
    hclust() %>%
    plot()
  
  # Close the pdf file
  dev.off()
}

# For creating unique gene names
unique_names_from_recurring <- function(name_vec) {
  
  # Number of entries
  n <- length(name_vec)
  
  # Duplicate the name vec - this will hold the transformed versions
  duplicate_vec <- name_vec
  
  for (i in 1:n) {
    n_occurences <- sum(name_vec == name_vec[i])
    
    if (n_occurences > 1) {
      duplicate_vec[i] <- paste0(
        duplicate_vec[i],
        ".",
        sum(name_vec[1:i] == name_vec[i])
      )
    }
  }
  duplicate_vec
}

# User inputs from command line
input_arguments <- function() {
  option_list <- list(
    
    # Convert all files in target destination (default is FALSE)
    optparse::make_option(c("--datasets"),
                          type = "character",
                          default = NA,
                          help = "names of dataset files [default= %default]",
                          metavar = "character"
    ),
    
    # Directory to read from
    optparse::make_option(c("-d", "--dir"),
                          type = "character",
                          default = ".",
                          help = "directory to read files from [default= %default]",
                          metavar = "character"
    ),
    
    # Number of genes in each dataset
    optparse::make_option(c("--n_samples"),
                          type = "integer",
                          default = NA,
                          help = "Number of samples in the datasets (if NA is given this will be calculated - this is included to avoid repetitive calculations) [default= %default]",
                          metavar = "integer"
    ),
    
    # Number of iterations in MCMC step of MDI
    optparse::make_option(c("-n", "--n_iter"),
                          type = "integer",
                          default = NA,
                          help = "Number of iterations in the MCMC step of MDI (if NA is given will calcualte as the number of rows in the MDI output file) [default= %default]",
                          metavar = "integer"
    ),
    
    # Thinning factor in MCMC step of MDI
    optparse::make_option(c("-t", "--thin"),
                          type = "integer",
                          default = 1,
                          help = "Thinnging factor in the MCMC step of MDI [default= %default]",
                          metavar = "integer"
    ),
    
    # Burn-in for MCMC step of MDI
    optparse::make_option(c("-b", "--burn"),
                          type = "integer",
                          default = 0,
                          help = "Burn-in for the MCMC step of MDI [default= %default]",
                          metavar = "integer"
    ),
    
    # .csv file describing which probes are not present in which datasets
    optparse::make_option(c("--samples_present"),
                          type = "character",
                          default = NA, # "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/samples_present_per_dataset.csv",
                          help = "Name of .csv files describing which samples are missing in which dataset [default= %default]",
                          metavar = "character"
    ),
    
    ### MOVE THIS TO A PRE-PROCESSING STAGE
    # .csv connecting probe IDs to genes
    optparse::make_option(c("--probe_key"),
                          type = "character",
                          default = NA, # "/home/MINTS/sdc56/Desktop/ra_chris_wallace/Analysis/probe_key_unique_names.csv",
                          help = "Name of .csv files describing which probes are associated with which gene [default= %default]",
                          metavar = "character"
    ),
    
    # Plot types to generate
    optparse::make_option(c("--plot_type"),
                          type = "character",
                          default = ".png",
                          help = "Plot types, one of ``.png'' or ``.pdf'' [default= %default]",
                          metavar = "character"
    ),
    
    # Instruction to plot trees
    optparse::make_option(c("--plot_trees"),
                          type = "logical",
                          default = FALSE,
                          help = "Instruction to plot dendrogram from clustering on PSMs [default= %default]",
                          metavar = "logical"
    ),
    
    # Instruction to plot adjusted rand index
    optparse::make_option(c("--plot_rand_index"),
                          type = "logical",
                          default = TRUE,
                          help = "Instruction to plot adjusted rand index comparing clusterings [default= %default]",
                          metavar = "logical"
    ),
    
    # Instruction to plot trees
    optparse::make_option(c("--plot_similarity_matrices"),
                          type = "logical",
                          default = TRUE,
                          help = "Instruction to plot heatmaps of the clustering over iterations [default= %default]",
                          metavar = "logical"
    ),
    
    # Instruction to plot phi scatter plots
    optparse::make_option(c("--plot_phi_series"),
                          type = "logical",
                          default = TRUE,
                          help = "Instruction to plot phi parameters over MCMC iterations for each dataset [default= %default]",
                          metavar = "logical"
    ),
    
    # Instruction to plot phi densities
    optparse::make_option(c("--plot_phi_densities"),
                          type = "logical",
                          default = TRUE,
                          help = "Instruction to plot density of phi parameters for each dataset [default= %default]",
                          metavar = "logical"
    ),
    
    # Instruction to plot phi densities
    optparse::make_option(c("--plot_phi_histograms"),
                          type = "logical",
                          default = TRUE,
                          help = "Instruction to plot histogram of phi parameters for each dataset [default= %default]",
                          metavar = "logical"
    ),
    
    # Instruction to plot phi densities
    optparse::make_option(c("--plot_mass_parameters"),
                          type = "logical",
                          default = TRUE,
                          help = "Instruction to plot mass parameter values across iterations for each dataset [default= %default]",
                          metavar = "logical"
    ),
    
    # Location for raw data
    optparse::make_option(c("--raw_data_dir"),
                          type = "character",
                          default = NA,
                          help = "Path to the directory containing the raw data .csv files (if NA no expression heatmaps can be generated) [default= %default]",
                          metavar = "character"
    ),
    
    
    # Instruction to plot the raw data
    optparse::make_option(c("--plot_raw_data"),
                          type = "logical",
                          default = TRUE,
                          help = "Instruction to plot heatmaps of the raw data [default= %default]",
                          metavar = "logical"
    ),
    
    # Set the random seed
    optparse::make_option(c("--seed"),
                          type = "integer",
                          default = 1,
                          help = "Random seed for script [default= %default]",
                          metavar = "integer"
    ),
    
    
    # Threshold for proportion of iterations we require samples to have the same
    # label to be considered ''fused''
    optparse::make_option(c("--fusion_threshold"),
                          type = "double",
                          default = 0.5,
                          help = "Threshold for proportion of iterations we require samples to have the same label to be considered ''fused'' [default= %default]",
                          metavar = "double"
    ),
    
    
    # Instruction to plot raw data heatmaps for fused samples
    optparse::make_option(c("--plot_fused_samples"),
                          type = "logical",
                          default = TRUE,
                          help = "Instruction to plot heatmaps of the expression data for fused genes [default= %default]",
                          metavar = "logical"
    ),
    
    # Instruction to plot heatmaps of the expression data side-by-side with the similarity matrix
    optparse::make_option(c("--plot_comparison"),
                          type = "logical",
                          default = TRUE,
                          help = "Instruction to plot heatmaps of the expression data side-by-side with the similarity matrix [default= %default]",
                          metavar = "logical"
    ),
    
    # Instruction to plot heatmap of the average adjusted rand index between datasets
    optparse::make_option(c("--plot_arandi_matrix"),
                          type = "logical",
                          default = TRUE,
                          help = "Instruction to plot heatmap of the average adjusted rand index between datasets [default= %default]",
                          metavar = "logical"
    ),
    
    # Instruction to plot the number of clusters present per iteration
    optparse::make_option(c("--plot_clusters"),
                          type = "logical",
                          default = TRUE,
                          help = "Instruction to plot the number of clusters present per iteration [default= %default]",
                          metavar = "logical"
    ),
    
    # Instruction to time programme
    optparse::make_option(c("--time"),
                          type = "logical",
                          default = FALSE,
                          help = "instruciton to record runtime of function [default= %default]",
                          metavar = "logical"
    )
  )
  opt_parser <- optparse::OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
}


# === Main script ==============================================================

# Set the ggplot2 theme
theme_set(theme_bw())

no_missingness <- TRUE

# print("In script")

args <- input_arguments()
save_plots <- T
stm_i <- Sys.time()

seed <- args$seed
set.seed(seed)

# print("Inputs received?")

# Colour palette for PSMs
col_pal <- colorRampPalette(c("white", "#146EB4"))(100)
breaks_0_1 <- c(
  seq(0, 1, length.out = ceiling(length(col_pal)) + 1)
)

# If there are missing entries than we include -1 in the PSM to indicate these
if(no_missingness){
  col_pal_sim <- col_pal
  my_breaks <- breaks_0_1
} else {
  col_pal_sim <- colorRampPalette(c("#FF9900", "white", "#146EB4"))(100)
  
  my_breaks <- c(
    seq(-1, 0, length.out = ceiling(palette_length_expr / 2) + 1),
    seq(1 / palette_length_expr, 1, length.out = floor(palette_length_expr / 2))
  )
}

col_pal_expr <- colorRampPalette(c("#146EB4", "white", "#FF9900"))(100)
palette_length_expr <- length(col_pal_expr)

# All possible datasets (and the names of the probes present columns)
# all_datasets <- c(
#   "CD14",
#   "CD15",
#   "CD19",
#   "CD4",
#   "CD8",
#   "IL",
#   "PLA",
#   "RE",
#   "TR"
# )

# Directory holding the raw data files
data_dir <- args$raw_data_dir
do_no_raw_data_plots <- F

if (is.na(data_dir)) {
  do_no_raw_data_plots <- T
}

# print("Me")

# Read in the file relating the probe IDs to the related gene
# SHOULD THIS BE PRE-PROCESSING?
# probe_key <- fread(args$probe_key)

# Read in the MDI output file
file_path <- args$dir
# file_path <- "Analysis/MDI_runs/vsn_many_seeds/"

# MDI call specific values
n_iter <- args$n_iter
thin <- args$thin
burn <- args$burn

# Number of samples present in datasets
# The number of samples is the number of columns in the mcmc_output
# exclusing the columns containing information on the phi and mass parameters
# divided by the number of datasets present
# (i.e. (ncol(MDI_OUTPUT) - n_datasets - choose(n_datasets, 2)) / n_datasets)
n_samples <- args$n_samples

# The files / tissues used in the MDI
# Have to have the order from the call
# CD4, CD8, CVD19, CD14, CD15, platelets, ileonic, colonic and rectal biopsies
# transcriptome data for six circulating immune cell types (CD4+ T lymphocytes,
# CD8+ T lymphocytes, CD19+ B lymphocytes, CD14+ monocytes, CD15+ granulocytes,
# platelets) as well as ileal, colonic, and rectal biopsies (IL, TR, RE)
# 323 healthy Europeans
files_present <- args$datasets %>%
  strsplit(., " ") %>%
  unlist()

n_datasets <- length(files_present)
col_names <- paste0("D", 1:n_datasets)

# Remove the file extension
dataset_names <- tools::file_path_sans_ext(files_present)

# Instructions for different plots and calculations
do_dendrograms_ie_trees <- args$plot_trees
do_rand_plot <- args$plot_rand_index
do_similarity_matrices_plot <- args$plot_similarity_matrices
do_phis_series <- args$plot_phi_series
do_phis_densities <- args$plot_phi_densities
do_phis_histograms <- args$plot_phi_histograms
do_mass_parameter_plots <- args$plot_mass_parameters
do_arandi_matrices <- args$plot_arandi_matrix
do_clusters_series <- args$plot_clusters

# Plot based upon raw data - do not allow to be TRUE if no raw data is given
if (do_no_raw_data_plots) {
  cat(paste0(
    "\nNo directory given for raw data.\nPlots incorporating raw data ",
    "are not possible regardless of user instruction. \nDefaulting to FALSE.\n."
  ))
  do_raw_data_heatmap <- FALSE
  do_fused_samples_heatmaps <- FALSE
  do_comparison_plots <- FALSE
} else {
  do_raw_data_heatmap <- args$plot_raw_data
  do_fused_samples_heatmaps <- args$plot_fused_samples
  do_comparison_plots <- args$plot_comparison
}

# Plot type
plot_type <- args$plot_type

# Fusion threshold
fusion_threshold <- args$fusion_threshold

# Create the common save path
save_path <- file_path

# For Rand index plots
# Save the common part of each title
generic_title <- "MDI: Adjusted Rand index for"
generic_save_name <- file_path


mdi_output_files <- list.files(
  path = file_path,
  full.names = T,
  include.dirs = F,
  pattern = ".csv$"
) # %>%
# grep("csv", ., value = TRUE)

num_files <- length(mdi_output_files)
file_names <- basename(tools::file_path_sans_ext(mdi_output_files))

# print(mdi_output_files)

# For plotting phis
phis <- list()
count <- 0
start_index <- ceiling(burn / thin) + 1 # Consider letting this equal "burn"

# Convert this into useable output using functions provided by Sam Mason
mcmc_out_lst <- list() # vector("list", num_files)

# Large palette of qualitative colours
# (see https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r)
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
col_vector <- brewer.pal(n = 12, name = "Paired") # unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# print("Hi")
for (i in 1:num_files) {
  curr_name <- file_names[[i]]
  mcmc_out_lst[[curr_name]] <- fread(mdi_output_files[[i]])
}
# print("ho")
if (is.na(n_iter)) {
  n_iter <- nrow(mcmc_out_lst[[1]]) * thin
}

# The effective number of iterations saved
eff_n_iter <- n_iter / thin # - burn

# Note that as each mdi output is on the same data, we only need to count one dataset
if (is.na(n_samples)) {
  dataset_1_ind <- colnames(mcmc_out_lst[[1]]) %>%
    grep("Dataset1", .)
  
  sample_names <- colnames(mcmc_out_lst[[1]])[dataset_1_ind] %>%
    stringr::str_remove_all("Dataset1_")
  # %>%
  # stringr::str_remove_all("X")
  
  n_samples <- sample_names %>%
    length()
}


# Boolean instructing labels to be included in heatmaps
show_heatmap_labels <- TRUE

if (n_samples > 50) {
  show_heatmap_labels <- FALSE
}

# Read in the probes present - this is a matrix of bools with the first column
# as the probe IDs and the remaining columns corrresponding to the cell types
# with TRUE indicating the probe is present in this cell type (i.e. not added
# manually with an imputed value) and FALSE indicates we added it in.
samples_present <- args$samples_present
check_missingness <- !is.na(samples_present)

# If we have not given the meta-data to check this, assume all samples are
# present in all datasets
if (check_missingness) {
  samples_present_dt <- fread(samples_present)
  all_datasets <- colnames(samples_present_dt)[2:ncol(samples_present_dt)]
} else {
  samples_present_dt <- matrix(data = T, ncol = n_datasets, nrow = n_samples) %>%
    as.data.frame() %>%
    set_rownames(sample_names) %>%
    set_colnames(dataset_names)
  
  all_datasets <- dataset_names
}

# If using datasets with missing entries, keep only the meta data for the relevant files
cols_to_keep <- all_datasets %in% dataset_names

# curr_viable_datasets <- list(
#   c(
#     "CD14",
#     "CD15",
#     "CD19",
#     "CD4",
#     "CD8",
#     "IL",
#     "PLA",
#     "RE",
#     "TR"
#   ),
#   c(
#     "MDItestdata1",
#     "MDItestdata2",
#     "MDItestdata3",
#     "MDItestdata4",
#     "MDItestdata5",
#     "MDItestdata6"
#   )
# )
#
# # Check columns are as expected.
# if(length(all_datasets) == length(curr_viable_datasets[[1]])){
#   if(! all.equal(all_datasets[match(curr_viable_datasets[[1]], all_datasets)], curr_viable_datasets[[1]])){
#     cat("\nDatasets appears to be CEDAR (based on samples_present_per_dataset.csv), but datasets in columns do not match expected.")
#     stop(paste("Expected", curr_viable_datasets[[1]]))
#
#   }
# }
#
# if(length(all_datasets) == length(curr_viable_datasets[[2]])){
#   if(! all.equal(all_datasets[match(curr_viable_datasets[[2]], all_datasets)], curr_viable_datasets[[2]])){
#     cat("\nDatasets appears to be Yeast (based on samples_present_per_dataset.csv), but datasets in columns do not match expected.")
#     stop(paste("Expected", curr_viable_datasets[[2]]))
#   }
# }

mdi_input_files <- list.files(
  path = data_dir,
  full.names = T,
  include.dirs = F,
  pattern = ".csv$"
) # %>%
# grep("csv", ., value = TRUE)

lower_files <- stringr::str_to_lower(files_present) %>%
  tools::file_path_sans_ext()

mdi_input_files_lower <- stringr::str_to_lower(mdi_input_files) %>%
  tools::file_path_sans_ext() %>%
  basename()

expr_data_order <- match(lower_files, mdi_input_files_lower)

# === Plotting phis ==========================================================

# Plot phi value between each dataset combination across iterations
if (do_phis_series) {
  if (n_datasets > 1) {
    # print(str(mcmc_out_lst))
    # print(file_path)
    
    cat("\n\nPlotting phi values across iterations.\n")
    mdiHelpR::plotPhiSeries(mcmc_out_lst,
                            file_path,
                            num_files,
                            n_datasets,
                            start_index,
                            eff_n_iter,
                            save_plots = T,
                            burn = burn,
                            thin = thin
    )
  } else {
    cat("\nOnly one dataset. No phi parameters.")
  }
}


# === Prepare the tibble =======================================================

# Create an empty dataframe with column names corresponding to dataset numbers
compare_df <- data.frame(matrix(ncol = n_datasets, nrow = n_samples))

# Set the cell type to the column names (make sure order is as per MDI call)
colnames(compare_df) <- dataset_names

cat("\nConstructing tibble.\n")

# Now put everything in a tibble
compare_tibble <- tibble(
  mdi = rep(file_names, n_datasets),
  dataset = rep(dataset_names, num_files), # unlist(lapply(dataset_names, rep, num_files))
  seed = unlist(lapply(1:num_files, rep, n_datasets)), # rep(1:num_files, num_datasets),
  phis = rep(vector("list", num_files), n_datasets),
  pred_allocation = list(compare_df),
  mdi_allocation = rep(vector("list", num_files), n_datasets),
  n_clust = rep(vector("list", num_files), n_datasets),
  similarity_matrix = list(data.frame(matrix(ncol = n_samples, nrow = n_samples))),
  correlation_matrix = list(data.frame(matrix(ncol = n_samples, nrow = n_samples))),
  expression_data = rep(vector("list", num_files), n_datasets),
  non_zero_probes_ind = rep(vector("list", num_files), n_datasets),
  non_zero_probes = rep(vector("list", num_files), n_datasets),
  fused_probes = rep(vector("list", num_files), n_datasets),
  fused_non_zero_probes = rep(vector("list", num_files), n_datasets),
  mass_parameter = rep(vector("list", num_files), n_datasets)
)

# compare_tibble$phis <- phis

# === MDI output ===============================================================

# Global similarity between datasets
phis <- list()

# Collect the number of clusters present in each iteration
n_clust_list <- vector("list", length = n_datasets) # list()

for (j in 1:num_files) {
  
  # mdi_pos_sim_mat <- list()
  
  # Capture the allocation information in the named lists and the predicted
  # allocation in the dataframe
  phis[[j]] <- mcmc_out_lst[[j]] %>%
    dplyr::select(contains("Phi"))
  
  
  for (i in 1:n_datasets) {
    dataset_name <- paste0("Dataset", i)
    
    # Extract mass parameters
    curr_mass_parameter <- "MassParameter_" %>% paste0(i)
    .mass_parameters <- mcmc_out_lst[[j]][[curr_mass_parameter]]
    
    # Get the allocation and drop the burn in
    .mdi_alloc <- mcmc_out_lst[[j]] %>%
      dplyr::select(contains(dataset_name)) %>%
      magrittr::extract(start_index:eff_n_iter, )
    
    # cat(paste("\n\n", dataset_name, " currently\n\n"))
    # print(str(.mdi_alloc))
    
    # print(nrow(.mdi_alloc))
    # print(ncol(.mdi_alloc))
    
    # If there is the same number of clusters in each iteration apply(., 1, unique)
    # retrns a matrix rather than a list with each column representing the unique
    # entries in one iteration (or in a given row of .mdi_alloc).
    # Due to this we must check if we have a list or not and treat the object
    # appropriately.
    unique_labels <- .mdi_alloc %>%
      apply(1, unique)
    
    # Find the number of clusters in each iteration
    if (typeof(unique_labels) == "list") {
      n_clust_list[[i]] <- unique_labels %>%
        lapply(length) %>%
        unlist()
    } else {
      if (!is.null(nrow(unique_labels))) {
        n_clust_list[[i]] <- unique_labels %>%
          nrow() %>%
          rep((eff_n_iter - (start_index - 1))) # as we include start_index
      } else {
        n_clust_list[[i]] <- rep(1, (eff_n_iter - (start_index - 1)))
      }
    }
    
    
    
    # if(i > 1){
    # print(str(unique_labels))
    # print(str(n_clust_list[[i]]))
    # }
    # This step does not need to be part of the for loop
    # if (i == 1) {
    #
    #   # Find the Probe IDs
    #   probe_names <- colnames(.mdi_alloc) %>%
    #     stringr::str_remove_all(paste0(dataset_name, "_")) %>%
    #     stringr::str_remove_all("X")
    #
    #   # Find the sample names
    #   sample_names <- colnames(.mdi_alloc) %>%
    #       stringr::str_remove_all(paste0(dataset_name, "_")) %>%
    #       stringr::str_remove_all("X")
    #
    #   # I think this should be a pre-processing step now.
    #   if(FALSE){
    #   # Find the relevant part of the key
    #   probe_key_rel <- probe_key[probe_key$ProbeID %in% probe_names, ]
    #
    #   # Pull out the Gene IDs in the correct order
    #   gene_id <- probe_key_rel %>%
    #     .[match(probe_names, .$ProbeID)] %>%
    #     .$Unique_gene_name
    #
    #   }
    #
    #   samples_present_dt <- samples_present_dt[samples_present_dt$V1 %in% probe_names, ]
    #
    #   # Add the gene names to the samples_present dataframe
    #   samples_present_dt$Gene_names <- gene_id # [match(probe_key_rel$ProbeID, samples_present_dt$V1)]
    #
    #   # Make sure the order is as in the MDI data
    #   # samples_present_dt <- samples_present_dt[match(probe_names, samples_present_dt$V1)]
    #
    #   # unique_gene_id <- gene_id
    # }
    
    curr_dataset <- dataset_names[[i]]
    
    # Create and save the posterior similarity matrix (PSM) for the current
    # allocation
    # We transpose as we are interested in how the genes cluster rahter than the
    # people
    # .sim_mat <- make_psm(.mdi_alloc) %>%
    #   set_colnames(gene_id) %>%
    #   set_rownames(gene_id)
    
    .sim_mat <- mdiHelpR::makePSM(.mdi_alloc) %>%
      set_colnames(sample_names) %>%
      set_rownames(sample_names)
    
    # Use the maxpear() function from mcclust to interpret the PSM as a clustering
    .pred_alloc <- .sim_mat %>%
      mcclust::maxpear()
    
    if (check_missingness) {
      # extract the empty genes
      rel_empty_samples <- samples_present_dt %>%
        dplyr::select(dplyr::one_of(c(curr_dataset, "Gene_names")))
      
      # Create a diagonal matrix of genes x genes with -2 on the diagonal entries
      # corresponding to empty genes and 0's elsewhere
      empty_gene_mat <- diag(rel_empty_samples[[curr_dataset]]) * (-2)
      
      # Now highlight the empty probes by setting their diagonal to -1 in the PSM
      .sim_mat <- .sim_mat %>%
        add(empty_gene_mat)
    }
    
    
    
    # Set the row nad column names of the PSM
    # row.names(.sim_mat) <- colnames(.sim_mat) <- gene_id
    
    # Set the column names of the MDI to the Probe IDs (consider using gene IDs)
    # colnames(.mdi_alloc) <- probe_names
    colnames(.mdi_alloc) <- sample_names
    
    # Save these objects to the tibble
    compare_tibble$mdi_allocation[i + (j - 1) * num_files][[1]] <- .mdi_alloc
    compare_tibble$similarity_matrix[i + (j - 1) * num_files][[1]] <- .sim_mat
    
    # Record the numebr of clusters per iteration
    compare_tibble$n_clust[i + (j - 1) * num_files][[1]] <- n_clust_list[[i]]
    
    # Record this in the tibble
    compare_tibble$pred_allocation[i + (j - 1) * num_files][[1]] <- .pred_alloc$cl
    
    # Record mass parameters
    compare_tibble$mass_parameter[i + (j - 1) * num_files][[1]] <- .mass_parameters
  }
}

# === Probes present ===========================================================

# print("Finding probes present.")

# # Find which probes are relevant from the full set
# probes_actually_present_ind <- samples_present_dt %>%
#   magrittr::use_series("V1") %>%
#   magrittr::is_in(probe_names)
#
# # Select these
# probes_actually_present <- samples_present_dt[probes_actually_present_ind, ]
#
# # Find the appropriate order
# probes_order <- match(probe_names, probes_actually_present$V1)
#
# # Order the probes so comparable to allocation data frame
# samples_present_ordered <- probes_actually_present[probes_order, ] %>%
#   as.data.frame() %>%
#   set_colnames(all_datasets)
#
# # Remove the irrelevant columns and set row names
# samples_present_final <- samples_present_ordered %>%
#   magrittr::set_rownames(probes_actually_present$V1) %>%
#   magrittr::extract(, cols_to_keep)
#
# colnames(samples_present_final) <- dataset_names

# === Phi denisty plots ================================================================

# plotPhiDensities <- function(phis, file_path, start_index, eff_n_iter,
#          burn = 0,
#          thin = 1) {
#   loc_dir <- paste0(file_path, "Phi_density_plots/")
#   dir.create(loc_dir, showWarnings = FALSE)
#
#   for (i in 1:length(phis)) {
#     for (j in 1:ncol(phis[[i]])) {
#       curr_phi <- colnames(phis[[i]])[j]
#
#       density_plot_name <- paste0(loc_dir, curr_phi, "_density_plot", plot_type)
#
#       density_title <- paste(curr_phi, ": density plot")
#
#       # Find which datasets are used here
#       dataset_indices <- readr::parse_number(curr_phi)
#
#       # The following steps do not work if we have more than 9 datasets (which is not relevant to me)
#       # This extracts the first dataset index (i.e. from 67 takes 6)
#       dataset_index_1 <- dataset_indices[1] / 10 %>%
#         floor(.)
#
#       # This extracts the second dataset index
#       dataset_index_2 <- dataset_indices[1] %% (10)
#
#       dataset_1 <- dataset_names[dataset_index_1]
#       dataset_2 <- dataset_names[dataset_index_2]
#
#       density_subtitle <- paste(
#         "Iterations",
#         (burn + thin),
#         "through",
#         n_iter,
#         "for",
#         dataset_1,
#         "and",
#         dataset_2
#       )
#
#       ggplot2::ggplot(data = phis[[i]][start_index:eff_n_iter, ], ggplot2::aes_string(x = curr_phi)) +
#         ggplot2::geom_density() +
#         ggplot2::labs(
#           title = density_title,
#           subtitle = density_subtitle
#         )
#
#       ggplot2::ggsave(density_plot_name)
#     }
#   }
# }

if (do_phis_densities) {
  if (n_datasets > 1) {
    cat("\nSaving phi density plots.\n")
    
    mdiHelpR::plotPhiDensities(phis, file_path, start_index, eff_n_iter,
                               burn = burn,
                               thin = thin
    )
  } else {
    cat("\nOnly one dataset. No phi parameters.")
  }
}

# === Phi histograms ===========================================================

if (do_phis_histograms) {
  if (n_datasets > 1) {
    cat("\nSaving phi histogram plots.\n")
    plotPhiHistograms(phis, file_path, start_index, eff_n_iter,
                      burn = burn,
                      thin = thin
    )
  } else {
    cat("\nOnly one dataset. No phi parameters.")
  }
}

# === Plot posterior similarity matrices =======================================

# If making heatplots of the clusterings across iterations
if (do_similarity_matrices_plot) {
  cat("\nSaving heatmaps of PSMs.\n")
  
  plot_similarity_matrices(
    compare_tibble$similarity_matrix,
    # samples_present_final,
    dataset_names,
    num_files,
    n_datasets,
    file_path,
    col_pal = col_pal_sim,
    breaks = my_breaks,
    show_labels = show_heatmap_labels
  )
}

# === Plot clusters per iteration ==============================================

# Plot the number of clusters present per iteration
if (do_clusters_series) {
  cat("\nPlotting the number of clusters present per iteration\n")
  
  # n_clust_list <- compare_tibble$n_clust
  
  plot_clusters_present(n_clust_list,
                        dataset_names,
                        n_datasets,
                        start_index,
                        eff_n_iter,
                        thin,
                        file_path,
                        # gen_main_title = "Number of clusters present per iteration",
                        # gen_save_name = NULL,
                        plot_type = plot_type
  )
}


# === Plot adjusted rand index==================================================
# If instructed to make Rand index plots
if (do_rand_plot) {
  cat("\nSaving scatter plots of adjusted rand index comparing final clustering to clustering at each iteration.\n")
  plot_rand_index(
    compare_tibble$mdi_allocation,
    file_path,
    dataset_names,
    num_files,
    n_datasets,
    eff_n_iter,
    burn = burn,
    thin = thin
  )
}

# If instructed do Rand index matrices and the heatmap
if (do_arandi_matrices) {
  if (n_datasets > 1) {
    arandi_matrices <- makeARIMatrices(compare_tibble$mdi_allocation, n_datasets)
    # arandi_matrices <- arandi_matrices(compare_tibble$mdi_allocation, n_datasets)
    
    avg_arandi_matrix <- averageMatrix(arandi_matrices) %>%
      set_colnames(dataset_names) %>%
      set_rownames(dataset_names)
    
    # avg_arandi_matrix <- average_matrix(arandi_matrices) %>%
    #   set_colnames(dataset_names) %>%
    #   set_rownames(dataset_names)
    
    arandi_pheatmap_title <- "Heatmap comparing adjusted rand index across datasets"
    arandi_pheatmap_file_name <- paste0(save_path, "Arandi_heatmap", plot_type)
    
    pheatmap(avg_arandi_matrix,
             cluster_rows = F,
             cluster_cols = F,
             main = arandi_pheatmap_title,
             filename = arandi_pheatmap_file_name,
             color = col_pal # ,
             # breaks = breaks_0_1
    )
  } else {
    cat("\nOnly one dataset. No comparison between datasets possible.")
  }
}

# === Plot mass parameters =====================================================

if (do_mass_parameter_plots) {
  cat("\nPlotting mass parameters over iterations.\n")
  plot_mass_parameters(compare_tibble,
                       thin = thin,
                       file_path = file_path,
                       plot_type = plot_type
  )
}

# === Heatmap expression data ==================================================
# If no raw data given, do not attempt this section.
if (!do_no_raw_data_plots) {
  loc_dir <- paste0(save_path, "Expression_heatmaps/")
  
  if (do_raw_data_heatmap) {
    cat("\nSaving gene expression heatmaps.\n")
    
    
    dir.create(loc_dir, showWarnings = FALSE)
  }
  
  # Generic title and filename for pheatmap
  gen_ph_title <- ": heatmap of expression data"
  gen_ph_file_name <- paste0(loc_dir, "pheatmap_")
  
  # Find the expression data files
  mdi_input_files <- list.files(
    path = data_dir,
    full.names = T,
    include.dirs = F,
    pattern = ".csv$"
  ) # %>%
  # grep("csv", ., value = TRUE)
  
  input_file_names <- basename(tools::file_path_sans_ext(mdi_input_files))
  data_to_keep <- match(files_present, input_file_names)
  
  input_file_names <- input_file_names[data_to_keep]
  mdi_input_files <- mdi_input_files[data_to_keep]
  
  datasets_relevant <- input_file_names
  relevant_input_files <- mdi_input_files
  
  expression_datasets <- input_file_names # %>%
  # gsub("([^\\_]+)\\_.*", "\\1", .)
  # stringr::str_replace("_sma_mat_nv", "")
  
  expression_datasets <- dataset_names
  # stop()
  
  # datasets_relevant_indices <- files_present %>%
  #   tools::file_path_sans_ext() %>%
  #   match(expression_datasets)
  
  # print("Here's a print statement")
  # print(datasets_relevant_indices)
  # print(mdi_input_files)
  
  # datasets_relevant <- expression_datasets[datasets_relevant_indices]
  # relevant_input_files <- mdi_input_files[datasets_relevant_indices]
  
  # print("Dog")
  # print(expression_datasets)
  # print("God")
  # print(datasets_relevant_indices)
  # print("Soon")
  # print(datasets_relevant)
  # print("Hi")
  # print(expr_data_order)
  
  # datasets_relevant <- datasets_relevant[expr_data_order]
  
  # print("MAmmy")
  # print(datasets_relevant)
  # print(relevant_input_files)
  
  # relevant_input_files <- relevant_input_files[expr_data_order]
  # print(relevant_input_files)
  
  mega_df <- data.frame(matrix(nrow = n_samples, ncol = 0)) %>%
    magrittr::set_rownames(sample_names) # magrittr::set_rownames(gene_id)
  
  big_annotation <- data.frame(matrix(nrow = n_samples, ncol = n_datasets)) %>%
    magrittr::set_rownames(sample_names) %>% # magrittr::set_rownames(gene_id) %>%
    magrittr::set_colnames(datasets_relevant)
  
  n_total_clusters <- 0
  
  # Vector of columns in each expression dataset (controls showing column names in
  # heatmaps)
  n_people <- c()
  
  data_files <- list()
  
  for (i in 1:n_datasets) {
    curr_dataset <- datasets_relevant[[i]]
    # print(curr_dataset)
    
    file_name <- gen_ph_file_name %>%
      paste0(datasets_relevant[[i]], plot_type)
    
    ph_title <- datasets_relevant[[i]] %>%
      paste0(gen_ph_title)
    
    # Prepare the clustering information as an annotation row for pheatmap
    pred_clustering <- compare_tibble$pred_allocation[compare_tibble$dataset == curr_dataset][[1]] %>%
      as.factor() %>%
      as.data.frame() %>%
      magrittr::set_rownames(sample_names) %>% # magrittr::set_rownames(gene_id) %>%
      magrittr::set_colnames(c("Cluster"))
    
    # Read in the expression data
    f <- relevant_input_files[[i]]
    # print(f)
    
    expression_data <- fread(f, header = T)
    
    # Convert from probe ids to gene ids if necessary
    # if(sum(expression_data[[1]] %in% probe_names) > 0){
    if (sum(expression_data[[1]] %in% sample_names) > 0) {
      expression_data[, 1] <- sample_names # gene_id[match(expression_data[[1]], probe_names)]
    }
    
    if (!all.equal(unname(unlist(expression_data[[1]])), sample_names)) { # gene_id)){
      stop("Gene ids not matching in expression data.")
    }
    
    # Tidy (remove NAs and row name column) and convert to the appropriate format
    # for pheatmap
    expression_data_tidy <- expression_data %>%
      magrittr::extract(, -1) %>%
      as.matrix() %>%
      magrittr::set_rownames(sample_names) # gene_id)
    
    expression_data_tidy[is.na(expression_data_tidy)] <- 0
    
    # Add the expression data to our tibble
    # print(compare_tibble$dataset == curr_dataset)
    curr_ind <- which(compare_tibble$dataset == curr_dataset)
    # print(curr_ind)
    num_occurences_dataset <- length(compare_tibble$expression_data[compare_tibble$dataset == curr_dataset])
    for (k in 1:num_occurences_dataset) {
      # print(k)
      # print(compare_tibble$dataset[curr_ind][[k]])
      compare_tibble$expression_data[curr_ind][[k]] <- expression_data_tidy
      compare_tibble$non_zero_probes_ind[curr_ind][[k]] <- .non_zero_probes <- rowSums(expression_data_tidy) != 0
      compare_tibble$non_zero_probes[curr_ind][[k]] <- names(.non_zero_probes)[.non_zero_probes]
      compare_tibble$correlation_matrix[curr_ind][[k]] <- expression_data_tidy %>%
        t() %>%
        cor()
      
      # pheatmap::pheatmap(cor(t(expression_data_tidy)))
    }
    
    data_files[[i]] <- expression_data
    
    # Specify colors based on cluster labels
    cluster_labels <- levels(pred_clustering$Cluster)
    n_clusters <- length(cluster_labels)
    
    n_total_clusters <- max(n_clusters, n_total_clusters)
    
    expr_min <- min(expression_data_tidy)
    expr_max <- max(expression_data_tidy)
    
    # expr_breaks <- define_breaks(col_pal_expr, lb = expr_min, ub = expr_max) %>% unique()
    expr_breaks <- mdiHelpR::defineBreaks(col_pal_expr,
                                          lb = expr_min,
                                          ub = expr_max
    ) %>%
      unique()
    
    n_people <- c(n_people, ncol(expression_data_tidy))
    show_expr_col_names <- TRUE
    if (n_people[i] > 50) {
      show_expr_col_names <- FALSE
    }
    
    if (do_raw_data_heatmap) {
      if (n_clusters > 12) {
        cat("\nToo many clusters. Cannot include annotation row.\n")
        
        
        # Pheatmap
        expression_data_tidy %>%
          pheatmap(
            filename = file_name,
            main = ph_title,
            color = col_pal_expr,
            breaks = expr_breaks,
            show_rownames = show_heatmap_labels,
            show_colnames = show_expr_col_names
          )
      } else {
        # col_pal <- col_vector[1:n_clusters] %>%
        # sample(col_vector, n_clusters) %>%
        # magrittr::set_names(cluster_labels)
        
        # Use distinct, colour-blind friendly annotation colours
        col_pal <- viridis(n_clusters) %>% 
          magrittr::set_names(cluster_labels)
        
        annotation_colors <- list(Cluster = col_pal)
        
        row_order <- pred_clustering$Cluster %>% order()
        
        # Pheatmap
        expression_data_tidy[row_order, ] %>%
          pheatmap(
            filename = file_name,
            main = ph_title,
            cluster_rows = F,
            annotation_row = pred_clustering,
            annotation_colors = annotation_colors,
            color = col_pal_expr,
            breaks = expr_breaks,
            show_rownames = show_heatmap_labels,
            show_colnames = show_expr_col_names
          )
      }
    }
    
    # Record the clustering for the current dataset for annotation purposes
    big_annotation[[curr_dataset]] <- as.factor(pred_clustering$Cluster)
    
    # Record the expression data from the current dataset in the compound dataset
    mega_df <- mega_df %>%
      dplyr::bind_cols(as.data.frame(expression_data_tidy))
  }
  
  # row.names(mega_df) <- probe_names
  
  big_file_name <- gen_ph_file_name %>%
    paste0("all_datasets", plot_type)
  
  big_ph_title <- "All datasets" %>%
    paste0(gen_ph_title)
  
  mega_matrix <- mega_df %>%
    as.matrix() %>%
    magrittr::set_rownames(sample_names) # gene_id)
  
  expr_min <- min(mega_matrix)
  expr_max <- max(mega_matrix)
  
  # expr_breaks <- define_breaks(col_pal_expr, lb = expr_min, ub = expr_max) %>% unique()
  expr_breaks <- mdiHelpR::defineBreaks(col_pal_expr,
                                        lb = expr_min,
                                        ub = expr_max
  ) %>%
    unique()
  
  if (do_raw_data_heatmap) {
    if (TRUE) { # n_total_clusters > 20) {
      pheatmap(mega_matrix,
               filename = big_file_name,
               main = big_ph_title,
               color = col_pal_expr,
               breaks = expr_breaks,
               show_colnames = F,
               show_rownames = show_heatmap_labels
      )
    } else {
      col_pal <- sample(col_vector, n_total_clusters) %>%
        magrittr::set_names(cluster_labels)
      
      annotation_colors <- list(Cluster = col_pal)
      
      pheatmap(mega_matrix,
               filename = big_file_name,
               main = big_ph_title,
               annotation_row = big_annotation,
               annotation_colors = annotation_colors,
               color = col_pal_expr,
               breaks = expr_breaks,
               cluster_cols = F,
               show_colnames = F,
               show_rownames = show_heatmap_labels
      )
    }
  }
}
# === Fused probes =============================================================

# Find which probes are ''fused'' across datasets
# We save this as a named list to the tibble. Each entry in the list corresponds
# to a dataset and records the fused probes between the current dataset and the
# entry name
if (n_datasets > 1) {
  cat("\nFinding ''fused'' probes.\n")
  for (k in 1:num_files) {
    for (i in 1:n_datasets) {
      dataset_i <- dataset_names[i]
      
      fused_probes_curr <- vector("list", n_datasets) %>%
        magrittr::set_names(dataset_names)
      
      fused_non_zero_probes_curr <- vector("list", n_datasets) %>%
        magrittr::set_names(dataset_names)
      
      # fused_probes_curr[dataset_i] <- 1.0
      
      for (j in 1:n_datasets) {
        dataset_j <- dataset_names[j]
        
        count <- count + 1
        
        # Find the unnormalised count for the numebr of times each probe has the same label across iterations
        .fusion_count <- compare_tibble$mdi_allocation[compare_tibble$dataset == dataset_i][[k]] %>%
          magrittr::equals(compare_tibble$mdi_allocation[compare_tibble$dataset == dataset_j][[k]])
        
        # convert this to a propotion
        .fusion_prob <- (1 / nrow(.fusion_count)) * colSums(.fusion_count)
        
        # Record as ''fused'' those probes for which the proportion of times they
        # have a common labelling exceeds some user-defined threshold (defualt of 0.5)
        
        # compare_tibble$fused_probes[
        #   compare_tibble$dataset == dataset_i
        # ][[k]] <-
        
        # This is rather ugly - we are recording a named list associating each
        # dataset with each other one
        fused_probes_curr[[dataset_j]] <- .fused_probes_ind_i <- .fusion_prob > fusion_threshold
        
        fused_non_zero_probes_curr[[dataset_j]] <- .fused_probes_ind_i & compare_tibble$non_zero_probes_ind[compare_tibble$dataset == dataset_i][[k]]
        
        # compare_tibble$fused_non_zero_probes[
        #   compare_tibble$dataset == dataset_j
        # ][[k]] <- .fused_probes_ind_j & compare_tibble$non_zero_probes_ind[compare_tibble$dataset == dataset_j][[k]]
      }
      
      # Record the named lists within the tibble
      compare_tibble$fused_probes[
        compare_tibble$dataset == dataset_i
        ][[k]] <- fused_probes_curr
      
      compare_tibble$fused_non_zero_probes[
        compare_tibble$dataset == dataset_i
        ][[k]] <- fused_non_zero_probes_curr
    }
  }
  
  # === Gene expression for fused and unfused genes ==============================
  
  if (do_fused_samples_heatmaps) {
    cat("\nSaving heatmaps of pairwise fused gene expression across datasets.\n")
    
    # print(str(compare_tibble$fused_probes))
    # print(str(samples_present_dt))
    
    row.names(samples_present_dt) <- row.names(compare_tibble$expression_data[[1]])
    # print(str(samples_present_dt))
    
    fused_gene_heatmaps(
      compare_tibble$expression_data,
      compare_tibble$fused_probes,
      sample_names, # gene_id,
      dataset_names,
      file_path,
      n_datasets,
      plot_type,
      samples_present_dt # ,
      # show_row_labels = TRUE
    )
  }
} else {
  cat("\nOnly one dataset. No fusion possible.")
}
# === Compare similarity and gene expression ===================================

if (do_comparison_plots) {
  cat("\nSaving comparison of similarity matrix and expression data.\n")
  # plot_comparison_expression_to_clustering(
  mdiHelpR::plotComparisonData2Clustering(
    compare_tibble,
    dataset_names,
    n_datasets,
    file_path,
    plot_type,
    col_pal_sim = col_pal_sim,
    sim_breaks = my_breaks,
    show_row_labels = show_heatmap_labels,
    col_pal_expr = col_pal_expr,
    expr_breaks = expr_breaks
  )
  
  mdiHelpR::plotComparisonCorrelationPSMData(
    # plot_comparison_corr_sim_expr(
    compare_tibble,
    dataset_names,
    n_datasets,
    file_path,
    plot_type,
    col_pal_sim = col_pal_sim,
    sim_breaks = my_breaks,
    show_row_labels = show_heatmap_labels,
    col_pal_expr = col_pal_expr,
    expr_breaks = expr_breaks
  )
}

# === Timing ===================================================================

saveRDS(
  compare_tibble,
  paste0(save_path, "compare_tibble.rds")
)

if (args$time) {
  cat("\n")
  print(Sys.time() - stm_i)
  cat("\n")
}

# Well done
praise()