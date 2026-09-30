##########################################################################################################################

# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
# Created: December 2025; Last edited: September 2026
# In this script, we create an empirical and 5000 synthetic networks for each network generator algorithm. 

# IMPORTANT NOTES:
# 1. This script relies on a .sh script with identical name for execution. 
# 2. The script sources "00_constructNetworksFunction.R". Ensure you have this script in your working directory. 
# 3. The script requires path changing. The paths correspond to my HPC environment, NOT GitHub!
# 4. The script uses several default function arguments. Please check these. 
# 5. The script pre-defines parameter ranges for kappa, tau, and lambda for SENCA. If you wish to change these ranges,
#    this needs to be done manually. 
# 6. The script produces five replicates of each graph. This is computationally intensive and could be improved, but was 
#    a necessary initial addition to check if seeding was performed correctly and has not been altered. You may improve
#    this per your preferences. It should have no further effect on any results. 
# 7. The corresponding SLURM script saves output files in a designated folder. You can do this or choose not to, up to you. 
#    Ensure that you check this in the 01_NetworkGeneratorComparison.sh script (lines 3-4). 
# 8. You can also adapt the memory, CPUs and allocated time in the .sh script, primarily depending on the your array size
#    and the parameter grid size, as the latter requires a significant amount of memory. 

# For any issues, feel free to report a new issue on this GitHub repo or contact inez.derkx@swisstph.ch

##########################################################################################################################

### SET UP R ENVIRONMENT ###

# Empty list
rm(list = ls())

# I ran these scripts on a server and had issues with package versions.
# This ensures that all necessary packages used in the script were the ones I installed locally. 
.libPaths(c("/scicore/home/chitnis/derkx0000/R/ubuntu/4.2.1-foss-2022a", .libPaths()))
cat("Does path exist?", dir.exists("/scicore/home/chitnis/derkx0000/R/ubuntu/4.2.1-foss-2022a/igraph"), "\n")

# Load required libraries
library(igraph, lib.loc)
library(lhs, lib.loc)
library(dplyr, lib.loc)
library(randnet, lib.loc)
library(tidyr, lib.loc)
library(ggplot2, lib.loc)

# Set local root directory 
LOCAL_ROOT_DIR <- "/scicore/home/chitnis/derkx0000/GraphComparison/"
setwd(LOCAL_ROOT_DIR)

# Generate a new output folder to store all networks
ifelse(!dir.exists(file.path(LOCAL_ROOT_DIR, "Net_Sens")),
       dir.create(file.path(LOCAL_ROOT_DIR, "Net_Sens")), FALSE)

# Set path for storing output
out_path <- file.path(LOCAL_ROOT_DIR, "Net_Sens")

# Set path for storing graphs
out_path_graphs <- file.path(out_path, "Graphs")

##########################################################################################################################

### Load all helper functions (1-9)

# All functions are stored in this R script and corresponding .Rdata file: "HelperFunctions.RData"
source("00_constructNetworkFunctions.R")

##########################################################################################################################

### Execute code

## 1. Get the task_id and seed from SLURM environment
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])
base_seed <- 1000 + (task_id * 10000) # Give each task 10,000 "room"
max_seed = 5000 # adapt this where needed 

## 2. Load data  

# The edgelist is stored in ServerData.rds as a dataframe called 'dog_contact_zone_1'
# If you have your data stored differently or elsewhere, adapt the following:
ServerData <- readRDS("~/BaseData/Chad/ServerData.rds")
if (!exists("ServerData")) {
  stop("Error: The 'ServerData' object was not loaded from ServerData.rds.")
}

## 3. Runs and replicates

# Generate empty replicates list
# This is to check whether the seeding has been performed correctly throughout the pipeline
replicate_results <- list() 

# Generate empty graphs list > to store all output graphs
total_graphs_list <- list()

# Create graph replicates for each seed
for (r in 1:5) {
  
  # Generate the five networks
  results_from_func <- construct_four_networks(
    empirical_edgelist = ServerData$dog_contact_zone_1,
    kappa_grid = seq(28, 30, by = 1),  # adapt kappa grid as needed
    tau_grid = seq(0.5, 0.95, by = 0.05), # adapt tau grid as needed
    lambda_grid = seq(28, 33, by = 1), # adapt lambda grid as needed
    base_seed = base_seed)
  
  # Extract all objects from output list
  graphs_to_process <- results_from_func$graph_list
  best_p <- results_from_func$best_params
  coms <- results_from_func$communities
  
  # Store the graph in the nested list
  total_graphs_list[[paste0("replicate_", r)]] <- graphs_to_process
  
  # Calculate metrics for each graph in the list
  for (graph_name in names(graphs_to_process)) {
    
    # select each graph
    g <- graphs_to_process[[graph_name]]
    
    metrics <- data.frame(
      Seed = base_seed,
      Replicate = r,
      Graph = graph_name,
      Opt_Lambda = best_p$lambda,
      Opt_Tau    = best_p$tau,
      Opt_Kappa  = best_p$kappa,
      Nodes = vcount(g),
      Edges = ecount(g),
      Density = edge_density(g),
      Avg_Clustering_Coefficient = transitivity(g, type = "global"),
      Avg_Path_Length = mean_distance(g, directed = FALSE),
      Avg_Degree = mean(degree(g)),
      Med_Degree = median(degree(g)),
      Max_Degree = max(degree(g)),
      Avg_Betweenness = mean(betweenness(g)),
      com_sbm = coms$sbm_com,
      com_dcsbm = coms$dcsbm_com,
      stringsAsFactors = FALSE
    )
    replicate_results[[length(replicate_results) + 1]] <- metrics
  }
}


## 4. Bind summary statistics for evaluating the graph generators

# Bind results in a df from replicate_results
if (length(replicate_results) > 0) {
  output_df <- do.call(rbind, replicate_results)
  
  # Create the full path string
  file_name <- paste0("results_seed_", base_seed, "_seed_", max_seed, ".csv")
  final_out_path <- file.path(out_path, file_name)
  
  # Save dataframe > see path set at top of script
  write.csv(output_df, final_out_path, row.names = FALSE)
  cat("\nSuccess: Results saved to", final_out_path, "\n")
} else {
  
  # Report if no results generated
  cat("\nError: No results were generated. Check construct_four_networks logic.\n")
}


## 5. Check validity of results

# Stability Check: Variance across replicates should be 0
stability_check <- output_df %>%
  dplyr::group_by(Seed, Graph) %>%
  dplyr::summarise(SD_Degree = sd(Avg_Degree),
            SD_Bet = sd(Avg_Betweenness),
            SD_Edges = sd(Edges), 
            SD_optL = sd(Opt_Lambda),
            SD_optK = sd(Opt_Kappa),
            SD_optT = sd(Opt_Tau))

# Report stability check
if(sum(stability_check$SD_Degree, na.rm = TRUE) == 0 &&
   sum(stability_check$SD_Bet, na.rm = TRUE) == 0 &&
   sum(stability_check$SD_Edges, na.rm = TRUE) == 0 &&
   sum(stability_check$SD_optL, na.rm = TRUE) == 0 &&
   sum(stability_check$SD_optK, na.rm = TRUE) == 0 &&
   sum(stability_check$SD_optT, na.rm = TRUE) == 0) {
  cat("SUCCESS: Each seed produces identical networks every time.\n")
} else {
  cat("WARNING: Stochastic leakage detected. Check seed offsets.\n")
}

# For conditional use:
stability_sum <- sum(stability_check$SD_Degree, stability_check$SD_Bet, 
                     stability_check$SD_Edges, stability_check$SD_optL, 
                     stability_check$SD_optK, stability_check$SD_optT, 
                     na.rm = TRUE)

# Save stability check for needing replicates or not 
is_stable <- (stability_sum == 0)
if(is_stable){
  
  cat("SUCCESS: Each seed produces identical networks. Proceeding to extraction.\n")
  
  # Extract empirical graph
  empirical_anchor <- total_graphs_list[[1]]$empirical_graph
  
  # Because replicates are identical -> only take replicate 1
  task_candidates <- list(
    spatial = total_graphs_list[[1]]$spatial_graph,
    sbm     = total_graphs_list[[1]]$sbm_graph,
    dcsbm   = total_graphs_list[[1]]$dcsbm_graph,
    random  = total_graphs_list[[1]]$random_graph)
  
  # Output of five graphs (per seed)
  vetted_output <- list(
    anchor = empirical_anchor,
    candidates = task_candidates)
  
  # Save output
  output_rds_path <- file.path(out_path_graphs, paste0("vetted_graphs_task_", task_id, "_seed_", max_seed, ".rds"))
  saveRDS(vetted_output, output_rds_path)
  cat(paste("Vetted graphs for task", task_id, "saved successfully.\n"))
  
} else {
  # WARNING: If stability fails, we do NOT save the graphs. 
  # IF fail, inspect the stability results to see where leakage occurs. 
  cat("CRITICAL WARNING: Stochastic leakage detected (SD > 0).\n")
  cat("The replicates are NOT identical. Extraction aborted.\n")
}
  


### End of script

