##########################################################################################################################

# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
# Created: December 2025; Last edited: September 2026
# In this script, we create an empirical and 5000 synthetic networks for each network generator algorithm for all 
# additional locations.


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
# 9. The data is stored within country and location specific folders. It would be easiest to follow this set-up for 
#    script execution. They are within a BaseData folder stored as "Country_Location", e.g. "Indonesia_Habi".


##########################################################################################################################

### SET UP R ENVIRONMENT ###

# Empty list
rm(list = ls())

# I ran these scripts on a server and had issues with package versions.
# This ensures that all necessary packages used in the script were the ones I installed locally. 
.libPaths(c("/scicore/home/chitnis/derkx0000/R/ubuntu/4.2.1-foss-2022a", .libPaths()))
print(.libPaths())
cat("Does path exist?", dir.exists("/scicore/home/chitnis/derkx0000/R/ubuntu/4.2.1-foss-2022a/igraph"), "\n")

# Load required libraries
library(igraph, lib.loc)
library(lhs)
library(dplyr, lib.loc)
library(randnet, lib.loc)
library(tidyr, lib.loc)
library(ggplot2, lib.loc)
library(purrr)
library(readr)

# Set local root directory 
LOCAL_ROOT_DIR <- "/scicore/home/chitnis/derkx0000/GraphComparison"
setwd(LOCAL_ROOT_DIR)

# Generate a new output folder to store all networks
ifelse(!dir.exists(file.path(LOCAL_ROOT_DIR, "Net_Sens")),
       dir.create(file.path(LOCAL_ROOT_DIR, "Net_Sens")), FALSE)

# Set path for storing output
out_path <- file.path(LOCAL_ROOT_DIR, "Net_Sens")

# Set path for storing summary output of these other locations graphs
out_path_results <- file.path(out_path, "OtherLocations")

# Create output folder for other locations output
ifelse(!dir.exists(out_path_results),
       dir.create(out_path_results), FALSE)

# Set path for storing graphs
out_path_graphs <- file.path(out_path, "Graphs/OtherLocations")

# Create output folder for other locations graphs
ifelse(!dir.exists(out_path_graphs),
       dir.create(out_path_graphs), FALSE)


##########################################################################################################################

### Loading helper functions

#-------------------------------------------------------------------------------------------------------------------------
# Helper function 1: open and reads files in data folder
read_and_validate_file <- function(file_path, col_names) {

  # This function reads files in file path
  #' @param file_path: path of files storage
  #' @param col_names: the names of columns in csv file
  
  # Get basename of files
  file_name <- basename(file_path)
    
  # Read the CSV, skipping the header (row 1). 
  # Using show_col_types = FALSE suppresses messages
  data <- read_csv(file_path, col_names = col_names, skip = 1, show_col_types = FALSE)
    
  # Check if the resulting tibble has any data rows
  if (nrow(data) > 0) {
    
    # If data exists, add the source file name and return it
    return(mutate(data, source_file = file_name))
  } 
  else {
    # If only the header was present (0 rows of data), return NULL and print a message
    message(paste("Warning: File", file_name, "was skipped as it contained only a header (0 data rows)."))
    return(NULL)
  }
}
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------  
# Helper function 2: Open data frames
open_data_create_edgelist <- function(data_folder, folder_dir, col_names, file_pattern){

  # This function takes all csv files with data and uses it to create an edge list
  #' @param data_folder: name of folder where the .csv files are stored
  #' @param folder_dir: the directory of the folder
  #' @param col_names: the names of columns in csv file
  #' @param file_pattern: the pattern to recognize the file for selecting only desired files
  #' @return a dataframe edge list

  ## Step 1: identify files

  # Set file path
  filepath <- file.path(folder_dir, data_folder)
  
  # List all files matching the pattern "D<number>_clean.csv"
  all_files <- list.files(
    path = filepath,
    pattern = file_pattern,
    full.names = TRUE)
  
  # Read all files, filter out empty ones, and combine (using helper function 1)
  combined_data_list <- map(all_files, read_and_validate_file, col_names = col_names)
  
  # Remove all NULL elements (the skipped files)
  data_list_cleaned <- purrr::compact(combined_data_list)
  if (length(data_list_cleaned) == 0) {
    stop("No valid data rows found across all files after filtering. Check file contents.")
  }
  
  # bind_rows combines the list of data frames into a single one.
  combined_data_raw <- bind_rows(data_list_cleaned)
  
  # Turn into edgelist
  edge_list_df <- combined_data_raw %>%
    select(name, peerId2) %>%
    filter(name != peerId2) %>% # remove self loops
    unique() %>%
    rename(dog = name, peer = peerId2)
  
  # Return df with edgelist 
  return(edge_list_df)
  
}
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# Helper function 3: construct five networks (NOTE: incl. NCRG)
construct_five_networks <- function(empirical_edgelist, 
                                    square_edge_size = 1,
                                    kappa_grid, 
                                    tau_grid, 
                                    lambda_grid,
                                    K_max = 20,
                                    base_seed){

    # This function relies on all previous helper functions to fulfill network construction
  #' @param empirical_edgelist: an edgelist as detailed in helper function 1
  #' @param square_edge_size: edge size for grid. Default = 1
  #' @param kappa_grid: parameter range for kappa
  #' @param tau_grid: parameter range for tau
  #' @param lambda_grid: parameter range for lambda
  #' @param K_max: maximum K for BIC
  #' @param base_seed: seed for reproducibility
  #' @return list of lists with graphs, optimal SENCA parameters and SBM/DCSBM communities
  
  ## Step 1: Create empirical network from edgelist
  
  # Create empirical network with helper function 1
  empirical_network <- empirical_net(empirical_edgelist)

  # Save the graph
  empirical_graph <- empirical_network$graph

  # Calculate number of nodes and edges
  empirical_nodes <- length(degree(empirical_graph))
  empirical_edges <- ecount(empirical_graph)
  
 ## Step 2: network generation using SENCA

    cat("\nStarting SENCA graph construction.\n")
  
  # Set square edge sizes and x_coordinates for network construction
  N_nodes <- square_edge_size^2 * empirical_nodes
  
  set.seed(base_seed + 10) # Specific block for X
  x_coordinates <- square_edge_size * randomLHS(N_nodes, 1)
  
  set.seed(base_seed + 20) # Specific block for Y
  y_coordinates <- square_edge_size * randomLHS(N_nodes, 1)
  
  # Initiate grid search
  network_grid <- grid_optimization(empirical_net = empirical_graph,
                                    kappa_grid = kappa_grid,
                                    tau_grid = tau_grid,
                                    lambda_grid = lambda_grid,
                                    base_seed = base_seed + 30)

  # Save optimized parameter combination based on KS distance
  optimal_net <- network_grid[which.min(network_grid$ks_dist), ]
  
  # Get optimal grid parameters
  optimal_parameters <- c(optimal_net$lambda, 
                          optimal_net$tau,
                          optimal_net$kappa)
  
  # Construct final network with optimized parameters
  spatial_graph <- construct_network(parameters = optimal_parameters,
                                     N_nodes_large = N_nodes, 
                                     x_coordinates,
                                     y_coordinates,
                                     seed = base_seed + 40)
  
  if(exists("spatial_graph")){
    cat("\nSpatial graph constructed. Moving on to next algorithm.\n")
  } else {
    stop("\nError: spatial graph not constructed. Aborting process...\n")
  }
  

  ## Step 2: Generate network with SBM and DCSBM
  
  cat("\nStarting SBM and DCSBM graph construction.\n")
  
  # Turn empirical network into adjacency matrix
  A <- as.matrix(as_adjacency_matrix(empirical_graph, sparse = FALSE))
  
  # Re-set seed because of RNG
  set.seed(base_seed + 100) 

  # Bayesian Information Criterion (BIC) likelihood using 'K_max' as set in function arguments
  # and adjacency matrix 'A'. Test both models (SBM and DCSBM).
  BIC_emp <- randnet::LRBIC(A, Kmax = K_max, model = "both") 

  # Save number of communities for SBM and DCSBM
  k_hat_sbm = BIC_emp$SBM.K 
  k_hat_dcsbm <- BIC_emp$DCSBM.K 
  
  # We re-perform this if the number of optimal communities equals the number of maximum communities
  # Increasing K should not affect the graph's optimal K if it does not currently exceed K_max
  # E.g., if K_max = 10, SBM.K is 10 and DCSBM.K is 6, increasing K_max may change SBM.K to 11 (or not)
  # but should not change DCSBM.K from 6 to 7 or higher. 
  while(k_hat_sbm == K_max | k_hat_dcsbm == K_max){
    
    # Increase by increments of 1
    cat("Maximum K has been reached. Increasing K by 1.\n")
    K_max <- K_max + 1

    # Recalculate BIC and community number
    BIC_emp <- randnet::LRBIC(A, Kmax = K_max, model = "both") 
    k_hat_sbm <- BIC_emp$SBM.K 
    k_hat_dcsbm <- BIC_emp$DCSBM.K
  }
  
  # Print communities 
  cat(paste0("\nSBM has ", k_hat_sbm, " communities and DCSBM has ", k_hat_dcsbm, " communities.\n"))
  
  # Save for later inspection
  communities <- list(sbm_com = k_hat_sbm,
                      dcsbm_com = k_hat_dcsbm)
  
  # We use the number of communities to assign clusters with RSSP
  # We do NOT use the Laplacian matrix here!
  set.seed(base_seed + 110) 
  if(k_hat_sbm > 1){
    communities_sbm <- reg.SSP(A,K=k_hat_sbm,lap=FALSE) 
    sbm_generated <- SBM.estimate(A, communities_sbm$cluster)
  
  # If no communities were detected, all have membership to the same community
  } else {
    communities_sbm <- rep(1, empirical_nodes)
    sbm_generated <- SBM.estimate(A, communities_sbm)
  }
  
  # Repeat for DCSBM: now we DO use the Laplacian matrix
  set.seed(base_seed + 120)
  if(k_hat_dcsbm > 1){
    communities_dcsbm <- reg.SSP(A,K=k_hat_dcsbm,lap=TRUE) 
    dcsbm_generated <- DCSBM.estimate(A, communities_dcsbm$cluster)
    
  } else {
    communities_dcsbm <- rep(1, empirical_nodes)
    dcsbm_generated <- DCSBM.estimate(A, communities_dcsbm)
  }
  
  # Use g and B from sbm_generated to create graph
  sbm_graph <- generate_sbm_graph(num_nodes = empirical_nodes, 
                                  c = sbm_generated$g, 
                                  B = sbm_generated$B,
                                  seed = base_seed + 130)
  
  # Use Phat from dcsbm_generated to create graph
  dcsbm_graph <- generate_dcsbm_graph(num_nodes = empirical_nodes, 
                                      Phat = dcsbm_generated$Phat,
                                      seed = base_seed + 140)
  
  if(exists("sbm_graph") && exists("dcsbm_graph")){
    cat("\nBoth SBM graphs constructed. Moving on to next algorithm.\n")
  } else {
    stop("\nError: either SBM or DCSBM not constructed. Aborting process...\n")
  }
  
  ## Step 3: Generate Erdos-Renyi (random) graph (ERM)
  
  cat("\nStarting random graph construction.\n")
  
  # Create random Erdos & Renyi graph: edges = empirical edges
  set.seed(base_seed + 400)
  random_graph = igraph::sample_gnm(empirical_nodes, 
                                          empirical_edges,
                                          FALSE, 
                                          FALSE)
  
  if(exists("random_graph")){
    cat("\nRandom graph constructed. Merging graphs in list now.\n")
  } else {
    stop("\nError: random graph not constructed. Aborting process...\n")
  }

  ## Step 3: Generate Newman-Clustered Random Graph (NCRG)
  
  cat("\nStarting Newman clustered random graph construction.\n")
  
  # First, let's examine the graph characteristics
  diagnose_newman_fit(empirical_graph)
  
  # Make joint degree df (this returns 3 columns: node, s_i, t_i)
  joint_degrees <- calculate_joint_degree_sequence(empirical_graph)
  
  # FIX: Strip out the 'node' column so it only passes exactly 2 columns (s_i, t_i)
  ncr_graph <- random_clustered_graph(joint_degree_sequence = joint_degrees[, c("s_i", "t_i")],
                                      simplify_graph = TRUE,
                                      seed = base_seed + 250)
  
  if(exists("ncr_graph")){
    cat("\nNewman clustered random graph constructed. Merging graphs in list now.\n")
  } else {
    stop("\nError: Newman clustered random graph not constructed. Aborting process...\n")
  }
  
  ## Step 4: Generate Erdos-Renyi (random) graph
  
  cat("\nStarting random graph construction.\n")
  
  # Create random Erdos & Renyi graph: edges = empirical edges
  set.seed(base_seed + 400)
  random_graph = igraph::sample_gnm(empirical_nodes, 
                                    empirical_edges,
                                    FALSE, 
                                    FALSE)
  
  if(exists("random_graph")){
    cat("\nRandom graph constructed. Merging graphs in list now.\n")
  } else {
    stop("\nError: random graph not constructed. Aborting process...\n")
  }
  

  ## Step 5: list and return results 


  # Create list of graphs
  graph_list <- list(empirical_graph = empirical_graph, 
                     spatial_graph = spatial_graph, 
                     sbm_graph = sbm_graph, 
                     dcsbm_graph = dcsbm_graph, 
                     ncr_graph = ncr_graph,
                     random_graph = random_graph)
  cat(paste0("\nUsed the following parameters for grid optimization. Kappa = ", optimal_net$kappa, 
             "; Tau = ", optimal_net$tau, "; Lambda = ", optimal_net$lambda, ".\n"))
  
  
  # Return graph list
  list(graph_list = graph_list,
       best_params = optimal_net,
       communities = communities)
  
}
#-------------------------------------------------------------------------------------------------------------------------

# Source 00 and 01a scripts for remaining helper functions
source("00_constructNetworkFunctions.R")
source("01a_Run_NCRG.R")

##########################################################################################################################

### Execute code

## 1. Set up arguments and seeding from .sh script

# Set arguments
args <- commandArgs(trailingOnly = TRUE)

# The SLURM script passes: SEED, COUNTRY, LOCATION, FILE_DIR
if (length(args) < 4) {
  stop("Usage: Rscript 9_Graphs_Other_Locations.R <SEED> <COUNTRY> <LOCATION> <FILE_DIR>", call. = FALSE)
}

# Assign arguments
TASK_ID  <- as.numeric(args[1]) # This is the SEED/Array ID
COUNTRY  <- as.character(args[2])
LOCATION <- as.character(args[3])
FILE_DIR <- as.character(args[4])

# Set seed using task id
base_seed <- 1000 + (TASK_ID * 10000) # Give each task 10,000 "room"
max_seed = 5000 # adapt


## 2. Load data


# Set general col names and file pattern
col_names_contact_file <- c("name","deviceNr","deviceId","timestamp","peerId","rssi","TIME","peerId2")

# Set file pattern
file_pattern <- "^D[0-9]+_clean\\.csv$"

# Set data folder with country and location
data_folder <- paste0(COUNTRY, "_", LOCATION)

# Create edgelist using helper function 1
edgelist <- open_data_create_edgelist(data_folder, FILE_DIR, col_names_contact_file, file_pattern)
  

## 3. Runs and replicates


# Generate empty lists
replicate_results <- list()
total_graphs_list <- list()

# Create graph replicates for each seed
for (r in 1:5) {
  
  # Script deviation from 01 Script! 
  # Given that we are testing multiple locations in one script, with graphs of different sizes and densities,
  # we're using a floor and ceiling of lambda (for grid optimization in SENCA) and iterate over it. 
  # This takes a little longer 
  
  # Set max lambda
  current_lambda_max <- 18

  # Don't let it drop below the start of your grid
  floor_lambda <- 3  
  success <- FALSE
  
  # As long as we're not successful, we keep going 
  while (!success && current_lambda_max >= floor_lambda) {
    
    # Create the grid for this attempt
    temp_lambda_grid <- seq(3, current_lambda_max, by = 1)
    
    results_from_func <- tryCatch({
      
      message(paste("Attempting construction with Lambda max:", current_lambda_max))
      
      # Try the construction
      res <- construct_five_networks(
        empirical_edgelist = edgelist,
        kappa_grid  = seq(3, 15, by = 1),
        tau_grid    = seq(0.5, 0.95, by = 0.05),
        lambda_grid = temp_lambda_grid,
        base_seed   = base_seed)
      
      success <- TRUE  # If we reach here, it worked!
      res              # Return the result to the variable
      
    }, error = function(e) {
      message(paste("Failed at Lambda", current_lambda_max, ":", e$message))
      current_lambda_max <<- current_lambda_max - 1  # Decrement
      return(NULL)
    })
  }
  
  if (!success) {
    stop("Network construction failed even after scaling Lambda down to minimum.")
  }

  
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


# Bind results (df contains all replicates)
if (length(replicate_results) > 0) {
  output_df <- do.call(rbind, replicate_results)
  
  # Create the full path string
  file_name <- paste0(LOCATION, "_results_seed_", base_seed, "_seed_", max_seed, ".csv")
  final_out_path <- file.path(out_path_results, file_name)
  
  # Save dataframe
  write.csv(output_df, final_out_path, row.names = FALSE)
  cat("\nSuccess: Results saved to", final_out_path, "\n")
} else {
  
  # Report if not results generated
  cat("\nError: No results were generated. Check construct_five_networks logic.\n")
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
    spatial   = total_graphs_list[[1]]$spatial_graph,
    sbm       = total_graphs_list[[1]]$sbm_graph,
    dcsbm     = total_graphs_list[[1]]$dcsbm_graph,
    ncrg      = total_graphs_list[[1]]$ncr_graph, 
    random    = total_graphs_list[[1]]$random_graph)
  
  # Output of five graphs (per seed)
  vetted_output <- list(
    anchor = empirical_anchor,
    candidates = task_candidates)
  
  # Save output
  output_rds_path <- file.path(out_path_graphs, paste0(LOCATION, "_vetted_graphs_task_",TASK_ID, "_seed_", max_seed, ".rds"))
  saveRDS(vetted_output, output_rds_path)
  cat(paste("Vetted graphs for task", TASK_ID, "saved successfully.\n"))
  
} else {
  # WARNING: If stability fails, we do NOT save the graphs.
  # This prevents "leaky" data from reaching your SEIR simulations.
  cat("CRITICAL WARNING: Stochastic leakage detected (SD > 0).\n")
  cat("The replicates are NOT identical. Extraction aborted.\n")
}



### End of script

