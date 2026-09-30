###########################################################################################################################

# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
# Created: December 2025; Last edited: September 2026
# In this script, we select which graphs will be used for downstream modeling. 

# IMPORTANT NOTES:
# 0. The is a standalone script. There is no corresponding .sh script. 
# 1. This script requires the output of the 01 script. 

# For any issues, feel free to report a new issue on this GitHub repo or contact inez.derkx@swisstph.ch

##########################################################################################################################

### SET UP R ENVIRONMENT ###

# Empty list
rm(list = ls())

# Load required libraries
library(igraph)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(readr)
library(data.table)

# Set local root directory 
LOCAL_ROOT_DIR <- "/scicore/home/chitnis/derkx0000/GraphComparison"
setwd(LOCAL_ROOT_DIR)

# Output folders and paths. 
out_path <- file.path(LOCAL_ROOT_DIR, "Net_Sens")
ifelse(!dir.exists(file.path(out_path, "Plots")),
       dir.create(file.path(out_path, "Plots")), FALSE)

##########################################################################################################################

### Helper functions

#-------------------------------------------------------------------------------------------------------------------------
# 1. Function to summarize graphs
summarize_graphs <- function(path,
                             csv_pattern,
                             seed_nr,
                             max_task_id = NULL){

  # This function creates summaries of all generated graphs from script 01. 
  #' @param path: path where results df is stored ("results_seed_", base_seed, "_seed_", max_seed, ".csv")
  #' @param csv_pattern: pattern for file name recognition, e.g. "results_seed_(\\d+)_seed_.*\\.csv" from above.
  #' @param seed_nr: should correspond to the number in path ("_seed_") for file recognition, aka total array size. 
  #' @param max_task_id: number of array tasks that will be summarized
  #' @return A summary dataframe. The function saves and plots this file automatically

  ## Step 1: Source and gather files

  # Adapt path to seed number to locate csv files 
  seed <- paste0("Seed_", seed_nr)
  csv_path <- file.path(path, seed)
  
  # Add plot path
  plots_path <- file.path(path, "Plots")
  
  # output path 
  out_seed <- paste0("Seed_", max_task_id)
  output_path <- file.path(path, out_seed)
  ifelse(!dir.exists(output_path),
         dir.create(output_path), FALSE)
  
  # list all ouptut files 
  file_list <- list.files(path = csv_path, 
                          pattern = csv_pattern, 
                          full.names = TRUE)
  
  # Optionally filter by task ID
  if (!is.null(max_task_id)) {
    pattern <- paste0("^results_seed_(\\d+)_seed_", seed_nr, "\\.csv$")
    base_seeds <- as.integer(sub(pattern, "\\1", basename(file_list)))
    max_base_seed <- (max_task_id * 10000) + 1000
    file_list <- file_list[!is.na(base_seeds) & base_seeds <= max_base_seed]
  }
  
  # Combine all csv files in one dataframe
  sensitivity_data <- do.call(rbind, lapply(file_list, read.csv))

  ## Step 2: Transform and save summary data
  
  # Save in new summary file
  summary_check <- sensitivity_data %>%
    dplyr::select(Seed, Replicate, Graph, Avg_Degree, Med_Degree, Max_Degree, Avg_Betweenness, Density, 
                  Avg_clust = Avg_Clustering_Coefficient, Opt_Lambda, Opt_Tau, Opt_Kappa) %>%
    tidyr::pivot_longer(cols = Avg_Degree:Opt_Kappa, 
                        names_to = c("Params"),
                        values_to = c("Values")) %>%

    # Save only one replicate > in previous script we checked that all should be identical.
    dplyr::filter(Replicate == 1) %>%
    dplyr::mutate(SeedTotal = max_task_id)
  
  # Set summary file name — include task ID subset in filename if filtered
  summary_filename <- if (!is.null(max_task_id)) {
    paste0("summary_check_top", max_task_id, ".csv")
  } else {
    "summary_check.csv"
  }
  
  # Set path to save
  summary_path <- file.path(csv_path, summary_filename)
  write.csv(summary_check, summary_path, row.names = FALSE)

  ## Step 3: Plot summary data for diagnostics
  
  # Help plot: this plot is for examination of data. It is not a final plot!
  params_summary <- ggplot(summary_check, aes(x = Graph, 
                                              y = Values, 
                                              group = Graph,
                                              fill = Graph)) +
    geom_boxplot() +
    theme_minimal() +
    scale_fill_manual(values = palette_paper) +
    facet_wrap(~Params, scales = 'free', ncol = 3, axis.labels = "margins") +
    labs(title = "Network Sensitivity over Random Seeds",
         subtitle = "Variance shows the impact of randomness on the generator",
         y = "Network metric value",
         x = "Graph type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
  
  # Plot file nam — include task ID subset in filename if filtere
  plot_name <- if (!is.null(max_task_id)) {
    paste0("params_summary_seed_", seed_nr, "_top", max_task_id, ".png")
  } else {
    paste0("params_summary_seed_", seed_nr, ".png")
  }

  # Save path and plot
  plot1_path <- file.path(plots_path, plot_name)
  ggsave(plot1_path, params_summary, width = 10, height = 10)
  cat(paste0("Parameter summary has been saved to: ", plot1_path, ".\n"))
  
  # Return summary file
  return(summary_check)
}
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# 2. Function to select graph matches 
select_matches <- function(path, 
                           graph_types, 
                           graph_pattern, 
                           n_to_keep = 100,
                           mypalette,
                           seed_nr,
                           max_task_id){
  
  
  # The main goal of this function is to find the 'n_to_keep' number of files (default = 100) that would best represent
  # the average degree distribution of a specific graph type based on all 5000 graphs. For more details on this logic, 
  # read the methods section of the corresponding manuscript. 
  #' @param path: base path containing folder with graphs 
  #' @param graph_types string of graph types (e.g. ("spatial", "sbm", "dcsbm", "random", "newclust_graph"))
  #' @paramgraph_pattern: pattern string (e.g. "vetted_graphs_task_.*\\.rds") for file recognition
  #' @param n_to_keep: how many graphs to keep for downstream analysis? Default = 100
  #' @param mypalette: what colour palette are you using for your diagnostics graphs?
  #' @param seed_nr: should correspond to the number in path ("_seed_") for file recognition, aka total array size. 
  #' @param max_task_id: number of array tasks that will be summarized
  #' @return final ensemble of graphs to use for SEIR/SIS simulations

  ## Step 1: Set paths

  # First adapt out path to seed
  seed <- paste0("Graphs/Seed_", seed_nr)
  graphs_path <- file.path(path, seed)
  
  # Add plot path
  plots_path <- file.path(path, "Plots")
  
  # output path 
  out_seed <- paste0("Graphs/Seed_", max_task_id)
  output_path <- file.path(path, out_seed)
  ifelse(!dir.exists(output_path),
         dir.create(output_path), FALSE)
  
  
  ## Step 2: create helper functions
  
  # Helper function 1: get degree summary stats for comparison
  get_degree_stats <- function(g) {

    # This function takes a graph and calculates six degree statistics:
    # The minimum, 0.25 quantile, mean, median, 0.75 quantile, and max,
    # as to capture the full shape of the graph's degree distribution. 
    #' @param 'g': an igraph graph object
    #' @return the summary stats

    # Get the degrees of graph g
    d <- as.numeric(igraph::degree(g))

    # Calculate stats
    stats <- c(
      min  = min(d),
      q25  = quantile(d, 0.25, names = FALSE),
      mean = mean(d),
      median = median(d),
      q75  = quantile(d, 0.75, names = FALSE),
      max  = max(d))
    
    # Return stats
    return(stats)
  }
  
  # Helper function 2: Kolmogorov-Smirnov distance
  get_ks_dist <- function(g_candidate, g_empirical) {

    # This function computes the KS distance between two degree distributions
    # ks.test is part of the stats packagea
    #' @param 'g_candidate': candidate graph (igraph object)
    #' @param 'g_empirical': empirical graph (igraph object)
    #' @return KS distance statistic

    # Ensure both graphs exist
    if (is.null(g_candidate) || is.null(g_empirical)) return(NA)

    # Calculate degrees
    d1 <- as.numeric(igraph::degree(g_candidate))
    d2 <- as.numeric(igraph::degree(g_empirical))

    # ks.test returns a list; 'statistic' is the D value (distance)
    return(as.numeric(ks.test(d1, d2)$statistic))
  }
  
  # Optionally filter by task ID
  if (!is.null(max_task_id)) {
    
    # Obtain graphs from rds files
    file_list <- list.files(path = graphs_path, 
                            pattern = graph_pattern, 
                            full.names = TRUE)
    
    # Set pattern of files
    pattern <- paste0("^vetted_graphs_task_(\\d+)_seed_", seed_nr, "\\.rds$")
    
    # Extract task ids
    task_ids <- as.integer(sub(pattern, "\\1", basename(file_list)))
    
    # Retain only those with desired task ids
    file_list <- file_list[!is.na(task_ids) & task_ids <= max_task_id]
    
  } else{
    
    # Obtain graphs from rds files
    file_list <- list.files(path = graphs_path, 
                            pattern = graph_pattern, 
                            full.names = TRUE)
  }
  
  # Extract empirical graph
  graph1 <- readRDS(file_list[1])
  empirical <- graph1$anchor
  emp_degrees <- table(igraph::degree(empirical))
  emp_degree_df <- data.frame(task_id = 1,
                              degree  = as.numeric(names(emp_degrees)),
                              frequency = as.numeric(emp_degrees),
                              type = "empirical")
  cat("Empirical degree stats:")
  print(get_degree_stats(empirical))
  
  # Betweenness
  emp_between <- table(round(igraph::betweenness(empirical), 0))
  emp_between_df <- data.frame(task_id = 1,
                               between  = as.numeric(names(emp_between)),
                               frequency = as.numeric(emp_between),
                               type = "empirical")
  
  # Map seed/task id to the actual graph objects
  all_task_data <- lapply(file_list, function(f) {
    
    # Extract Task ID from filename for linking
    task_id <- as.numeric(gsub(".*?task_([0-9]+).*", "\\1", basename(f)))
    data <- readRDS(f)
    
    # Return a list where each graph is tagged with its origin task_id
    # This ensures 'linkage' to the specific seed it was made with
    return(list(
      task_id        = task_id,
      anchor         = data$anchor,
      spatial        = data$candidates$spatial,
      sbm            = data$candidates$sbm,
      dcsbm          = data$candidates$dcsbm,
      random         = data$candidates$random,
      newclust_graph = data$candidates$newclust_graph))
  })
  
  # Initialize lists
  final_ensemble <- list()
  global_distribution_list <- list()
  betweenness_distribution_list <- list()
  all_centroids <- list()
  best_graphs_indices_list <- list()
  ks_distances_list <- list()
  
  ### Go through all graphs
  
  # For each graph type ...
  for (type in graph_types) {
    
    # Process them separately
    cat(paste("\n--- Processing:", type, "---\n"))
    
    # Create a list of length equal to the number of tasks
    candidates <- lapply(all_task_data, function(x) {
      list(task_id = x$task_id, 
           graph = x[[type]]) # If 'type' isn't found, this becomes NULL
    })
    
    
    # Check that number of graphs found = number of seeds
    cat(sprintf("Found %d candidate graphs for %s\n", length(candidates), type))
    
    
    ### Obtain degree distributions
    
    # Extract degree distribution for all graphs
    type_distributions <- lapply(candidates, function(c) {
      deg_counts <- table(round(igraph::degree(c$graph), 0))
      data.frame(task_id   = c$task_id,
                 degree    = as.numeric(names(deg_counts)),
                 frequency = as.numeric(deg_counts),
                 type      = type)
    })
    
    # Bind all seeds for this specific type
    global_distribution_list[[type]] <- do.call(rbind, type_distributions)
    
    
    ### Obtain betweenness distributions
    bet_distributions <- lapply(candidates, function(c) {
      bet_counts <- table(round(igraph::betweenness(c$graph),0))
      data.frame(task_id = c$task_id,
                 between  = as.numeric(names(bet_counts)),
                 frequency = as.numeric(bet_counts),
                 type = type)
    })
    betweenness_distribution_list[[type]] <- do.call(rbind, bet_distributions)
    
    
    ### Obtain summary stats and rank graphs
    
    # Calculate errors for all summary points across ALL 100 or 1000 candidates
    summary_matrix <- sapply(candidates, function(c) {
      
      # Calculate stats for the graph object inside the list
      sim_stats <- get_degree_stats(c$graph)
    }) %>% t() 
    print(head(summary_matrix))
    
    # KS values
    ks_values <- sapply(candidates, function(c) get_ks_dist(c$graph, empirical))
    
    # Compute the centroid of all summary stats
    model_centroid <- colMeans(summary_matrix, na.rm = TRUE)
    cat("Computed Model Centroid (Average Stats):\n")
    print(model_centroid)
    
    # Transform for saving
    centroid_row <- as.data.frame(t(model_centroid))
    centroid_row$type <- type
    all_centroids[[type]] <- centroid_row
    
    # Create error matrix
    error_matrix <- abs(sweep(summary_matrix, 2, model_centroid, "-"))
    
    # Rank the errors (1 is best fit)
    rank_matrix <- apply(error_matrix, 2, rank, ties.method = "min")
    cat("Errors of best fitting graphs (summary stats - avg summary stats):\n")
    print(head(rank_matrix, 20))
    
    # Sum the ranks for each graph (Total score out of a possible 5 * 100)
    rank_sums <- rowSums(rank_matrix)
    
    # Select the top 20 best-fitting graphs based on Rank Sum
    best_indices <- order(rank_sums)[1:min(n_to_keep, length(rank_sums))]
    best_graphs_indices_list[[type]] <- best_indices
    
    
    # Store the 'Elite' graphs + their origin metadata
    final_ensemble[[type]] <- lapply(best_indices, function(idx) {
      list(
        origin_task_id = candidates[[idx]]$task_id,
        graph = candidates[[idx]]$graph,
        rank_score = rank_sums[idx])
    })
    
    # Store all KS distances for later analysis if needed
    ks_distances_list[[type]] <- data.frame(
      task_id = sapply(candidates, function(x) x$task_id),
      ks_dist = ks_values,
      graph_type = type
    )
    
    cat(sprintf("Successfully selected the top %d %s graphs.\n", n_to_keep, type))
  }
  
  # Bind the rows together
  centroids_df <- do.call(rbind, all_centroids)
  cat("Full summary of degree stats for each graph type:\n")
  print(centroids_df)
  all_ks_df <- rbindlist(ks_distances_list)
  
  
  ### Save betweenness distribution
  all_graphs_bet_df <- do.call(rbind, betweenness_distribution_list)
  all_graphs_bet_df <- rbind(all_graphs_bet_df, emp_between_df)
  
  ### Plot degree distribution
  
  # Put all degree distributions together
  all_graphs_dist_df <- do.call(rbind, global_distribution_list)
  
  # Add empirical df
  all_graphs_dist_df <- rbind(all_graphs_dist_df, emp_degree_df)
  
  
  # Obtain the mean and quantiles average over all seeds
  degree_summarized <- all_graphs_dist_df %>%
    group_by(type, degree) %>%
    summarise(q_25 = quantile(frequency, 0.25, names = FALSE),
              mean = mean(frequency),
              q_75 = quantile(frequency, 0.75, names = FALSE))
  
  # Plot average degree distributions
  degree_plot <- ggplot(degree_summarized, aes(x = degree, 
                                               y = mean, 
                                               group = type,
                                               color = type,
                                               fill = type)) +
    geom_line() +
    geom_ribbon(aes(ymin=q_25, ymax=q_75, colour = type), alpha = 0.3) +
    scale_fill_manual(values = mypalette) +
    scale_color_manual(values = mypalette) +
    labs(title = "Degree distribution by graphp type",
         subtitle = "Graph shows average and quantiles over all seeds",
         y = "Frequency",
         x = "Degree") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # Save plot
  plot_name <- paste0("degree_plot_seed", seed_nr, ".png")
  plot_name_path <- file.path(plots_path, plot_name)
  ggsave(plot_name_path, degree_plot, width = 10, height = 10)
  
  # Save best graphs
  best_indices_df <- do.call(rbind, lapply(names(best_graphs_indices_list), function(type) {
    data.frame(
      graphs = best_graphs_indices_list[[type]],
      graph_type = type,
      stringsAsFactors = FALSE)
  }))
  
  ### Save all output files
  
  # Best graphs
  write.csv(best_indices_df, file.path(output_path, "best_graphs_indices.csv"), row.names = FALSE)
  
  # KS distances
  write.csv(all_ks_df, file.path(output_path, "ks_distances_all_sims.csv"), row.names=FALSE)
  
  # Final ensemble
  ensemble_path <- file.path(output_path, "MASTER_ENSEMBLE_FOR_DM.rds")
  saveRDS(final_ensemble, ensemble_path)
  
  # Full degree distribution
  degree_path <- file.path(output_path, "all_seeds_degree_distributions.csv")
  write.csv(all_graphs_dist_df, degree_path, row.names = FALSE)
  
  # Full betweenness distribution
  betweenness_path <- file.path(output_path, "all_seeds_between_distributions.csv")
  write.csv(all_graphs_bet_df, betweenness_path, row.names = FALSE)
  
  # Degree summary
  baseline_path <- file.path(output_path, "degree_model_baselines.csv")
  write.csv(centroids_df, baseline_path, row.names = FALSE)
  
  # Return the ensemble file
  return(final_ensemble)
}
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# 3. Function to compute Mahalanobis distance
mh_distance <- function(path, graph_types, graph_pattern){
  
  #' Compare measure of synthetic networks to empirical network
  #' @param path General path to where relevant folders are stored
  #' @param graph_types Character vector of models (e.g., c("ER"))
  #' @param graph_pattern Pattern of files that store the graph
  #' @return A summary table ranking models by their proximity to the empirical data
  
  
  # Helper function: network topology
  get_network_topology <- function(g) {
    if (is.null(g) || !igraph::is_igraph(g)) return(rep(NA, 5))
    
    # Calculate key structural metrics
    # Using 'mean' of distributions to represent the global state
    stats <- c(
      avg_degree     = mean(igraph::degree(g)),
      med_degree     = median(igraph::degree(g)),
      max_degree     = max(igraph::degree(g)),
      avg_density    = mean(igraph::edge_density(g)),
      avg_clust      = igraph::transitivity(g, type = "global"),
      avg_between    = mean(igraph::betweenness(g)))
    return(stats)
  }
  
  # Obtain graphs from rds files
  file_list <- list.files(path = path, 
                          pattern = graph_pattern, 
                          full.names = TRUE)
  
  # Get empirical network
  graph1 <- readRDS(file_list[1])
  empirical_graph <- graph1$anchor
  
  # Obtain summary stats for empirical network
  cat("Extracting empirical network stats...\n")
  emp_vector <- get_network_topology(empirical_graph)
  
  # Map seed/task id to the actual graph objects
  cat("Extracting all synthetic graphs...\n")
  all_task_data <- lapply(file_list, function(f) {
    
    # Extract Task ID from filename for linking
    task_id <- as.numeric(gsub(".*?task_([0-9]+).*", "\\1", basename(f)))
    data <- readRDS(f)
    
    # Return a list where each graph is tagged with its origin task_id
    # This ensures 'linkage' to the specific seed it was made with
    return(list(
      task_id = task_id,
      anchor  = data$anchor,
      spatial = data$candidates$spatial,
      sbm     = data$candidates$sbm,
      dcsbm   = data$candidates$dcsbm,
      random  = data$candidates$random,
      newclust_graph = data$candidates$newclust_graph))
  })
  
  # Create empty list
  performance_list <- list()
  
  # Loop through the four graph types
  cat("Extracting synthetic network stats...\n")
  
  for (type in graph_types) {
    cat(sprintf("Processing population: %s (%d graphs)...\n", type, length(all_task_data)))
    
    # Extract the specific graph type from all tasks
    type_graphs <- lapply(all_task_data, function(x) x[[type]])
    
    # Calculate fingerprints for all 5000 simulations
    sim_matrix <- do.call(rbind, lapply(type_graphs, get_network_topology))
    
    # Clean data (remove rows with Inf or NA)
    sim_matrix <- sim_matrix[complete.cases(sim_matrix) & !rowSums(is.infinite(sim_matrix)), ]
    
    if (nrow(sim_matrix) < 100) {
      warning(paste("Insufficient valid simulations for type:", type))
      next
    }
    
    # Mahalanobis distance 
    mu_model <- colMeans(sim_matrix) # The average version of this model
    sigma_model <- cov(sim_matrix) #The internal variance/correlation of this model's features
    
    # Regularization to prevent singular matrix errors
    diag(sigma_model) <- diag(sigma_model) + 1e-6
    
    # Calculate D^2 distance from empirical point to the model cloud
    d2 <- mahalanobis(x = matrix(emp_vector, nrow = 1), 
                      center = mu_model, 
                      cov = sigma_model)
    
    # Calculate actual Mahalanobis distance
    dist_val <- sqrt(d2)
    print(dist_val)
    
    performance_list[[type]] <- data.frame(
      Graph_Type = type,
      Mahalanobis_Dist = dist_val,
      Valid_Sims = nrow(sim_matrix))
  }
  
  # Combine results for function output
  results <- do.call(rbind, performance_list) %>%
    arrange(Mahalanobis_Dist) %>%
    mutate(
      Rank = row_number(),
      Similarity_Score = 1 / (1 + Mahalanobis_Dist)) # Normalize 0 to 1
  return(results)
  cat("\n--- Analysis Complete. ---\n")
}
#-------------------------------------------------------------------------------------------------------------------------

##########################################################################################################################

### Select graphs


## Step 1: set up for graph examination


# Colour palette
palette_paper <- c("#3A405A", "#FF8465", "#99B2DD", "#FCD2A2", "#A26F39", "#5C8A6F")

# Set patterns
graphs_pattern = "vetted_graphs_task_.*\\.rds"
summary_pattern = "results_seed_(\\d+)_seed_.*\\.csv"

# Set graph types
five_graphs = c("spatial", "sbm", "dcsbm", "random", "newclust_graph")


## Step 2: examine first 100 graphs


# Summary
summary_check_100 <- summarize_graphs(path = out_path,
                                      csv_pattern = summary_pattern,
                                      seed_nr = 5000, #this should be 100
                                      max_task_id = 100)

# Save summary
write.csv(summary_check_100, "Net_Sens/Seed_100/summary_check.csv", row.names = F)

# Select graphs to keep
matches_seed100 <- select_matches(path = out_path, 
                                  graph_types = five_graphs,
                                  graph_pattern = graphs_pattern,
                                  n_to_keep = 100,
                                  mypalette = palette_paper,
                                  seed_nr = 5000,
                                  max_task_id = 100)


## Step 3: examine first 1000 graphs


# Summary
summary_check_1000 <- summarize_graphs(path = out_path,
                                       csv_pattern = summary_pattern,
                                       seed_nr = 5000, #this should be 100
                                       max_task_id = 1000) 
# Save summary
write.csv(summary_check_1000, "Net_Sens/Seed_1000/summary_check.csv", row.names = F)

# Select graphs
matches_seed1000 <- select_matches(path = out_path, 
                                   graph_types = five_graphs,
                                   graph_pattern = graphs_pattern,
                                   n_to_keep = 100,
                                   mypalette = palette_paper,
                                   seed_nr = 5000,
                                   max_task_id = 1000)


## Step 4: examine first 2000 graphs


# Summary
summary_check_2000 <- summarize_graphs(path = out_path,
                                       csv_pattern = summary_pattern,
                                       seed_nr = 5000, 
                                       max_task_id = 2000) 

# Save summary
write.csv(summary_check_2000, "Net_Sens/Seed_2000/summary_check.csv", row.names = F)

# Select graphs
matches_seed2000 <- select_matches(path = out_path, 
                                   graph_types = five_graphs,
                                   graph_pattern = graphs_pattern,
                                   n_to_keep = 100,
                                   mypalette = palette_paper,
                                   seed_nr = 5000,
                                   max_task_id = 2000)


## Step 5: examine first 5000 graphs


# Seed 5000
summary_check_5000 <- summarize_graphs(path = out_path,
                                       csv_pattern = summary_pattern,
                                       seed_nr = 5000, #this should be 100
                                       max_task_id = 5000) 
write.csv(summary_check_5000, "Net_Sens/Seed_5000/summary_check.csv", row.names = F)

matches_seed5000 <- select_matches(path = out_path, 
                                   graph_types = five_graphs,
                                   graph_pattern = graphs_pattern,
                                   n_to_keep = 100,
                                   mypalette = palette_paper,
                                   seed_nr = 5000,
                                   max_task_id = NULL)

graph_path = file.path(out_path, paste0("Graphs/Seed_5000"))
mahala_dist5000 <- mh_distance(path = graph_path, 
                               graph_types = five_graphs,
                               graph_pattern = graphs_pattern)
write.csv(mahala_dist5000, "Net_Sens/Seed_5000/mahala_dist5000", row.names=F)



### End of script ###