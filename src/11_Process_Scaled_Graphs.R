###########################################################################################################################

# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
# Created: 22 December 2025; Last edited: 13 January 2026
# Script compares different network construction algorithms

# Legend of script:

### IS FOR NEW SECTIONS IN CAPITAL ###
### Is for headings (e.g., a new function)
# is for 'small' commands (e.g., rename or merge)

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

# Output folders
out_path <- file.path(LOCAL_ROOT_DIR, "Scaling")

##########################################################################################################################


### Set functions

# 1. Function to summarize graphs
summarize_graphs <- function(path,
                             csv_pattern,
                             seed_nr){
  
  # Set paths
  plots_path <- file.path(path, "Plots")
  
  # Aggregate results from all graphs
  file_list <- list.files(path = path, 
                          pattern = csv_pattern, 
                          full.names = TRUE)
  
  # Combine all csv files in one df
  sensitivity_data <- do.call(rbind, lapply(file_list, read.csv))
  summary(sensitivity_data)
  
  # Save in new summary file
  summary_check <- sensitivity_data %>%
    dplyr::select(Seed, Replicate, Graph, Edges, Avg_Degree, Med_Degree, Max_Degree, Avg_Betweenness, Density, 
                  Avg_clust = Avg_Clustering_Coefficient, Opt_Lambda, Opt_Tau, Opt_Kappa) %>%
    tidyr::pivot_longer(cols = Edges: Opt_Kappa, 
                        names_to = c("Params"),
                        values_to = c("Values")) %>%
    dplyr::filter(Replicate == 1) %>%
    dplyr::mutate(SeedTotal = seed_nr)
  
  # Save summary file
  summary_path <- file.path(path, "scaled_summary_check.csv")
  write.csv(summary_check, summary_path, row.names = FALSE)
  
  # Plot summary file
  params_summary <- ggplot(summary_check, aes(x = Graph, 
                                              y = Values, 
                                              group = Graph,
                                              fill = Graph)) +
    geom_boxplot() +
    theme_minimal() +
    scale_fill_manual(values = palette_paper) +
    facet_wrap(~Params, scales = 'free', ncol = 4, axis.labels = "margins") +
    labs(title = "Network Sensitivity over Random Seeds",
         subtitle = "Variance shows the impact of randomness on the generator",
         y = "Network metric value",
         x = "Graph type") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom")
  
  # Save plot
  plot_name <- paste0("params_summary_seed_", seed_nr, ".png")
  plot1_path <- file.path(plots_path, plot_name)
  ggsave(plot1_path, params_summary, width = 10, height = 10)
  
  # Print that plot is made
  cat(paste0("Parameter summary has been saved to: ", plot1_path, ".\n"))
  
  # Return summary file
  return(summary_check)
  
}

# 2. Function to select graph matches 
select_matches <- function(path, 
                           graph_types, 
                           graph_pattern, 
                           n_to_keep,
                           mypalette,
                           seed_nr){
  
  
  # Function to examine and select graphs for disease modelling:
  #' @param path
  #' @param graph_types 
  #' @param graph_pattern
  #' @param n_to_keep
  #' @param mypalette
  #' @param seed_nr
  
  # Set paths
  graphs_path <- file.path(path, "Graphs")
  if (!dir.exists(graphs_path)) dir.create(graphs_path)
  plots_path  <- file.path(path, "Plots")
  if (!dir.exists(plots_path)) dir.create(plots_path)
  
  ### Function set up: helper functions and files
  
  # 1. Helper function: get degree summary stats for comparison
  get_degree_stats <- function(d) {
    
    if (length(d) == 0) return(rep(NA, 6))
    stats <- c(
      min  = min(d),
      q25  = quantile(d, 0.25, names = FALSE),
      mean = mean(d),
      median = median(d),
      q75  = quantile(d, 0.75, names = FALSE),
      max  = max(d))
    return(stats)
  }
  
  # 2. Helper function: compute KS distance between two degree distributions
  get_ks_dist <- function(d_candidate, d_empirical) {
    
    if (length(d_candidate) == 0 || length(d_empirical) == 0) return(NA)  
    
    # ks.test returns a list; 'statistic' is the D value (distance)
    return(as.numeric(ks.test(d_candidate, d_empirical)$statistic))
  }
  
  # Obtain graphs from rds files
  file_list <- list.files(path = graphs_path, 
                          pattern = graph_pattern, 
                          full.names = TRUE)
  if (length(file_list) == 0) stop("No files found matching the pattern.")
  
  ## Empirical baseline 
  
  # Extract empirical graph and degrees
  emp_data <- readRDS(file_list[1])
  empirical <- emp_data$anchor
  emp_degree_vector <- as.numeric(igraph::degree(empirical))
  emp_degrees <- table(emp_degree_vector)
  emp_degree_df <- data.frame(task_id = 1,
                              degree  = as.numeric(names(emp_degrees)),
                              frequency = as.numeric(emp_degrees),
                              type = "empirical")
  cat("Empirical degree stats:")
  print(get_degree_stats(emp_degree_vector))
  rm(emp_data)
  
  # Calculate degree fractions
  total_nodes <- vcount(empirical)
  degree_frac_vector <- emp_degree_vector / total_nodes
  degree_fractions <- table(degree_frac_vector)
  emp_degree_frac_df <- data.frame(task_id = 1,
                              degree  = as.numeric(names(degree_fractions)),
                              frequency = as.numeric(degree_fractions),
                              type = "empirical")
  
  
  ## Start with data extraction
  
  cat("Start with degree stats extraction from all graphs.\n")
  
  # Initialize lists
  meta_stats_list <- list() # Stores task_id, type, and the 6 summary stats
  dist_plot_list <- list() # Stores the degree frequency tables (small)
  frac_plot_list <- list() # Store degree fraction frequency tables (small)
  
  # Loop throuh file list to get stats and distances
  for (i in seq_along(file_list)) {
    f <- file_list[i]
    if (i %% 50 == 0) cat(sprintf("Processing file %d of %d...\n", i, length(file_list)))
    
    # Extract Task ID from filename
    task_id <- as.numeric(gsub("Scaled_graphs_task_([0-9]+).*", "\\1", basename(f)))
    
    # Load file
    data_obj <- readRDS(f)
    graph_candidates <- data_obj$candidates
    rm(data_obj)
    
    # Go through all four graphs types
    for (type in graph_types) {
      
      g <- graph_candidates[[type]]
      if (is.null(g)) next
      
      # Extract degree vector
      d_vec <- as.numeric(igraph::degree(g))
      
      # Calculate degree fractions
      n <- vcount(g)
      d_frac <- d_vec / n
      
      # Store tiny bits of data
      meta_stats_list[[paste0(task_id, "_", type)]] <- list(
        task_id = task_id,
        type    = type,
        stats   = get_degree_stats(d_vec),
        ks      = get_ks_dist(d_vec, emp_degree_vector),
        ks_frac = get_ks_dist(d_frac, degree_frac_vector),
        f_path  = f)
      
      # Store frequency table for plotting
      counts <- table(d_vec)
      dist_plot_list[[paste0(task_id, "_", type)]] <- data.table(
        task_id   = task_id,
        degree    = as.numeric(names(counts)),
        frequency = as.numeric(counts),
        type      = type)
      
      # Store fraction frequency table for plotting
      fracs <- table(d_frac)
      frac_plot_list[[paste0(task_id, "_", type)]] <- data.table(
        task_id   = task_id,
        degree    = as.numeric(names(fracs)),
        frequency = as.numeric(fracs),
        type      = type)
    }
    
    # Explicitly remove object to be safe
    rm(data_obj)
  }
  
  
  ## Ranking graphs
  
  cat("Continue with ranking graphs according to degree similarity.\n")
  
  # Save empty lists
  final_ensemble     <- list()
  best_indices_list   <- list()
  ks_distances_list   <- list()
  centroids_list      <- list()
  
  # Go through graph types
  for (type in graph_types) {
    
    # Process them separately
    cat(paste("\n--- Processing:", type, "---\n"))
    
    # Filter through the degree stats list for this type
    type_meta <- Filter(function(x) x$type == type, meta_stats_list)
    if (length(type_meta) == 0) next
    
    # Convert to matrix for centroid calculation
    summary_matrix <- do.call(rbind, lapply(type_meta, function(x) x$stats))
    t_ids          <- sapply(type_meta, function(x) x$task_id)
    
    # Compute Model Centroid
    model_centroid <- colMeans(summary_matrix, na.rm = TRUE)
    centroids_list[[type]] <- as.data.frame(t(model_centroid)) %>% mutate(type = type)
    
    # Ranking Logic (Absolute distance from centroid)
    error_matrix <- abs(sweep(summary_matrix, 2, model_centroid, "-"))
    rank_matrix  <- apply(error_matrix, 2, rank, ties.method = "min")
    rank_sums    <- rowSums(rank_matrix)
    
    # Get top N task IDs
    best_indices      <- order(rank_sums)[1:min(n_to_keep, length(rank_sums))]
    best_task_ids <- t_ids[best_indices]
    
    # Get the actual graphs for these IDs only
    cat(sprintf("Retrieving %d graphs for %s ensemble...\n", length(best_task_ids), type))
    
    # Get ensemble of graphs
    final_ensemble[[type]] <- lapply(best_task_ids, function(tid) {
      meta_entry <- Filter(function(x) x$task_id == tid, type_meta)[[1]]
      tmp_data   <- readRDS(meta_entry$f_path)
      
      # Save in list
      list(
        origin_task_id = tid,
        graph          = tmp_data$candidates[[type]],
        rank_score     = rank_sums[which(t_ids == tid)])
    })
    
    # Log metadata for saving
    best_indices_list[[type]] <- data.frame(task_id = best_task_ids, graph_type = type)
    ks_distances_list[[type]] <- data.frame(
      task_id = t_ids,
      ks_dist = sapply(type_meta, function(x) x$ks),
      ks_frac_dist = sapply(type_meta, function(x) x$ks_frac),
      graph_type = type)
  }
  
  ## Last bit: saving and plotting!
  
  cat("\nFinalizing plots and files...\n")
  
  # Combine all distributions
  all_graphs_dist_df <- rbindlist(dist_plot_list)
  all_graphs_dist_df <- rbind(all_graphs_dist_df, emp_degree_df, fill = TRUE)
  
  # Summary for ribbon plot
  degree_summarized <- all_graphs_dist_df %>%
    group_by(type, degree) %>%
    summarise(q_25 = quantile(frequency, 0.25, na.rm = TRUE),
              mean = mean(frequency, na.rm = TRUE),
              q_75 = quantile(frequency, 0.75, na.rm = TRUE), .groups = 'drop')
  
  # Plot
  degree_plot <- ggplot(degree_summarized, aes(x = degree, y = mean, group = type, color = type, fill = type)) +
    geom_line(linewidth = 0.8) +
    geom_ribbon(aes(ymin = q_25, ymax = q_75), alpha = 0.2, color = NA) +
    scale_fill_manual(values = mypalette) +
    scale_color_manual(values = mypalette) +
    theme_bw() + 
    theme(legend.position = "bottom") +
    labs(title = paste("Degree Distribution - Seed", seed_nr),
         subtitle = "Solid line: Mean | Ribbon: IQR",
         x = "Degree", y = "Frequency")
  ggsave(file.path(plots_path, paste0("degree_plot_seed", seed_nr, ".png")), degree_plot, width = 10, height = 8)
  
  # Degree fractions!
  all_graphs_frac_df <- rbindlist(frac_plot_list)
  all_graphs_frac_df <- rbind(all_graphs_frac_df, emp_degree_frac_df, fill = TRUE)
  
  # Summary for ribbon plot
  degree_frac_summarized <- all_graphs_frac_df %>%
    group_by(type, degree) %>%
    summarise(q_25 = quantile(frequency, 0.25, na.rm = TRUE),
              mean = mean(frequency, na.rm = TRUE),
              q_75 = quantile(frequency, 0.75, na.rm = TRUE), .groups = 'drop')
  
  # Plot
  degree_frac_plot <- ggplot(degree_frac_summarized, aes(x = degree, y = mean, group = type, color = type, fill = type)) +
    geom_line(linewidth = 0.8) +
    geom_ribbon(aes(ymin = q_25, ymax = q_75), alpha = 0.2, color = NA) +
    scale_fill_manual(values = mypalette) +
    scale_color_manual(values = mypalette) +
    theme_bw() + 
    theme(legend.position = "bottom") +
    labs(title = paste("Degree (fraction) Distribution - Seed", seed_nr),
         subtitle = "Solid line: Mean | Ribbon: IQR",
         x = "Degree", y = "Frequency")
  ggsave(file.path(plots_path, paste0("degree_frac_plot_seed", seed_nr, ".png")), degree_frac_plot, width = 10, height = 8)
  
  
  # Save CSVs
  write.csv(rbindlist(best_indices_list), file.path(graphs_path, "best_graphs_indices.csv"), row.names = FALSE)
  write.csv(rbindlist(ks_distances_list), file.path(graphs_path, "ks_distances_all_sims.csv"), row.names = FALSE)
  write.csv(rbindlist(centroids_list), file.path(graphs_path, "degree_model_baselines.csv"), row.names = FALSE)
  write.csv(all_graphs_dist_df, file.path(graphs_path, "all_seeds_degree_distributions.csv"), row.names = FALSE)
  write.csv(all_graphs_frac_df, file.path(graphs_path, "all_seeds_degree_frac_distributions.csv"), row.names = FALSE)
  
  # Save Master Ensemble
  saveRDS(final_ensemble, file.path(graphs_path, "MASTER_ENSEMBLE_FOR_DM.rds"))
  
  cat("Done! Selection complete.\n")
  return(final_ensemble)
}


### Apply functions

# Colour palette
palette_paper <- c("#3A405A", "#FF8465", "#99B2DD", "#FCD2A2", "#A26F39")

# Set patterns
graphs_pattern = "Scaled_graphs_task_.*\\.rds"
summary_pattern = "scaling_results_seed_.*\\.csv"

# Set graph types
four_graphs = c("spatial", "sbm", "dcsbm", "random")

## Seed 5000
summary_check_5000 <- summarize_graphs(path = out_path,
                                       csv_pattern = summary_pattern,
                                       seed_nr = 5000) 

matches_seed5000 <- select_matches(path = out_path, 
                                   graph_types = four_graphs,
                                   graph_pattern = graphs_pattern,
                                   n_to_keep = 20,
                                   mypalette = palette_paper,
                                   seed_nr = 5000)


### Plot total degree distributions
degree_frac_summarized <- all_seeds_degree_frac_distributions %>%
  group_by(type, degree) %>%
  summarise(min = min(frequency),
            q_25 = quantile(frequency, 0.25, na.rm = TRUE),
            mean = mean(frequency, na.rm = TRUE),
            median = median(frequency, na.rm=TRUE),
            q_75 = quantile(frequency, 0.75, na.rm = TRUE), 
            max = max(frequency), .groups = 'drop')
ggplot(degree_frac_summarized, aes(x = degree, y = median, group = type, color = type, fill = type)) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = q_25, ymax = q_75), alpha = 0.2, color = NA) +
  scale_fill_manual(values = palette_paper) +
  scale_color_manual(values = palette_paper) +
  theme_bw() + 
  theme(legend.position = "bottom") +
  labs(title = "Degree fraction distribution comparison",
       subtitle = "Solid line: Mean | Ribbon: IQR",
       x = "Degree", y = "Frequency")

### Lastly, plot the KS distances by graph type

# Import
ks_distances_5000 <- read_csv("Scaling/Graphs/ks_distances_all_sims.csv") %>%
  mutate(Seed_Total = '5000')

# Plot
ks_summary <- ggplot(ks_distances_5000, aes(x = as.factor(graph_type), 
                                           y = ks_dist, 
                                           fill = as.factor(graph_type))) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = palette_paper) +
  labs(title = "Kolmogorov-Smirnov distance distributions per graph type",
       x = "Graph type",
       y = "K-S distance",
       fill = "Seed Total") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
plotpath <- file.path(out_path, "Plots/ks_distances_plot_scaled.png")
ggsave(plotpath, ks_summary, width = 10, height = 10) 

# Plot
ks_frac_summary <- ggplot(ks_distances_5000, aes(x = as.factor(graph_type), 
                                            y = ks_frac_dist, 
                                            fill = as.factor(graph_type))) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = palette_paper) +
  labs(title = "Kolmogorov-Smirnov distance distributions per graph type",
       subtitle = 'For degree fractions rather than absolutes',
       x = "Graph type",
       y = "K-S distance",
       fill = "Seed Total") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
plotpath <- file.path(out_path, "Plots/ks_distances_frac_plot_scaled.png")
ggsave(plotpath, ks_frac_summary, width = 10, height = 10) 

### End of script ###