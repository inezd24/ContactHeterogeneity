###########################################################################################################################

# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
# Created: 01 December 2025; Last edited: 03 December 2025
# Script compares different network construction algorithms

# Legend of script:

### IS FOR NEW SECTIONS IN CAPITAL ###
### Is for headings (e.g., a new function)
# is for 'small' commands (e.g., rename or merge)

##########################################################################################################################

### SET UP R ENVIRONMENT ###

# Empty environment
rm(list = ls())

# Set working directory
setwd("/scicore/home/chitnis/derkx0000/GraphComparison")

# Load required libraries
library(igraph) 
library(randnet)
library(ggplot2)
library(ggraph)
library(patchwork)
library(ggcorrplot)
library(lhs)


##########################################################################################################################

### Helper functions

# 1. Make empirical network: this is the basis of our comparisons and is zone-based data
empirical_net <- function(df){
  
  # Function generates a simple network from a dataframe with variables 'dog' and 'peer' as edgelist components
  # @param: 'df': dataframe that can be used as edgelist and has the columns 'dog' and 'peer'
  # Necessary packages: igraph
  
  cat("\nStarting with construction of empirical graph.\n")
  
  # Check that the columns dog and peer exist in the dataframe
  if (!all(c("dog", "peer") %in% names(df))) 
    stop("Input data frame must contain columns named 'dog' and 'peer'.\n")
  
  # Create full network
  graph <- graph_from_data_frame(df[,c("dog", "peer")], directed = FALSE)
  
  # If it is is not a simple (undirected, unweighted) graph, turn it into one 
  if (!is_simple(graph)){ 
    cat("Converting graph to simple graph.\n")
    graph <- igraph::simplify(graph)
  } 
  
  # Save important network parameters
  degree_g <- degree(graph)
  transitivity_g <- transitivity(graph, type = 'global')
  edge_density_g <- edge_density(graph)
  nodes_g <- vcount(graph)
  
  cat("Finished with construction of empirical graph.\n")
  
  # Return graph object
  list(
    graph = graph, 
    network_parameters = 
      list(degree_g = degree_g,
           transitivity_g = transitivity_g,
           edge_density_g = edge_density_g,
           nodes_g = nodes_g))
}

# 2. Helper function to compute the Degree Distribution
calculate_degree_distribution <- function(graph) {
  degrees <- degree(graph)
  
  # Calculate the frequency (Probability Mass Function, PMF)
  dist_table <- table(degrees)
  dist_pmf <- dist_table / sum(dist_table)
  
  # Return a dataframe/vector mapping degree value to probability
  return(as.data.frame(dist_pmf))
}

# 3. Distance statistic 1: Kolmogorov-Smirnov (KS) Distance
ks_distance <- function(synthetic_degrees, empirical_degrees) {
  
  # Convert degrees to ECDFs (Empirical Cumulative Distribution Function)
  ecdf_synth <- ecdf(synthetic_degrees)
  ecdf_emp <- ecdf(empirical_degrees)
  
  # Get all unique degree values from both sets
  all_degrees <- sort(unique(c(synthetic_degrees, empirical_degrees)))
  
  # Calculate the maximum vertical difference between the two CDFs
  ks_dist <- max(abs(ecdf_synth(all_degrees) - ecdf_emp(all_degrees)))
  return(ks_dist)
}

# 4. Distance statistic 2: Chi-Squared (χ²) Distance
chi2_distance <- function(synth_df, emp_df) {
  
  # Standardize column names
  colnames(synth_df) <- c("degree", "P_synth")
  colnames(emp_df)   <- c("degree", "P_emp")
  
  # Merge on the union of all degrees
  merged <- merge(synth_df, emp_df, by = "degree", all = TRUE)
  
  # Replace missing probabilities with 0
  merged[is.na(merged)] <- 0
  
  # Small epsilon to avoid division by zero
  eps <- 1e-10
  
  # Symmetric Chi-squared distance:
  chi2 <- sum((merged$P_synth - merged$P_emp)^2 /
                (merged$P_synth + merged$P_emp + eps))
  
  return(chi2)
}

# 5. Network construction
construct_network <- function(parameters, N_nodes_large, 
                              x_coordinates, 
                              y_coordinates, 
                              seed = NULL) {
  
  #' @param parameters A numeric vector c(m, p, s), where m is the Poisson mean for 
  #'   preferential attachment, p is the proportion of local dogs, and s is the spatial decay parameter.
  #' @param N_nodes_large The total number of nodes in the network.
  #' @param x_coordinates A numeric vector of x-coordinates for all nodes.
  #' @param y_coordinates A numeric vector of y-coordinates for all nodes.
  #' @param seed a numeric vector of seed to set
  #' @return An igraph object (a simulated network).
  
  # Parameter Initialization
  l <- parameters[1] # Mean number of peers for preferential attachment
  t <- parameters[2] # Proportion of local dogs
  k <- parameters[3] # Spatial decay parameter
  
  # Initialize the Adjacency Matrix
  adjacency <- matrix(0, nrow = N_nodes_large, ncol = N_nodes_large)
  
  # How many local dogs?
  number_of_local_dogs <- round(N_nodes_large * t)
  if (number_of_local_dogs > N_nodes_large)
    stop(paste('Warning: Number of local dogs exceeds total sample'))
  
  # Set seed (optional)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  ### Spatial edge creation
  
  # We use the raw seed for this phase
  if (!is.null(seed)) set.seed(seed)
  
  # Loop only for k < j to avoid processing the same dyad twice and ensure only one check per pair.
  for (i in 1:(N_nodes_large - 1)) { # prevents self-loops
    for (j in (i + 1):N_nodes_large) { # this prevents double dyads (e.g., 5-2 and 2-5)
      
      # Calculate Euclidean distances between the two nodes
      d_ij <- sqrt((x_coordinates[i] - x_coordinates[j])^2 +
                     (y_coordinates[i] - y_coordinates[j])^2)
      
      # Calculate the probability of attachment (r) based on distance (p_ij)
      p_ij <- exp(-k * d_ij)
      
      # Create edge if probability is higher than random
      if (runif(1) < p_ij) {
        adjacency[i, j] <- 1
        adjacency[j, i] <- 1
      }
    }
  }
  
  
  ### Preferential attachment
  
  # We use seed + 10^6 to ensure this stream never overlaps with Phase 1
  if (!is.null(seed)) set.seed(seed + 1000000)
  
  
  # Calculate degrees
  degrees <- colSums(adjacency)
  
  # Prepare the list of nodes for sampling (indices repeated by degree)
  indices_repeated <- unlist(lapply(1:length(degrees),
                                    function(j) rep(j, degrees[j])))
  
  # Check that we'd have more than 0 pref_nodes
  if(N_nodes_large > number_of_local_dogs){
    
    # Loop over preferred nodes
    pref_nodes <- 1:(N_nodes_large - number_of_local_dogs)
    for (k in pref_nodes) {
      
      # Generate the desired number of peers
      desired_peers <- rpois(1, l)
      
      # Check the size of the available pool (total existing links)
      pool_size <- length(indices_repeated)
      
      # *** STRICT CHECK IMPLEMENTATION ***
      if (desired_peers > pool_size) {
        
        # If the desired sample size exceeds the available population, stop the execution.
        stop(paste0(
          "ERROR: Parameter 'L' (mean Poisson connections) is too high. ",
          "Node ", k, " requires ", desired_peers, " connections, but only ", 
          pool_size, " unique peer links are currently available in the network. ",
          "Please reduce the value of 'L'."))
      }
      # *** End of strict check *** 
      
      number_of_peers <- desired_peers
      if (number_of_peers > 0 && length(indices_repeated) > 0) {
        
        # Sampling without replacement as according to Mirjam's python script
        # Deterministic for the chosen seed
        indices_peers <- sample(indices_repeated, number_of_peers, replace = FALSE)
        
        for (peer in indices_peers) {
          adjacency[peer, k] <- 1
          adjacency[k, peer] <- 1
        }
        
        indices_repeated <- c(indices_repeated, indices_peers)
      }
    }
    
  }
  
  ### Clean up and output
  
  # Remove self edges (if any)
  diag(adjacency) <- 0
  
  # Convert the adjacency matrix to an igraph object (H)
  H <- graph_from_adjacency_matrix(adjacency, mode = "undirected")
  
  # Add spatial coordinates as vertex attributes
  V(H)$x_coord <- x_coordinates
  V(H)$y_coord <- y_coordinates
  
  return(H)
}

# 6. Generate empty SBM graph
generate_sbm_graph <- function(num_nodes, # number of nodes
                               c,         # block membership
                               B,         # block probability matrix
                               seed = NULL){       
  
  # Force a seed
  if (!is.null(seed)) set.seed(seed)
  
  # Make an empty graph with the identical number of nodes and add the blocks
  g_sim <- make_empty_graph(n = num_nodes, directed = FALSE)
  V(g_sim)$community <- factor(c)
  
  # Create full B matrix
  full_B_matrix <- B[c, c]
  
  # Iterate through all unique pairs of nodes (i, j) where i < j
  for (i in 1:(num_nodes - 1)) {
    for (j in (i + 1):num_nodes) {
      
      # Get the community IDs of the two nodes
      
      # Get the probability of connection from the matrix
      prob <- full_B_matrix[i,j]
      
      # Draw a random number; if it's less than the probability, add an edge
      # runif(1) is deterministic based on 'seed'
      if (runif(1) < prob) { # use independent Bernoulli trials
        
        # If success: add the edges to the new graph
        g_sim <- add_edges(g_sim, c(i, j))
      }
    }
  }
  return(g_sim)
}

# 7. Helper function: Generate graph for DCSBM 
generate_dcsbm_graph <- function(num_nodes, # number of nodes
                                 Phat, # edge probability matrix
                                 seed = NULL){
  
  # Force a seed
  if (!is.null(seed)) set.seed(seed)
  
  # Make an empty graph with the identical number of nodes and add the blocks
  g_sim <- igraph::make_empty_graph(n = num_nodes, directed = FALSE)
  
  # Iterate through all unique pairs of nodes (i, j) where i < j
  for (i in 1:(num_nodes - 1)) {
    for (j in (i + 1):num_nodes) {
      
      print(paste0("Testing node", i, " and node ", j))
      
      # Get the probability of connection
      # This is the key change: we multiply the block parameter B by the
      # degree propensities of the two nodes.
      prob <- Phat[i,j]
      
      # Draw a random number; if it's less than the probability, add an edge
      # runif(1) is deterministic based on 'seed'     
      if (runif(1) < prob) { # use independent Bernoulli trials
        
        # If success: add the edges to the new graph
        g_sim <- add_edges(g_sim, c(i, j))
      }
    }
  }
  return(g_sim)
}

# 8. Function to optimize grid
grid_optimization <- function(empirical_net,
                              kappa_grid,
                              tau_grid,
                              lambda_grid,
                              square_size = 1,
                              base_seed = 42){
  
  ### Define input parameters for construction
  
  # Define degree vectors and frequency table
  empirical_degrees_vector <- degree(empirical_net)
  empirical_dist_df <- calculate_degree_distribution(empirical_net)
  
  # Set nr. of empirical nodes
  empirical_nodes <- length(empirical_degrees_vector)
  
  # Set total runs: e.g. 3 * 10 * 10 = 300
  total_runs <- length(kappa_grid) * length(tau_grid) * length(lambda_grid)
  
  # Initialize results table
  results_grid <- data.frame(
    kappa = numeric(total_runs),
    tau = numeric(total_runs),
    lambda = numeric(total_runs),
    ks_dist = numeric(total_runs),
    chi2_dist = numeric(total_runs), 
    av_degree = numeric(total_runs),
    max_degree = numeric(total_runs),
    edges = numeric(total_runs),
    edge_density = numeric(total_runs))
  run_counter <- 1
  
  # Set seed to ensure coordinates are fixed
  set.seed(base_seed)
  
  # Set square edge sizes and x_coordinates
  square_edge_size <- square_size
  N_nodes <- square_edge_size^2 * empirical_nodes
  x_coor <- square_edge_size * runif(N_nodes)
  y_coor <- square_edge_size * runif(N_nodes)
  
  
  ### Find optimized kappa, tau, lambda
  
  # Initiate grid search
  print(paste("Starting Grid Search with", total_runs, "runs..."))
  
  for (kappa_val in kappa_grid) {
    for (tau_val in tau_grid) {
      for (lambda_val in lambda_grid) {
        
        # Set distinct run seed
        run_seed <- base_seed + run_counter + 5000
        
        # Construct the Network
        temp_parameters <- c(kappa_val, tau_val, lambda_val) # [m, p, s] in construct_network code
        synthetic_graph <- construct_network(parameters = temp_parameters, 
                                             N_nodes_large = N_nodes, 
                                             x_coordinates = x_coor, 
                                             y_coordinates = y_coor, 
                                             seed = run_seed)
        
        edges <- ecount(synthetic_graph)
        density <- edge_density(synthetic_graph)
        
        # Get Synthetic Degree Distribution
        synthetic_degrees_vector <- degree(synthetic_graph)
        synthetic_dist_df <- calculate_degree_distribution(synthetic_graph)
        
        # Calculate Distances
        current_ks_dist <- ks_distance(synthetic_degrees_vector, empirical_degrees_vector)
        current_chi2_dist <- chi2_distance(synthetic_dist_df, empirical_dist_df)
        
        # Store Results
        results_grid[run_counter, ] <- c(
          kappa = kappa_val,
          tau = tau_val,
          lambda = lambda_val,
          ks_dist = current_ks_dist,
          chi2_dist = current_chi2_dist,
          av_degree = mean(synthetic_degrees_vector), 
          med_degree = median(synthetic_degrees_vector),
          max_degree = max(synthetic_degrees_vector),
          edges = edges,
          edge_density = density)
        
        # Change for seed
        run_counter <- run_counter + 1
        
        # Print results
        print(paste("Run", run_counter - 1, ": κ=", kappa_val, ", τ=", tau_val, ", λ=", lambda_val, 
                    "| KS Dist:", round(current_ks_dist, 4)))
      }
    }
  }
  print("Grid Search complete.")
  
  ### Identify Optimal Parameters of p, s, m
  
  # Find the row corresponding to the minimum KS distance
  min_ks_row <- results_grid[which.min(results_grid$ks_dist), ]
  
  # Find the row corresponding to the minimum Chi-squared distance
  min_chi2_row <- results_grid[which.min(results_grid$chi2_dist), ]
  
  # Combine and display the best results
  optimal_params <- list(
    KS_Min_Result = min_ks_row,
    Chi2_Min_Result = min_chi2_row)
  
  print("--- Optimization Results ---")
  print("Parameters that MINIMIZE Kolmogorov Distance (KS):")
  print(min_ks_row)
  print("Parameters that MINIMIZE Chi-Squared Distance (χ²):")
  print(min_chi2_row)
  
  return(results_grid)
}

# 9. ReCoN graph generator by Staudt et al. (2017)

### Final function to construct all five networks
construct_five_networks <- function(empirical_edgelist, 
                                    square_edge_size = 1,
                                    kappa_grid, 
                                    tau_grid, 
                                    lambda_grid,
                                    K_max = 20,
                                    base_seed){
  
  set.seed(base_seed)
  
  ### Empirical network
  
  # Create empirical network and degree distribution
  empirical_network <- empirical_net(empirical_edgelist)
  empirical_graph <- empirical_network$graph
  empirical_nodes <- length(degree(empirical_graph))
  empirical_edges <- ecount(empirical_graph)
  
  ### Generator 1: Mirjam's network generator
  
  # Set square edge sizes and x_coordinates
  N_nodes <- square_edge_size^2 * empirical_nodes
  
  set.seed(base_seed + 10) # Specific block for X
  x_coordinates <- square_edge_size * randomLHS(N_nodes, 1)
  
  set.seed(base_seed + 20) # Specific block for Y
  y_coordinates <- square_edge_size * randomLHS(N_nodes, 1)
  
  # Optimize grid
  network_grid <- grid_optimization(empirical_net = empirical_graph,
                                    kappa_grid = kappa_grid,
                                    tau_grid = tau_grid,
                                    lambda_grid = lambda_grid,
                                    base_seed = base_seed + 30)
  optimal_net <- network_grid[which.min(network_grid$ks_dist), ]
  
  # Get optimal grid parameters
  optimal_parameters <- c(optimal_net$lambda, 
                          optimal_net$tau,
                          optimal_net$kappa)

  # Construct network with these parameters
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
  
  ### Generators 2 and 3: Stochastic Block Model
  
  cat("\nStarting SBM and DCSBM graph construction.\n")
  
  # Turn empirical network into adjacency matrix
  A <- as.matrix(as_adjacency_matrix(empirical_graph, sparse = FALSE))
  
  # Define the range of K values to test
  set.seed(base_seed + 100) 
  BIC_emp <- randnet::LRBIC(A, Kmax = K_max, model = "both") # Bayesian Information Criterion (BIC) likelihood 
  k_hat_sbm = BIC_emp$SBM.K # 13 communities
  k_hat_dcsbm <- BIC_emp$DCSBM.K # 6 communities

  # We re-perform this if the number of optimal communities equals the number of maximum communities
  if(k_hat_sbm == K_max | k_hat_dcsbm == K_max){
    new_K_max = K_max + 1
    BIC_emp <- randnet::LRBIC(A, Kmax = new_K_max, model = "both") 
    k_hat_sbm = BIC_emp$SBM.K # 13 communities
    BIC_emp$SBM.BIC # BIC values
    k_hat_dcsbm <- BIC_emp$DCSBM.K # 6 communities
    BIC_emp$DCSBM.BIC # BIC values
  } 
  
  # Print communities 
  cat(paste0("\nSBM has ", k_hat_sbm, " communities and DCSBM has ", k_hat_dcsbm, " communities.\n"))
  
  # We use the number of communities to assign clusters and estimate the models
  set.seed(base_seed + 110) 
  if(k_hat_sbm > 1){
    communities_sbm <- reg.SSP(A,K=k_hat_sbm,lap=FALSE) # because NOT degree-corrected
    sbm_generated <- SBM.estimate(A, communities_sbm$cluster)
    
  } else {
    communities_sbm <- rep(1, empirical_nodes)
    sbm_generated <- SBM.estimate(A, communities_sbm)
  }
    
  # Repeat for DCSBM
  set.seed(base_seed + 120)
  if(k_hat_dcsbm > 1){
    communities_dcsbm <- reg.SSP(A,K=k_hat_dcsbm,lap=TRUE) # because degree-corrected
    dcsbm_generated <- DCSBM.estimate(A, communities_dcsbm$cluster)
    
  } else {
    communities_dcsbm <- rep(1, empirical_nodes)
    dcsbm_generated <- DCSBM.estimate(A, communities_dcsbm)
  }
  
  # Generate graph from model
  sbm_graph <- generate_sbm_graph(num_nodes = empirical_nodes, 
                                  c = sbm_generated$g, 
                                  B = sbm_generated$B,
                                  seed = base_seed + 130)
  
  dcsbm_graph <- generate_dcsbm_graph(num_nodes = empirical_nodes, 
                                      Phat = dcsbm_generated$Phat,
                                      seed = base_seed + 140)
  
  if(exists("sbm_graph") && exists("dcsbm_graph")){
    cat("\nBoth SBM graphs constructed. Moving on to next algorithm.\n")
  } else {
    stop("\nError: either SBM or DCSBM not constructed. Aborting process...\n")
  }
  
  ### Generator 4: ReCon
  
  #cat("\nStarting adapted ReCoN algorithm graph construction.\n")
  
  #recon_graph <- recon_generator(empirical_graph, x, steps)
  
  # if(exists("recon_graph")){
  #   cat("\nRecon_graph constructed. Moving on to next algorithm.\n")
  # } else {
  #   stop("\nError: either ReCoN graph not constructed. Aborting process...\n")
  # }
  
  ### Generator 5: random network
  
  cat("\nStarting random graph construction.\n")
  
  # Create random Erdos & Renyi graph: edges = empirical edges
  set.seed(base_seed + 300)
  random_graph = igraph::erdos.renyi.game(empirical_nodes, 
                                          empirical_edges,
                                          "gnm",
                                          FALSE, 
                                          FALSE)
  
  if(exists("random_graph")){
    cat("\nRandom graph constructed. Merging graphs in list now.\n")
  } else {
    stop("\nError: random graph not constructed. Aborting process...\n")
  }
  
  # Create list of graphs
  graph_list <- list(empirical_graph = empirical_graph, 
                     spatial_graph = spatial_graph, 
                     sbm_graph = sbm_graph, 
                     dcsbm_graph = dcsbm_graph, 
                     random_graph = random_graph)
  cat(paste0("\nUsed the following parameters for grid optimization. Kappa = ", optimal_net$kappa, 
             "; Tau = ", optimal_net$tau, "; Lambda = ", optimal_net$lambda, ".\n"))
  
  
  # Return graph list
  return(graph_list)
  
}


### Function to summarize these graphs
summarize_graph_list <- function(graph_list){
  
  # Iterate over the list and calculate metrics for each graph
  results_list <- lapply(seq_along(graph_list), function(i) {
    
    # Get each graph
    g <- graph_list[[i]]
    graph_name <- names(graph_list)[i]
    
    # Calculate degree 
    degrees <- igraph::degree(g)
    
    # Betweenness
    betweenness <- igraph::betweenness(g)
    
    # Calculate metrics
    nodes <- igraph::vcount(g)
    edges <- igraph::ecount(g)
    density <- igraph::edge_density(g)
    
    # Average Local Clustering Coefficient: calculated as the mean of the local clustering
    # coefficients of all vertices.
    clust_coeff <- igraph::transitivity(g, type = "average")
    
    # Average Path Length: 'unconnected = TRUE' ensures a result is returned even
    # if the graph is disconnected (mean of finite distances).
    path_length <- igraph::mean_distance(g, unconnected = TRUE)
    
    # Degree metrics
    avg_degree <- mean(degrees)
    med_degree <- median(degrees)
    max_degree <- max(degrees)
    
    # Eigen centrality
    eigen <- eigen_centrality(g)
    
    # Create a single-row data frame
    data.frame(
      Graph = graph_name,
      Nodes = nodes,
      Edges = edges,
      Density = density,
      Avg_Clustering_Coefficient = clust_coeff,
      Avg_Path_Length = path_length,
      Avg_Degree = avg_degree,
      Med_Degree = med_degree,
      Max_Degree = max_degree,
      Avg_Betweenness = mean(betweenness),
      Med_Betweenness = median(betweenness),
      Max_Betweenness = max(betweenness),
      stringsAsFactors = FALSE)
  })
  
  # Combine all single-row data frames into one
  network_metrics <- do.call(rbind, results_list)
  
  return(network_metrics)
  
}


### Function to plot over graph_list
plot_consistent_networks <- function(graph_list, layout_type = "fr") {
  
  # We use the first graph as the reference for layout calculation
  ref_graph <- graph_list[[1]]
  
  # Calculate the fixed layout coordinates ONCE
  cat(paste0("Calculating fixed layout using '", layout_type, "' from the reference graph...\n"))
  
  # Generate the plots by applying the fixed layout to all graphs
  plots <- lapply(seq_along(graph_list), function(i) {
    g <- graph_list[[i]]
    graph_name <- names(graph_list)[i]
    
    cat(paste0("Generating plot for: ", graph_name, "\n"))
    
    # Create the ggraph object, explicitly using the pre-calculated layout
    p <- ggraph(g, layout = layout_type) +
      geom_edge_fan(
        alpha = 0.5, 
        edge_width = 0.3, 
        edge_colour = "gray50") +
      geom_node_point(
        size = 3, 
        color = "#1e3a8a",  # Dark Blue color
        alpha = 0.8) +
      theme_graph() + # A clean theme designed for ggraph
      labs(title = paste("Network:", graph_name)) +
      
      # Make sure the plot area is square and the coordinates are fixed
      coord_fixed()
    
    return(p)
  })
  
  # Filter out any NULL results
  plots <- plots[!sapply(plots, is.null)]
  
  # Combine in one plot
  combined_plot <- patchwork::wrap_plots(plots, ncol = 3, nrow = 2) +
    patchwork::plot_annotation(
      title = "Graph Generator Comparison",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
  
  return(combined_plot)
}


### Function to compare KS distances
compare_ks_distances <- function(graph_list, country_name, location_name){
  
  # Create empty df for ks distances
  ks_distances <- data.frame(
    graph1 = as.character(), 
    graph2 = as.character(), 
    ks_dist = as.numeric())
  
  # Get names of all graphs
  graph_names <- names(graph_list)
  
  # Create all pairwise combinations (including self-comparison)
  graph_pairs <- expand.grid(graph1 = graph_names, graph2 = graph_names, stringsAsFactors = FALSE)
  
  # Get KS distance between all pairs
  calculate_distance <- function(name1, name2) {
    g1 <- graph_list[[name1]]
    g2 <- graph_list[[name2]]
    
    # Calculate degree distributions
    degree_g1 <- igraph::degree(g1)
    degree_g2 <- igraph::degree(g2)
    
    # Get the Kolmogorov-Smirnov D statistic
    ks_dist <- ks_distance(degree_g1, degree_g2)
    return(ks_dist)
  }
  
  # Get final comparison dataframe
  ks_distances <- graph_pairs
  ks_distances$ks_dist <- mapply(
    calculate_distance, 
    ks_distances$graph1, 
    ks_distances$graph2)
  
  # Create levels for correct plotting
  ks_distances <- ks_distances %>%
    mutate(graph1 = factor(graph1, levels = c("empirical_graph",
                                              "spatial_graph",
                                              "sbm_graph",
                                              "dcsbm_graph",
                                              "random_graph")),
           graph2 = factor(graph2, levels = c("empirical_graph",
                                              "spatial_graph",
                                              "sbm_graph",
                                              "dcsbm_graph",
                                              "random_graph"))) %>%
    # Add a column for the country and the location
    mutate(country = country_name, 
           location = location_name)
  
  # Correlation plot
  ks_heat <- ggplot(ks_distances, aes(x = graph1, y = graph2, fill = ks_dist)) +
    geom_tile() +
    theme_minimal() +
    theme(text = element_text(size = 20),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    labs(title = location_name)
  
  # Return df and plot
  list(plot = ks_heat,
       df = ks_distances)
  
}


### End of script ###

save.image(file = "constructNetworksFunction.RData")
