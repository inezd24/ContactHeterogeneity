###########################################################################################################################

# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
# Created: 01 December 2025; Last edited: 17 December 2025
# Script compares different network construction algorithms

# Legend of script:

### IS FOR NEW SECTIONS IN CAPITAL ###
### Is for headings (e.g., a new function)
# is for 'small' commands (e.g., rename or merge)

##########################################################################################################################

### SET UP R ENVIRONMENT ###

# Empty list
rm(list = ls())

.libPaths(c("/scicore/home/chitnis/derkx0000/R/ubuntu/4.2.1-foss-2022a", .libPaths()))
print(.libPaths())
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

# Output folders
ifelse(!dir.exists(file.path(LOCAL_ROOT_DIR, "Net_Sens")),
       dir.create(file.path(LOCAL_ROOT_DIR, "Net_Sens")), FALSE)
out_path <- file.path(LOCAL_ROOT_DIR, "Net_Sens")
out_path_graphs <- file.path(out_path, "Graphs")
plot_path <- file.path(LOCAL_ROOT_DIR, "Plots")

##########################################################################################################################


### Loading helper functions


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
      
      # print(paste0("Testing node", i, " and node ", j))
      
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
  
  # Set seed to ensure coordinates are fixed
  set.seed(base_seed)
  
  # Set square edge sizes and x_coordinates
  square_edge_size <- square_size
  N_nodes <- square_edge_size^2 * empirical_nodes
  x_coor <- square_edge_size * runif(N_nodes)
  y_coor <- square_edge_size * runif(N_nodes)
  
  # Initialize results table
  results_grid <- data.frame(matrix(ncol = 10, nrow = total_runs))
  colnames(results_grid) <- c("kappa","tau","lambda","ks_dist","chi2_dist",
                              "av_degree","med_degree", "max_degree","edges","edge_density")
  
  # Set run counter
  run_counter <- 1
  
  
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
        
        # Print results > not for now 
        #print(paste("Run", run_counter - 1, ": κ=", kappa_val, ", τ=", tau_val, ", λ=", lambda_val, 
        #            "| KS Dist:", round(current_ks_dist, 4)))
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

# 09. Joint degree sequence
#' Compute the joint degree sequence of a graph for the Newman clustered model
#'
#' @param graph An undirected igraph object.
#' @return A data.frame with columns: node, s_i (independent-edge stubs),
#'   t_i (triangle stubs). Suitable for passing directly to random_clustered_graph().
#' @details t_i is the number of triangles node i participates in, via
#'   igraph::count_triangles(). s_i = total_degree - 2 * t_i, since each triangle
#'   contributes 2 edges to a node's degree. Negative s_i indicates edges shared
#'   across multiple triangles, violating the Newman model's edge-disjoint assumption.
calculate_joint_degree_sequence <- function(graph) {
  N <- vcount(graph)
  
  # count_triangles() returns per-node triangle counts directly.
  # Note: sum(t_vector) / 3 gives total triangles in the graph (each triangle
  # is counted once per node, so divide by 3 at the graph level only).
  t_vector <- count_triangles(graph)
  
  # Each triangle contributes 2 edges to a node's total degree (one to each
  # partner), so s_i = total_degree - 2 * t_i.
  total_degrees <- degree(graph)
  s_vector <- total_degrees - 2 * t_vector
  
  # Negative s_i means edges are shared across multiple triangles, violating
  # the Newman model's edge-disjoint assumption. Clamp to 0 and warn.
  if (any(s_vector < 0)) {
    n_affected <- sum(s_vector < 0)
    mean_triangles_per_edge <- sum(t_vector) / (3 * ecount(graph))
    warning(sprintf(
      "%d / %d nodes (%.1f%%) have negative independent-edge stubs.\n  Mean triangles per edge: %.1f (Newman model assumes ~1).\n  The graph's clique structure violates the edge-disjoint triangle assumption;\n  consider an alternative null model (ERGM, degree-preserving rewiring).",
      n_affected, N, 100 * n_affected / N, mean_triangles_per_edge
    ))
    s_vector <- pmax(s_vector, 0L)
  }
  
  return(data.frame(
    node = seq_len(N),
    s_i  = s_vector,
    t_i  = t_vector
  ))
}

# 10. Helper function: Newman clustering
#' Generate a Newman clustered random graph (configuration model with triangles)
#'
#' @param joint_degree_sequence A matrix, data.frame, or list of 2-element vectors
#'   c(s_i, t_i), where s_i is the number of independent edge stubs and t_i is the
#'   number of triangle stubs for node i. Columns must be in order (independent, triangle).
#'   The sum of all s_i must be even; the sum of all t_i must be divisible by 3.
#' @param simplify_graph Logical. If TRUE (default), remove self-loops and
#'   multi-edges after construction. Set to FALSE to keep the raw configuration
#'   model output.
#' @param seed Optional integer seed for reproducibility.
#' @return An undirected igraph object.
random_clustered_graph <- function(joint_degree_sequence,
                                   simplify_graph = TRUE,
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # --- Normalise input to a list of c(s_i, t_i) vectors -------------------
  if (is.matrix(joint_degree_sequence) || is.data.frame(joint_degree_sequence)) {
    m <- as.matrix(joint_degree_sequence)
    if (ncol(m) != 2L)
      stop("Matrix/data.frame input must have exactly 2 columns: (independent_stubs, triangle_stubs).")
    joint_degree_sequence <- lapply(seq_len(nrow(m)), function(i) unname(m[i, ]))
  }
  
  N <- length(joint_degree_sequence)
  if (N == 0L) return(make_empty_graph(n = 0L, directed = FALSE))
  
  # --- Build stub pools (vectorised) ---------------------------------------
  s_vector <- vapply(joint_degree_sequence, `[[`, numeric(1), 1)
  t_vector <- vapply(joint_degree_sequence, `[[`, numeric(1), 2)
  
  # Validate: NA/NaN first (comparisons on NA produce NA, not FALSE)
  if (anyNA(s_vector) || anyNA(t_vector))
    stop("Degree sequence contains NA or NaN values.")
  if (any(s_vector < 0 | t_vector < 0 |
          s_vector != floor(s_vector) | t_vector != floor(t_vector)))
    stop("All degrees must be non-negative integers.")
  
  # --- Validate stub-pool sizes --------------------------------------------
  if (sum(s_vector) %% 2L != 0L)
    stop("Sum of independent-edge degrees must be even.")
  if (sum(t_vector) %% 3L != 0L)
    stop("Sum of triangle degrees must be divisible by 3.")
  
  # --- Build and shuffle stub pools ----------------------------------------
  # rep(seq_len(N), times = ...) avoids loop-with-c() memory reallocation
  ilist <- rep(seq_len(N), times = s_vector)
  tlist <- rep(seq_len(N), times = t_vector)
  
  if (length(ilist) > 0) ilist <- sample(ilist)
  if (length(tlist) > 0) tlist <- sample(tlist)
  
  # --- Build edge list -----------------------------------------------------
  edges_to_add <- integer(0)
  
  # Independent edges: after shuffling, consecutive pairs (ilist[1], ilist[2]),
  # (ilist[3], ilist[4]), ... form a random perfect matching on the stubs.
  # igraph interprets a flat integer vector as consecutive pairs.
  if (length(ilist) > 0) {
    edges_to_add <- ilist
  }
  
  # Triangle edges: reshape tlist into a 3 x (T/3) matrix; each column is one
  # triangle (n1, n2, n3). rbind each pair of rows then c() concatenates into
  # the interleaved (n1[1],n2[1], n1[2],n2[2], ...) layout igraph expects.
  if (length(tlist) > 0) {
    triplets <- matrix(tlist, nrow = 3)
    n1 <- triplets[1, ]
    n2 <- triplets[2, ]
    n3 <- triplets[3, ]
    
    t_edges <- c(rbind(n1, n2), rbind(n1, n3), rbind(n2, n3))
    edges_to_add <- c(edges_to_add, t_edges)
  }
  
  # --- Construct graph -----------------------------------------------------
  G <- make_empty_graph(n = N, directed = FALSE)
  if (length(edges_to_add) > 0) {
    G <- add_edges(G, edges_to_add)
  }
  
  if (simplify_graph) {
    G <- igraph::simplify(G, remove.multiple = TRUE, remove.loops = TRUE)
  }
  
  return(G)
}

# 10a. extra: get output for ncrg
diagnose_newman_fit <- function(graph, jds = NULL) {
  if (is.null(jds)) jds <- calculate_joint_degree_sequence(graph)
  
  total_triangles     <- sum(jds$t_i) / 3
  implied_edges       <- total_triangles * 3
  actual_edges        <- ecount(graph)
  triangles_per_edge  <- implied_edges / actual_edges
  pct_negative        <- 100 * mean(jds$s_i == 0 & degree(graph) > 0)
  largest_clique_size <- length(largest_cliques(graph)[[1]])
  
  cat(sprintf("Nodes:                  %d\n",   vcount(graph)))
  cat(sprintf("Edges:                  %d\n",   actual_edges))
  cat(sprintf("Triangles:              %.0f\n", total_triangles))
  cat(sprintf("Implied edges (triags): %.0f\n", implied_edges))
  cat(sprintf("Mean triangles/edge:    %.1f\n", triangles_per_edge))
  cat(sprintf("Largest clique:         %d\n",   largest_clique_size))
  cat(sprintf("Nodes with s_i = 0:     %.1f%%\n", pct_negative))
  cat(sprintf("Newman model suitable:  %s\n",
              if (triangles_per_edge < 2) "Likely" else "Unlikely"))
  
  invisible(list(
    total_triangles    = total_triangles,
    triangles_per_edge = triangles_per_edge,
    largest_clique     = largest_clique_size,
    pct_si_zero        = pct_negative
  ))
}


# 11. Final integration
construct_five_networks <- function(empirical_edgelist, 
                                    square_edge_size = 1,
                                    kappa_grid, 
                                    tau_grid, 
                                    lambda_grid,
                                    K_max = 20,
                                    base_seed){
  
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
    
    cat("maximum K has been reached. Need to adjust K through K + 1")
    new_K_max = K_max + 1
    BIC_emp <- randnet::LRBIC(A, Kmax = new_K_max, model = "both") 
    k_hat_sbm = BIC_emp$SBM.K # 13 communities
    BIC_emp$SBM.BIC # BIC values
    k_hat_dcsbm <- BIC_emp$DCSBM.K # 6 communities
    BIC_emp$DCSBM.BIC # BIC values
  } 
  
  # Print communities 
  cat(paste0("\nSBM has ", k_hat_sbm, " communities and DCSBM has ", k_hat_dcsbm, " communities.\n"))
  
  # Save for later inspection
  communities <- list(sbm_com = k_hat_sbm,
                      dcsbm_com = k_hat_dcsbm)
  
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
  
  ### Generator 4: NCRG
  
  cat("\nStarting Newman clustered random graph construction.\n")
  
  # First, let's examine the graph characteristics
  diagnose_newman_fit(empirical_graph)
  
  # Make joint degree df
  joint_degrees <- calculate_joint_degree_sequence(empirical_graph)
  
  # Make graph
  ncr_graph <- random_clustered_graph(joint_degree_sequence = joint_degrees,
                                      simplify_graph = TRUE,
                                      seed = base_seed + 250)
  
  ### Generator 5: random network
  
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
  
  # Create list of graphs
  graph_list <- list(empirical_graph = empirical_graph, 
                     spatial_graph = spatial_graph, 
                     sbm_graph = sbm_graph, 
                     dcsbm_graph = dcsbm_graph, 
                     newclust_graph = ncr_graph,
                     random_graph = random_graph)
  cat(paste0("\nUsed the following parameters for grid optimization. Kappa = ", optimal_net$kappa, 
             "; Tau = ", optimal_net$tau, "; Lambda = ", optimal_net$lambda, ".\n"))
  
  
  # Return graph list
  list(graph_list = graph_list,
       best_params = optimal_net,
       communities = communities)
  
}


### Setting up arguments


## 1. Get the seed from SLURM environment
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])
base_seed <- 1000 + (task_id * 10000) # Give each task 10,000 "room"
max_seed = 5000 # adapt

## 2. Load data
ServerData <- readRDS("~/BaseData/Chad/ServerData.rds")
if (!exists("ServerData")) {
  stop("Error: The 'ServerData' object was not loaded from ServerData.rds.")
}


## 3. Runs and replicates

# Generate empty lists
replicate_results <- list()
total_graphs_list <- list()

# Create graph replicates for each seed
for (r in 1:5) {
  
  # Generate the five networks
  results_from_func <- construct_five_networks(
    empirical_edgelist = ServerData$dog_contact_zone_1,
    kappa_grid = seq(28, 30, by = 1), 
    tau_grid = seq(0.5, 0.95, by = 0.05),
    lambda_grid = seq(28, 33, by = 1),
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

# Bind results (df contains all replicates)
if (length(replicate_results) > 0) {
  output_df <- do.call(rbind, replicate_results)
  
  # Create the full path string
  file_name <- paste0("results_seed_", base_seed, "_seed_", max_seed, ".csv")
  final_out_path <- file.path(out_path, file_name)
  
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
  # This prevents "leaky" data from reaching your SEIR simulations.
  cat("CRITICAL WARNING: Stochastic leakage detected (SD > 0).\n")
  cat("The replicates are NOT identical. Extraction aborted.\n")
}
  


### End of script

