##########################################################################################################################

# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
# Created: December 2025; Last edited: September 2026
# This script acts as a repository for helper function used throughout analyses. In their current form, the remaining 
# scripts do not use 'source()' to call these functions, but call them locally. You may adapt this as desired. Please 
# note that the functions are well described here, but not elsewhere.

##########################################################################################################################

### SET UP R ENVIRONMENT ###

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

#-------------------------------------------------------------------------------------------------------------------------
# 1. Make empirical network: this function takes an edgelist and turns it into a basic igraph graph object. 
empirical_net <- function(df){
  
  # Function generates a simple network from a dataframe with variables 'dog' and 'peer' as edgelist components
  # It does not care whether interactions are weighted, repeated, or else, as it will turn the edgelist into 
  # an undirected, unweighted graph. 
  #' @param df: dataframe that can be used as edgelist and has the columns 'dog' and 'peer'. 
  #' @return list of graph and parameters
  # Necessary packages: igraph
  
  cat("\nStarting with construction of empirical graph.\n")
  
  # Check that the columns dog and peer exist in the dataframe
  if (!all(c("dog", "peer") %in% names(df))) 
    stop("Input data frame must contain columns named 'dog' and 'peer'.\n")
  
  # Create full network
  graph <- graph_from_data_frame(df[,c("dog", "peer")], directed = FALSE)
  
  # If it is is not a simple (undirected, unweighted) graph, turn it into one 
  if (!is_simple(graph)){ 
    cat("Graph is not simple. Converting graph to simple graph.\n")
    cat("Ensure that you check graph for any potential issues after simplifying.\n")
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
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# 2. Helper function to compute the Degree Distribution
calculate_degree_distribution <- function(graph) {
  
  # Function returns a dataframe with empirical degree PMF
  #' @param graph: any igraph grpah object. If you use a different graph 
  # type, make sure to convert it to an igraph object. 
  # Necessary packages: igraph

  # Save the degrees of the nodes in the graph
  degrees <- degree(graph)
  
  # Calculate the frequency of degree categories
  dist_table <- table(degrees)

  # Calculate probability
  dist_pmf <- dist_table / sum(dist_table)
  
  # Return a dataframe mapping degree value to probability
  return(as.data.frame(dist_pmf))
}
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# 3. Distance statistic 1: Kolmogorov-Smirnov (KS) Distance
ks_distance <- function(synthetic_degrees, empirical_degrees) {

  # Function calculates the Kolmogorov-Smirnov distance between two 
  # degree distributions
  #' @param synthetic_degrees: degree distribution of other graph
  #' @param empirical_degrees: degree distribution of reference graph
  # Necessary packages: igraph, stats (for ecdf)
  
  # Convert degrees to ECDFs (Empirical Cumulative Distribution Function)
  ecdf_synth <- ecdf(synthetic_degrees)
  ecdf_emp <- ecdf(empirical_degrees)
  
  # Get all unique degree values from both sets
  all_degrees <- sort(unique(c(synthetic_degrees, empirical_degrees)))
  
  # Calculate the maximum vertical difference between the two CDFs
  ks_dist <- max(abs(ecdf_synth(all_degrees) - ecdf_emp(all_degrees)))

  # Return KS distance
  return(ks_dist)
}
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# 4. Distance statistic 2: Chi-Squared (χ²) Distance
chi2_distance <- function(synth_df, emp_df) {

  # Function calculates the Chi-Squared (χ²) distance between two 
  # degree distributions. Requires output from helper function 2. 
  #' @param synth_df: df of degrees and probability of other graph
  #' @param emp_df: df of degrees and probability of reference graph
  # Necessary packages: /
  
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
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# 5. Spatially Explicit Network Construction Algorithm (SENCA)
construct_network <- function(parameters, N_nodes_large, 
                              x_coordinates, 
                              y_coordinates, 
                              seed = NULL) {
  
  # This function is a spatially-explicit network construction algorithm to generate the SENCA network
  # It is adapted from code by Laager et al. (2018), who used it to create degree-optimized networks
  #' @param parameters A numeric vector c(m, p, s), where m is the Poisson mean for 
  #'   preferential attachment, p is the proportion of local dogs, and s is the spatial decay parameter.
  #'   See Laager et al. (2018) for detailed descriptions of these parameters .
  #' @param N_nodes_large The total number of nodes in the network.
  #' @param x_coordinates A numeric vector of x-coordinates for all nodes.
  #' @param y_coordinates A numeric vector of y-coordinates for all nodes.
  #' @param seed a numeric vector of seed to set
  #' @return An igraph object (a simulated network, 'H').

  ## Step 1: initialization
  
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
  
  ## Step 2: spatial edge creation
  
  # Loop only for k < j to avoid processing the same dyad twice and ensure only one check per pair.
  for (i in 1:(N_nodes_large - 1)) { # prevents self-loops
    for (j in (i + 1):N_nodes_large) { # this prevents double dyads (e.g., 5-2 and 2-5)
      
      # Calculate Euclidean distances between the two nodes
      d_ij <- sqrt((x_coordinates[i] - x_coordinates[j])^2 +
                     (y_coordinates[i] - y_coordinates[j])^2)
      
      # Calculate the probability of attachment (p_ij) based on distance (d_ij) and spatial decay (k)
      p_ij <- exp(-k * d_ij)
      
      # Create edge if probability is higher than random
      if (runif(1) < p_ij) {
        adjacency[i, j] <- 1
        adjacency[j, i] <- 1
      }
    }
  }
  
  ## Step 3: preferential attachment
  
  # As we will use R's P-RNG again, we use seed + 10^6 to ensure 
  # this stream never overlaps with step 2
  if (!is.null(seed)) set.seed(seed + 1000000)
  
  # Calculate degrees
  degrees <- colSums(adjacency)
  
  # Prepare the list of nodes for sampling (indices repeated by degree)
  indices_repeated <- unlist(lapply(1:length(degrees),
                                    function(j) rep(j, degrees[j])))
  
  # Check: we should have more than 0 pref_nodes
  if(N_nodes_large > number_of_local_dogs){
    
    # Loop over preferred nodes
    pref_nodes <- 1:(N_nodes_large - number_of_local_dogs)
    for (k in pref_nodes) {
      
      # Generate the desired number of peers
      desired_peers <- rpois(1, l)
      
      # Check: size of the available pool (total existing links)
      pool_size <- length(indices_repeated)
      
      # Check: desired sample size should not exceed total population
      if (desired_peers > pool_size) {
        
        # If the desired sample size exceeds the available population, stop the execution.
        stop(paste0(
          "ERROR: Parameter 'L' (mean Poisson connections) is too high. ",
          "Node ", k, " requires ", desired_peers, " connections, but only ", 
          pool_size, " unique peer links are currently available in the network. ",
          "Please reduce the value of 'L'."))
      }
      
      # Set number of peers as desired nr of peers
      number_of_peers <- desired_peers

      # Check that both are not empty
      if (number_of_peers > 0 && length(indices_repeated) > 0) {
        
        # Sampling without replacement as according to the original python code
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
  
  ## Step 4: Clean up and output
  
  # Remove self edges (if any; same as simplification)
  diag(adjacency) <- 0
  
  # Convert the adjacency matrix to an undirected igraph object (H)
  H <- graph_from_adjacency_matrix(adjacency, mode = "undirected")
  
  # Add spatial coordinates as vertex attributes
  V(H)$x_coord <- x_coordinates
  V(H)$y_coord <- y_coordinates
  
  # Return graph
  return(H)
}
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# 6. Function to optimize parameter grid (this function calls function 5)
grid_optimization <- function(empirical_net,
                              kappa_grid,
                              tau_grid,
                              lambda_grid,
                              square_size = 1,
                              base_seed = 42){

  # This function performs grid optimization for the kappa, tau, and lambda parameters,
  # given a certain parameter range for each parameter and a certain square size. 
  # This function is also fully adapted from Laager et al. (2018)'s python code. 
  #' @param empirical_net: the empirical network as an igraph object
  #' @param: 'kappa_grid': parameter range for kappa
  #' @param: 'tau_grid': parameter range for tau
  #' @param: 'lambda_grid': parameter range for lambda
  #' @param: 'square_size': size of each square for grid
  #' @param: 'base_seed': seed (this is set)
  #' @return: full results table from which optimized parameter values can be obtained
  
  ## Step 1: Define input parameters for construction
  
  # Degree
  empirical_degrees_vector <- degree(empirical_net)

  # Degree probability (helper function 2 output)
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
  x_coor <- square_edge_size * runif(N_nodes) # gives point between 0 and 1
  y_coor <- square_edge_size * runif(N_nodes) # gives point between 0 and 1
  
  # Initialize results table
  results_grid <- data.frame(matrix(ncol = 10, nrow = total_runs))
  colnames(results_grid) <- c("kappa","tau","lambda","ks_dist","chi2_dist",
                              "av_degree","med_degree", "max_degree","edges","edge_density")
  
  # Set run counter
  run_counter <- 1
  
  ## Step 2: Find optimized kappa, tau, lambda
  
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
        
        # Count edges and get edge density
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
        
        # Uncomment below if you want to print the results (will significantly slow the process)
        #print(paste("Run", run_counter - 1, ": κ=", kappa_val, ", τ=", tau_val, ", λ=", lambda_val, 
        #            "| KS Dist:", round(current_ks_dist, 4)))
      }
    }
  }
  print("Grid Search complete.")
  
  ## Step 3: Identify Optimal Parameters of p, s, m
  
  # Find the row corresponding to the minimum KS distance
  min_ks_row <- results_grid[which.min(results_grid$ks_dist), ]
  
  # Find the row corresponding to the minimum Chi-squared distance
  min_chi2_row <- results_grid[which.min(results_grid$chi2_dist), ]
  
  # Combine and display the best results > choose the parameters based on these results
  # High preference given to the KS distance over the X2 distance in manuscript
  optimal_params <- list(
    KS_Min_Result = min_ks_row,
    Chi2_Min_Result = min_chi2_row)
  
  print("--- Optimization Results ---")
  print("Parameters that MINIMIZE Kolmogorov Distance (KS):")
  print(min_ks_row)
  print("Parameters that MINIMIZE Chi-Squared Distance (χ²):")
  print(min_chi2_row)
  
  # Return results table
  return(results_grid)
}
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# 7. Generate graph with SBM
generate_sbm_graph <- function(num_nodes, 
                               c,     
                               B,        
                               seed = NULL){     

  # Function generates a network using a stochastic block model
  # This function already requires you to know c and B for a given empirical network
  #' @param num_nodes: number of nodes (of your empirical graph)
  #' @param c: block membership
  #' @param B: block probability matrix 
  #' @param seed: optional seed
  #' @return an igraph object ('g_sim')
  
  # Set optional seed (per function argument)
  if (!is.null(seed)) set.seed(seed)
  
  # Make an empty graph with the identical number of nodes and add the blocks
  g_sim <- igraph::make_empty_graph(n = num_nodes, directed = FALSE)
  V(g_sim)$community <- factor(c)
  
  # Create full B matrix
  full_B_matrix <- B[c, c]
  
  # Iterate through all unique pairs of nodes (i, j) where i < j
  for (i in 1:(num_nodes - 1)) {
    for (j in (i + 1):num_nodes) {
            
      # Get the probability of connection from the matrix
      prob <- full_B_matrix[i,j]
      
      # Draw a random number; if it's less than the probability, add an edge
      # runif(1) is deterministic based on 'seed'
      if (runif(1) < prob) { 
        
        # If success: add the edges to the new graph
        g_sim <- add_edges(g_sim, c(i, j))
      }
    }
  }
  return(g_sim)
}
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# 8. Generate graph with DCSBM
generate_dcsbm_graph <- function(num_nodes, 
                                 Phat, 
                                 seed = NULL){

  # Function generates a network using a degree-corrected stochastic block model
  # This function already requires you to know Phat for a given empirical network
  #' @param num_nodes: number of nodes (of your empirical graph)
  #' @param Phat: edge probability matrix
  #' @param seed: optional seed
  #' @return an igraph object ('g_sim')
  
  # Force a seed
  if (!is.null(seed)) set.seed(seed)
  
  # Make an empty graph with the identical number of nodes and add the blocks
  g_sim <- igraph::make_empty_graph(n = num_nodes, directed = FALSE)
  
  # Iterate through all unique pairs of nodes (i, j) where i < j
  for (i in 1:(num_nodes - 1)) {
    for (j in (i + 1):num_nodes) {
      
      # Uncomment this if you want to print every node combination
      # print(paste0("Testing node", i, " and node ", j)) 
      
      # Get the probability of connection
      prob <- Phat[i,j]
      
      # Draw a random number; if it's less than the probability, add an edge
      # runif(1) is deterministic based on 'seed'     
      if (runif(1) < prob) { 
        
        # If success: add the edges to the new graph
        g_sim <- add_edges(g_sim, c(i, j))
      }
    }
  }
  return(g_sim)
}
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# 9. Pipeline to construct all networks
construct_four_networks <- function(empirical_edgelist, 
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
  
  ## Step 3: Generate Erdos-Renyi (random) graph
  
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

  ## Step 4: list and return results 
  
  # Create list of graphs, including empirical graph
  graph_list <- list(empirical_graph = empirical_graph, 
                     spatial_graph = spatial_graph, 
                     sbm_graph = sbm_graph, 
                     dcsbm_graph = dcsbm_graph, 
                     random_graph = random_graph)
  # Report the parameters used for SENCA
  cat(paste0("\nUsed the following parameters for graph optimization. Kappa = ", optimal_net$kappa, 
             "; Tau = ", optimal_net$tau, "; Lambda = ", optimal_net$lambda, ".\n"))
  
  # Return graph list, including the optimal parameters for SENCA and communities for SBM/DCSBM
  list(graph_list = graph_list,
       best_params = optimal_net,
       communities = communities)
  
}
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# Helper function 10: Joint degree sequence (without correction)
calculate_joint_degree_sequence <- function(graph) {

  # This function calculates the joint degree sequence of a given graph,
  # with edge degree and triangle degree
  #' @param graph: an igraph graph object
  #' @return a list dataframe with nr of nodes, triangle count and independent-stub count

  # Vertex count
  N <- vcount(graph)
  
  # count_triangles() returns per-node triangle counts directly.
  # Note: sum(t_vector) / 3 gives total triangles in the graph (each triangle
  # is counted once per node, so divide by 3 at the graph level only).
  t_vector <- count_triangles(graph)
  total_degrees <- degree(graph)

  # Number of independent-edge stubs
  # Each triangle contributes 2 edges to a node's total degree (one to each
  # partner), so s_i = total_degree - 2 * t_i.
  s_vector <- total_degrees - 2 * t_vector
  
  # Negative s_i means edges are shared across multiple triangles, violating
  # the Newman model's edge-disjoint assumption. Clamp to 0 and warn.
  if (any(s_vector < 0)) {
    n_affected <- sum(s_vector < 0)
    mean_triangles_per_edge <- sum(t_vector) / (3 * ecount(graph))
    warning(sprintf(
      "%d / %d nodes (%.1f%%) have negative independent-edge stubs.\n  Mean triangles per edge: %.1f (Newman model assumes ~1).\n  
      The graph's clique structure violates the edge-disjoint triangle assumption;\n  consider an alternative null model (ERGM, degree-preserving rewiring).",
      n_affected, N, 100 * n_affected / N, mean_triangles_per_edge
    ))
    s_vector <- pmax(s_vector, 0L)
    
    # Ensure even sum for valid configuration model
    if (sum(s_vector) %% 2 != 0) {
      s_vector[which.max(s_vector)] <- s_vector[which.max(s_vector)] + 1
    }
  }
  
  # Return dataframe with results 
  return(data.frame(
    node = seq_len(N),
    s_i  = s_vector, # number of independent-edge stubs
    t_i  = t_vector # number of triangles the node participates in.
  ))
}
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# Helper function 11: create a NCRG graph
random_clustered_graph <- function(joint_degree_sequence,
                                   simplify_graph = TRUE,
                                   seed = NULL) {

  # This function creates an NCRG from the joint degree sequence 
  #' @param joint_degree_sequence: output from helper function 1
  #' @param simplify_graph: whether to simplify or not. Default = TRUE
  #' @param seed: whether a seed is used or not. Default = NULL. 
  #' @return an igraph graph object  
  
  # Set optional seed
  if (!is.null(seed)) set.seed(seed)
  
  ## Step 1: Normalise input to a list of c(s_i, t_i) vectors

  # Turn joint degree sequence into matrix
  if (is.matrix(joint_degree_sequence) || is.data.frame(joint_degree_sequence)) {
    m <- as.matrix(joint_degree_sequence)

    # Check if only 2 columns
    if (ncol(m) != 2L)
      stop("Matrix/data.frame input must have exactly 2 columns: (independent_stubs, triangle_stubs).")
    joint_degree_sequence <- lapply(seq_len(nrow(m)), function(i) unname(m[i, ]))
  }
  
  # Get number of nodes
  N <- length(joint_degree_sequence)
  if (N == 0L) return(make_empty_graph(n = 0L, directed = FALSE))
  
  ## Step 2: Build stub pools and check vectors
  s_vector <- vapply(joint_degree_sequence, `[[`, numeric(1), 1)
  t_vector <- vapply(joint_degree_sequence, `[[`, numeric(1), 2)
  
  # Check for NA values
  if (anyNA(s_vector) || anyNA(t_vector))
    stop("Degree sequence contains NA or NaN values.")

  # Check that all degrees are non-negative
  if (any(s_vector < 0 | t_vector < 0 |
          s_vector != floor(s_vector) | t_vector != floor(t_vector)))
    stop("All degrees must be non-negative integers.")
  
  # Check independent-edge degrees are even
  if (sum(s_vector) %% 2L != 0L)
    stop("Sum of independent-edge degrees must be even.")

  # Check triangles are divisible by 3
  if (sum(t_vector) %% 3L != 0L)
    stop("Sum of triangle degrees must be divisible by 3.")
  
  # Build and shuffle stub pools 
  ilist <- rep(seq_len(N), times = s_vector)
  tlist <- rep(seq_len(N), times = t_vector)
  if (length(ilist) > 0) ilist <- sample(ilist)
  if (length(tlist) > 0) tlist <- sample(tlist)
  
  ## Step 3: Build edge list 
  edges_to_add <- integer(0)
  
  # Independent edges: 
  # The consecutive pairs (ilist[1], ilist[2], etc.) after shuffling form a random 
  # perfect matching on the stubs (because of how igraph reads it)
  if (length(ilist) > 0) {
    edges_to_add <- ilist
  }
  
  # Triangle edges: 
  # Reshape tlist into a 3 x (T/3) matrix; each column is one triangle (n1, n2, n3).
  if (length(tlist) > 0) {
    triplets <- matrix(tlist, nrow = 3)
    n1 <- triplets[1, ]
    n2 <- triplets[2, ]
    n3 <- triplets[3, ]
    
    # Bind each pair of rows and concatenate
    t_edges <- c(rbind(n1, n2), rbind(n1, n3), rbind(n2, n3))

    # Add edges
    edges_to_add <- c(edges_to_add, t_edges)
  }
  
  ## Step 4: construct graph 
  G <- make_empty_graph(n = N, directed = FALSE)

  # Check if we need to add edges
  if (length(edges_to_add) > 0) {
    G <- add_edges(G, edges_to_add)
  }
  
  # Optional: simplify (based on function argument)
  if (simplify_graph) {
    G <- igraph::simplify(G, remove.multiple = TRUE, remove.loops = TRUE)
  }
  
  # Return graph
  return(G)
}
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# Helper function 12: check the fit
diagnose_newman_fit <- function(graph, jds = NULL) {

   # This checks whether the NCRG is a good representation of the empirical network. This
   # function is purely there to examine the Newman fit, not to generate files. 
  #' @param graph: NCRG from helper function 2
  #' @param jds: joint degree sequence of graph (helper function 1)

  # Check if JDS is included. Else, calculate
  if (is.null(jds)) jds <- calculate_joint_degree_sequence(graph)
  
  # Important diagnostics. Names should be informative. 
  total_triangles     <- sum(jds$t_i) / 3
  implied_edges       <- total_triangles * 3
  actual_edges        <- ecount(graph)
  triangles_per_edge  <- implied_edges / actual_edges
  pct_negative        <- 100 * mean(jds$s_i == 0 & degree(graph) > 0)
  largest_clique_size <- length(largest_cliques(graph)[[1]])
  
  # Print diagnostics for evaluation of suitability
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
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# Helper function 13: summarize graphs
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
# Helper function 14: select matches
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
# Helper function 15: compute Mahalanobis distance
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

### End of script ###


