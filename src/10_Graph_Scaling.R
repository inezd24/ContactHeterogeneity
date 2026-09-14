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
LOCAL_ROOT_DIR <- "/scicore/home/chitnis/derkx0000/GraphComparison"
setwd(LOCAL_ROOT_DIR)

# Output folders
ifelse(!dir.exists(file.path(LOCAL_ROOT_DIR, "Scaling")),
       dir.create(file.path(LOCAL_ROOT_DIR, "Scaling")), FALSE)
out_path <- file.path(LOCAL_ROOT_DIR, "Scaling")
ifelse(!dir.exists(file.path(out_path, "Graphs")),
       dir.create(file.path(out_path, "Graphs")), FALSE)
out_path_graphs <- file.path(out_path, "Graphs")
ifelse(!dir.exists(file.path(out_path, "Plots")),
       dir.create(file.path(out_path, "Plots")), FALSE)
plot_path <- file.path(out_path, "Plots")

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

# 5. Spatial network construction
construct_network <- function(parameters, 
                              N_nodes_large, 
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
                              nodes, 
                              x_coor, 
                              y_coor,
                              empirical_degrees, 
                              empirial_dist,
                              base_seed){

  
  # Set total runs: e.g. 3 * 10 * 10 = 300
  total_runs <- length(kappa_grid) * length(tau_grid) * length(lambda_grid)
  
  # Initialize results table
  results_grid <- data.frame(matrix(ncol = 10, nrow = total_runs))
  colnames(results_grid) <- c("kappa","tau","lambda","ks_dist","chi2_dist",
                              "av_degree","med_degree", "max_degree","edges","edge_density")
  
  # Set run counter
  run_counter <- 1
  
  # Set seed to ensure coordinates are fixed
  set.seed(base_seed)
  
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
                                             N_nodes_large = nodes, 
                                             x_coordinates = x_coor, 
                                             y_coordinates = y_coor, 
                                             seed = run_seed)
        
        # Calculate edges and edge density
        edges <- ecount(synthetic_graph)
        density <- edge_density(synthetic_graph)
        
        # Get Synthetic Degree Distribution
        synthetic_degrees_vector <- degree(synthetic_graph)
        synthetic_dist_df <- calculate_degree_distribution(synthetic_graph)
        
        # Calculate Distances
        current_ks_dist <- ks_distance(synthetic_degrees_vector, empirical_degrees)
        current_chi2_dist <- chi2_distance(synthetic_dist_df, empirical_dist)
        
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

# 9. Function to generate adjacency matrix with dcsbm
generate_scaled_dcsbm_network <- function(num_nodes_new, c_orig, B, theta_orig, seed = NULL){       
  
  # Force a seed
  if (!is.null(seed)) set.seed(seed)
  
  ### Generate new community structure
  
  # Calculate community proportions from the original graph
  prop_orig <- prop.table(table(c_orig))
  print(prop_orig)
  
  # Generate a new community vector (c_new) for the new nodes
  c_new <- sample(x = names(prop_orig), size = num_nodes_new, replace = TRUE, prob = prop_orig)
  print(table(c_new))
  
  # Generate a new theta vector (theta_new) by sampling with replacement
  theta_new <- sample(x = theta_orig, size = num_nodes_new, replace = TRUE)
  
  
  ### Construct the new, scaled Phat matrix (Phat_scaled)
  
  # Create a matrix of theta products using outer() for efficiency
  theta_prod_matrix <- outer(theta_new, theta_new)
  
  # Use matrix subsetting to create the full B matrix for the new nodes
  # We convert the new community vector to integers for correct indexing
  full_B_matrix <- B[as.integer(c_new), as.integer(c_new)]
  
  # Calculate the scaled Phat matrix
  phat_scaled <- theta_prod_matrix * full_B_matrix
  
  ### Generate the final adjacency matrix from the new Phat matrix.
  
  # Generate a matrix of random uniform numbers
  random_matrix <- matrix(runif(num_nodes_new * num_nodes_new), nrow = num_nodes_new)
  
  # Create a new adjacency matrix based on the probabilities
  adj_matrix <- ifelse(random_matrix < phat_scaled, 1, 0)
  
  # Symmetrize the matrix
  adj_matrix[lower.tri(adj_matrix, diag = TRUE)] <- 0
  adj_matrix <- adj_matrix + t(adj_matrix)
  
  return(adj_matrix)
}

# 10. Function to generate adjacency matrix with sbm
generate_scaled_sbm_network <- function(num_nodes_new, c_orig, B, seed = NULL){       
  
  # Force a seed
  if (!is.null(seed)) set.seed(seed)
  
  ### Generate new community structure
  
  # Calculate community proportions from the original graph
  prop_orig <- prop.table(table(c_orig))
  print(prop_orig)
  
  # Generate a new community vector (c_new) for the new nodes
  c_new <- sample(x = names(prop_orig), size = num_nodes_new, replace = TRUE, prob = prop_orig)
  print(table(c_new))
  
  # Use matrix subsetting to create the full B matrix for the new nodes
  # We convert the new community vector to integers for correct indexing
  full_B_matrix <- B[as.integer(c_new), as.integer(c_new)]
  
  
  ### Generate the final adjacency matrix from the new Phat matrix.
  
  # Generate a matrix of random uniform numbers
  random_matrix <- matrix(runif(num_nodes_new * num_nodes_new), nrow = num_nodes_new)
  
  # Create a new adjacency matrix based on the probabilities
  adj_matrix <- ifelse(random_matrix < full_B_matrix, 1, 0)
  
  # Symmetrize the matrix
  adj_matrix[lower.tri(adj_matrix, diag = TRUE)] <- 0
  adj_matrix <- adj_matrix + t(adj_matrix)
  
  return(adj_matrix)
}

# 11. Formula to calculate edges from density and nodes
density_formula <- function(D, n){
  
  # Given that D = (2*E) / (n*(n-1))
  E = (D*n*(n-1))/2
  return(E)
}

# 12. Function to scale all networks
scale_networks <- function(empirical_edgelist, 
                           square_edge_size = 1, # default was 1, but we want to scale it to a larger size
                           kappa_grid, 
                           tau_grid, 
                           lambda_grid,
                           K_max = 20,
                           base_seed){
  
  # This function takes an empirical network and scales it using four different scaling algorithms.
  # It requires the following arguments:
  #' @param empirical_edgelist A dataframe with an unweighted edgelist (2 columns)
  #' @param square_edge_size The population density/size scaler. Default is 1 (no scaling)
  #' @param kappa_grid Grid to determine kappa for the spatial network 
  #' @param tau_grid Grid to determine tau for the spatial network
  #' @param lambda_grid Grid to determine lambda for the spatial network
  #' @param K_max Maximum number of communities allowed in generative models. Default is 20. 
  #' @param base_seed The base seed from the shell script for reproducibility of findings
  
  
  ### Empirical network
  
  # Create empirical network and degree distribution
  empirical_network <- empirical_net(empirical_edgelist)
  empirical_graph <- empirical_network$graph
  empirical_nodes <- length(degree(empirical_graph))
  empirical_edges <- ecount(empirical_graph)
  empirical_density <- edge_density(empirical_graph)
  empirical_degrees_vector <- degree(empirical_graph)
  empirical_dist_df <- calculate_degree_distribution(empirical_graph)
  
  
  ### Set up scaling parameters
  
  # Set square edge sizes and x_coordinates
  N_nodes <- square_edge_size^2 * empirical_nodes
  
  
  ### Generator 1: spatial model
  
  set.seed(base_seed + 10) # Specific block for X
  x_coordinates <- square_edge_size * randomLHS(N_nodes, 1)
  
  set.seed(base_seed + 20) # Specific block for Y
  y_coordinates <- square_edge_size * randomLHS(N_nodes, 1)
  
  # Optimize grid
  network_grid <- grid_optimization(empirical_net = empirical_graph,
                                    kappa_grid = kappa_grid,
                                    tau_grid = tau_grid,
                                    lambda_grid = lambda_grid,
                                    nodes = N_nodes, 
                                    x_coor = x_coordinates, 
                                    y_coor = y_coordinates,
                                    empirical_degrees = empirical_degrees_vector, 
                                    empirial_dist = empirical_dist_df,
                                    base_seed = base_seed + 30)
  
  # Choose grid that yields lowest ks distance with empirical network
  optimal_net <- network_grid[which.min(network_grid$ks_dist), ]
  
  # Obtain optimal grid parameters
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
  
  ### Generators 2 and 3: Stochastic Block Models
  
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
  
  # Generate graphs from each generative model
  
  #SBM
  sbm_A <- generate_scaled_sbm_network(num_nodes_new = N_nodes, 
                                           c_orig = sbm_generated$g, 
                                           B = sbm_generated$B,
                                           seed = base_seed + 130)
  sbm_graph <- graph_from_adjacency_matrix(sbm_A, mode = "undirected")
  
  #DCSBM
  dcsbm_A <- generate_scaled_dcsbm_network(num_nodes_new = N_nodes, 
                                               c_orig = dcsbm_generated$g, 
                                               B = dcsbm_generated$B,
                                               theta_orig = dcsbm_generated$Psi,
                                               seed = base_seed + 140)
  dcsbm_graph <- graph_from_adjacency_matrix(dcsbm_A, mode = "undirected")
  
  
  # Give error if graphs failed to be made
  if(exists("sbm_graph") && exists("dcsbm_graph")){
    cat("\nBoth SBM graphs constructed. Moving on to next algorithm.\n")
  } else {
    stop("\nError: either SBM or DCSBM not constructed. Aborting process...\n")
  }

  
  ### Generator 4: random network
  
  cat("\nStarting random graph construction.\n")
  
  # Use density of empirical network to calculate new number of edges
  N_edges = round(density_formula(empirical_density, N_nodes),0)

  # Create random Erdos & Renyi graph based on the same edge density as the empirical graph
  set.seed(base_seed + 400)
  random_graph = igraph::sample_gnm(N_nodes, 
                                    N_edges,
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
  list(graph_list = graph_list,
       best_params = optimal_net,
       communities = communities)
  
}


#### EXECUTION ####


## 1. Setting up arguments
args <- commandArgs(trailingOnly = TRUE)

# Set arguments
TASK_ID <- as.numeric(args[1])
SQUARE_EDGE <- as.numeric(args[2])
K_MAX <- as.numeric(args[3])
MAX_SEED <- as.numeric(args[4])

# Seeding
BASE_SEED <- 1000 + (TASK_ID * 10000) # Give each task 10,000 "room"


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
  results_from_func <- scale_networks(empirical_edgelist = ServerData$dog_contact_zone_1[,c(2,3)], 
                                      square_edge_size = SQUARE_EDGE, # default was 1, but we want to scale it to a larger size
                                      kappa_grid= seq(25, 35, by = 1), 
                                      tau_grid = seq(0.5, 0.95, by = 0.05), 
                                      lambda_grid = seq(22, 30, by = 1),
                                      K_max = K_MAX,
                                      base_seed = BASE_SEED) 
  
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
      Seed = BASE_SEED,
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
      stringsAsFactors = FALSE)
    replicate_results[[length(replicate_results) + 1]] <- metrics
  }
}


## 4. Bind summary statistics for evaluating the graph generators

# Bind results (df contains all replicates)
if (length(replicate_results) > 0) {
  output_df <- do.call(rbind, replicate_results)
  
  # Create the full path string
  file_name <- paste0("scaling_results_seed_", BASE_SEED, "_seed_", MAX_SEED, ".csv")
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
                   SD_optT = sd(Opt_Tau)) %>%
  ungroup()

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
  output_rds_path <- file.path(out_path_graphs, paste0("Scaled_graphs_task_", TASK_ID, "_seed_", MAX_SEED, ".rds"))
  saveRDS(vetted_output, output_rds_path)
  cat(paste("Scaled graphs for task", TASK_ID, "saved successfully.\n"))
  
} else {
  
  # WARNING: If stability fails, the graphs are NOT saved. 
  cat("CRITICAL WARNING: Stochastic leakage detected (SD > 0).\n")
  cat("The replicates are NOT identical. Extraction aborted.\n")
}



### End of script