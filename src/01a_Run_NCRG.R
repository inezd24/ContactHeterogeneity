##########################################################################################################################
# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
# Created: May 2026; Last edited: September 2026
# This script serves as an addendum for script 01. It constructs the NCRG separately and adds the graphs to the original
# .rds files. 

# IMPORTANT NOTES:
# 0. 

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

# Output folders and paths
ifelse(!dir.exists(file.path(LOCAL_ROOT_DIR, "Net_Sens")),
       dir.create(file.path(LOCAL_ROOT_DIR, "Net_Sens")), FALSE)
net_dir <- file.path(LOCAL_ROOT_DIR, "Net_Sens")
out_path_results <- file.path(net_dir, "Seed_5000")
out_path_graphs <- file.path(net_dir, "Graphs/Seed_5000")
plot_path <- file.path(LOCAL_ROOT_DIR, "Plots")

##########################################################################################################################

### Helper functions

# First get functions from original helper function script
source("00_constructNetworkFunctions.R")

#-------------------------------------------------------------------------------------------------------------------------
# Helper function 1: Joint degree sequence (without correction)
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
# Helper function 2: create a NCRG graph
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
# Helper function 3: check the fit
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

##########################################################################################################################


## Step 1: set-up

# Get task id and seed
args <- commandArgs(trailingOnly = TRUE)
task_id  <- as.numeric(args[1])
base_seed <- 1000 + (task_id * 10000)
max_seed  <- 5000

# Get path of graphs generated in script 01
graphs_path <- file.path(out_path_graphs, 
                         paste0("vetted_graphs_task_", task_id, 
                                "_seed_", max_seed, ".rds"))
if (!file.exists(graphs_path)) {
  stop(sprintf("File not found: %s", graphs_path))
}

# For each task, load the four synthetic graphs and empirical graph separately 
vetted <- readRDS(graphs_path)
empirical_graph <- vetted$anchor

cat(sprintf("Task %d: loaded graph list. Candidates: %s\n",
            task_id, paste(names(vetted$candidates), collapse = ", ")))

# Check if there is already a newclust_graph as a safety measure
if ("newclust_graph" %in% names(vetted$candidates)) {
  cat("Newman graph already present — nothing to do.\n")
  quit(save = "no", status = 0)
}

## Step 2: load data

# Load serverdata. Again, change according to where your data is stored. 
ServerData <- readRDS("~/BaseData/Chad/ServerData.rds")

# Build empirical graph
empirical_network <- empirical_net(ServerData$dog_contact_zone_1)
empirical_graph   <- empirical_network$graph

# Free up space
rm(ServerData)
rm(empirical_network)

## Step 3: Generate Newman graph

# How well would Newman fit empirical graph?
diagnose_newman_fit(empirical_graph)

# Calculate and save joint degree sequence 
joint_degrees <- calculate_joint_degree_sequence(empirical_graph)

# Generate NCR graph 
ncr_graph <- random_clustered_graph(
  joint_degree_sequence = joint_degrees[, c("s_i", "t_i")],
  simplify_graph        = TRUE,
  seed                  = base_seed + 250
)

cat(sprintf("Newman graph: %d nodes, %d edges.\n", 
            vcount(ncr_graph), ecount(ncr_graph)))

## Step 4: Add graph to original vetted list

# Add to candidates list
vetted$candidates[["newclust_graph"]] <- ncr_graph

# Save rds
saveRDS(vetted, graphs_path)
cat(sprintf("Saved updated graph list to %s\n", graphs_path))

# Get path of results csv
results_path <- file.path(out_path_results,
                          paste0("results_seed_", base_seed, 
                                 "_seed_", max_seed, ".csv"))

# Read results csv
if (file.exists(results_path)) {
  existing_df <- read.csv(results_path)
  
  # Skip if already appended
  if (any(existing_df$Graph == "newclust_graph")) {
    cat("Metrics already present in CSV — skipping append.\n")
    quit(save = "no", status = 0)
  }
  
  # Pull shared fields from an existing row for this task
  ref_row <- existing_df[existing_df$Seed == base_seed &
                           existing_df$Graph == "empirical_graph", ][1, ]
  
  # Create new row
  new_row <- data.frame(
    Seed       = base_seed,
    Replicate  = ref_row$Replicate,
    Graph      = "newclust_graph",
    Opt_Lambda = ref_row$Opt_Lambda,
    Opt_Tau    = ref_row$Opt_Tau,
    Opt_Kappa  = ref_row$Opt_Kappa,
    Nodes      = vcount(ncr_graph),
    Edges      = ecount(ncr_graph),
    Density    = edge_density(ncr_graph),
    Avg_Clustering_Coefficient = transitivity(ncr_graph, type = "global"),
    Avg_Path_Length            = mean_distance(ncr_graph, directed = FALSE),
    Avg_Degree                 = mean(degree(ncr_graph)),
    Med_Degree                 = median(degree(ncr_graph)),
    Max_Degree                 = max(degree(ncr_graph)),
    Avg_Betweenness            = mean(betweenness(ncr_graph)),
    com_sbm    = ref_row$com_sbm,
    com_dcsbm  = ref_row$com_dcsbm,
    stringsAsFactors = FALSE
  )
  
  # Overwrite csv 
  write.csv(rbind(existing_df, new_row), results_path, row.names = FALSE)
  cat(sprintf("Metrics appended to %s\n", results_path))
  
} else {
  cat("No existing results CSV found — skipping metrics append.\n")
}


### END OF SCRIPT