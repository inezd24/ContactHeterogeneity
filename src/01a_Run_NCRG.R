###########################################################################################################################

# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
# Created: May 2026
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
ifelse(!dir.exists(file.path(LOCAL_ROOT_DIR, "Net_Sens")),
       dir.create(file.path(LOCAL_ROOT_DIR, "Net_Sens")), FALSE)
net_dir <- file.path(LOCAL_ROOT_DIR, "Net_Sens")
out_path_results <- file.path(net_dir, "Seed_5000")
out_path_graphs <- file.path(net_dir, "Graphs/Seed_5000")
plot_path <- file.path(LOCAL_ROOT_DIR, "Plots")

##########################################################################################################################

# Write arguments

args <- commandArgs(trailingOnly = TRUE)
task_id  <- as.numeric(args[1])
base_seed <- 1000 + (task_id * 10000)
max_seed  <- 5000

# --- Load existing graph list --------------------------------------------
graphs_path <- file.path(out_path_graphs, 
                         paste0("vetted_graphs_task_", task_id, 
                                "_seed_", max_seed, ".rds"))

if (!file.exists(graphs_path)) {
  stop(sprintf("File not found: %s", graphs_path))
}

vetted <- readRDS(graphs_path)
empirical_graph <- vetted$anchor

cat(sprintf("Task %d: loaded graph list. Candidates: %s\n",
            task_id, paste(names(vetted$candidates), collapse = ", ")))

# --- Skip if already patched (idempotent) --------------------------------
if ("newclust_graph" %in% names(vetted$candidates)) {
  cat("Newman graph already present — nothing to do.\n")
  quit(save = "no", status = 0)
}

# ------------ Helper functions ------------------------------------------------

# Function to load empirical graph
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

# Joint degree sequence (without correction)
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
    
    # Ensure even sum for valid configuration model
    if (sum(s_vector) %% 2 != 0) {
      s_vector[which.max(s_vector)] <- s_vector[which.max(s_vector)] + 1
    }
  }
  
  return(data.frame(
    node = seq_len(N),
    s_i  = s_vector,
    t_i  = t_vector
  ))
}

# Joint degree sequence WITH DEGREE CONSERVATION AND CLIQUE CORRECTION
calculate_joint_degree_sequence_correction <- function(graph) {
  N <- vcount(graph)
  
  # 1. Gather true raw metrics
  t_vector <- count_triangles(graph)
  total_degrees <- degree(graph)
  
  # 2. Vectorized decomposition with edge conservation
  # Initialize vectors to match true degree bounds
  s_vector <- integer(N)
  t_corrected <- integer(N)
  
  for (i in seq_len(N)) {
    ki <- total_degrees[i]
    Ti <- t_vector[i]
    
    if (ki - (2 * Ti) >= 0) {
      # No structural overflow: Node can accommodate all its triangles safely
      t_corrected[i] <- Ti
      s_vector[i]     <- ki - (2 * Ti)
    } else {
      # Clique/Overlap Correction: Triangles exceed degree budget.
      # Max out triangles safely to conserve true edge count constraint (ki = s_i + 2t_i)
      t_corrected[i] <- floor(ki / 2)
      s_vector[i]     <- ki %% 2
    }
  }
  
  # 3. Handle Configuration Handshakes (Parity Requirements)
  
  # Parity Rule A: Sum of independent stubs (s_i) must be even
  if (sum(s_vector) %% 2 != 0) {
    # Adjust a non-zero single stub to fix matching parity safely
    pos_s <- which(s_vector > 0)
    if (length(pos_s) > 0) {
      idx <- pos_s[1]
      s_vector[idx] <- s_vector[idx] - 1
    } else {
      # If all s_vector elements are 0, swap 1 triangle corner for a matching stub pair
      pos_t <- which(t_corrected > 0)
      if (length(pos_t) > 0) {
        idx <- pos_t[1]
        t_corrected[idx] <- t_corrected[idx] - 1
        s_vector[idx]    <- s_vector[idx] + 2
      }
    }
  }
  
  # Parity Rule B: Sum of triangle stubs (t_i) must be divisible by 3
  rem_t <- sum(t_corrected) %% 3
  if (rem_t != 0) {
    # Step down triangle stubs by the remainder, transferring budget safely to single edges
    for (step in seq_len(rem_t)) {
      pos_t <- which(t_corrected > 0)
      if (length(pos_t) > 0) {
        idx <- pos_t[1]
        t_corrected[idx] <- t_corrected[idx] - 1
        s_vector[idx]    <- s_vector[idx] + 2
      }
    }
  }
  
  return(data.frame(
    node = seq_len(N),
    s_i  = as.integer(s_vector),
    t_i  = as.integer(t_corrected)
  ))
}

# NCRG graph
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

# Function to determine fit:
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


# --------- Load data ----------------------------------------------------------

ServerData <- readRDS("~/BaseData/Chad/ServerData.rds")

# Build empirical graph
empirical_network <- empirical_net(ServerData$dog_contact_zone_1)
empirical_graph   <- empirical_network$graph

# Free up space
rm(ServerData)
rm(empirical_network)

# --- Generate Newman graph -----------------------------------------------

# Apply function
diagnose_newman_fit(empirical_graph)

joint_degrees <- calculate_joint_degree_sequence(empirical_graph)
joint_degrees_corrected <- calculate_joint_degree_sequence_correction(empirical_graph)

ncr_graph <- random_clustered_graph(
  joint_degree_sequence = joint_degrees[, c("s_i", "t_i")],
  simplify_graph        = TRUE,
  seed                  = base_seed + 250
)

# Correction
ncr_graph_corrected <- random_clustered_graph(
  joint_degree_sequence = joint_degrees_corrected[, c("s_i", "t_i")],
  simplify_graph        = TRUE,
  seed                  = base_seed + 250
)

cat(sprintf("Newman graph: %d nodes, %d edges.\n", 
            vcount(ncr_graph), ecount(ncr_graph)))

# --- Patch into candidates and resave ------------------------------------
vetted$candidates[["newclust_graph"]] <- ncr_graph
saveRDS(vetted, graphs_path)
cat(sprintf("Saved updated graph list to %s\n", graphs_path))

# --- Append metrics to existing results csv ------------------------------
results_path <- file.path(out_path_results,
                          paste0("results_seed_", base_seed, 
                                 "_seed_", max_seed, ".csv"))

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
  
  write.csv(rbind(existing_df, new_row), results_path, row.names = FALSE)
  cat(sprintf("Metrics appended to %s\n", results_path))
  
} else {
  cat("No existing results CSV found — skipping metrics append.\n")
}