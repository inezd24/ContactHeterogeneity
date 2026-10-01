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
LOCAL_ROOT_DIR <- "/scicore/home/chitnis/derkx0000/Derkx2026_publication"
setwd(LOCAL_ROOT_DIR)

# Output folders and paths
ifelse(!dir.exists(file.path(LOCAL_ROOT_DIR, "Net_Sens")),
       dir.create(file.path(LOCAL_ROOT_DIR, "Net_Sens")), FALSE)
net_dir <- file.path(LOCAL_ROOT_DIR, "Net_Sens")
out_path_results <- file.path(net_dir, "Seed_5000")
out_path_graphs <- file.path(net_dir, "Graphs/Seed_5000")
plot_path <- file.path(LOCAL_ROOT_DIR, "Plots")

##########################################################################################################################

### Load all helper functions (1-15)

# All functions are stored in this R script and corresponding .Rdata file: "HelperFunctions.RData"
source("00_constructNetworksFunction.R")

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