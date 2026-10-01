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
LOCAL_ROOT_DIR <- "/scicore/home/chitnis/derkx0000/Derkx2026_publication"
setwd(LOCAL_ROOT_DIR)

# Output folders and paths. 
out_path <- file.path(LOCAL_ROOT_DIR, "Net_Sens")
ifelse(!dir.exists(file.path(out_path, "Plots")),
       dir.create(file.path(out_path, "Plots")), FALSE)

##########################################################################################################################

### Load all helper functions (1-15)

# All functions are stored in this R script and corresponding .Rdata file: "HelperFunctions.RData"
source("00_constructNetworksFunction.R")

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