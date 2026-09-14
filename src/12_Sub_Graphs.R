###########################################################################################################################

# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
# Created: 09 February 2025; Last edited: 09 February 2025
# Script compares different network construction algorithms

# Legend of script:

### IS FOR NEW SECTIONS IN CAPITAL ###
### Is for headings (e.g., a new function)
# is for 'small' commands (e.g., rename or merge)

##########################################################################################################################

### SET UP R ENVIRONMENT ###

# Empty list
rm(list = ls())

# Set library path
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
library(readr)

# Set local root directory 
LOCAL_ROOT_DIR <- "/scicore/home/chitnis/derkx0000/GraphComparison"
setwd(LOCAL_ROOT_DIR)

# Output folders
plot_path <- file.path(LOCAL_ROOT_DIR, "Plots")

##########################################################################################################################

### Helper functions


# 1. Generate empirical_graph
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

# 2. General Function to calculate metrics for graphs coming out of scaling algorithm
calculate_network_metrics <- function(g, label) {
  if (is.null(g)) return(NULL)
  
  # Create desired output
  data.frame(
    graph_type = label,
    nodes = vcount(g),
    edges = ecount(g),
    density = edge_density(g),
    avg_degree = mean(degree(g)),
    med_degree = median(degree(g)),
    max_degree = max(degree(g)),
    transitivity = transitivity(g, type = "global"),
    avg_path = mean_distance(g, directed = FALSE),
    stringsAsFactors = FALSE)
}

# 3. Haversine distance
haversine_dist <- function(lon1, lat1, lon2, lat2) {
  R <- 6371000 # Earth radius in meters
  phi1 <- lat1 * pi / 180
  phi2 <- lat2 * pi / 180
  delta_phi <- (lat2 - lat1) * pi / 180
  delta_lambda <- (lon2 - lon1) * pi / 180
  
  a <- sin(delta_phi/2)^2 + cos(phi1) * cos(phi2) * sin(delta_lambda/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  return(R * c)
}

# 4. Distance statistic 1: Kolmogorov-Smirnov (KS) Distance
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

# 5. Inducing sub-graphs
get_subgraph_metrics <- function(full_graph, 
                                 n_samples = 5,
                                 sample_perc = 0.25, 
                                 method = 'random',
                                 households, 
                                 devices,
                                 seed = NULL) {
  
  # Function goal: automate obtaining sub-graphs from empirical network
  # Function parameters:
  #' @param full_graph is the actual empirical graph 
  #' @param n_samples How many graphs do we want? Default is set at 5 graphs
  #' @param sample_size How big each graph? Default is set at 50 as approx 25% of real graph (235)
  #' @param method Sampling method used. Default is 'random', but can also be 'neighbourhood' or 'spatial'. 
  #' @param households Data frame with all households of the dogs
  #' @param devices The devices that were actually employed
  #' @param seed Setting a seed for reproducibility. Default is set to no seed. 
  
  # Set a seed if desired
  if (!is.null(seed)) set.seed(seed)
  
  # Get vertices of full graph
  emp_v <- V(full_graph)
  
  # Get degrees of full graph
  degrees_emp <- degree(full_graph)
  
  # Set sample size
  sample_size <- round(sample_perc * length(emp_v), 0)
  
  # For each required sub-graph, sample using one of two methods:
  map_dfr(1:n_samples, function(i) {
    
    ## Random sampling
    if (method == "random") {
      
      # Standard random sampling
      sub_v <- sample(emp_v, size = sample_size, replace = FALSE)
      
      
      ## Breadth-first search sampling
    } else if (method == "neighbourhood") {
      
      # Pick a random starting node and take its nearest neighbors (Breadth-First Search)
      start_node <- sample(emp_v, 1)
      
      # Get nodes ordered by distance from the start_node
      bfs_res <- igraph::bfs(full_graph, root = start_node, unreachable = TRUE, order = TRUE) # unreachable = TRUE bc we always want 50 nodes. 
      
      # Take nodes up to the set sample size 
      sub_v <- bfs_res$order[1:sample_size]
      
      
      ## Spatial sampling
    } else if (method == 'spatial') {
      
      # Get overview of used devices 
      gps <- households %>%
        rename(x = CoordN, y = CoordE) %>%
        inner_join(devices %>% select(gcs = V1), by = 'gcs')
      
      # Calculate Bounding Box
      lon_min <- min(gps$x)
      lon_max <- max(gps$x)
      lat_min <- min(gps$y)
      lat_max <- max(gps$y)
      
      # Calculate Dimensions in Meters
      avg_lat <- (lat_min + lat_max) / 2
      width_m <- haversine_dist(lon_min, avg_lat, lon_max, avg_lat)
      height_m <- haversine_dist(lon_min, lat_min, lon_min, lat_max)
      
      # Calculate area
      total_area_km2 <- (width_m * height_m) / 1e6 # 1.32 km2, seems about right
      
      # Print the metrics
      cat("--- Bounding Box Metrics ---\n")
      cat(sprintf("Width:  %.2f meters\n", width_m))
      cat(sprintf("Height: %.2f meters\n", height_m))
      cat(sprintf("Area:   %.2f square km\n\n", total_area_km2))
      
      # Sample 25% of area (so 0.5 of width and 0.5 of height)
      sample_width_deg <- (lon_max - lon_min) * 0.5
      sample_height_deg <- (lat_max - lat_min) * 0.5
      
      # Randomly sample the lower part of the bounding box
      rand_lon <- runif(1, lon_min, lon_max - sample_width_deg)
      rand_lat <- runif(1, lat_min, lat_max - sample_height_deg)
      
      # Get starting and ending coords
      lon_start = rand_lon
      lon_end   = rand_lon + sample_width_deg
      lat_start = rand_lat
      lat_end   = rand_lat + sample_height_deg
      
      # Select nodes that are within the coords
      nodes <- gps %>%
        filter(x >= lon_start & x <= lon_end ) %>%
        filter(y >= lat_start & y <= lat_end ) %>%
        select(gcs)
      
      # Select only those nodes
      sub_v <- emp_v[name %in% nodes$gcs]
      
    }
    
    # Make sub graph
    sub_g <- induced_subgraph(full_graph, sub_v)
    
    # Calculate the desired network metrics
    summary <- calculate_network_metrics(sub_g, paste0("empirical_sub_", i))
    
    # Compute ks distance
    degrees_sub <- degree(sub_g)
    ks <- ks_distance(degrees_sub, degrees_emp)
    
    # Add to summary
    summary <- summary %>%
      mutate(ks_dist = ks, 
             avg_degree_perc = avg_degree / nodes, 
             max_degree_perc = max_degree / nodes)
  })
}


### Execution

# Get contact data
ServerData <- readRDS("~/BaseData/Chad/ServerData.rds")
if (!exists("ServerData")) {
  stop("Error: The 'ServerData' object was not loaded from ServerData.rds.")
}

# Get GPS data
Household_Locations <- read_delim("~/BaseData/Chad/Household_Locations_6.csv", 
                                    delim = ";", escape_double = FALSE, trim_ws = TRUE)

# And device data
deployed_devices <- read.table("~/BaseData/Chad/Deployed_Devices/deployed_devices_6.txt", quote="\"", comment.char="")

# Transform data to network
empirical_graph <- empirical_net(ServerData$dog_contact_zone_1[,c(2,3)])
empirical_graph <- empirical_graph$graph

# Apply random sampling to get sub graphs
sub_sample_random <- get_subgraph_metrics(empirical_graph, 
                                          n_samples = 1000, 
                                          sample_perc = 0.25, 
                                          method = 'random',
                                          households = NULL,
                                          devices = NULL,
                                          seed = 42)
mean(sub_sample_random$density)# 0.063 compared to 0.063 full graph
mean(sub_sample_random$transitivity) # 0.391 compared to 0.404 full graph
mean(sub_sample_random$avg_path) # 3.142 compared to 2.870 full graph

# Also test neighbourhood sampling
sub_sample_neigh <- get_subgraph_metrics(empirical_graph, 
                                         n_samples = 1000, 
                                         sample_perc = 0.25, 
                                         method = 'neighbourhood',
                                         households = NULL,
                                         devices = NULL,
                                         seed = 42)
mean(sub_sample_neigh$density)# 0.254 compared to 0.063 full graph
mean(sub_sample_neigh$transitivity) # 0.518 compared to 0.404 full graph
mean(sub_sample_neigh$avg_path) # 1.931 compared to 2.870 full graph

# Also test neighbourhood sampling
sub_sample_spat <- get_subgraph_metrics(empirical_graph, 
                                         n_samples = 1000, 
                                         sample_perc = 0.25, 
                                         method = 'spatial',
                                         households = Household_Locations,
                                         devices = deployed_devices,
                                         seed = 42)
mean(sub_sample_spat$density)# 0.132 compared to 0.063 full graph
mean(sub_sample_spat$transitivity) # 0.472 compared to 0.404 full graph
mean(sub_sample_spat$avg_path) # 2.570 compared to 2.870 full graph

# Combine
sub_sample_metrics <- sub_sample_random %>%
  mutate(method = 'random') %>%
  bind_rows(sub_sample_neigh %>% mutate(method = 'neighbourhood')) %>%
  bind_rows(sub_sample_spat %>% mutate(method = 'spatial'))
write.csv(sub_sample_metrics, file.path(LOCAL_ROOT_DIR, "Outfiles/sub_sample_metrics.csv"))

# Colour palette
palette_paper <- c("#3A405A", "#FF8465", "#99B2DD", "#FCD2A2", "#A26F39")

# Visualisation of nodes and edges
nodes_edges <- ggplot(sub_sample_metrics, aes(x = nodes, y = edges, group = method)) +
  geom_point(aes(color = method)) +
  scale_color_manual(values = palette_paper) +
  theme_minimal() +
  labs(title = "Node and Edge comparison", y = "Edges", x = "Nodes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot_path, "nodes_edges_subgraphs.png"), nodes_edges, width = 10, height = 10)

# Visualization of transitivity and density
trans_dens <- ggplot(sub_sample_metrics, aes(x = density, y = transitivity, group = method)) +
  geom_point(aes(color = method)) +
  geom_point(aes(x=0.063247863, y=0.40411693), color = 'blue') +
  scale_color_manual(values = palette_paper) +
  theme_minimal() +
  labs(title = "Transitivity and density comparison", y = "Global Transitivity", x = "Edge density") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot_path, "trans_dens_subgraphs.png"), trans_dens, width = 10, height = 10)

# Visualisation of max and average degree
avg_max_degree <- ggplot(sub_sample_metrics, aes(x = avg_degree, y = max_degree, group = method)) +
  geom_point(aes(color = method)) +
  scale_color_manual(values = palette_paper) +
  theme_minimal() +
  labs(title = "Degree metric comparison averaged over all runs", y = "(avg) Maximum degree", x = "(avg) Average degree") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot_path, "avg_max_degree_subgraphs.png"), avg_max_degree, width = 10, height = 10)

# Visualisation of max and average degree fractions
xv <- mean(degree(empirical_graph))/vcount(empirical_graph)
yv <- max(degree(empirical_graph))/vcount(empirical_graph)
avg_max_frac_degree <- ggplot(sub_sample_metrics, aes(x = avg_degree_perc, y = max_degree_perc, group = method)) +
  geom_point(aes(color = method)) +
  scale_color_manual(values = palette_paper) +
  geom_point(aes(x=xv, y=yv), color = 'blue') +
  theme_minimal() +
  labs(title = "Degree metric comparison averaged over all runs", y = "(avg) Maximum degree perc", x = "(avg) Average degree perc") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot_path, "avg_max_degree_frac_subgraphs.png"), avg_max_degree, width = 10, height = 10)


### End of script ###