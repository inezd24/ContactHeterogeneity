###########################################################################################################################


# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
# Created: 19 January 2026; Last edited: 19 January 2026
# SEIR Simulation script updated for Randomized Sampling & SLURM Integration
# Non-factorial approach


####################################################################################################
### SET-UP R ENVIRONMENT
####################################################################################################

# Empty list
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(furrr)
  library(data.table)
})

# Set local root directory where the simulation folders reside.
LOCAL_ROOT_DIR <- "/scicore/home/chitnis/derkx0000/GraphComparison"
setwd(LOCAL_ROOT_DIR)

# Set cores
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 2))


####################################################################################################
### HELPER FUNCTIONS
####################################################################################################

# 1. Neighbor list: computes all neighbors of each node
neighbor_list <- function(g) {
  
  # This function has the following arguments: 
  #' @g This is an igraph graph object
  
  # calculate the neighbors and turn into integer
  lapply(adjacent_vertices(g, V(g), mode = "all"), as.integer)
}

# 2. S, I, S function 
sis_discrete <- function(nbs, start, delta, beta) {
  
  # This function has the following arguments: 
  #' @nbs The precomputed list of neighbors from neighbor_list().
  #' @start The ID of the node where the outbreak begins.
  #' @delta The mean duration a node is infectious (in days), for a Poisson distribution.
  #' @beta The transmission probability per edge per day.
  
  ### Initialize vectors
  
  # N: Number of individuals (taken from list of neighbours)
  N <- length(nbs)
  
  # S: A vector to track the state of each node (0=S, 1=E, 2=I, 3=D).
  # IMPORTANT: S here stands for state, not for susceptible!!
  S <- integer(N)    
  
  # Tleft: A vector of integers to act as timers for exposed and infectious nodes.
  # It stores the number of days left until they change state.
  Tleft <- integer(N)  
  
  # Set the initial state for the 'start' node > 1L = infected
  S[start] <- 1L
  
  # How many days are left until the initial node changes state?
  # Use delta as mean for Poisson distribution
  Tleft[start] <- max(1L, rpois(1, delta)) 
  
  # Initialize lists to track nodes in each state.
  infectious <- start
  
  # Initiate the time step (day) of the simulation
  # time_step <- 1L
  time_step <- 0L # start at day 0 
  
  # The total infection events (S -> E)
  total_infections_event_count <- 0L
  
  # Daily state counts
  infectious_per_day <- length(infectious)
  
  
  ### The actual simulation loop
  
  
  # The simulation continues as long as there are infectious nodes.
  # A max of 1000 days is set, but this is a number that is very unlikely to be reached
  while (length(infectious) > 0L && time_step < 1000) {
    
    
    ### Recovery (I > S)
    
    
    # Check if there is anyone in the infectious class
    if (length(infectious) > 0L) {
      
      # For every infectious node, decrease their timer by 1 day.
      Tleft[infectious] <- Tleft[infectious] - 1L
      
      # Identify nodes whose timer has reached 0 or less.
      becoming_S <- infectious[Tleft[infectious] <= 0L]
      
      # If there is at least one person ready to transition to the removed state:
      if (length(becoming_S) > 0L) {
        
        # Update their state to back to 0 ("S")
        S[becoming_S] <- 0L
        
        # Remove these individuals from the 'infectious' list.
        infectious <- setdiff(infectious, becoming_S)
      }
    }
    
    
    ### New infections (S > I)
    
    # New infections are based on nodes that were Infectious at the start of this step
    if (length(infectious) > 0L) {
      
      # Identify neighbours of infectious individuals
      all_neighbors <- unique(unlist(nbs[infectious]))
      
      #From that list of neighbors, keep only those who are Susceptible.
      susceptible_neighbors <- all_neighbors[S[all_neighbors] == 0L]
      
      # only proceed if there's actually a susceptible person nearby.
      if (length(susceptible_neighbors) > 0L) {
        
        # For every susceptible neighbor, roll a random number (0 to 1).
        # If the number is less than 'beta' (your transmission probability), they get infected.
        # This is a vectorized operation: it rolls the dice for the whole group at once.
        newI <- susceptible_neighbors[runif(length(susceptible_neighbors)) < beta]
        
        # If the dice roll resulted in at least one new infection:
        if (length(newI) > 0L) {
          
          # Change their state in the Master State Vector 'S' to 1 (Infectious).
          S[newI] <- 1L
          
          # Assign how many days they will stay in the 'Infectious' state.
          Tleft[newI] <- rpois(length(newI), delta) + 1L
          
          # Add the newly infected IDs to the 'exposed' tracking list.
          infectious <- c(infectious, newI)
          
          # Tracker: add to total infection number
          total_infections_event_count <- total_infections_event_count + length(newI)
        }
      }
    }
    
    # Add time stamp
    time_step <- time_step + 1L
    
    # Store the count of nodes at the end of the day.
    infectious_per_day <- c(infectious_per_day, length(infectious))
    
  }
  
  
  ### Recovery (I > S)
  
  # We check which infected nodes have a timer of 1 or less (time to recover).
  # if (length(infectious) > 0L) {
  #  keep_infected <- Tleft[infected] > 1L
  #  recovered_nodes <- infected[!keep_infected]
  
  # If there are recovered nodes
  #  if (length(recovered_nodes) > 0L) {
  
  # Move back to susceptible state
  #   S[recovered_nodes] <- 0L
  # }
  
  # Update the list of currently infected nodes and decrement their timers.
  # infected <- infected[keep_infected]
  # Tleft[infected] <- Tleft[infected] - 1L
  # }
  
  # Store the count of infected nodes at the end of the day.
  # infected_per_day <- c(infected_per_day, length(infected))
  #  }
  
  # Return all results for comprehensive analysis
  return(list(
    infectious_per_day = infectious_per_day,
    total_infections = total_infections_event_count,
    N = N,
    simulation_duration = time_step)) # Exclude the initial day because day 0
}

# 4 selecting graph type
get_graph <- function(graph_type, graph_idx) {
  
  # Ensure we are dealing with a single string, not a vector
  type_val <- as.character(graph_type[1])
  
  # 1. Handle Empirical
  if (type_val == "empirical") {
    
    # Assuming your lookup has a specific structure for empirical
    return(graph_lookup$empirical$empirical)
  }
  
  # 2. Handle Synthetic
  if (type_val %in% names(graph_lookup$synthetic)) {
    
    # Pull the specific instance. Note the [[1]] to ensure we have a single index.
    idx_val <- as.integer(graph_idx[1])
    return(igraph::upgrade_graph(graph_lookup$synthetic[[type_val]][[idx_val]][[2]]))  
  }
  
  stop("Unknown graph type: ", type_val)
}


####################################################################################################
### EXECUTION BLOCK
####################################################################################################

# Set script arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript script.R <graph_file> <master_seed> <total_tasks> <sims_per_task> <n_graphs>")
}

# Detail six arguments
GRAPH_FILE     <- args[1]
MASTER_SEED    <- as.integer(args[2])
TOTAL_TASKS    <- as.integer(args[3])
SIMS_PER_TASK  <- as.integer(args[4])
N_GRAPHS       <- as.numeric(args[5])

# Obtain task id from system
TASK_ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))

# Set total number of simulations
N_TOTAL_SIMS <- TOTAL_TASKS * SIMS_PER_TASK
cat("Task:", TASK_ID, "; Master seed:", MASTER_SEED, "\n")


####################################################################################################
### LOAD GRAPHS
####################################################################################################

# Read graphs from graph file
synthetic_graphs <- readRDS(GRAPH_FILE)

# Load empirical graph
vetted_graphs <- readRDS(sprintf("~/GraphComparison/Net_Sens/Graphs/Seed_5000/vetted_graphs_task_710_seed_5000.rds"))
emp_graph <- igraph::upgrade_graph(vetted_graphs$anchor)
rm(vetted_graphs)

# Add to graphs
graph_lookup <- list(
  synthetic = synthetic_graphs,
  empirical = list(empirical = emp_graph))



################################################################################
### BALANCED SAMPLING ACROSS GRAPH TYPES
################################################################################


# Set graph indices
graph_indices <- 1:N_GRAPHS

# Set beta levels
beta_levels <- seq(0.01, 0.20, by = 0.005)

# Set delta levels
delta_levels <- seq(3,9, by = 1)

# For empirical and synthetic graphs
graph_registry <- list(
  "dcsbm" = list(
    source = "synthetic",
    indices = 1:100),
  "random" = list(
    source = "synthetic",
    indices = 1:100),
  "sbm" = list(
    source = "synthetic",
    indices = 1:100),
  "spatial" = list(
    source = "synthetic",
    indices = 1:100),
  "newclust_graph" = list(
    source = "synthetic",
    indices = 1:100),
  "empirical" = list(
    source = "empirical",
    indices = NA_integer_))

# Get graph types from graph registry
graph_types <- names(graph_registry)

# E.g. if total simulations is 20000, we run 5000 per graph type
sims_per_graph_type <- N_TOTAL_SIMS / length(graph_types)

# Make sure it's an integer
if(sims_per_graph_type != floor(sims_per_graph_type)) {
  stop("N_TOTAL_SIMS must be divisible by number of graph types")
}

task_grids <- lapply(graph_types, function(gtype) {
  
  # Set a unique seed per graph type for reproducibility
  set.seed(MASTER_SEED + which(graph_types == gtype))
  
  # A. Generate the Random Parameters
  # We use replace = TRUE to ensure true independent random draws
  sampled <- data.frame(
    beta = sample(beta_levels, size = sims_per_graph_type, replace = TRUE),
    delta = sample(delta_levels, size = sims_per_graph_type, replace = TRUE),
    graph_type = gtype
  )
  
  # B. Assign Graph Indices randomly (Synthetic) or NA (Empirical)
  if (!is.na(graph_registry[[gtype]]$indices[1])) {
    sampled$graph_idx <- sample(graph_indices, size = sims_per_graph_type, replace = TRUE)
  } else {
    sampled$graph_idx <- NA_integer_
  }
  
  # C. Generate unique simulation seeds
  # This is the "ID" of the specific stochastic outbreak
  sampled$sim_seed <- sample.int(.Machine$integer.max, size = sims_per_graph_type)
  
  return(sampled)
})

# --- 2. Final Assembly ---
task_grid <- dplyr::bind_rows(task_grids)

# --- 3. The Global Shuffle ---
# Still important for SLURM load balancing
set.seed(MASTER_SEED + 123)
task_grid <- task_grid[sample(nrow(task_grid)), ]

# --- 4. Slicing for SLURM ---
start_idx <- (TASK_ID) * SIMS_PER_TASK + 1
end_idx   <- min(start_idx + SIMS_PER_TASK - 1, nrow(task_grid))
task_grid <- task_grid[start_idx:end_idx, ]


####################################################################################################
### PARALLEL PLAN
####################################################################################################

# Set up parallelization with n cores
future::plan(multisession, workers = N_CORES)

# Print number of cores used > check if uses all available cores
cat(
  sprintf(
    "Parallelization enabled using %d cores (future multisession).\n",
    future::nbrOfWorkers()
  )
)

####################################################################################################
### RUN SIMULATIONS
####################################################################################################


# Run a single simulation
run_single_sim <- function(p) {
  
  # Get the right graph from the graph file
  g <- get_graph(
    graph_type = p$graph_type,
    graph_idx  = p$graph_idx)
  
  # Checks that the graph is a valid igraph object
  if (!is_igraph(g)) {
    stop("Error in run_simulation: 'graph' object is not a valid igraph object.")
  }
  
  # Check that the graph has actual vertices
  if (vcount(g) < 1) {
    stop(sprintf("Error in run_simulation: Graph '%s' has zero vertices.", 
                 deparse(substitute(g))))
  }
  
  # Identify neighbors using helper function 1.   
  nbs <- neighbor_list(g)
  
  # Save total network size
  N <- length(nbs)
  
  # Use seed from grid to set starting node
  set.seed(p$sim_seed)
  start_node <- sample.int(N, 1)
  
  # Run discrete SEIR model
  sim_result <- sis_discrete(
    nbs = nbs,
    start = start_node,
    delta = p$delta,
    beta = p$beta)
  
  # Calculate daily active cases and deaths
  active_cases_daily <- sim_result$infectious_per_day
  
  # Record the total duration of the simulation  (time_step - 1L)
  duration <- sim_result$simulation_duration
  
  # Sum of counts and average count
  sum_I <- sum(sim_result$infectious_per_day)
  avg_I_fraction <- (sum_I / (duration + 1)) / N
  
  # Average susceptible fraction over time
  avg_S_fraction <- 1 - avg_I_fraction
  
  # Create summary data
  summary <- data.frame(
    
    # Parameters used
    graph_type = p$graph_type,
    graph_idx  = p$graph_idx,
    beta       = p$beta,
    delta      = p$delta,
    seed       = p$sim_seed,
    
    # Max case number 
    peak_active_cases = max(active_cases_daily, na.rm = TRUE), # Max of E + I
    
    # Average case number
    avg_active_cases = mean(active_cases_daily, na.rm = TRUE), # Avg of E + I
    
    # Duration
    duration_days = duration,    
    
    # Average time-averaged compartment fractions
    avg_susceptible_fraction = avg_S_fraction,
    avg_infectious_fraction = avg_I_fraction,
    
    # Total S-> E events (Total cases)
    total_infections = sim_result$total_infections)
  
  # Output
  list(summary = summary,
       I = sim_result$infectious_per_day,
       N = sim_result$N)
}

# Run the single simulation over the whole task grid
results <- furrr::future_map(
  split(task_grid, seq_len(nrow(task_grid))),
  run_single_sim,
  .options = furrr_options(seed = TRUE))

# Add all summary results into one data frame
summary_df <- rbindlist(lapply(results, `[[`, "summary"))

# Turn counts (I, E, R) into time series
timeseries_df <- rbindlist(lapply(seq_along(results), function(i) {
  data.table(
    sim_id = i,
    day = seq_along(results[[i]]$I) - 1,
    I = results[[i]]$I,
    graph_type = summary_df$graph_type[i],
    graph_idx  = summary_df$graph_idx[i],
    beta       = summary_df$beta[i],
    delta      = summary_df$delta[i],
    seed       = summary_df$seed[i],
    N          = results[[i]]$N)
}))


####################################################################################################
### OUTPUT
####################################################################################################


# Set new directory for output
out_name <- paste0("Outfiles/SIS_Randomized_", N_TOTAL_SIMS)
out_dir <- file.path(LOCAL_ROOT_DIR, out_name)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# Save summary output
saveRDS(summary_df, file.path(out_dir,sprintf("Summary_Task%03d_Seed%d_TotalSimulations%d.rds", TASK_ID, MASTER_SEED, N_TOTAL_SIMS)))
saveRDS(timeseries_df, file.path(out_dir,sprintf("TimeSeries_Task%03d_Seed%d_TotalSimulations%d.rds", TASK_ID, MASTER_SEED, N_TOTAL_SIMS)))

cat("Task", TASK_ID, "finished successfully\n")



### End of script ###