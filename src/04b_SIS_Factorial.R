##########################################################################################################################

# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
# Created: January 2026; Last edited: September 2026
# This script produces S-I-S modeled disease outbreak simulations using a factorial sampling approach. This approach
# is discussed in the methods of the corresponding Derkx et al. (2026) article, but not used in the main analyses. 

# IMPORTANT NOTES:
# 1. This script relies on a .sh script with identical name for execution. 
# 2. Very important: you HAVE to check the SLURM script and adapt it. It has instructions written in it. 
# 3. Note the script is nearly identical to the 03a script, apart from its sampling approach 
# 3. Als note that this script is not performed for all locations, as it was just used as a methods check
# 4. The beta values (transmission rates) are pre-defined in this script, as they are not varied throughout the 
#    manuscript. If you want to change this, do so manually.
# 5. You have to change the path of where the data is stored in the .sh file.  

##########################################################################################################################

### SET-UP R ENVIRONMENT

# Set-up R environment
rm(list = ls())

suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(furrr)
  library(data.table)
})

# Set local root directory
LOCAL_ROOT_DIR <- "/scicore/home/chitnis/derkx0000/GraphComparison"
setwd(LOCAL_ROOT_DIR)

# Set cores for parallel processing
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 2))

##########################################################################################################################

### Helper functions 

#-------------------------------------------------------------------------------------------------------------------------
# Helper function 1: Compute neighbour list
neighbor_list <- function(g) {
  
  # This function takes a graph and computes a list of the neighbours of all individuals
  #' @param g: This is an igraph graph object
  #' @return list of integers with neighbours

  # calculate the neighbors and turn into integer
  lapply(adjacent_vertices(g, V(g), mode = "all"), as.integer)
}
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# Helper function 2: S-I-S simulation
sis_discrete <- function(nbs, start, delta, beta) {
  
  # This function has the following arguments: 
  #' @param nbs: The precomputed list of neighbors from neighbor_list().
  #' @param start: The ID of the node where the outbreak begins.
  #' @param delta: The mean duration a node is infectious (in days), for a Poisson distribution.
  #' @param beta: The transmission probability per edge per day.
  #' @return list with outbreak summary

  ## Step 1: Initialize vectors
  
  # N: Number of individuals (taken from list of neighbours)
  N <- length(nbs)
  
  # S: A vector to track the state of each node (0=S, 1=E, 2=I, 3=D).
  # IMPORTANT: S here stands for state, not for susceptible!
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
  time_step <- 0L # start at day 0 
  
  # The total infection events (S -> E)
  total_infections_event_count <- 0L
  
  # Daily state counts
  infectious_per_day <- length(infectious)
  
  
  ## Step 2: initiate the simulation loop
  
  
  # The simulation continues as long as there are infectious nodes.
  # A max of 1000 days is set, but this is a number that is very unlikely to be reached
  while (length(infectious) > 0L && time_step < 1000) {
    
    ## Step 2.1 Transitions: I > S
    
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
    
    ## Step 2.2 Transitions: S > I
    
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
  
  ## Step 3: Return list

  # Return all results for comprehensive analysis
  return(list(
    infectious_per_day = infectious_per_day,
    total_infections = total_infections_event_count,
    N = N,
    simulation_duration = time_step)) # Exclude the initial day because day 0
}
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# Helper function 3: Selecting graphs
get_graph <- function(graph_type, graph_idx) {

  # This function handles looking up a graph type and id
  #' @param graph_type: which of five graph types you are looking up
  #' @param graph_idx: the index of the graph you are looking up
  #' @return graph 

  # Ensure we are dealing with a single string, not a vector
  type_val <- as.character(graph_type[1])
  
  ## Step 1. Handle Empirical
  if (type_val == "empirical") {
    
    # Assuming your lookup has a specific structure for empirical
    return(graph_lookup$empirical$empirical)
  }
  
  ## Step 2. Handle Synthetic
  if (type_val %in% names(graph_lookup$synthetic)) {
    
    # Pull the specific instance. 
    # If you are uncertain, check the hierarchical listing structure where graphs are stored
    idx_val <- as.integer(graph_idx[1])
    return(igraph::upgrade_graph(graph_lookup$synthetic[[type_val]][[idx_val]][[2]]))  
  }
  
  stop("Unknown graph type: ", type_val)
}
#-------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------
# Helper function 4: Run a single S-I-S simulation
run_single_sim <- function(p) {
  
  # This function runs a single simulation of the SEIR model
  #' @param p: the task grid 
  #' @return list with summary output

  ## Step 1: Set-up
  
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

  ## Step 2: Run simulation
  
  # Run discrete SEIR model
  sim_result <- sis_discrete(
    nbs = nbs,
    start = start_node,
    delta = p$delta,
    beta = p$beta)
  
  ## Step 3: summarize output

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
  
  # Return output
  list(summary = summary,
       N = sim_result$N)
}
#-------------------------------------------------------------------------------------------------------------------------


##########################################################################################################################


## Step 1: set-up

# Set script arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: Rscript script.R <graph_file> <master_seed> <total_tasks> <n_graphs> <n_reps>")

# Set five arguments (see SLURM script for details)
GRAPH_FILE    <- args[1]
MASTER_SEED   <- as.integer(args[2])
TOTAL_TASKS   <- as.integer(args[3])
N_GRAPHS      <- as.integer(args[4])
N_REPS        <- as.integer(args[5])

# Obtain task id from system
TASK_ID       <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))


## Step 2: Load graphs


# Synthetic
synthetic_graphs <- readRDS(GRAPH_FILE)

# Empirical
vetted_graphs <- readRDS(
  "~/GraphComparison/Net_Sens/Graphs/Seed_5000/vetted_graphs_task_710_seed_5000.rds"
)

# Combine in list
graph_lookup <- list(
  synthetic = synthetic_graphs, 
  empirical = list(empirical = igraph::upgrade_graph(vetted_graphs$anchor)
  )
)

# Remove for memory
rm(vetted_graphs, synthetic_graphs)
gc()


## Step 3: Build factorial grid


# Set beta levels
beta_levels <- seq(0.01, 0.20, by = 0.005)

# Set delta levels
delta_levels <- seq(3,9, by = 1)

# Get graph names
syn_types <- names(graph_lookup$synthetic)

# Synthetic combinations
syn_grid <- expand.grid(
  beta = beta_levels,
  delta = delta_levels,
  graph_idx = 1:N_GRAPHS,
  graph_type = syn_types,
  replicate = 1:N_REPS,
  stringsAsFactors = FALSE
)

# Set emp reps
EMP_REPS = N_REPS * N_GRAPHS

# Empirical combinations
emp_grid <- expand.grid(
  beta = beta_levels,
  delta = delta_levels,
  graph_idx = NA_integer_,
  graph_type = "empirical",
  replicate = 1:EMP_REPS,
  stringsAsFactors = FALSE
)

# Combine
base_grid <- rbind(syn_grid, emp_grid)

# Global shuffle
set.seed(MASTER_SEED + 123)
base_grid <- base_grid[sample(nrow(base_grid)), ]

# Add simulation seeds
set.seed(MASTER_SEED) 
base_grid$sim_seed <- sample.int(.Machine$integer.max, 
                                 size = nrow(base_grid)
)


## Step 4: SLURM slicing


# Set total number of simulations
TOTAL_SIMS <- nrow(base_grid) 

# Set total number of simulations per task
SIMS_PER_TASK <- ceiling(TOTAL_SIMS / TOTAL_TASKS) 

# Set start and end id of task
start_idx <- (TASK_ID * SIMS_PER_TASK) + 1 
end_idx <- min(start_idx + SIMS_PER_TASK - 1, TOTAL_SIMS) 

# Create task grid
task_grid <- base_grid[start_idx:end_idx, ] 

cat(sprintf( "Task %d running %d simulations\n", TASK_ID, nrow(task_grid) ))


## Step 5: Parallelization execution


# Set up parallelization with n cores
future::plan(multisession, workers = N_CORES)

# Print number of cores used > check if uses all available cores
cat(
  sprintf(
    "Parallelization enabled using %d cores (future multisession).\n",
    future::nbrOfWorkers()
  )
)


## Step 6: Run simulations


# Run the single simulation over the whole task grid
results <- furrr::future_map(
  split(task_grid, seq_len(nrow(task_grid))),
  run_single_sim,
  .options = furrr_options(seed = TRUE))

# Add all summary results into one data frame
summary_df <- rbindlist(lapply(results, `[[`, "summary"))


## Step 7: Save output


# Set new directory for output
out_name <- paste0("Outfiles/SIS_Factorial_", TOTAL_SIMS)
out_dir <- file.path(LOCAL_ROOT_DIR, out_name)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# Save summary output
saveRDS(summary_df, file.path(out_dir,sprintf("Summary_Task%03d_Seed%d_TotalSimulations%d.rds", TASK_ID, MASTER_SEED, TOTAL_SIMS)))

cat("Task", TASK_ID, "finished successfully\n")



### End of script ###