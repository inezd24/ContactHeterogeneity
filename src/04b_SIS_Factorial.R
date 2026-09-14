###########################################################################################################################


# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
# Created: 14 January 2026; Last edited: 19 January 2026
# SIS Simulation script updated for Randomized Sampling & SLURM Integration
# Non-factorial approach


####################################################################################################
### SET-UP R ENVIRONMENT
####################################################################################################

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

# 2. S, E, I, R function without birth or death rate (for now)
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

# 3. selecting graph type
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


################################################################################
### EXECUTION BLOCK
################################################################################


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: Rscript script.R <graph_file> <master_seed> <total_tasks> <n_graphs> <n_reps>")

GRAPH_FILE    <- args[1]
MASTER_SEED   <- as.integer(args[2])
TOTAL_TASKS   <- as.integer(args[3])
N_GRAPHS      <- as.integer(args[4])
N_REPS        <- as.integer(args[5])

TASK_ID       <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))


################################################################################
### LOAD GRAPHS
################################################################################


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
rm(vetted_graphs, synthetic_graphs); gc()
gc()


################################################################################
### BUILD FACTORIAL GRID
################################################################################


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


################################################################################
### SLURM SLICING
################################################################################


TOTAL_SIMS <- nrow(base_grid) 
SIMS_PER_TASK <- ceiling(TOTAL_SIMS / TOTAL_TASKS) 
start_idx <- (TASK_ID * SIMS_PER_TASK) + 1 
end_idx <- min(start_idx + SIMS_PER_TASK - 1, TOTAL_SIMS) 
task_grid <- base_grid[start_idx:end_idx, ] 

cat(sprintf( "Task %d running %d simulations\n", TASK_ID, nrow(task_grid) ))


################################################################################
### PARALLELIZATION
################################################################################


future::plan(multisession, workers = N_CORES)

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
  
  # Run discrete SIS model
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
  list(summary = summary)
}


################################################################################
### RUN SIMULATION
################################################################################


# Run the single simulation over the whole task grid
results <- furrr::future_map(
  split(task_grid, seq_len(nrow(task_grid))),
  run_single_sim,
  .options = furrr_options(seed = TRUE))

# Add all summary results into one data frame
summary_df <- rbindlist(lapply(results, `[[`, "summary"))


################################################################################
### OUTPUT
################################################################################

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