##########################################################################################################################

# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
# Created: January 2026; Last edited: September 2026
# This script produces S-E-I-R modeled disease outbreak simulations using a randomized sampling approach, as per the main
# methodology of the corresponding Derkx et al. (2026) manuscript. 

# IMPORTANT NOTES:
# 1. This script relies on a .sh script with identical name for execution. 
# 2. Very important: you HAVE to check the SLURM script and adapt it. It has instructions written in it. 
# 3. In 'Step 2: Load graphs', please note that there are two options: the code will either search for the ensemble file
#    and empirical graph in the four additional network locations ("Sabaneta", "Habi", "Hepang", "Romana") or it will do so
#    in the Seed folders ("Seed_100", "Seed_1000", "Seed_2000", "Seed_5000"). The latter refer to the Chad network and 
#    corresponding graphs. For the simulations in the manuscript, "Seed_5000" is used as we generated 5000 graphs for 
#    each graph type. Which of these you will actually simulate is determined in the accompanying SLURM script. 
# 4. Make sure that in those lines 406-416, the paths perfectly correspond to your own paths! Otherwise, R won't be able
#    to find the graphs. 
# 5. The beta values (transmission rates) are pre-defined in this script, as they are not varied throughout the 
#    manuscript. If you want to change this, do so manually. 


# For any issues, feel free to report a new issue on this GitHub repo or contact inez.derkx@swisstph.ch

##########################################################################################################################

### SET UP R ENVIRONMENT ###

# Empty environment
rm(list = ls())

# Set working directory
setwd("/scicore/home/chitnis/derkx0000/Derkx2026_publication")

# Load required libraries
suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(furrr)
  library(data.table)
})

# Set local root directory where the simulation folders reside.
LOCAL_ROOT_DIR <- "/scicore/home/chitnis/derkx0000/Derkx2026_publication"
setwd(LOCAL_ROOT_DIR)

# Set cores
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
# Helper function 2: S-E-I-R simulation
seir_discrete <- function(nbs, start, delta = 2L, sigma = 30L, beta) {
  
  # This function simulates an SEIR model using the following arguments: 
  #' @param nbs The precomputed list of neighbors from neighbor_list().
  #' @param start The ID of the initially infected node (patient zero)
  #' @param delta The mean duration a node is infectious (in days)
  #' @param sigma The mean duration a node is latently infected
  #' @param beta The transmission probability per edge per day.
  #' @return list with outbreak summary
  
  ## Step 1: Initialize vectors
  
  # N: Number of individuals
  N <- length(nbs)
  
  # S: A vector to track the state of each node (0=S, 1=E, 2=I, 3=D).
  # IMPORTANT: S here stands for state, not for susceptible!!
  S <- integer(N) 
  
  # Tleft: A vector of integers to act as timers for exposed and infectious nodes.
  # It stores the number of days left until they change state.
  Tleft <- integer(N)  
  
  # Set the initial state for the 'start' node.
  # 2 = I, a.k.a. infectious (exposed individuals cannot infect yet)
  S[start] <- 2L
  
  # Sample a timer from a Poisson distribution with delta as mean
  Tleft[start] <- max(1L, rpois(1, delta))
  
  # Initialize lists to track nodes in each state.
  infectious <- start
  exposed  <- integer(0)
  removed  <- integer(0)
  
  # Initiate the time step (day) of the simulation
  time_step <- 0L # start at day 0 
  
  # The total infection events (S -> E)
  total_infections_event_count <- 0L
  
  # Daily state counts
  infectious_per_day <- length(infectious)
  exposed_per_day <- length(exposed)
  removed_per_day <- length(removed)
  
  
  ## Step 2: initiate the simulation loop

  
  # The simulation continues as long as there are exposed or infectious nodes.
  # A max of 1000 days is set as upper limit
  while ((length(infectious) > 0L || length(exposed) > 0L)  && time_step < 1000) {
    
    ## Step 2.1 Transitions: E -> I
    
    # Check if there is anyone in the exposed class
    if (length(exposed) > 0L) {
      
      # For every exposed node, decrease their timer by 1 day.
      Tleft[exposed] <- Tleft[exposed] - 1L
      
      # Identify nodes whose timer has reached 0 or less.
      becoming_I <- exposed[Tleft[exposed] <= 0L]
      
      # If there is at least one person ready to transition to the infectious state:
      if (length(becoming_I) > 0L) {
        
        # Update their state to 2 ('infectious')
        S[becoming_I] <- 2L
        
        # Assign a new timer for how long they will remain infectious.
        Tleft[becoming_I] <- rpois(length(becoming_I), delta) + 1L # Minimum of 1 day 
        
        # Add these newly sick individuals to the 'infectious' tracking list.
        infectious <- c(infectious, becoming_I)
        
        # Remove these individuals from the 'exposed' list.
        exposed <- setdiff(exposed, becoming_I)
      }
    }
    
    ## Step 2.2 Transitions: I -> R
    
    # Check if there is anyone in the infectious class
    if (length(infectious) > 0L) {
      
      # For every infectious node, decrease their timer by 1 day.
      Tleft[infectious] <- Tleft[infectious] - 1L
      
      # Identify nodes whose timer has reached 0 or less.
      becoming_R <- infectious[Tleft[infectious] <= 0L]
      
      # If there is at least one person ready to transition to the removed state:
      if (length(becoming_R) > 0L) {
        
        # Update their state to 3 ('removed')
        S[becoming_R] <- 3L
        
        # Add these  individuals to the 'removed' tracking list.
        removed <- c(removed, becoming_R)
        
        # Remove these individuals from the 'infectious' list.
        infectious <- setdiff(infectious, becoming_R)
      }
    }
    
    ## Step 2.3 Transitions: S -> E (New Infections)
    
    # New infections are based on nodes that were infectious at the start of this step
    if (length(infectious) > 0L) {
      
      # Identify neighbours of infectious individuals
      all_neighbors <- unique(unlist(nbs[infectious]))
      
      # From that list of neighbors, keep only those who are Susceptible.
      susceptible_neighbors <- all_neighbors[S[all_neighbors] == 0L]
      
      # only proceed if there's actually a susceptible person nearby.
      if (length(susceptible_neighbors) > 0L) {
        
        # For every susceptible neighbor, roll a random number (0 to 1).
        # If the number is less than 'beta' (your transmission probability), they get infected.
        # This is a vectorized operation: it rolls the dice for the whole group at once.
        newE <- susceptible_neighbors[runif(length(susceptible_neighbors)) < beta]
        
        # If the dice roll resulted in at least one new infection:
        if (length(newE) > 0L) {
          
          # Change their state in the Master State Vector 'S' to 1 (Exposed/Latent).
          S[newE] <- 1L
          
          # Assign how many days they will stay in the 'Exposed' state.
          Tleft[newE] <- rpois(length(newE), sigma) + 1L
          
          # Add the newly infected IDs to the 'exposed' tracking list.
          exposed <- c(exposed, newE)
          
          # Tracker: add to total infection number
          total_infections_event_count <- total_infections_event_count + length(newE)
        }
      }
    }
    
    # Add timestemp
    time_step <- time_step + 1L
    
    # Store the count of nodes at the end of the day.
    infectious_per_day <- c(infectious_per_day, length(infectious))
    exposed_per_day <- c(exposed_per_day, length(exposed))
    removed_per_day <- c(removed_per_day, length(removed))
  }

  ## Step 3: Return list 
  
  # Return all results, including the daily counts for E, I, R
  return(list(
    infectious_per_day = infectious_per_day,
    exposed_per_day = exposed_per_day,
    removed_per_day = removed_per_day,
    total_infections = total_infections_event_count, # Total S -> E events
    N = N,
    simulation_duration = time_step)) 
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
# Helper function 4: Run one single simulation 
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
  
  # Run one round of the discrete SEIR model
  sim_result <- seir_discrete(
    nbs = nbs,
    start = start_node,
    beta = p$beta,
    delta = p$delta,
    sigma = p$sigma)
  
  ## Step 3: summarize output
  
  # Calculate daily active cases and deaths
  active_cases_daily <- sim_result$infectious_per_day + sim_result$exposed_per_day
  deaths_daily <- sim_result$removed_per_day
  
  # Record the total duration of the simulation 
  duration <- sim_result$simulation_duration
  
  # Sum of counts across all days 
  sum_I <- sum(sim_result$infectious_per_day)
  sum_E <- sum(sim_result$exposed_per_day)
  sum_R <- sum(sim_result$removed_per_day)
  
  # Average compartment count over time 
  # (Divided by the total number of time points: duration + 1)
  avg_I_fraction <- (sum_I / (duration + 1)) / N
  avg_E_fraction <- (sum_E / (duration + 1)) / N
  avg_R_fraction <- (sum_R / (duration + 1)) / N
  
  # Average susceptible fraction over time
  avg_S_fraction <- 1 - (avg_I_fraction + avg_E_fraction + avg_R_fraction)
  
  # Create summary data
  summary <- data.frame(
    
    # Parameters used
    graph_type = p$graph_type,
    graph_idx  = p$graph_idx,
    beta       = p$beta,
    delta      = p$delta,
    sigma      = p$sigma,
    seed       = p$sim_seed,
    
    # Max case number 
    peak_active_cases = max(active_cases_daily, na.rm = TRUE), # Max of E + I
    peak_exposed_cases = max(sim_result$exposed_per_day, na.rm=TRUE), # Max of E only
    peak_infectious_cases = max(sim_result$infectious_per_day, na.rm = TRUE), # Max of I only
    peak_deaths = max(sim_result$removed_per_day, na.rm=TRUE), # Max of R only
    
    # Average case number
    avg_active_cases = mean(active_cases_daily, na.rm = TRUE), # Avg of E + I
    avg_exposed_cases = mean(sim_result$exposed_per_day, na.rm=TRUE), # Avg of E only
    avg_infectious_cases = mean(sim_result$infectious_per_day, na.rm = TRUE), # Avg of I only
    avg_deaths = mean(sim_result$removed_per_day, na.rm=TRUE), # Avg of R 
    
    # Duration
    duration_days = duration,    
    
    # Time-averaged compartment fractions
    avg_susceptible_fraction = avg_S_fraction,
    avg_exposed_fraction = avg_E_fraction,
    avg_infectious_fraction = avg_I_fraction,
    avg_removed_fraction = avg_R_fraction, 
    
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
if (length(args) < 7) {
  stop("Usage: Rscript script.R <master_seed> <total_tasks> <sims_per_task> <delta> <sigma> <n_graphs> <curr_type>")
}

# Detail seven arguments (see SLURM script)
MASTER_SEED    <- as.integer(args[1])
TOTAL_TASKS    <- as.integer(args[2])
SIMS_PER_TASK  <- as.integer(args[3])
DELTA_VAL      <- as.numeric(args[4])
SIGMA_VAL      <- as.numeric(args[5])
N_GRAPHS       <- as.numeric(args[6])
CURR_TYPE      <- as.character(args[7])


# Obtain task id from system
TASK_ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))

# Set total number of simulations
N_TOTAL_SIMS <- TOTAL_TASKS * SIMS_PER_TASK
cat("Task:", TASK_ID, "; Master seed:", MASTER_SEED, "\n")


## Step 2: Load graphs


# Read file
# IMPORTANT: see notes at top of script
if (CURR_TYPE %in% c("Sabaneta", "Habi", "Hepang", "Romana")) {
  GRAPH_FILE <- paste0("/scicore/home/chitnis/derkx0000/Derkx2026_publication/Net_Sens/Graphs/OtherLocations/", CURR_TYPE, "_MASTER_ENSEMBLE_FOR_DM.rds")
  EMP_FILE <- paste0("/scicore/home/chitnis/derkx0000/Derkx2026_publication/Net_Sens/Graphs/OtherLocations/", CURR_TYPE, "_vetted_graphs_task_710_seed_5000.rds")
} else if (CURR_TYPE %in% c("Seed_100", "Seed_1000", "Seed_2000", "Seed_5000")) {
  GRAPH_FILE <- paste0("/scicore/home/chitnis/derkx0000/Derkx2026_publication/Net_Sens/Graphs/", CURR_TYPE, "/MASTER_ENSEMBLE_FOR_DM.rds")
  EMP_FILE <- paste0("/scicore/home/chitnis/derkx0000/Derkx2026_publication/Net_Sens/Graphs/", CURR_TYPE, "/vetted_graphs_task_710_seed_5000.rds")
} else {
  stop("CURR_TYPE not correct or not given. Please check. Run aborted.\n")
}

# Read synthetic and empirical graphs
synthetic_graphs <- readRDS(GRAPH_FILE)
vetted_graphs <- readRDS(EMP_FILE)
emp_graph <- igraph::upgrade_graph(vetted_graphs$anchor)
rm(vetted_graphs)

# Add graphs to list 
graph_lookup <- list(
  synthetic = synthetic_graphs,
  empirical = list(empirical = emp_graph))


## Step 3: balanced sampling across graph types


# Set graph indices
graph_indices <- 1:N_GRAPHS

# Set beta levels according to manuscript
beta_levels <- seq(0.01, 0.20, by = 0.005)

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
  "ncrg" = list(
    source = "synthetic",
    indices = 1:100),
  "empirical" = list(
    source = "empirical",
    indices = NA_integer_))

# Get graph type names from graph registry
graph_types <- names(graph_registry)

# E.g. if total simulations is 20000, we run 4000 per graph for 5 graph types
sims_per_graph_type <- N_TOTAL_SIMS / length(graph_types)

# Make sure it's an integer
if(sims_per_graph_type != floor(sims_per_graph_type)) {
  stop("N_TOTAL_SIMS must be divisible by number of graph types")
}

# Create task grids
task_grids <- lapply(graph_types, function(gtype) {
  
  # Set a unique seed per graph type for reproducibility
  set.seed(MASTER_SEED + which(graph_types == gtype))
  
  # Generate the Random Parameters
  # We use replace = TRUE to ensure true independent random draws
  sampled <- data.frame(
    beta = sample(beta_levels, size = sims_per_graph_type, replace = TRUE),
    graph_type = gtype
  )
  
  # Assign Graph Indices randomly (Synthetic) or NA (Empirical)
  if (!is.na(graph_registry[[gtype]]$indices[1])) {
    sampled$graph_idx <- sample(graph_indices, size = sims_per_graph_type, replace = TRUE)
  } else {
    sampled$graph_idx <- NA_integer_
  }
  
  # Generate unique simulation seeds
  # This is the "ID" of the specific stochastic outbreak
  sampled$sim_seed <- sample.int(.Machine$integer.max, size = sims_per_graph_type)
  
  # Return grid
  return(sampled)
})

# Assemble task grid
task_grid <- dplyr::bind_rows(task_grids)

# Add static parameters > these don't change for the SEIR model 
task_grid$delta <- DELTA_VAL
task_grid$sigma <- SIGMA_VAL

# We shuffle all rows from task grid 
# This ensures that SLURM Task #1 isn't doing e.g. 500 DCSBM runs, 
# but a mix of all graph types and betas.
set.seed(MASTER_SEED + 123)
task_grid <- task_grid[sample(nrow(task_grid)), ]

# Slice for SLURM
start_idx <- (TASK_ID) * SIMS_PER_TASK + 1
end_idx   <- min(start_idx + SIMS_PER_TASK - 1, nrow(task_grid))

# Final grid for a specific SLURM task
task_grid <- task_grid[start_idx:end_idx, ]


## Step 4: Set up parallelization 


# Set up parallelization with n cores
future::plan(multisession, workers = N_CORES)

# Print number of cores used > check if uses all available cores
cat(
  sprintf(
    "Parallelization enabled using %d cores (future multisession).\n",
    future::nbrOfWorkers()
  )
)


## Step 5: Run simulations


# Run a single simulation over the whole task grid
results <- furrr::future_map(
  split(task_grid, seq_len(nrow(task_grid))),
  run_single_sim,
  .options = furrr_options(seed = TRUE))

# Add all summary results into one data frame
summary_df <- rbindlist(lapply(results, `[[`, "summary"))


## Step 6: Save output


# Set new directory for output
out_name <- paste0("Outfiles/", CURR_TYPE, "_Randomized_", N_TOTAL_SIMS)
out_dir <- file.path(LOCAL_ROOT_DIR, out_name)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# Save summary output
saveRDS(summary_df, file.path(out_dir,sprintf("Summary_Task%03d_Seed%d_TotalSimulations%d.rds", TASK_ID, MASTER_SEED, N_TOTAL_SIMS)))

cat("Task", TASK_ID, "finished successfully\n")



### End of script ###