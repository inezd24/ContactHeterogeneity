###########################################################################################################################

# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
# Created: 16 January 2026; Last edited: 19 January 2026
# Based on script 7 (plot_simulations) but adapted for the randomized sampling 

# Legend of script:

### IS FOR NEW SECTIONS IN CAPITAL ###
### Is for headings (e.g., a new function)
# is for 'small' commands (e.g., rename or merge)

# Start with empty list
rm(list = ls(all = TRUE))

##########################################################################################################################

##### SCRIPT ARGUMENTS #####

# Get arguments from the command line (passed by Slurm)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript 05_Graph_Plotting.sh <model> <type> <master_seed> <total_sims>", call. = FALSE)
}

# Set arguments
curr_model <- args[1] # e.g., "SEIR"
curr_type <- args[2] #e.g., 'Randomized" or "Hepang"
master_seed  <- args[3] # e.g., 42
total_sims <- args[4] # e.g., 250000, TOTAL_TASKS * SIMS_PER_TASK
message(paste("Processing Model:", curr_model, " for ", curr_type, " with Seed:", master_seed, ", yielding total simulations: ", total_sims))


##########################################################################################################################

### SET UP R ENVIRONMENT ###

# Load necessary packages
library(ggplot2)
library(ggpubr)
library(data.table)
library(stringr)
library(purrr)
library(dplyr)
library(tidyr)

# Set directories: working directory, plots, out files
setwd("/scicore/home/chitnis/derkx0000/GraphComparison/")
files_directory <- "/scicore/home/chitnis/derkx0000/GraphComparison/Outfiles/"

##########################################################################################################################


##### LOAD DATA #####

# Function to open summary files
open_rds <- function(file_dir, 
                     file_pattern, 
                     model, 
                     type, 
                     total_simulations, 
                     seed){
  
  # File path for locations
  files_path_new <- file.path(file_dir, paste0(type, "_", model, "_Randomized_", total_simulations))
  
  # File path normal
  files_path_old <- file.path(file_dir, paste0(model, "_", type, "_", total_simulations))
  
  # Determine which path actually exists on the cluster
  if (dir.exists(files_path_new)) {
    files_path <- files_path_new
  } else if (dir.exists(files_path_old)) {
    files_path <- files_path_old
  } else {
    warning(sprintf("Neither path found:\n  [New]: %s\n  [Old]: %s", files_path_new, files_path_old))
    return(NULL)
  }
  
  # Get list of all files
  files <- list.files(path = files_path, 
                      pattern = file_pattern, 
                      full.names = TRUE)
  
  if (length(files) == 0) return(NULL)
  
  # Open all files 
  files_list <- lapply(files, function(f) {
    dt <- as.data.table(readRDS(f))
    if (nrow(dt) == 0) return(NULL)
    return(dt) 
  })
  
  # Bind all data frames together
  summary_df <- rbindlist(files_list, fill = TRUE)
  
  cat(sprintf("Successfully merged %d files from %s\n", length(files), files_path))
  return(summary_df)
}

#### SUMMARY DATA #####

# Open files
summary_data <- open_rds(file_dir = files_directory, 
                         file_pattern = "Summary", 
                         model = curr_model, 
                         type = curr_type,
                         total_simulations = total_sims, 
                         seed = master_seed)

# Check
if(is.null(summary_data) || nrow(summary_data) == 0) stop("Summary: no files found for this combination!")

# Save output
date_suffix <- format(Sys.Date(), "_%d%B%y")
summary_name = paste0(curr_type, "_Summary_", curr_model,"_", total_sims,"_mseed_", master_seed, date_suffix, ".csv")
write.csv(summary_data, file.path(files_directory, summary_name), row.names = FALSE)



# Disease cases (from summary)
cases <- summary_data %>%
  # Filtering out failed outbreaks (peak > 1) to focus on established epidemics
  #filter(peak_active_cases > 1) %>%
  pivot_longer(cols = c(peak_active_cases, 
                        avg_active_cases, 
                        total_infections), 
               names_to = "measure", 
               values_to = 'value') %>%
  group_by(graph_type, beta, delta, measure) %>%
  summarise(min = min(value), 
            q_25 = quantile(value, 0.25, names = FALSE),
            median = median(value),
            q_75 = quantile(value, 0.75, names = FALSE),
            max = max(value), .groups = 'drop') %>%
  mutate(delta_fac = as.factor(as.character(delta)),
         measure_fac = as.factor(measure))

# Dynamically select columns that exist in the dataframe (handles Exposed/Recovered presence)
comp_cols <- names(summary_data)[grep("fraction", names(summary_data))]

# Set order
target_order <- c(
  "avg_removed_fraction",
  "avg_infectious_fraction", 
  "avg_exposed_fraction",
  "avg_susceptible_fraction")

# Compartment fractions (from summary)
fractions <- summary_data %>%
  select(graph_type, beta, delta, all_of(comp_cols)) %>%
  pivot_longer(cols = all_of(comp_cols), names_to = "compartment", values_to = 'fraction') %>%
  group_by(graph_type, beta, delta, compartment) %>%
  summarise(mean = mean(fraction), .groups = 'drop') %>%
  mutate(delta_fac = as.factor(as.character(delta))) %>%
  mutate(compartment = factor(compartment, levels = intersect(target_order, unique(compartment))))

# remove mother file for memory
rm(summary_data)

# Colour palette
scenario_colors <-  c("#3A405A", "#FF8465", "#99B2DD", "#FCD2A2", "#A26F39")

# Set labels
new_cases <- c("peak_active_cases" = "Peak Active Cases", "avg_active_cases" = "Avg Prev.", "total_infections" = "Total Infections")

### Cases

# Plot the peak and average number of cases as well as total infections by beta
p1 <- ggplot(cases, aes(x = beta, color = graph_type)) +
  geom_ribbon(aes(ymin=q_25, ymax=q_75, fill = graph_type), alpha = 0.3, color = NA) +
  geom_line(aes(y = median), linewidth =1) +
  labs(title = paste0(curr_model, ": Simulation Summary"),
       subtitle = paste0("Total simulations: ", total_sims),
       x = "Transmission Rate (β)", 
       y = "Cases") + 
  theme_bw() +
  scale_color_manual(values = scenario_colors) +
  scale_fill_manual(values = scenario_colors) +
  theme(legend.position = "bottom",
        text = element_text(size = 15))

# Facet by measure (and delta if SIS)
if(length(unique(cases$delta)) > 1){
  p1 <- p1 + facet_grid(measure_fac ~ delta_fac, scales = "free", labeller = labeller(measure_fac = new_cases))
} else {
  p1 <- p1 + facet_wrap(~measure_fac, scales = "free", labeller = labeller(measure_fac = new_cases))
}

### Compartments

# For labelling the fractions
all_labels <- c(
  "avg_susceptible_fraction" = "Susceptible fraction",
  "avg_exposed_fraction"     = "Exposed fraction",
  "avg_infectious_fraction"  = "Infectious fraction",
  "avg_removed_fraction"     = "Removed fraction") # For SIS

# Create a filtered label vector for this specific model
current_labels <- all_labels[levels(fractions$compartment)] 

# Plot fraction of compartments in area by beta
p2 <- ggplot(fractions, aes(x = beta, y = mean, fill = compartment)) +
  geom_area(alpha = 0.8) +
  labs(title = paste0(curr_model, ": Compartment Fractions across Beta"),
       subtitle = paste0("Total simulations: ", total_sims),
       y = "Fraction") +
  theme_bw() +
  scale_fill_manual(values = scenario_colors, labels = current_labels) +
  theme(legend.position = "bottom",
        text = element_text(size = 15))

# Facet by delta ONLY if there is more than one delta value (Static in SEIR)
if(length(unique(fractions$delta)) > 1){
  p2 <- p2 + facet_grid(graph_type ~ delta_fac)
} else {
  p2 <- p2 + facet_wrap(~graph_type)
}


### Save output

ggsave(filename = paste0("Plots/", curr_type, "_", curr_model,"_", total_sims,"_mseed_", master_seed, "_SummaryPlots1.png"), plot = p1, width = 14, height = 8)
ggsave(filename = paste0("Plots/", curr_type, "_", curr_model,"_", total_sims,"_mseed_", master_seed, "_Fractions_beta1.png"), plot = p2, width = 14, height = 8)

cat("Summary Plots saved to working directory /Plots\n")

rm(p1)
rm(p2)
rm(cases)
rm(fractions)
gc()



#### TIMESERIES DATA #####


# Open files
time_data <- open_rds(file_dir = files_directory, 
                      file_pattern = "TimeSeries", 
                      model = curr_model, 
                      type = curr_type,
                      total_simulations = total_sims, 
                      seed = master_seed)

# Check
if(is.null(time_data) || nrow(time_data) == 0) stop("Time series: no files found for this combination!")

# Get total number of individuals
N <- unique(time_data$N) # or set N <- 1000 for example

# Create fractions
if (curr_model == 'SEIR') {
  ts_frac <- time_data %>%
    mutate(
      avg_susceptible_fraction = (N - (E + I + R))/N, 
      avg_exposed_fraction = E / N,
      avg_removed_fraction = R / N,
      avg_infectious_fraction = I / N) %>%
    select(day, starts_with("frac"), everything())  # optional: reorder
} else {
  ts_frac <- time_data %>%
    mutate(
      avg_susceptible_fraction = (N - I)/N, 
      avg_infectious_fraction = I / N ) %>%
    select(day, starts_with("frac"), everything())  # optional: reorder
}


# Time series preparation
ts_cols <- names(ts_frac)[grep("frac", names(ts_frac))]

# Compartment fractions (over time, from time series)
time_fractions <- ts_frac %>%
  select(graph_type, beta, delta, day, all_of(ts_cols)) %>%
  pivot_longer(cols = all_of(ts_cols), names_to = "compartment", values_to = 'fraction') %>%
  group_by(graph_type, beta, delta, day, compartment) %>%
  summarise(q_25 = quantile(fraction, 0.25, names = FALSE),
            median = median(fraction),
            q_75 = quantile(fraction, 0.75, names = FALSE),
            mean = mean(fraction),
            N = n(), 
            .groups = 'drop') %>%
  filter(beta == 0.01 | beta == 0.2) %>%
  mutate(delta_fac = as.factor(as.character(delta)),
         beta_fac = as.factor(as.character(beta)))

# Plot compartment fraction by time series 
p3 <- ggplot(time_fractions, aes(x = day, y = mean, fill = compartment)) +
  geom_area(alpha = 0.8) +
  labs(title = paste0(curr_model, ": Compartment Fractions over time"), 
       subtitle = paste0("Total simulations: ", total_sims),
       y = "Fraction") +
  theme_bw() +
  scale_fill_manual(values = scenario_colors, labels = current_labels) +
  theme(legend.position = "bottom",
        text = element_text(size = 15))

# Facet by delta ONLY if there is more than one delta value (Static in SEIR)
if(length(unique(time_fractions$delta)) > 1){
  p3 <- p3 + facet_grid(graph_type ~ delta_fac)
} else {
  p3 <- p3 + facet_grid(beta_fac ~graph_type)
}

#Save
ggsave(filename = paste0("Plots/", curr_type, "_", curr_model, "_", total_sims, "_Fractions_time.png"), plot = p3, width = 10, height = 8)

# Save output
time_name = paste0(curr_type, "_TimeSeries_",curr_model,"_", total_sims, "_8May26.csv")
write.csv(time_data, file.path(files_directory, time_name), row.names = FALSE)

cat("Time Series Plots saved to working directory /Plots\n")
cat("Done!\n")




