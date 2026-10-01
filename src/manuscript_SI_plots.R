###########################################################################################################################

# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
# Created: March 2026; Last edited: September 2026
# This script creates the exact plots and tables used in the manuscript SI

# Important notes:
# 1. The images are manually saved by date, e.g. "Figure1_27_05_26.png". If you want, you can automate this for some form
#    of version control. I prefer doing it manually.
# 2. The ordering and naming of the graphs is very important here. It appears it may not always be consistent across graphs,
#    but the result is that all images use the same graph type naming and alphabetical order. Please do not change this. 
# 3. After saving each plot, the script removes large input files and the plots themselves for efficient memory use. Comment
#    those script lines if you want to inspect the plots within the environment. 

##########################################################################################################################


### SET UP R ENVIRONMENT ###

# Empty list
rm(list = ls())

# Load required libraries
library(dplyr) # data wrangling
library(tidyr) # data wrangling
library(purrr) # data opening
library(ggplot2) # for plotting
library(readr) # for opening the csv files
library(ggh4x) # for the faceting
library(scales) # for comma in 30,000
library(egg) # for faceting
library(ggtext) # for font labels
library(igraph) # graph object
library(ggraph) # for igraph plotting
library(cowplot) # for putting graphs together
library(data.table)
library(ggpubr)

# Set local root directory 
LOCAL_ROOT_DIR <- "/scicore/home/chitnis/derkx0000/Derkx2026_publicationn"
setwd(LOCAL_ROOT_DIR)

# Output folders
out_path <- file.path(LOCAL_ROOT_DIR, "Manuscript")
ifelse(!dir.exists(out_path),
       dir.create(out_path), FALSE)

##########################################################################################################################

### Helper functions

#-------------------------------------------------------------------------------------------------------------------------
# Helper function 1: find local maxima in outbreak data
find_outbreak_threshold_sensitive <- function(values, adjust_bw) {

  # This function finds the local maxima in our outbreak data. 
  # See the 'optimal_bandwidth.R' script for a thorough analysis
  #' @param 'values'
  #' @param 'adjust_bw': value to adjust bandwidth for sensitivity
  
  # 1. Basic cleaning
  values <- na.omit(values)
  if(length(unique(values)) < 5) return(0.05)
  
  # 2. Set sensitivity using a small 'adjust' value
  dens <- density(values, adjust = adjust_bw, from = 0, to = 1)
  
  # 3. Find peaks (local maxima)
  #    A point is a peak if it's higher than its two neighbors
  is_peak <- diff(sign(diff(dens$y))) == -2
  peaks <- which(is_peak) + 1
  
  # 4. Filter peaks: Ignore tiny "wiggles" by requiring peaks to be a certain 
  #    percentage of the maximum density height
  significant_peaks <- peaks[dens$y[peaks] > (max(dens$y) * 0.02)]
  
  if(length(significant_peaks) >= 2) {
    # The valley is the lowest point between the FIRST peak and the LAST peak
    first_p <- significant_peaks[1]
    last_p  <- significant_peaks[length(significant_peaks)]
    
    # Slice the density to only the area between peaks
    between_y <- dens$y[first_p:last_p]
    valley_idx <- which.min(between_y) + first_p - 1
    threshold  <- dens$x[valley_idx]
    if (!exists("threshold")) cat("Threshold not calculate. Revert to quantile-based threshold.\n")
    if (exists("threshold")) cat("Threshold calculated.\n")
    
  } else {
    # If bimodal logic fails (common at very low or very high beta), we use a quantile-based fallback. 
    # For bimodal epidemic data, the 10th percentile is often a safe 'floor' for major outbreaks 
    # NOTE: this is only suitable if the distribution is wide.
    if(max(values) > 0.1) {
      threshold <- 0.05 # A standard "5% prevalence" definition for epidemic establishment
    } else {
      threshold <- 1.0 # Everything is minor if the max is tiny
    }
  }
  
  return(threshold)
}
#-------------------------------------------------------------------------------------------------------------------------

##########################################################################################################################

### ALL MANUSCRIPT SI IMAGES

# Colour palette 1
scenario_colors <- c(
  "ERM"       = "#FAA76B",
  "NCRG"      = "#5C8A6F",
  "SBM"       = "#FCD2A2",
  "SENCA"     = "#AB4B54",
  "DCSBM"     = "#99B2DD",
  "Empirical" = "#495AA5")

# Colour palette 2
seed_colors <- c(
  "100"      = "#CAD2C5",
  "1000"      = "#84A98C",
  "2000"    = "#52796F",
  "5000"    = "#354F52")


#-------------------------------------------------------------------------------------------------------------------------
## SUPPLEMENTARY FIGURE 1  ----
#-------------------------------------------------------------------------------------------------------------------------

# Topic: six graph parameters by graph type and seed

# Open all summaries
files <- list.files(
  path = "Net_Sens", 
  pattern = "^summary_check\\.csv$", 
  full.names = TRUE, 
  recursive = TRUE
)
summary_checks <- do.call(rbind, lapply(files, read.csv))

# Leave out spatial optimization parameters, select only 1 replicate since identical
summary_check_no_opt <- summary_checks %>% filter(!grepl('Opt', Params), 
                                                  Replicate == 1)
summary_check_no_opt$Graph_type <- ifelse(summary_check_no_opt$Graph == 'dcsbm_graph', 
                                          'DCSBM',
                                          ifelse(summary_check_no_opt$Graph == 'empirical_graph', 
                                                 'Empirical',
                                                 ifelse(summary_check_no_opt$Graph == 'random_graph', 
                                                        "ERM",
                                                        ifelse(summary_check_no_opt$Graph == 'sbm_graph', 
                                                               'SBM', 
                                                               ifelse(summary_check_no_opt$Graph == 'newclust_graph', "NCRG", "SENCA"
                                                                      )))))

summary_check_no_opt$Graph_type <- factor(summary_check_no_opt$Graph_type, 
                                          levels = c("DCSBM", 'Empirical', "ERM","NCRG", "SBM", "SENCA"))

summary_check_no_opt$Parameters <- ifelse(summary_check_no_opt$Params == 'Density', 
                                          'Density', 
                                          ifelse(summary_check_no_opt$Params == 'Avg_clust', 
                                                 'Global Clustering',
                                                 ifelse(summary_check_no_opt$Params ==  'Avg_Betweenness', 
                                                        'Mean Betweenness',
                                                        ifelse(summary_check_no_opt$Params == "Avg_Degree", 
                                                               'Mean Degree',
                                                               ifelse(summary_check_no_opt$Params ==  'Med_Degree', 
                                                                      'Median Degree', 
                                                                      "Maximum Degree"
                                                               )
                                                        )
                                                 )
                                          )
)
summary_check_no_opt$Parameters <- factor(summary_check_no_opt$Parameters, 
                                          levels = c('Density', 
                                                     'Global Clustering', 
                                                     "Mean Betweenness", 
                                                     "Mean Degree", 
                                                     'Median Degree', 
                                                     "Maximum Degree")
)


# Plot 1: network topology overview per seed number
summary_plot_seed <- ggplot(summary_check_no_opt, aes(x = Graph_type, 
                                                      y = Values, 
                                                      fill = as.character(SeedTotal))) +
  geom_boxplot(position=position_dodge(1)) +
  theme_bw() +
  scale_fill_manual(values = seed_colors) +
  facet_wrap(~Parameters, scales = 'free', ncol = 3, axis.labels = "margins") +
  labs(y = "Network metric value",
       x = "Graph type",
       fill = "Seed Total") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    text = element_text(size = 20, family = "sans"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    strip.background = element_rect(fill = '#EEEBD3', color = "black"),
    strip.text = element_text(size = 20),
    panel.grid.major = element_line(color = "grey90"), # Visible major grid
    panel.grid.minor = element_blank(),               # Remove distracting minor lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.spacing = unit(1.5, "lines"))
summary_plot_seed
plotpath <- file.path(out_path, "Supp_Figure_1_27_05_26.png")
ggsave(plotpath, summary_plot_seed, width = 14, height = 14)           


#-------------------------------------------------------------------------------------------------------------------------
## SUPPLEMENTARY FIGURE 2  ----
#-------------------------------------------------------------------------------------------------------------------------

# Topic: running mean of betweenness and degree

# Select only degree and betweenness
plot_data <- summary_check_no_opt %>%
  filter(Parameters == "Mean Degree" | Parameters == 'Mean Betweenness') %>%
  filter(Graph_type != 'Empirical') %>%
  select(-(Replicate)) %>%
  group_by(Graph_type, SeedTotal, Parameters) %>%
  arrange(Seed) %>%
  mutate(n = row_number()) %>%
  
  # Calculate cumulative statistics
  mutate(
    running_mean = cumsum(Values) / n,
    running_var = (cumsum(Values^2) / n) - (running_mean^2),
    running_sd = sqrt(pmax(running_var, 0)),
    running_se = running_sd / sqrt(n),
    ci_lower = running_mean - (1.96 * running_se),
    ci_upper = running_mean + (1.96 * running_se)) %>%
  ungroup() %>%
  mutate(SeedSize = factor(as.character(SeedTotal), levels = c("5000", "2000", "1000", "100")))

# Plot running mean
running_mean_plot <- ggplot(plot_data, aes(x = n, y = running_mean, group = SeedSize)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = SeedSize), alpha = 0.2) +
  geom_line(aes(color = SeedSize), linewidth = 1.5) +
  facet_grid(Parameters~Graph_type, scales = "free") +
  labs(
    x = "Number of Seeds Sampled (n)",
    y = "Cumulative average of measure",
    fill = 'Total Seeds',
    color = "Total Seeds") +
  theme_bw() +
  scale_fill_manual(values = seed_colors) +
  scale_color_manual(values = seed_colors) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    text = element_text(size = 20, family = "sans"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    strip.background = element_rect(fill = '#EEEBD3', color = "black"),
    strip.text = element_text(size = 20),
    panel.grid.major = element_line(color = "grey90"), # Visible major grid
    panel.grid.minor = element_blank(),               # Remove distracting minor lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.spacing = unit(1.5, "lines"))
running_mean_plot
plotpath <- file.path(out_path, "Supp_Figure_2_27_05_26.png")
ggsave(plotpath, running_mean_plot, width = 14, height = 10)      


#-------------------------------------------------------------------------------------------------------------------------
## SUPPLEMENTARY FIGURE 3 ----
#-------------------------------------------------------------------------------------------------------------------------

# Topic: Kolmogorov-Smirnov distances by graph type and seed

# Import
ks_distances_100 <- read_csv("Net_Sens/Graphs/Seed_100/ks_distances_all_sims.csv") %>%
  mutate(Seed_Total = '100')
ks_distances_1000 <- read_csv("Net_Sens/Graphs/Seed_1000/ks_distances_all_sims.csv") %>%
  mutate(Seed_Total = '1000')
ks_distances_2000 <- read_csv("Net_Sens/Graphs/Seed_2000/ks_distances_all_sims.csv") %>%
  mutate(Seed_Total = '2000')
ks_distances_5000 <- read_csv("Net_Sens/Graphs/Seed_5000/ks_distances_all_sims.csv") %>%
  mutate(Seed_Total = '5000')

# Combine 
all_ks_distances <- ks_distances_100 %>%
  bind_rows(ks_distances_1000, ks_distances_2000, ks_distances_5000) %>%
  mutate(Graph = ifelse(graph_type == 'dcsbm', 'DCSBM',
                        ifelse(graph_type == 'sbm', 'SBM',
                               ifelse(graph_type == 'spatial', 'SENCA',
                                      ifelse(graph_type == 'newclust_graph', 'NCRG', "ERM"))))) %>%
  mutate(Graph_Type = factor(Graph, levels = c("DCSBM", 'Empirical', "ERM","NCRG", "SBM", "SENCA"))) %>%
  mutate(SeedSize = factor(as.character(Seed_Total), levels = c("5000", "2000", "1000", "100")))


# Plot
ks_summary <- ggplot(all_ks_distances, aes(x = as.factor(Graph_Type), 
                                           y = ks_dist, 
                                           fill = SeedSize)) +
  geom_boxplot() +
  labs(x = "Graph type",
       y = "K-S distance",
       fill = "Seed Total") +
  theme_bw() +
  scale_fill_manual(values = seed_colors) +
  scale_color_manual(values = seed_colors) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    text = element_text(size = 20, family = "sans"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    strip.background = element_rect(fill = '#EEEBD3', color = "black"),
    strip.text = element_text(size = 20),
    panel.grid.major = element_line(color = "grey90"), # Visible major grid
    panel.grid.minor = element_blank(),               # Remove distracting minor lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.spacing = unit(1.5, "lines"))
ks_summary
plotpath <- file.path(out_path, "Supp_Figure_3_27_05_26.png")
ggsave(plotpath, ks_summary, width = 10, height = 10)    


#-------------------------------------------------------------------------------------------------------------------------
## SUPPLEMENTARY FIGURE 4 ----
#-------------------------------------------------------------------------------------------------------------------------

# Load data
MASTER_ENSEMBLE_FOR_DM <- readRDS("~/Derkx2026_publication/Net_Sens/Graphs/Seed_5000/MASTER_ENSEMBLE_FOR_DM.rds")
all_seeds_between_distributions_5000 <- read_csv("Net_Sens/Graphs/Seed_5000/all_seeds_between_distributions.csv")
empirical <- readRDS("~/Derkx2026_publication/Net_Sens/Graphs/Seed_5000/vetted_graphs_task_3981_seed_5000.rds")
emp_between <- table(round(igraph::betweenness(empirical$anchor), 0))
emp_between_df <- data.frame(task_id = 1,
                             between  = as.numeric(names(emp_between)),
                             frequency = as.numeric(emp_between),
                             type = "empirical")

# Combine
between_summarized <- all_seeds_between_distributions_5000 %>%
  bind_rows(emp_between_df) %>%
  rename(Type = type) %>%
  group_by(Type, between) %>%
  summarise(q_25 = quantile(frequency, 0.25, names = FALSE),
            median = median(frequency),
            q_75 = quantile(frequency, 0.75, names = FALSE))
between_summarized$Graph <- ifelse(between_summarized$Type == 'dcsbm', 'DCSBM',
                                   ifelse(between_summarized$Type == 'sbm', 'SBM',
                                          ifelse(between_summarized$Type == 'empirical', 'Empirical',
                                                 ifelse(between_summarized$Type == 'spatial', 'SENCA', 
                                                        ifelse(between_summarized$Type == 'newclust_graph', 'NCRG', "ERM")))))
between_summarized$Graph  <- factor(between_summarized$Graph, levels = c('ERM', 'SBM','SENCA','NCRG', 'DCSBM','Empirical'))

# Plot
between_plot <- ggplot(between_summarized, aes(x = between, y = median, group = Graph, color = Graph, fill = Graph)) +
  
  # Ribbon first, then Line so the line sits on top
  geom_ribbon(aes(ymin = q_25, ymax = q_75), alpha = 0.2, color = NA) + 
  geom_line(linewidth = 1) +
  
  # Labels
  labs(
    y = "Frequency",
    x = "Betweenness",
    color = "Network Type",
    fill = "Network Type") +
  
  # Scale
  # Color Palette
  scale_color_manual(
    values = scenario_colors,
    breaks = c("DCSBM", "Empirical", "ERM",'NCRG', "SBM", "SENCA")
  ) +
  scale_fill_manual(
    values = scenario_colors,
    breaks = c("DCSBM", "Empirical", "ERM",'NCRG', "SBM", "SENCA")
  ) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  scale_x_log10() +
  
  # Theme
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    text = element_text(size = 20),
    # Full Grid Look as requested
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    # Polish legend
    legend.title = element_text(),
    legend.background = element_blank(),
    plot.subtitle = element_text(size = 12, color = "grey30", margin = margin(b = 10)))
between_plot
plotpath <- file.path(out_path, "Supp_Figure_4_27_05_26.png")
ggsave(plotpath, between_plot, width = 14, height = 10)


#-------------------------------------------------------------------------------------------------------------------------
## SUPPLEMENTARY FIGURE 5  ----
#-------------------------------------------------------------------------------------------------------------------------

# Topic: Remainder of SIS delta values (4 to 8)

Summary_SIS_ran <- read_csv("Outfiles/Summary_SIS_6000000_mseed_100_01June26.csv")
SIS_cases <- Summary_SIS_ran %>%
  
  # Filtering out failed outbreaks (peak > 1) to focus on established epidemics
  filter(total_infections > 0) %>%
  filter(delta >= 4 & delta <= 8) %>%
  pivot_longer(cols = c(peak_active_cases, avg_active_cases, total_infections), names_to = "measure", values_to = 'value') %>%
  group_by(graph_type, beta, delta, measure) %>%
  summarise(min = min(value), 
            q_25 = quantile(value, 0.25, names = FALSE),
            median = median(value),
            mean = mean(value),
            q_75 = quantile(value, 0.75, names = FALSE),
            max = max(value), .groups = 'drop') %>%
  mutate(delta_fac = as.factor(as.character(delta)))

# Adapt graph type names
SIS_cases$Graph <- ifelse(SIS_cases$graph_type == 'dcsbm', 'DCSBM',
                          ifelse(SIS_cases$graph_type == 'sbm', 'SBM',
                                 ifelse(SIS_cases$graph_type == 'empirical', 'Empirical',
                                        ifelse(SIS_cases$graph_type == 'spatial', 'SENCA',
                                               ifelse(SIS_cases$graph_type == 'newclust_graph', 'NCRG', "ERM")))))
SIS_cases$Graph <- factor(SIS_cases$Graph, 
                          levels = c("ERM",
                                     'SBM',
                                     'SENCA', 
                                     'NCRG',
                                     'DCSBM', 
                                     'Empirical'))
SIS_cases$delta_fac <- paste0("\u03B4 = ", SIS_cases$delta)

# Plot all cases
SIS_cases_plot <- ggplot(SIS_cases, aes(x = beta, color = Graph, fill = Graph)) +
  
  # Use a slightly thinner ribbon and a distinct mean line
  geom_ribbon(aes(ymin = q_25, ymax = q_75), alpha = 0.2, color = NA) +
  geom_line(aes(y = median), linewidth = 1.1) +
  
  # Labels
  labs(
    x = "Transmission Rate (\u03B2)", 
    y = NULL, # Remove the generic 'Cases' label; the facet strips usually explain the Y
    color = "Graph Type", 
    fill = "Graph Type") + 
  
  # Faceting
  facet_grid(
    measure ~ delta_fac, 
    #nest_line = element_line(colour = "black"),
    scales = "free_y", 
    labeller = labeller(
      measure = c(
        avg_active_cases = "Average Prevalence",
        peak_active_cases = "Peak Prevalence",
        total_infections = "Total Infections"))) +
  
  # Color Palette
  scale_color_manual(
    values = scenario_colors,
    breaks = c("DCSBM", "Empirical", "ERM", "SBM", "SENCA")
  ) +
  scale_fill_manual(
    values = scenario_colors,
    breaks = c("DCSBM", "Empirical", "ERM", "SBM", "SENCA")
  ) +
  
  # Refined Theme
  theme_bw() + # Base theme
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    text = element_text(size = 20, family = "sans"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    strip.background = element_rect(fill = '#EEEBD3', color = "black"),
    strip.text = element_text(size = 20),
    panel.grid.major = element_line(color = "grey90"), # Visible major grid
    panel.grid.minor = element_blank(),               # Remove distracting minor lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.spacing = unit(1.5, "lines")) +
  
  # Use commas for thousands (30,000 instead of 30000)
  scale_y_continuous(label=comma)
SIS_cases_plot
plotpath <- file.path(out_path, "Supp_Figure_5_27_05_26.png")
ggsave(plotpath, SIS_cases_plot, width = 16, height = 16) 


#-------------------------------------------------------------------------------------------------------------------------
## SUPPLEMENTARY FIGURE 6 ----
#-------------------------------------------------------------------------------------------------------------------------

# Topic: factorial sampling of SIS and SEIR graphs

# Import
Summary_SIS_ran <- read_csv("Outfiles/Factorial_Summary_SIS_163800000_mseed_100_05June26.csv")
Summary_SEIR_ran <- read_csv("Outfiles/Factorial_Summary_SEIR_23400000_mseed_100_05June26.csv")
empirical <- readRDS("~/Derkx2026_publication/Net_Sens/Graphs/Seed_5000/vetted_graphs_task_3981_seed_5000.rds")
nodes <- length(V(empirical$anchor))

# SEIR summary
SEIR_cases <- Summary_SEIR_ran %>%
  
  # Filtering out failed outbreaks (peak > 1) to focus on established epidemics
  filter(total_infections > 0) %>%
  mutate(peak_active_cases_perc = peak_active_cases/nodes, 
         avg_active_cases_perc = avg_active_cases/nodes,
         total_infections = total_infections/nodes) %>%
  pivot_longer(cols = c(peak_active_cases_perc, 
                        avg_active_cases_perc, 
                        total_infections), 
               names_to = "measure", 
               values_to = 'value'
  ) %>%
  group_by(graph_type, beta, delta, measure
  ) %>%
  summarise(min = min(value), 
            q_25 = quantile(value, 0.25, names = FALSE),
            median = median(value),
            mean = mean(value),
            q_75 = quantile(value, 0.75, names = FALSE),
            max = max(value), .groups = 'drop')
SEIR_cases$delta_fac <- "bold('(A)') ~ italic('SEIR model, \u03B4 = 2 days')"

# SIS summary
SIS_cases <- Summary_SIS_ran %>%
  
  # Filtering out failed outbreaks (peak > 1) to focus on established epidemics
  filter(total_infections > 0) %>%
  filter(delta == 3 | delta == 9
  )%>%
  mutate(peak_active_cases_perc = peak_active_cases/nodes, 
         avg_active_cases_perc = avg_active_cases/nodes
  ) %>%
  pivot_longer(cols = c(peak_active_cases_perc, 
                        avg_active_cases_perc, 
                        total_infections), 
               names_to = "measure", 
               values_to = 'value'
  ) %>%
  group_by(graph_type, beta, delta, measure
  ) %>%
  summarise(min = min(value), 
            q_25 = quantile(value, 0.25, names = FALSE),
            median = median(value),
            mean = mean(value),
            q_75 = quantile(value, 0.75, names = FALSE),
            max = max(value), .groups = 'drop'
  ) %>%
  mutate(delta_fac = as.factor(as.character(delta)))
SIS_cases$delta_fac <- ifelse(SIS_cases$delta_fac == '3', 
                              "bold('(B)') ~ italic('SIS model, \u03B4 = 3 days')", 
                              "bold('(C)') ~ italic('SIS model, \u03B4 = 9 days')")
# Join
model_cases <- rbind(SEIR_cases, SIS_cases)

# Adapt graph type names
model_cases$Graph <- ifelse(model_cases$graph_type == 'dcsbm', 
                            'DCSBM',
                            ifelse(model_cases$graph_type == 'sbm', 
                                   'SBM',
                                   ifelse(model_cases$graph_type == 'empirical', 
                                          'Empirical',
                                          ifelse(model_cases$graph_type == 'spatial', 
                                                 'SENCA', 
                                                 ifelse(model_cases$graph_type == 'newclust_graph', "NCRG",
                                                 "ERM")))))
model_cases$Graph <- factor(model_cases$Graph, 
                            levels = c("ERM",
                                       'SBM',
                                       'SENCA', 
                                       'NCRG',
                                       'DCSBM', 
                                       'Empirical'))
model_cases$Empty <- " "
model_cases$Empty2 <- " "

# Colour palette
scenario_colors <- c(
  "ERM"       = "#FAA76B",
  "NCRG"      = "#5C8A6F",
  "SBM"       = "#FCD2A2",
  "SENCA"     = "#AB4B54",
  "DCSBM"     = "#99B2DD",
  "Empirical" = "#495AA5")

# Plot all cases
cases_plot <- ggplot(model_cases, aes(x = beta, 
                                      color = Graph, 
                                      fill = Graph)) +
  
  # Use a slightly thinner ribbon and a distinct mean line
  geom_ribbon(aes(ymin = q_25, ymax = q_75), 
              alpha = 0.2, 
              color = NA) +
  geom_line(aes(y = median), 
            linewidth = 1.1) +
  
  # Labels
  labs(
    x = "Transmission Rate (\u03B2)", 
    y = "Fraction or total number",
    color = "Graph Type", 
    fill = "Graph Type") + 
  
  # Faceting
  facet_nested_wrap(
    ~ Empty + delta_fac + Empty2 + measure, 
    scales = "free", 
    ncol = 3,
    
    # Apply styling level-by-level
    strip = strip_nested(
      background_x = list(
        element_blank(),
        element_rect(fill = '#EEEBD3'),             # Level 1: delta_fac 
        element_blank(),                            # Level 2: empty spacer
        element_rect(fill = '#FFFFFF')),            # Level 3: measure
      text_x = list(
        element_blank(),
        element_text(hjust = 0),                    # Level 1: Letter left, Delta italic
        element_blank(),                            # Level 2: hidden
        element_text(face = "plain", hjust = 0.5)), # Level 3: centered measure
      by_layer_x = TRUE),
    labeller = labeller(
      delta_fac = label_parsed, 
      measure = c(
        avg_active_cases_perc = "Average Prevalence",
        peak_active_cases_perc = "Peak Prevalence",
        total_infections = "Total Infections")
    )) +
  
  # Color Palette
  scale_color_manual(
    values = scenario_colors,
    breaks = c("DCSBM", "Empirical", "ERM", "SBM", "SENCA")
  ) +
  scale_fill_manual(
    values = scenario_colors,
    breaks = c("DCSBM", "Empirical", "ERM", "SBM", "SENCA")
  ) +
  
  # Refined Theme
  theme_bw() + # Base theme
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 15),
    text = element_text(size = 15),
    # Clean up the strips (headers)
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    # This makes the gap between Level 2 and 3 so small the line can't draw
    panel.spacing.y = unit(0, "lines"),
    # Ensure the 'Full Grid' look
    panel.grid.major = element_line(color = "grey90"), # Visible major grid
    panel.grid.minor = element_blank(),  # Remove distracting minor lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    # Spacing between facets
    panel.spacing = unit(1.5, "lines")) +
  
  # Use commas for thousands (30,000 instead of 30000)
  scale_y_continuous(label=comma)
cases_plot
plotpath <- file.path(out_path, "Supp_Figure_6_27_05_26.png")
ggsave(plotpath, cases_plot, width = 12, height = 16)   


#-------------------------------------------------------------------------------------------------------------------------
## SUPPLEMENTARY FIGURE 7 ----
#-------------------------------------------------------------------------------------------------------------------------

# Topic: factorial sampling outbreak types

# Import SEIR data
size_and_duration <-read_csv("Outfiles/Factorial_Summary_SEIR_23400000_mseed_100_05June26.csv") %>%
  #filter(total_infections > 0) %>%
  select(graph_type, graph_idx, beta, delta, sigma, seed, peak_deaths, duration_days, total_infections) %>%
  mutate(Graph = ifelse(graph_type == 'dcsbm', 'DCSBM',
                        ifelse(graph_type == 'empirical', 'Empirical',
                               ifelse(graph_type == 'random', 'ERM',
                                      ifelse(graph_type == 'sbm', 'SBM', 
                                             ifelse(graph_type == 'newclust_graph', 'NCRG', 'SENCA'))))),
         final_size = peak_deaths/235) %>%
  filter(beta == 0.01 | beta == 0.1 | beta == 0.2) %>%
  rename_with(~c("duration"), c(duration_days))
size_and_duration$Graph <- factor(size_and_duration$Graph, 
                                  levels = c("DCSBM", "Empirical", "ERM",'NCRG', 'SBM', "SENCA"))

# Apply to your data table by group
setDT(size_and_duration)
size_and_duration[, dynamic_threshold := find_outbreak_threshold_sensitive(values = final_size,
                                                                           adjust_bw = 0.3), 
                  by = .(Graph, beta)]

# Classify based on the specific threshold found for that group
size_and_duration[, outbreak_type := ifelse(final_size <= dynamic_threshold, 
                                            "Minor", 
                                            "Major")]
table(size_and_duration$outbreak_type)

# Custom colour palette
graph_levels <- c("DCSBM", "Empirical", "ERM",'NCRG', 'SBM', "SENCA")
custom_pal <- c(
  "ERM.Minor" = "#FAA76B", "ERM.Major" = "#CB8350", # Orange 
  "SBM.Minor"         = "#FCD2A2", "SBM.Major"         = "#C4A27B",
  "SENCA.Minor"     = "#AB4B54", "SENCA.Major"     = "#78343A", 
  "DCSBM.Minor"       = "#99B2DD", "DCSBM.Major"       = "#3E5787",
  "Empirical.Minor"   = "#495AA5", "Empirical.Major"   = "#273160",
  "NCRG.Minor"        = "#5C8A6F", "NCRG.Major"        = '#3B604B') 

# Plot
comparison_plot <- ggplot(size_and_duration, aes(x = duration, y = final_size)) +
  geom_point(
    aes(color = interaction(Graph, outbreak_type, sep = ".")),
    alpha = 0.8) +
  scale_color_manual(values = custom_pal,
                     # This tells ggplot: "Only show these 5 entries in the legend"
                     breaks = paste0(graph_levels, ".Major"),
                     # This renames them from "DCSBM.Major" back to just "DCSBM"
                     labels = graph_levels) +
  facet_grid(beta ~ Graph) +
  labs(
    color = "Network Type (Major outbreak = Darker Shading)",
    x = "Duration (Days)",
    y = "Final Outbreak Size") +
  theme_bw() + 
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15),
    text = element_text(size = 20, family = "sans"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    strip.background = element_rect(fill = '#EEEBD3', color = "black"),
    strip.text = element_text(size = 20),
    panel.grid.major = element_line(color = "grey90"), # Visible major grid
    panel.grid.minor = element_blank(),               # Remove distracting minor lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.spacing = unit(1.5, "lines"))
comparison_plot
plotpath <- file.path(out_path, "Supp_Figure_7_27_05_26.png")
ggsave(plotpath, comparison_plot, width = 18, height = 12)  


#-------------------------------------------------------------------------------------------------------------------------
## SUPPLEMENTARY FIGURE 8  ----
#-------------------------------------------------------------------------------------------------------------------------

# Topic: Factorial sampling SIS remaining deltas (4 to 8)

Summary_SIS_ran <- read_csv("Outfiles/Factorial_Summary_SIS_163800000_mseed_100_05June26.csv")
setDT(Summary_SIS_ran)
SIS_cases <- Summary_SIS_ran %>%
  
  # Filtering out failed outbreaks (peak > 1) to focus on established epidemics
  filter(total_infections > 0) %>%
  filter(delta >= 4 & delta <= 8) %>%
  pivot_longer(cols = c(peak_active_cases, avg_active_cases, total_infections), names_to = "measure", values_to = 'value') %>%
  group_by(graph_type, beta, delta, measure) %>%
  summarise(min = min(value), 
            q_25 = quantile(value, 0.25, names = FALSE),
            median = median(value),
            mean = mean(value),
            q_75 = quantile(value, 0.75, names = FALSE),
            max = max(value), .groups = 'drop') %>%
  mutate(delta_fac = as.factor(as.character(delta)))

# Adapt graph type names
SIS_cases$Graph <- ifelse(SIS_cases$graph_type == 'dcsbm', 'DCSBM',
                          ifelse(SIS_cases$graph_type == 'sbm', 'SBM',
                                 ifelse(SIS_cases$graph_type == 'empirical', 'Empirical',
                                        ifelse(SIS_cases$graph_type == 'spatial', 'SENCA', 
                                               ifelse(SIS_cases$graph_type == "newclust_graph", "NCRG", "ERM")))))
SIS_cases$Graph <- factor(SIS_cases$Graph, 
                          levels = c("ERM",
                                     'SBM',
                                     'SENCA',
                                     'NCRG',
                                     'DCSBM', 
                                     'Empirical'))
SIS_cases$delta_fac <- paste0("\u03B4 = ", SIS_cases$delta)

# Colour palette
scenario_colors <- c(
  "ERM"       = "#FAA76B",
  "NCRG"      = "#5C8A6F",
  "SBM"       = "#FCD2A2",
  "SENCA"     = "#AB4B54",
  "DCSBM"     = "#99B2DD",
  "Empirical" = "#495AA5")

# Plot all cases
SIS_cases_plot <- ggplot(SIS_cases, aes(x = beta, color = Graph, fill = Graph)) +
  
  # Use a slightly thinner ribbon and a distinct mean line
  geom_ribbon(aes(ymin = q_25, ymax = q_75), alpha = 0.2, color = NA) +
  geom_line(aes(y = median), linewidth = 1.1) +
  
  # Labels
  labs(
    x = "Transmission Rate (\u03B2)", 
    y = NULL, # Remove the generic 'Cases' label; the facet strips usually explain the Y
    color = "Graph Type", 
    fill = "Graph Type") + 
  
  # Faceting
  facet_grid(
    measure ~ delta_fac, 
    #nest_line = element_line(colour = "black"),
    scales = "free_y", 
    labeller = labeller(
      measure = c(
        avg_active_cases = "Average Prevalence",
        peak_active_cases = "Peak Prevalence",
        total_infections = "Total Infections"))) +
  
  # Color Palette
  scale_color_manual(
    values = scenario_colors,
    breaks = c("DCSBM", "Empirical", "ERM","NCRG", "SBM", "SENCA")
  ) +
  scale_fill_manual(
    values = scenario_colors,
    breaks = c("DCSBM", "Empirical", "ERM","NCRG", "SBM", "SENCA")
  ) +
  
  # Refined Theme
  theme_bw() + # Base theme
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    text = element_text(size = 20, family = "sans"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    strip.background = element_rect(fill = '#EEEBD3', color = "black"),
    strip.text = element_text(size = 20),
    panel.grid.major = element_line(color = "grey90"), # Visible major grid
    panel.grid.minor = element_blank(),               # Remove distracting minor lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.spacing = unit(1.5, "lines")) +
  
  # Use commas for thousands (30,000 instead of 30000)
  scale_y_continuous(label=comma)
SIS_cases_plot
plotpath <- file.path(out_path, "Supp_Figure_8_27_05_26.png")
ggsave(plotpath, SIS_cases_plot, width = 16, height = 16) 


#-------------------------------------------------------------------------------------------------------------------------
## SUPPLEMENTARY FIGURE 9  ----
#-------------------------------------------------------------------------------------------------------------------------

# Topic: outbreak type for SEIR randomized

# Adapt df
size_and_duration <-  read_csv("size_and_duration_SEIR_3000000.csv")
size_duration_longer <- size_and_duration %>%
  select(Graph, beta, outbreak_type, final_size, duration) %>%
  pivot_longer(
    cols = duration:final_size,
    names_to = "metric",
    values_to = 'values') %>%
  mutate(empty1 = "",
         empty2 = "",
         metric_type = ifelse(metric == 'final_size', 'Final Outbreak Size', 'Outbreak Duration'))

# Plot
outbreak_type_plot <- ggplot(size_duration_longer, aes(x = as.factor(beta), y = values, fill = outbreak_type)) +
  geom_boxplot(
    alpha = 0.8, 
    width = 0.7, 
    outlier.alpha = 0.2, 
    outlier.size = 0.5,
    linewidth = 0.6) +
  facet_grid(
    metric_type ~ Graph, 
    scales = "free_y") +
  scale_fill_manual(values = c("Minor" = "#99B2DD", "Major" = "#4F6DA0")) +
  labs(x = "Transmission Rate (\u03B2)",
       y = "Values",
       fill = "Outbreak Classification") +
  theme_bw() + # Base theme
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15),
    text = element_text(size = 20, family = "sans"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    strip.background = element_rect(fill = '#EEEBD3', color = "black"),
    strip.text = element_text(size = 20),
    panel.grid.major = element_line(color = "grey90"), # Visible major grid
    panel.grid.minor = element_blank(),               # Remove distracting minor lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.spacing = unit(1.5, "lines")) 
outbreak_type_plot
plotpath <- file.path(out_path, "Supp_Figure_9_27_05_26.png")
ggsave(plotpath, outbreak_type_plot, width = 16, height = 8)


#-------------------------------------------------------------------------------------------------------------------------
## SUPPLEMENTARY FIGURE 11 ----
#-------------------------------------------------------------------------------------------------------------------------

# Topic: Graph examples for all graph types and four locations

# Open ensemble of graphs for seed = 5000
locations <- c("Habi", "Hepang", "Sabaneta", "Romana")

# 1. Alphabetical titles (The order your final plots will display)
titles <- c("DCSBM", "Empirical", "ERM", "NCRG", "SBM", "SENCA") 

# 2. Named palettes (mapped precisely to the model identifiers)
model_names_alpha <- c("dcsbm", "empirical", "random", "ncrg", "sbm", "spatial") 

node_pal <- c("#99B2DD","#495AA5", "#FAA76B","#5C8A6F" ,"#FCD2A2", "#AB4B54") %>% 
  purrr::set_names(model_names_alpha)

stroke_pal <- c("#4F6DA0","#344078","#CB8350", "#3B604B", "#C4A27B", "#78343A") %>% 
  purrr::set_names(model_names_alpha)

# Empirical graphs
empirical_graphs <- list(
  Sabaneta = readRDS("~/Derkx2026_publication/Net_Sens/Graphs/OtherLocations/Sabaneta_vetted_graphs_task_4933_seed_5000.rds")[[1]],
  Hepang = readRDS("~/Derkx2026_publication/Net_Sens/Graphs/OtherLocations/Hepang_vetted_graphs_task_4933_seed_5000.rds")[[1]],
  Habi = readRDS("~/Derkx2026_publication/Net_Sens/Graphs/OtherLocations/Habi_vetted_graphs_task_4933_seed_5000.rds")[[1]],
  Romana = readRDS("~/Derkx2026_publication/Net_Sens/Graphs/OtherLocations/Romana_vetted_graphs_task_4933_seed_5000.rds")[[1]])

# Plot the data
final_plots <- locations %>% 
  purrr::set_names(locations) %>% 
  map(function(loc) {
    
    # Load RDS files with graphs per location
    file_pattern <- paste0("^", loc, ".*_MASTER_ENSEMBLE_FOR_DM\\.rds$")
    file_path <- list.files(path = file.path(LOCAL_ROOT_DIR, "Net_Sens/Graphs/OtherLocations"), 
                            pattern = file_pattern, full.names = TRUE)
    
    if (length(file_path) == 0) {
      warning(paste("No file found for:", loc))
      return(NULL)
    }
    
    rds_data <- readRDS(file_path[1])
    
    # 3. Safely build a named list of ALL graphs
    all_graphs_raw <- list(
      dcsbm     = list(graph = rds_data[["dcsbm"]][[1]][[2]],     colour = node_pal["dcsbm"],     stroke_colour = stroke_pal["dcsbm"]),
      random    = list(graph = rds_data[["random"]][[1]][[2]],    colour = node_pal["random"],    stroke_colour = stroke_pal["random"]),
      ncrg      = list(graph = rds_data[["ncrg"]][[1]][[2]],      colour = node_pal["ncrg"],      stroke_colour = stroke_pal["ncrg"]),
      sbm       = list(graph = rds_data[["sbm"]][[1]][[2]],       colour = node_pal["sbm"],       stroke_colour = stroke_pal["sbm"]),
      spatial   = list(graph = rds_data[["spatial"]][[1]][[2]],   colour = node_pal["spatial"],   stroke_colour = stroke_pal["spatial"]),
      empirical = list(graph = empirical_graphs[[loc]],           colour = node_pal["empirical"], stroke_colour = stroke_pal["empirical"])
    )
    
    # 4. Enforce strict sequential order matching the alphabetical 'titles' array
    loc_graphs <- all_graphs_raw[c("dcsbm", "empirical", "random", "ncrg", "sbm", "spatial")]
    
    # Plot graphs 
    styled_plots <- Map(function(item, ttl) {
      
      network <- item$graph
      set.seed(42)
      d <- igraph::degree(network)
      node_sizes <- rescale(sqrt(d), to = c(1, 10))
      
      # Use unname() to strip structure and pass purely the hex string value to ggraph
      plot_fill   <- unname(item$colour)
      plot_stroke <- unname(item$stroke_colour)
      
      ggraph(network, layout = "nicely") +
        geom_edge_link(edge_color = "gray70", alpha = 0.7) +
        geom_node_point(fill = plot_fill, shape = 21, color = plot_stroke, 
                        stroke = 0.8, size = node_sizes, alpha = 0.8) +
        labs(title = ttl) + 
        theme_graph() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 18, family = 'sans'),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))
    }, loc_graphs, titles)
    
    # Combine into the 1x6 composite for this location
    plot_grid(
      plotlist = styled_plots, 
      ncol = 6, 
      label_size = 18) + 
      labs(caption = paste("Location:", loc)) +
      theme(plot.background = element_rect(fill = "white", color = NA))
  })

# Put graphs in single arrangement
graphs <- ggarrange(final_plots$Habi,
                    final_plots$Hepang,
                    final_plots$Sabaneta,
                    final_plots$Romana, 
                    ncol = 1,
                    labels = c("(A)", "(B)", "(C)", "(D)"))
graphs
plotpath <- file.path(out_path, "Supp_Figure_11_27_05_26.png")
ggsave(plotpath, graphs, width = 20, height = 16) 


#-------------------------------------------------------------------------------------------------------------------------
## SUPPLEMENTARY FIGURE 12 ----
#-------------------------------------------------------------------------------------------------------------------------

# Topic: Graph parameters for all 4 locations and models

# Set locations and colours
locations <- c("Habi", "Hepang", "Sabaneta", "Romana")
scenario_colors <- c(
  "ERM"       = "#FAA76B",
  "NCRG"      = "#5C8A6F",
  "SBM"       = "#FCD2A2",
  "SENCA"     = "#AB4B54",
  "DCSBM"     = "#99B2DD",
  "Empirical" = "#495AA5")

# Load and clean summary files
all_summary_data <- map_dfr(locations, function(loc) {
  
  # Look for the specific summary_check file for each location
  file_path <- list.files(path = file.path(LOCAL_ROOT_DIR, "Net_Sens/OtherLocations"), 
                          pattern = paste0(loc, "_summary_check\\.csv$"), full.names = TRUE)
  
  if(length(file_path) == 0) return(NULL)
  
  read_csv(file_path) %>%
    filter(!grepl('Opt', Params), Replicate == 1) %>%
    mutate(Location = loc) # Keep track of the site
})

# Apply new factor levels
all_summary_data <- all_summary_data %>%
  mutate(
    Graph_type = case_when(
      Graph == 'dcsbm_graph' ~ 'DCSBM',
      Graph == 'empirical_graph' ~ 'Empirical',
      Graph == 'random_graph' ~ "ERM",
      Graph == 'ncr_graph' ~ "NCRG",
      Graph == 'sbm_graph' ~ 'SBM',
      TRUE ~ "SENCA"
    ),
    Graph_type = factor(Graph_type, levels = c("DCSBM", "Empirical", "ERM","NCRG", 'SBM', "SENCA")),
    
    Parameters = case_when(
      Params == 'Density' ~ 'Density', 
      Params == 'Avg_clust' ~ 'Global Clustering',
      Params == 'Avg_Betweenness' ~ 'Mean Betweenness',
      Params == 'Avg_Degree' ~ 'Mean Degree',
      Params == 'Med_Degree' ~ 'Median Degree', 
      TRUE ~ "Maximum Degree"
    ),
    Parameters = factor(Parameters, levels = c('Density', 'Global Clustering', "Mean Betweenness", "Mean Degree", 'Median Degree', "Maximum Degree"))
  )

# Plot
combined_summary_plot <- ggplot(all_summary_data, aes(x = Graph_type, y = Values, fill = Graph_type)) +
  
  # Standardized Boxplots
  geom_boxplot(
    width = 0.6, 
    outlier.size = 0.8, 
    outlier.alpha = 0.3, 
    linewidth = 0.6) +
  
  # FACET GRID: Rows = Metrics, Columns = Locations
  # This allows for easy vertical comparison of metrics across sites
  facet_grid(Parameters ~ Location, 
             scales = 'free_y') +
  # Styling
  scale_fill_manual(values = scenario_colors) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.05, 0.1))) +
  
  labs(
    y = "Network Metric Value",
    x = NULL,
  ) +
  
  theme_bw() + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 30),
    text = element_text(size = 30, family = "sans"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 20),
    strip.background = element_rect(fill = '#EEEBD3', color = "black"),
    strip.text = element_text(size = 20, face = "bold"),
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines"))
combined_summary_plot
plotpath <- file.path(out_path, "Supp_Figure_12_27_05_26.png")
ggsave(plotpath, combined_summary_plot, width = 20, height = 24)   


#-------------------------------------------------------------------------------------------------------------------------
## SUPPLEMENTARY FIGURE 13 ----
#-------------------------------------------------------------------------------------------------------------------------

# Topic: degree distributions for all four locations

# Open degree distributions of models for both seed numbers
all_degrees_Habi <- read_csv("Net_Sens/Graphs/OtherLocations/Habi_all_seeds_degree_distributions.csv") %>%
  mutate(Network = 'Habi')
all_degrees_Hepang <- read_csv("Net_Sens/Graphs/OtherLocations/Hepang_all_seeds_degree_distributions.csv") %>%
  mutate(Network = 'Hepang')
all_degrees_Sabaneta <- read_csv("Net_Sens/Graphs/OtherLocations/Sabaneta_all_seeds_degree_distributions.csv") %>%
  mutate(Network = 'Sabaneta')
all_degrees_Romana <- read_csv("Net_Sens/Graphs/OtherLocations/Romana_all_seeds_degree_distributions.csv") %>%
  mutate(Network = 'Romana')

# Combine 
all_seeds_degree_distributions <- all_degrees_Habi %>%
  bind_rows(all_degrees_Hepang, all_degrees_Romana, all_degrees_Sabaneta) %>%
  mutate(Network = trimws(Network)) %>%
  mutate(Graph_Type = case_when(
    type == "dcsbm"     ~ "DCSBM",
    type == "sbm"       ~ "SBM",
    type == "empirical" ~ "Empirical",
    type == 'ncrg'      ~ "NCRG",
    type == "spatial"   ~ "SENCA",
    TRUE                ~ "ERM"
  )) %>%
  mutate(Graph_Type = factor(Graph_Type, levels = c("SBM",
                                                    'ERM',
                                                    'SENCA',
                                                    'NCRG',
                                                    'DCSBM', 
                                                    'Empirical')))

# Summarise 
degree_summarized <- all_seeds_degree_distributions %>%
  group_by(Network, Graph_Type, degree) %>%
  summarise(
    q_25 = quantile(frequency, 0.25, names = FALSE, na.rm = TRUE),
    median = median(frequency, na.rm = TRUE),
    mean = mean(frequency, na.rm=TRUE),
    q_75 = quantile(frequency, 0.75, names = FALSE, na.rm = TRUE),
    .groups = "drop" 
  )

# Colour palette
scenario_colors <- c(
  "ERM"       = "#FAA76B",
  "NCRG"      = "#5C8A6F",
  "SBM"       = "#FCD2A2",
  "SENCA"     = "#AB4B54",
  "DCSBM"     = "#99B2DD",
  "Empirical" = "#495AA5")

# Plot 2: degree distributions per network
degree_plot <- ggplot(degree_summarized, aes(x = degree, 
                                             y = median, 
                                             group = Graph_Type, 
                                             color = Graph_Type, 
                                             fill = Graph_Type)) +
  
  # Ribbon first, then Line so the line sits on top
  geom_ribbon(aes(ymin = q_25, ymax = q_75), alpha = 0.2, color = NA) + 
  geom_line(linewidth = 1) +
  
  # Labels
  labs(
    y = "Frequency",
    x = "Degree (\u03BA)",
    color = "Graph Type",
    fill = "Graph Type") +
  
  # Scale
  scale_color_manual(
    values = scenario_colors,
    breaks = c("DCSBM", "Empirical", "ERM", "SBM", "SENCA")
  ) +
  scale_fill_manual(
    values = scenario_colors,
    breaks = c("DCSBM", "Empirical", "ERM", "SBM", "SENCA")
  ) +
  
  #scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  # scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  #scale_x_log10() +
  
  # Facet
  facet_wrap(~Network, scales = 'free', ncol = 2, axis.labels = "margins") +
  
  # Theme
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    text = element_text(size = 20, family = "sans"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    strip.background = element_rect(fill = '#EEEBD3', color = "black"),
    strip.text = element_text(size = 20),
    panel.grid.major = element_line(color = "grey90"), # Visible major grid
    panel.grid.minor = element_blank(),               # Remove distracting minor lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.spacing = unit(1.5, "lines"))
degree_plot
plotpath <- file.path(out_path, "Supp_Figure_13_27_05_26.png")
ggsave(plotpath, degree_plot, width = 12, height = 6)  


#-------------------------------------------------------------------------------------------------------------------------
## SUPPLEMENTARY FIGURE 14 ----
#-------------------------------------------------------------------------------------------------------------------------


# Topic: KS distances for all graph types and four locations

# Import
ks_distances_Habi <- read_csv("Net_Sens/Graphs/OtherLocations/Habi_ks_distances_all_sims.csv") %>%
  mutate(Network = 'Habi')
ks_distances_Hepang <- read_csv("Net_Sens/Graphs/OtherLocations/Hepang_ks_distances_all_sims.csv") %>%
  mutate(Network = 'Hepang')
ks_distances_Sabaneta <- read_csv("Net_Sens/Graphs/OtherLocations/Sabaneta_ks_distances_all_sims.csv") %>%
  mutate(Network = 'Sabaneta')
ks_distances_Romana <- read_csv("Net_Sens/Graphs/OtherLocations/Romana_ks_distances_all_sims.csv") %>%
  mutate(Network = 'Romana')

# Combine 
all_ks_distances <- ks_distances_Habi %>%
  bind_rows(ks_distances_Hepang, ks_distances_Sabaneta, ks_distances_Romana) %>%
  mutate(Graph = ifelse(graph_type == 'dcsbm', 'DCSBM',
                        ifelse(graph_type == 'sbm', 'SBM',
                               ifelse(graph_type == 'spatial', 'SENCA',
                                      ifelse(graph_type == 'ncrg', "NCRG", "ERM"))))) %>%
  mutate(Graph_Type = factor(Graph, levels = c("ERM",
                                               'SBM',
                                               'SENCA', 
                                               'NCRG',
                                               'DCSBM', 
                                               'Empirical')))
ks_summary <- ggplot(all_ks_distances, aes(x = Network, 
                                           y = ks_dist, 
                                           fill = Graph_Type)) +
  geom_boxplot() +
  labs(x = "Graph Type",
       y = "K-S distance") +
  theme_bw() +
  
  # Color Palette
  scale_color_manual(
    values = scenario_colors,
    breaks = c("DCSBM", "Empirical", "ERM",'NCRG', "SBM", "SENCA")
  ) +
  scale_fill_manual(
    values = scenario_colors,
    breaks = c("DCSBM", "Empirical", "ERM", 'NCRG', "SBM", "SENCA")
  ) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    text = element_text(size = 20, family = "sans"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    strip.background = element_rect(fill = '#EEEBD3', color = "black"),
    strip.text = element_text(size = 20),
    panel.grid.major = element_line(color = "grey90"), # Visible major grid
    panel.grid.minor = element_blank(),               # Remove distracting minor lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.spacing = unit(1.5, "lines"))
ks_summary
plotpath <- file.path(out_path, "Supp_Figure_14_27_05_26.png")
ggsave(plotpath, ks_summary, width = 12, height = 6)  


##########################################################################################################################

### ALL MANUSCRIPT SI TABLES

#-------------------------------------------------------------------------------------------------------------------------
## SUPPLEMENTARY TABLE 1  ----
#-------------------------------------------------------------------------------------------------------------------------

files <- list.files(
  path = "Net_Sens", 
  pattern = "^summary_check\\.csv$", 
  full.names = TRUE, 
  recursive = TRUE
) 
summary_checks <- do.call(rbind, lapply(files, read.csv))

# Leave out spatial optimization parameters, select only 1 replicate since identical
summary_checks_no_opt <- summary_checks %>% filter(!grepl('Opt', Params), 
                                                   Replicate == 1)

# Summarize convergence
convergence_summary <- summary_checks_no_opt %>%
  filter(Params %in% c("Avg_Degree", "Avg_Betweenness")) %>%
  filter(Graph != "empirical_graph") %>%
  group_by(Graph, Params, SeedTotal) %>%
  summarise(
    mean_val = mean(Values),
    se       = sd(Values) / sqrt(n()),
    rse_pct  = (se / mean_val) * 100,
    ci_width = 1.96 * se,
    .groups  = "drop") %>%
  
  # Sort to ensure the lag() function calculates the diff between the right rows
  arrange(Graph, Params, SeedTotal) %>%
  group_by(Graph, Params) %>%
  mutate(
    mean_drift_pct = (abs(mean_val - lag(mean_val)) / lag(mean_val)) * 100,
    precision_gain = (1 - (se / lag(se))) * 100)


#-------------------------------------------------------------------------------------------------------------------------
## SUPPLEMENTARY TABLE 3  ----
#-------------------------------------------------------------------------------------------------------------------------

size_and_duration <- read_csv("Outfiles/Randomized_Summary_SEIR_3000000_mseed_100_01June26.csv") %>%
  #filter(total_infections > 0) %>%
  select(graph_type, graph_idx, beta, delta, sigma, seed, peak_deaths, duration_days, total_infections) %>%
  mutate(Graph = ifelse(graph_type == 'dcsbm', 'DCSBM',
                        ifelse(graph_type == 'empirical', 'Empirical',
                               ifelse(graph_type == 'random', 'ERM',
                                      ifelse(graph_type == 'sbm', 'SBM', 
                                             ifelse(graph_type == 'newclust_graph', 'NCRG', 'SENCA'))))),
         final_size = peak_deaths/235) %>%
  filter(beta == 0.01 | beta == 0.1 | beta == 0.2) %>%
  rename_with(~c("duration"), c(duration_days))
size_and_duration$Graph <- factor(size_and_duration$Graph, 
                                  levels = c("DCSBM", "Empirical", "ERM", 'NCRG', 'SBM', "SENCA"))

# Apply to your data table by group
setDT(size_and_duration)
size_and_duration[, dynamic_threshold := find_outbreak_threshold_sensitive(values = final_size,
                                                                           adjust_bw = 0.3), 
                  by = .(Graph, beta)]

# Classify based on the specific threshold found for that group
size_and_duration[, outbreak_type := ifelse(final_size <= dynamic_threshold, 
                                            "Minor", 
                                            "Major")]
table(size_and_duration$outbreak_type)


supp_table3 <- size_and_duration %>% 
  group_by(Graph, outbreak_type, beta) %>%
  summarize(min = min(final_size),
            mean = mean(final_size),
            median = median(final_size),
            sd = sd(final_size),
            max = max(final_size),
            min_team = min(duration),
            mean_time = mean(duration),
            median_time = median(duration),
            sd_time = sd(duration),
            max_time = max(duration))


### END OF SCRIPT