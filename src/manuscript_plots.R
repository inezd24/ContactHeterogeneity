###########################################################################################################################

# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
# Created: March 2026; Last edited: September 2026
# This script creates the exact plots and tables used in the manuscript

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
library(dplyr)  
library(tidyr) 
library(ggplot2)  
library(readr)   
library(ggh4x)   
library(scales)   
library(igraph)  
library(ggraph)  
library(cowplot)  
library(data.table)
library(forcats)

# Set local root directory 
LOCAL_ROOT_DIR <- "/scicore/home/chitnis/derkx0000/GraphComparison"
setwd(LOCAL_ROOT_DIR)

# Output folders
out_path <- file.path(LOCAL_ROOT_DIR, "Manuscript")
ifelse(!dir.exists(out_path),
       dir.create(out_path), FALSE)

# Directory of files
files_directory <- "/scicore/home/chitnis/derkx0000/GraphComparison/Outfiles/"

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


### ALL MANUSCRIPT IMAGES

# ------------------------------------------------------------------------------
## FIGURE 1 ----
# ------------------------------------------------------------------------------

# Open ensemble of graphs for seed = 5000
MASTER_ENSEMBLE_FOR_DM <- readRDS("~/GraphComparison/Net_Sens/Graphs/Seed_5000/MASTER_ENSEMBLE_FOR_DM.rds")

# All empirical graphs are identical. Choose any task. 
empirical <- readRDS("~/GraphComparison/Net_Sens/Graphs/Seed_5000/vetted_graphs_task_3981_seed_5000.rds")

# Save total node number 
nodes <- length(V(empirical$anchor))

# Set models and colors
models <- c("spatial", "sbm", 'dcsbm', 'random', 'newclust_graph')
colors <- c("#FAA76B","#FCD2A2", "#AB4B54", "#99B2DD", "#5C8A6F")
stroke_colors <- c("#CB8350", "#C4A27B", "#78343A", "#4F6DA0", "#3B604B")

# Get list of graph types
graph_list <- Map(function(model, col, stroke_col) {
  list(graph = MASTER_ENSEMBLE_FOR_DM[[model]][[3]]$graph, 
       colour = col, 
       stroke_colour = stroke_col)
}, models, colors, stroke_colors)

# Add empirical graph and rename
graph_list[[6]] <- list(graph = empirical$anchor, 
                        colour = "#495AA5", 
                        stroke_colour = "#344078")
names(graph_list) <- c(models, "empirical")

# Reorder
graph_list <- graph_list[c("dcsbm","empirical", "random", "newclust_graph", "sbm", "spatial")]

# Function to plot a single graph
plot_graph <- function(list_item, title_text) {
  
  # Extract graph and colour from list to match publication style
  network <- list_item$graph
  node_colour <- list_item$colour
  stroke_colour <- list_item$stroke_colour
  
  # Set seed for identical graph each time
  set.seed(42)
  
  # Calculate degree
  d <- igraph::degree(network)
  
  # Rescale degree for node size in graph
  node_sizes <-  rescale(sqrt(d), to = c(1, 10))
  
  # Plot
  ggraph(network, layout = "nicely") +
    geom_edge_link(edge_color = "gray70", alpha = 0.7) +
    geom_node_point(fill = node_colour, shape = 21, color = stroke_colour, 
                    stroke = 0.8, size = node_sizes, alpha = 0.8) +
    
    # Faceting (not really, but kind of)
    labs(title = title_text) + 
    theme_graph() +
    theme(
      plot.title = element_text(
        hjust = 0.5, 
        size = 20, 
        face = "plain",
        family = 'sans',
        margin = margin(b = 10),
        background = ggplot2::element_rect(fill = "white", 
                                           colour = "black", 
                                           linewidth = 0.5)),
      plot.margin = margin(15, 15, 15, 15),
      panel.border = element_rect(colour = "black", 
                                  fill = NA, 
                                  linewidth = 0.5))
}

# Titles for plots
titles <- c('DCSBM', 'Empirical', "ERM","NCRG", 'SBM', "SENCA")

# Apply function to all plots
styled_plots <- Map(plot_graph, graph_list, titles)

# 3. Use cowplot::plot_grid with horizontal and vertical alignment ('vh')
all_graphs <- cowplot::plot_grid(plot_grid(
  plotlist = styled_plots, 
  ncol = 3, 
  labels = c('(A)', '(B)', '(C)', '(D)', '(E)'),
  label_size = 20,
  label_fontface = "bold",
  label_fontfamily = "sans",
  align = 'vh')) +      
  theme(plot.background = element_rect(fill = "white", colour = NA)) 
all_graphs
plotpath <- file.path(out_path, "Figure1_27_05_26.png")
ggsave(plotpath, all_graphs, width = 18, height = 12) 

# Remove files
rm(MASTER_ENSEMBLE_FOR_DM)
rm(empirical)
rm(all_graphs)


# ------------------------------------------------------------------------------
## FIGURE 2 ----
# ------------------------------------------------------------------------------

# Only use 5000 seeds
summary_check <- read_csv("Net_Sens/Seed_5000/summary_check.csv")

# Leave out spatial optimization parameters, select only 1 replicate since identical
summary_check_no_opt <- summary_check %>% filter(!grepl('Opt', Params), 
                                                 Replicate == 1)

# Rename graph type variable levels
summary_check_no_opt$Graph_type <- ifelse(summary_check_no_opt$Graph == 'dcsbm_graph', 
                                          'DCSBM',
                                          ifelse(summary_check_no_opt$Graph == 'empirical_graph', 
                                                 'Empirical',
                                                 ifelse(summary_check_no_opt$Graph == 'random_graph', 
                                                        "ERM",
                                                        ifelse(summary_check_no_opt$Graph == 'sbm_graph', 
                                                               'SBM', 
                                                               ifelse(summary_check_no_opt$Graph == 'newclust_graph', "NCRG",
                                                                      "SENCA")))))

# Relevel the graph type variable
summary_check_no_opt$Graph_type <- factor(summary_check_no_opt$Graph_type, 
                                          levels = c("DCSBM", 'Empirical', "ERM","NCRG", "SBM", "SENCA"))

# Rename degree parameters for plotting
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

# Relevel Parameters variable for plotting 
summary_check_no_opt$Parameters <- factor(summary_check_no_opt$Parameters, 
                                          levels = c('Density', 
                                                     'Global Clustering', 
                                                     "Mean Betweenness", 
                                                     "Mean Degree", 
                                                     'Median Degree', 
                                                     "Maximum Degree")
                                          )

# Summarize by mean, max, median, min
parameter_summary <- summary_check_no_opt %>% 
  group_by(Graph_type, Params) %>% 
  summarise(mean = mean(Values),
            max = max(Values),
            median = median(Values),
            min = min(Values))

# Colour palette
scenario_colors <- c(
  "ERM"       = "#FAA76B",
  "SBM"       = "#FCD2A2",
  "SENCA"     = "#AB4B54",
  "DCSBM"     = "#99B2DD",
  "Empirical" = "#495AA5",
  "NCRG"      = "#5C8A6F")

# Plot 2: network topology overview per seed number
summary_plot <- ggplot(summary_check_no_opt, aes(x = Graph_type, 
                                                 y = Values, 
                                                 fill = Graph_type)) +
  
  # Boxplots: Slightly narrower with subtle outliers
  geom_boxplot(
    width = 0.6, 
    outlier.size = 1, 
    outlier.alpha = 0.4, 
    linewidth = 0.7) +
  
  # Faceting 
  facet_wrap(
    ~ Parameters, 
    scales = 'free_y', 
    ncol = 3) +
  
  # Colors and Scales
  scale_fill_manual(
    values = scenario_colors,
    breaks = c("DCSBM", "Empirical","ERM","NCRG", "SBM", "SENCA")
  ) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.05, 0.1))) +
  
  # Labels 
  labs(
    y = "Network Metric Value",
    x = NULL, # Removed as the X-axis labels and legend make this redundant
    fill = 'Network Type') +
  
  # Theme
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
    panel.spacing = unit(1.5, "lines"))
summary_plot
plotpath <- file.path(out_path, "Figure2_27_05_26.png")
ggsave(plotpath, summary_plot, width = 18, height = 16)  

# Remove large files
rm(summary_check)
rm(summary_check_no_opt)
rm(summary_plot)
# ------------------------------------------------------------------------------
## TABLE 2 ----
# ------------------------------------------------------------------------------

# Open data
all_seeds_degree_distributions_5000 <- read_csv("Net_Sens/Graphs/Seed_5000/all_seeds_degree_distributions.csv") %>%
  mutate(Seed_run = '5000_seeds')

# Summarize 
degree_summary <- all_seeds_degree_distributions_5000 %>%
  uncount(weights = frequency) %>%
  group_by(type, task_id) %>%
  summarise(
    min_degree = min(degree),
    q25 = quantile(degree, 0.25),
    median = median(degree),
    mean = mean(degree),
    q75 = quantile(degree, 0.75),
    max_degree = max(degree),
    .groups = "drop"
  ) %>%
  group_by(type) %>%
  summarise(
    min = mean(min_degree),
    q25 = mean(q25),
    median = mean(median),
    mean = mean(mean),
    q75 = mean(q75),
    max = mean(max_degree),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ sprintf("%.3f", .x)))

# These results are summarized in table 2 of the manuscript
degree_summary

# Remove large files
rm(all_seeds_degree_distributions_5000)


# ------------------------------------------------------------------------------
## FIGURE 3 ----
# ------------------------------------------------------------------------------

# Open data
all_seeds_degree_distributions_5000 <- read_csv("Net_Sens/Graphs/Seed_5000/all_seeds_degree_distributions.csv") %>%
  mutate(Seed_run = '5000_seeds')

# Summarize data
degree_summarized <- all_seeds_degree_distributions_5000 %>%
  rename(Type = type) %>%
  group_by(Type, degree) %>%
  summarise(
    q_25   = quantile(frequency, 0.25, names = FALSE),
    median = median(frequency),
    q_75   = quantile(frequency, 0.75, names = FALSE),
    .groups = "drop"                                      # drop grouping cleanly
  ) %>%
  mutate(Graph = case_when(
    Type == "dcsbm"    ~ "DCSBM",
    Type == "sbm"      ~ "SBM",
    Type == "empirical" ~ "Empirical",
    Type == "spatial"  ~ "SENCA",
    Type == 'newclust_graph' ~ "NCRG",
    TRUE               ~ "ERM"
  )) %>%
  mutate(Graph = fct_relevel(Graph, "Empirical", after = Inf))
  
# Colour palette
scenario_colors <- c(
  "ERM"       = "#FAA76B",
  "SBM"       = "#FCD2A2",
  "SENCA"     = "#AB4B54",
  "DCSBM"     = "#99B2DD",
  "Empirical" = "#495AA5",
  "NCRG"      = "#5C8A6F")

# Plot
degree_plot <- ggplot(degree_summarized, aes(x = degree, 
                                             y = median, 
                                             group = Graph, 
                                             color = Graph, 
                                             fill = Graph)) +
  
  # Ribbon first, then Line so the line sits on top
  geom_ribbon(aes(ymin = q_25, ymax = q_75), alpha = 0.2, color = NA) + 
  geom_line(aes(size=Graph)) +
  scale_size_manual(values=c(1, 1, 1, 1, 1, 1.5), guide = "none") +
  
  # Labels
  labs(
    y = "Frequency",
    x = "Degree (\u03BA)",
    color = "Network Type",
    fill = "Network Type") +
  
  # Scale
  scale_color_manual(
    values = scenario_colors,
    breaks = c("DCSBM", "Empirical", "ERM","NCRG", "SBM", "SENCA")
  ) +
  scale_fill_manual(
    values = scenario_colors,
    breaks = c("DCSBM", "Empirical", "ERM","NCRG", "SBM", "SENCA")
  ) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  #scale_x_log10() +
  
  # Theme
  theme_bw() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 20),
    text = element_text(size = 20),
    # Full Grid Look as requested
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    # Polish legend
    legend.title = element_text(),
    legend.background = element_blank(),
    plot.subtitle = element_text(size = 12, color = "grey30", margin = margin(b = 10))) +
  guides(
    color = guide_legend(reverse = FALSE),
    fill  = guide_legend(reverse = FALSE)
  )
degree_plot

library(patchwork)
# 1. Create the inset plot zoomed to x = 50
inset_plot <- degree_plot +
  coord_cartesian(xlim = c(0, 30)) +
  labs(x = NULL, y = NULL) +  # Removes axis titles to keep the inset clean
  theme(
    legend.position = "none", # Hide the redundant legend
    text = element_text(size = 11), # Shrink text size for the inset
    plot.background = element_rect(fill = "white", color = NA) # Keep background opaque
  )
inset_plot

# 2. Place the inset in the top-right corner of the main panel
final_plot <- ggdraw(degree_plot) +
  draw_plot(
    inset_plot, 
    x = 0.45,    # X-position of the inset's bottom-left corner
    y = 0.55,    # Y-position of the inset's bottom-left corner
    width = 0.40, # Width of the inset relative to canvas
    height = 0.40 # Height of the inset relative to canvas
  )
final_plot
plotpath <- file.path(out_path, "Figure3_27_05_26.png")
ggsave(plotpath, final_plot, width = 16, height = 10)

# Remove large files
rm(all_seeds_degree_distributions_5000)
rm(final_plot)
rm(degree_plot)


# ------------------------------------------------------------------------------
## FIGURE 4 ----
# ------------------------------------------------------------------------------

# Import all files: SIS and SEIR output
Summary_SIS_ran <- read_csv("Outfiles/Randomized_Summary_SIS_6000000_mseed_100_01June26.csv")
Summary_SEIR_ran <- read_csv("Outfiles/Randomized_Summary_SEIR_3000000_mseed_100_01June26.csv")

# Import empirical graph and calculate node total
empirical <- readRDS("~/GraphComparison/Net_Sens/Graphs/Seed_5000/vetted_graphs_task_3981_seed_5000.rds")
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
                                                 ifelse(model_cases$graph_type == 'newclust_graph','NCRG', 'ERM')))))
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
  "ERM"      = "#FAA76B",
  "SBM"      = "#FCD2A2",
  "SENCA"    = "#AB4B54",
  "DCSBM"    = "#99B2DD",
  "Empirical" = "#495AA5",
  "NCRG"      = "#5C8A6F")

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
    breaks = c("DCSBM", "Empirical", "ERM", "NCRG", "SBM", "SENCA")
  ) +
  scale_fill_manual(
    values = scenario_colors,
    breaks = c("DCSBM", "Empirical", "ERM", "NCRG", "SBM", "SENCA")
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
plotpath <- file.path(out_path, "Figure4_27_05_26.png")
ggsave(plotpath, cases_plot, width = 12, height = 16)   

# Remove large files
rm(Summary_SIS_ran)
rm(Summary_SEIR_ran)
rm(SIS_cases)
rm(SEIR_cases)
rm(model_cases)
rm(cases_plot)
  
  
# ------------------------------------------------------------------------------
## FIGURE 5 ----
# ------------------------------------------------------------------------------

# Import SEIR data and mutate names
size_and_duration <- read_csv("Outfiles/Randomized_Summary_SEIR_3000000_mseed_100_01June26.csv") %>%
  select(graph_type, graph_idx, beta, delta, sigma, seed, peak_deaths, duration_days, total_infections) %>%
  mutate(Graph = ifelse(graph_type == 'dcsbm', 'DCSBM',
                        ifelse(graph_type == 'empirical', 'Empirical',
                               ifelse(graph_type == 'random', 'ERM',
                                      ifelse(graph_type == 'sbm', 'SBM', 
                                             ifelse(graph_type == 'newclust_graph', 'NCRG', 'SENCA'))))),
         final_size = peak_deaths/235) %>%
  
  # Retain only three transmission rate values for efficient plotting
  filter(beta == 0.01 | beta == 0.1 | beta == 0.2) %>%
  rename_with(~c("duration"), c(duration_days))

# Relevel graph factor
size_and_duration$Graph <- factor(size_and_duration$Graph, 
                                  levels = c("DCSBM", "Empirical", "ERM", 'NCRG', 'SBM', "SENCA"))

# Set df as data table
setDT(size_and_duration)

# Use helper function 'find_outbreak_threshold_sensitive' to identify outbreak thresholds
# We use 0.3 for bandwidth adjustment. See 'optimal_bandwidth.R' and methods of corresponding
# manuscript for explanation of this. 
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
  # (GraphType) . (outbreak_type)
  "ERM.Minor"         = "#FAA76B", "ERM.Major"         = "#CB8350",
  "SBM.Minor"         = "#FCD2A2", "SBM.Major"         = "#C4A27B",
  "SENCA.Minor"       = "#AB4B54", "SENCA.Major"       = "#78343A", 
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
    panel.spacing = unit(1.5, "lines")) #+
guides(size=FALSE)+guides(color = guide_legend(override.aes = list(size = 4)))
comparison_plot
plotpath <- file.path(out_path, "Figure5_27_05_26.png")
ggsave(plotpath, comparison_plot, width = 18, height = 12)  
write.csv(size_and_duration, "/Outfiles/size_and_duration_SEIR_3000000.csv", row.names = FALSE)

# Remove large files
rm(size_and_duration)
rm(comparison_plot)


# -------------------------------------------------------
## FIGURE 6 ----
# -------------------------------------------------------

# Load data
Hepang_SEIR <- read_csv("Outfiles/Hepang_Summary_SEIR_3000000_mseed_100_04June26.csv") %>%
  mutate(Location = 'Hepang',
         final_size = peak_deaths/61)
Habi_SEIR <- read_csv("Outfiles/Habi_Summary_SEIR_3000000_mseed_100_04June26.csv") %>%
  mutate(Location = 'Habi',
         final_size = peak_deaths/81)
Sabaneta_SEIR <- read_csv("Outfiles/Sabaneta_Summary_SEIR_3000000_mseed_100_04June26.csv") %>%
  mutate(Location = 'Sabaneta',
         final_size = peak_deaths/121)
Romana_SEIR <- read_csv("Outfiles/Romana_Summary_SEIR_3000000_mseed_100_04June26.csv") %>%
  mutate(Location = 'Romana',
         final_size = peak_deaths/57)

# Combine and convert
otherlocations_SEIR <- rbind(Hepang_SEIR, Habi_SEIR, Sabaneta_SEIR, Romana_SEIR)
setDT(otherlocations_SEIR)

# Apply outbreak threshold identification to data table by group
# Requires helper function 1
otherlocations_SEIR[, dynamic_threshold := find_outbreak_threshold_sensitive(
  values = final_size, 
  adjust_bw = 0.3
), 
by = .(Location, graph_type, beta)]

# Classify based on the specific threshold found for that group
otherlocations_SEIR[, outbreak_type := ifelse(final_size <= dynamic_threshold, 
                                              "Minor", 
                                              "Major")]
table(otherlocations_SEIR$outbreak_type)

# SEIR summary
SEIR_cases <- otherlocations_SEIR %>%
  
  # Filtering out failed outbreaks (peak > 1) to focus on established epidemics
  filter(peak_active_cases > 0) %>%
  mutate(N = ifelse(Location == 'Habi', 81,
                    ifelse(Location == 'Hepang', 61, 
                           ifelse(Location == 'Sabaneta', 121, 57)))) %>%
  pivot_longer(cols = c(peak_active_cases, avg_active_cases, total_infections), names_to = "measure", values_to = 'value') %>%
  group_by(Location, graph_type, outbreak_type, beta, delta, measure) %>%
  summarise(min = min(value), 
            q_25 = quantile(value/N, 0.25, names = FALSE),
            median = median(value/N),
            mean = mean(value/N),
            q_75 = quantile(value/N, 0.75, names = FALSE),
            max = max(value/N), .groups = 'drop')

# Adapt graph type names
SEIR_cases$Graph <- ifelse(SEIR_cases$graph_type == 'dcsbm', 'DCSBM',
                           ifelse(SEIR_cases$graph_type == 'sbm', 'SBM',
                                  ifelse(SEIR_cases$graph_type == 'empirical', 'Empirical',
                                         ifelse(SEIR_cases$graph_type == 'spatial', 'SENCA', 
                                                ifelse(SEIR_cases$graph_type == 'ncrg', "NCRG", "ERM")))))

# Relevel factor
SEIR_cases$Graph <- factor(SEIR_cases$Graph, 
                            levels = c("ERM",
                                       'SBM',
                                       'SENCA',
                                       'NCRG',
                                       'DCSBM', 
                                       'Empirical'))

# Colour palette
scenario_colors <- c(
  "ERM"      = "#FAA76B",
  "SBM"      = "#FCD2A2",
  "SENCA"    = "#AB4B54",
  "DCSBM"    = "#99B2DD",
  "Empirical" = "#495AA5",
  "NCRG"      = "#5C8A6F")

# Plot all cases
otherlocations_prevalence <- ggplot(SEIR_cases, aes(x = beta, color = Graph, fill = Graph, group = interaction(Graph, outbreak_type))) +
  
  # Ribbon for uncertainty (25th to 75th percentile)
  # Color is NA to avoid borders around the fill
  geom_ribbon(aes(ymin = q_25, ymax = q_75), alpha = 0.2, color = NA) +
  
  # Line for the median, mapped to linetype
  geom_line(aes(y = median, linetype = outbreak_type), linewidth = 1.1) +
  
  # Labels and Legends
  labs(
    x = "Transmission Rate (\u03B2)", 
    y = NULL, 
    color = "Graph Type", 
    fill = "Graph Type",
    linetype = "Outbreak Type") + 
  
  # Faceting
  facet_grid(measure ~ Location,
             scales = 'free',
             labeller = labeller(
               measure = c(
                 avg_active_cases = "Average Prevalence",
                 peak_active_cases = "Peak Prevalence",
                 total_infections = "Total Infections"))) + 
  
  # Manual Scales for Color and Fill
  scale_color_manual(
    values = scenario_colors,
    breaks = c("DCSBM", "Empirical", "ERM","NCRG", "SBM", "SENCA")
  ) +
  scale_fill_manual(
    values = scenario_colors,
    breaks = c("DCSBM", "Empirical", "ERM","NCRG", "SBM", "SENCA")
  ) +
  
  # Optional: Clean up linetypes (e.g., solid for primary, dashed for secondary)
  scale_linetype_manual(values = c("solid", "dashed")) +

  # Refined Theme
  theme_bw() + 
  theme(
    legend.position = "bottom",
    legend.box = "vertical", # Stack the Color and Linetype legends for clarity
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    text = element_text(size = 18, family = "sans"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    strip.background = element_rect(fill = '#EEEBD3', color = "black"),
    strip.text = element_text(size = 18),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.spacing = unit(1.2, "lines")
  )
otherlocations_prevalence
plotpath <- file.path(out_path, "Figure6_27_05_26.png")
ggsave(plotpath, otherlocations_prevalence, width = 12, height = 12)   

# Remove large files
rm(SEIR_cases)
rm(otherlocations_SEIR)
rm(Romana_SEIR)
rm(Sabaneta_SEIR)
rm(Habi_SEIR)
rm(Hepang_SEIR)
rm(otherlocations_prevalence)


# ------------------------------------------------------------------------------
## END OF SCRIPT ----
# ------------------------------------------------------------------------------
