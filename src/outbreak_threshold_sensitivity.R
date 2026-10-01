##########################################################################################################################

# Script by Inez Derkx, contact: inez.derkx@swisstph.ch
# In this script, we classify SEIR outbreaks as minor or major based on the distribution of final outbreak sizes, and
# run the bandwidth sensitivity analysis reported in the manuscript (Methods, outbreak classification).

# IMPORTANT NOTES:
# 1. The threshold is the lowest point of a kernel density estimate between the first and last significant peak of the
#    final-size distribution, calculated separately for each Graph x beta combination.
# 2. The main analysis uses a bandwidth adjustment of H_MAIN = 0.3. The sensitivity analysis repeats the classification
#    for every value in H_GRID and compares classifications across them.
# 3. 'adjust' in density() MULTIPLIES the default bandwidth (Silverman's rule of thumb, bw.nrd0); it is not the bandwidth
#    itself. The absolute bandwidths used per group are saved in the threshold table (column bw_abs).

##########################################################################################################################

### SET UP R ENVIRONMENT ###

# Set working directory
setwd("/scicore/home/chitnis/derkx0000/Derkx2026_publication")

# Load libraries
library(dplyr)
library(readr)
library(data.table)

# Settings
H_MAIN         <- 0.3                          # bandwidth adjustment for the main analysis
H_GRID         <- c(0.1, 0.2, 0.3, 0.4, 0.5)   # bandwidth adjustments for the sensitivity analysis
PEAK_THRESHOLD <- 0.02                         # peaks must exceed this fraction of the maximum density
N_POP          <- 235                          # population size used to convert deaths into a fraction
BETA_KEEP      <- c(0.01, 0.1, 0.2)            # transmission rates included in the analysis

IN_FILE  <- "Outfiles/Randomized_Summary_SEIR_3000000_mseed_100_01June26.csv"
OUT_DIR  <- "Outfiles"

stopifnot(H_MAIN %in% H_GRID)

##########################################################################################################################

### Helper functions

#-------------------------------------------------------------------------------------------------------------------------
# Helper function 1: find outbreak threshold
find_outbreak_threshold <- function(values, adjust_bw, peak_threshold = PEAK_THRESHOLD) {
  
  # This function finds the valley between minor and major outbreaks in a final-size distribution
  #' @param values: final outbreak sizes (fraction of population)
  #' @param adjust_bw: multiplier on the default bandwidth passed to density(adjust = ...)
  #' @param peak_threshold: minimum peak height relative to the maximum density
  #' @return list with the threshold, the method used, and the absolute bandwidth
  
  values <- values[!is.na(values)]
  if (length(unique(values)) < 5) {
    return(list(threshold = 0.05, method = "fallback_few_values", bw_abs = NA_real_))
  }
  
  # Kernel density estimate restricted to the 0-1 range of population fraction
  dens <- density(values, adjust = adjust_bw, from = 0, to = 1)
  
  # Local maxima: points higher than both neighbours
  peaks <- which(diff(sign(diff(dens$y))) == -2) + 1
  
  # Ignore small wiggles
  significant_peaks <- peaks[dens$y[peaks] > max(dens$y) * peak_threshold]
  
  if (length(significant_peaks) >= 2) {
    
    # The valley is the lowest point between the first and the last significant peak
    first_p    <- significant_peaks[1]
    last_p     <- significant_peaks[length(significant_peaks)]
    valley_idx <- which.min(dens$y[first_p:last_p]) + first_p - 1
    
    list(threshold = dens$x[valley_idx], method = "valley", bw_abs = dens$bw)
    
  } else {
    
    # Fallback when the distribution is not bimodal (common at very low or very high beta)
    # 5% prevalence if there are large outbreaks, otherwise everything is minor
    list(threshold = if (max(values) > 0.1) 0.05 else 1.0,
         method    = "fallback_unimodal",
         bw_abs    = dens$bw)
  }
}
#-------------------------------------------------------------------------------------------------------------------------

##########################################################################################################################

### Preparation and execution

## Step 1: load data

# Read csv
size_and_duration <- read_csv(IN_FILE, show_col_types = FALSE) %>%
  select(graph_type, graph_idx, beta, delta, sigma, seed, peak_deaths, duration_days, total_infections) %>%
  mutate(Graph = case_when(graph_type == "dcsbm"          ~ "DCSBM",
                           graph_type == "empirical"      ~ "Empirical",
                           graph_type == "random"         ~ "ERM",
                           graph_type == "sbm"            ~ "SBM",
                           graph_type == "newclust_graph" ~ "NCRG",
                           TRUE                           ~ "SENCA"),
         final_size = peak_deaths / N_POP) %>%
  filter(beta %in% BETA_KEEP) %>%
  rename(duration = duration_days)

# Change graph names
size_and_duration$Graph <- factor(size_and_duration$Graph,
                                  levels = c("DCSBM", "Empirical", "ERM", "NCRG", "SBM", "SENCA"))

# Turn into data table
setDT(size_and_duration)
cat(sprintf("Simulations included: %d\n", nrow(size_and_duration)))


## Step 2: Thresholds for every Graph x beta x bandwidth 


# Find outbreak thresholds 
thresholds <- size_and_duration[, {
  res <- lapply(H_GRID, function(h) find_outbreak_threshold(final_size, adjust_bw = h))
  .(h         = H_GRID,
    threshold = vapply(res, `[[`, numeric(1),   "threshold"),
    method    = vapply(res, `[[`, character(1), "method"),
    bw_abs    = vapply(res, `[[`, numeric(1),   "bw_abs"))
}, by = .(Graph, beta)]

# Report groups where the valley could not be found
fallbacks <- thresholds[method != "valley"]
if (nrow(fallbacks) > 0) {
  cat("\nGroups using a fallback threshold:\n")
  print(fallbacks)
} else {
  cat("\nValley-based threshold found for all Graph x beta x h combinations.\n")
}


## Step 3: Classify outbreaks for every bandwidth ###

                
# Set threshold and column names              
h_label    <- function(h) sprintf("h%.1f", h)
thr_cols   <- paste0("thr_",   h_label(H_GRID))
class_cols <- paste0("class_", h_label(H_GRID))

# Wide threshold table: one column per bandwidth, joined to each simulation by Graph x beta
thr_wide <- dcast(thresholds, Graph + beta ~ paste0("thr_", h_label(h)), value.var = "threshold")
size_and_duration <- thr_wide[size_and_duration, on = .(Graph, beta)]

# Loop over grid
for (i in seq_along(H_GRID)) {
  size_and_duration[, (class_cols[i]) := fifelse(final_size <= get(thr_cols[i]), "Minor", "Major")]
}

# Main analysis classification
main_col <- paste0("class_", h_label(H_MAIN))
size_and_duration[, dynamic_threshold := get(paste0("thr_", h_label(H_MAIN)))]
size_and_duration[, outbreak_type     := get(main_col)]

# Print output
cat("\nMain analysis (h =", H_MAIN, "):\n")
print(table(size_and_duration$outbreak_type)) # This is the classification used in the manuscript

                
## Step 4: Sensitivity analysis 

                
# A simulation is stable if it receives the same classification at every bandwidth in H_GRID
size_and_duration[, stable_all := Reduce(`&`, lapply(.SD, function(x) x == .SD[[1]])),
                  .SDcols = class_cols]

# Overall stability across the full bandwidth range
pct_stable_overall <- 100 * mean(size_and_duration$stable_all)

# Stability per Graph x beta combination
stability_by_group <- size_and_duration[, .(n          = .N,
                                            n_unstable = sum(!stable_all),
                                            pct_stable = 100 * mean(stable_all)),
                                        by = .(Graph, beta)][order(beta, Graph)]

# Where are the differences concentrated?
unstable_by_beta <- size_and_duration[, .(n_unstable = sum(!stable_all)), by = beta][order(beta)]
unstable_by_beta[, pct_of_all_unstable := 100 * n_unstable / sum(n_unstable)]
lowest_beta       <- min(BETA_KEEP)
min_stable_other  <- stability_by_group[beta != lowest_beta, min(pct_stable)]

# Pairwise agreement between bandwidths (% of simulations with identical classification)
pairwise <- CJ(h1 = H_GRID, h2 = H_GRID)[h1 < h2]
pairwise[, pct_agree := mapply(function(a, b) {
  100 * mean(size_and_duration[[paste0("class_", h_label(a))]] ==
               size_and_duration[[paste0("class_", h_label(b))]])
}, h1, h2)]
pairwise[, n_differ := mapply(function(a, b) {
  sum(size_and_duration[[paste0("class_", h_label(a))]] !=
        size_and_duration[[paste0("class_", h_label(b))]])
}, h1, h2)]

# Are the classifications at h = 0.2 and h = 0.3 identical?
n_diff_02_03 <- pairwise[h1 == 0.2 & h2 == 0.3, n_differ]

# Truncate (not round) percentages reported with ">=" so the claim is never overstated
floor_to <- function(x, digits) floor(x * 10^digits) / 10^digits


## Step 5:  report numbers
                                                     

cat("\nReport sensitivity analysis output:\n")
cat(sprintf("Bandwidth adjustments tested: %s (main analysis: %.1f)\n",
            paste(H_GRID, collapse = ", "), H_MAIN))
cat(sprintf("Simulations with identical classification across h = %.1f-%.1f: %.1f%% (%d of %d)\n",
            min(H_GRID), max(H_GRID), pct_stable_overall,
            sum(size_and_duration$stable_all), nrow(size_and_duration)))
cat(sprintf("Minimum stability across all Graph x beta combinations with beta != %s: >= %.2f%%\n",
            lowest_beta, floor_to(min_stable_other, 2)))
cat(sprintf("Classifications at h = 0.2 and h = 0.3: %s (%d simulations differ)\n",
            if (n_diff_02_03 == 0) "IDENTICAL" else "NOT identical", n_diff_02_03))

cat("\nPairwise agreement between bandwidths:\n")
print(pairwise)



### END OF SCRIPT ###
