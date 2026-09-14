# Function to find optimal bandwidth
find_optimal_bw <- function(values, 
                            bw_range = seq(0.1, 0.8, by = 0.05),
                            peak_threshold = 0.02) {
  
  results <- lapply(bw_range, function(bw) {
    
    # 
    dens <- density(values, adjust = bw, from = 0, to = 1)
    
    is_peak <- diff(sign(diff(dens$y))) == -2
    peaks <- which(is_peak) + 1
    sig_peaks <- peaks[dens$y[peaks] > max(dens$y) * peak_threshold]
    n_peaks <- length(sig_peaks)
    
    # Valley depth ratio: how pronounced is the dip between the two main peaks?
    # Higher = clearer separation = better for threshold-finding
    valley_ratio <- NA_real_
    if (n_peaks >= 2) {
      first_p <- sig_peaks[1]
      last_p  <- sig_peaks[length(sig_peaks)]
      valley_y <- min(dens$y[first_p:last_p])
      peak_mean <- mean(dens$y[c(first_p, last_p)])
      valley_ratio <- 1 - (valley_y / peak_mean)  # 0 = flat, 1 = deep valley
    }
    
    list(bw = bw, 
         n_peaks = n_peaks, 
         valley_ratio = valley_ratio)
  })
  
  results_df <- do.call(rbind, lapply(results, as.data.frame))
  
  # Prefer bandwidths that give exactly 2 peaks, then maximise valley depth
  bimodal <- results_df[!is.na(results_df$valley_ratio) & results_df$n_peaks == 2, ]
  
  if (nrow(bimodal) > 0) {
    optimal_bw <- bimodal$bw[which.max(bimodal$valley_ratio)]
  } else {
    # Fall back to narrowest bw that gives >= 2 peaks, or middle of range
    multi <- results_df[!is.na(results_df$valley_ratio), ]
    optimal_bw <- if (nrow(multi) > 0) multi$bw[which.max(multi$valley_ratio)] else median(bw_range)
  }
  
  optimal_bw
}

# Function to find outbreak type threshold
find_outbreak_threshold_sensitive <- function(values, 
                                              bw_range = seq(0.1, 0.8, by = 0.05),
                                              peak_threshold = 0.02) {
  values <- na.omit(values)
  if (length(unique(values)) < 5) return(0.05)
  
  adjust_bw <- find_optimal_bw(values, bw_range, peak_threshold)
  cat(sprintf("Optimal bandwidth: %.2f\n", adjust_bw))
  
  dens <- density(values, adjust = adjust_bw, from = 0, to = 1)
  
  is_peak <- diff(sign(diff(dens$y))) == -2
  peaks <- which(is_peak) + 1
  significant_peaks <- peaks[dens$y[peaks] > max(dens$y) * peak_threshold]
  
  if (length(significant_peaks) >= 2) {
    first_p <- significant_peaks[1]
    last_p  <- significant_peaks[length(significant_peaks)]
    between_y <- dens$y[first_p:last_p]
    valley_idx <- which.min(between_y) + first_p - 1
    threshold  <- dens$x[valley_idx]
    cat("Threshold calculated.\n")
  } else {
    threshold <- if (max(values) > 0.1) 0.05 else 1.0
    cat("Threshold not calculated. Reverting to fallback.\n")
  }
  
  return(threshold)
}

# Alter data
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

# Apply function 
setDT(size_and_duration)
size_and_duration[, dynamic_threshold := find_outbreak_threshold_sensitive(values = final_size),
                  by = .(Graph, beta)]
hist(size_and_duration$dynamic_threshold)

# Classify based on the specific threshold found for that group
size_and_duration[, outbreak_type := ifelse(final_size <= dynamic_threshold, 
                                            "Minor", 
                                            "Major")]
table(size_and_duration$outbreak_type)
79718 / 151695


79500 / 151913