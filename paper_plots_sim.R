#!/usr/bin/env Rscript
# ============================================================================
# Paper Plots for Simulation Study: Split-Group Elbow Analysis
# ============================================================================
# Creates publication-quality plots showing RMSE and variance elbow curves
# for two randomly-split groups of replicates across multiple H_true values.
#
# Output: One plot per location with 2 columns (RMSE, Variance) × 7 rows (H_true)
# ============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(patchwork)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Location to analyze
LOCATION <- "NY"  # Options: "GOM" (Gulf of Mexico) or "NY" (New York)

# List of H_true values to analyze
H_TRUE_VALUES <- 2:7

# Random seed for group splitting (for reproducibility)
SPLIT_SEED <- 42

# Number of replicates per group
N_GROUP_A <- 50
N_GROUP_B <- 50

# Whether to include titles and subtitles on plots
SHOW_TITLES <- TRUE

# Output directory
OUTPUT_DIR <- file.path(paste0("results/sim_posterior_",LOCATION),"585")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("\n", strrep("=", 70), "\n")
cat("PAPER PLOTS: Split-Group Elbow Analysis\n")
cat(strrep("=", 70), "\n\n")
cat("Location:", LOCATION, "\n")
cat("H_true values:", paste(H_TRUE_VALUES, collapse = ", "), "\n")
cat("Group A size:", N_GROUP_A, "\n")
cat("Group B size:", N_GROUP_B, "\n")
cat("Random seed:", SPLIT_SEED, "\n\n")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Load consolidated results for a given H_true
load_results <- function(H_true, location) {
    base_dir <- file.path("results", "sim_posterior_GOM","585")
    if (location == "NY") {
        base_dir <- file.path("results", "sim_posterior_NY","585")
    }
    
    results_dir <- file.path(base_dir, paste0("Htrue_", H_true))
    
    # Check if directory exists
    if (!dir.exists(results_dir)) {
        cat("Warning: Directory not found:", results_dir, "\n")
        return(NULL)
    }
    
    # Load consolidated files
    metrics_file <- file.path(results_dir, "consolidated_metrics.rds")
    curves_file <- file.path(results_dir, "timewise_posterior_summaries.rds")
    per_rep_file <- file.path(results_dir, "per_replicate_posterior_stats.rds")
    
    if (!file.exists(metrics_file) || !file.exists(curves_file) || !file.exists(per_rep_file)) {
        cat("Warning: Missing consolidated files for H_true =", H_true, "\n")
        return(NULL)
    }
    
    metrics <- readRDS(metrics_file)
    curves <- readRDS(curves_file)
    per_replicate <- readRDS(per_rep_file)
    
    list(
        H_true = H_true,
        metrics = metrics,
        curves = curves,
        per_replicate = per_replicate
    )
}

#' Split replicates into two random groups
split_replicates <- function(all_reps, n_group_a, seed) {
    set.seed(seed)
    group_A <- sample(all_reps, size = n_group_a, replace = FALSE)
    group_B <- setdiff(all_reps, group_A)
    list(A = group_A, B = group_B)
}

#' Analyze split group for RMSE and variance
analyze_split_group <- function(metrics_df, curves_df, per_replicate_df, group_reps, group_name) {
    
    # Filter to group replicates
    metrics_group <- metrics_df %>% filter(repl %in% group_reps)
    curves_group <- curves_df %>% filter(repl %in% group_reps)
    per_replicate_df_group <- per_replicate_df %>% filter(repl_num %in% group_reps)
    
    # 1) Average RMSE by H_fit
    rmse_by_H <- per_replicate_df_group %>%
        group_by(H_fit) %>%
        summarise(
            mean_rmse = mean(diff_mean_vs_mu, na.rm = TRUE),
            median_rmse = median(diff_mean_vs_mu, na.rm = TRUE),
            .groups = "drop"
        )
    
    # 2) Across-replicate variance by H_fit
    curves_mean <- curves_group %>%
        select(H_fit, repl, post_mean_t) %>%
        unnest(post_mean_t) %>%
        group_by(H_fit, repl) %>%
        mutate(year_idx = row_number()) %>%
        ungroup()
    
    var_by_H <- curves_mean %>%
        group_by(H_fit, year_idx) %>%
        summarize(var_across_repl = var(post_mean_t, na.rm = TRUE), .groups = "drop") %>%
        group_by(H_fit) %>%
        summarize(
            mean_var = mean(var_across_repl, na.rm = TRUE),
            max_var = max(var_across_repl, na.rm = TRUE),
            .groups = "drop"
        )
   
    # Combine RMSE and variance
    results <- rmse_by_H %>%
        left_join(var_by_H, by = "H_fit") %>%
        mutate(group = group_name)
    
    results
}

#' Detect elbow using PCA method (perpendicular distance)
find_elbow_pca <- function(x, y) {
    if (length(x) < 3) return(x[1])
    
    # Find first decreasing point
    first_dec <- which(diff(y) < 0)[1]
    if (is.na(first_dec)) return(x[1])
    
    # Find overall minimum
    min_idx <- which.min(y)
    
    # Search range: first decrease to min+1
    search_range <- first_dec:min(min_idx + 1, length(y))
    if (length(search_range) < 2) return(x[first_dec])
    
    x_valid <- x[search_range]
    y_valid <- y[search_range]
    
    if (length(x_valid) < 2) return(x_valid[1])
    
    # Normalize to [0, 1]
    x_norm <- (x_valid - min(x_valid)) / (max(x_valid) - min(x_valid))
    y_norm <- (y_valid - min(y_valid)) / (max(y_valid) - min(y_valid))
    
    # Line from first to last point
    x1 <- x_norm[1]
    y1 <- y_norm[1]
    x2 <- x_norm[length(x_norm)]
    y2 <- y_norm[length(y_norm)]
    
    # Perpendicular distance from each point to line
    distances <- abs((y2 - y1) * x_norm - (x2 - x1) * y_norm + x2 * y1 - y2 * x1) /
                 sqrt((y2 - y1)^2 + (x2 - x1)^2)
    
    x_valid[which.max(distances)]
}

#' Detect elbow using slope change method
find_elbow_slope <- function(x, y) {
    if (length(x) < 3) return(x[1])
    
    # Find first decreasing point
    first_dec <- which(diff(y) < 0)[1]
    if (is.na(first_dec)) return(x[1])
    
    # Find overall minimum
    min_idx <- which.min(y)
    
    # Search range: first decrease to min+1
    search_range <- first_dec:min(min_idx + 1, length(y))
    if (length(search_range) < 3) return(x[first_dec])
    
    x_valid <- x[search_range]
    y_valid <- y[search_range]
    
    # Normalize
    x_norm <- (x_valid - min(x_valid)) / (max(x_valid) - min(x_valid))
    y_norm <- (y_valid - min(y_valid)) / (max(y_valid) - min(y_valid))
    
    # Calculate slopes between consecutive points
    slopes <- diff(y_norm) / diff(x_norm)
    
    # Find largest change in slope
    if (length(slopes) < 2) return(x_valid[1])
    slope_changes <- abs(diff(slopes))
    
    x_valid[which.max(slope_changes) + 1]
}

#' Detect elbow using angle method (closest to 90 degrees)
find_elbow_angle <- function(x, y) {
    if (length(x) < 3) return(x[1])
    
    # Find first decreasing point
    first_dec <- which(diff(y) < 0)[1]
    if (is.na(first_dec)) return(x[1])
    
    # Find overall minimum
    min_idx <- which.min(y)
    
    # Search range: first decrease to min+1
    search_range <- first_dec:min(min_idx + 1, length(y))
    if (length(search_range) < 3) return(x[first_dec])
    
    x_valid <- x[search_range]
    y_valid <- y[search_range]
    
    # Normalize
    x_norm <- (x_valid - min(x_valid)) / (max(x_valid) - min(x_valid))
    y_norm <- (y_valid - min(y_valid)) / (max(y_valid) - min(y_valid))
    
    # Calculate angles at each point
    angles <- numeric(length(x_norm) - 2)
    for (i in 2:(length(x_norm) - 1)) {
        v1 <- c(x_norm[i] - x_norm[i-1], y_norm[i] - y_norm[i-1])
        v2 <- c(x_norm[i+1] - x_norm[i], y_norm[i+1] - y_norm[i])
        
        dot_prod <- sum(v1 * v2)
        norm_v1 <- sqrt(sum(v1^2))
        norm_v2 <- sqrt(sum(v2^2))
        
        if (norm_v1 > 0 && norm_v2 > 0) {
            angle_rad <- acos(pmax(-1, pmin(1, dot_prod / (norm_v1 * norm_v2))))
            angles[i-1] <- angle_rad * 180 / pi
        }
    }
    
    # Find angle closest to 90 degrees
    x_valid[which.min(abs(angles - 90)) + 1]
}

#' Identify elbows for all H_true, metrics, and groups using specified method
identify_elbows <- function(results_df, method = "pca") {
    
    # Select elbow detection function
    elbow_func <- switch(method,
        "pca" = find_elbow_pca,
        "slope" = find_elbow_slope,
        "angle" = find_elbow_angle,
        find_elbow_pca  # default to PCA
    )
    
    elbows <- list()
    
    for (h_true in unique(results_df$H_true)) {
        for (metric_name in c("mean_rmse", "max_var")) {
            subset <- results_df %>% filter(H_true == h_true)
            
            # Group A
            ga_data <- subset %>% filter(group == "Group A") %>% arrange(H_fit)
            if (nrow(ga_data) > 0 && metric_name %in% names(ga_data)) {
                elbow_a <- elbow_func(ga_data$H_fit, ga_data[[metric_name]])
                elbows[[length(elbows) + 1]] <- data.frame(
                    H_true = h_true,
                    metric = ifelse(metric_name == "mean_rmse", "RMSE", "Max Variance"),
                    group = "Group A",
                    elbow = elbow_a
                )
            }
            
            # Group B
            gb_data <- subset %>% filter(group == "Group B") %>% arrange(H_fit)
            if (nrow(gb_data) > 0 && metric_name %in% names(gb_data)) {
                elbow_b <- elbow_func(gb_data$H_fit, gb_data[[metric_name]])
                elbows[[length(elbows) + 1]] <- data.frame(
                    H_true = h_true,
                    metric = ifelse(metric_name == "mean_rmse", "RMSE", "Max Variance"),
                    group = "Group B",
                    elbow = elbow_b
                )
            }
        }
    }
    
    bind_rows(elbows)
}

# ============================================================================
# LOAD DATA FOR ALL H_TRUE VALUES
# ============================================================================

cat("Loading data for all H_true values...\n")

all_data <- list()
for (H_true in H_TRUE_VALUES) {
    cat("  Loading H_true =", H_true, "... ")
    data <- load_results(H_true, LOCATION)
    if (!is.null(data)) {
        all_data[[as.character(H_true)]] <- data
        cat("✓\n")
    } else {
        cat("✗ (skipped)\n")
    }
}

if (length(all_data) == 0) {
    stop("No data loaded. Check directory paths and file locations.")
}

cat("\nLoaded data for", length(all_data), "H_true values\n\n")

# ============================================================================
# SPLIT-GROUP ANALYSIS FOR EACH H_TRUE
# ============================================================================

cat("Performing split-group analysis...\n\n")

split_results_list <- list()

for (H_true_str in names(all_data)) {
    H_true <- as.integer(H_true_str)
    cat("H_true =", H_true, "\n")
    
    data <- all_data[[H_true_str]]
    metrics <- data$metrics
    curves <- data$curves
    per_replicate <- data$per_replicate

    # Get all replicates
    all_reps <- sort(unique(metrics$repl))
    n_reps <- length(all_reps)
    
    cat("  Total replicates:", n_reps, "\n")
    
    # Split into groups
    groups <- split_replicates(all_reps, N_GROUP_A, SPLIT_SEED)
    
    cat("  Group A:", length(groups$A), "replicates\n")
    cat("  Group B:", length(groups$B), "replicates\n")
    
    # Analyze both groups
    results_A <- analyze_split_group(metrics, curves, per_replicate, groups$A, "Group A")
    results_B <- analyze_split_group(metrics, curves, per_replicate, groups$B, "Group B")
    
    # Combine and add H_true
    combined <- bind_rows(results_A, results_B) %>%
        mutate(H_true = H_true)
    
    split_results_list[[H_true_str]] <- combined
    cat("  ✓ Analysis complete\n\n")
}

# Combine all results
split_results_all <- bind_rows(split_results_list)

cat("Split-group analysis complete for all H_true values\n\n")

# ============================================================================
# CREATE PLOTS
# ============================================================================

cat("Creating split-group plots...\n\n")

# Prepare data for plotting
plot_data_rmse <- split_results_all %>%
    select(H_true, H_fit, group, mean_rmse) %>%
    rename(value = mean_rmse) %>%
    mutate(metric = "RMSE")

plot_data_var <- split_results_all %>%
    select(H_true, H_fit, group, max_var) %>%
    rename(value = max_var) %>%
    mutate(metric = "Max Variance")

plot_data <- bind_rows(plot_data_rmse, plot_data_var) %>%
    mutate(
        H_true_label = paste0("H = ", H_true),
        metric = factor(metric, levels = c("RMSE", "Max Variance"))
    )

# Identify elbows using all three methods
elbows_pca <- identify_elbows(split_results_all, method = "pca") %>% mutate(method = "PCA")
elbows_slope <- identify_elbows(split_results_all, method = "slope") %>% mutate(method = "Slope")
elbows_angle <- identify_elbows(split_results_all, method = "angle") %>% mutate(method = "Angle")

# Combine all elbows
all_elbows <- bind_rows(elbows_pca, elbows_slope, elbows_angle) %>%
    mutate(
        H_true_label = paste0("H = ", H_true),
        elbow_label = paste(method, group, sep = " - ")
    ) %>%
    distinct(H_true, metric, group, method, elbow, H_true_label, elbow_label, .keep_all = TRUE)

# Create single combined plot with all three methods
p_combined <- ggplot(plot_data, aes(x = H_fit, y = value, color = group, group = group)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    # Add vertical lines for all elbows, colored by method
    geom_vline(data = all_elbows,
               aes(xintercept = elbow, linetype = group, alpha = method),
               linewidth = 0.7) +
    facet_grid(metric ~ H_true_label, scales = "free_y") +
    scale_x_continuous(breaks = seq(2, 20, by = 1)) +
    scale_color_manual(
        values = c("Group A" = "#E41A1C", "Group B" = "#377EB8"),
        name = "Replicate Group"
    ) +
    scale_linetype_manual(
        values = c("Group A" = "dashed", "Group B" = "dotted"),
        name = "Group"
    ) +
    scale_alpha_manual(
        values = c("PCA" = 1.0, "Slope" = 0.6, "Angle" = 0.3),
        name = "Elbow Method"
    ) +
    labs(
        title = if (SHOW_TITLES) paste0("Split-Group Elbow Analysis: ", LOCATION, " (All Methods)") else NULL,
        subtitle = if (SHOW_TITLES) paste0(
            "Each group has ", N_GROUP_A, " randomly sampled replicates | ",
            "Vertical lines show elbows from PCA (dark), Slope (medium), and Angle (light) methods"
        ) else NULL,
        x = "h",
        y = NULL
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 9),
        strip.text = element_text(face = "bold", size = 9),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom",
        legend.box = "horizontal",
        panel.spacing = unit(0.8, "lines")
    )

# Save combined plot
plot_file_combined <- file.path(OUTPUT_DIR, paste0("split_group_elbow_", LOCATION, "_all_methods.png"))
ggsave(plot_file_combined, p_combined, width = 12, height = 14, dpi = 300)
cat("✓ Saved combined plot:", plot_file_combined, "\n\n")

print(p_combined)

# Also save individual method plots
make_method_plot <- function(method_label, elbows_df) {
    elbows_df <- elbows_df %>%
        mutate(H_true_label = paste0("H = ", H_true)) %>%
        distinct(H_true, metric, group, elbow, H_true_label, .keep_all = TRUE)
    
    ggplot(plot_data, aes(x = H_fit, y = value, color = group, group = group)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        geom_vline(data = elbows_df,
                   aes(xintercept = elbow, linetype = group),
                   linewidth = 0.7, alpha = 0.7) +
        facet_grid(metric ~ H_true_label, scales = "free_y") +
        scale_x_continuous(breaks = seq(2, 20, by = 1)) +
        scale_color_manual(
            values = c("Group A" = "#E41A1C", "Group B" = "#377EB8"),
            name = "Replicate Group"
        ) +
        scale_linetype_manual(
            values = c("Group A" = "dashed", "Group B" = "dotted"),
            name = "Group"
        ) +
        labs(
            title = if (SHOW_TITLES) paste0("Split-Group Elbow Analysis: ", LOCATION, " (", method_label, ")") else NULL,
            subtitle = if (SHOW_TITLES) paste0("Each group has ", N_GROUP_A, " randomly sampled replicates") else NULL,
            x = "h",
            y = NULL
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 10),
            strip.text = element_text(face = "bold", size = 9),
            strip.text.y = element_text(angle = 0),
            legend.position = "bottom",
            legend.box = "horizontal",
            panel.spacing = unit(0.8, "lines")
        )
}

p_pca <- make_method_plot("PCA (Perpendicular Distance)", elbows_pca)
p_slope <- make_method_plot("Slope Change", elbows_slope)
p_angle <- make_method_plot("Angle (90°)", elbows_angle)

plot_file_pca <- file.path(OUTPUT_DIR, paste0("split_group_elbow_", LOCATION, "_pca.png"))
plot_file_slope <- file.path(OUTPUT_DIR, paste0("split_group_elbow_", LOCATION, "_slope.png"))
plot_file_angle <- file.path(OUTPUT_DIR, paste0("split_group_elbow_", LOCATION, "_angle.png"))

ggsave(plot_file_pca, p_pca, width = 12, height = 14, dpi = 300)
ggsave(plot_file_slope, p_slope, width = 12, height = 14, dpi = 300)
ggsave(plot_file_angle, p_angle, width = 12, height = 14, dpi = 300)

cat("✓ Saved PCA plot:", plot_file_pca, "\n")
cat("✓ Saved Slope plot:", plot_file_slope, "\n")
cat("✓ Saved Angle plot:", plot_file_angle, "\n\n")

# ============================================================================
# MONTE CARLO SPLIT ANALYSIS
# ============================================================================

cat(strrep("=", 70), "\n")
cat("MONTE CARLO ANALYSIS: Robustness Across Random Splits\n")
cat(strrep("=", 70), "\n\n")

N_SPLITS <- 100
cat("Running", N_SPLITS, "random splits...\n\n")

# Store results from all splits
mc_results <- list()

for (split_i in 1:N_SPLITS) {
    if (split_i %% 10 == 0) cat("  Split", split_i, "/", N_SPLITS, "\n")
    
    # Use different seed for each split
    split_seed <- SPLIT_SEED + split_i
    
    # Analyze each H_true with this split
    for (H_true_str in names(all_data)) {
        H_true <- as.integer(H_true_str)
        data <- all_data[[H_true_str]]
        
        # Get all replicates and split
        all_reps <- sort(unique(data$metrics$repl))
        groups <- split_replicates(all_reps, N_GROUP_A, split_seed)
        
        # Analyze both groups
        results_A <- analyze_split_group(data$metrics, data$curves, data$per_replicate, 
                                         groups$A, "Group A")
        results_B <- analyze_split_group(data$metrics, data$curves, data$per_replicate, 
                                         groups$B, "Group B")
        
        combined <- bind_rows(results_A, results_B) %>% mutate(H_true = H_true)
        
        # Compute elbows for all three methods
        for (method_name in c("pca", "slope", "angle")) {
            elbows <- identify_elbows(combined, method = method_name)
            
            # Store results
            for (i in 1:nrow(elbows)) {
                mc_results[[length(mc_results) + 1]] <- data.frame(
                    split_id = split_i,
                    seed = split_seed,
                    H_true = elbows$H_true[i],
                    metric = elbows$metric[i],
                    group = elbows$group[i],
                    method = method_name,
                    elbow = elbows$elbow[i]
                )
            }
        }
    }
}

mc_df <- bind_rows(mc_results)

cat("\n✓ Completed", N_SPLITS, "splits\n\n")

# ============================================================================
# ANALYZE AGREEMENT BETWEEN GROUPS
# ============================================================================

cat("Analyzing agreement between groups across splits...\n\n")

# For each split, method, H_true, and metric, check if Group A and Group B agree
agreement_analysis <- mc_df %>%
    filter(group != "Overall") %>%
    select(split_id, H_true, metric, method, group, elbow) %>%
    pivot_wider(names_from = group, values_from = elbow, names_prefix = "elbow_") %>%
    mutate(agree = `elbow_Group A` == `elbow_Group B`) %>%
    group_by(H_true, metric, method) %>%
    summarize(
        n_splits = n(),
        n_agree = sum(agree, na.rm = TRUE),
        agreement_rate = mean(agree, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(method, H_true, metric)

cat("Agreement rates (% of splits where Group A and Group B chose same elbow):\n\n")
print(agreement_analysis, n = Inf)

# Save agreement analysis
agreement_file <- file.path(OUTPUT_DIR, paste0("mc_agreement_", LOCATION, ".csv"))
write.csv(agreement_analysis, agreement_file, row.names = FALSE)
cat("\n✓ Saved agreement analysis:", agreement_file, "\n\n")

# ============================================================================
# ANALYZE DISTRIBUTION OF ELBOW LOCATIONS
# ============================================================================

cat("Analyzing distribution of elbow locations across splits...\n\n")

# For each method, H_true, metric, and group, get distribution of elbows
elbow_distributions <- mc_df %>%
    group_by(method, H_true, metric, group, elbow) %>%
    summarize(count = n(), .groups = "drop") %>%
    group_by(method, H_true, metric, group) %>%
    mutate(
        total = sum(count),
        proportion = count / total
    ) %>%
    arrange(method, H_true, metric, group, desc(count))

# Summary statistics for each combination
elbow_summary <- mc_df %>%
    group_by(method, H_true, metric, group) %>%
    summarize(
        mean_elbow = mean(elbow, na.rm = TRUE),
        median_elbow = median(elbow, na.rm = TRUE),
        sd_elbow = sd(elbow, na.rm = TRUE),
        min_elbow = min(elbow, na.rm = TRUE),
        max_elbow = max(elbow, na.rm = TRUE),
        n_unique = n_distinct(elbow),
        .groups = "drop"
    ) %>%
    arrange(method, H_true, metric, group)

cat("Summary statistics for elbow locations:\n\n")
print(elbow_summary, n = Inf)

# Save distributions
dist_file <- file.path(OUTPUT_DIR, paste0("mc_elbow_distributions_", LOCATION, ".csv"))
write.csv(elbow_distributions, dist_file, row.names = FALSE)
cat("\n✓ Saved elbow distributions:", dist_file, "\n\n")

summary_file <- file.path(OUTPUT_DIR, paste0("mc_elbow_summary_", LOCATION, ".csv"))
write.csv(elbow_summary, summary_file, row.names = FALSE)
cat("✓ Saved elbow summary:", summary_file, "\n\n")

# ============================================================================
# VISUALIZE MONTE CARLO RESULTS
# ============================================================================

cat("Creating Monte Carlo visualization...\n\n")

# Create histogram plots for elbow distributions
mc_plot_data <- mc_df %>%
    filter(group != "Overall") %>%
    mutate(
        H_true_label = paste0("H = ", H_true),
        method_label = case_when(
            method == "pca" ~ "PCA",
            method == "slope" ~ "Slope",
            method == "angle" ~ "Angle"
        )
    ) %>% filter(method == "pca")  # Focus on PCA for visualization

p_mc_dist <- ggplot(mc_plot_data, aes(x = elbow, fill = group)) +
    geom_histogram(position = "dodge", bins = 15, alpha = 0.7) +
    facet_grid(metric ~ H_true_label, scales = "free_y")+# + method_label, scales = "free_y") +
    scale_fill_manual(
        values = c("Group A" = "#E41A1C", "Group B" = "#377EB8"),
        name = "Group"
    ) +
    scale_x_continuous(breaks = seq(2, 20, by = 1))+
    labs(
        title = if (SHOW_TITLES) paste0("Monte Carlo Analysis: Elbow Distributions (", N_SPLITS, " Splits) - ", LOCATION) else NULL,
        subtitle = if (SHOW_TITLES) "Distribution of elbow locations across random 50/50 splits" else NULL,
        x = "Elbow location (h)",
        y = "Count"
    ) +
    theme_bw() +
    theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10),
        strip.text = element_text(face = "bold", size = 7),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom",
        panel.spacing = unit(0.5, "lines")
    )

mc_dist_file <- file.path(OUTPUT_DIR, paste0("mc_elbow_distributions_", LOCATION, ".png"))
ggsave(mc_dist_file, p_mc_dist, width = 16, height = 8, dpi = 300)
cat("✓ Saved MC distribution plot:", mc_dist_file, "\n\n")

# Create agreement rate visualization
agreement_plot_data <- agreement_analysis %>%
    mutate(
        H_true_label = paste0("H = ", H_true),
        method_label = case_when(
            method == "pca" ~ "PCA",
            method == "slope" ~ "Slope",
            method == "angle" ~ "Angle"
        )
    ) %>% filter(method == "pca")  # Focus on PCA for visualization

p_agreement <- ggplot(agreement_plot_data, aes(x = H_true_label, y = agreement_rate * 100)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,fill='navy') +
    geom_text(aes(label = sprintf("%.0f%%", agreement_rate * 100)),
              position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
    facet_wrap(~ metric, ncol = 2) +
    labs(
        title = if (SHOW_TITLES) paste0("Group Agreement Rates Across ", N_SPLITS, " Random Splits - ", LOCATION) else NULL,
        subtitle = if (SHOW_TITLES) "Percentage of splits where Group A and Group B selected the same elbow" else NULL,
        x = "True H",
        y = "Agreement Rate (%)"
    ) +
    ylim(0, 105) +
    theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "bottom"
    )

agreement_plot_file <- file.path(OUTPUT_DIR, paste0("mc_agreement_rates_", LOCATION, ".png"))
ggsave(agreement_plot_file, p_agreement, width = 10, height = 6, dpi = 300)
cat("✓ Saved agreement rate plot:", agreement_plot_file, "\n\n")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat(strrep("=", 70), "\n")
cat("SUMMARY: Elbow Detection Methods Comparison\n")
cat(strrep("=", 70), "\n\n")

# Compare elbows across methods
method_comparison <- data.frame(
    H_true = elbows_pca$H_true,
    group = elbows_pca$group,
    metric = elbows_pca$metric,
    elbow_pca = elbows_pca$elbow,
    elbow_slope = elbows_slope$elbow,
    elbow_angle = elbows_angle$elbow
) %>%
    mutate(
        agree_all = (elbow_pca == elbow_slope) & (elbow_slope == elbow_angle),
        agree_pca_slope = elbow_pca == elbow_slope,
        agree_pca_angle = elbow_pca == elbow_angle,
        agree_slope_angle = elbow_slope == elbow_angle
    )

cat("Agreement between methods:\n\n")
cat("All three methods agree:", sum(method_comparison$agree_all), "/", nrow(method_comparison), 
    "cases (", round(100 * mean(method_comparison$agree_all), 1), "%)\n", sep = "")
cat("PCA and Slope agree:", sum(method_comparison$agree_pca_slope), "/", nrow(method_comparison), 
    "(", round(100 * mean(method_comparison$agree_pca_slope), 1), "%)\n", sep = "")
cat("PCA and Angle agree:", sum(method_comparison$agree_pca_angle), "/", nrow(method_comparison), 
    "(", round(100 * mean(method_comparison$agree_pca_angle), 1), "%)\n", sep = "")
cat("Slope and Angle agree:", sum(method_comparison$agree_slope_angle), "/", nrow(method_comparison), 
    "(", round(100 * mean(method_comparison$agree_slope_angle), 1), "%)\n\n", sep = "")

# Save method comparison
method_comp_file <- file.path(OUTPUT_DIR, paste0("elbow_methods_comparison_", LOCATION, ".csv"))
write.csv(method_comparison, method_comp_file, row.names = FALSE)
cat("✓ Saved method comparison:", method_comp_file, "\n\n")

cat(strrep("=", 70), "\n")
cat("✓ PAPER PLOTS COMPLETE\n")
cat(strrep("=", 70), "\n\n")
