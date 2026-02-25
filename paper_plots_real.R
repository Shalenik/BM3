# =========================
# CONSOLIDATED ELBOW PLOTS FOR REAL DATA (GOM & NY, SSP126/245/585)
# =========================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

real_scenarios <- c("ssp126_real", "ssp245_real", "ssp585_real")
real_regions <- c("Gulf_of_Mexico", "New_York")
H_range <- 2:7

get_real_elbow_data <- function(region, scenario) {
    root <- file.path("results", scenario, region)
    out <- list()
    for (H in H_range) {
        H_dir <- file.path(root, paste0("H_", H))
        rep_dirs <- list.files(H_dir, pattern = "^rep_", full.names = TRUE)
        if (length(rep_dirs) == 0) next
        # Use first replicate with a fit file
        fitfile <- list.files(rep_dirs[1], pattern = "^fit_.*\\.rds$", full.names = TRUE)
        if (length(fitfile) == 0) next
        fit <- readRDS(fitfile[1])
        # Extract posterior means for each time point
        ysim <- fit$fit$BM3$ysim
        n_iter <- nrow(ysim)
        burnin_idx <- ceiling(n_iter / 2)
        ysim_post <- ysim[burnin_idx:n_iter, 2:ncol(ysim)]
        # Variance across posterior draws (mean and max over time)
        var_post <- apply(ysim_post, 2, var, na.rm = TRUE)
        out[[length(out)+1]] <- data.frame(
            H = H,
            mean_var = mean(var_post, na.rm = TRUE),
            max_var = max(var_post, na.rm = TRUE)
        )
    }
    bind_rows(out)
}

elbow_plot <- function(df, region, scenario) {
    if (nrow(df) < 3) {
        return(ggplot() + annotate("text", x = 1, y = 0, label = "Insufficient data", size = 4) + theme_void())
    }
    # Perpendicular elbow for mean_var
    elbow_mean <- find_elbow(df$H, df$mean_var, minimize = TRUE)
    elbow_max <- find_elbow(df$H, df$max_var, minimize = TRUE)
    p <- ggplot(df, aes(x = H)) +
        geom_line(aes(y = mean_var, color = "Mean Variance"), size = 1) +
        geom_line(aes(y = max_var, color = "Max Variance"), size = 1, linetype = "dashed") +
        geom_vline(xintercept = elbow_mean, color = "#2ca02c", linetype = "dotted", size = 0.7) +
        geom_vline(xintercept = elbow_max, color = "#ff7f0e", linetype = "dotted", size = 0.7) +
        scale_color_manual(values = c("Mean Variance" = "#2ca02c", "Max Variance" = "#ff7f0e")) +
        labs(
            title = paste(region, scenario),
            x = "H",
            y = "Variance",
            color = NULL
        ) +
        theme_bw(base_size = 10) +
        theme(plot.title = element_text(size = 9, face = "bold"), legend.position = "bottom")
    return(p)
}

cat("\nGenerating consolidated elbow plot grid for real data...\n")
elbow_grid <- list()
for (i in seq_along(real_regions)) {
    for (j in seq_along(real_scenarios)) {
        region <- real_regions[i]
        scenario <- real_scenarios[j]
        df <- get_real_elbow_data(region, scenario)
        elbow_grid[[length(elbow_grid)+1]] <- elbow_plot(df, region, scenario)
    }
}

elbow_grid_plot <- (
    elbow_grid[[1]] | elbow_grid[[2]] | elbow_grid[[3]]) /
    (elbow_grid[[4]] | elbow_grid[[5]] | elbow_grid[[6]])

elbow_grid_plot <- elbow_grid_plot +
    plot_annotation(
        title = "Elbow Plots for Real Data (Mean/Max Variance by H)",
        subtitle = "Rows: Gulf of Mexico, New York | Columns: SSP126, SSP245, SSP585\nDotted lines: Perpendicular elbow for each variance curve",
        theme = theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 10, hjust = 0.5)
        )
    )

ggsave("results/ElbowGrid_RealData.png", elbow_grid_plot, width = 14, height = 7, dpi = 300)
cat("Saved: results/ElbowGrid_RealData.png\n")
#!/usr/bin/env Rscript
# Paper plots: Comprehensive 6×3 grid comparing NY and GOM across H_true values

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(ggplot2)
    library(patchwork)
})

# ============================================================================
# CONFIGURATION
# ============================================================================

scenario <- "585"
H_true_values <- c(2, 4, 7)
locations <- c("NY", "GOM")
show_overall_title <- FALSE  # Set to TRUE to show overall title and subtitle

# Color palette for H_fit values
h_colors <- c(
    "2" = "#1f77b4",
    "3" = "#ff7f0e", 
    "4" = "#2ca02c",
    "5" = "#d62728",
    "6" = "#9467bd",
    "7" = "#8c564b",
    "8" = "#e377c2",
    "10" = "#7f7f7f",
    "12" = "#bcbd22",
    "16" = "#17becf",
    "20" = "#aec7e8"
)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Load data for a specific location and H_true
load_htrue_data <- function(loc, h_true) {
    base_dir <- paste0("/work/users/s/h/shaleni/GRL_private/results/sim_posterior_", loc)
    root <- file.path(base_dir, scenario, paste0("Htrue_", h_true))
    
    # Load curves data
    curves_file <- file.path(root, "timewise_posterior_summaries.rds")
    if (!file.exists(curves_file)) {
        warning("Curves file not found: ", curves_file)
        return(NULL)
    }
    curves_data <- readRDS(curves_file)
    
    # Load averaged Y matrix
    avg_Y_file <- file.path(root, "averaged_Y_matrix.rds")
    if (!file.exists(avg_Y_file)) {
        warning("Averaged Y matrix not found: ", avg_Y_file)
        return(NULL)
    }
    avg_Y_data <- readRDS(avg_Y_file)
    
    # Load per-replicate stats for determining "disagreement"
    repl_stats_file <- file.path(root, "per_replicate_posterior_stats.rds")
    if (!file.exists(repl_stats_file)) {
        warning("Replicate stats not found: ", repl_stats_file)
        return(NULL)
    }
    repl_stats <- readRDS(repl_stats_file)
    
    # Load sim for true mean
    rep_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
    rep_dirs <- rep_dirs[grepl("rep_[0-9]+$", rep_dirs)]
    if (length(rep_dirs) == 0) return(NULL)
    
    sim_file <- file.path(rep_dirs[1], "sim.rds")
    if (!file.exists(sim_file)) return(NULL)
    sim <- readRDS(sim_file)
    
    list(
        curves = curves_data,
        Y_avg = avg_Y_data$Y_avg,
        years = avg_Y_data$years,
        true_mean = avg_Y_data$true_mean,
        repl_stats = repl_stats,
        cutoff_year = sim$years[max(which(!is.na(sim$y_masked)))],
        loc = loc,
        h_true = h_true
    )
}

#' Find elbow point using perpendicular distance method
find_elbow <- function(x, y, minimize = TRUE) {
    if (length(x) < 3) return(x[1])
    
    # If maximizing, flip the curve
    if (!minimize) {
        y <- -y
    }
    
    # Skip any initial increase - find where decrease starts
    start_idx <- 1
    for (i in 1:(length(y)-1)) {
        if (y[i+1] < y[i]) {
            start_idx <- i
            break
        }
    }
    
    # Find the minimum from the decreasing region onwards
    y_subset <- y[start_idx:length(y)]
    min_idx_rel <- which.min(y_subset)
    min_idx <- start_idx + min_idx_rel
    
    # Only analyze points from start of decrease to minimum
    x_mono <- x[start_idx:min_idx]
    y_mono <- y[start_idx:min_idx]
    
    if (length(x_mono) < 3) return(x_mono[1])
    
    # Normalize x and y to 0-1 scale for fair distance calculation
    x_norm <- (x_mono - min(x_mono)) / (max(x_mono) - min(x_mono))
    y_norm <- (y_mono - min(y_mono)) / (max(y_mono) - min(y_mono))
    
    # Line from first to last point in the decreasing region
    x1 <- x_norm[1]
    y1 <- y_norm[1]
    x2 <- x_norm[length(x_norm)]
    y2 <- y_norm[length(y_norm)]
    
    # Calculate perpendicular distance from each point to the line
    distances <- abs((y2 - y1) * x_norm - (x2 - x1) * y_norm + x2 * y1 - y2 * x1) /
                 sqrt((y2 - y1)^2 + (x2 - x1)^2)
    
    # The elbow is the point with maximum perpendicular distance
    max_dist <- max(distances)
    candidates <- which(distances >= max_dist)
    elbow_idx <- candidates[1]  # Take the first (earliest H) among close candidates
    
    return(x_mono[elbow_idx])
}

#' Determine which H_fit values "disagree" based on RMSE from true mean
#' Returns vector of H_fit values that should be dashed
get_disagreeing_hfits <- function(curves_avg, true_mean, years, threshold_percentile = 0.75) {
    # Calculate RMSE for each H_fit between averaged posterior mean and true mean
    rmse_by_h <- curves_avg %>%
        mutate(true_mean_val = rep(true_mean, length(unique(H_fit)))) %>%
        group_by(H_fit) %>%
        summarize(
            rmse = sqrt(mean((post_mean - true_mean_val)^2, na.rm = TRUE)),
            .groups = "drop"
        )
    
    # Determine threshold: H_fits with RMSE above the threshold_percentile are "disagreeing"
    rmse_threshold <- quantile(rmse_by_h$rmse, threshold_percentile, na.rm = TRUE)
    
    disagreeing <- rmse_by_h %>%
        filter(rmse > rmse_threshold) %>%
        pull(H_fit)
    
    return(disagreeing)
}

#' Create Column 1: Convergence plot with Y_avg, posterior means, and true mean
create_convergence_plot <- function(data, show_title = TRUE, show_legend = FALSE) {
    # Prepare curves data
    curves_long <- data$curves %>%
        mutate(
            year = map(post_median_t, ~ data$years),
            post = post_median_t
        ) %>%
        select(H_fit, repl, year, post) %>%
        unnest(cols = c(year, post))
    
    curves_avg <- curves_long %>%
        group_by(H_fit, year) %>%
        summarize(post_mean = mean(post, na.rm = TRUE), .groups = "drop")
    
    # Determine disagreeing H_fit values based on RMSE from true mean
    disagreeing <- get_disagreeing_hfits(curves_avg, data$true_mean, data$years)
    
    curves_avg <- curves_avg %>%
        mutate(
            linetype = ifelse(H_fit %in% disagreeing, "dashed", "solid"),
            H_fit_char = as.character(H_fit)
        )
    
    # Prepare Y_avg data
    Y_avg_long <- as.data.frame(data$Y_avg) %>%
        mutate(year = data$years) %>%
        pivot_longer(cols = -year, names_to = "H_fit_str", values_to = "Y_avg_value") %>%
        mutate(H_fit = as.numeric(gsub("V", "", H_fit_str)))
    
    # True mean data
    truth_df <- data.frame(year = data$years, true_mean = data$true_mean)
    
    # Shading for forecast period
    shade_df <- data.frame(
        xmin = data$cutoff_year,
        xmax = max(data$years),
        ymin = -Inf,
        ymax = Inf
    )
    
    # Create plot
    p <- ggplot() +
        geom_rect(
            data = shade_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightgray", alpha = 0.3
        ) +
        geom_line(
            data = Y_avg_long,
            aes(x = year, y = Y_avg_value, group = factor(H_fit)),
            color = "darkgrey", linewidth = 0.5, alpha = 0.7
        ) +
        geom_point(
            data = truth_df,
            aes(x = year, y = true_mean),
            color = "red", size = .8, alpha = 0.6
        ) +
        geom_line(
            data = curves_avg,
            aes(x = year, y = post_mean, color = H_fit_char, linetype = linetype),
            linewidth = 0.8
        ) +
        scale_color_manual(
            values = h_colors, 
            name = "H",
            breaks = as.character(sort(as.numeric(names(h_colors))))
        ) +
        scale_linetype_identity() +
        guides(color = guide_legend(nrow = 2)) +
        labs(
            title = paste0(data$loc, " (H=", data$h_true, ")"),
            x = NULL,
            y = "Temp Anomaly"
        ) +
        theme_bw(base_size = 10) +
        theme(
            plot.title = element_text(face = "bold", size = 10),
            legend.position = if (show_legend) "bottom" else "none",
            legend.text = element_text(size = 6),
            legend.key.size = unit(0.3, "cm"),
            legend.title = element_text(size = 7),
            axis.title.y = element_text(size = 8),
            axis.text = element_text(size = 7)
        )
    
    return(p)
}

#' Create Column 2: Normalized overlay of RMSE and variance metrics
create_dual_elbow_plot <- function(data, show_title = TRUE, show_legend = FALSE) {
    # Calculate average performance metrics
    perf_summary <- data$repl_stats %>%
        filter(H_true == data$h_true) %>%
        group_by(H_fit) %>%
        summarize(
            mean_rmse = mean(diff_mean_vs_mu_full, na.rm = TRUE),
            se_rmse = sd(diff_mean_vs_mu_full, na.rm = TRUE) / sqrt(n()),
            mean_var = mean(var_mean, na.rm = TRUE),
            se_var = sd(var_mean, na.rm = TRUE) / sqrt(n()),
            max_var = mean(var_max, na.rm = TRUE),
            se_max_var = sd(var_max, na.rm = TRUE) / sqrt(n()),
            .groups = "drop"
        ) %>%
        arrange(H_fit)
    
    # Normalize all metrics and errors to 0-1 scale
    perf_normalized <- perf_summary %>%
        mutate(
            rmse_norm = (mean_rmse - min(mean_rmse)) / (max(mean_rmse) - min(mean_rmse)),
            se_rmse_norm = se_rmse / (max(mean_rmse) - min(mean_rmse)),
            mean_var_norm = (mean_var - min(mean_var)) / (max(mean_var) - min(mean_var)),
            se_var_norm = se_var / (max(mean_var) - min(mean_var)),
            max_var_norm = (max_var - min(max_var)) / (max(max_var) - min(max_var)),
            se_max_var_norm = se_max_var / (max(max_var) - min(max_var))
        ) %>%
        pivot_longer(
            cols = c(rmse_norm, mean_var_norm, max_var_norm),
            names_to = "metric",
            values_to = "normalized_value"
        ) %>%
        mutate(
            se = case_when(
                metric == "rmse_norm" ~ se_rmse_norm,
                metric == "mean_var_norm" ~ se_var_norm,
                metric == "max_var_norm" ~ se_max_var_norm
            ),
            metric_label = case_when(
                metric == "rmse_norm" ~ "RMSE",
                metric == "mean_var_norm" ~ "Mean Variance",
                metric == "max_var_norm" ~ "Max Variance"
            )
        )
    
    # Find elbows for normalized data
    elbow_rmse <- find_elbow(perf_normalized$H_fit[perf_normalized$metric == "rmse_norm"], 
                             perf_normalized$normalized_value[perf_normalized$metric == "rmse_norm"], 
                             minimize = TRUE)
    elbow_mean_var <- find_elbow(perf_normalized$H_fit[perf_normalized$metric == "mean_var_norm"], 
                                 perf_normalized$normalized_value[perf_normalized$metric == "mean_var_norm"], 
                                 minimize = TRUE)
    elbow_max_var <- find_elbow(perf_normalized$H_fit[perf_normalized$metric == "max_var_norm"], 
                                perf_normalized$normalized_value[perf_normalized$metric == "max_var_norm"], 
                                minimize = TRUE)
    
    # Create overlay plot with error bars
    p <- ggplot(perf_normalized, aes(x = H_fit, y = normalized_value, color = metric_label)) +
        geom_line(linewidth = 1) +
        geom_errorbar(aes(ymin = normalized_value - se, ymax = normalized_value + se), width = 0.3, alpha = 0.5) +
        geom_vline(xintercept = elbow_rmse, color = "#1f77b4", 
                  linetype = "dotted", linewidth = 0.6, alpha = 0.7) +
        geom_vline(xintercept = elbow_mean_var, color = "#2ca02c", 
                  linetype = "dotted", linewidth = 0.6, alpha = 0.7) +
        geom_vline(xintercept = elbow_max_var, color = "#ff7f0e", 
                  linetype = "dotted", linewidth = 0.6, alpha = 0.7) +
        scale_color_manual(
            values = c("RMSE" = "#1f77b4", "Mean Variance" = "#2ca02c", "Max Variance" = "#ff7f0e"),
            name = NULL
        ) +
        labs(
            title = if (show_title) "Mean and Max Posterior Variance (avg over r), RMSE (avg over r)" else NULL,
            x = NULL,
            y = "Normalized Value"
        ) +
        theme_bw(base_size = 9) +
        theme(
            plot.title = element_text(size = 9),
            axis.title = element_text(size = 7),
            axis.text = element_text(size = 6),
            legend.position = if (show_legend) "bottom" else "none",
            legend.text = element_text(size = 6),
            legend.key.size = unit(0.3, "cm")
        )
    
    return(p)
}

#' Create Column 3: Mean across-replicate variance by H
create_variance_elbow_plot <- function(data, show_title = TRUE) {
    # Calculate across-replicate variance of posterior means
    # For each H_fit and time point, get variance of posterior means across replicates
    
    curves_long <- data$curves %>%
        mutate(
            year = map(post_median_t, ~ data$years),
            post_mean = post_mean_t  # posterior mean for this replicate
        ) %>%
        select(H_fit, repl, year, post_mean) %>%
        unnest(cols = c(year, post_mean))
    
    # Calculate variance across replicates for each H_fit and time point
    across_rep_var <- curves_long %>%
        group_by(H_fit, year) %>%
        summarize(
            var_post_mean = var(post_mean, na.rm = TRUE),
            .groups = "drop"
        )
    
    # Summarize: mean and SD of variance over time for each H_fit
    var_summary <- across_rep_var %>%
        group_by(H_fit) %>%
        summarize(
            mean_var = mean(var_post_mean, na.rm = TRUE),
            sd_var = sd(var_post_mean, na.rm = TRUE),
            se_var = sd_var / sqrt(n()),
            .groups = "drop"
        ) %>%
        arrange(H_fit)
    
    # Create plot
    p <- ggplot(var_summary, aes(x = H_fit, y = mean_var)) +
        geom_line(color = "#2ca02c", linewidth = 1) +
        geom_errorbar(aes(ymin = mean_var - se_var, ymax = mean_var + se_var),
                     width = 0.3, color = "#2ca02c", alpha = 0.6) +
        labs(
            title = if (show_title) "Variance of Posterior Means (avg over t)" else NULL,
            x = "H",
            y = "Variance"
        ) +
        theme_bw(base_size = 9) +
        theme(
            plot.title = element_text(size = 9),
            axis.title = element_text(size = 7),
            axis.text = element_text(size = 6)
        )
    
    return(p)
}

# ============================================================================
# MAIN PLOTTING ROUTINE
# ============================================================================

cat("Creating paper plots (6×3 grid)...\n\n")

# Storage for all plots
all_plots <- list()
plot_idx <- 1

# Loop through locations and H_true values
for (loc in locations) {
    cat("Processing", loc, "...\n")
    
    for (h in H_true_values) {
        cat("  H_true =", h, "...\n")
        
        # Load data
        data <- load_htrue_data(loc, h)
    
        if (is.null(data)) {
            cat("    Skipping (data not found)\n")
            # Create placeholder plots
            all_plots[[plot_idx]] <- ggplot() + 
                annotate("text", x = 0.5, y = 0.5, label = "Data not found", size = 5) +
                theme_void()
            all_plots[[plot_idx + 1]] <- all_plots[[plot_idx]]
            all_plots[[plot_idx + 2]] <- all_plots[[plot_idx]]
            plot_idx <- plot_idx + 3
            next
        }
        
        # Create the three columns for this row
        show_title <- (plot_idx <= 3)  # Only show titles in first row
        show_legend_col1 <- (plot_idx == 16)  # Only show legend at very bottom (last GOM row)
        show_legend_col2 <- (plot_idx == 16)  # Only show legend at very bottom (last GOM row)
        
        p1 <- create_convergence_plot(data, show_title = show_title, show_legend = show_legend_col1)
        p2 <- create_dual_elbow_plot(data, show_title = show_title, show_legend = show_legend_col2)
        p3 <- create_variance_elbow_plot(data, show_title = show_title)
        
        all_plots[[plot_idx]] <- p1
        all_plots[[plot_idx + 1]] <- p2
        all_plots[[plot_idx + 2]] <- p3
        
        plot_idx <- plot_idx + 3
    }
}

cat("\nAssembling final grid...\n")

# Ensure we have 18 plots (6 rows × 3 columns)
if (length(all_plots) != 18) {
    warning("Expected 18 plots but got ", length(all_plots))
}

# Create the 6×3 grid with proper layout
final_plot <- wrap_plots(all_plots, ncol = 3, nrow = 6)

if (show_overall_title) {
    final_plot <- final_plot +
        plot_annotation(
            title = "Model Selection Analysis: NY (top) vs GOM (bottom) for H = 2, 4, 7",
            subtitle = "Column 1: Black=climate data, Colors=posterior means (dashed=low selection rate), Red=true mean | Column 2: RMSE & Variance | Column 3: Variance metrics",
            theme = theme(
                plot.title = element_text(size = 14, face = "bold"),
                plot.subtitle = element_text(size = 10)
            )
        )
}

# Save the plot
print(final_plot)
output_file <- "/work/users/s/h/shaleni/GRL_private/results/paper_plots_6x3_grid.png"
cat("Saving to:", output_file, "\n")

ggsave(
    output_file,
    final_plot,
    width = 14,
    height = 16,
    dpi = 300,
    bg = "white"
)

cat("\nDone! Plot saved to:", output_file, "\n")

# =========================
# DEBUGGING: GOM H_true=4, H_fit=2,3 outlier analysis
# =========================

gom_h4_data <- load_htrue_data("GOM", 4)
if (!is.null(gom_h4_data)) {
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(purrr)
    set.seed(123)
    # Extract curves for H_fit=2 and H_fit=3
    curves_long <- gom_h4_data$curves %>%
        filter(H_fit %in% c(2, 3)) %>%
        mutate(year = map(post_median_t, ~ gom_h4_data$years),
               post = post_median_t) %>%
        select(H_fit, repl, year, post) %>%
        unnest(cols = c(year, post))

    # True mean for comparison
    truth_df <- data.frame(year = gom_h4_data$years, true_mean = gom_h4_data$true_mean)

    # Calculate RMSE for each replicate
    rmse_by_repl <- curves_long %>%
        left_join(truth_df, by = "year") %>%
        group_by(H_fit, repl) %>%
        summarize(rmse = sqrt(mean((post - true_mean)^2, na.rm = TRUE)), .groups = "drop")

    # Identify outliers (top 5% RMSE for each H_fit)
    outlier_df <- rmse_by_repl %>%
        group_by(H_fit) %>%
        mutate(outlier = rmse > quantile(rmse, 0.95)) %>%
        filter(outlier)

    print("Replicates with high RMSE (potential outliers):")
    print(outlier_df)

    # For each outlier, plot simulated data and a sample of posterior curves
    for (i in seq_len(nrow(outlier_df))) {
        hfit <- outlier_df$H_fit[i]
        repl <- outlier_df$repl[i]
        cat("\nPlotting outlier: H_fit=", hfit, ", repl=", repl, "\n")

        # Find the rep directory
        rep_dir <- file.path("/work/users/s/h/shaleni/GRL_private/results/sim_posterior_GOM/585/Htrue_4",
                            sprintf("rep_%03d", repl))
        sim_file <- file.path(rep_dir, "sim.rds")
        if (!file.exists(sim_file)) {
            cat("Sim file not found for repl", repl, "\n")
            next
        }
        sim <- readRDS(sim_file)

        # Check for valid y and y_masked
        n_years <- length(sim$years)
        valid_y <- !is.null(sim$y) && length(sim$y) == n_years
        valid_y_masked <- !is.null(sim$y_masked) && length(sim$y_masked) == n_years
        if (!valid_y || !valid_y_masked) {
            cat("Skipping repl", repl, ": sim$y or sim$y_masked missing or wrong length\n")
            next
        }
        sim_df <- data.frame(year = sim$years, y = sim$y, y_masked = sim$y_masked)

        # Posterior samples: try to get a sample of posterior curves for this replicate and H_fit
        # Find the corresponding posterior samples if available
        # If not available, just plot the median curve
        posterior_sample <- NULL
        if (!is.null(gom_h4_data$curves$post_samples)) {
            posterior_sample <- gom_h4_data$curves %>%
                filter(H_fit == hfit, repl == repl) %>%
                pull(post_samples)
            if (length(posterior_sample) > 0 && !is.null(posterior_sample[[1]])) {
                # Sample up to 20 posterior draws
                post_mat <- posterior_sample[[1]]
                n_draws <- min(20, ncol(post_mat))
                draw_idx <- sample(seq_len(ncol(post_mat)), n_draws)
                post_draws <- as.data.frame(post_mat[, draw_idx, drop = FALSE])
                post_draws$year <- gom_h4_data$years
                post_draws_long <- post_draws %>%
                    pivot_longer(-year, names_to = "draw", values_to = "value")
            } else {
                post_draws_long <- NULL
            }
        } else {
            post_draws_long <- NULL
        }

        # Median posterior curve for this replicate
        median_curve <- curves_long %>%
            filter(H_fit == hfit, repl == repl) %>%
            select(year, post) %>%
            mutate(year = as.numeric(year))

        # Plot grid: Simulated data (masked/unmasked), true mean, posterior curves
        p_sim <- ggplot(sim_df, aes(x = year)) +
            geom_line(aes(y = y), color = "black", size = 0.7, alpha = 0.7) +
            geom_point(aes(y = y_masked), color = "blue", size = 1, alpha = 0.7) +
            labs(title = paste0("Simulated data (H_fit=", hfit, ", repl=", repl, ")"),
                 y = "Simulated y") +
            theme_bw()

        p_post <- ggplot() +
            {if (!is.null(post_draws_long)) geom_line(data = post_draws_long, aes(x = year, y = value, group = draw), color = "grey", alpha = 0.5)} +
            geom_line(data = median_curve, aes(x = year, y = post), color = "darkgreen", size = 1) +
            geom_line(data = truth_df, aes(x = year, y = true_mean), color = "red", size = 1, linetype = "dashed") +
            labs(title = "Posterior curves (sample) & median", y = "Posterior mean") +
            theme_bw()

        gridExtra::grid.arrange(p_sim, p_post, ncol = 2)
    }
}


# =========================
