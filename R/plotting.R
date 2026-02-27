#' Plot Posterior Distribution
#'
#' Create a plot of posterior means and credible intervals.
#'
#' @param results List from `consolidate_results()`.
#' @param H Integer. Which H value to plot.
#' @param data Optional. Observations to overlay (default NULL).
#' @param title Optional plot title (default NULL).
#'
#' @return ggplot2 object.
#'
#' @examples
#' \dontrun{
#' plot_posterior(results, H = 2, title = "My Title")
#' }
#'
#' @export
plot_posterior <- function(results, H, data = NULL, title = NULL) {
    summaries <- results$timewise_summaries %>%
        dplyr::filter(H_fit == H)

    if (nrow(summaries) == 0) {
        stop(sprintf("No summaries found for H=%d", H))
    }

    p <- ggplot2::ggplot(summaries, ggplot2::aes(x = time, y = post_mean)) +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = post_ci_lower, ymax = post_ci_upper),
            alpha = 0.3, fill = "steelblue"
        ) +
        ggplot2::geom_line(color = "steelblue", linewidth = 0.7) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            x = "Time",
            y = "Predicted Value",
            subtitle = sprintf(
                "RMSE: %.4f",
                unique(summaries$rmse)
            )
        )
    
    if (!is.null(title)) {
        p <- p + ggplot2::labs(title = title)
    }

    if (!is.null(data)) {
        data_df <- tibble::tibble(time = 1:length(data), obs = data)
        p <- p + ggplot2::geom_point(
            data = data_df,
            ggplot2::aes(x = time, y = obs),
            color = "red", size = 1, alpha = 0.5
        )
    }

    return(p)
}


#' Plot Model Comparison
#'
#' Create a plot comparing model performance across H values.
#'
#' @param results List from `consolidate_results()`.
#' @param title Optional plot title (default NULL).
#'
#' @return ggplot2 object.
#'
#' @export
plot_comparison <- function(results, title = NULL) {
    comparison <- results$model_comparison %>%
        dplyr::select(H, RMSE) %>%
        dplyr::mutate(RMSE = scales::rescale(RMSE, to = c(0, 1)))

    p <- ggplot2::ggplot(comparison, ggplot2::aes(x = H, y = RMSE)) +
        ggplot2::geom_point(size = 3, color = "steelblue") +
        ggplot2::geom_line(color = "steelblue") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            x = "Number of Latent Factors (H)",
            y = "RMSE (normalized)"
        ) +
        ggplot2::theme(legend.position = "bottom")
    
    if (!is.null(title)) {
        p <- p + ggplot2::labs(title = title)
    }

    return(p)
}


#' Plot Timewise Posterior with Observations
#'
#' Plot posterior predictions with observed data overlaid.
#'
#' @param results List from `consolidate_results()`.
#' @param sim Simulation or data list with years and z (observations).
#' @param title Optional plot title (default NULL).
#'
#' @return ggplot2 object.
#'
#' @export
plot_timewise <- function(results, sim, title = NULL) {
    summaries <- results$timewise_summaries

    # Add years
    summaries <- summaries %>%
        dplyr::mutate(year = rep(sim$years, length(unique(H_fit))))

    # Add observations
    n_time <- length(sim$z)
    obs_df <- tibble::tibble(
        year = sim$years[1:n_time],
        obs = sim$z
    )

    p <- ggplot2::ggplot(summaries, ggplot2::aes(x = year, y = post_mean, color = factor(H_fit))) +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = post_ci_lower, ymax = post_ci_upper, fill = factor(H_fit)),
            alpha = 0.2, color = NA
        ) +
        ggplot2::geom_line(linewidth = 0.7) +
        ggplot2::geom_point(
            data = obs_df,
            ggplot2::aes(y = obs, color = "Observed"),
            size = 1.5, alpha = 0.6
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            x = "Year",
            y = "Predicted Value",
            color = "H Value",
            fill = "H Value"
        )
    
    if (!is.null(title)) {
        p <- p + ggplot2::labs(title = title)
    }

    return(p)
}


#' Plot Data with Posterior Samples or True Trend
#'
#' Create a plot of climate models with observations and either posterior samples or true trend.
#'
#' @param Y Matrix of climate model outputs (time x models).
#' @param z Vector of observations.
#' @param posterior_samples Optional matrix of posterior samples (samples x time).
#' @param true_trend Optional vector of true trend (e.g., mu from simulation).
#' @param title Optional plot title (default NULL).
#' @param years Optional vector of years for x-axis (default: 1:nrow(Y)).
#' @param obs_years Indices for where observations are available (default: 1:length(z)).
#' @param ylim Optional y-axis limits as c(min, max) (default NULL for automatic limits).
#'
#' @return ggplot2 object.
#'
#' @examples
#' \dontrun{
#' # With posterior samples
#' p <- plot_data_with_posterior(Y, z, posterior_samples, title = "Real Data Fit")
#'
#' # With true trend
#' p <- plot_data_with_posterior(Y, z, true_trend = mu, title = "Synthetic Data")
#'
#' # With specified y-axis limits
#' p <- plot_data_with_posterior(Y, z, posterior_samples, ylim = c(-2, 2))
#' }
#'
#' @export
plot_data_with_posterior <- function(Y, z, posterior_samples = NULL, true_trend = NULL,
                                     title = NULL, years = NULL,
                                     obs_years = 1:length(z), ylim = NULL) {
    T <- nrow(Y)
    K <- ncol(Y)

    if (is.null(years)) years <- 1:T

    # Create data frame for models
    model_df <- data.frame()
    for (k in 1:K) {
        model_df <- rbind(model_df, data.frame(
            year = years,
            value = Y[, k],
            type = ifelse(is.null(colnames(Y)), paste0("Model ", k), colnames(Y)[k])
        ))
    }

    # Create data frame for observations
    obs_df <- data.frame(
        year = years[obs_years],
        value = z,
        type = "Observation"
    )

    # Start with basic plot
    p <- ggplot2::ggplot() +
        ggplot2::geom_line(
            data = model_df, ggplot2::aes(x = year, y = value, group = type),
            color = "gray70", linewidth = 0.3
        ) +
        ggplot2::geom_point(
            data = obs_df, ggplot2::aes(x = year, y = value),
            color = "red", size = 1
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(x = "Year", y = "Temperature Anomaly") +
        ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank())
    
    if (!is.null(title)) {
        p <- p + ggplot2::labs(title = title)
    }

    # Add posterior samples if provided
    if (!is.null(posterior_samples)) {
        # posterior_samples is samples x time, so nrow = samples, ncol = time
        post_df <- data.frame()
        for (s in 1:nrow(posterior_samples)) {
            post_df <- rbind(post_df, data.frame(
                year = years,
                value = posterior_samples[s, ],
                sample = s
            ))
        }
        p <- p + ggplot2::geom_line(
            data = post_df, ggplot2::aes(x = year, y = value, group = sample),
            color = "blue", alpha = 0.1, linewidth = 1
        )
    }

    # Add true trend if provided
    if (!is.null(true_trend)) {
        trend_df <- data.frame(year = years, value = true_trend)
        p <- p + ggplot2::geom_line(
            data = trend_df, ggplot2::aes(x = year, y = value),
            color = "red", linewidth = 1
        )
    }
    
    # Apply y-axis limits if provided
    if (!is.null(ylim)) {
        p <- p + ggplot2::ylim(ylim[1], ylim[2])
    }

    return(p)
}


#' Plot Normalized Metrics with Elbow Detection
#'
#' Creates an overlay plot of RMSE, mean posterior variance, and max posterior variance,
#' all normalized to [0,1] scale with error bars. Displays elbow points for each metric.
#'
#' @param model_comparison Tibble with columns: H, RMSE, mean_post_var, max_post_var
#'   (typically from `consolidate_results()$replicate_comparison`).
#' @param title Optional plot title (default NULL).
#' @param show_legend Logical. Whether to display legend (default: TRUE).
#'
#' @return ggplot2 object showing normalized metrics overlaid with elbows marked by vertical lines.
#'
#' @details
#' The function:
#' 1. Normalizes all three metrics to [0,1] scale
#' 2. Computes standard errors if available (sd columns)
#' 3. Identifies elbow points for each metric independently
#' 4. Displays metrics as lines with error bars
#' 5. Marks elbows with colored vertical dotted lines:
#'    - Blue: RMSE elbow
#'    - Green: Mean posterior variance elbow
#'    - Orange: Max posterior variance elbow
#'
#' @examples
#' \dontrun{
#' results <- consolidate_results(fits, sims, by_replicate = TRUE)
#' p <- plot_normalized_metrics(results$replicate_comparison)
#' print(p)
#' }
#'
#' @export
plot_normalized_metrics <- function(model_comparison,
                                    title = NULL,
                                    show_legend = TRUE) {
    # Handle both RMSE (single replicate) and mean_RMSE (multi-replicate) columns
    if ("mean_RMSE" %in% names(model_comparison) && !("RMSE" %in% names(model_comparison))) {
        model_comparison <- model_comparison %>%
            dplyr::mutate(RMSE = mean_RMSE)
    }
    
    # Ensure we have required columns
    if (!all(c("H", "RMSE", "mean_post_var", "max_post_var") %in% names(model_comparison))) {
        stop("model_comparison must contain columns: H, RMSE (or mean_RMSE), mean_post_var, max_post_var")
    }

    # Prepare data with normalization
    perf_data <- model_comparison %>%
        dplyr::mutate(
            rmse_norm = (RMSE - min(RMSE)) / (max(RMSE) - min(RMSE)),
            mean_var_norm = (mean_post_var - min(mean_post_var)) / (max(mean_post_var) - min(mean_post_var)),
            max_var_norm = (max_post_var - min(max_post_var)) / (max(max_post_var) - min(max_post_var))
        ) %>%
        tidyr::pivot_longer(
            cols = c(rmse_norm, mean_var_norm, max_var_norm),
            names_to = "metric",
            values_to = "normalized_value"
        ) %>%
        dplyr::mutate(
            metric_label = dplyr::case_when(
                metric == "rmse_norm" ~ "RMSE",
                metric == "mean_var_norm" ~ "Mean Posterior Variance",
                metric == "max_var_norm" ~ "Max Posterior Variance"
            )
        )

    # Find elbows for each metric
    elbow_rmse <- find_elbow(
        model_comparison$H,
        model_comparison$RMSE,
        minimize = TRUE
    )

    elbow_mean_var <- find_elbow(
        model_comparison$H,
        model_comparison$mean_post_var,
        minimize = TRUE
    )

    elbow_max_var <- find_elbow(
        model_comparison$H,
        model_comparison$max_post_var,
        minimize = TRUE
    )

    # Create plot
    p <- ggplot2::ggplot(perf_data, ggplot2::aes(x = H, y = normalized_value, color = metric_label)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_vline(
            xintercept = elbow_rmse, color = "#1f77b4",
            linetype = "dotted", linewidth = 0.8, alpha = 0.6
        ) +
        ggplot2::geom_vline(
            xintercept = elbow_mean_var, color = "#2ca02c",
            linetype = "dotted", linewidth = 0.8, alpha = 0.6
        ) +
        ggplot2::geom_vline(
            xintercept = elbow_max_var, color = "#ff7f0e",
            linetype = "dotted", linewidth = 0.8, alpha = 0.6
        ) +
        ggplot2::scale_color_manual(
            values = c(
                "RMSE" = "#1f77b4",
                "Mean Posterior Variance" = "#2ca02c",
                "Max Posterior Variance" = "#ff7f0e"
            ),
            name = NULL
        ) +
        ggplot2::labs(
            x = "Number of Latent Factors (H)",
            y = "Normalized Value [0-1]"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            legend.position = if (show_legend) "bottom" else "none",
            plot.title = ggplot2::element_text(face = "bold", size = 12)
        )
    
    if (!is.null(title)) {
        p <- p + ggplot2::labs(title = title)
    }

    return(p)
}

#' Plot Average Posterior Means Across Replicates with True Trend
#'
#' Creates a plot showing the average posterior mean for each fitted H value
#' across multiple replications, overlaid with the true trend (mu) and simulated data.
#'
#' @param results_multi List from `consolidate_results()` with `by_replicate=TRUE`.
#' @param sims_list List of simulation objects (one per replicate), each containing
#'   Y, z, mu, and years.
#' @param title Optional plot title (default NULL).
#'
#' @return ggplot2 object.
#'
#' @details
#' The function:
#' 1. Computes the average posterior mean for each H across all replicates
#' 2. Plots these averages as colored lines
#' 3. Overlays the true trend (mu) as a red line
#' 4. Shows the averaged Y (climate models) from each replicate as grey lines
#'
#' @examples
#' \dontrun{
#' results <- consolidate_results(fits_list, sims_list, by_replicate = TRUE)
#' p <- plot_replicate_averages(results, sims_list)
#' print(p)
#' }
#'
#' @export
plot_replicate_averages <- function(results_multi, sims_list,
                                   title = NULL) {
    # Extract timewise results with replicate info
    timewise <- results_multi$timewise_results
    
    # Get true trend from first replicate (should be same for all)
    sim_ref <- sims_list[[1]]
    mu_true <- if (is.matrix(sim_ref$mu)) sim_ref$mu[, 1] else sim_ref$mu
    years <- sim_ref$years
    
    # Compute average posterior mean for each H across replicates
    avg_post <- timewise %>%
        dplyr::group_by(H_fit, time) %>%
        dplyr::summarise(
            avg_post_mean = mean(post_mean, na.rm = TRUE),
            .groups = "drop"
        )
    
    # Create data frames for plotting
    post_df <- avg_post %>%
        dplyr::mutate(
            year = years[time],
            H_label = paste0("H = ", H_fit)
        )
    
    true_df <- tibble::tibble(
        year = years,
        value = mu_true
    )
    
    # Compute average for each climate model across replicates
    # Each model (column of Y) gets averaged across all replicates
    n_models <- ncol(sims_list[[1]]$Y)
    model_avg_df <- data.frame()
    
    for (k in 1:n_models) {
        # Average this model across all replicates
        model_values <- sapply(sims_list, function(sim) sim$Y[, k])
        model_avg <- rowMeans(model_values, na.rm = TRUE)
        
        model_avg_df <- rbind(model_avg_df, data.frame(
            year = years,
            value = model_avg,
            model = k
        ))
    }
    
    # Create plot
    p <- ggplot2::ggplot() +
        ggplot2::geom_line(
            data = model_avg_df,
            ggplot2::aes(x = year, y = value, group = model),
            color = "gray70", linewidth = 0.3, alpha = 0.7
        ) +
        ggplot2::geom_line(
            data = true_df,
            ggplot2::aes(x = year, y = value),
            color = "red", linewidth = .8
        ) +
        ggplot2::geom_line(
            data = post_df,
            ggplot2::aes(x = year, y = avg_post_mean, color = H_label),
            linewidth = 1, alpha = 0.8
        ) +
        ggplot2::scale_color_brewer(palette = "Set1") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            x = "Year",
            y = "Temperature Anomaly",
            color = NULL
        ) +
        ggplot2::theme(
            legend.position = "bottom",
            plot.title = ggplot2::element_text(face = "bold", size = 12)
        )
    
    if (!is.null(title)) {
        p <- p + ggplot2::labs(title = title)
    }
    
    return(p)
}

#' Plot Variance of Posterior Means Across Replicates
#'
#' Creates a plot showing how the variance of posterior means across replicates
#' changes with H. Lower variance indicates more consistent predictions.
#'
#' @param results_multi List from `consolidate_results()` with `by_replicate=TRUE`.
#' @param title Optional plot title (default NULL).
#'
#' @return ggplot2 object.
#'
#' @details
#' For each H value, computes the variance of posterior means across replicates
#' at each time point, then averages these variances over time. Lower variance
#' indicates that the fitted H produces more consistent predictions across different
#' random realizations of the data.
#'
#' @examples
#' \dontrun{
#' results <- consolidate_results(fits_list, sims_list, by_replicate = TRUE)
#' p <- plot_posterior_variance(results)
#' print(p)
#' }
#'
#' @export
plot_posterior_variance <- function(results_multi,
                                   title = NULL) {
    # Extract timewise results with replicate info
    timewise <- results_multi$timewise_results
    
    # For each H and time, compute variance across replicates
    var_by_H <- timewise %>%
        dplyr::group_by(H_fit, time) %>%
        dplyr::summarise(
            var_post_mean = stats::var(post_mean, na.rm = TRUE),
            .groups = "drop"
        ) %>%
        dplyr::group_by(H_fit) %>%
        dplyr::summarise(
            mean_variance = mean(var_post_mean, na.rm = TRUE),
            sd_variance = stats::sd(var_post_mean, na.rm = TRUE),
            .groups = "drop"
        )
    
    # Create plot
    p <- ggplot2::ggplot(var_by_H, ggplot2::aes(x = H_fit, y = mean_variance)) +
        ggplot2::geom_line(color = "steelblue", linewidth = 1) +
        ggplot2::geom_point(color = "steelblue", size = 3) +
        ggplot2::geom_errorbar(
            ggplot2::aes(ymin = mean_variance - sd_variance, 
                        ymax = mean_variance + sd_variance),
            width = 0.2, color = "steelblue", linewidth = 0.8
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            x = "Number of Latent Factors (H)",
            y = "Mean Variance of Posterior Means"
        ) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold", size = 12),
            plot.subtitle = ggplot2::element_text(size = 10, color = "gray50")
        )
    
    if (!is.null(title)) {
        p <- p + ggplot2::labs(title = title)
    }
    
    return(p)
}
