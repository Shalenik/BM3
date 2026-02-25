#' Plot Posterior Distribution
#'
#' Create a plot of posterior means and credible intervals.
#'
#' @param results List from `consolidate_results()`.
#' @param H Integer. Which H value to plot.
#' @param data Optional. Observations to overlay (default NULL).
#'
#' @return ggplot2 object.
#'
#' @examples
#' \dontrun{
#' plot_posterior(results, H = 2)
#' }
#'
#' @export
plot_posterior <- function(results, H, data = NULL) {
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
            title = sprintf("Posterior Distribution (H=%d)", H),
            x = "Time",
            y = "Predicted Value",
            subtitle = sprintf(
                "RMSE: %.4f",
                unique(summaries$rmse)
            )
        )

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
#'
#' @return ggplot2 object.
#'
#' @export
plot_comparison <- function(results) {
    comparison <- results$model_comparison %>%
        dplyr::select(H, RMSE) %>%
        dplyr::mutate(RMSE = scales::rescale(RMSE, to = c(0, 1)))

    p <- ggplot2::ggplot(comparison, ggplot2::aes(x = H, y = RMSE)) +
        ggplot2::geom_point(size = 3, color = "steelblue") +
        ggplot2::geom_line(color = "steelblue") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = "Model Comparison Across H Values",
            x = "Number of Latent Factors (H)",
            y = "RMSE (normalized)"
        ) +
        ggplot2::theme(legend.position = "bottom")

    return(p)
}


#' Plot Timewise Posterior with Observations
#'
#' Plot posterior predictions with observed data overlaid.
#'
#' @param results List from `consolidate_results()`.
#' @param sim Simulation or data list with years and z (observations).
#'
#' @return ggplot2 object.
#'
#' @export
plot_timewise <- function(results, sim) {
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
            title = "Posterior Predictions Across H Values",
            x = "Year",
            y = "Predicted Value",
            color = "H Value",
            fill = "H Value"
        )

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
#' @param title Plot title.
#' @param years Optional vector of years for x-axis (default: 1:nrow(Y)).
#' @param obs_years Indices for where observations are available (default: 1:length(z)).
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
#' }
#'
#' @export
plot_data_with_posterior <- function(Y, z, posterior_samples = NULL, true_trend = NULL,
                                     title = "Data and Posterior", years = NULL,
                                     obs_years = 1:length(z)) {
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
            color = "red", size = 3
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = title, x = "Year", y = "Temperature Anomaly") +
        ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank())

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
            color = "blue", alpha = 0.1, linewidth = 1.5
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

    return(p)
}


#' Plot Normalized Metrics with Elbow Detection
#'
#' Creates an overlay plot of RMSE, mean posterior variance, and max posterior variance,
#' all normalized to [0,1] scale with error bars. Displays elbow points for each metric.
#'
#' @param model_comparison Tibble with columns: H, RMSE, mean_post_var, max_post_var
#'   (typically from `consolidate_results()$replicate_comparison`).
#' @param title Plot title (default: "Normalized Metrics with Elbow Detection").
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
                                    title = "Normalized Metrics with Elbow Detection",
                                    show_legend = TRUE) {
    # Ensure we have required columns
    if (!all(c("H", "RMSE", "mean_post_var", "max_post_var") %in% names(model_comparison))) {
        stop("model_comparison must contain columns: H, RMSE, mean_post_var, max_post_var")
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
            title = title,
            x = "Number of Latent Factors (H)",
            y = "Normalized Value [0-1]"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            legend.position = if (show_legend) "bottom" else "none",
            plot.title = ggplot2::element_text(face = "bold", size = 12)
        )

    return(p)
}
