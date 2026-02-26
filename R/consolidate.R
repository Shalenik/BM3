#' Consolidate Fitting Results
#'
#' Create summary statistics from fitted BFACT models across multiple H values.
#' Combines individual fit objects into posterior summaries and model comparison metrics.
#'
#' @param fits List of fit objects (output from `fit_bfact_model()`), or a list of such lists
#'   for multiple replications. See details.
#' @param sim List containing simulation data (Y, z, years), or list of such lists
#'   for multiple replications.
#' @param by_replicate Logical. If TRUE, compute per-replicate statistics (default FALSE).
#'   If fits and sim are lists of lists, this will aggregate across replications.
#'
#' @return A list containing:
#'   - timewise_summaries: tibble with posterior means/CIs per H
#'   - model_comparison: tibble comparing metrics across H values
#'   - per_replicate_stats: per-replicate RMSE if by_replicate=TRUE
#'   - replicate_comparison: tibble with mean/sd RMSE by H across replications (if multiple reps)
#'
#' @details
#' Consolidates multiple H fits into a single summary table with columns:
#'   - H_fit: fitted H value
#'   - time: time index
#'   - post_mean: posterior mean prediction
#'   - post_ci_lower, post_ci_upper: 95% credible interval
#'   - se: squared error at each time point
#'   - rmse: overall RMSE for that H
#'
#' For multiple replications, pass:
#'   - fits: list of lists, where each inner list is fits for one replicate
#'   - sim: list of simulation lists, one per replicate
#'
#' @examples
#' \dontrun{
#' sim <- simulate_data(H_true = 2, K = 20)
#' fit_H2 <- fit_bfact_model(sim, H = 2)
#' fit_H3 <- fit_bfact_model(sim, H = 3)
#'
#' results <- consolidate_results(
#'     fits = list(fit_H2, fit_H3),
#'     sim = sim
#' )
#' head(results$timewise_summaries)
#' }
#'
#' @export
consolidate_results <- function(fits, sim, by_replicate = FALSE) {
    # Check if this is multi-replicate setup (fits[[1]] is a list of fits)
    is_multi_rep <- is.list(fits[[1]]) && is.list(sim[[1]])

    if (is_multi_rep && by_replicate) {
        # Multi-replicate case: aggregate across replications
        message(sprintf("Consolidating %d replicates with multiple H fits...", length(fits)))

        # Helper function to compute posterior variance from fit object
        compute_post_var <- function(fit) {
            if (!is.null(fit$ysim)) {
                ysim_mat <- as.matrix(fit$ysim)
                if (ncol(ysim_mat) < 2) {
                    return(NA_real_)
                }
                # Use only second half of posterior samples (after burn-in)
                burn_start <- floor(nrow(ysim_mat) / 2) + 1
                pred_mat <- ysim_mat[burn_start:nrow(ysim_mat), -1, drop = FALSE]
                # Variance at each time point (across samples), then average over time
                post_vars <- apply(pred_mat, 2, stats::var, na.rm = TRUE)
                return(mean(post_vars, na.rm = TRUE))
            }
            return(NA_real_)
        }

        rep_results_list <- mapply(function(fits_rep, sim_rep, rep_idx) {
            # Call consolidate_results on each replicate
            results_rep <- consolidate_results(fits_rep, sim_rep, by_replicate = FALSE)
            results_rep$timewise_summaries %>%
                dplyr::mutate(replicate = rep_idx, .before = 1)
        }, fits, sim, seq_along(fits), SIMPLIFY = FALSE)

        timewise_results <- dplyr::bind_rows(rep_results_list)

        # Extract model comparison metrics
        rep_results <- purrr::map_df(seq_along(fits), function(rep_idx) {
            results_rep <- consolidate_results(fits[[rep_idx]], sim[[rep_idx]], by_replicate = FALSE)
            results_rep$model_comparison %>%
                dplyr::mutate(replicate = rep_idx, .before = 1)
        })

        # Compute posterior variances for each fit
        post_vars_by_fit <- list()
        for (rep_idx in seq_along(fits)) {
            for (h_idx in seq_along(fits[[rep_idx]])) {
                fit <- fits[[rep_idx]][[h_idx]]
                H <- if (!is.null(fit$H)) fit$H else dim(fit$U_samples)[3]
                post_vars_by_fit[[paste0("r", rep_idx, "_H", H)]] <- compute_post_var(fit)
            }
        }

        # Create a tibble of posterior variances
        post_vars_tbl <- tibble::tibble(
            replicate = rep(seq_along(fits), each = length(fits[[1]])),
            H = rep(sapply(fits[[1]], function(f) {
                if (!is.null(f$H)) f$H else dim(f$U_samples)[3]
            }), length(fits)),
            post_var = unlist(post_vars_by_fit)
        )

        # Aggregate metrics across replicates
        replicate_comparison <- rep_results %>%
            dplyr::group_by(H) %>%
            dplyr::summarise(
                mean_RMSE = mean(RMSE, na.rm = TRUE),
                sd_RMSE = stats::sd(RMSE, na.rm = TRUE),
                min_RMSE = min(RMSE, na.rm = TRUE),
                max_RMSE = max(RMSE, na.rm = TRUE),
                n_reps = dplyr::n(),
                .groups = "drop"
            ) %>%
            dplyr::left_join(
                post_vars_tbl %>%
                    dplyr::group_by(H) %>%
                    dplyr::summarise(
                        mean_post_var = mean(post_var, na.rm = TRUE),
                        max_post_var = max(post_var, na.rm = TRUE),
                        .groups = "drop"
                    ),
                by = "H"
            ) %>%
            dplyr::arrange(mean_RMSE)

        message("  Replicate comparison (aggregated across all replicates):")
        print(replicate_comparison)

        result <- list(
            replicate_comparison = replicate_comparison,
            replicate_results = rep_results,
            timewise_results = timewise_results,
            post_var_by_rep = post_vars_tbl,
            fits_all = fits,
            sims_all = sim
        )
        return(result)
    }

    # Single replicate case: proceed with original logic
    if (!is.list(fits[[1]])) fits <- list(fits)

    message(sprintf("Consolidating %d fit(s)...", length(fits)))

    get_H <- function(fit) {
        if (!is.null(fit$H)) {
            return(fit$H)
        }
        if (!is.null(fit$U_samples)) {
            return(dim(fit$U_samples)[3])
        }
        return(NA_integer_)
    }

    get_nsim <- function(fit) {
        if (!is.null(fit$nsim)) {
            return(fit$nsim)
        }
        if (!is.null(fit$beta_samples)) {
            return(dim(fit$beta_samples)[1])
        }
        return(NA_integer_)
    }

    get_posterior_summary <- function(fit) {
        if (!is.null(fit$posterior_summary)) {
            return(tibble::as_tibble(fit$posterior_summary))
        }

        if (!is.null(fit$ysim)) {
            ysim_mat <- as.matrix(fit$ysim)
            if (ncol(ysim_mat) < 2) {
                stop("fit$ysim must have at least two columns (precision + mean curve)")
            }
            # Use only second half of posterior samples (after burn-in)
            burn_start <- floor(nrow(ysim_mat) / 2) + 1
            pred_mat <- ysim_mat[burn_start:nrow(ysim_mat), -1, drop = FALSE]

            post_mean <- colMeans(pred_mat, na.rm = TRUE)
            post_ci_lower <- apply(pred_mat, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
            post_ci_upper <- apply(pred_mat, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
            return(tibble::tibble(
                time = seq_along(post_mean),
                post_mean = post_mean,
                post_ci_lower = post_ci_lower,
                post_ci_upper = post_ci_upper
            ))
        }

        stop("No posterior summary available on fit object.")
    }

    compute_rmse <- function(post_summary) {
        # Prefer full series (mu) if available, otherwise use observations (z)
        if (!is.null(sim$mu)) {
            # Use full series
            mu_full <- if (is.matrix(sim$mu)) sim$mu[, 1] else sim$mu
            if (length(mu_full) != nrow(post_summary)) {
                return(NA_real_)
            }
            return(sqrt(mean((post_summary$post_mean - mu_full)^2, na.rm = TRUE)))
        }

        if (is.null(sim$z)) {
            return(NA_real_)
        }

        # Fall back to observation window
        obs_idx <- if (!is.null(sim$T1) && !is.null(sim$T2)) {
            sim$T1:sim$T2
        } else {
            seq_len(length(sim$z))
        }
        obs_idx <- obs_idx[obs_idx <= nrow(post_summary)]
        if (length(obs_idx) == 0) {
            return(NA_real_)
        }
        mu_obs <- post_summary$post_mean[obs_idx]
        z_obs <- sim$z[seq_len(length(mu_obs))]
        sqrt(mean((z_obs - mu_obs)^2, na.rm = TRUE))
    }

    # Extract H values
    H_values <- sapply(fits, get_H)

    # Combine posterior summaries
    timewise_summaries <- purrr::map_df(fits, function(fit) {
        summary_tbl <- get_posterior_summary(fit)
        rmse <- compute_rmse(summary_tbl)

        summary_tbl %>%
            dplyr::mutate(
                H_fit = get_H(fit),
                rmse = rmse,
                .before = 1
            )
    })


    # Compute posterior variance metrics for each fit (after burn-in)
    calc_post_vars <- function(fit) {
        if (!is.null(fit$ysim)) {
            ysim_mat <- as.matrix(fit$ysim)
            burn_start <- floor(nrow(ysim_mat) / 2) + 1
            pred_mat <- ysim_mat[burn_start:nrow(ysim_mat), -1, drop = FALSE]
            vars <- apply(pred_mat, 2, stats::var, na.rm = TRUE)
            return(list(
                mean_post_var = mean(vars, na.rm = TRUE),
                max_post_var = max(vars, na.rm = TRUE)
            ))
        }
        list(mean_post_var = NA_real_, max_post_var = NA_real_)
    }

    post_var_tbl <- lapply(fits, calc_post_vars)

    model_comparison <- tibble::tibble(
        H = H_values,
        RMSE = sapply(fits, function(f) compute_rmse(get_posterior_summary(f))),
        mean_post_var = vapply(post_var_tbl, function(x) x$mean_post_var, numeric(1)),
        max_post_var = vapply(post_var_tbl, function(x) x$max_post_var, numeric(1))
    ) %>%
        dplyr::arrange(H)

    result <- list(
        timewise_summaries = timewise_summaries,
        model_comparison = model_comparison,
        fits = fits
    )

    if (by_replicate) {
        # Extract per-replicate statistics
        per_rep_stats <- purrr::map_df(fits, function(fit) {
            tibble::tibble(
                H_fit = get_H(fit),
                n_samples = get_nsim(fit),
                rmse = compute_rmse(get_posterior_summary(fit))
            )
        })
        result$per_replicate_stats <- per_rep_stats
    }

    return(result)
}
