# Plot Normalized Metrics with Elbow Detection

Creates an overlay plot of RMSE, mean posterior variance, and max
posterior variance, all normalized to 0,1 scale with error bars.
Displays elbow points for each metric.

## Usage

``` r
plot_normalized_metrics(model_comparison, title = NULL, show_legend = TRUE)
```

## Arguments

- model_comparison:

  Tibble with columns: H, RMSE, mean_post_var, max_post_var (typically
  from `consolidate_results()$replicate_comparison`).

- title:

  Optional plot title (default NULL).

- show_legend:

  Logical. Whether to display legend (default: TRUE).

## Value

ggplot2 object showing normalized metrics overlaid with elbows marked by
vertical lines.

## Details

The function:

1.  Normalizes all three metrics to 0,1 scale

2.  Computes standard errors if available (sd columns)

3.  Identifies elbow points for each metric independently

4.  Displays metrics as lines with error bars

5.  Marks elbows with colored vertical dotted lines:

    - Blue: RMSE elbow

    - Green: Mean posterior variance elbow

    - Orange: Max posterior variance elbow

## Examples

``` r
if (FALSE) { # \dontrun{
results <- consolidate_results(fits, sims, by_replicate = TRUE)
p <- plot_normalized_metrics(results$replicate_comparison)
print(p)
} # }
```
