# Plot Timewise Posterior with Observations

Plot posterior predictions with observed data overlaid.

## Usage

``` r
plot_timewise(results, sim, title = NULL)
```

## Arguments

- results:

  List from
  [`consolidate_results()`](https://shalenik.github.io/BM3/reference/consolidate_results.md).

- sim:

  Simulation or data list with years and z (observations).

- title:

  Optional plot title (default NULL).

## Value

ggplot2 object.
