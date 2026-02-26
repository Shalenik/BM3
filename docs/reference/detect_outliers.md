# Detect Outlier Climate Models Using RMSE Thresholding

Iteratively removes climate model columns whose RMSE from the ensemble
mean exceeds a user-defined threshold. Assumes input is a matrix or data
frame with time series in rows and models in columns.

## Usage

``` r
detect_outliers(
  data,
  start_index,
  end_index,
  threshold_fn = mean_sd_threshold,
  threshold_dial = 2.5
)
```

## Arguments

- data:

  A numeric matrix or data frame with model outputs.

- start_index:

  Integer, starting row index (e.g., corresponding to 2015).

- end_index:

  Integer, ending row index (e.g., corresponding to 2100).

- threshold_fn:

  A function taking a numeric vector of RMSE values and returning an
  RMSE threshold. Default is
  [`mean_sd_threshold()`](https://shalenik.github.io/BM3/reference/mean_sd_threshold.md).

- threshold_dial:

  Numeric multiplier passed to `threshold_fn`. Default is 2.5.

## Value

A character vector of column names identified as outliers.
