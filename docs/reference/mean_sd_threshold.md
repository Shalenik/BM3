# RMSE Threshold Function (Mean + k SD)

Computes a threshold for outlier detection as the mean plus `dial`
standard deviations.

## Usage

``` r
mean_sd_threshold(rmse_vals, dial = 2.5)
```

## Arguments

- rmse_vals:

  A numeric vector of RMSE values.

- dial:

  Numeric multiplier for the standard deviation. Default is 2.5.

## Value

A numeric threshold value.
