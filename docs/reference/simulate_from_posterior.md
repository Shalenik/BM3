# Simulate Data from a Posterior Draw of a Fitted BFACT Object

Generates a simulated dataset from a single posterior draw of a fitted
BFACT object.

## Usage

``` r
simulate_from_posterior(
  fit,
  s = NULL,
  years_start = 1850,
  years_end = 2100,
  T2 = 173,
  T1 = 1,
  J = 6,
  OBS_COL = 1,
  seed = NULL
)
```

## Arguments

- fit:

  A fitted BFACT object (as returned by BFACT()).

- s:

  Integer index of the posterior draw to use (default: random draw).

- years_start:

  First year (default: 1850).

- years_end:

  Last year (default: 2100).

- T2:

  Last observed time index (default: 173).

- T1:

  First observed time index (default: 1).

- J:

  Number of spline basis functions (default: 6).

- OBS_COL:

  Which column of Y is the observed series (default: 1).

- seed:

  Random seed for reproducibility (default: NULL).

## Value

A list with simulated data and metadata.
