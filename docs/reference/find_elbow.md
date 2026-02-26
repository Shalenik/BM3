# Find Elbow Point in a Curve (Perpendicular Distance Method)

Identifies the "elbow" point in a metric-vs-H curve using the
perpendicular distance method. Useful for model selection when fitting
multiple H values (e.g., identifying the optimal H by finding where RMSE
or variance improvements plateau).

## Usage

``` r
find_elbow(x, y, minimize = TRUE)
```

## Arguments

- x:

  Numeric vector of x-values (typically H values).

- y:

  Numeric vector of y-values (typically RMSE or variance).

- minimize:

  Logical. If TRUE (default), assumes lower y is better (e.g., RMSE). If
  FALSE, assumes higher y is better (e.g., likelihood).

## Value

Numeric. The x-value (H) at the identified elbow point.

## Details

The elbow method finds where the curve transitions from steep to flat.
The algorithm:

1.  Identifies the decreasing region of the curve

2.  Normalizes x and y to 0,1 scale

3.  Draws a line from the start to the end of the decreasing region

4.  Calculates perpendicular distance from each point to this line

5.  Returns the x-value (H) with maximum perpendicular distance

This is particularly useful for model selection when the optimal H is
unclear from the metric alone.

## Examples

``` r
if (FALSE) { # \dontrun{
# Find elbow in RMSE curve
H_values <- 2:7
rmse_values <- c(0.85, 0.42, 0.38, 0.36, 0.35, 0.35)
elbow_h <- find_elbow(H_values, rmse_values, minimize = TRUE)
print(elbow_h) # Likely returns H around 3-4
} # }
```
