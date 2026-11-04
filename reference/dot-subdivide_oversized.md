# Subdivide oversized intervals

Identifies intervals between consecutive points that exceed the size
limit and generates interior boundaries to split them. Returns all extra
boundaries needed.

## Usage

``` r
.subdivide_oversized(points, limit)
```

## Arguments

- points:

  Integer vector: all boundary points (including endpoints 0 and
  chr_end) in sorted order.

- limit:

  Integer: maximum allowed interval width.

## Value

Integer vector of extra interior boundaries needed. Empty if no
intervals exceed the limit.

## Examples

``` r
if (FALSE) { # \dontrun{
  extra <- .subdivide_oversized(c(0L, 100L, 500L), limit = 100L)
  # Interval [100, 500] is 400bp, needs splitting
} # }
```
