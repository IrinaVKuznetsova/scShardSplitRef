# Generate split boundaries for an interval

Calculates interior boundary points that divide a coordinate interval
\[start, end) into chunks of approximately equal size, each at most
`limit` base pairs wide.

## Usage

``` r
.build_boundaries(start, end, limit)
```

## Arguments

- start:

  Integer: interval start (inclusive).

- end:

  Integer: interval end (exclusive).

- limit:

  Integer: maximum width per chunk in base pairs.

## Value

Integer vector of boundary positions. Empty if interval width \<= limit.

## Details

The function calculates the number of chunks needed, then generates
evenly-spaced boundaries. Interior boundaries are rounded to the nearest
integer. The first and last boundaries are excluded (they are the
interval endpoints).

## Examples

``` r
if (FALSE) { # \dontrun{
  boundaries <- .build_boundaries(0L, 1000L, 250L)
  # Returns something like: 250, 500, 750
} # }
```
