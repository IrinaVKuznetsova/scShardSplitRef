# Convert boundary points to interval data frame

Transforms a vector of boundary positions into a data frame representing
intervals \[start, end) between consecutive boundaries.

## Usage

``` r
.boundaries_to_df(chr, start, end, boundaries)
```

## Arguments

- chr:

  Character: chromosome name for all rows.

- start:

  Integer: the minimum start coordinate (first boundary).

- end:

  Integer: the maximum end coordinate (last boundary).

- boundaries:

  Integer vector: interior boundary points (may be empty).

## Value

Data frame with columns 'chr', 'start', 'end'. Each row represents an
interval between consecutive boundary points.

## Examples

``` r
if (FALSE) { # \dontrun{
  df <- .boundaries_to_df("chr1", 0L, 100L, c(30L, 60L))
  # Returns:
  #   chr start end
  # 1 chr1    0  30
  # 2 chr1   30  60
  # 3 chr1   60 100
} # }
```
