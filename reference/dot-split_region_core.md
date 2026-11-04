# Split a single coordinate interval

Splits a coordinate interval \[start, end) into chunks of at most
`limit` base pairs. Returns a data frame with the resulting intervals.

## Usage

``` r
.split_region_core(chr, start, end, limit)
```

## Arguments

- chr:

  Character: chromosome name.

- start:

  Integer: interval start (inclusive).

- end:

  Integer: interval end (exclusive).

- limit:

  Integer: maximum width per interval in base pairs.

## Value

Data frame with columns 'chr', 'start', 'end' representing the split
intervals. If the interval width \<= limit, a single row is returned.

## Examples

``` r
if (FALSE) { # \dontrun{
  result <- .split_region_core("chr1", 0L, 1000L, 250L)
} # }
```
