# Format a coordinate value as a plain fixed-point integer string

`sprintf("%s", x)` and other implicit string coercions of a *numeric*
(double) vector go through
[`as.character()`](https://rdrr.io/r/base/character.html), which follows
`getOption("scipen")`/`getOption("digits")` and will happily render
large, round genomic coordinates in scientific notation (e.g. `1e+08`
instead of `100000000`). Integers never do this, but
[`read.table()`](https://rdrr.io/r/utils/read.table.html)-derived or
arithmetic-derived coordinate columns are not guaranteed to stay integer
(e.g. `matched$start - matched$RegionStart + 1` promotes to double).

## Usage

``` r
.format_coord(x)
```

## Arguments

- x:

  Numeric (integer or double) vector of coordinates. `NA` values are
  passed through as `NA` (not the string `"NA"`).

## Value

Character vector, same length as `x`, with each value rendered as a
plain integer string (no scientific notation, no decimal point).

## Details

This helper forces fixed-point, non-scientific formatting regardless of
whether `x` is integer or double, and is safe for coordinates well
within double's exact-integer range (+/- 2^53), which comfortably covers
any genomic coordinate.

## Examples

``` r
if (FALSE) { # \dontrun{
  .format_coord(1e8)      # "100000000", not "1e+08"
  .format_coord(100000000L) # "100000000"
} # }
```
