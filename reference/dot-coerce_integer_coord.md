# Coerce a coordinate column to integer, safely

Genomic coordinates must always be whole-number `integer`, never
`double`. Plain [`as.integer()`](https://rdrr.io/r/base/integer.html) is
not safe to use directly on trusted input because it silently returns
`NA` (with only a warning) for values outside R's 32-bit integer range,
and silently truncates fractional values. This helper fails loudly
instead in both cases so a type problem is caught immediately at the
point of reading, rather than surfacing much later as a
scientific-notation string or a mysteriously missing row.

## Usage

``` r
.coerce_integer_coord(x, label = "coordinate")
```

## Arguments

- x:

  Numeric vector (integer, double, or character coercible to numeric) of
  coordinate values.

- label:

  Character scalar used in error messages to identify which column
  failed (e.g. `"RegionStart"`).

## Value

Integer vector, same length as `x`.

## Examples

``` r
if (FALSE) { # \dontrun{
  .coerce_integer_coord(c(0, 100000000), "RegionStart")
  .coerce_integer_coord(3.5, "RegionStart") # errors: not whole numbers
} # }
```
