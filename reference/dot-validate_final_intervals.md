# Validate final split intervals

Checks that all final intervals satisfy the size limit constraint.
Raises an error if any interval exceeds the limit.

## Usage

``` r
.validate_final_intervals(boundaries, chr_end, limit, chr)
```

## Arguments

- boundaries:

  Integer vector: boundary positions (interior points only).

- chr_end:

  Integer: chromosome end coordinate.

- limit:

  Integer: maximum allowed interval width.

- chr:

  Character: chromosome name for error messages.

## Value

Invisible `TRUE` if all intervals satisfy the limit. Raises error
otherwise.

## Examples

``` r
if (FALSE) { # \dontrun{
  .validate_final_intervals(c(100L, 200L), 300L, 50L, "chr1")
} # }
```
