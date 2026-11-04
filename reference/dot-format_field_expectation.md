# Format field requirement description

Creates a human-readable string describing the allowed field count
range.

## Usage

``` r
.format_field_expectation(min_fields, max_fields)
```

## Arguments

- min_fields:

  Integer: minimum fields.

- max_fields:

  Integer or NA: maximum fields. If NA, no upper limit.

## Value

Character: description such as "3-9 TAB-separated fields" or "at least 9
TAB-separated fields".

## Examples

``` r
if (FALSE) { # \dontrun{
  .format_field_expectation(3L, 9L)
  # Returns: "3-9 TAB-separated fields"
  .format_field_expectation(9L, NA_integer_)
  # Returns: "at least 9 TAB-separated fields"
} # }
```
