# Find lines with invalid field counts

Identifies line indices where the field count violates min/max
constraints.

## Usage

``` r
.find_invalid_field_counts(field_counts, min_fields, max_fields = NA_integer_)
```

## Arguments

- field_counts:

  Integer vector: number of fields per line.

- min_fields:

  Integer: minimum required fields.

- max_fields:

  Integer or NA: maximum allowed fields.

## Value

Integer vector: indices of lines with invalid field counts, sorted and
deduplicated. Empty vector if all lines are valid.

## Examples

``` r
if (FALSE) { # \dontrun{
  bad <- .find_invalid_field_counts(c(2, 3, 2), 3L, NA_integer_)
} # }
```
