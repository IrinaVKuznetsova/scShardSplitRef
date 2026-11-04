# Validate a TAB-delimited file and return its data lines

Performs structural validation (i.e., extension check,
comment/blank-line filtering, field-count range check), and returns the
validated, filtered lines.

## Usage

``` r
.validate_tab_file_lines(path, min_fields, label, max_fields = NA_integer_)
```

## Arguments

- path:

  Character: file path to validate.

- min_fields:

  Integer: minimum number of required fields per line.

- label:

  Character: human-readable file type label (*e.g.*, "GTF", "BED") for
  error messages.

- max_fields:

  Integer or `NA`: maximum number of allowed fields. If `NA` (default),
  no upper limit is enforced.

## Value

Character vector of validated data lines (comments and blank lines
removed). Empty character vector if the file has no data lines.
