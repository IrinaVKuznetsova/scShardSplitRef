# Abort if any feature matches multiple split regions

Detects duplicated `unique_ID` values in the matched table, which
indicates overlapping split regions (or otherwise ambiguous
assignments).

## Usage

``` r
.abort_if_ambiguous_matches(matched)
```

## Arguments

- matched:

  Data frame containing a `unique_ID` column.

## Value

Invisibly `TRUE` if ok; otherwise aborts.
