# Abort if any GTF features match no split region

Abort if any GTF features match no split region

## Usage

``` r
.abort_if_unmatched_features(matched, gtf)
```

## Arguments

- matched:

  Data frame containing `unique_ID`.

- gtf:

  Original GTF data frame containing `unique_ID`, `Chr`, `feature`,
  `start`, and `end`.

## Value

Invisibly `TRUE` if ok; otherwise aborts.
