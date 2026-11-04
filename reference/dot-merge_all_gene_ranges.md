# Merge gene ranges for every chromosome, once

Groups gene annotations by chromosome and merges overlapping ranges
within each group, producing a lookup used by both the upfront
feasibility check and per-region boundary fixing – so the merge work
happens exactly once per chromosome, rather than being recomputed by
every downstream consumer.

## Usage

``` r
.merge_all_gene_ranges(genes)
```

## Arguments

- genes:

  Data frame with gene annotations (columns 'chr', 'start', 'end').

## Value

Named list, keyed by chromosome, of data frames with columns 'start',
'end' representing merged, non-overlapping gene ranges. Empty list if
`genes` has 0 rows.

## Examples

``` r
if (FALSE) { # \dontrun{
  merged_gene_ranges <- .merge_all_gene_ranges(genes)
  merged_gene_ranges[["chr1"]]
} # }
```
