# Fix split boundaries to avoid genes

Adjusts interval boundaries to ensure none fall within gene annotations.
Iteratively shifts boundaries past genes and subdivides oversized
intervals until all constraints are satisfied.

## Usage

``` r
.fix_split_boundaries(
  split_df,
  gene_ranges,
  limit,
  shift_by = 1L,
  clearance = 1L
)
```

## Arguments

- split_df:

  Data frame with columns 'chr', 'start', 'end' representing initial
  split intervals.

- gene_ranges:

  Data frame with already-merged, non-overlapping gene ranges for the
  same chromosome (columns 'start', 'end'), as produced by
  [`.merge_all_gene_ranges()`](https://IrinaVKuznetsova.github.io/scShardSplitRef/reference/dot-merge_all_gene_ranges.md).
  Passed in pre-merged so this function does not redo the sort-and-merge
  work for every region on the chromosome.

- limit:

  Integer: maximum allowed interval width in base pairs.

- shift_by:

  Integer: step size (bp) used to walk boundaries rightward when
  clearance is violated. Defaults to 1L.

- clearance:

  Integer: minimum required clearance (bp) from nearby genes on both
  sides. Defaults to 1L.

## Value

Data frame with same structure as `split_df` but with adjusted
boundaries that satisfy the size and gene-avoidance constraints.

## Details

The algorithm iteratively:

1.  Shifts boundaries rightward in `shift_by` bp increments until they
    satisfy the `clearance` gene distance on both sides

2.  Subdivides intervals exceeding the limit

3.  Shifts again after subdivision

4.  Repeats until convergence

Raises an error if constraints cannot be satisfied.

## Examples

``` r
if (FALSE) { # \dontrun{
  fixed_df <- fix_split_boundaries(
    split_df, gene_ranges, limit = 50000, shift_by = 1, clearance = 1
  )
} # }
```
