# Merge overlapping gene ranges

Combines overlapping or adjacent gene intervals on a single chromosome
into merged ranges. Gaps between genes are preserved.

## Usage

``` r
.merge_overlapping_gene_ranges(chr_genes)
```

## Arguments

- chr_genes:

  Data frame with gene coordinates for one chromosome (columns 'start',
  'end'). May contain duplicates and overlaps.

## Value

Data frame with columns 'start', 'end' representing merged gene ranges
with no overlaps. Rows are sorted by start position. Empty data frame if
input has 0 rows.

## Details

Algorithm:

1.  Sorts genes by start position

2.  Iterates through genes, merging overlapping ranges

3.  Returns list of disjoint merged intervals

## Examples

``` r
if (FALSE) { # \dontrun{
  merged <- .merge_overlapping_gene_ranges(chr_genes)
} # }
```
