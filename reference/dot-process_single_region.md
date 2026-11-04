# Process a single genomic region

Splits a single genomic region according to the size limit and adjusts
boundaries to avoid genes. If the region has a non-zero start
(centromere), it is split into two parts: before and after the
centromere.

## Usage

``` r
.process_single_region(
  region,
  merged_gene_ranges,
  limit,
  shift_by,
  clearance = 1L
)
```

## Arguments

- region:

  Data frame row with columns 'chr', 'start', 'end'.

- merged_gene_ranges:

  Named list, keyed by chromosome, of data frames with columns 'start',
  'end' – gene ranges already merged via
  [`.merge_all_gene_ranges()`](https://IrinaVKuznetsova.github.io/scShardSplitRef/reference/dot-merge_all_gene_ranges.md).
  Passed in rather than recomputed here so the merge work happens once
  per chromosome regardless of how many regions share it.

- limit:

  Integer: maximum allowed interval width in base pairs.

- shift_by:

  Integer: step size (bp) used to walk boundaries rightward when
  clearance is violated.

- clearance:

  Integer: minimum required clearance (bp) from nearby genes on both
  sides.

## Value

Data frame with columns 'chr', 'start', 'end' representing the processed
split intervals for this region.

## Details

If the region start is \> 0, the function assumes this marks a
centromere and splits the region into \[0, start) and \[start, end).
Each part is independently split according to the size limit. Boundaries
are then adjusted to avoid genes in the region.

## Examples

``` r
if (FALSE) { # \dontrun{
  region <- data.frame(chr = "chr1", start = 100000, end = 200000)
  genes_subset <- genes[genes$chr == "chr1", ]
  result <- .process_single_region(region, genes_subset,
                                   limit = 50000,
                                   shift_by = 1,
                                   clearance = 1)
} # }
```
