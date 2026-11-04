# Determine split regions for reference sharding

Builds BED-like split intervals from genomic regions and GTF gene
annotations. If the start position in a genomic region is not 0, it is
assumed to mark the start of a centromere and the region is split around
that point. Proposed boundaries are shifted rightward when they fall too
close to gene intervals, and oversized chunks are repeatedly re-split
until all segments satisfy the requested size limit.

## Usage

``` r
determine_split_regions(
  bed,
  gtf,
  output_bed,
  limit = 2L^29L,
  shift_by = 1L,
  clearance = 1L
)
```

## Arguments

- bed:

  Path to a 3-column BED-like file (chr, start, end). Columns should be
  TAB-delimited with chromosome name, region start (0-based), and region
  end coordinates.

- gtf:

  Path to a GTF file containing gene annotations. Must contain standard
  9 GTF columns. Comments starting with '#' are ignored.

- output_bed:

  Full file path and filename of the BED file to be written. Directory
  must exist.

- limit:

  Maximum allowed width (bp) for output intervals. Defaults to 2^29
  (536,870,912 bp). Must be larger than any single gene plus the shift
  distance.

- shift_by:

  Step size (bp) used to iteratively walk boundaries to the right when
  they violate gene clearance rules. Defaults to 1 bp.

- clearance:

  Minimum required clearance (bp) between split boundaries and gene
  intervals on both sides. Defaults to 1 bp.

## Value

Invisible character: path to the written BED file.

## Details

The function performs the following steps:

1.  Reads genomic regions and gene annotations from input files.

2.  Checks that the provided `limit` is feasible given gene positions
    and requested two-sided boundary clearance.

3.  For each region, splits it at the centromere if required and applies
    the size limit.

4.  Adjusts boundaries so each is at least `clearance` bp away from
    genes on both sides, moving by `shift_by` bp per iteration when
    needed.

5.  Iteratively subdivides intervals exceeding the limit until all
    intervals satisfy the constraint.

6.  Writes final split regions to BED file.

## Examples

``` r
  chromosomes <- system.file(
    "extdata",
    "AlgorithmToy.bed",
    package = "scShardSplitRef",
    mustWork = TRUE
  )

  genes <- system.file(
    "extdata",
    "AlgorithmToy.gtf",
    package = "scShardSplitRef",
    mustWork = TRUE
  )

  determine_split_regions(
    bed = chromosomes,
    gtf = genes,
    output_bed = file.path(tempdir(), "split_regions.bed"),
    limit = 2L^29L,
    shift_by = 1L,
    clearance = 1L
  )
#> Preparing split regions from
#> /home/runner/work/_temp/Library/scShardSplitRef/extdata/AlgorithmToy.bed and
#> /home/runner/work/_temp/Library/scShardSplitRef/extdata/AlgorithmToy.gtf into
#> /tmp/RtmpGo9qKK/split_regions.bed (`limit` = 536870912, `shift_by` = 1,
#> `clearance` = 1).
#> Built 10 split intervals across 7 chromosomes and wrote them to
#> /tmp/RtmpGo9qKK/split_regions.bed.
```
