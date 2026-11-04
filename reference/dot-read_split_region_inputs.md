# Read and prepare split-region inputs

Read and prepare split-region inputs

## Usage

``` r
.read_split_region_inputs(bed, gtf, limit, clearance)
```

## Arguments

- bed:

  BED input path.

- gtf:

  GTF input path.

- limit:

  Maximum interval width.

- clearance:

  Boundary clearance distance.

## Value

Named list containing `regions` and `merged_gene_ranges`.
