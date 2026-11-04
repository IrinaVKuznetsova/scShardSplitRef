# Check that files exist

Validates that both BED and GTF file paths point to existing files.

## Usage

``` r
.check_files_exist(bed_path, gtf_path)
```

## Arguments

- bed_path:

  Character: BED file path.

- gtf_path:

  Character: GTF file path.

## Value

Invisible TRUE if both files exist. Raises error otherwise.

## Examples

``` r
if (FALSE) { # \dontrun{
  .check_files_exist("regions.bed", "genes.gtf")
} # }
```
