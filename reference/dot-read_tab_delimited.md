# Shared `read.table()` wrapper for TAB-delimited genomic files

All of the package's BED/GTF readers
([`.read_bed_file()`](https://IrinaVKuznetsova.github.io/scShardSplitRef/reference/dot-read_bed_file.md),
[`.read_gtf_file()`](https://IrinaVKuznetsova.github.io/scShardSplitRef/reference/dot-read_gtf_file.md)
in `determine_split_regions.R`;
[`.read_and_validate_bed()`](https://IrinaVKuznetsova.github.io/scShardSplitRef/reference/dot-read_and_validate_bed.md),
[`.read_and_validate_gtf()`](https://IrinaVKuznetsova.github.io/scShardSplitRef/reference/dot-read_and_validate_gtf.md)
here) call
[`utils::read.table()`](https://rdrr.io/r/utils/read.table.html) with
the same core options (`sep = "\t"`, `header = FALSE`, `quote = ""`,
`stringsAsFactors = FALSE`), differing only in whether they read from a
file path or already-read lines, whether comments are stripped, and
whether ragged rows are padded. This factors out just that repeated
call.

## Usage

``` r
.read_tab_delimited(file = NULL, text = NULL, comment_char = "#", fill = FALSE)
```

## Arguments

- file:

  Character: path to read from. Exactly one of `file`/`text` must be
  supplied.

- text:

  Character vector of already-read lines to parse via a text connection.
  Exactly one of `file`/`text` must be supplied.

- comment_char:

  Character: comment-line prefix, passed to
  [`read.table()`](https://rdrr.io/r/utils/read.table.html)'s
  `comment.char`. Defaults to `"#"`, matching
  [`read.table()`](https://rdrr.io/r/utils/read.table.html)'s own
  built-in default (several callers below rely on this default
  implicitly, exactly as they relied on
  [`read.table()`](https://rdrr.io/r/utils/read.table.html)'s default
  before this refactor). Use `""` to disable comment stripping.

- fill:

  Logical: passed to
  [`read.table()`](https://rdrr.io/r/utils/read.table.html)'s `fill`.
  `TRUE` pads ragged (short) rows with `NA` instead of erroring.
  Defaults to `FALSE`.

## Value

Data frame, exactly as returned by
[`utils::read.table()`](https://rdrr.io/r/utils/read.table.html).

## Details

Deliberately narrow in scope: each caller keeps its own file-exists
check, column-count check, column naming, and error messages exactly as
before, so this changes no external behavior, message, or signature –
only the duplicated
[`read.table()`](https://rdrr.io/r/utils/read.table.html) invocation
itself.
