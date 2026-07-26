#' Process a GTF file by assigning split-region sequence names
#'
#' Reads a BED-like "split regions" file and a GTF, assigns each GTF feature to
#' exactly one split-region interval, and replaces the GTF `Chr` field with a
#' split-region name of the form `Chr-RegionStart-RegionEnd`.
#'
#' A feature is considered inside a region when:
#' - `start >= RegionStart`
#' - `start <  RegionEnd`
#' - `end   <= RegionEnd`
#'
#' The function aborts if:
#' - any GTF feature matches more than one region (overlapping regions), or
#' - any GTF feature matches no region (regions don't fully cover the GTF).
#'
#' @section Output file:
#' The function will automatically name the file based on the input file and
#' prepends the following string onto the filename:
#' "B1_FINAL_MODIFIED_GTF_*genome_name*_*genome_version*.gtf", placing
#' the file in the directory specified in the `out_path`.
#'
#' @param split_regions_bed Character: File path to a BED-like file
#'   containing split intervals as produced by [determine_split_regions()].
#'   The first three columns must be: `Chr`, `RegionStart`, `RegionEnd`.
#' @inheritParams determine_split_regions
#' @param keep_attributes Character vector of attribute keys to keep (e.g.,
#'   `c("gene_id", "transcript_id", "gene_name")`). If provided, the GTF
#'   attribute column will be filtered to keep only these fields. Defaults to
#'   `NULL` (no filtering).
#' @param genome_name Genome identifier used in the output filename.
#' @param genome_version Genome version string used in the output filename.
#' @param out_path Directory where the output GTF will be written. Defaults to
#'   the current working directory (`"."`).
#'
#' @returns The function writes the filtered GTF file to disk.
#'
#' @examples
#' reg_file <- system.file(
#'   "extdata/AlgorithmToy_DetermineSplitReg.bed",
#'   package = "scShardSplitRef",
#'   mustWork = TRUE
#' )
#' gtf_file <- system.file(
#'   "extdata/AlgorithmToy.gtf",
#'   package = "scShardSplitRef",
#'   mustWork = TRUE
#' )
#'
#' process_gtf(
#'   split_regions_bed = reg_file,
#'   gtf = gtf_file,
#'   genome_name = "Example",
#'   genome_version = "v1.0",
#'   out_path = tempdir()
#' )
#'
#' @export
process_gtf <- function(
  split_regions_bed,
  gtf,
  genome_name,
  genome_version = NULL,
  keep_attributes = NULL,
  out_path = "."
) {
  cli::cli_alert_info("Validating inputs and reading files...")

  # Validate that genome_name and genome_version are provided
  if (missing(genome_name) || is.null(genome_name)) {
    cli::cli_abort(
      call = rlang::caller_env(),
      c(x = "Parameter {.arg genome_name} is required.")
    )
  }

  if (is.null(genome_version)) {
    cli::cli_abort(
      call = rlang::caller_env(),
      c(x = "Parameter {.arg genome_version} is required.")
    )
  }
  validated <- .validate_inputs(
    bed = split_regions_bed,
    gtf = gtf
  )
  regions <- validated$bed
  gtf <- validated$gtf

  gtf$unique_ID <- sprintf("ID%i", seq_len(nrow(gtf)))
  regions <- .prepare_split_regions(regions)

  if (.regions_overlap(regions)) {
    # Overlapping regions mean a feature really can match more than one
    # region, which the fast single-candidate lookup below cannot detect
    # (it assumes non-overlapping regions to avoid a full cross join). This
    # is the invalid-input path, so falling back to the original full
    # Chr-only join here is fine -- it only runs when we're about to abort
    # anyway, and it preserves the exact original ambiguous-match error
    # (message, duplicated IDs, etc).
    merged <- merge(
      gtf,
      regions[, c("Chr", "RegionStart", "RegionEnd", "NEW"), drop = FALSE],
      by = "Chr",
      all.x = TRUE,
      sort = FALSE
    )
    matched <- merged[.rows_within_region(merged), , drop = FALSE]
  } else {
    matched_gtf <- .match_features_to_regions(gtf, regions)
    matched <- matched_gtf[!is.na(matched_gtf$NEW), , drop = FALSE]
  }

  .abort_if_ambiguous_matches(matched)
  .abort_if_unmatched_features(matched, gtf)

  matched$start <- matched$start - matched$RegionStart + 1L
  matched$end <- matched$end - matched$RegionStart + 1L

  out <- .build_gtf_output(matched)

  # Guard against scientific notation: `start`/`end` are numeric (possibly
  # double, after the RegionStart-relative shift arithmetic) and must always
  # render as plain fixed-point integers in the written GTF.
  out[[4]] <- .format_coord(out[[4]])
  out[[5]] <- .format_coord(out[[5]])

  # Handle attribute filtering if parameters are provided
  if (!is.null(keep_attributes)) {
    filtered <- .filter_attributes(out$attribute, keep_attributes)
    .check_filtered_completeness(filtered, out$attribute)
    out$attribute <- filtered
  }

  # Write to file
  if (!dir.exists(out_path)) {
    cli::cli_abort(
      call = rlang::caller_env(),
      "Output directory does not exist: {.file {out_path}}"
    )
  }
  output_file <- .build_output_filename(out_path, genome_name, genome_version)
  .write_filtered_gtf(out, output_file)
  cli::cli_alert_info(c(
    i = "{.bold Finished writing processed GTF: {output_file}}"
  ))

  invisible(NULL)
}

#' Prepare split regions by adding `NEW` region names
#'
#' Adds a `NEW` column to the split-regions data frame, formatted as
#' `Chr-RegionStart-RegionEnd`.
#'
#' @param regions Data frame with columns `Chr`, `RegionStart`, `RegionEnd`.
#'
#' @returns The input data frame with an added `NEW` column.
#' @dev
.prepare_split_regions <- function(regions) {
  regions$NEW <- sprintf(
    "%s-%s-%s",
    regions$Chr,
    .format_coord(regions$RegionStart),
    .format_coord(regions$RegionEnd)
  )
  regions
}

#' Format a coordinate value as a plain fixed-point integer string
#'
#' `sprintf("%s", x)` and other implicit string coercions of a *numeric*
#' (double) vector go through `as.character()`, which follows
#' `getOption("scipen")`/`getOption("digits")` and will happily render large,
#' round genomic coordinates in scientific notation (e.g. `1e+08` instead of
#' `100000000`). Integers never do this, but `read.table()`-derived or
#' arithmetic-derived coordinate columns are not guaranteed to stay integer
#' (e.g. `matched$start - matched$RegionStart + 1` promotes to double).
#'
#' This helper forces fixed-point, non-scientific formatting regardless of
#' whether `x` is integer or double, and is safe for coordinates well within
#' double's exact-integer range (+/- 2^53), which comfortably covers any
#' genomic coordinate.
#'
#' @param x Numeric (integer or double) vector of coordinates. `NA` values
#'   are passed through as `NA` (not the string `"NA"`).
#'
#' @returns Character vector, same length as `x`, with each value rendered as
#'   a plain integer string (no scientific notation, no decimal point).
#'
#' @examples
#' \dontrun{
#'   .format_coord(1e8)      # "100000000", not "1e+08"
#'   .format_coord(100000000L) # "100000000"
#' }
#'
#' @dev
.format_coord <- function(x) {
  if (!is.numeric(x)) {
    cli::cli_abort(
      call = rlang::caller_env(),
      "{.arg x} must be numeric (integer or double), not {.cls {class(x)}}."
    )
  }

  out <- sprintf("%.0f", x)
  out[is.na(x)] <- NA_character_
  out
}

#' Check whether any split regions overlap within a chromosome
#'
#' The fast interval-assignment used to match GTF features to split regions
#' (see [.match_features_to_regions()]) assumes at most one split region can
#' contain any given position on a chromosome. This checks that invariant
#' upfront, cheaply (region count is typically small relative to feature
#' count), so [process_gtf()] can choose the fast path when it holds and
#' fall back to the original, slower-but-exhaustive matching when it
#' doesn't -- which is also the only case where the ambiguous-match error
#' actually fires, so falling back there costs nothing in practice.
#'
#' @param regions Data frame with columns `Chr`, `RegionStart`, `RegionEnd`.
#'
#' @returns Logical scalar: `TRUE` if any chromosome has overlapping region
#'   intervals, `FALSE` otherwise.
#'
#' @dev
.regions_overlap <- function(regions) {
  by_chr <- split(regions[, c("RegionStart", "RegionEnd")], regions$Chr)

  any(vapply(
    by_chr,
    FUN.VALUE = logical(1L),
    FUN = function(r) {
      if (nrow(r) <= 1L) {
        return(FALSE)
      }
      ord <- order(r$RegionStart)
      starts <- r$RegionStart[ord]
      ends <- r$RegionEnd[ord]
      any(starts[-1L] < ends[-length(ends)])
    }
  ))
}

#' Match GTF features to split regions
#'
#' Assigns each GTF feature to the single split-region interval that fully
#' contains it, per chromosome.
#'
#' @details
#' Rather than joining every feature against every split region on its
#' chromosome (a full cross join, which duplicates each feature once per
#' shard on that chromosome before being filtered back down -- expensive
#' for genome-scale GTFs split into many shards), this does a direct,
#' vectorized per-chromosome lookup: because split regions are asserted to
#' be non-overlapping (see [.regions_overlap()]), each feature
#' start position can belong to at most one region, found via a single
#' sorted-boundary [findInterval()] call instead of a full join. This turns
#' an O(features * regions_per_chromosome) join into O(features log
#' regions_per_chromosome).
#'
#' A feature is matched to a region only when it is fully contained
#' (`start >= RegionStart` and `end <= RegionEnd`); a feature whose start
#' falls in one region but whose end crosses into the next (or that falls
#' entirely before/after all regions on its chromosome) is left unmatched,
#' exactly as under the original cross-join-and-filter approach.
#'
#' @param gtf Data frame with columns `Chr`, `start`, `end`, `unique_ID`
#'   (plus other standard GTF columns).
#' @param regions Data frame with columns `Chr`, `RegionStart`, `RegionEnd`,
#'   `NEW`.
#'
#' @returns `gtf` with three additional columns, `NEW`, `RegionStart`,
#'   `RegionEnd`, populated for matched rows and `NA` for unmatched rows.
#'   Row order is identical to the input `gtf`.
#'
#' @dev
.match_features_to_regions <- function(gtf, regions) {
  gtf_by_chr <- split(seq_len(nrow(gtf)), gtf$Chr)
  region_by_chr <- split(seq_len(nrow(regions)), regions$Chr)

  new_col <- rep(NA_character_, nrow(gtf))
  region_start_col <- rep(NA_integer_, nrow(gtf))
  region_end_col <- rep(NA_integer_, nrow(gtf))

  for (chr in names(gtf_by_chr)) {
    r_idx <- region_by_chr[[chr]]
    if (is.null(r_idx)) {
      next # no split regions for this chromosome; features stay unmatched
    }

    g_idx <- gtf_by_chr[[chr]]

    ord <- order(regions$RegionStart[r_idx])
    r_idx <- r_idx[ord]
    starts <- regions$RegionStart[r_idx]

    # For each feature start, the candidate region is the one whose
    # RegionStart is the largest value <= the feature start (0 = before
    # the first region on this chromosome, i.e. no candidate).
    pos <- findInterval(gtf$start[g_idx], starts)
    has_candidate <- pos > 0L

    cand_g <- g_idx[has_candidate]
    cand_r <- r_idx[pos[has_candidate]]

    # Full containment also requires the feature end to fit inside the same
    # candidate region; otherwise it's left unmatched (e.g. it crosses a
    # shard boundary).
    fits <- gtf$end[cand_g] <= regions$RegionEnd[cand_r]

    new_col[cand_g[fits]] <- regions$NEW[cand_r[fits]]
    region_start_col[cand_g[fits]] <- regions$RegionStart[cand_r[fits]]
    region_end_col[cand_g[fits]] <- regions$RegionEnd[cand_r[fits]]
  }

  gtf$NEW <- new_col
  gtf$RegionStart <- region_start_col
  gtf$RegionEnd <- region_end_col
  gtf
}

#' Identify rows where a feature lies fully within a split region
#'
#' Applies the matching rule used by [process_gtf()]. A feature `start`, `end`
#' should be fully inside a region:
#'  `RegionStart`, `RegionEnd`:
#'  `start >= RegionStart`
#'  `end <= RegionEnd`
#'
#' Retained for reference/tests; the main pipeline now enforces this rule
#' directly inside [.match_features_to_regions()] without a full cross join.
#'
#' @param df Data frame containing `start`, `end`, `RegionStart`, `RegionEnd`.
#'
#' @returns Logical vector of length `nrow(df)`.
#' @dev
.rows_within_region <- function(df) {
  df$start >= df$RegionStart & df$end <= df$RegionEnd
}

#' Abort if any feature matches multiple split regions
#'
#' Detects duplicated `unique_ID` values in the matched table, which indicates
#' overlapping split regions (or otherwise ambiguous assignments).
#'
#' @param matched Data frame containing a `unique_ID` column.
#'
#' @returns Invisibly `TRUE` if ok; otherwise aborts.
#' @dev

.abort_if_ambiguous_matches <- function(matched) {
  if (nrow(matched) == 0L) {
    return(invisible(TRUE))
  }

  dup_ids <- unique(matched$unique_ID[duplicated(matched$unique_ID)])
  n_dup <- length(dup_ids)

  if (n_dup > 0L) {
    cli::cli_abort(
      call = rlang::caller_env(),
      c(
        "{n_dup} ambiguous split-region matches detected.",
        x = "{n_dup} feature IDs matched multiple intervals.",
        i = "Affected IDs:",
        i = paste0("* {.val ", dup_ids, "}"),
        i = "Ensure the split-regions BED has no overlapping intervals per chromosome."
      ),
      .envir = environment()
    )
  }

  invisible(TRUE)
}

#' Abort if any GTF features match no split region
#'
#' @param matched Data frame containing `unique_ID`.
#' @param gtf Original GTF data frame containing `unique_ID`, `Chr`, `feature`,
#'   `start`, and `end`.
#'
#' @returns Invisibly `TRUE` if ok; otherwise aborts.
#' @dev

.abort_if_unmatched_features <- function(matched, gtf) {
  missing_ids <- setdiff(gtf$unique_ID, matched$unique_ID)
  n_missing <- length(missing_ids)

  if (n_missing == 0L) {
    return(invisible(TRUE))
  }

  missing_rows <- gtf[
    gtf$unique_ID %in% missing_ids,
    c("Chr", "feature", "start", "end"),
    drop = FALSE
  ]

  preview_n <- min(5L, nrow(missing_rows))
  preview_items <- sprintf(
    "{.code %s:%s:%s-%s}",
    missing_rows$Chr[seq_len(preview_n)],
    missing_rows$feature[seq_len(preview_n)],
    missing_rows$start[seq_len(preview_n)],
    missing_rows$end[seq_len(preview_n)]
  )

  cli::cli_abort(
    call = rlang::caller_env(),
    c(
      "{n_missing} GTF features do not fall within any split-region interval.",
      x = "{n_missing} unmatched feature{?s} found.",
      i = "Example:",
      i = paste0("* ", preview_items),
      i = "Ensure split regions fully cover the annotation coordinates.",
      i = "Handle boundary-crossing features before calling {.fn process_gtf}."
    ),
    .envir = environment()
  )
}

#' Build the output GTF data frame
#'
#' Selects and renames columns so the output matches GTF layout, replacing
#' `Chr` with the split-region name stored in `NEW`.
#'
#' @param matched Data frame containing `NEW` and standard GTF columns.
#'
#' @returns A data frame with GTF columns where the first column is `Chr`.
#' @dev
.build_gtf_output <- function(matched) {
  out <- matched[, c(
    "NEW",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute"
  )]

  colnames(out)[1] <- "Chr"
  out
}

#' Filter GTF attributes
#'
#' Extract and keep only specified attributes from GTF attribute column
#'
#' @param attributes Character vector of GTF attribute strings
#' @param keep_attributes Character vector of attribute names to keep
#'
#' @returns Character vector of filtered attribute strings
#'
#' @details
#' Extracts key-value pairs from GTF attribute column and keeps only
#' specified attributes in the original order.
#'
#' @dev
.filter_attributes <- function(attributes, keep_attributes) {
  pattern <- paste0("^(?:", paste(keep_attributes, collapse = "|"), ")\\s+")

  str_list <- stringi::stri_split_regex(attributes, ";\\s*")

  vapply(
    str_list,
    FUN.VALUE = character(1L),
    FUN = function(pairs) {
      pairs <- pairs[nzchar(pairs)]
      key <- sub("\\s.*$", "", pairs)
      keep <- key %in% keep_attributes
      if (!any(keep)) {
        return("")
      }
      paste(pairs[keep], collapse = "; ")
    }
  )
}

#' Check filtered attributes completeness
#'
#' Verify that filtering didn't result in empty attributes
#'
#' @param filtered Character vector of filtered attributes
#' @param original Character vector of original attributes
#'
#' @returns Invisible TRUE if valid
#'
#' @dev
.check_filtered_completeness <- function(filtered, original) {
  is_empty_filtered <- !nzchar(filtered)
  is_nonempty_original <- nzchar(original)

  empty_idx <- which(is_empty_filtered & is_nonempty_original)

  if (length(empty_idx) > 0L) {
    cli::cli_warn(
      c(
        "{length(empty_idx)} lines had no matching attributes.",
        i = "Check your keep_attributes list."
      )
    )
  }

  invisible(TRUE)
}

#' Build output filename for filtered GTF
#'
#' @param out_path Character: output directory
#' @param genome_name Character: genome name
#' @param genome_version Character: genome version
#'
#' @returns Character: full file path
#'
#' @dev
.build_output_filename <- function(out_path, genome_name, genome_version) {
  filename <- sprintf(
    "B1_FINAL_MODIFIED_GTF_%s_%s.gtf",
    genome_name,
    genome_version
  )
  file.path(out_path, filename)
}

#' Write filtered GTF file
#'
#' @param result Data frame with filtered GTF data
#' @param output_file Character: output file path
#'
#' @returns Called for its side-effects of writing a file to disk.
#'
#' @dev
.write_filtered_gtf <- function(result, output_file) {
  invisible(utils::write.table(
    result,
    file = output_file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  ))
}
