#' Determine split regions for reference sharding
#'
#' Builds BED-like split intervals from genomic regions and GTF gene
#' annotations. If the start position in a genomic region is not 0, it is
#' assumed to mark the start of a centromere and the region is split around
#' that point. Proposed boundaries are shifted rightward when they fall too
#' close to gene intervals, and oversized chunks are repeatedly re-split until
#' all segments satisfy the requested size limit.
#'
#' @param bed Path to a 3-column BED-like file (chr, start, end).
#'   Columns should be TAB-delimited with chromosome name, region start
#'   (0-based), and region end coordinates.
#' @param gtf Path to a GTF file containing gene annotations.
#'   Must contain standard 9 GTF columns. Comments starting with '#' are
#'   ignored.
#' @param output_bed Path to the BED file to be written. Directory must exist.
#' @param limit Maximum allowed width (bp) for output intervals.
#'   Defaults to 2^29 (536,870,912 bp). Must be larger than any single gene
#'   plus the shift distance.
#' @param shift_by Step size (bp) used to iteratively walk boundaries to the
#'   right when they violate gene clearance rules. Defaults to 1 bp.
#' @param clearance Minimum required clearance (bp) between split boundaries
#'   and gene intervals on both sides. Defaults to 1 bp.
#'
#' @returns Invisible character: path to the written BED file.
#'
#' @details
#' The function performs the following steps:
#' 1. Reads genomic regions and gene annotations from input files.
#' 2. Checks that the provided `limit` is feasible given gene positions and
#'    requested two-sided boundary clearance.
#' 3. For each region, splits it at the centromere if required and applies
#'    the size limit.
#' 4. Adjusts boundaries so each is at least `clearance` bp away from genes on
#'    both sides, moving by `shift_by` bp per iteration when needed.
#' 5. Iteratively subdivides intervals exceeding the limit until all intervals
#'    satisfy the constraint.
#' 6. Writes final split regions to BED file.
#'
#' @examples
#'   chromosomes <- system.file(
#'     "extdata",
#'     "AlgorithmToy.bed",
#'     package = "scShardSplitRef",
#'     mustWork = TRUE
#'   )
#'
#'   genes <- system.file(
#'     "extdata",
#'     "AlgorithmToy.gtf",
#'     package = "scShardSplitRef",
#'     mustWork = TRUE
#'   )
#'
#'   determine_split_regions(
#'     bed = chromosomes,
#'     gtf = genes,
#'     output_bed = file.path(tempdir(), "split_regions.bed"),
#'     limit = 2L^29L,
#'     shift_by = 1L,
#'     clearance = 1L
#'   )
#' @export
determine_split_regions <- function(
  bed,
  gtf,
  output_bed,
  limit = 2L^29L,
  shift_by = 1L,
  clearance = 1L
) {
  args <- .validate_split_region_args(
    shift_by = shift_by,
    clearance = clearance
  )

  shift_by <- args$shift_by
  clearance <- args$clearance

  .inform_split_region_start(
    bed = bed,
    gtf = gtf,
    output_bed = output_bed,
    limit = limit,
    shift_by = shift_by,
    clearance = clearance
  )

  inputs <- .read_split_region_inputs(
    bed = bed,
    gtf = gtf,
    limit = limit,
    clearance = clearance
  )

  out <- .build_split_region_output(
    regions = inputs$regions,
    merged_gene_ranges = inputs$merged_gene_ranges,
    limit = limit,
    shift_by = shift_by,
    clearance = clearance
  )

  .write_bed_file(out, output_bed)

  .inform_split_region_done(
    out = out,
    output_bed = output_bed
  )

  invisible(output_bed)
}

#' Validate scalar arguments for split-region generation
#'
#' @param shift_by Step size for boundary shifting.
#' @param clearance Minimum required boundary clearance around genes.
#'
#' @returns Named list containing integer `shift_by` and `clearance`.
#'
#' @dev
.validate_split_region_args <- function(shift_by, clearance) {
  list(
    shift_by = .validate_scalar_integerish(
      x = shift_by,
      arg = "shift_by",
      min_value = 1L
    ),
    clearance = .validate_scalar_integerish(
      x = clearance,
      arg = "clearance",
      min_value = 0L
    )
  )
}

#' Validate and coerce a scalar integer-like numeric argument
#'
#' @param x Value to validate.
#' @param arg Argument name used in error messages.
#' @param min_value Minimum permitted value.
#'
#' @returns Integer scalar.
#'
#' @dev
.validate_scalar_integerish <- function(x, arg, min_value) {
  valid <- is.numeric(x) &&
    length(x) == 1L &&
    !is.na(x) &&
    x >= min_value

  if (!valid) {
    cli::cli_abort(
      call = rlang::caller_env(),
      "{.arg {arg}} must be a single numeric value >= {min_value}."
    )
  }

  as.integer(x)
}

#' Inform user that split-region generation has started
#'
#' @param bed BED input path.
#' @param gtf GTF input path.
#' @param output_bed Output BED path.
#' @param limit Maximum interval width.
#' @param shift_by Boundary shift step size.
#' @param clearance Boundary clearance distance.
#'
#' @returns Invisible `TRUE`.
#'
#' @dev
.inform_split_region_start <- function(
  bed,
  gtf,
  output_bed,
  limit,
  shift_by,
  clearance
) {
  cli::cli_inform(
    paste(
      "Preparing split regions from {.file {bed}} and {.file {gtf}}",
      "into {.file {output_bed}}",
      "({.var limit} = {limit},",
      "{.var shift_by} = {shift_by},",
      "{.var clearance} = {clearance})."
    )
  )

  invisible(TRUE)
}

#' Read and prepare split-region inputs
#'
#' @param bed BED input path.
#' @param gtf GTF input path.
#' @param limit Maximum interval width.
#' @param clearance Boundary clearance distance.
#'
#' @returns Named list containing `regions` and `merged_gene_ranges`.
#'
#' @dev
.read_split_region_inputs <- function(bed, gtf, limit, clearance) {
  regions <- .read_bed_file(bed, 3L)
  genes <- .read_gene_annotations(gtf)

  .assert_feasible_limit(genes, limit, clearance)

  list(
    regions = regions,
    merged_gene_ranges = .merge_all_gene_ranges(genes)
  )
}

#' Read gene annotations from a GTF file
#'
#' @param gtf GTF input path.
#'
#' @returns Data frame containing gene feature rows only.
#'
#' @dev
.read_gene_annotations <- function(gtf) {
  genes <- .read_gtf_file(gtf)
  genes[genes$feature == "gene", , drop = FALSE]
}

#' Build final split-region output data frame
#'
#' @param regions BED-like regions data frame.
#' @param merged_gene_ranges Named list of merged gene ranges by chromosome.
#' @param limit Maximum interval width.
#' @param shift_by Boundary shift step size.
#' @param clearance Boundary clearance distance.
#'
#' @returns Data frame of final split intervals.
#'
#' @dev
.build_split_region_output <- function(
  regions,
  merged_gene_ranges,
  limit,
  shift_by,
  clearance
) {
  split_regions <- lapply(
    seq_len(nrow(regions)),
    function(i) {
      .process_single_region(
        region = regions[i, ],
        merged_gene_ranges = merged_gene_ranges,
        limit = limit,
        shift_by = shift_by,
        clearance = clearance
      )
    }
  )

  do.call(rbind, split_regions)
}

#' Inform user that split-region generation has finished
#'
#' @param out Final split-region data frame.
#' @param output_bed Output BED path.
#'
#' @returns Invisible `TRUE`.
#'
#' @dev
.inform_split_region_done <- function(out, output_bed) {
  cli::cli_inform(
    "{.strong Built {nrow(out)} split intervals across {length(unique(out$chr))} chromosome{?s} and wrote them to {.file {output_bed}}.}"
  )

  invisible(TRUE)
}


#' Read BED file with validation
#'
#' Reads a BED-like file and ensures it has the required minimum number of
#' columns. The first three columns are named 'chr', 'start', 'end'.
#'
#' @param path Character: file path to BED file.
#' @param min_cols Integer: minimum required number of columns. Defaults to 3L.
#'
#' @return Data frame with columns 'chr', 'start', 'end' plus any additional
#'   columns from the input file.
#'
#' @details
#' Raises an error if the file does not exist or has fewer than `min_cols`
#' columns. Assumes TAB-delimited format with no header.
#'
#' @examples
#' \dontrun{
#'   bed_df <- .read_bed_file("regions.bed", min_cols = 3L)
#' }
#'
#' @dev
.read_bed_file <- function(path, min_cols = 3L) {
  if (!file.exists(path)) {
    cli::cli_abort(
      call = rlang::caller_env(),
      "BED file not found: {.file {path}}"
    )
  }

  d <- .read_tab_delimited(file = path, fill = TRUE)

  ncol_d <- ncol(d)
  if (ncol_d < min_cols) {
    cli::cli_abort(
      call = rlang::caller_env(),
      "BED file must have at least {min_cols} columns, found {ncol_d}."
    )
  }

  colnames(d)[seq_len(3L)] <- c("chr", "start", "end")
  d$start <- .coerce_integer_coord(d$start, "start")
  d$end <- .coerce_integer_coord(d$end, "end")
  return(d)
}

#' Read GTF file with validation
#'
#' Reads a GTF file and ensures it has the required 9 columns. Comment lines
#' starting with '#' are automatically skipped.
#'
#' @param path Character: file path to GTF file.
#'
#' @return Data frame with standard 9 GTF columns: chr, source, feature,
#'   start, end, score, strand, frame, attribute.
#'
#' @details
#' Raises an error if the file does not exist or has fewer than 9 columns.
#' Assumes TAB-delimited format with no header.
#'
#' @examples
#' \dontrun{
#'   gtf_df <- .read_gtf_file("annotations.gtf")
#' }
#'
#' @dev
.read_gtf_file <- function(path) {
  if (!file.exists(path)) {
    cli::cli_abort(
      call = rlang::caller_env(),
      "GTF file not found: {.file {path}}"
    )
  }

  d <- .read_tab_delimited(file = path, comment_char = "#", fill = TRUE)

  ncol_d <- ncol(d)

  if (ncol_d < 9L) {
    cli::cli_abort(
      call = rlang::caller_env(),
      "GTF file must have at least 9 columns, found {ncol_d}."
    )
  }

  colnames(d) <- c(
    "chr",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute"
  )
  d$start <- .coerce_integer_coord(d$start, "start")
  d$end <- .coerce_integer_coord(d$end, "end")
  return(d)
}

#' Write BED file
#'
#' Writes a data frame to a TAB-delimited BED format file without quotes
#' or row names.
#'
#' @param df Data frame with 'chr', 'start', 'end' columns.
#' @param path Character: output file path.
#'
#' @returns Called for its side effects, writes a file to local disk.
#'
#' @examples
#' \dontrun{
#'   .write_bed_file(split_regions_df, "output.bed")
#' }
#'
#' @dev
.write_bed_file <- function(df, path) {
  invisible(
    utils::write.table(
      df,
      file = path,
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t"
    )
  )
}

#' Process a single genomic region
#'
#' Splits a single genomic region according to the size limit and adjusts
#' boundaries to avoid genes. If the region has a non-zero start (centromere),
#' it is split into two parts: before and after the centromere.
#'
#' @param region Data frame row with columns 'chr', 'start', 'end'.
#' @param merged_gene_ranges Named list, keyed by chromosome, of data frames
#'   with columns 'start', 'end' -- gene ranges already merged via
#'   [.merge_all_gene_ranges()]. Passed in rather than recomputed here so the
#'   merge work happens once per chromosome regardless of how many regions
#'   share it.
#' @param limit Integer: maximum allowed interval width in base pairs.
#' @param shift_by Integer: step size (bp) used to walk boundaries rightward
#'  when clearance is violated.
#' @param clearance Integer: minimum required clearance (bp) from nearby genes
#'  on both sides.
#'
#' @return Data frame with columns 'chr', 'start', 'end' representing the
#'   processed split intervals for this region.
#'
#' @details
#' If the region start is > 0, the function assumes this marks a centromere
#' and splits the region into [0, start) and [start, end). Each part is
#' independently split according to the size limit. Boundaries are then
#' adjusted to avoid genes in the region.
#'
#' @examples
#' \dontrun{
#'   region <- data.frame(chr = "chr1", start = 100000, end = 200000)
#'   genes_subset <- genes[genes$chr == "chr1", ]
#'   result <- .process_single_region(region, genes_subset,
#'                                    limit = 50000,
#'                                    shift_by = 1,
#'                                    clearance = 1)
#' }
#'
#' @dev
.process_single_region <- function(
  region,
  merged_gene_ranges,
  limit,
  shift_by,
  clearance = 1L
) {
  chr <- region$chr
  region_start <- as.integer(region$start)
  region_end <- as.integer(region$end)

  .validate_bed_coordinates(region_start, region_end, chr)

  # Split around centromere if start > 0
  if (region_start > 0L) {
    pre_df <- .split_region_core(chr, 0L, region_start, limit)
    post_df <- .split_region_core(chr, region_start, region_end, limit)
    split_df <- rbind(pre_df, post_df)
  } else {
    split_df <- .split_region_core(chr, 0L, region_end, limit)
  }

  # Fix boundaries that fall inside genes, using the already-merged ranges
  # for this chromosome (empty data frame if this chromosome has no genes).
  chr_gene_ranges <- merged_gene_ranges[[chr]]
  if (is.null(chr_gene_ranges)) {
    chr_gene_ranges <- data.frame(start = integer(0L), end = integer(0L))
  }
  .fix_split_boundaries(split_df, chr_gene_ranges, limit, shift_by, clearance)
}

#' Validate BED coordinate format
#'
#' Checks that start and end coordinates satisfy BED format requirements:
#' start >= 0 and start < end.
#'
#' @param start Integer: region start position (0-based).
#' @param end Integer: region end position (exclusive).
#' @param chr Character: chromosome name for error messages.
#'
#' @return Invisible TRUE if valid. Raises error if invalid.
#'
#' @examples
#' \dontrun{
#'   .validate_bed_coordinates(0L, 100L, "chr1")
#' }
#'
#' @dev
#'

.validate_bed_coordinates <- function(start, end, chr) {
  # Extract scalar values if they're from a data frame row
  if (length(start) > 1L) {
    start <- start[1L]
  }
  if (length(end) > 1L) {
    end <- end[1L]
  }
  if (length(chr) > 1L) {
    chr <- chr[1L]
  }

  if (anyNA(start) || anyNA(end) || anyNA(chr)) {
    cli::cli_abort(
      call = rlang::caller_env(),
      "Missing values in BED region coordinates."
    )
  }

  if (start < 0L) {
    cli::cli_abort(
      call = rlang::caller_env(),
      "Region start must be >= 0 for BED coordinates ({.var chr}={chr})."
    )
  }

  if (start >= end) {
    cli::cli_abort(
      call = rlang::caller_env(),
      "Invalid BED region: {.var start} ({start}) must be < {.var end}
      ({end}) for {.var chr}={chr}."
    )
  }

  return(invisible(TRUE))
}

#' Split a single coordinate interval
#'
#' Splits a coordinate interval [start, end) into chunks of at most `limit`
#' base pairs. Returns a data frame with the resulting intervals.
#'
#' @param chr Character: chromosome name.
#' @param start Integer: interval start (inclusive).
#' @param end Integer: interval end (exclusive).
#' @param limit Integer: maximum width per interval in base pairs.
#'
#' @return Data frame with columns 'chr', 'start', 'end' representing the
#'   split intervals. If the interval width <= limit, a single row is returned.
#'
#' @examples
#' \dontrun{
#'   result <- .split_region_core("chr1", 0L, 1000L, 250L)
#' }
#'
#' @dev
.split_region_core <- function(chr, start, end, limit) {
  width <- end - start
  boundaries <- if (width <= limit) {
    integer(0)
  } else {
    .build_boundaries(start, end, limit)
  }
  .boundaries_to_df(chr, start, end, boundaries)
}

#' Generate split boundaries for an interval
#'
#' Calculates interior boundary points that divide a coordinate interval
#' [start, end) into chunks of approximately equal size, each at most `limit`
#' base pairs wide.
#'
#' @param start Integer: interval start (inclusive).
#' @param end Integer: interval end (exclusive).
#' @param limit Integer: maximum width per chunk in base pairs.
#'
#' @return Integer vector of boundary positions. Empty if interval width <=
#'  limit.
#'
#' @details
#' The function calculates the number of chunks needed, then generates
#' evenly-spaced boundaries. Interior boundaries are rounded to the nearest
#' integer. The first and last boundaries are excluded (they are the interval
#' endpoints).
#'
#' @examples
#' \dontrun{
#'   boundaries <- .build_boundaries(0L, 1000L, 250L)
#'   # Returns something like: 250, 500, 750
#' }
#'
#' @dev
.build_boundaries <- function(start, end, limit) {
  width <- end - start
  if (width <= limit) {
    return(integer(0))
  }

  n_chunks <- ceiling(width / limit)
  step <- width / n_chunks

  as.integer(round(seq(
    start + step,
    end - step,
    by = step
  )))
}

#' Convert boundary points to interval data frame
#'
#' Transforms a vector of boundary positions into a data frame representing
#' intervals [start, end) between consecutive boundaries.
#'
#' @param chr Character: chromosome name for all rows.
#' @param start Integer: the minimum start coordinate (first boundary).
#' @param end Integer: the maximum end coordinate (last boundary).
#' @param boundaries Integer vector: interior boundary points (may be empty).
#'
#' @return Data frame with columns 'chr', 'start', 'end'. Each row represents
#'   an interval between consecutive boundary points.
#'
#' @examples
#' \dontrun{
#'   df <- .boundaries_to_df("chr1", 0L, 100L, c(30L, 60L))
#'   # Returns:
#'   #   chr start end
#'   # 1 chr1    0  30
#'   # 2 chr1   30  60
#'   # 3 chr1   60 100
#' }
#'
#' @dev
.boundaries_to_df <- function(chr, start, end, boundaries) {
  points <- c(start, boundaries, end)
  data.frame(
    chr = chr,
    start = points[-length(points)],
    end = points[-1L],
    stringsAsFactors = FALSE
  )
}

# ============================================================================
# Boundary Fixing (Reduced Complexity)
# ============================================================================

#' Fix split boundaries to avoid genes
#'
#' Adjusts interval boundaries to ensure none fall within gene annotations.
#' Iteratively shifts boundaries past genes and subdivides oversized intervals
#' until all constraints are satisfied.
#'
#' @param split_df Data frame with columns 'chr', 'start', 'end' representing
#'   initial split intervals.
#' @param gene_ranges Data frame with already-merged, non-overlapping gene
#'   ranges for the same chromosome (columns 'start', 'end'), as produced by
#'   [.merge_all_gene_ranges()]. Passed in pre-merged so this function does
#'   not redo the sort-and-merge work for every region on the chromosome.
#' @param limit Integer: maximum allowed interval width in base pairs.
#' @param shift_by Integer: step size (bp) used to walk boundaries rightward
#'   when clearance is violated. Defaults to 1L.
#' @param clearance Integer: minimum required clearance (bp) from nearby genes
#'   on both sides. Defaults to 1L.
#'
#' @return Data frame with same structure as `split_df` but with adjusted
#'   boundaries that satisfy the size and gene-avoidance constraints.
#'
#' @details
#' The algorithm iteratively:
#' 1. Shifts boundaries rightward in `shift_by` bp increments until they
#'    satisfy the `clearance` gene distance on both sides
#' 2. Subdivides intervals exceeding the limit
#' 3. Shifts again after subdivision
#' 4. Repeats until convergence
#'
#' Raises an error if constraints cannot be satisfied.
#'
#' @examples
#' \dontrun{
#'   fixed_df <- fix_split_boundaries(
#'     split_df, gene_ranges, limit = 50000, shift_by = 1, clearance = 1
#'   )
#' }
#'
#' @dev
.fix_split_boundaries <- function(
  split_df,
  gene_ranges,
  limit,
  shift_by = 1L,
  clearance = 1L
) {
  if (nrow(split_df) <= 1L || nrow(gene_ranges) == 0L) {
    return(split_df)
  }

  chr <- split_df$chr[1L]
  chr_end <- as.integer(max(split_df$end))
  boundaries <- as.integer(split_df$end[-nrow(split_df)])

  boundaries <- .refine_boundaries(
    boundaries,
    gene_ranges,
    chr_end,
    limit,
    shift_by,
    clearance
  )

  .validate_final_intervals(boundaries, chr_end, limit, chr)
  .boundaries_to_df(chr, 0L, chr_end, boundaries)
}

#' Iteratively refine boundaries
#'
#' Core algorithm that repeatedly adjusts boundaries by shifting past genes
#' and subdividing oversized intervals until convergence.
#'
#' @param boundaries Integer vector: initial boundary positions.
#' @param gene_ranges Data frame with merged gene ranges (columns 'start',
#'  'end').
#' @param chr_end Integer: chromosome end coordinate.
#' @param limit Integer: maximum allowed interval width.
#' @param shift_by Integer: step size (bp) for rightward boundary updates.
#' @param clearance Integer: minimum required clearance (bp) from nearby genes
#'   on both sides.
#' @param max_iter Integer: maximum number of iterations. Defaults to 100L.
#'
#' @return Integer vector: refined boundary positions after convergence or
#'   max iterations.
#'
#' @details
#' Each iteration:
#' 1. Shifts all boundaries past overlapping genes
#' 2. Collects extra boundaries from oversized intervals
#' 3. Deduplicates and filters boundaries
#' 4. Shifts again
#' 5. Checks for convergence (identical to previous iteration)
#'
#' @examples
#' \dontrun{
#'   refined <- .refine_boundaries(
#'     c(100L, 200L), gene_ranges, 300L, 50L, 1L, 1L
#'   )
#' }
#'
#' @dev
.refine_boundaries <- function(
  boundaries,
  gene_ranges,
  chr_end,
  limit,
  shift_by,
  clearance = 1L,
  max_iter = 100L
) {
  for (iter in seq_len(max_iter)) {
    old_boundaries <- boundaries

    # Shift past genes
    boundaries <- .shift_boundaries_past_genes(
      boundaries,
      gene_ranges,
      chr_end,
      shift_by,
      clearance
    )

    # Subdivide oversized intervals
    points <- c(0L, boundaries, chr_end)
    extra_boundaries <- .subdivide_oversized(points, limit)
    boundaries <- sort(unique(c(boundaries, extra_boundaries)))
    boundaries <- boundaries[boundaries > 0L & boundaries < chr_end]

    # Shift again after subdivision
    boundaries <- .shift_boundaries_past_genes(
      boundaries,
      gene_ranges,
      chr_end,
      shift_by,
      clearance
    )

    # Check convergence
    if (identical(boundaries, old_boundaries)) {
      return(boundaries)
    }
  }

  ## Max iterations exceeded
  stop(
    sprintf(
      ".refine_boundaries() failed to converge after %d iterations",
      max_iter
    ),
    call. = FALSE
  )
}

#' Validate final split intervals
#'
#' Checks that all final intervals satisfy the size limit constraint.
#' Raises an error if any interval exceeds the limit.
#'
#' @param boundaries Integer vector: boundary positions (interior points only).
#' @param chr_end Integer: chromosome end coordinate.
#' @param limit Integer: maximum allowed interval width.
#' @param chr Character: chromosome name for error messages.
#'
#' @return Invisible `TRUE` if all intervals satisfy the limit. Raises error
#'  otherwise.
#'
#' @examples
#' \dontrun{
#'   .validate_final_intervals(c(100L, 200L), 300L, 50L, "chr1")
#' }
#'
#' @dev
.validate_final_intervals <- function(boundaries, chr_end, limit, chr) {
  final <- c(0L, boundaries, chr_end)
  widths <- diff(final)

  if (any(widths > limit)) {
    bad_idx <- which(widths > limit)
    cli::cli_abort(
      call = rlang::caller_env(),
      "Unable to split chromosome {chr} to satisfy limit={limit}.
       Oversized intervals at positions: {bad_idx}."
    )
  }

  invisible(TRUE)
}

#' Subdivide oversized intervals
#'
#' Identifies intervals between consecutive points that exceed the size limit
#' and generates interior boundaries to split them. Returns all extra
#' boundaries needed.
#'
#' @param points Integer vector: all boundary points (including endpoints 0
#'   and chr_end) in sorted order.
#' @param limit Integer: maximum allowed interval width.
#'
#' @return Integer vector of extra interior boundaries needed. Empty if no
#'   intervals exceed the limit.
#'
#' @examples
#' \dontrun{
#'   extra <- .subdivide_oversized(c(0L, 100L, 500L), limit = 100L)
#'   # Interval [100, 500] is 400bp, needs splitting
#' }
#'
#' @dev
.subdivide_oversized <- function(points, limit) {
  widths <- diff(points)
  idx <- which(widths > limit)

  if (length(idx) == 0L) {
    return(integer(0L))
  }

  n <- ceiling(widths[idx] / limit)

  extra <- unlist(
    Map(
      function(i, k) {
        width <- widths[i]
        step <- width / k
        as.integer(round(seq(
          points[i] + step,
          points[i + 1] - step,
          by = step
        )))
      },
      idx,
      n
    ),
    use.names = FALSE
  )

  extra
}

#' Shift boundaries past genes
#'
#' Adjusts a vector of boundaries that fall too close to gene intervals by
#' shifting them rightward. Each boundary is independently shifted until it
#' satisfies the requested clearance around genes. Results are deduplicated and
#' filtered.
#'
#' @param boundaries Integer vector: boundary positions to adjust.
#' @param gene_ranges Data frame with gene ranges (columns 'start', 'end').
#' @param chr_end Integer: chromosome end coordinate (upper bound).
#' @param shift_by Integer: step size (bp) for rightward boundary updates.
#' @param clearance Integer: minimum required clearance (bp) from nearby genes
#'   on both sides.
#'   Defaults to 1L.
#'
#' @return Integer vector: adjusted, deduplicated, and filtered boundaries.
#'   Values outside (0, chr_end) are removed.
#'
#' @examples
#' \dontrun{
#'   shifted <- .shift_boundaries_past_genes(
#'     c(60L, 150L), gene_ranges, 300L, shift_by = 1L, clearance = 1L
#'   )
#' }
#'
#' @dev
.shift_boundaries_past_genes <- function(
  boundaries,
  gene_ranges,
  chr_end,
  shift_by = 1L,
  clearance = 1L
) {
  if (length(boundaries) == 0L || nrow(gene_ranges) == 0L) {
    return(boundaries)
  }

  shifted <- vapply(
    boundaries,
    FUN = .shift_single_boundary,
    FUN.VALUE = integer(1L),
    gene_ranges = gene_ranges,
    chr_end = chr_end,
    shift_by = shift_by,
    clearance = clearance
  )

  shifted <- sort(unique(shifted))
  shifted[shifted > 0L & shifted < chr_end]
}

#' Shift a single boundary past genes
#'
#' Adjusts one boundary position that falls inside or too close to gene
#' intervals by shifting it rightward until it satisfies the requested
#' clearance around genes. Each shift jumps directly to the nearest
#' shift_by-reachable position that clears the currently-violating range(s),
#' rather than walking forward one shift_by step at a time, so the result is
#' identical to a naive per-step walk but computed in O(number of blocking
#' gene clusters) instead of O(distance / shift_by).
#'
#' Candidate gene ranges are also found via two binary searches
#' (`findInterval()` against `gene_ranges$start` and `gene_ranges$end`,
#' both sorted ascending since `gene_ranges` is pre-merged and
#' non-overlapping) rather than a linear scan of every range on the
#' chromosome, turning an O(n_genes) check into O(log n_genes + k), where
#' k is the (typically small) number of ranges actually within `clearance`
#' of `boundary`.
#'
#' @param boundary Integer: initial boundary position.
#' @param gene_ranges Data frame with gene ranges (columns 'start', 'end'),
#'   already merged and sorted ascending by start (and therefore also by
#'   end, since ranges are non-overlapping).
#' @param chr_end Integer: chromosome end coordinate (upper bound).
#' @param shift_by Integer: step size (bp) for rightward boundary updates.
#' @param clearance Integer: minimum required clearance (bp) from nearby genes
#'   on both sides.
#'
#' @return Integer: adjusted boundary position. May be the original position
#'   if it did not overlap any gene, or identical to its start position if
#'   the shift would exceed chr_end.
#'
#' @examples
#' \dontrun{
#'   adjusted <- .shift_single_boundary(60L, gene_ranges, 300L, 1L, 1L)
#' }
#'
#' @dev
.shift_single_boundary <- function(
  boundary,
  gene_ranges,
  chr_end,
  shift_by,
  clearance
) {
  if (shift_by < 1L) {
    cli::cli_abort(
      call = rlang::caller_env(),
      "{.arg shift_by} must be >= 1 for iterative boundary updates."
    )
  }

  step <- as.integer(shift_by)
  clearance <- as.integer(clearance)
  n_ranges <- nrow(gene_ranges)

  repeat {
    # Ranges that could have `boundary` within their clearance zone are
    # exactly those with start < boundary + clearance AND end > boundary -
    # clearance. Both bounds are found via binary search instead of
    # scanning every range:
    #  - hi: last range index with start < boundary + clearance
    #  - lo: first range index with end > boundary - clearance
    hi <- findInterval(boundary + clearance - 1L, gene_ranges$start)
    lo <- findInterval(boundary - clearance, gene_ranges$end) + 1L

    if (lo > hi || hi < 1L) {
      break
    }

    candidates <- lo:min(hi, n_ranges)

    inside_clearance <-
      boundary > (gene_ranges$start[candidates] - clearance) &
      boundary < (gene_ranges$end[candidates] + clearance)

    if (!any(inside_clearance)) {
      break
    }

    # Jump directly to the first shift_by-reachable position that clears
    # every currently-violating gene range, instead of walking forward one
    # shift_by step at a time. With the default shift_by = 1L, a boundary
    # landing inside e.g. a megabase-scale gene previously required on the
    # order of a million single-bp iterations (each rescanning gene_ranges)
    # to clear it; this reaches the same final position in O(1), and the
    # outer repeat loop still runs again in case the new position lands
    # inside a different gene's clearance zone.
    threshold <- max(gene_ranges$end[candidates[inside_clearance]]) +
      clearance
    n_steps <- max(1L, as.integer(ceiling((threshold - boundary) / step)))
    updated <- boundary + n_steps * step

    if (updated >= chr_end || updated <= boundary) {
      break
    }

    boundary <- updated
  }

  return(boundary)
}

# ============================================================================
# Gene Range Processing
# ============================================================================

#' Merge gene ranges for every chromosome, once
#'
#' Groups gene annotations by chromosome and merges overlapping ranges
#' within each group, producing a lookup used by both the upfront
#' feasibility check and per-region boundary fixing -- so the merge work
#' happens exactly once per chromosome, rather than being recomputed by
#' every downstream consumer.
#'
#' @param genes Data frame with gene annotations (columns 'chr', 'start',
#'  'end').
#'
#' @return Named list, keyed by chromosome, of data frames with columns
#'   'start', 'end' representing merged, non-overlapping gene ranges. Empty
#'   list if `genes` has 0 rows.
#'
#' @examples
#' \dontrun{
#'   merged_gene_ranges <- .merge_all_gene_ranges(genes)
#'   merged_gene_ranges[["chr1"]]
#' }
#'
#' @dev
.merge_all_gene_ranges <- function(genes) {
  if (nrow(genes) == 0L) {
    return(list())
  }

  lapply(
    split(genes[, c("start", "end")], genes$chr),
    .merge_overlapping_gene_ranges
  )
}

#' Assert that the size limit is feasible
#'
#' Checks that the requested size `limit` is larger than any single gene
#' with two-sided boundary clearance. Raises an informative error if not
#' feasible.
#'
#' @param genes Data frame with gene annotations (columns 'chr', 'start',
#'  'end').
#' @param limit Integer: requested maximum interval width.
#' @param clearance Integer: minimum required clearance (bp) from nearby genes
#'   on both sides.
#'
#' @return Invisible TRUE if feasible. Raises error with remediation advice
#'   if not.
#'
#' @details
#' Merges overlapping genes per chromosome via [.merge_all_gene_ranges()],
#' then checks if any merged range plus two-sided clearance (2 * `clearance`)
#' exceeds the limit. Reports the worst (largest) blocking interval and
#' suggests an appropriate minimum limit value.
#'
#' @examples
#' \dontrun{
#'   .assert_feasible_limit(genes, limit = 50000, clearance = 1)
#' }
#'
#' @dev
.assert_feasible_limit <- function(genes, limit, clearance) {
  if (nrow(genes) == 0L) {
    return(invisible(TRUE))
  }

  merged_gene_ranges <- .merge_all_gene_ranges(genes)

  merged <- do.call(
    rbind,
    Map(
      function(chr, m) {
        m$chr <- chr
        m
      },
      names(merged_gene_ranges),
      merged_gene_ranges
    )
  )

  blocked_width <- merged$end - merged$start + (2L * as.integer(clearance))
  max_blocked <- max(blocked_width, na.rm = TRUE)

  if (limit < max_blocked) {
    worst <- merged[which.max(blocked_width), , drop = FALSE]
    cli::cli_abort(
      call = rlang::caller_env(),
      c(
        x = "Limit {limit} is smaller than largest blocked gene+clearance \\
               interval ({max_blocked} bp).",
        i = "Worst interval: {worst$chr}:{worst$start}-{worst$end}",
        i = "Increase limit to at least {max_blocked} (gene width + 2*clearance)."
      )
    )
  }

  return(invisible(TRUE))
}

#' Merge overlapping gene ranges
#'
#' Combines overlapping or adjacent gene intervals on a single chromosome
#' into merged ranges. Gaps between genes are preserved.
#'
#' @param chr_genes Data frame with gene coordinates for one chromosome
#'   (columns 'start', 'end'). May contain duplicates and overlaps.
#'
#' @return Data frame with columns 'start', 'end' representing merged gene
#'   ranges with no overlaps. Rows are sorted by start position.
#'   Empty data frame if input has 0 rows.
#'
#' @details
#' Algorithm:
#' 1. Sorts genes by start position
#' 2. Iterates through genes, merging overlapping ranges
#' 3. Returns list of disjoint merged intervals
#'
#' @examples
#' \dontrun{
#'   merged <- .merge_overlapping_gene_ranges(chr_genes)
#' }
#'
#' @dev
.merge_overlapping_gene_ranges <- function(chr_genes) {
  if (nrow(chr_genes) == 0L) {
    return(data.frame(start = integer(0L), end = integer(0L)))
  }

  ranges <- chr_genes[, c("start", "end")]
  ranges <- ranges[order(ranges$start, ranges$end), , drop = FALSE]

  s <- ranges$start
  e <- ranges$end

  # A new block starts whenever the next start is > current end
  new_block <- c(TRUE, s[-1L] > cummax(e)[-length(e)])

  block_id <- cumsum(new_block)

  data.frame(
    start = tapply(s, block_id, min),
    end = tapply(e, block_id, max),
    row.names = NULL
  )
}
