#' Determine split regions for reference sharding
#'
#' Builds BED-like split intervals from genomic regions and GTF gene
#' annotations. If the start position in a genomic region is not 0, it is
#' assumed to mark the start of a centromere and the region is split around
#' that point. Proposed boundaries are shifted rightward when they fall inside
#' gene intervals, and oversized chunks are repeatedly re-split until all
#' segments satisfy the requested size limit.
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
#' @param shift_by Number of base pairs to shift a boundary right when it falls
#'   inside a gene. Defaults to 1 bp.
#'
#' @return Invisible character: path to the written BED file.
#'
#' @details
#' The function performs the following steps:
#' 1. Reads genomic regions and gene annotations from input files.
#' 2. Checks that the provided `limit` is feasible given gene positions.
#' 3. For each region, splits it at the centromere (if start > 0) and applies
#'    the size limit.
#' 4. Adjusts boundaries that fall within genes by shifting them right.
#' 5. Iteratively subdivides intervals exceeding the limit until all intervals
#'    satisfy the constraint.
#' 6. Writes final split regions to BED file in the active \R session's
#'  `tempdir()`.
#'
#' @examples
#'   # load example files from this package and write the output bed file to
#'   # R's tempdir()
#'
#'   chromosomes <- system.file("extdata/A3_toy_centromeres.bed",
#'                    package = "scShardSplitRef",
#'                    mustWork = TRUE)
#'
#'   genes <- system.file("extdata/A3_toy_all_scenarios_2chr.gtf",
#'              package = "scShardSplitRef",
#'              mustWork = TRUE)
#'
#'   determine_split_regions(
#'     bed = chromosomes,
#'     gtf = genes,
#'     output_bed = file.path(tempdir(), "A3_toy_split_regions.bed"),
#'     limit = 2L^29L,
#'     shift_by = 1L
#'   )
#' @autoglobal
#' @export
determine_split_regions <- function(
  bed,
  gtf,
  output_bed,
  limit = 2L^29L,
  shift_by = 1L
) {
  cli::cli_inform(
    "Preparing split regions from {.file bed} and \\
     {.file gtf} into {.file output_bed} \\
     ({.var limit} = {limit}, {.var shift_by} = {shift_by}). 
    \n
    \n"
  )

  regions <- .read_bed_file(bed, 3L)
  genes <- .read_gtf_file(gtf)
  genes <- genes[genes$feature == "gene", , drop = FALSE]

  .assert_feasible_limit(genes, limit, shift_by)

  split_regions <- lapply(
    seq_len(nrow(regions)),
    function(i) {
      .process_single_region(regions[i, ], genes, limit, shift_by)
    }
  )

  out <- do.call(rbind, split_regions)
  .write_bed_file(out, output_bed)

  cli::cli_inform(
    "{.strong Built {nrow(out)} split intervals across \\
     {length(unique(out$chr))} chromosome{?s} and wrote them to \\
     {.file output_bed}}."
  )

  return(invisible(output_bed))
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
      "BED file not found: {.file {path}}",
    )
  }

  d <- utils::read.table(
    path,
    sep = "\t",
    header = FALSE,
    quote = "",
    stringsAsFactors = FALSE,
    fill = TRUE
  )

  ncol_d <- ncol(d)
  if (ncol_d < min_cols) {
    cli::cli_abort(
      call = rlang::caller_env(),
      "BED file must have at least {min_cols} columns, found {ncol_d}."
    )
  }

  colnames(d)[seq_len(3L)] <- c("chr", "start", "end")
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

  d <- utils::read.table(
    path,
    sep = "\t",
    header = FALSE,
    quote = "",
    comment.char = "#",
    stringsAsFactors = FALSE,
    fill = TRUE
  )

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
#' @return Invisible TRUE.
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
  invisible(TRUE)
}

#' Process a single genomic region
#'
#' Splits a single genomic region according to the size limit and adjusts
#' boundaries to avoid genes. If the region has a non-zero start (centromere),
#' it is split into two parts: before and after the centromere.
#'
#' @param region Data frame row with columns 'chr', 'start', 'end'.
#' @param genes Data frame with gene annotations including 'chr', 'start',
#'  'end'.
#' @param limit Integer: maximum allowed interval width in base pairs.
#' @param shift_by Integer: base pairs to shift a boundary right when inside a
#'  gene.
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
#'                                    shift_by = 1)
#' }
#'
#' @dev
.process_single_region <- function(region, genes, limit, shift_by) {
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

  # Fix boundaries that fall inside genes
  chr_genes <- genes[genes$chr == chr, , drop = FALSE]
  .fix_split_boundaries(split_df, chr_genes, limit, shift_by)
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
#' @param chr_genes Data frame with gene annotations for the same chromosome,
#'   including columns 'start', 'end'.
#' @param limit Integer: maximum allowed interval width in base pairs.
#' @param shift_by Integer: base pairs to shift a boundary right when inside
#'   a gene. Defaults to 1L.
#'
#' @return Data frame with same structure as `split_df` but with adjusted
#'   boundaries that satisfy the size and gene-avoidance constraints.
#'
#' @details
#' The algorithm iteratively:
#' 1. Shifts boundaries that fall inside genes rightward by shift_by bp
#' 2. Subdivides intervals exceeding the limit
#' 3. Shifts again after subdivision
#' 4. Repeats until convergence
#'
#' Raises an error if constraints cannot be satisfied.
#'
#' @examples
#' \dontrun{
#'   fixed_df <- fix_split_boundaries(
#'     split_df, chr_genes, limit = 50000, shift_by = 1
#'   )
#' }
#'
#' @dev
.fix_split_boundaries <- function(split_df, chr_genes, limit, shift_by = 1L) {
  if (nrow(split_df) <= 1L || nrow(chr_genes) == 0L) {
    return(split_df)
  }

  chr <- split_df$chr[1L]
  chr_end <- as.integer(max(split_df$end))
  gene_ranges <- .merge_overlapping_gene_ranges(chr_genes)
  boundaries <- as.integer(split_df$end[-nrow(split_df)])

  boundaries <- .refine_boundaries(
    boundaries,
    gene_ranges,
    chr_end,
    limit,
    shift_by
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
#' @param shift_by Integer: base pairs to shift right when inside a gene.
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
#'     c(100L, 200L), gene_ranges, 300L, 50L, 1L
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
  max_iter = 100L
) {
  for (iter in seq_len(max_iter)) {
    old_boundaries <- boundaries

    # Shift past genes
    boundaries <- .shift_boundaries_past_genes(
      boundaries,
      gene_ranges,
      chr_end,
      shift_by
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
      shift_by
    )

    # Check convergence
    if (identical(boundaries, old_boundaries)) {
      return(boundaries)
    }
  }

  # Max iterations exceeded
  return(boundaries)
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
#' Adjusts a vector of boundaries that fall inside gene intervals by shifting
#' them rightward. Each boundary is independently shifted until it no longer
#' overlaps any gene. Results are deduplicated and filtered.
#'
#' @param boundaries Integer vector: boundary positions to adjust.
#' @param gene_ranges Data frame with gene ranges (columns 'start', 'end').
#' @param chr_end Integer: chromosome end coordinate (upper bound).
#' @param shift_by Integer: base pairs to shift right per iteration.
#'   Defaults to 1L.
#'
#' @return Integer vector: adjusted, deduplicated, and filtered boundaries.
#'   Values outside (0, chr_end) are removed.
#'
#' @examples
#' \dontrun{
#'   shifted <- .shift_boundaries_past_genes(
#'     c(60L, 150L), gene_ranges, 300L, shift_by = 1L
#'   )
#' }
#'
#' @dev
.shift_boundaries_past_genes <- function(
  boundaries,
  gene_ranges,
  chr_end,
  shift_by = 1L
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
    shift_by = shift_by
  )

  shifted <- sort(unique(shifted))
  shifted[shifted > 0L & shifted < chr_end]
}

#' Shift a single boundary past genes
#'
#' Adjusts one boundary position that falls inside gene intervals by
#' iteratively shifting it rightward until it no longer overlaps any gene.
#'
#' @param boundary Integer: initial boundary position.
#' @param gene_ranges Data frame with gene ranges (columns 'start', 'end').
#' @param chr_end Integer: chromosome end coordinate (upper bound).
#' @param shift_by Integer: base pairs to shift right per iteration.
#'
#' @return Integer: adjusted boundary position. May be the original position
#'   if it did not overlap any gene, or identical to its start position if
#'   the shift would exceed chr_end.
#'
#' @examples
#' \dontrun{
#'   adjusted <- .shift_single_boundary(60L, gene_ranges, 300L, 1L)
#' }
#'
#' @dev
.shift_single_boundary <- function(boundary, gene_ranges, chr_end, shift_by) {
  repeat {
    inside <- gene_ranges$start < boundary & gene_ranges$end >= boundary

    if (!any(inside)) {
      break
    }

    updated <- as.integer(max(gene_ranges$end[inside])) + as.integer(shift_by)

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

#' Assert that the size limit is feasible
#'
#' Checks that the requested size `limit` is larger than any single gene
#' (plus the shift distance). Raises an informative error if not feasible.
#'
#' @param genes Data frame with gene annotations (columns 'chr', 'start',
#'  'end').
#' @param limit Integer: requested maximum interval width.
#' @param shift_by Integer: base pairs to shift right when inside a gene.
#'
#' @return Invisible TRUE if feasible. Raises error with remediation advice
#'   if not.
#'
#' @details
#' First merges overlapping genes per chromosome, then checks if any merged
#' range (plus shift_by) exceeds the limit. Reports the worst (largest) blocking
#' interval and suggests an appropriate minimum limit value.
#'
#' @examples
#' \dontrun{
#'   .assert_feasible_limit(genes, limit = 50000, shift_by = 1)
#' }
#'
#' @dev
.assert_feasible_limit <- function(genes, limit, shift_by) {
  if (nrow(genes) == 0L) {
    return(invisible(TRUE))
  }

  merged <- do.call(
    rbind,
    lapply(split(genes, genes$chr), function(chr_genes) {
      m <- .merge_overlapping_gene_ranges(chr_genes)
      m$chr <- unique(chr_genes$chr)[1L]
      m
    })
  )

  blocked_width <- merged$end - merged$start + as.integer(shift_by)
  max_blocked <- max(blocked_width, na.rm = TRUE)

  if (limit < max_blocked) {
    worst <- merged[which.max(blocked_width), , drop = FALSE]
    cli::cli_abort(
      call = rlang::caller_env(),
      c(
        x = "Limit {limit} is smaller than largest blocked gene \\
               interval ({max_blocked} bp).",
        i = "Worst interval: {worst$chr}:{worst$start}-{worst$end}",
        i = "Increase limit to at least {max_blocked}."
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
