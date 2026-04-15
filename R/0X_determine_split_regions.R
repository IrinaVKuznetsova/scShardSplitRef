#' Determine split regions for reference sharding
#'
#' Builds BED-like split intervals from genomic regions and GTF gene
#' annotations. If the start position in a genomic region is not 0, it is
#' assumed to mark the start of a centromere and the region is split around
#' that point. Proposed boundaries are shifted rightward when they fall inside
#' gene intervals, and oversized chunks are repeatedly re-split until all
#' segments satisfy the requested size limit.
#'
#' @param regions_bed Path to a 3-column BED-like file (chr, start, end).
#' @param genes_gft Path to a GTF file containing gene annotations.
#' @param output_bed Path to the BED file to be written.
#' @param limit Maximum allowed width (bp) for output intervals.
#' @param shift_by Number of base pairs to shift a boundary to the right when
#'   it falls inside a gene.
#'
#' @return Invisible path to the written BED file.
#' @export
determine_split_regions <- function(
  regions_bed,
  genes_gft,
  output_bed,
  limit = 2L^29L,
  shift_by = 1L
) {
  cli::cli_alert_info(
    "Preparing split regions from {.file {regions_bed}} and {.file {genes_gft}}
     into {.file {output_bed}} (limit = {limit}, shift_by = {shift_by})."
  )

  regions <- utils::read.table(
    regions_bed,
    sep = "\t",
    header = FALSE,
    quote = "",
    stringsAsFactors = FALSE,
    fill = TRUE,
    col.names = c("chr", "start", "end")
  )

  genes <- utils::read.table(
    genes_gft,
    sep = "\t",
    header = FALSE,
    quote = "",
    stringsAsFactors = FALSE,
    fill = TRUE,
    col.names = c(
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
  )
  genes <- genes[genes$feature == "gene", , drop = FALSE]

  assert_feasible_limit(genes, limit, shift_by)

  split_regions <- vector("list", nrow(regions))

  for (i in seq_len(nrow(regions))) {
    region <- regions[i, ]
    chr_genes <- genes[genes$chr == region$chr, , drop = FALSE]

    split_df <- split_region_core(region, limit)

    split_df <- fix_split_boundaries(
      split_df,
      chr_genes,
      limit,
      shift_by
    )

    split_regions[[i]] <- split_df
  }

  out <- do.call(rbind, split_regions)

  utils::write.table(
    out,
    file = output_bed,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    sep = "\t"
  )

  cli::cli_alert_success(
    "Built {nrow(out)} split intervals across
     {length(unique(out$chr))} chromosome(s) and wrote them to
     {.file {output_bed}}."
  )

  invisible(output_bed)
}

# -------------------------------------------------------------------------
# Core splitting helpers
# -------------------------------------------------------------------------

split_region_core <- function(region, limit) {
  chr <- region$chr
  start <- as.integer(region$start)
  end <- as.integer(region$end)

  if (start < 0L) {
    cli::cli_abort("RegionStart must be >= 0 for BED coordinates.")
  }
  if (start >= end) {
    cli::cli_abort("BED regions must satisfy {.field RegionStart < RegionEnd}.")
  }

  if (start == 0L) {
    boundaries <- build_boundaries(0L, end, limit)
    return(boundaries_to_df(chr, 0L, end, boundaries))
  }

  pre_df <- boundaries_to_df(
    chr,
    0L,
    start,
    build_boundaries(0L, start, limit)
  )

  post_df <- boundaries_to_df(
    chr,
    start,
    end,
    build_boundaries(start, end, limit)
  )

  rbind(pre_df, post_df)
}

build_boundaries <- function(start, end, limit) {
  width <- end - start
  if (width <= limit) {
    return(integer(0L))
  }

  n_chunks <- ceiling(width / limit)
  step <- width / n_chunks

  as.integer(round(seq(
    start + step,
    end - step,
    by = step
  )))
}

boundaries_to_df <- function(chr, start, end, boundaries) {
  points <- c(start, boundaries, end)
  data.frame(
    chr = chr,
    start = points[-length(points)],
    end = points[-1L],
    stringsAsFactors = FALSE
  )
}

# -------------------------------------------------------------------------
# Feasibility assertion
# -------------------------------------------------------------------------

assert_feasible_limit <- function(genes, limit, shift_by) {
  if (nrow(genes) == 0L) {
    return(invisible(TRUE))
  }

  merged <- do.call(
    rbind,
    lapply(split(genes, genes$chr), function(chr_genes) {
      m <- merge_overlapping_gene_ranges(chr_genes)
      m$chr <- unique(chr_genes$chr)[1L]
      m
    })
  )

  blocked_width <- merged$end - merged$start + as.integer(shift_by)
  max_blocked <- max(blocked_width, na.rm = TRUE)

  if (limit < max_blocked) {
    worst <- merged[which.max(blocked_width), , drop = FALSE]
    cli::cli_abort(
      c(
        "x" = "{.val limit} {limit} is smaller than the largest merged blocked
               gene interval after shifting ({max_blocked} bp).",
        "i" = "Largest blocking interval: {worst$chr}:
               {worst$start}-{worst$end}",
        "i" = "Increase {.var limit} to at least {max_blocked}."
      )
    )
  }

  invisible(TRUE)
}

# -------------------------------------------------------------------------
# Boundary fixing logic
# -------------------------------------------------------------------------

fix_split_boundaries <- function(split_df, chr_genes, limit, shift_by = 1L) {
  if (nrow(split_df) <= 1L || nrow(chr_genes) == 0L) {
    return(split_df)
  }

  chr <- split_df$chr[1L]
  chr_end <- as.integer(max(split_df$end))
  gene_ranges <- merge_overlapping_gene_ranges(chr_genes)

  boundaries <- as.integer(split_df$end[-nrow(split_df)])

  repeat {
    old <- boundaries

    boundaries <- shift_boundaries_past_genes(
      boundaries,
      gene_ranges,
      chr_end,
      shift_by
    )

    points <- c(0L, boundaries, chr_end)
    boundaries <- sort(unique(c(
      boundaries,
      subdivide_oversized(points, limit)
    )))
    boundaries <- boundaries[boundaries > 0L & boundaries < chr_end]

    boundaries <- shift_boundaries_past_genes(
      boundaries,
      gene_ranges,
      chr_end,
      shift_by
    )

    final <- c(0L, boundaries, chr_end)
    if (any(diff(final) > limit)) {
      cli::cli_abort(
        "Unable to build valid split intervals for chromosome {chr}
         with {.var limit} = {limit}."
      )
    }

    if (identical(boundaries, old)) break
  }

  boundaries_to_df(chr, 0L, chr_end, boundaries)
}

subdivide_oversized <- function(points, limit) {
  extra <- integer(0L)
  for (i in seq_len(length(points) - 1L)) {
    width <- points[i + 1L] - points[i]
    if (width > limit) {
      n <- ceiling(width / limit)
      step <- width / n
      extra <- c(
        extra,
        as.integer(round(seq(
          points[i] + step,
          points[i + 1L] - step,
          by = step
        )))
      )
    }
  }
  extra
}

shift_boundaries_past_genes <- function(
  boundaries,
  gene_ranges,
  chr_end,
  shift_by = 1L
) {
  if (length(boundaries) == 0L || nrow(gene_ranges) == 0L) {
    return(boundaries)
  }

  for (i in seq_along(boundaries)) {
    boundary <- boundaries[i]
    repeat {
      inside <- gene_ranges$start < boundary &
        gene_ranges$end >= boundary
      if (!any(inside)) {
        break
      }

      updated <- as.integer(max(gene_ranges$end[inside])) +
        as.integer(shift_by)
      if (updated >= chr_end || updated <= boundary) {
        break
      }
      boundary <- updated
    }
    boundaries[i] <- boundary
  }

  boundaries <- sort(unique(boundaries))
  boundaries[boundaries > 0L & boundaries < chr_end]
}

merge_overlapping_gene_ranges <- function(chr_genes) {
  if (nrow(chr_genes) == 0L) {
    return(data.frame(start = integer(0), end = integer(0)))
  }

  ranges <- data.frame(
    start = as.integer(chr_genes$start),
    end = as.integer(chr_genes$end),
    stringsAsFactors = FALSE
  )
  ranges <- ranges[order(ranges$start, ranges$end), , drop = FALSE]

  merged_start <- ranges$start[1L]
  merged_end <- ranges$end[1L]
  merged <- list()
  idx <- 1L

  for (i in seq_len(nrow(ranges))[-1L]) {
    if (ranges$start[i] <= merged_end) {
      merged_end <- max(merged_end, ranges$end[i])
    } else {
      merged[[idx]] <- data.frame(
        start = merged_start,
        end = merged_end
      )
      idx <- idx + 1L
      merged_start <- ranges$start[i]
      merged_end <- ranges$end[i]
    }
  }

  merged[[idx]] <- data.frame(start = merged_start, end = merged_end)
  do.call(rbind, merged)
}

# -------------------------------------------------------------------------
# Unit tests (move to tests/testthat/ in package)
# -------------------------------------------------------------------------

if (getRversion() >= "3.5.0") {
  testthat::test_that("build_boundaries respects limit", {
    testthat::expect_equal(
      build_boundaries(0L, 100L, 200L),
      integer(0L)
    )
    testthat::expect_true(length(build_boundaries(0L, 100L, 30L)) > 0L)
  })

  testthat::test_that("split_region_core rejects invalid BED", {
    bad <- data.frame(chr = "chr1", start = 10, end = 10)
    testthat::expect_error(split_region_core(bad, 100))
  })

  testthat::test_that("split_region_core produces valid partitions", {
    r <- data.frame(chr = "chr1", start = 0L, end = 100L)
    out <- split_region_core(r, limit = 30L)
    testthat::expect_true(all(out$end > out$start))
    testthat::expect_true(max(out$end - out$start) <= 30L)
  })

  testthat::test_that("merge_overlapping_gene_ranges collapses overlaps", {
    genes <- data.frame(
      chr = "chr1",
      start = c(10, 15, 30),
      end = c(20, 25, 40)
    )
    merged <- merge_overlapping_gene_ranges(genes)
    testthat::expect_equal(nrow(merged), 2L)
  })
}
