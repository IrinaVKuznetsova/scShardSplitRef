#' Determine split regions for reference sharding
#'
#' Builds BED-like split intervals from genomic regions and GTF gene
#' annotations. If the start position in the genomic region is not 0, it is
#' assumed to be the start of a centromere, and splitting is initiated from
#' there. Proposed boundaries are shifted rightward when they fall inside
#' a gene, and oversized chunks are re-split until boundaries are valid and
#' segment widths are within the requested limit.
#'
#' @param regions_bed Path to a 3-column BED-like file with chromosome,
#'   start, and end coordinates. If the start position is not 0, it is assumed
#'   to be the start of a centromere.
#' @param genes_gft Path to a GTF file containing gene annotations.
#' @param output_bed Path to the BED file that will be written with the final
#'   split intervals.
#' @param limit Maximum allowed width (in bp) for each output interval. Default
#'   is 2^29 (512 million bp) as per Cell Ranger docs.
#' @param shift_by Number of base pairs to shift a boundary to the right after
#'   the end of an overlapping gene.
#'
#' @return The path to the written BED file.
#'
#' @examples
#' regions_bed <- "inst/extdata/IN0_toy_centromeres_for_gtf.bed"
#' genes_gft <- "inst/extdata/A3_toy_all_scenarios_2chr.gtf"
#' output_bed <- tempfile(fileext = ".bed")
#' determine_split_regions(
#'   regions_bed,
#'   genes_gft,
#'   output_bed,
#'   limit = 92,
#'   shift_by = 1
#' )
#' @autoglobal
#' @export
determine_split_regions <- function(
  regions_bed,
  genes_gft,
  output_bed,
  limit = 2^29,
  shift_by = 1L
) {
  cli::cli_alert_info(
    "Preparing split regions from {.file {regions_bed}} and {.file {genes_gft}} into {.file {output_bed}} (limit = {limit}, shift_by = {shift_by})."
  )

  split_regions <- list()

  # TODO: split this out as seperate function that includes checks for file
  # existence and format, and handle errors gracefully.
  regions <- utils::read.table(
    file = regions_bed,
    sep = "\t",
    header = FALSE,
    quote = "",
    stringsAsFactors = FALSE,
    fill = TRUE,
    col.names = c("chr", "start", "end")
  )

  # TODO: split this out as seperate function that includes checks for file
  # existence and format, and handle errors gracefully.
  genes <- utils::read.table(
    file = genes_gft,
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

  genes <- genes[genes$feature == "gene", ]

  # A valid split is impossible if limit is smaller than any merged blocked
  # gene interval once the boundary shift is applied.
  if (nrow(genes) > 0L) {
    # The later boundary-fixing logic works on merged blocked intervals rather
    # than individual genes, so the up-front feasibility check must use the
    # same representation or it can miss impossible cases.
    merged_gene_ranges <- do.call(
      rbind,
      lapply(split(genes, genes$chr), function(chr_genes) {
        merged_ranges <- merge_overlapping_gene_ranges(chr_genes)
        merged_ranges$chr <- unique(chr_genes$chr)[1L]
        merged_ranges
      })
    )
    # A boundary shifted off the right edge of a blocked interval needs
    # `shift_by` extra bases of clearance, so that extra width is part of the
    # minimum feasible chunk size.
    blocked_widths <-
      merged_gene_ranges$end - merged_gene_ranges$start + as.integer(shift_by)
    max_blocked_width <- max(blocked_widths, na.rm = TRUE)

    if (limit < max_blocked_width) {
      largest_block <- merged_gene_ranges[
        which.max(blocked_widths),
        , drop = FALSE
      ]
      cli::cli_abort(c(
        "x" = paste(
          "`limit` (", limit, ") is smaller than the largest merged",
          "blocked gene interval after shifting (", max_blocked_width,
          " bp)."
        ),
        "i" = paste(
          "Largest blocking interval:",
          sprintf(
            "%s:%s-%s",
            largest_block$chr,
            largest_block$start,
            largest_block$end
          )
        ),
        "i" = paste(
          "Increase `limit` to at least", max_blocked_width,
          "to produce valid splits."
        )
      ))
    }
  }

  for (i in seq_len(nrow(regions))) {
    region <- regions[i, ]
    chr_genes <- genes[genes$chr == region$chr, ]

    # If start is not 0, assume it is a centromere start and split around it.
    if (region$start != 0L) {
      cli::cli_alert_warning(
        "Region {region$chr}:{region$start}-{region$end} does not start at 0; assuming {region$start} is centromere start."
      )

      # Build boundaries for the pre-centromere region (0 to region$start).
      if (region$start <= limit) {
        # No need to split the pre-centromere region.
        pre_boundaries <- integer(0)
      } else {
        n_chunks <- ceiling(region$start / limit)
        chunk_size <- region$start / n_chunks
        pre_boundaries <- as.integer(
          round(seq(chunk_size, region$start - chunk_size, by = chunk_size))
        )
      }

      # Split coordinates for the current pre-centromere region.
      pre_split_df <- data.frame(
        chr = region$chr,
        start = c(0L, pre_boundaries),
        end = c(pre_boundaries, as.integer(region$start))
      )

      # Build boundaries for the post-centromere region.
      if (region$end - region$start <= limit) {
        # No need to split post-centromere; use region end as boundary.
        post_boundaries <- as.integer(region$end)
      } else {
        n_chunks <- ceiling((region$end - region$start) / limit)
        chunk_size <- (region$end - region$start) / n_chunks
        post_boundaries <- as.integer(
          round(seq(region$start + chunk_size, region$end, by = chunk_size))
        )
      }

      post_split_df <- data.frame(
        chr = region$chr,
        start = c(as.integer(region$start), post_boundaries),
        end = c(post_boundaries, as.integer(region$end))
      )

      split_df <- rbind(pre_split_df, post_split_df)
    } else {
      # Region starts at 0, so split only by size limit.
      if (region$end <= limit) {
        # No split needed; use whole region.
        split_df <- data.frame(
          chr = region$chr,
          start = 0L,
          end = as.integer(region$end)
        )
      } else {
        n_chunks <- ceiling(region$end / limit)
        chunk_size <- region$end / n_chunks
        boundaries <- as.integer(
          round(seq(chunk_size, region$end - chunk_size, by = chunk_size))
        )
        split_df <- data.frame(
          chr = region$chr,
          start = c(0L, boundaries),
          end = c(boundaries, as.integer(region$end))
        )
      }
    }

    # Shift boundaries away from genes, then re-split oversize chunks.
    split_df <- fix_split_boundaries(split_df, chr_genes, limit, shift_by)

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
    "Built {nrow(out)} split intervals across {length(unique(out$chr))} chromosome(s) and wrote them to {.file {output_bed}}."
  )
  output_bed
}

# Check every internal boundary in split_df against chr_genes. Any boundary
# that falls inside a blocked gene interval is shifted rightward to just past
# that interval. If the shift causes an oversized chunk, the outer loop tries
# to add more boundaries and re-check them.
fix_split_boundaries <- function(split_df, chr_genes, limit, shift_by = 1L) {
  if (nrow(split_df) <= 1L || nrow(chr_genes) == 0L) {
    return(split_df)
  }

  chr <- split_df$chr[1L]
  chr_end <- as.integer(max(split_df$end))
  # Treat overlapping genes as one blocked interval so a boundary cannot be
  # shifted out of one gene and directly into another overlapping gene.
  gene_ranges <- merge_overlapping_gene_ranges(chr_genes)

  # Internal boundaries are the shared points between consecutive rows.
  boundaries <- as.integer(split_df$end[-nrow(split_df)])

  repeat {
    # Keep looping until shifting/re-splitting stops changing the boundaries.
    previous_boundaries <- boundaries

    # First push any existing boundary that still lands inside a blocked gene
    # interval onto a safe edge of that interval.
    boundaries <- shift_boundaries_past_genes(
      boundaries = boundaries,
      gene_ranges = gene_ranges,
      chr_end = chr_end,
      shift_by = shift_by
    )

    # Re-subdivide any chunks that became oversized after shifting. These new
    # boundaries are checked again on the next pass through the outer loop.
    all_points <- c(0L, boundaries, chr_end)
    extra_boundaries <- integer(0)

    for (i in seq_len(length(all_points) - 1L)) {
      seg_start <- all_points[i]
      seg_end <- all_points[i + 1L]
      if ((seg_end - seg_start) > limit) {
        # Even subdivision gives a neutral first guess; the new inner
        # boundaries are validated against blocked genes immediately below.
        n_sub_chunks <- ceiling((seg_end - seg_start) / limit)
        sub_chunk_size <- (seg_end - seg_start) / n_sub_chunks
        inner <- as.integer(round(seq(
          seg_start + sub_chunk_size,
          seg_end - sub_chunk_size,
          by = sub_chunk_size
        )))
        extra_boundaries <- c(extra_boundaries, inner)
      }
    }

    boundaries <- sort(unique(c(boundaries, extra_boundaries)))
    boundaries <- boundaries[boundaries > 0L & boundaries < chr_end]
    # Re-check the newly added boundaries because subdivision can place them
    # back inside blocked gene intervals.
    boundaries <- shift_boundaries_past_genes(
      boundaries = boundaries,
      gene_ranges = gene_ranges,
      chr_end = chr_end,
      shift_by = shift_by
    )

    final_points <- c(0L, boundaries, chr_end)
    oversized_segments <- (final_points[-1L] - final_points[-length(final_points)]) > limit
    if (any(oversized_segments)) {
      # If an oversized segment survives after both passes, the current limit
      # and shift rules cannot produce a valid set of split points.
      bad_idx <- which(oversized_segments)[1L]
      cli::cli_abort(c(
        "x" = paste(
          "Unable to build valid split intervals for chromosome", chr,
          "with `limit` =", limit, "."
        ),
        "i" = paste(
          "Blocked segment remains larger than `limit` after shifting:",
          sprintf(
            "%s:%s-%s",
            chr,
            final_points[bad_idx],
            final_points[bad_idx + 1L]
          )
        ),
        "i" = "Increase `limit` or reduce `shift_by` to allow a valid split."
      ))
    }

    if (identical(boundaries, previous_boundaries)) break
  }

  split_points <- c(0L, boundaries, chr_end)
  data.frame(
    chr = chr,
    start = split_points[-length(split_points)],
    end = split_points[-1L],
    stringsAsFactors = FALSE
  )
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
      inside <- gene_ranges$start < boundary & gene_ranges$end >= boundary
      if (!any(inside)) {
        break
      }

      # When a boundary lands inside one or more overlapping blocked gene
      # intervals, move it to just after the furthest right end so it is no
      # longer inside any of them.
      updated_boundary <-
        as.integer(max(gene_ranges$end[inside])) + as.integer(shift_by)
      if (updated_boundary >= chr_end || updated_boundary <= boundary) {
        break
      }
      boundary <- updated_boundary
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

  # Sort genes by position so overlapping intervals can be collapsed with one
  # pass.
  ranges <- data.frame(
    start = as.integer(chr_genes$start),
    end = as.integer(chr_genes$end),
    stringsAsFactors = FALSE
  )
  ranges <- ranges[order(ranges$start, ranges$end), , drop = FALSE]

  merged_start <- ranges$start[1L]
  merged_end <- ranges$end[1L]
  merged_ranges <- vector("list", nrow(ranges))
  out_idx <- 1L

  if (nrow(ranges) > 1L) {
    for (i in 2:nrow(ranges)) {
      # GTF coordinates are inclusive, so ranges that touch at a boundary are
      # still part of the same blocked interval for splitting.
      if (ranges$start[i] <= merged_end) {
        merged_end <- max(merged_end, ranges$end[i])
      } else {
        # Store the finished interval and start a new one.
        merged_ranges[[out_idx]] <- data.frame(
          start = merged_start,
          end = merged_end
        )
        out_idx <- out_idx + 1L
        merged_start <- ranges$start[i]
        merged_end <- ranges$end[i]
      }
    }
  }

  merged_ranges[[out_idx]] <- data.frame(start = merged_start, end = merged_end)
  do.call(rbind, merged_ranges[seq_len(out_idx)])
}
