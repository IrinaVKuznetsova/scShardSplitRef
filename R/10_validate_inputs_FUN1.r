#' @import dplyr
#' @import stringr
#' @importFrom utils read.table write.table













# 16 Feb 2026
# IK
# License:  GPL-3.0 



# FUNCTION-1 
# Sets formatting checks of GTF and BED files and raises error


# mamba activate 00_R_env
#library(tidyverse)
validate_inputs <- function(gtf_path, centromere_path) {

    # --- ERROR 0: Check file paths ----------------------------------------------------
    if (!file.exists(gtf_path)) {
        stop("ERROR 0: GTF file does not exist. Check the provided directory:", gtf_path)
    }
    if (!file.exists(centromere_path)) {
        stop("ERROR 1: Centromere BED file does not exist. Check the provided directory:", centromere_path)
    }

    # --- Helper: check TAB-delimited structure ---------------------------------------
    check_tab_file <- function(path, min_fields, max_fields = NULL, label, error_num) {

        lines <- suppressWarnings(readLines(path))
        lines <- lines[!grepl("^#", lines)]   # remove comment lines

        # Split strictly on TAB
        split_lines <- strsplit(lines, "\t")

        # Count fields per line
        field_counts <- sapply(split_lines, length)

        # Check minimum fields
        bad_min <- which(field_counts < min_fields)

        # Check maximum fields (optional)
        bad_max <- if (!is.null(max_fields)) which(field_counts > max_fields) else integer(0)

        # Combine bad rows
        bad <- sort(unique(c(bad_min, bad_max)))

        if (length(bad) > 0) {
            stop(
                paste0(
                    "ERROR ", error_num, ": Invalid ", label, " formatting.\n",
                    "Problematic line numbers: ", paste(bad, collapse = ", "), "\n\n",
                    "Example problematic line:\n",
                    lines[bad[1]], "\n\n",
                    "Expected: ",
                    if (!is.null(max_fields))
                        paste0(min_fields, "-", max_fields, " TAB-separated fields")
                    else
                        paste0("=", min_fields, " TAB-separated fields"),
                    "\nThis file must be TAB-delimited.\n"
                )
            )
        }

        #return(split_lines)
    }

    # --- ERROR 2: Validate GTF --------------------------------------------------------
    # GTF must have at least 9 TAB-separated fields
    # (attributes column may contain spaces but must be in field 9)
    check_tab_file(
        path = gtf_path,
        min_fields = 9,
        label = "GTF",
        error_num = 2
    )

    # --- ERROR 3: Validate BED --------------------------------------------------------
    # BED must have exactly 3 TAB-separated fields (Chr, start, end)
    check_tab_file(
        path = centromere_path,
        min_fields = 3,
        max_fields = 3,
        label = "BED",
        error_num = 3
    )

    # --- If all checks passed, read files safely --------------------------------------
    gtf <- read.table(gtf_path, sep = "\t", header = FALSE, quote = "", comment.char = "#")
    colnames(gtf) <- c("Chr", "source", "feature", "start", "end",
                       "score", "strand", "frame", "attribute")

    bed <- read.table(centromere_path, sep = "\t", header = FALSE)
    colnames(bed) <- c("Chr", "CentrStart", "ChrEND")

    return(list(gtf = gtf, bed = bed))
}
