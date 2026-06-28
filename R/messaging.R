#' Get scShardSplitRef verbosity (internal)
#'
#' @return Character scalar: "quiet", "minimal", or "verbose".
#' @dev
.scs_get_verbosity <- function() {
  v <- as.character(getOption("scShardSplitRef.verbosity") %||% "verbose")
  if (!v %in% c("quiet", "minimal", "verbose")) {
    v <- "verbose"
  }
  v
}

#' Should scShardSplitRef emit informational chatter? (internal)
#'
#' Policy:
#' - quiet:   no
#' - minimal: no
#' - verbose: yes
#'
#' @return Logical scalar.
#' @dev
.scs_emit_info <- function() {
  .scs_get_verbosity() == "verbose"
}

#' Should scShardSplitRef emit success messages? (internal)
#'
#' Policy matches info by default (only in verbose).
#' If you prefer success in minimal too, change to:
#' `.scs_get_verbosity() %in% c("minimal", "verbose")`.
#'
#' @return Logical scalar.
#' @dev
.scs_emit_success <- function() {
  .scs_get_verbosity() == "verbose"
}

#' scShardSplitRef cli info alert (internal)
#'
#' Wrapper for [cli::cli_alert_info()] that respects
#' `options(scShardSplitRef.verbosity = "quiet"|"minimal"|"verbose")`.
#'
#' @param text Character scalar. A cli glue string (supports `{}` and `{?s}`).
#' @param ... Named values made available for `{}` interpolation.
#' @param .envir Environment used for interpolation. Defaults to caller env.
#'
#' @return Invisibly `TRUE`.
#' @dev
scs_cli_alert_info <- function(text, ..., .envir = parent.frame()) {
  if (!.scs_emit_info()) {
    return(invisible(TRUE))
  }

  dots <- list(...)
  if (length(dots) > 0L) {
    e <- rlang::env(.parent = .envir, !!!dots)
    cli::cli_alert_info(text, .envir = e)
  } else {
    cli::cli_alert_info(text, .envir = .envir)
  }

  invisible(TRUE)
}
