#' Get scShardSplitRef verbosity
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

#' Should scShardSplitRef emit informational chatter?
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
