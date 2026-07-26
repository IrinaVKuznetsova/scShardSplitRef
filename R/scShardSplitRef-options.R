#' Get or Set scShardSplitRef Options
#'
#' A convenience function to get or set options used by \pkg{scShardSplitRef}.
#'
#' @param ... Named options to set, or no arguments to retrieve current values.
#'
#' @returns A list of current option values.
#'
#' @examples
#' # See currently set options for scShardSplitRef
#' scShardSplitRef_options()
#'
#' # Set verbosity to "quiet" to suppress info/success messages
#' scShardSplitRef_options(scShardSplitRef.verbosity = "quiet")
#' scShardSplitRef_options()
#'
#' @export
#' @family scShardSplitRef-options
#'
scShardSplitRef_options <- function(...) {
  dots <- list(...)
  if (length(dots) == 0L) {
    return(options()[grep("^scShardSplitRef\\.", names(options()))])
  }
  do.call(options, dots)
}
