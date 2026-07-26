#' Map package verbosity to derived options (internal)
#'
#' @param verbosity Character scalar: "quiet", "minimal", or "verbose".
#'
#' @returns Named list of options to be set.
#' @keywords internal
.map_verbosity <- function(verbosity) {
  v <- as.character(verbosity %||% "verbose")
  if (!v %in% c("quiet", "minimal", "verbose")) {
    v <- "verbose"
  }

  list(
    rlib_message_verbosity = switch(
      v,
      quiet = "quiet",
      minimal = "minimal",
      verbose = "verbose"
    ),
    rlib_warning_verbosity = switch(
      v,
      quiet = "quiet",
      minimal = "verbose",
      verbose = "verbose"
    ),
    warn = switch(
      v,
      quiet = -1L,
      minimal = 0L,
      verbose = 0L
    ),
    datatable.showProgress = switch(
      v,
      quiet = FALSE,
      minimal = FALSE,
      verbose = TRUE
    )
  )
}

#' Initialize package options (internal)
#'
#' Extracted from `.onLoad()` to allow direct testing.
#'
#' @returns Invisibly NULL.
#' @keywords internal
.init_package_options <- function() {
  op <- options()

  op.defaults <- list(
    scShardSplitRef.verbosity = "verbose"
  )

  # withr::local_options(.local_envir = topenv()) registers automatic
  # restoration when the package namespace is unloaded; no explicit
  # .onUnload()/withr::deferred_run() is needed (or wanted -- calling
  # deferred_run() manually on top of this double-registers/races with
  # withr's own automatic restoration and produces a spurious "No deferred
  # expressions to run." message, notably during repeated
  # devtools::document()/pkgload::load_all() cycles).
  toset <- !(names(op.defaults) %in% names(op))
  if (any(toset)) {
    withr::local_options(op.defaults[toset], .local_envir = topenv())
  }

  verbosity <- getOption("scShardSplitRef.verbosity")
  mapped <- .map_verbosity(verbosity)

  withr::local_options(mapped, .local_envir = topenv())

  invisible(NULL)
}

# nocov start
.onLoad <- function(libname, pkgname) {
  .init_package_options()
}
# nocov end
