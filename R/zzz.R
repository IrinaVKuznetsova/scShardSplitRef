#' Map package verbosity to derived options (internal)
#'
#' @param verbosity Character scalar: "quiet", "minimal", or "verbose".
#'
#' @return Named list of options to be set.
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
#' @return Invisibly NULL.
#' @keywords internal
.init_package_options <- function() {
  penv <- parent.env(environment())

  op <- options()
  saved <- op[
    names(op) %in%
      c(
        "rlib_message_verbosity",
        "rlib_warning_verbosity",
        "warn"
      )
  ]

  pkg_env <- new.env(parent = emptyenv())
  pkg_env$old_options <- saved

  if (!exists(".pkg_env", envir = penv, inherits = FALSE)) {
    assign(".pkg_env", pkg_env, envir = penv)
  }

  op.defaults <- list(
    scShardSplitRef.verbosity = "verbose"
  )

  toset <- !(names(op.defaults) %in% names(op))
  if (any(toset)) {
    withr::local_options(op.defaults[toset], .local_envir = penv)
  }

  verbosity <- getOption("scShardSplitRef.verbosity")
  mapped <- .map_verbosity(verbosity)

  withr::local_options(mapped, .local_envir = penv)

  invisible(NULL)
}

# nocov start
.onLoad <- function(libname, pkgname) {
  .init_package_options()
}

.onUnload <- function(libpath) {
  penv <- parent.env(environment())
  withr::deferred_run(penv)
  invisible()
}
# nocov end
