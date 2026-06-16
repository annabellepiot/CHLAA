# Internal helpers

.check_named_list <- function(x, nm = "x") {
  if (!is.list(x) || is.null(names(x)) || any(names(x) == "")) {
    stop(nm, " must be a named list", call. = FALSE)
  }
  invisible(TRUE)
}

.require_suggested <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for this functionality.", call. = FALSE)
  }
  invisible(TRUE)
}
