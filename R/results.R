#' @export
print.et_fit <- function(x, ...) {
  cat("<et_fit> ", x$n_sweeps, " sweeps (burn ", x$n_burn, ", seed ", x$seed,
      ") in ", sprintf("%.1f", x$elapsed), " s\n", sep = "")
  cat("  ", sprintf("%.3f", x$elapsed / max(1, x$n_sweeps + x$n_burn)),
      " s/sweep\n", sep = "")
  cat("  blocks:\n")
  for (b in x$blocks) {
    cat("    ", b$kind, ": ", paste(b$vars, collapse = ", "), "\n", sep = "")
  }
  for (term in names(x$depends)) {
    d <- x$depends[[term]]
    cat("  depends[", term, "]: ",
        if (is.null(d)) "(not derived - an et_julia() body is opaque)"
        else paste(d, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' @export
summary.et_fit <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
  df <- as.data.frame(object)
  stats <- t(vapply(df, function(v) {
    c(mean = mean(v), sd = stats::sd(v), stats::quantile(v, probs, names = FALSE))
  }, numeric(2 + length(probs))))
  colnames(stats) <- c("mean", "sd", paste0(probs * 100, "%"))
  structure(list(stats = as.data.frame(stats), fit = object),
            class = "summary.et_fit")
}

#' @export
print.summary.et_fit <- function(x, digits = 4, ...) {
  cat("Posterior summary (", nrow(as.data.frame(x$fit)), " draws)\n\n", sep = "")
  print(round(x$stats, digits))
  invisible(x)
}

#' Flatten a fit's draws to a data frame, one column per scalar quantity.
#'
#' A vector parameter `alpha` of length 3 becomes `alpha[1]`, `alpha[2]`,
#' `alpha[3]`; a matrix parameter is flattened column-major with `[i,j]` names.
#'
#' @param x An `et_fit`.
#' @param ... Ignored.
#' @return A data frame with one row per draw.
#' @export
as.data.frame.et_fit <- function(x, ...) {
  cols <- list()
  for (nm in names(x$draws)) {
    v <- x$draws[[nm]]
    if (is.null(dim(v))) {
      cols[[nm]] <- as.numeric(v)
      next
    }
    if (length(dim(v)) == 2L) {
      for (k in seq_len(ncol(v))) cols[[sprintf("%s[%d]", nm, k)]] <- v[, k]
      next
    }
    d <- dim(v)[-1]
    idx <- expand.grid(lapply(d, seq_len))
    flat <- matrix(v, nrow = dim(v)[1])
    for (k in seq_len(ncol(flat))) {
      cols[[sprintf("%s[%s]", nm, paste(idx[k, ], collapse = ","))]] <- flat[, k]
    }
  }
  as.data.frame(cols, check.names = FALSE)
}

#' Convert a fit to a `coda::mcmc` object.
#'
#' Opens the fit to R's MCMC diagnostics ecosystem (`coda::effectiveSize`,
#' `coda::gelman.diag` across chains, trace plots).
#'
#' @param x An `et_fit`.
#' @param ... Ignored.
#' @return A `coda::mcmc` object.
#' @export
as_mcmc <- function(x, ...) UseMethod("as_mcmc")

#' @rdname as_mcmc
#' @export
as_mcmc.et_fit <- function(x, ...) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("as_mcmc() needs the coda package.", call. = FALSE)
  }
  coda::mcmc(as.matrix(as.data.frame(x)), start = x$n_burn + 1L, thin = 1L)
}

#' The generated Julia source behind a fit.
#' @param x An `et_fit`.
#' @return The Julia source, invisibly; printed by default.
#' @export
et_source <- function(x) {
  stopifnot(inherits(x, "et_fit"))
  cat(x$source$src)
  invisible(x$source$src)
}
