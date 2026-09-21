# Validator for the derived `depends=` annotations.
#
# Under-declaring `depends=` is silent: a real gradient contribution is
# dropped with no error and no warning. The annotation is derived from the
# transpiled bodies rather than hand-written, and this is the check on that:
# move each parameter, confirm no term it was excluded from moves.

#' Check the derived `depends=` annotations against the model's behaviour.
#'
#' For each likelihood term, perturbs every parameter **not** in that term's
#' declared dependency set and asserts the term's value does not change. A term
#' that does change has an under-declared `depends=`.
#'
#' This is a probe, not a proof: a dependency that happens to be inactive at the
#' tested parameter values would pass. Run it at more than one `factor`, and at a
#' realistic `x_init`, if the model has regime-dependent structure.
#'
#' @param model An [et_model()].
#' @param blocks The sampler blocks, as for [et_loglik()].
#' @param factor Multiplicative perturbation. Shrinks toward zero by default so
#'   probability parameters stay inside their support.
#' @param x_init Optional starting trajectory to evaluate at.
#' @param tol Absolute tolerance in log-density units.
#' @return A data frame with one row per (parameter, term): the declared
#'   membership, the observed change, and a verdict. Invisibly, plus a message
#'   summarising. Rows with `verdict == "UNDER-DECLARED"` are bugs.
#' @export
et_check_depends <- function(model, blocks = list(), factor = 0.7, x_init = NULL,
                             tol = 1e-8) {
  et_require_session()
  gen <- et_julia_source(model, blocks)
  mod <- et_load_module(gen)

  arg <- ""
  if (!is.null(x_init)) {
    x_init <- as.matrix(x_init); storage.mode(x_init) <- "integer"
    check_x_init(x_init, model)
    nm <- paste0(mod, "_depx")
    JuliaCall::julia_assign(nm, x_init)
    arg <- sprintf(", x_init=Matrix{Int}(Main.%s)", nm)
  }
  raw <- JuliaCall::julia_eval(sprintf("%s.et_check_depends(; factor=%s%s)",
                                       mod, julia_float(factor), arg))

  terms <- c("epidemic", "observation")
  has_obs <- !is.null(model$data$observation_process) ||
             !is.null(model$data$observation_weight)
  if (!has_obs) terms <- "epidemic"

  rows <- list()
  for (par in names(raw)) {
    deltas <- as.numeric(raw[[par]])
    for (k in seq_along(terms)) {
      term <- terms[k]
      declared_set <- gen$depends[[term]]
      derived <- !is.null(declared_set)
      declared <- derived && (par %in% declared_set)
      moved <- abs(deltas[k]) > tol
      verdict <- if (!derived) "not derived"
                 else if (moved && !declared) "UNDER-DECLARED"
                 else if (!moved && declared) "over-declared (harmless)"
                 else "ok"
      rows[[length(rows) + 1L]] <- data.frame(
        parameter = par, term = term, declared = declared,
        delta = deltas[k], verdict = verdict, stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  out <- out[order(out$term, out$parameter), , drop = FALSE]
  rownames(out) <- NULL

  bad <- sum(out$verdict == "UNDER-DECLARED")
  if (bad) {
    warning(bad, " under-declared dependency(ies) found. The affected ",
            "parameters' gradients are silently WRONG. See the returned table.",
            call. = FALSE)
  } else {
    message("depends= check passed: no term moves under a parameter it excludes.")
  }
  out
}
