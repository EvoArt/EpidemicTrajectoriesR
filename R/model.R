# Parameters, priors, deterministics, and the two likelihood terms.

#' Declare a model parameter and its prior.
#'
#' @param dist The prior distribution: [beta_dist()], [gamma_dist()] and
#'   friends, or [custom_dist()] for anything else. Omit only for
#'   `kind = "latent"`. For a vector parameter whose elements have different
#'   priors, a list of distributions, one per element (emitted as
#'   `PracticalBayes.arraydist`).
#' @param init Initial value. Length must match `n` (or `prod(dim)`).
#' @param n Length of a vector parameter; `1` for a scalar. A vector parameter is
#'   emitted as `PracticalBayes.filldist(dist, n)`. With a list `dist`, it
#'   defaults to the list's length.
#' @param dim Dimensions of a matrix-valued parameter (`kind = "latent"` only).
#' @param kind `"sampled"` (the default) for a parameter with a prior, or
#'   `"latent"` for one **owned by a conjugate kernel**: emitted as a
#'   placeholder distribution whose density is constant, such as an initial
#'   infection probability updated by its own kernel.
#' @param ... For the deprecated `et_par()`, passed on to `prior()`.
#' @return An object of class `prior`.
#' @export
#' @examples
#' parameters <- list(
#'   alpha = prior(gamma_dist(2, 0.005), init = 0.005),
#'   p     = prior(beta_dist(2, 2),      init = 0.6))
prior <- function(dist = NULL, init, n = 1L, dim = NULL,
                  kind = c("sampled", "latent")) {
  kind <- match.arg(kind)
  # One distribution per element: `list(beta_dist(..), beta_dist(..))`.
  if (is.list(dist) && !inherits(dist, "custom_dist")) {
    if (!length(dist) || !all(vapply(dist, inherits, logical(1), "custom_dist"))) {
      stop("prior(): a list `dist` must hold one distribution per element.",
           call. = FALSE)
    }
    if (missing(n)) n <- length(dist)
    if (n != length(dist)) {
      stop("prior(): `dist` has ", length(dist), " distributions but n = ", n, ".",
           call. = FALSE)
    }
    if (n < 2L) {
      stop("prior(): a list of distributions is for a vector parameter; give ",
           "one distribution for a scalar.", call. = FALSE)
    }
    if (kind != "sampled") {
      stop("prior(): a list of distributions is for a sampled parameter.", call. = FALSE)
    }
    dist <- structure(dist, class = c("dist_list", "custom_dist"))
  } else if (kind == "sampled" && !inherits(dist, "custom_dist")) {
    stop("prior(): a sampled parameter needs a distribution from ",
         "normal_dist(), gamma_dist(), beta_dist(), ... (or custom_dist() ",
         "for anything else).", call. = FALSE)
  }
  if (kind == "latent" && is.null(dim)) dim <- length(init)
  init <- as.numeric(init)
  expected <- if (!is.null(dim)) prod(dim) else n
  if (length(init) != expected) {
    stop("prior(): `init` has length ", length(init), " but the parameter has ",
         expected, " element(s). A wrong length here becomes a wrong Julia ",
         "indexing convention later, so it is checked now.", call. = FALSE)
  }
  structure(list(prior = dist, init = init, n = as.integer(n),
                 dim = if (is.null(dim)) NULL else as.integer(dim),
                 kind = kind), class = "prior")
}

#' @export
print.prior <- function(x, ...) {
  cat("<prior> ", x$kind,
      if (!is.null(x$prior)) paste0(" ~ ", dist_to_julia(x$prior)) else "",
      "  init=[", paste(signif(x$init, 4), collapse = ", "), "]\n", sep = "")
  invisible(x)
}

#' Assemble the full model: data, parameters, priors and deterministics.
#'
#' @param data An [et_data()] object.
#' @param parameters Named list of [prior()] declarations.
#' @param derived Named list of `quote()`d expressions over other parameter
#'   names, PracticalBayes deterministics (`:=`). E.g.
#'   `list(m = quote(m_tilde + 1))`.
#' @param entry_time Optional per-individual entry time: a numeric vector, or the
#'   name of an entry in `data`'s `extras`. Supplying it switches on ET's entry
#'   gate, which scores disease transitions before entry but not the survival
#'   factor. Requires the transitions to declare an [et_survival()]: the same
#'   function is passed through automatically, which is the only way the
#'   subtraction removes exactly what was multiplied in.
#' @param likelihood How the hidden trajectories are handled.
#'
#'   * `"augmented"` (the default): the trajectory is a latent block, resampled
#'     by iFFBS once per Gibbs sweep, and the parameters see the complete-data
#'     likelihood given it. Works for any model, coupled or not.
#'   * `"marginal"`: the trajectories are summed out exactly by the forward
#'     algorithm, and the parameters are fitted against the marginal likelihood
#'     alone, with no latent block, no `x_init` and no conjugate kernels. Valid
#'     only when individuals are independent given the parameters: give each
#'     individual its own group, or declare `coupled_transitions = list()` in
#'     [et_data()]. The rates must not read the aggregates, and the observation
#'     and starting-state functions must not read `X`; this is checked
#'     numerically when the module is first used.
#'
#'   Cheaper and exactly mixing where it applies, because there is no
#'   trajectory to mix over.
#'
#'   For an [et_latent_data()] model this is ignored: its chains are always
#'   summed out, and its latent variables are the latent block.
#' @param log_density Extra terms of the log density over the parameters alone:
#'   a `function(model)`, or a list of them, each returning a number that is
#'   added to the log density (PracticalBayes' `@addlogprob!`). For what a prior
#'   cannot say: a penalty (towards `Se + Sp > 1`, say), or the log-Jacobian of
#'   a parameter sampled through a transformation, with its carrier's prior
#'   cancelled. The bodies are transpiled like the rest of the model, may call
#'   the data's helpers, and may use [et_logpdf()]; `model$x` reads the
#'   parameters and the `derived` values.
#' @return An object of class `et_model`.
#' @export
et_model <- function(data, parameters, derived = list(), entry_time = NULL,
                     likelihood = c("augmented", "marginal"), log_density = list()) {
  likelihood <- match.arg(likelihood)
  latent <- inherits(data, "et_latent_data")
  if (!inherits(data, "et_data") && !latent) {
    stop("et_model(): `data` must come from et_data() or et_latent_data().",
         call. = FALSE)
  }
  if (latent) {
    likelihood <- "latent"
    if (!is.null(entry_time)) {
      stop("et_model(): `entry_time` does not apply to et_latent_data(); write ",
           "entry into the weights.", call. = FALSE)
    }
  }
  if (!is.list(parameters) || !length(parameters) || is.null(names(parameters))) {
    stop("et_model(): `parameters` must be a non-empty NAMED list of prior().",
         call. = FALSE)
  }
  check_julia_name(names(parameters), "parameter name")
  for (nm in names(parameters)) {
    if (!inherits(parameters[[nm]], "prior")) {
      stop("et_model(): parameter '", nm, "' must be created with prior().",
           call. = FALSE)
    }
  }
  if ("X" %in% names(parameters)) {
    stop("et_model(): 'X' is reserved for the latent trajectory block.",
         call. = FALSE)
  }
  if (latent) {
    clash <- intersect(names(parameters), c("latents", names(data$latents)))
    if (length(clash)) {
      stop("et_model(): parameter name(s) clash with the latent block: ",
           paste(clash, collapse = ", "), ".", call. = FALSE)
    }
    owned <- names(parameters)[vapply(parameters, function(p) p$kind == "latent",
                                      logical(1))]
    if (length(owned)) {
      stop("et_model(): an et_latent_data() model has no conjugate kernels, so ",
           "kind = \"latent\" parameter(s) ", paste(owned, collapse = ", "),
           " would never be informed. Give them a prior.", call. = FALSE)
    }
  }
  if (!is.list(derived)) stop("et_model(): `derived` must be a list.", call. = FALSE)
  if (length(derived)) {
    if (is.null(names(derived)) || any(!nzchar(names(derived)))) {
      stop("et_model(): every entry of `derived` must be named.", call. = FALSE)
    }
    check_julia_name(names(derived), "derived name")
    clash <- intersect(names(derived), names(parameters))
    if (length(clash)) {
      stop("et_model(): derived name(s) clash with parameters: ",
           paste(clash, collapse = ", "), ".", call. = FALSE)
    }
    known <- c(names(parameters), names(derived))
    for (nm in names(derived)) {
      unknown <- missing_from(expr_symbols(derived[[nm]]), known)
      unknown <- setdiff(unknown, c(names(.et_builtin_map), "T", "F", "TRUE", "FALSE"))
      if (length(unknown)) {
        stop("et_model(): derived '", nm, "' refers to unknown name(s): ",
             paste(unknown, collapse = ", "),
             ". Only parameters and other derived values are in scope.",
             call. = FALSE)
      }
    }
  }

  if (!is.null(entry_time)) {
    if (is.null(data$transitions$survival)) {
      stop("et_model(): `entry_time` requires the transitions to declare an ",
           "et_survival(). The entry gate removes the survival factor before ",
           "entry, so it has to know which factor that is.", call. = FALSE)
    }
    if (is.character(entry_time)) {
      if (!(entry_time %in% names(data$extras))) {
        stop("et_model(): entry_time '", entry_time, "' is not an entry of ",
             "`extras`. Available: ", paste(names(data$extras), collapse = ", "),
             ".", call. = FALSE)
      }
    } else {
      entry_time <- as.integer(entry_time)
      if (length(entry_time) != data$n_individuals) {
        stop("et_model(): `entry_time` has length ", length(entry_time),
             " but there are ", data$n_individuals, " individuals.", call. = FALSE)
      }
    }
  }

  if (likelihood == "marginal") check_marginal(data, parameters, entry_time)

  if (is.function(log_density) || inherits(log_density, "et_julia")) {
    log_density <- list(log_density)
  }
  if (!is.list(log_density)) {
    stop("et_model(): `log_density` must be a function(model), or a list of them.",
         call. = FALSE)
  }
  for (k in seq_along(log_density)) {
    f <- log_density[[k]]
    if (inherits(f, "et_julia")) next
    if (!is.function(f) || !identical(names(formals(f)), "model")) {
      stop("et_model(): log_density term ", k, " must be a function(model).",
           call. = FALSE)
    }
  }

  # The reference set every downstream check uses: which names a parameter
  # expression may legally mention.
  structure(list(data = data, parameters = parameters, derived = derived,
                 entry_time = entry_time, likelihood = likelihood,
                 log_density = unname(log_density),
                 par_names = names(parameters),
                 all_names = c(names(parameters), names(derived))),
            class = "et_model")
}

# The log_density terms: their functions, and for each, what it reads. A term
# that reads a name that is neither a parameter nor a derived value would fail
# at run time with a missing field, so it is refused here.
log_density_terms <- function(model) {
  helpers <- model$data$helpers %||% list()
  lapply(seq_along(model$log_density), function(k) {
    nm <- sprintf("et_log_density_%d", k)
    r <- et_transpile_role(model$log_density[[k]], "log_density", nm, helpers)
    unknown <- setdiff(r$reads, model$all_names)
    if (length(unknown)) {
      stop("et_model(): log_density term ", k, " reads model$",
           paste(unknown, collapse = ", model$"), ", which is not a parameter or ",
           "a derived value.", call. = FALSE)
    }
    list(name = nm, src = r$src,
         depends = if (r$opaque) NULL else expand_to_sampled(r$reads, model))
  })
}

log_density_depends <- function(model) {
  terms <- log_density_terms(model)
  stats::setNames(lapply(terms, `[[`, "depends"),
                  vapply(terms, `[[`, character(1), "name"))
}

# The model-body lines adding each term; `depends=` where the term's reads are
# known and the model kind takes it.
log_density_lines <- function(model, with_depends = TRUE) {
  vapply(log_density_terms(model), function(t) {
    dep <- if (with_depends && !is.null(t$depends))
      paste0("  depends=", julia_symbol_tuple(t$depends)) else ""
    sprintf("@addlogprob! %s(pars)%s", t$name, dep)
  }, character(1))
}

log_density_src <- function(model) {
  vapply(log_density_terms(model), `[[`, character(1), "src")
}

is_marginal <- function(model) identical(model$likelihood, "marginal")

# For the entry points that read a sampled trajectory, which a collapsed fit
# does not have.
require_augmented <- function(model, caller) {
  if (is_marginal(model) || is_latent(model)) {
    stop(caller, " needs a sampled trajectory, and ",
         if (is_latent(model)) "an et_latent_data() model"
         else "a model with likelihood = \"marginal\"",
         " has none: its trajectories are summed out.", call. = FALSE)
  }
  invisible(TRUE)
}

# Both classes whose trajectories are summed out: no `X`, no iFFBS block, no
# conjugate kernel.
is_collapsed <- function(model) is_marginal(model) || is_latent(model)

# What a collapsed model cannot carry, refused at declaration rather than
# discovered inside Julia. ET's own `require_independent` makes the coupling
# check again on the built data, and its numeric tripwire covers what no
# declaration can: a rate or callback that reads the latent state.
check_marginal <- function(data, parameters, entry_time) {
  if (!is.null(data$latents)) {
    stop("et_model(): likelihood = \"marginal\" does not take et_data() latents; ",
         "their update conditions on a sampled trajectory. For latents with the ",
         "chains summed out, use et_latent_data().", call. = FALSE)
  }
  if (!is.null(entry_time)) {
    stop("et_model(): likelihood = \"marginal\" does not support `entry_time`. ",
         "The entry gate changes the per-step factor the forward algorithm ",
         "sums over.", call. = FALSE)
  }
  if (!is.null(data$likelihood_weight)) {
    stop("et_model(): likelihood = \"marginal\" uses the FULL observation ",
         "model. `likelihood_weight` splits off a factor for a conjugate ",
         "kernel, and a collapsed model has no trajectory for such a kernel to ",
         "condition on. Drop it and sample those parameters in the marginal ",
         "likelihood with their actual priors.", call. = FALSE)
  }
  latent <- names(parameters)[vapply(parameters, function(p) p$kind == "latent",
                                     logical(1))]
  if (length(latent)) {
    stop("et_model(): likelihood = \"marginal\" has no conjugate kernels, so ",
         "parameter(s) ", paste(latent, collapse = ", "), " declared with ",
         "kind = \"latent\" would never be informed. Give them a prior.",
         call. = FALSE)
  }
  grouped <- !is.null(data$group) && anyDuplicated(data$group) > 0L
  shared <- is.null(data$group) || grouped ||
    !is.null(data$affected_individuals)
  declared_none <- !is.null(data$coupled_transitions) &&
    length(data$coupled_transitions) == 0L
  if (shared && !declared_none) {
    stop("et_model(): likelihood = \"marginal\" needs individuals that are ",
         "independent given the parameters, but they share groups or ",
         "neighbours and nothing declares that none of their transitions ",
         "depend on each other. Give each individual its own group, or pass ",
         "coupled_transitions = list() to et_data() if the model has no ",
         "between-individual coupling.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @export
print.et_model <- function(x, ...) {
  cat("<et_model>", if (is_marginal(x)) " marginal likelihood (trajectories summed out)",
      if (is_latent(x)) " weighted chains summed out, latents sampled", "\n", sep = "")
  print(x$data)
  cat("  parameters:\n")
  for (nm in x$par_names) {
    p <- x$parameters[[nm]]
    sz <- if (!is.null(p$dim)) paste0("[", paste(p$dim, collapse = "x"), "]")
          else if (p$n > 1L) paste0("[", p$n, "]") else ""
    cat("    ", nm, sz, if (p$kind == "latent") "  (conjugate-owned)"
        else paste0(" ~ ", dist_to_julia(p$prior)), "\n", sep = "")
  }
  if (length(x$derived)) {
    for (nm in names(x$derived)) {
      cat("    ", nm, " := ", deparse(x$derived[[nm]]), "\n", sep = "")
    }
  }
  invisible(x)
}

# Every symbol appearing in an R expression.
expr_symbols <- function(e) {
  if (is.symbol(e)) return(as.character(e))
  if (!is.call(e)) return(character())
  unique(unlist(lapply(as.list(e)[-1L], expr_symbols)))
}

# Expand a set of names to the sampled parameters behind them: a derived value
# resolves to the parameters its expression reads, transitively. This is what
# keeps `depends=` honest when a rate function reads `model$m` and `m` is
# `m_tilde + 1`.
expand_to_sampled <- function(names_in, model) {
  out <- character(); seen <- character()
  frontier <- names_in
  while (length(frontier)) {
    nm <- frontier[1]; frontier <- frontier[-1]
    if (nm %in% seen) next
    seen <- c(seen, nm)
    if (nm %in% model$par_names) { out <- c(out, nm); next }
    if (nm %in% names(model$derived)) {
      frontier <- c(frontier, expr_symbols(model$derived[[nm]]))
    }
    # A name that is neither is not a model variable (a local, a builtin); it
    # contributes no dependency.
  }
  unique(out)
}

#' @rdname prior
#' @export
et_par <- function(...) {
  warning("et_par() is deprecated; use prior()", call. = FALSE)
  prior(...)
}
