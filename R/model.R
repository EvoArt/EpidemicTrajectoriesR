# Parameters, priors, deterministics, and the two likelihood terms.

#' Declare a model parameter and its prior.
#'
#' @param dist The prior distribution: [beta_dist()], [gamma_dist()] and
#'   friends, or [custom_dist()] for anything else. Omit only for
#'   `kind = "latent"`.
#' @param init Initial value. Length must match `n` (or `prod(dim)`).
#' @param n Length of a vector parameter; `1` for a scalar. A vector parameter is
#'   emitted as `PracticalBayes.filldist(dist, n)`.
#' @param dim Dimensions of a matrix-valued parameter (`kind = "latent"` only).
#' @param kind `"sampled"` (the default) for a parameter with a prior, or
#'   `"latent"` for one **owned by a conjugate kernel**: emitted as a
#'   placeholder distribution whose density is constant, exactly as the badger
#'   model's `nu` is.
#' @return An object of class `prior`.
#' @export
#' @examples
#' parameters <- list(
#'   alpha = prior(gamma_dist(2, 0.005), init = 0.005),
#'   p     = prior(beta_dist(2, 2),      init = 0.6))
prior <- function(dist = NULL, init, n = 1L, dim = NULL,
                  kind = c("sampled", "latent")) {
  kind <- match.arg(kind)
  if (kind == "sampled" && !inherits(dist, "custom_dist")) {
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
#'   trajectory to mix over. It is also what makes [et_lfo_spec()]'s
#'   `method = "exact_hmm"` score available without latent draws.
#' @return An object of class `et_model`.
#' @export
et_model <- function(data, parameters, derived = list(), entry_time = NULL,
                     likelihood = c("augmented", "marginal")) {
  likelihood <- match.arg(likelihood)
  if (!inherits(data, "et_data")) {
    stop("et_model(): `data` must come from et_data().", call. = FALSE)
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

  # The reference set every downstream check uses: which names a parameter
  # expression may legally mention.
  structure(list(data = data, parameters = parameters, derived = derived,
                 entry_time = entry_time, likelihood = likelihood,
                 par_names = names(parameters),
                 all_names = c(names(parameters), names(derived))),
            class = "et_model")
}

is_marginal <- function(model) identical(model$likelihood, "marginal")

# For the entry points that read a sampled trajectory, which a collapsed fit
# does not have.
require_augmented <- function(model, caller) {
  if (is_marginal(model)) {
    stop(caller, " needs a sampled trajectory, and a model with likelihood = ",
         "\"marginal\" has none: its trajectories are summed out.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# What a collapsed model cannot carry, refused at declaration rather than
# discovered inside Julia. ET's own `require_independent` makes the coupling
# check again on the built data, and its numeric tripwire covers what no
# declaration can: a rate or callback that reads the latent state.
check_marginal <- function(data, parameters, entry_time) {
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
      "\n", sep = "")
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
