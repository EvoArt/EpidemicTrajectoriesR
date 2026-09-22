# Leave-future-out cross-validation: a wrapper over ET's truncation /
# LFOSpec / lfo_cv / compare. Nothing statistical is reimplemented here.
#
# et_run() closes over the module-level `const DATA`, since a normal fit only
# sees the whole series. LFO needs the same compiled model refit against a
# truncated copy at every cutoff, so et_lfo_run_src() emits a second entry
# point taking `data` as an argument. Injected with Base.include_string(), so
# it stays out of et_julia_source() and never perturbs the module-name hash.

#' Declare how each `extras` entry behaves under truncation.
#'
#' R mirror of ET's `truncation()`. Every entry of the model's `extras` (and
#' `group`, if time-varying information is smuggled into it) must be accounted
#' for; an undeclared entry that looks time-indexed is a Julia-side error rather
#' than a silent leak of post-cutoff information into a training fit.
#'
#' @param clamp Character vector of `extras` names to clamp at the cutoff
#'   (`min(v, cutoff)`) -- for "last known present" vectors.
#' @param filter Character vector of `extras` names holding `(id, time)` event
#'   lists, to drop entries after the cutoff.
#' @param copy Character vector of `extras` names to copy unchanged, for
#'   arrays a sampler mutates in place.
#' @param keep Character vector of `extras` names carried through untouched: an
#'   explicit "this is a covariate, not a leak".
#' @param strict Error on an undeclared, time-shaped extra (the default and the
#'   safe choice).
#' @return An object of class `et_lfo_truncation`.
#' @export
et_lfo_truncation <- function(clamp = character(), filter = character(),
                              copy = character(), keep = character(),
                              strict = TRUE) {
  structure(list(clamp = clamp, filter = filter, copy = copy, keep = keep,
                 strict = isTRUE(strict)),
            class = "et_lfo_truncation")
}

#' @export
print.et_lfo_truncation <- function(x, ...) {
  cat("<et_lfo_truncation>\n")
  for (r in c("clamp", "filter", "copy", "keep")) {
    if (length(x[[r]])) cat("  ", r, ": ", paste(x[[r]], collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

# `truncation(...)` call text.
truncation_src <- function(plan) {
  arg <- function(nm) {
    v <- plan[[nm]]
    if (!length(v)) return(NULL)
    sprintf("%s=%s", nm, julia_symbol_tuple(v))
  }
  parts <- Filter(Negate(is.null), lapply(c("clamp", "filter", "copy", "keep"), arg))
  parts <- c(parts, sprintf("strict=%s", if (plan$strict) "true" else "false"))
  sprintf("truncation(%s)", paste(parts, collapse = ", "))
}

#' Score a capture-recapture window.
#'
#' The observation model a mark-recapture study has: a dead individual is never
#' seen, a live one is caught with probability `p`, and if `test` is given,
#' a caught individual is also tested, so the observation carries information
#' about its infection state and not merely that it was alive.
#'
#' **Supply `test` whenever the comparison is about transmission.** Capture
#' alone separates alive from dead and nothing else, so every infection state
#' has the same observation weight; the held-out score then barely depends on
#' the transmission model, and an LFO comparison between transmission models
#' has almost nothing to work with.
#'
#' Supplying this instead of a hand-written [et_julia()] also fixes the two
#' things easiest to get wrong: the density is `-Inf` only where the data
#' actually contradict the trajectory, and the survival constraint is derived
#' from the same capture matrix, so the two can never disagree about who was
#' alive.
#'
#' @param caught Name of an `extras` entry holding an
#'   `n_timepoints x n_individuals` 0/1 capture matrix.
#' @param p Name of the model parameter giving the capture probability of a live
#'   individual.
#' @param unobservable_states States in which an individual cannot be caught, so
#'   that a capture contradicts the trajectory outright. Either integer state
#'   codes, or `"absorbing"` (the default) to use whichever state the model's
#'   declared transitions make absorbing: death, in the usual case. `NULL` for
#'   a model where every state is observable, when the capture term alone is the
#'   whole density.
#'
#'   `"absorbing"` is derived, not declared, so it cannot fall out of step with
#'   the transitions. Note that the *constraint* can only ever forbid the
#'   absorbing state (see `constrain_survival` in [et_lfo_spec()]), while this
#'   argument may name any number of states: a trajectory can be contradicted
#'   by a state the proposal has no way to avoid.
#' @param tests Optional diagnostic tests applied to a caught individual, as a
#'   list of [et_test()]. A bare [et_test()] is accepted for the single-test
#'   case. Results are read only where `caught` is 1.
#'
#'   More than one test is usual rather than exotic: with a single test,
#'   sensitivity, specificity and prevalence all compete to explain the same
#'   positives, and none of them is identified without fixing one.
#' @param infected_states Integer codes of the states a test can detect.
#' @return An object of class `et_capture_recapture`, accepted by
#'   [et_lfo_spec()] as `cell_logdensity`.
#' @export
et_capture_recapture <- function(caught, p, unobservable_states = "absorbing",
                                 tests = NULL, infected_states = NULL) {
  stopifnot(is.character(caught), length(caught) == 1L)
  stopifnot(is.character(p), length(p) == 1L)
  check_julia_name(c(caught, p), "name")
  if (!is.null(unobservable_states) &&
      !identical(unobservable_states, "absorbing")) {
    unobservable_states <- as.integer(unobservable_states)
    if (anyNA(unobservable_states) || !length(unobservable_states)) {
      stop("et_capture_recapture(): `unobservable_states` must be state codes, ",
           "\"absorbing\", or NULL.", call. = FALSE)
    }
  }
  if (inherits(tests, "et_test")) tests <- list(tests)
  if (!is.null(tests)) {
    if (!is.list(tests) || !length(tests) ||
        !all(vapply(tests, inherits, logical(1), "et_test"))) {
      stop("et_capture_recapture(): `tests` must be an et_test(), or a list ",
           "of them.", call. = FALSE)
    }
    if (is.null(infected_states)) {
      stop("et_capture_recapture(): scoring a test needs `infected_states` -- ",
           "which states it is meant to detect.", call. = FALSE)
    }
  }
  structure(list(caught = caught, p = p,
                 unobservable = unobservable_states,
                 tests = tests,
                 infected_states = if (is.null(infected_states)) NULL
                                   else as.integer(infected_states)),
            class = "et_capture_recapture")
}

#' One diagnostic test applied to caught individuals.
#'
#' @param result Name of an `extras` entry holding the result matrix, indexed
#'   `[t, i]`: `1` positive, `0` negative, anything else "not tested".
#' @param sensitivity Name of the model parameter for P(positive | infected).
#' @param specificity Name of the model parameter for P(negative | uninfected).
#'   `NULL` for a test that never gives a false positive, in which case a
#'   positive result rules the uninfected states out entirely.
#' @return An object of class `et_test`.
#' @export
et_test <- function(result, sensitivity, specificity = NULL) {
  stopifnot(is.character(result), length(result) == 1L)
  stopifnot(is.character(sensitivity), length(sensitivity) == 1L)
  check_julia_name(c(result, sensitivity, specificity), "name")
  structure(list(result = result, sensitivity = sensitivity,
                 specificity = specificity), class = "et_test")
}

#' @export
print.et_test <- function(x, ...) {
  cat("<et_test> ", x$result, "  se=", x$sensitivity,
      if (is.null(x$specificity)) "  sp=1 (no false positives)"
      else paste0("  sp=", x$specificity), "\n", sep = "")
  invisible(x)
}

#' @export
print.et_capture_recapture <- function(x, ...) {
  cat("<et_capture_recapture> ", x$caught, " ~ Bernoulli(", x$p,
      ") when alive", sep = "")
  if (is.null(x$tests)) {
    cat("\n  no test: the observation says alive vs dead only, so it carries",
        "\n  no information about infection -- see ?et_capture_recapture\n")
  } else {
    cat("\n  + ", length(x$tests), " test(s) on states {",
        paste(x$infected_states, collapse = ","), "}: ",
        paste(vapply(x$tests, function(tt) tt$result, character(1)),
              collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

# The observation log density. A capture proves the individual was alive, so a
# trajectory with it dead is contradicted (-Inf); otherwise it is the Bernoulli
# capture term, a dead individual explains a non-capture with certainty, and a
# caught individual additionally contributes its test result.
#
# With no `dead_state` the model has no mortality, the contradiction branch
# cannot arise, and the density is the capture term alone. Emitting the branch
# anyway would test the trajectory against a state code that does not exist.
cr_cell_logdensity <- function(cr) {
  body <- paste0(
    sprintf("lp = y == 1 ? log(model.%s) : log1p(-model.%s)\n", cr$p, cr$p),
    paste(vapply(cr$tests %||% list(), cr_test_src, character(1),
                 infected = cr$infected_states),
          collapse = ""),
    "lp")
  head <- sprintf("y = data.%s[t, i]\n", cr$caught)
  if (is.null(cr$unobservable)) return(et_julia(paste0(head, body)))
  # `_et_unobservable` is resolved from the model's own transitions where the
  # user said "absorbing", so it cannot drift out of step with them.
  et_julia(paste0(
    head,
    "if X[t, i] in _et_unobservable\n",
    "    y == 1 ? -Inf : 0.0\n",
    "else\n",
    indent(body), "\n",
    "end"))
}

# The unobservable set, bound in the module. "absorbing" asks ET for the state
# its declared transitions make absorbing: the same function `forward_simulate`
# and `survival_constrained` use, so the score and the proposal agree by
# construction rather than by the user restating a state code.
cr_unobservable_src <- function(cr) {
  if (identical(cr$unobservable, "absorbing")) {
    return(paste0(
      "const _et_unobservable = let a = EpidemicTrajectories._absorbing_state(DATA)\n",
      "    a === nothing ? Int[] : Int[a]\n",
      "end"))
  }
  sprintf("const _et_unobservable = %s",
          julia_vector(cr$unobservable, "Int"))
}

# The test factor, added only where the individual was actually caught. A state
# the test cannot detect contributes the false-positive term, so a positive on
# an uninfected individual is improbable rather than impossible, returning
# -Inf there would discard the draw over a test error.
cr_test_src <- function(test, infected) {
  inf_set <- sprintf("(%s)", paste(c(julia_int(infected), ""), collapse = ", "))
  spec_pos <- if (is.null(test$specificity)) "-Inf"
              else sprintf("log1p(-model.%s)", test$specificity)
  spec_neg <- if (is.null(test$specificity)) "0.0"
              else sprintf("log(model.%s)", test$specificity)
  sprintf(paste0(
    "if y == 1\n",
    "    r = data.%s[t, i]\n",
    "    if r == 1 || r == 0\n",
    "        if X[t, i] in %s\n",
    "            lp += r == 1 ? log(model.%s) : log1p(-model.%s)\n",
    "        else\n",
    "            lp += r == 1 ? %s : %s\n",
    "        end\n",
    "    end\n",
    "end\n"),
    test$result, inf_set, test$sensitivity, test$sensitivity,
    spec_pos, spec_neg)
}

# Known present: individual `i` is caught somewhere in [t, window end], so it
# must have been alive at `t`.
#
# The lookahead is bounded by the window, and that bound is not optional. A
# capture after `t + M` is outside the block being scored, and conditioning on
# it deletes a death branch that carries real probability mass: a lost
# positive contribution, which no reweighting restores. (Score one step whose
# observation is "no capture" with P(alive) = 0.4 and P(capture | alive) = 0.5:
# the honest answer is 0.4*0.5 + 0.6 = 0.8, and forcing survival on the strength
# of a later capture leaves 0.2.) So condition only on y_{(t+1):(t+M)}.
#
# `_et_seen_within` is an n_timepoints x n_individuals Bool: is there a capture
# at some u in [t, t + M - 1]? Precomputed per M, so the constraint is a lookup.
cr_known_present <- function(cr, M) {
  et_julia(sprintf("_et_seen_within_%s_M%s[t, i]", cr$caught, julia_int(M)))
}

cr_seen_within_src <- function(cr, M) {
  sprintf(paste0(
    "const _et_seen_within_%s_M%s = let y = DATA.%s, T = DATA.n_timepoints\n",
    "    [any(==(1), @view y[t:min(t + %s - 1, T), i])\n",
    "     for t in 1:T, i in 1:DATA.n_individuals]\n",
    "end"), cr$caught, julia_int(M), cr$caught, julia_int(M))
}

#' A leave-future-out scoring specification.
#'
#' R mirror of ET's `LFOSpec`. `cell_logdensity` (and `constrain`, if given) are
#' written directly in Julia via [et_julia()]: these run inside the forward
#' simulation, not inside the gradient, so there is no autodiff subset to
#' respect and nothing is gained by adding a transpiled R form. `model` in every
#' signature below is the fitted parameter NamedTuple (what `fit$draws` rows
#' become one at a time); `data` is the (possibly truncated) `DATA` object;
#' state codes are 1-based positions in the model's state vector.
#'
#' @param cell_logdensity Either an [et_capture_recapture()] (the usual case)
#'   which also derives `known_present` for you, or an [et_julia()] wrapping a
#'   Julia function `(model, data, X, i, t) -> Float64`: individual `i`'s
#'   observation log-density at `t` given the simulated states `X`. Return
#'   `-Inf` for an inadmissible trajectory; it is charged to that individual's
#'   own cell, never to the whole window.
#' @param truncation An [et_lfo_truncation()].
#' @param is_informative Optional [et_julia()] wrapping
#'   `(data, i, t) -> Bool`: whether that cell's density actually depends on the
#'   latent state. Strongly recommended, see `n_informative()` in ET's
#'   `lfo.jl`: on real data this can be a small fraction of the nominal cell
#'   count, which is the difference between "no granularity discriminates much"
#'   and a real comparison.
#' @param constrain_survival Simulate from the survival-constrained proposal
#'   instead of the raw dynamics: never enter the absorbing state when the data
#'   prove the individual was still present, and accumulate the importance
#'   weight that converts back to the true kernel. `TRUE` corrects the wasted
#'   draws; `FALSE` is the unconstrained estimator, kept so the two can be run
#'   on identical draws and compared.
#'
#'   **This addresses the absorbing state only, which is narrower than the
#'   problem it belongs to.** A draw scores `-Inf` whenever the simulated
#'   trajectory contradicts an observation, and an absorbing state is merely the
#'   worst case: once entered it cannot be left, so one bad step condemns the
#'   whole draw. A transient state that the observation also rules out — a
#'   recovered individual that later tests positive, say — produces the same
#'   `-Inf`, and this correction does nothing about it, because "do not enter
#'   state k at time t" is only a well-defined restriction with a computable
#'   weight when k is absorbing. Where transient states are the binding
#'   constraint, reach for a finer granularity (so the loss is charged to one
#'   cell rather than the window) or an exact forward recursion.
#'
#'   With an [et_capture_recapture()] score this needs nothing further: who was
#'   alive is read off the same capture matrix, looking ahead **only to the end
#'   of the scored window**. That bound matters. A capture after `t + M` is
#'   outside the block being scored, and conditioning on it deletes a death
#'   branch carrying real probability mass: a lost positive contribution that
#'   no reweighting restores. Otherwise supply `known_present`, and bound it
#'   yourself.
#' @param known_present Optional [et_julia()] whose body is a function of `i`
#'   and `t` only: individual `i` is known to be alive at `t`. **`data` is not
#'   in scope**: ET calls this as `known_present(i, t)`, so reach arrays through
#'   the generated module's own `DATA`, e.g. `t <= DATA.last_seen[i]`. Ignored
#'   unless `constrain_survival` is `TRUE`.
#' @param n_sim Forward trajectories per posterior draw per window.
#' @param seed Base RNG seed for the forward simulation.
#' @return An object of class `et_lfo_spec`.
#' @export
et_lfo_spec <- function(cell_logdensity, truncation, is_informative = NULL,
                        constrain_survival = FALSE, known_present = NULL,
                        n_sim = 1L, seed = 13L) {
  cr <- if (inherits(cell_logdensity, "et_capture_recapture")) cell_logdensity
  # A capture-recapture constraint depends on the window length, so it is built
  # by et_lfo_cv() once M is known -- see cr_known_present().
  if (!is.null(cr)) cell_logdensity <- cr_cell_logdensity(cr)
  if (!inherits(cell_logdensity, "et_julia")) {
    stop("et_lfo_spec(): `cell_logdensity` must come from et_capture_recapture() ",
         "or et_julia().", call. = FALSE)
  }
  if (isTRUE(constrain_survival) && is.null(cr) && is.null(known_present)) {
    stop("et_lfo_spec(): constrain_survival = TRUE needs to know who was alive. ",
         "Score with et_capture_recapture(), which derives it, or pass ",
         "`known_present`.", call. = FALSE)
  }
  if (!isTRUE(constrain_survival)) known_present <- NULL
  if (!inherits(truncation, "et_lfo_truncation")) {
    stop("et_lfo_spec(): `truncation` must come from et_lfo_truncation().", call. = FALSE)
  }
  if (!is.null(is_informative) && !inherits(is_informative, "et_julia")) {
    stop("et_lfo_spec(): `is_informative` must come from et_julia().", call. = FALSE)
  }
  if (!is.null(known_present) && !inherits(known_present, "et_julia")) {
    stop("et_lfo_spec(): `known_present` must come from et_julia().", call. = FALSE)
  }
  structure(list(cell_logdensity = cell_logdensity, truncation = truncation,
                 is_informative = is_informative, known_present = known_present,
                 capture_recapture = cr,
                 constrain_survival = isTRUE(constrain_survival),
                 n_sim = as.integer(n_sim), seed = as.integer(seed)),
            class = "et_lfo_spec")
}

#' @export
print.et_lfo_spec <- function(x, ...) {
  cat("<et_lfo_spec>  n_sim=", x$n_sim,
      if (isTRUE(x$constrain_survival)) "  survival-constrained"
      else "  naive forward simulation",
      "\n", sep = "")
  invisible(x)
}

# `LFOSpec(...)` call text. `known_present` expands to the co-gated
# (constrain, survival_weight) pair from `survival_constrained`, ET refuses
# one without the other (lfo_run.jl), so the R surface never offers them apart.
#
# `fit` is a placeholder: `LFOSpec` requires the field, but the module's own
# et_lfo_cv() (et_lfo_inject) never calls it: it rebuilds a second LFOSpec
# with the same cell_logdensity/plan/etc and its own `fit` closure over
# et_lfo_fit(), which is the one that actually refits this model per cutoff.
lfo_spec_src <- function(spec, plan_expr) {
  parts <- c(
    "fit = (train, t) -> error(\"unused placeholder: et_lfo_cv() supplies its own fit\")",
    sprintf("cell_logdensity = (model, data, X, i, t) -> begin\n%s\nend",
            indent(spec$cell_logdensity$src)),
    sprintf("plan = %s", plan_expr),
    sprintf("n_sim = %s", julia_int(spec$n_sim)),
    sprintf("seed = %s", julia_int(spec$seed))
  )
  if (!is.null(spec$is_informative)) {
    parts <- c(parts, sprintf(
      "is_informative = (data, i, t) -> begin\n%s\nend",
      indent(spec$is_informative$src)))
  }
  # The pair is destructured before the call: `a, b = f(...)` is an assignment,
  # and an assignment is not a legal keyword argument, so it cannot be inlined
  # into LFOSpec(...). Both halves then go in together, which is also the only
  # shape ET accepts, since it rejects one without the other.
  prelude <- ""
  if (!is.null(spec$known_present)) {
    prelude <- sprintf(paste0(
      "constrain, survival_weight = survival_constrained((i, t) -> begin\n",
      "%s\nend)\n"), indent(spec$known_present$src))
    parts <- c(parts, "constrain = constrain", "survival_weight = survival_weight")
  }
  paste0(prelude,
         sprintf("spec = LFOSpec(\n%s,\n)", indent(paste(parts, collapse = ",\n"))))
}

#' Run leave-future-out cross-validation.
#'
#' Refits the model at each cutoff `L:stride:(n_timepoints - M)`, forward-
#' simulates `M` steps and scores the held-out observations under every
#' requested granularity from the same fit and the same forward trajectories --
#' so a difference between granularities is never an artefact of separately-
#' drawn simulations. See ET's `lfo_cv` (lfo_run.jl) for the full contract.
#'
#' This always does EXACT refitting, with no importance-sampling reuse of an
#' earlier cutoff's posterior: every window's weights are uniform. The only
#' importance sampling here is the *inner* kind, correcting a guided forward
#' simulation, and it is switched on by `constrain_survival`. The distinction
#' matters: because the two are normalised differently, an outer weight ratio
#' is genuinely constant across draws and self-normalises, while an inner weight
#' must be divided by the number of draws. Self-normalising an inner weight
#' silently inflates the score.
#'
#' @param model An [et_model()], the same model passed to [et_sample()].
#' @param spec An [et_lfo_spec()].
#' @param L Minimum training window.
#' @param M Steps ahead scored per window.
#' @param granularity One or more of `"pointwise"`, `"joint"`, `"by_group"`.
#'   `"by_group"` uses the model data's own `group`.
#' @param stride Cutoff spacing (default every timepoint).
#' @param blocks Sampler blocks, as for [et_sample()].
#' @param n_sweeps,n_burn,n_adapts,adtype Passed to the refit at each cutoff, as
#'   for [et_sample()].
#' @param cache Optional directory to cache each cutoff's fit. Keyed on the
#'   cutoff only (not on granularity/scorer), so adding a scoring arm later
#'   costs seconds, not a re-run of every fit.
#' @param quiet Suppress progress messages.
#' @return An object of class `et_lfo_result`. `$granularities` names what was
#'   scored; use [et_lfo_elpd()] and [et_lfo_compare()] to read it.
#' @export
et_lfo_cv <- function(model, spec, L, M, granularity = "pointwise", stride = 1L,
                      blocks = list(), n_sweeps = 1000, n_burn = 0, n_adapts = 0,
                      adtype = "forwarddiff", cache = NULL, quiet = FALSE) {
  et_require_session()
  if (!inherits(model, "et_model")) {
    stop("et_lfo_cv(): `model` must come from et_model().", call. = FALSE)
  }
  if (!inherits(spec, "et_lfo_spec")) {
    stop("et_lfo_cv(): `spec` must come from et_lfo_spec().", call. = FALSE)
  }
  granularity <- match.arg(granularity, c("pointwise", "joint", "by_group"),
                           several.ok = TRUE)
  if ("by_group" %in% granularity && is.null(model$data$group)) {
    stop("et_lfo_cv(): granularity 'by_group' needs et_data(..., group=).",
         call. = FALSE)
  }

  gen <- et_julia_source(model, blocks)
  mod <- et_load_module(gen)
  JuliaCall::julia_command(sprintf(
    "Base.include_string(%s, \"using EpidemicTrajectories: truncation, LFOSpec, lfo_cv, Pointwise, Joint, ByGroup, survival_constrained\")",
    mod))

  et_lfo_inject(mod, model)   # idempotent: defines et_lfo_fit/et_lfo_cv once

  # `constrain_survival` on a model with no absorbing state would be a silent
  # no-op: ET's survival_weight derives the absorbing state itself and returns
  # an empty weight when there is none, so the run would look constrained and
  # correct nothing. Refuse instead.
  if (isTRUE(spec$constrain_survival)) {
    has_abs <- JuliaCall::julia_eval(sprintf(
      "Base.include_string(%s, \"EpidemicTrajectories._absorbing_state(DATA) !== nothing\")",
      mod))
    if (!isTRUE(has_abs)) {
      stop("et_lfo_cv(): constrain_survival = TRUE, but this model's ",
           "transitions make no state absorbing, so there is nothing the ",
           "constraint can forbid and it would silently correct nothing. Use ",
           "constrain_survival = FALSE, or check the transitions.",
           call. = FALSE)
    }
  }

  # The states a capture contradicts, bound once per module.
  if (!is.null(spec$capture_recapture) &&
      !is.null(spec$capture_recapture$unobservable)) {
    un_sym <- paste0(mod, "_lfo_unobs_src")
    JuliaCall::julia_assign(un_sym, sprintf(
      "if !isdefined(@__MODULE__, :_et_unobservable)\n%s\nend",
      cr_unobservable_src(spec$capture_recapture)))
    JuliaCall::julia_command(sprintf("Base.include_string(%s, Main.%s)",
                                     mod, un_sym))
  }

  gran_src <- paste(vapply(granularity, function(g) switch(g,
    pointwise = "Pointwise()", joint = "Joint()",
    by_group  = "ByGroup(DATA.group)"), character(1)), collapse = ", ")
  cache_src <- if (is.null(cache)) "nothing" else
    julia_string(julia_path(normalizePath(cache, mustWork = FALSE)))

  # A capture-recapture constraint is derived from the capture matrix rather
  # than asked for, so the constraint and the score cannot disagree about who
  # was alive -- and it depends on M, because the lookahead must stop at the end
  # of the scored block. Bound once per (matrix, M).
  cr <- spec$capture_recapture
  if (!is.null(cr) && isTRUE(spec$constrain_survival)) {
    spec$known_present <- cr_known_present(cr, M)
    sw_sym <- paste0(mod, "_lfo_seenwithin_src")
    JuliaCall::julia_assign(sw_sym, sprintf(
      "if !isdefined(@__MODULE__, :_et_seen_within_%s_M%s)\n%s\nend",
      cr$caught, julia_int(M), cr_seen_within_src(cr, M)))
    JuliaCall::julia_command(sprintf("Base.include_string(%s, Main.%s)",
                                     mod, sw_sym))
  }

  # one Julia expression: build the plan and the spec (both reference names --
  # `survival_constrained`, `LFOSpec`, that only resolve inside the module),
  # then call the module's own et_lfo_cv so the fit closure refits this model.
  # `spec_src` may be two statements (the survival pair is destructured before
  # the LFOSpec call), so it is spliced at statement level, not into `spec = `.
  body <- sprintf(paste0(
    "let\n",
    "    plan = %s\n",
    "%s\n",
    "    et_lfo_cv(spec; L=%s, M=%s, granularity=(%s,), stride=%s, cache=%s,\n",
    "              verbose=%s, n_sweeps=%s, n_burn=%s, n_adapts=%s, adtype=%s)\n",
    "end"),
    truncation_src(spec$truncation),
    indent(lfo_spec_src(spec, plan_expr = "plan")),
    julia_int(L), julia_int(M), gran_src, julia_int(stride), cache_src,
    if (quiet) "false" else "true",
    julia_int(n_sweeps), julia_int(n_burn), julia_int(n_adapts),
    adtype_to_julia(adtype))

  # The result stays in Julia, under a name of its own, not round-tripped
  # through R, which has no faithful representation of an LFOResult (a Dict of
  # WindowResult structs) to hand back to a later julia_assign(). et_lfo_elpd()
  # and et_lfo_compare() reference it there by name, exactly as et_residuals()
  # leaves the fitted module itself in Julia rather than returning it.
  res_sym <- paste0(mod, "_lfo_res_", length(.et_state$lfo_results %||% character()) + 1L)
  src_sym <- paste0(mod, "_lfo_call_src")
  JuliaCall::julia_assign(src_sym, body)
  JuliaCall::julia_command(sprintf(
    "Base.include_string(%s, \"const %s = \" * Main.%s)", mod, res_sym, src_sym))
  # `const %s = ...` above only works the FIRST time a given name is bound in a
  # module (Julia forbids re-binding a `const`); res_sym is suffixed by an
  # ever-growing counter precisely so a second et_lfo_cv() call in the same
  # session never collides with the first.
  .et_state$lfo_results <- c(.et_state$lfo_results %||% character(), res_sym)
  structure(list(sym = res_sym, module = mod, granularities = granularity,
                 L = L, M = M), class = "et_lfo_result")
}

# Inject, into the loaded module, a fit function usable as `LFOSpec.fit` and a
# thin `et_lfo_cv` wrapper that supplies it. Mirrors et_run() (codegen.R)
# exactly: same et_prepare!, same SPL, same INIT_PARS, except it takes
# `data` as an argument instead of closing over the module constant `DATA`,
# which is the one thing a normal fit never needs and LFO always does.
et_lfo_inject <- function(mod, model) {
  nm <- paste0(mod, "_lfo_inject_src")
  JuliaCall::julia_assign(nm, lfo_inject_src(model))
  JuliaCall::julia_command(sprintf("Base.include_string(%s, Main.%s)", mod, nm))
}

# The injected text, as a pure function of the model with no Julia session, so the
# generated source can be inspected and parse-checked on its own.
lfo_inject_src <- function(model) {
  has_obs <- !is.null(model$data$observation_process) ||
    !is.null(model$data$observation_weight)
  paste0(
    "if !isdefined(@__MODULE__, :et_lfo_fit)\n",
    "function et_lfo_fit(data; n_sweeps, n_burn=0, n_adapts=0, seed=1,\n",
    "                     adtype=ADTypes.AutoForwardDiff())\n",
    "    X0 = fill(1, data.n_timepoints, data.n_individuals)\n",
    "    reset_aggregates!(data)\n",
    "    apply_derived_summaries!(et_params(INIT_PARS), data, X0)\n",
    "    m = et_the_model(data, LOGLIK", if (has_obs) ", OBSLOGLIK" else "", ")\n",
    "    init = (; X=X0, INIT_PARS...)\n",
    "    # Stepped BY HAND, exactly as et_collect_run() is: `X` must be RETAINED\n",
    "    # here (LFO needs every draw's whole trajectory to forward-simulate from),\n",
    "    # which is the opposite of a normal fit's `save_states=(X=:buffer,)`. The\n",
    "    # bundled chain drops per-sweep `X` unconditionally, so the raw\n",
    "    # transitions from `step` are read directly instead.\n",
    "    rng = StableRNG(seed)\n",
    "    t, state = AbstractMCMC.step(rng, m, SPL; init=init, adtype=adtype,\n",
    "                                 n_adapts=n_adapts)\n",
    "    for _ in 1:n_burn\n",
    "        t, state = AbstractMCMC.step(rng, m, SPL, state; n_adapts=n_adapts)\n",
    "    end\n",
    "    draws = Vector{Any}(undef, n_sweeps)\n",
    "    Xs = Vector{Matrix{Int}}(undef, n_sweeps)\n",
    "    draws[1] = et_params(t); Xs[1] = copy(t.X)\n",
    "    for k in 2:n_sweeps\n",
    "        t, state = AbstractMCMC.step(rng, m, SPL, state; n_adapts=n_adapts)\n",
    "        draws[k] = et_params(t); Xs[k] = copy(t.X)\n",
    "    end\n",
    "    (draws, Xs)\n",
    "end\n",
    "# `data` truncated per cutoff, model refit fresh each time -- see truncate.jl\n",
    "# for why shortening n_timepoints alone is not enough.\n",
    "function et_lfo_cv(spec; L, M, granularity, stride=1, cache=nothing,\n",
    "                    verbose=true, n_sweeps, n_burn=0, n_adapts=0, adtype)\n",
    "    fitfn = (train, t) -> et_lfo_fit(train; n_sweeps=n_sweeps, n_burn=n_burn,\n",
    "                                      n_adapts=n_adapts, seed=1000 + t, adtype=adtype)\n",
    "    spec2 = LFOSpec(fit=fitfn, cell_logdensity=spec.cell_logdensity,\n",
    "                    plan=spec.plan, is_informative=spec.is_informative,\n",
    "                    constrain=spec.constrain, survival_weight=spec.survival_weight,\n",
    "                    n_sim=spec.n_sim, seed=spec.seed)\n",
    "    lfo_cv(spec2, DATA; L=L, M=M, granularity=granularity, stride=stride,\n",
    "          cache=cache, verbose=verbose)\n",
    "end\n",
    "end\n")
}

#' Total ELPD from a leave-future-out sweep.
#'
#' @param res An [et_lfo_cv()] result.
#' @param granularity Which granularity to total; defaults to the first one run.
#' @return A single number.
#' @export
et_lfo_elpd <- function(res, granularity = res$granularities[1]) {
  if (!inherits(res, "et_lfo_result")) {
    stop("et_lfo_elpd(): `res` must come from et_lfo_cv().", call. = FALSE)
  }
  granularity <- match.arg(granularity, res$granularities)
  JuliaCall::julia_eval(sprintf("Base.include_string(%s, \"elpd(%s, %s)\")",
                                res$module, res$sym, julia_symbol(granularity)))
}

#' Per-window diagnostics from a leave-future-out sweep.
#'
#' One row per cutoff: the window's score under each granularity, how many of
#' its posterior draws produced a finite score, and (if the spec supplied
#' `is_informative`) how many cells were informative. `n_finite / n_draws`
#' collapsing toward 0 as the horizon grows is the signature of the dead-
#' animal degeneracy this package's survival-constrained IS exists to fix --
#' plot it against [et_lfo_cv()]'s naive and constrained results side by side.
#'
#' @param res An [et_lfo_cv()] result.
#' @return A data frame: `cutoff`, `n_draws`, `n_informative`, and
#'   `elpd_<granularity>` / `n_finite_<granularity>` for each granularity run.
#' @export
et_lfo_windows <- function(res) {
  if (!inherits(res, "et_lfo_result")) {
    stop("et_lfo_windows(): `res` must come from et_lfo_cv().", call. = FALSE)
  }
  gnames <- julia_symbol_vector(res$granularities)
  body <- sprintf(paste0(
    "let\n",
    "    ws = %s.windows\n",
    "    Dict(\n",
    "        \"cutoff\" => [w.cutoff for w in ws],\n",
    "        \"n_draws\" => [w.n_draws for w in ws],\n",
    "        \"n_informative\" => [w.n_informative === nothing ? missing : w.n_informative for w in ws],\n",
    "        [string(\"elpd_\", g) => [w.elpd[g] for w in ws] for g in %s]...,\n",
    "        [string(\"n_finite_\", g) => [w.n_finite[g] for w in ws] for g in %s]...,\n",
    "    )\n",
    "end"), res$sym, gnames, gnames)
  src_sym <- paste0(res$module, "_lfo_windows_src")
  JuliaCall::julia_assign(src_sym, body)
  raw <- JuliaCall::julia_eval(sprintf("Base.include_string(%s, Main.%s)",
                                       res$module, src_sym))
  as.data.frame(raw, stringsAsFactors = FALSE)
}

#' Compare two leave-future-out sweeps.
#'
#' `a` minus `b`, on the cutoffs they share, under one granularity. Refuses
#' cross-granularity totals (a pointwise total sums `N*M` densities per window,
#' a joint total sums one; see ET's `compare`).
#'
#' @param a,b [et_lfo_cv()] results.
#' @param granularity Which granularity to compare.
#' @return A list with `diff`, `se` (windows treated as independent -- optimistic
#'   once `M > 1`, since windows then overlap), `se_indep` (every `M`-th window;
#'   the conservative reading), and `n_windows`.
#' @export
et_lfo_compare <- function(a, b, granularity = a$granularities[1]) {
  if (!inherits(a, "et_lfo_result") || !inherits(b, "et_lfo_result")) {
    stop("et_lfo_compare(): `a` and `b` must come from et_lfo_cv().", call. = FALSE)
  }
  if (a$module != b$module) {
    stop("et_lfo_compare(): both results must come from models generated in ",
         "the SAME Julia module -- compare results from the same et_model(), ",
         "or from two et_lfo_spec()s run through et_lfo_cv() on it.", call. = FALSE)
  }
  out <- JuliaCall::julia_eval(sprintf(
    "Base.include_string(%s, \"compare(%s, %s; granularity=%s)\")",
    a$module, a$sym, b$sym, julia_symbol(granularity)))
  out$granularity <- granularity
  out
}

#' @export
print.et_lfo_result <- function(x, ...) {
  cat("<et_lfo_result>  L=", x$L, " M=", x$M, "  granularity: ",
      paste(x$granularities, collapse = ", "), "\n", sep = "")
  for (g in x$granularities) {
    cat("  ", g, ": elpd = ", format(et_lfo_elpd(x, g), digits = 6), "\n", sep = "")
  }
  invisible(x)
}

#' Convergence diagnostics for a leave-future-out sweep's fits.
#'
#' One row per cutoff and parameter: effective sample size, its Monte Carlo
#' standard error, and rhat.
#'
#' **ESS is the number to read here, and rhat is the weaker of the two.** A
#' sweep fits one chain per cutoff, so there is no between-chain rhat; what is
#' reported is computed over the split halves of that single chain, which
#' catches a drifting chain but not one stuck in a single mode. ESS also
#' answers the question the caller actually has -- whether enough draws were
#' retained -- and it is sensitive to the autocorrelation a Gibbs sweep over a
#' latent state produces.
#'
#' Read the WORST cell rather than an average. One badly mixed parameter at one
#' cutoff invalidates that window's score, and a mean over cutoffs hides
#' exactly that.
#'
#' The draws are read from the fits `et_lfo_cv()` cached, not from its result:
#' a scored window keeps its scores and releases the parameter draws. So this
#' needs `cache=` to have been set, and to be called before that directory is
#' cleaned up.
#'
#' @param cache The cache directory given to [et_lfo_cv()].
#' @return A data frame: `cutoff`, `parameter`, `ess`, `mcse`, `rhat`,
#'   `n_draws`, ordered worst ESS first.
#' @export
et_lfo_diagnostics <- function(cache) {
  et_require_session()
  if (!dir.exists(cache)) {
    stop("et_lfo_diagnostics(): no such cache directory: ", cache, call. = FALSE)
  }
  files <- sort(list.files(cache, pattern = "^fit_t[0-9]+[.]jls$",
                           full.names = TRUE))
  if (!length(files)) {
    stop("et_lfo_diagnostics(): no fit_t*.jls in ", cache,
         " -- was et_lfo_cv() given cache=, and the directory kept?",
         call. = FALSE)
  }
  out <- do.call(rbind, lapply(files, function(f) {
    src <- sprintf(lfo_diagnostics_src(),
                   julia_string(julia_path(normalizePath(f, mustWork = TRUE))))
    d <- as.data.frame(JuliaCall::julia_eval(src), stringsAsFactors = FALSE)
    if (!nrow(d)) return(NULL)
    cbind(cutoff = as.integer(sub("^.*fit_t0*([0-9]+)[.]jls$", "\1", f)), d,
          stringsAsFactors = FALSE)
  }))
  out[order(out$ess), ]
}

# The Julia half, as a pure function of nothing so it can be parse-checked
# without a session. `%s` is the cache file path.
#
# FlexiChains is reached THROUGH PracticalBayes rather than imported. It is
# PracticalBayes's dependency, not this project's, so `import FlexiChains`
# fails with "not found in current path" even though the package is resolved
# and already loaded -- and adding a direct dependency for four functions would
# be a heavier change than asking the module that owns it.
#
# ess/rhat/mcse each return a FlexiSummary wrapping a 3-D array, hence `only`.
lfo_diagnostics_src <- function() paste0(
  "let FC = @eval(PracticalBayes, FlexiChains)\n",
  "    draws, _ = open(deserialize, %s)\n",
  "    S = length(draws)\n",
  # Scalar rates only: a vector parameter would need flattening into one column
  # per element, and nothing in these models has one.
  "    nms = [k for (k, v) in pairs(draws[1]) if v isa Real]\n",
  "    chain = FC.FlexiChain{Symbol}(S, 1, Dict(\n",
  "        FC.Parameter(n) => reshape(Float64[getproperty(d, n) for d in draws], S, 1)\n",
  "        for n in nms))\n",
  "    e = FC.ess(chain); m = FC.mcse(chain); r = FC.rhat(chain)\n",
  "    Dict(\n",
  "        \"parameter\" => String.(nms),\n",
  "        \"ess\"     => [Float64(only(e[FC.Parameter(n)])) for n in nms],\n",
  "        \"mcse\"    => [Float64(only(m[FC.Parameter(n)])) for n in nms],\n",
  "        \"rhat\"    => [Float64(only(r[FC.Parameter(n)])) for n in nms],\n",
  "        \"n_draws\" => fill(S, length(nms)),\n",
  "    )\n",
  "end")
