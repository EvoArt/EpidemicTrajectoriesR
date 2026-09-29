# Leave-future-out cross-validation: a wrapper over ET's truncation /
# LFOSpec / lfo_cv / compare. Nothing statistical is reimplemented here.
#
# et_run() closes over the module-level `const DATA`, since a normal fit only
# sees the whole series. LFO needs the same compiled model refit against a
# truncated copy at every cutoff, so lfo_inject_src() emits a second entry
# point taking `data` as an argument and rebuilding every data-dependent piece
# from it. Injected with Base.include_string(), so it stays out of
# et_julia_source() and never perturbs the module-name hash.

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
#' @param custom Named list of [et_julia()] rules, one per `extras` name, for
#'   anything the four fixed rules do not cover. Each body is a Julia function
#'   of `v` (the full-series value) and `cutoff`, returning the value the fit
#'   at that cutoff should see. The generated module's `DATA` is in scope.
#' @param strict Error on an undeclared, time-shaped extra (the default and the
#'   safe choice).
#' @param known_from Optional name of an `extras` vector giving, per individual,
#'   the first time the study knew it existed. Needed when a sampling period
#'   opens before the individual is first observed (an animal born into a
#'   watched group but not caught until later). An individual not yet known at
#'   a cutoff is left out of that cutoff's fit, and is neither scored nor
#'   simulated in its window. Defaults to the start of each sampling period.
#' @return An object of class `et_lfo_truncation`.
#' @export
et_lfo_truncation <- function(clamp = character(), filter = character(),
                              copy = character(), keep = character(),
                              custom = list(), strict = TRUE,
                              known_from = NULL) {
  if (length(custom)) {
    if (is.null(names(custom)) || any(!nzchar(names(custom))) ||
        !all(vapply(custom, inherits, logical(1), "et_julia"))) {
      stop("et_lfo_truncation(): `custom` must be a named list of et_julia() ",
           "rules, one per extras name.", call. = FALSE)
    }
    check_julia_name(names(custom), "extras name")
  }
  if (!is.null(known_from)) {
    if (!is.character(known_from) || length(known_from) != 1L) {
      stop("et_lfo_truncation(): `known_from` must be one extras name.",
           call. = FALSE)
    }
    check_julia_name(known_from, "extras name")
  }
  structure(list(clamp = clamp, filter = filter, copy = copy, keep = keep,
                 custom = custom, strict = isTRUE(strict),
                 known_from = known_from),
            class = "et_lfo_truncation")
}

#' @export
print.et_lfo_truncation <- function(x, ...) {
  cat("<et_lfo_truncation>\n")
  for (r in c("clamp", "filter", "copy", "keep")) {
    if (length(x[[r]])) cat("  ", r, ": ", paste(x[[r]], collapse = ", "), "\n", sep = "")
  }
  if (length(x$custom)) cat("  custom: ", paste(names(x$custom), collapse = ", "), "\n", sep = "")
  if (!is.null(x$known_from)) cat("  known from: ", x$known_from, "\n", sep = "")
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
  if (length(plan$custom)) {
    rules <- vapply(names(plan$custom), function(nm) sprintf(
      "%s => ((v, cutoff) -> begin\n%s\nend)", julia_symbol(nm),
      indent(plan$custom[[nm]]$src)), character(1))
    parts <- c(parts, sprintf("custom=(%s,)", paste(rules, collapse = ", ")))
  }
  parts <- c(parts, sprintf("strict=%s", if (plan$strict) "true" else "false"))
  if (!is.null(plan$known_from)) {
    parts <- c(parts, sprintf("known_from=%s", julia_symbol(plan$known_from)))
  }
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
# ET passes the window's last timepoint as `last_t`, so the bound holds whether
# or not windows overlap.
cr_known_present <- function(cr) {
  caught <- if (is.character(cr)) cr else cr$caught
  et_julia(sprintf("any(==(1), @view DATA.%s[t:last_t, i])", caught))
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
#' @param cell_logdensity How a simulated state is scored against an
#'   observation. `NULL` (the default) uses the model's own observation weight,
#'   which keeps the score and the fitted model on one observation model and is
#'   what `method = "exact_hmm"` scores with too. Otherwise an
#'   [et_capture_recapture()], which also derives `known_present` for you, or
#'   an [et_julia()] wrapping a Julia function `(model, data, X, i, t) ->
#'   Float64`: individual `i`'s observation log-density at `t` given the
#'   simulated states `X`. Return `-Inf` for an inadmissible trajectory; it is
#'   charged to that individual's own cell, never to the whole window.
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
#'   of the scored window**, whichever step of the window is being simulated.
#'   That bound matters. A capture after the window's last occasion is outside
#'   the block being scored, and conditioning on it deletes a death branch
#'   carrying real probability mass: a lost positive contribution that no
#'   reweighting restores. Otherwise supply `known_present`.
#' @param known_present Who is known to be alive, for `constrain_survival`:
#'   the name of a 0/1 capture matrix in `extras` (alive at `t` if caught at
#'   some occasion from `t` to the end of the scored window), or an
#'   [et_julia()] whose body is a function of `i`, `t` and `last_t`, the
#'   window's last occasion, e.g. `any(==(1), @view DATA.caught[t:last_t, i])`.
#'   `data` is not in scope, so reach arrays through the generated module's
#'   own `DATA`. Never look past `last_t`. Ignored unless `constrain_survival`
#'   is `TRUE`.
#' @param n_sim Forward trajectories per posterior draw per window.
#' @param seed Base RNG seed for the forward simulation.
#' @param observation_guided Simulate from the observation-guided proposal:
#'   each animal's next state is drawn in proportion to its move probability
#'   times how likely its own observations over the rest of the window are
#'   from that state, and the importance weight is carried with the score.
#'   The survival constraint is the special case that looks only at whether
#'   the animal is alive; this one also avoids transient states the data rule
#'   out, such as a seronegative animal converting before a negative test.
#'   Not with `constrain_survival`. For a coupled model it corrects only the
#'   joint score, so et_lfo_cv() then takes `granularity = "joint"` alone.
#' @param simulate_entrants Whether individuals entering after the cutoff are
#'   simulated into the forecast, drawn from the starting state at their entry
#'   time, so that they can affect the cohort's dynamics. That uses entry times
#'   only later data reveal; `FALSE` fixes the forecast population at the
#'   cutoff. It matters only for coupled models.
#' @param method How each posterior draw's predictive density is computed.
#'
#'   * `"simulate"` (the default): forward-simulate from the draw's sampled
#'     state at the cutoff and score with `cell_logdensity`. Works for any
#'     model.
#'   * `"exact_hmm"`: for individuals independent given the parameters, sum
#'     each one's state at the cutoff and its future path out exactly, by the
#'     forward algorithm, and score with the model's own observation process.
#'     No forward simulation, so no `cell_logdensity`, `n_sim` or survival
#'     constraint; the only Monte Carlo error left is the posterior draws'.
#'     Scores whole forecast blocks, so for `M > 1` it takes the `"joint"` and
#'     `"by_individual"` granularities.
#'
#'   The method and the model's `likelihood` are separate choices: an
#'   augmented (iFFBS) fit can be scored exactly, and a marginal fit can be
#'   scored by simulation, its trajectories then drawn afterwards from their
#'   exact posterior given the training data.
#' @return An object of class `et_lfo_spec`.
#' @export
et_lfo_spec <- function(cell_logdensity = NULL, truncation, is_informative = NULL,
                        constrain_survival = FALSE, known_present = NULL,
                        n_sim = 1L, seed = 13L, simulate_entrants = TRUE,
                        observation_guided = FALSE,
                        method = c("simulate", "exact_hmm")) {
  method <- match.arg(method)
  if (!inherits(truncation, "et_lfo_truncation")) {
    stop("et_lfo_spec(): `truncation` must come from et_lfo_truncation().", call. = FALSE)
  }
  if (!is.null(is_informative) && !inherits(is_informative, "et_julia")) {
    stop("et_lfo_spec(): `is_informative` must come from et_julia().", call. = FALSE)
  }
  if (method == "exact_hmm") {
    # Refuse rather than ignore: a run that looks corrected or simulated and is
    # neither is worse than an error.
    if (!is.null(cell_logdensity)) {
      stop("et_lfo_spec(): method = \"exact_hmm\" scores with the model's own ",
           "observation process; drop `cell_logdensity` so there is only one ",
           "observation model.", call. = FALSE)
    }
    if (isTRUE(constrain_survival) || !is.null(known_present) ||
        isTRUE(observation_guided)) {
      stop("et_lfo_spec(): method = \"exact_hmm\" integrates the future ",
           "exactly, so there is no forward proposal to constrain.", call. = FALSE)
    }
    if (as.integer(n_sim) != 1L) {
      stop("et_lfo_spec(): `n_sim` does not apply to method = \"exact_hmm\".",
           call. = FALSE)
    }
    return(structure(list(method = method, truncation = truncation,
                          is_informative = is_informative,
                          cell_logdensity = NULL, known_present = NULL,
                          capture_recapture = NULL, constrain_survival = FALSE,
                          n_sim = 1L, seed = as.integer(seed)),
                     class = "et_lfo_spec"))
  }
  cr <- if (inherits(cell_logdensity, "et_capture_recapture")) cell_logdensity
  # A capture-recapture score's constraint is built by et_lfo_cv() from the
  # same capture matrix -- see cr_known_present().
  if (!is.null(cr)) cell_logdensity <- cr_cell_logdensity(cr)
  if (is.character(known_present)) {
    check_julia_name(known_present, "extras name")
    known_present <- cr_known_present(known_present)
  }
  # NULL means the model's own observation model, filled in by et_lfo_cv()
  # once the model is known -- see own_cell_logdensity().
  if (!is.null(cell_logdensity) && !inherits(cell_logdensity, "et_julia")) {
    stop("et_lfo_spec(): `cell_logdensity` must come from et_capture_recapture() ",
         "or et_julia().", call. = FALSE)
  }
  if (isTRUE(constrain_survival) && is.null(cr) && is.null(known_present)) {
    stop("et_lfo_spec(): constrain_survival = TRUE needs to know who was alive. ",
         "Score with et_capture_recapture(), which derives it, or pass ",
         "`known_present`.", call. = FALSE)
  }
  if (isTRUE(constrain_survival) && isTRUE(observation_guided)) {
    stop("et_lfo_spec(): use observation_guided or constrain_survival, not both.",
         call. = FALSE)
  }
  if (!isTRUE(constrain_survival)) known_present <- NULL
  if (!is.null(known_present) && !inherits(known_present, "et_julia")) {
    stop("et_lfo_spec(): `known_present` must come from et_julia().", call. = FALSE)
  }
  structure(list(method = method,
                 cell_logdensity = cell_logdensity, truncation = truncation,
                 is_informative = is_informative, known_present = known_present,
                 capture_recapture = cr,
                 constrain_survival = isTRUE(constrain_survival),
                 simulate_entrants = isTRUE(simulate_entrants),
                 observation_guided = isTRUE(observation_guided),
                 n_sim = as.integer(n_sim), seed = as.integer(seed)),
            class = "et_lfo_spec")
}

#' @export
print.et_lfo_spec <- function(x, ...) {
  if (identical(x$method, "exact_hmm")) {
    cat("<et_lfo_spec>  exact: hidden states summed out by the forward algorithm\n")
    return(invisible(x))
  }
  cat("<et_lfo_spec>  n_sim=", x$n_sim,
      if (isTRUE(x$constrain_survival)) "  survival-constrained"
      else "  naive forward simulation",
      "\n", sep = "")
  invisible(x)
}

lfo_method <- function(spec) spec$method %||% "simulate"

# The model's own observation weight as a per-cell log density: the same
# generated function the fit and the exact scorer use.
own_cell_logdensity <- function(model) {
  d <- model$data
  if (!is.null(d$observation_weight)) {
    return(et_julia("log(data.observation_weight(model, data, X, i, t, X[t, i]))"))
  }
  if (!is.null(d$observation_process)) {
    return(et_julia("log(et_obs_process(model, data, X, i, t)[X[t, i]])"))
  }
  stop("et_lfo_cv(): the model has no observation process to score with; ",
       "give et_lfo_spec() a `cell_logdensity`.", call. = FALSE)
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
  exact <- lfo_method(spec) == "exact_hmm"
  parts <- c(
    "fit = (train, t) -> error(\"unused placeholder: et_lfo_cv() supplies its own fit\")",
    if (exact) "scorer = ExactHMM()"
    else sprintf("cell_logdensity = (model, data, X, i, t) -> begin\n%s\nend",
                 indent(spec$cell_logdensity$src)),
    sprintf("plan = %s", plan_expr),
    sprintf("n_sim = %s", julia_int(spec$n_sim)),
    sprintf("seed = %s", julia_int(spec$seed)),
    if (identical(spec$simulate_entrants, FALSE)) "entrants = false",
    if (isTRUE(spec$observation_guided)) "guide = true"
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
      "constrain, survival_weight = survival_constrained((i, t, last_t) -> begin\n",
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
#' @param granularity One or more of `"pointwise"`, `"by_individual"`, `"joint"`,
#'   `"by_group"`. `"by_individual"` scores each animal's full M-step history.
#'   `"by_group"` uses the model data's own `group`.
#' @param stride Cutoff spacing (default every timepoint).
#' @param blocks Sampler blocks, as for [et_sample()].
#' @param n_sweeps,n_burn,n_adapts,adtype Passed to the refit at each cutoff, as
#'   for [et_sample()].
#' @param cache Optional directory to cache each cutoff's fit. Keyed on the
#'   cutoff only (not on granularity/scorer), so adding a scoring arm later
#'   costs seconds, not a re-run of every fit.
#' @param x_init Optional `n_timepoints x n_individuals` integer matrix of state
#'   codes to start each cutoff's chain from. The default, all-susceptible, is a
#'   poor start for a transmission model: iFFBS resamples one individual at a
#'   time, so with nobody infected the only route in is the background hazard,
#'   and a small one leaves the chain effectively stuck. The transmission term
#'   is then multiplied by a zero infected count throughout, so its rate is
#'   unidentified — and two models differing only in that term fit identically,
#'   which looks like a working comparison that resolves nothing. Truncation
#'   keeps the full time dimension, so one matrix is the right shape at every
#'   cutoff.
#' @param cutoffs Optional explicit vector of cutoffs, overriding
#'   `L:stride:(n_timepoints - M)`. `L` and `stride` still describe the grid the
#'   cutoffs sit on, which a capture-recapture survival constraint uses to find
#'   the end of the window each step belongs to.
#' @param thin Keep every `thin`-th sweep after burn-in, so each cutoff's fit
#'   holds `n_sweeps %/% thin` draws. Every retained draw carries a whole
#'   `n_timepoints x n_individuals` trajectory, so on a large population this is
#'   what keeps the fit, and its cache, in memory.
#' @param fit_seed Base seed for the refits: the chain at cutoff `t` is seeded
#'   `fit_seed + t`. Independent chains at the same cutoffs need different
#'   values, and their own `cache` directories, since the cache is keyed on the
#'   cutoff alone.
#' @param quiet Suppress progress messages.
#' @return An object of class `et_lfo_result`. `$granularities` names what was
#'   scored; use [et_lfo_elpd()] and [et_lfo_compare()] to read it.
#' @export
et_lfo_cv <- function(model, spec, L, M, granularity = "pointwise", stride = 1L,
                      blocks = list(), n_sweeps = 1000, n_burn = 0, n_adapts = 0,
                      adtype = "forwarddiff", cache = NULL, x_init = NULL,
                      cutoffs = NULL, thin = 1L,
                      fit_seed = 1000L, quiet = FALSE) {
  et_require_session()
  if (!inherits(model, "et_model")) {
    stop("et_lfo_cv(): `model` must come from et_model().", call. = FALSE)
  }
  if (!inherits(spec, "et_lfo_spec")) {
    stop("et_lfo_cv(): `spec` must come from et_lfo_spec().", call. = FALSE)
  }
  granularity <- match.arg(granularity,
                           c("pointwise", "by_individual", "joint", "by_group"),
                           several.ok = TRUE)
  if ("by_group" %in% granularity && is.null(model$data$group)) {
    stop("et_lfo_cv(): granularity 'by_group' needs et_data(..., group=).",
         call. = FALSE)
  }
  if (lfo_method(spec) == "exact_hmm" && M > 1L &&
      any(!granularity %in% c("joint", "by_individual"))) {
    stop("et_lfo_cv(): method = \"exact_hmm\" scores whole forecast blocks, so ",
         "for M > 1 the granularity must be \"joint\" or \"by_individual\".",
         call. = FALSE)
  }
  if (is_marginal(model) && !is.null(x_init)) {
    stop("et_lfo_cv(): a marginal-likelihood model has no trajectory to start ",
         "from; drop `x_init`.", call. = FALSE)
  }
  if (lfo_method(spec) == "simulate" && is.null(spec$cell_logdensity)) {
    spec$cell_logdensity <- own_cell_logdensity(model)
  }

  gen <- et_julia_source(model, blocks)
  mod <- et_load_module(gen)
  JuliaCall::julia_command(sprintf(
    "Base.include_string(%s, \"using EpidemicTrajectories: truncation, LFOSpec, lfo_cv, Pointwise, ByIndividual, Joint, ByGroup, survival_constrained, ExactHMM, cell_elpd\")",
    mod))

  # Truncation keeps the full time dimension (it clamps sampling periods
  # instead), so one full-size matrix is the right shape at every cutoff.
  x_sym <- "nothing"
  if (!is.null(x_init)) {
    x_init <- as.matrix(x_init); storage.mode(x_init) <- "integer"
    check_x_init(x_init, model, "et_lfo_cv()")
    x_sym <- paste0("Main.", mod, "_lfo_xinit")
    JuliaCall::julia_assign(paste0(mod, "_lfo_xinit"), x_init)
    x_sym <- sprintf("Matrix{Int}(%s)", x_sym)
  }

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
    by_individual = "ByIndividual()",
    by_group  = "ByGroup(DATA.group)"), character(1)), collapse = ", ")
  cache_src <- if (is.null(cache)) "nothing" else
    julia_string(julia_path(normalizePath(cache, mustWork = FALSE)))

  # A capture-recapture constraint is derived from the capture matrix rather
  # than asked for, so the constraint and the score cannot disagree about who
  # was alive -- and it depends on M, because the lookahead must stop at the end
  # of the scored block. Bound once per (matrix, M).
  cr <- spec$capture_recapture
  if (!is.null(cr) && isTRUE(spec$constrain_survival)) {
    spec$known_present <- cr_known_present(cr)
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
    "              verbose=%s, n_sweeps=%s, n_burn=%s, n_adapts=%s, adtype=%s,\n",
    "              x_init=%s, fit_seed=%s, cutoffs=%s, thin=%s)\n",
    "end"),
    truncation_src(spec$truncation),
    indent(lfo_spec_src(spec, plan_expr = "plan")),
    julia_int(L), julia_int(M), gran_src, julia_int(stride), cache_src,
    if (quiet) "false" else "true",
    julia_int(n_sweeps), julia_int(n_burn), julia_int(n_adapts),
    adtype_to_julia(default_adtype(model, adtype)), x_sym, julia_int(fit_seed),
    if (is.null(cutoffs)) "nothing" else julia_vector(as.integer(cutoffs), "Int"),
    julia_int(thin))

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
                 L = L, M = M, method = lfo_method(spec),
                 likelihood = model$likelihood %||% "augmented"),
            class = "et_lfo_result")
}

# Inject, into the loaded module, a fit function usable as `LFOSpec.fit` and a
# thin `et_lfo_cv` wrapper that supplies it. Mirrors et_run() (codegen.R), with
# the same INIT_PARS and the same blocks, except that it takes `data` as an
# argument and builds the likelihoods and the sampler from it (et_loglik_for,
# et_latent_for, et_sampler_for) instead of using the constants built from
# `DATA`. Reusing those constants runs the latent update on the full series.
et_lfo_inject <- function(mod, model) {
  nm <- paste0(mod, "_lfo_inject_src")
  JuliaCall::julia_assign(nm, qualify_base(lfo_inject_src(model)))
  JuliaCall::julia_command(sprintf("Base.include_string(%s, Main.%s)", mod, nm))
}

# The injected text, as a pure function of the model with no Julia session, so the
# generated source can be inspected and parse-checked on its own.
lfo_inject_src <- function(model) {
  has_obs <- !is.null(model$data$observation_process) ||
    !is.null(model$data$observation_weight)
  paste0(
    "if !isdefined(@__MODULE__, :et_lfo_fit)\n",
    if (is_marginal(model)) lfo_fit_marginal_src() else lfo_fit_augmented_src(has_obs),
    "# `data` truncated per cutoff, model refit fresh each time -- see truncate.jl\n",
    "# for why shortening n_timepoints alone is not enough.\n",
    "function et_lfo_cv(spec; L, M, granularity, stride=1, cache=nothing,\n",
    "                    verbose=true, n_sweeps, n_burn=0, n_adapts=0, adtype,\n",
    "                    x_init=nothing, fit_seed=1000, cutoffs=nothing, thin=1)\n",
    "    fitfn = (train, t) -> et_lfo_fit(train; n_sweeps=n_sweeps, n_burn=n_burn,\n",
    "                                      n_adapts=n_adapts, seed=fit_seed + t, adtype=adtype,\n",
    "                                      x_init=x_init, thin=thin)\n",
    "    spec2 = LFOSpec(fit=fitfn, cell_logdensity=spec.cell_logdensity,\n",
    "                    plan=spec.plan, is_informative=spec.is_informative,\n",
    "                    constrain=spec.constrain, survival_weight=spec.survival_weight,\n",
    "                    n_sim=spec.n_sim, seed=spec.seed, scorer=spec.scorer,\n",
    "                    entrants=spec.entrants, guide=spec.guide)\n",
    "    lfo_cv(spec2, DATA; L=L, M=M, granularity=granularity, stride=stride,\n",
    "          cache=cache, cutoffs=cutoffs, verbose=verbose)\n",
    "end\n",
    "end\n",
    "nothing\n")
}

# A collapsed refit keeps parameter draws only: there is no trajectory. The
# cache then holds draws alone, which the exact scorer reads directly and the
# simulation scorer completes by backward sampling from the training data.
lfo_fit_marginal_src <- function() paste0(
  "function et_lfo_fit(data; n_sweeps, n_burn=0, n_adapts=0, seed=1,\n",
  "                     adtype=ADTypes.AutoForwardDiff(), x_init=nothing, thin=1)\n",
  "    x_init === nothing || error(\"a marginal-likelihood model takes no x_init\")\n",
  "    et_check_independent!()\n",
  "    m = et_the_model(data, et_marginal_for(data; threads=et_marginal_threads(adtype)))\n",
  "    spl = et_sampler_for(data, nothing)\n",
  "    rng = StableRNG(seed)\n",
  "    t, state = AbstractMCMC.step(rng, m, spl; init=INIT_PARS, adtype=adtype,\n",
  "                                 n_adapts=n_adapts)\n",
  "    for _ in 1:n_burn\n",
  "        t, state = AbstractMCMC.step(rng, m, spl, state; n_adapts=n_adapts)\n",
  "    end\n",
  "    n_keep = max(1, div(n_sweeps, thin))\n",
  "    draws = Vector{Any}(undef, n_keep)\n",
  "    k = 0\n",
  "    for sweep in 1:n_sweeps\n",
  "        sweep > 1 && ((t, state) = AbstractMCMC.step(rng, m, spl, state; n_adapts=n_adapts))\n",
  "        sweep % thin == 0 && k < n_keep || continue\n",
  "        k += 1\n",
  "        draws[k] = et_params(t)\n",
  "    end\n",
  "    draws\n",
  "end\n")

lfo_fit_augmented_src <- function(has_obs) {
  paste0(
    "function et_lfo_fit(data; n_sweeps, n_burn=0, n_adapts=0, seed=1,\n",
    "                     adtype=ADTypes.AutoForwardDiff(), x_init=nothing, thin=1)\n",
    "    # All-susceptible is a poor start for a model whose infection route is\n",
    "    # transmission. iFFBS resamples one individual at a time, so with nobody\n",
    "    # infected the only way in is the background hazard, and a small one\n",
    "    # leaves the chain effectively stuck. Worse, the transmission term is\n",
    "    # multiplied by a zero infected count throughout, so its rate is\n",
    "    # unidentified and two models differing only in that term fit\n",
    "    # identically. Pass `x_init` when a better configuration is available;\n",
    "    # truncation keeps the full time dimension, so one matrix is the right\n",
    "    # shape at every cutoff.\n",
    "    X0 = x_init === nothing ? fill(1, data.n_timepoints, data.n_individuals) :\n",
    "         Matrix{Int}(x_init)\n",
    "    reset_aggregates!(data)\n",
    "    apply_derived_summaries!(et_params(INIT_PARS), data, X0)\n",
    "    # Every artefact is rebuilt from `data`, not taken from the module\n",
    "    # constants. SPL's latent block closes over the data it was built from, so\n",
    "    # stepping SPL here would resample all T occasions against the full\n",
    "    # observations while the model saw only the truncated copy: a fit at\n",
    "    # cutoff t that has seen the future.\n",
    "    m = et_the_model(data, et_loglik_for(data)",
    if (has_obs) ", et_obsloglik_for(data)" else "", ")\n",
    "    spl = et_sampler_for(data, et_latent_for(data))\n",
    "    init = (; X=X0, INIT_PARS...)\n",
    "    # Stepped BY HAND, exactly as et_collect_run() is: `X` must be RETAINED\n",
    "    # here (LFO needs every draw's whole trajectory to forward-simulate from),\n",
    "    # which is the opposite of a normal fit's `save_states=(X=:buffer,)`. The\n",
    "    # bundled chain drops per-sweep `X` unconditionally, so the raw\n",
    "    # transitions from `step` are read directly instead.\n",
    "    rng = StableRNG(seed)\n",
    "    t, state = AbstractMCMC.step(rng, m, spl; init=init, adtype=adtype,\n",
    "                                 n_adapts=n_adapts)\n",
    "    for _ in 1:n_burn\n",
    "        t, state = AbstractMCMC.step(rng, m, spl, state; n_adapts=n_adapts)\n",
    "    end\n",
    "    # Every `thin`-th sweep is kept; each kept draw holds a whole trajectory.\n",
    "    n_keep = max(1, div(n_sweeps, thin))\n",
    "    draws = Vector{Any}(undef, n_keep)\n",
    "    Xs = Vector{Matrix{Int}}(undef, n_keep)\n",
    "    k = 0\n",
    "    for sweep in 1:n_sweeps\n",
    "        sweep > 1 && ((t, state) = AbstractMCMC.step(rng, m, spl, state; n_adapts=n_adapts))\n",
    "        sweep % thin == 0 && k < n_keep || continue\n",
    "        k += 1\n",
    "        draws[k] = et_params(t); Xs[k] = copy(t.X)\n",
    "    end\n",
    "    (draws, Xs)\n",
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
  JuliaCall::julia_assign(src_sym, qualify_base(body))
  raw <- JuliaCall::julia_eval(sprintf("Base.include_string(%s, Main.%s)",
                                       res$module, src_sym))
  as.data.frame(raw, stringsAsFactors = FALSE)
}

#' Per-cell scores from a leave-future-out sweep.
#'
#' The terms each window's score is the sum of. Under `"by_individual"` a cell
#' is one individual's whole forecast block, so this is one row per
#' (window, individual): what stacking over individual histories needs, and
#' where to look when a window scores `-Inf`. A window total cannot be split
#' back into these, because the log is taken per cell.
#'
#' @param res An [et_lfo_cv()] result.
#' @param granularity Which granularity; defaults to the first one run.
#' @return A data frame: `cutoff`, `cell` (the individual under
#'   `"by_individual"`, `1` under `"joint"`, `"i:m"` or `"g:m"` for the per-step
#'   granularities), `elpd`, and `n_finite`, the number of posterior draws that
#'   gave the cell a finite density.
#' @export
et_lfo_cells <- function(res, granularity = res$granularities[1]) {
  if (!inherits(res, "et_lfo_result")) {
    stop("et_lfo_cells(): `res` must come from et_lfo_cv().", call. = FALSE)
  }
  granularity <- match.arg(granularity, res$granularities)
  body <- sprintf(paste0(
    "let c = cell_elpd(%s, %s)\n",
    "    Dict(\"cutoff\" => c.cutoff,\n",
    "         \"cell\" => [k isa Tuple ? join(k, \":\") : string(k) for k in c.cell],\n",
    "         \"elpd\" => c.elpd, \"n_finite\" => c.n_finite)\n",
    "end"), res$sym, julia_symbol(granularity))
  src_sym <- paste0(res$module, "_lfo_cells_src")
  JuliaCall::julia_assign(src_sym, qualify_base(body))
  raw <- JuliaCall::julia_eval(sprintf("Base.include_string(%s, Main.%s)",
                                       res$module, src_sym))
  cell <- as.character(raw$cell)
  out <- data.frame(cutoff = as.integer(raw$cutoff),
                    cell = if (all(grepl("^[0-9]+$", cell))) as.integer(cell) else cell,
                    elpd = as.numeric(raw$elpd),
                    n_finite = as.integer(raw$n_finite),
                    stringsAsFactors = FALSE)
  out[order(out$cutoff, out$cell), , drop = FALSE]
}

#' A leave-future-out result as plain data.
#'
#' An [et_lfo_cv()] result is a handle into the Julia session that made it, so
#' it cannot be saved and reloaded. This collects everything a later gather
#' needs into ordinary R objects, ready for `saveRDS()`.
#'
#' @param res An [et_lfo_cv()] result.
#' @return A list: `windows` ([et_lfo_windows()]), `cells` (a named list of
#'   [et_lfo_cells()], one per granularity), and `method`, `likelihood`, `L`,
#'   `M`.
#' @export
et_lfo_tables <- function(res) {
  if (!inherits(res, "et_lfo_result")) {
    stop("et_lfo_tables(): `res` must come from et_lfo_cv().", call. = FALSE)
  }
  list(windows = et_lfo_windows(res),
       cells = stats::setNames(lapply(res$granularities,
                                      function(g) et_lfo_cells(res, g)),
                               res$granularities),
       method = res$method, likelihood = res$likelihood, L = res$L, M = res$M)
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

#' Model weights from a set of leave-future-out sweeps.
#'
#' A margin says which model won and by how much; a weight says how much of the
#' predictive job each model should be given, which is usually what a reader of
#' a model comparison wants to know.
#'
#' Two methods, answering different questions. `"stacking"` (the default)
#' chooses the weights that maximise the predictive density of the weighted
#' MIXTURE, window by window, so a model that predicts well where the others
#' predict badly earns weight even if its total is not the best.
#' `"pseudo_bma"` is a softmax of the totals; it ignores between-window
#' structure and collapses onto the single best model as the series lengthens,
#' and is offered because it is what most people mean by "model weight".
#'
#' Windows in which ANY candidate scored non-finite are dropped, and the count
#' is returned as `n_dropped`. Keeping them would let one model's degeneracy
#' set the weights -- which is the mortality failure the survival-constrained
#' proposal exists to fix, so it is reported rather than hidden.
#'
#' @param results A named list of [et_lfo_cv()] results, one per candidate. The
#'   names label the models in the output.
#' @param granularity Which granularity to weight on; defaults to the first one
#'   the results were run under. Mixing granularities is refused.
#' @param method `"stacking"` or `"pseudo_bma"`.
#' @return A data frame with `model` and `weight`, carrying `method`,
#'   `granularity`, `n_windows` and `n_dropped` as attributes.
#' @export
et_lfo_weights <- function(results, granularity = NULL, method = "stacking") {
  if (!is.list(results) || length(results) < 2L) {
    stop("et_lfo_weights(): `results` must be a list of at least two ",
         "et_lfo_cv() results.", call. = FALSE)
  }
  if (is.null(names(results)) || any(!nzchar(names(results)))) {
    stop("et_lfo_weights(): name the list, so the weights can be read back ",
         "against the models -- list(null = a, rate = b, suscept = c).",
         call. = FALSE)
  }
  for (r in results) {
    if (!inherits(r, "et_lfo_result")) {
      stop("et_lfo_weights(): every element must come from et_lfo_cv().",
           call. = FALSE)
    }
  }
  mods <- vapply(results, function(r) r$module, character(1))
  if (length(unique(mods)) != 1L) {
    stop("et_lfo_weights(): every result must come from models generated in ",
         "the SAME Julia module.", call. = FALSE)
  }
  method <- match.arg(method, c("stacking", "pseudo_bma"))
  if (is.null(granularity)) granularity <- results[[1]]$granularities[1]
  granularity <- match.arg(granularity, results[[1]]$granularities)

  # Symbol keys, not string keys: the call is built inside a Julia string that
  # is itself inside an R string, so a `"name"` here would close the outer
  # string early and emit invalid Julia. Symbols need no quotes at all.
  check_julia_name(names(results), "model name")
  pairs <- paste(sprintf("%s => %s",
                         vapply(names(results), julia_symbol, character(1)),
                         vapply(results, function(r) r$sym, character(1))),
                 collapse = ", ")
  out <- JuliaCall::julia_eval(sprintf(
    "Base.include_string(%s, \"model_weights(Dict(%s); granularity=%s, method=%s)\")",
    results[[1]]$module, pairs, julia_symbol(granularity),
    julia_symbol(method)))

  d <- data.frame(model = as.character(out$names),
                  weight = as.numeric(out$weights),
                  stringsAsFactors = FALSE)
  d <- d[order(-d$weight), , drop = FALSE]
  rownames(d) <- NULL
  attr(d, "method") <- method
  attr(d, "granularity") <- granularity
  attr(d, "n_windows") <- out$n_windows
  attr(d, "n_dropped") <- out$n_dropped
  d
}

#' Model weights from per-window scores.
#'
#' The same weighting as [et_lfo_weights()], but from a matrix of per-window
#' scores rather than from live [et_lfo_cv()] results.
#'
#' This is what a two-stage pipeline needs. An [et_lfo_cv()] result holds a
#' handle into a generated Julia module, so it does not survive the session
#' that made it -- a sweep that fits on a cluster and gathers afterwards has
#' the numbers but not the objects. Per-window scores, being a data frame,
#' do survive, and are enough to weight from.
#'
#' Windows with a non-finite score for ANY model are dropped and counted, for
#' the reason [et_lfo_weights()] documents.
#'
#' @param scores A numeric matrix, one row per window and one column per
#'   model, of per-window scores under a single granularity. Column names, if
#'   present, label the models.
#' @param method `"stacking"` or `"pseudo_bma"`.
#' @return A data frame with `model` and `weight`, carrying `method`,
#'   `n_windows` and `n_dropped` as attributes.
#' @export
et_weights_from_scores <- function(scores, method = "stacking") {
  et_require_session()
  scores <- as.matrix(scores)
  if (!is.numeric(scores) || ncol(scores) < 2L) {
    stop("et_weights_from_scores(): `scores` must be a numeric matrix with ",
         "one column per model and at least two models.", call. = FALSE)
  }
  nms <- colnames(scores)
  if (is.null(nms)) nms <- paste0("model", seq_len(ncol(scores)))
  check_julia_name(nms, "model name")
  method <- match.arg(method, c("stacking", "pseudo_bma"))

  JuliaCall::julia_assign("_et_w_scores", scores)
  JuliaCall::julia_command(
    "using EpidemicTrajectories: model_weights")
  out <- JuliaCall::julia_eval(sprintf(
    "model_weights(_et_w_scores; names=%s, method=%s)",
    julia_symbol_vector(nms), julia_symbol(method)))

  d <- data.frame(model = nms, weight = as.numeric(out$weights),
                  stringsAsFactors = FALSE)
  d <- d[order(-d$weight), , drop = FALSE]
  rownames(d) <- NULL
  attr(d, "method") <- method
  attr(d, "n_windows") <- out$n_windows
  attr(d, "n_dropped") <- out$n_dropped
  d
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
  JuliaCall::julia_command("import Serialization")
  out <- do.call(rbind, lapply(files, function(f) {
    src <- sprintf(lfo_diagnostics_src(),
                   julia_string(julia_path(normalizePath(f, mustWork = TRUE))))
    d <- as.data.frame(JuliaCall::julia_eval(src), stringsAsFactors = FALSE)
    if (!nrow(d)) return(NULL)
    cbind(cutoff = as.integer(sub("^.*fit_t0*([0-9]+)[.]jls$", "\\1", f)), d,
          stringsAsFactors = FALSE)
  }))
  out[order(out$ess), ]
}

#' The parameter draws a leave-future-out sweep cached.
#'
#' One row per (cutoff, draw), one column per scalar parameter, read from the
#' fits [et_lfo_cv()] wrote to `cache`. Plain data, so draws from separate
#' chains (separate caches) can be compared or pooled after the session that
#' made them has gone.
#'
#' @param cache The cache directory given to [et_lfo_cv()].
#' @return A data frame: `cutoff`, `draw`, then the parameters.
#' @export
et_lfo_draws <- function(cache) {
  et_require_session()
  files <- sort(list.files(cache, pattern = "^fit_t[0-9]+[.]jls$",
                           full.names = TRUE))
  if (!length(files)) {
    stop("et_lfo_draws(): no fit_t*.jls in ", cache, call. = FALSE)
  }
  JuliaCall::julia_command("import Serialization")
  do.call(rbind, lapply(files, function(f) {
    raw <- JuliaCall::julia_eval(sprintf(paste0(
      "let fit = open(Serialization.deserialize, %s)\n",
      "    draws = fit isa Tuple ? fit[1] : fit\n",
      "    nms = [k for (k, v) in pairs(draws[1]) if v isa Real]\n",
      "    Dict(String(n) => Float64[getproperty(d, n) for d in draws] for n in nms)\n",
      "end"), julia_string(julia_path(normalizePath(f, mustWork = TRUE)))))
    d <- as.data.frame(raw)
    cbind(cutoff = as.integer(sub("^.*fit_t0*([0-9]+)[.]jls$", "\\1", f)),
          draw = seq_len(nrow(d)), d)
  }))
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
  # An augmented fit is cached as (draws, trajectories), a collapsed one as
  # the draws alone; destructuring the latter would take its first two draws.
  "    fit = open(Serialization.deserialize, %s)\n",
  "    draws = fit isa Tuple ? fit[1] : fit\n",
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
