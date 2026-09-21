# The badger bTB model in R, tuned.
#
# The SAME model as badger_friendly.R -- same states, same rates, same
# observation process, same posterior -- with every performance seam the package
# offers turned on, plus post-hoc residuals. Read the friendly version first;
# this is a diff against it.
#
# What changes, and why each is exact rather than an approximation:
#
#   1. coupled_transitions = list(c("S", "E"))
#      The only way one badger changes another's fate is by contributing to its
#      force of infection. A neighbour whose realised move is not S -> E has the
#      same probability whatever the focal does, so it cancels on normalisation
#      and can be skipped.
#
#   2. likelihood_weight -- the scalar observation path
#      The likelihood needs exactly one entry of the weight vector, w[X[t,i]].
#      The vector form allocates a whole array per (individual, time) and throws
#      all but one element away. It also lets the LIKELIHOOD score a different
#      factor from the FILTER, which is what keeps `etas` conjugate.
#
#   3. et_summary() -- the O(n_states) coupling
#      The default coupling recomputes, per candidate focal state, how probable
#      every affected neighbour's realised move is. Keeping per-(group, time)
#      running totals of the S -> E and S -> S moves makes it a closed form in
#      those totals, with no neighbour loop. Those totals depend on where the
#      individual goes NEXT, which et_aggregate() cannot express -- hence
#      et_summary(), where you write the reverse yourself.
#
#   4. the depends= split -- DERIVED, not declared
#      Each Gibbs block skips likelihood terms its parameters cannot move. You do
#      not write the annotation: the package reads which model$x each of your
#      functions touches and emits it. et_check_depends() verifies it.
#
#   5. residuals -- survival, exposure and latent period
#      The trajectory is streamed to disc during the fit and the residuals are
#      computed from it afterwards, so they can be recomputed, or new ones added,
#      without refitting.
#
# Run:  Rscript inst/examples/showcase/badger_power.R
#       BADGER_SWEEPS=200 Rscript ...        (default is 20)

library(EpidemicTrajectoriesR)

# Locate the sibling loader whether this is Rscript-ed, source()d or pasted.
.here <- local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f) && file.exists(file.path(dirname(f[1]), "badger_data.R"))) {
    return(dirname(normalizePath(f[1])))
  }
  # sys.source() / source() set `ofile` on the calling frame instead.
  for (k in rev(sys.parents())) {
    of <- tryCatch(get("ofile", envir = sys.frame(k)), error = function(e) NULL)
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  for (d in c(getwd(), "inst/examples/showcase")) {
    if (file.exists(file.path(d, "badger_data.R"))) return(normalizePath(d))
  }
  stop("cannot locate badger_data.R; run from the showcase directory")
})
source(file.path(.here, "badger_data.R"))

states <- c("S", "E", "I", "D")

## ---------------------------------------------------------------------------
## Aggregates -- the two counts, plus two running totals for the coupling
## ---------------------------------------------------------------------------
# n_infectious and n_alive are declared exactly as in the friendly version: the
# package derives both directions of their update from the `A[i] <- A[i] + x`
# shape.
#
# nSE and nSS cannot be. Their contribution depends on where the individual goes
# NEXT (X[t + 1, i]), which that shape does not expose, so they are written with
# et_summary() -- you supply the reverse yourself.
#
# The contract is the same either way: applying the update and then reversing it
# must leave the array exactly as it was. That reversibility is what lets the
# sampler take an individual out, refilter, and put it back, and it is the whole
# basis of iFFBS being both correct and cheap.
#
# Note `X[t + 1, i] == "E"`: the trajectory holds state CODES, and the package
# resolves the name to its index. Comparing against a name your state space does
# not contain is an error rather than a comparison that never matches.

aggs <- et_aggregates(
  states,
  arrays = list(
    n_infectious = et_array("Int", c(raw$n_groups, raw$n_timepoints)),
    n_alive      = et_array("Int", c(raw$n_groups, raw$n_timepoints)),

    nSE = et_summary("Int", c(raw$n_groups, raw$n_timepoints),
      update = function(model, data, X, state, i, t, reverse) {
        g <- data$social_group[i, t]
        if (g > 0 && t < data$n_timepoints) {
          contrib <- (state == "S") && (X[t + 1, i] == "E")
          if (reverse) {
            nSE[g, t] <- nSE[g, t] - contrib
          } else {
            nSE[g, t] <- nSE[g, t] + contrib
          }
        }
      }),

    nSS = et_summary("Int", c(raw$n_groups, raw$n_timepoints),
      update = function(model, data, X, state, i, t, reverse) {
        g <- data$social_group[i, t]
        if (g > 0 && t < data$n_timepoints) {
          contrib <- (state == "S") && (X[t + 1, i] == "S")
          if (reverse) {
            nSS[g, t] <- nSS[g, t] - contrib
          } else {
            nSS[g, t] <- nSS[g, t] + contrib
          }
        }
      })),

  update = function(model, data, X, state, i, t) {
    if (data$social_group[i, t] > 0) {
      n_infectious[data$social_group[i, t], t] <-
        n_infectious[data$social_group[i, t], t] + (state == "I")
    }
    if (data$social_group[i, t] > 0) {
      n_alive[data$social_group[i, t], t] <-
        n_alive[data$social_group[i, t], t] + (state != "D")
    }
  })

## ---------------------------------------------------------------------------
## Rates -- identical to the friendly version
## ---------------------------------------------------------------------------

infection <- function(model, data, i, t) {
  g <- data$social_group[i, t]
  if (g == 0) return(0)
  infectious <- data$aggregates$n_infectious[g, t]
  alive      <- data$aggregates$n_alive[g, t]
  if (alive == 0) return(0)
  foi <- model$lambda * model$alpha[g] +
         model$beta * infectious / ((alive / data$K)^model$q)
  -expm1(-foi)
}

progression <- function(model, data, i, t) {
  x <- data$k / model$tau
  s <- 1
  term <- 1
  for (j in 1:(data$k - 1)) {
    term <- term * x / j
    s <- s + term
  }
  1 - s * exp(-x)
}

survival_fn <- function(model, data, i, t) {
  tt <- t + 1
  age <- if (tt <= data$n_timepoints) data$age[i, tt]
         else data$age[i, data$n_timepoints] + (tt - data$n_timepoints)
  if (age < 0) return(1)
  y1 <- model$b2 * (age - 1)
  y2 <- model$b2 * age
  exp(-model$c1 + (model$a2 / model$b2) * (-exp(y1) * expm1(y2 - y1)))
}

trans <- et_transitions(states,
  survival = et_survival(survival_fn, death = "D"),
  "S -> E" = infection,
  "E -> I" = progression)

starting_state <- function(model, data, X, i, t) {
  p <- numeric(data$n_states)
  st <- data$sampling_period[i][1]
  nuE <- 0
  nuI <- 0
  if (data$birth_time[i] < st) {
    idx <- 0
    for (j in 1:length(data$nu_times)) {
      if (data$nu_times[j] == st) idx <- j
    }
    if (idx > 0) {
      nuE <- model$nu[idx, 1]
      nuI <- model$nu[idx, 2]
    }
  }
  p[1] <- 1 - nuE - nuI
  p[2] <- nuE
  p[3] <- nuI
  p[4] <- 0
  p
}

## ---------------------------------------------------------------------------
## 3. The O(n_states) coupling
## ---------------------------------------------------------------------------
# This is the one place the R DSL hands over to Julia, and deliberately so. The
# closed form assumes a single coupled transition (S -> E) with a per-group
# frequency-dependent force of infection -- a modelling choice the generic core
# does not bake in, and one that cannot be inferred from the spec.
#
# et_julia() takes literal Julia anywhere a transpiled R function is accepted.
# It reads the same aggregates the R code declared above.

coupling <- et_julia('
    t == data.n_timepoints && return ones(n_states)
    g = data.social_group[i, t]
    g == 0 && return ones(n_states)
    I_minus = data.aggregates.n_infectious[g, t]
    M_minus = data.aggregates.n_alive[g, t]
    nSE, nSS = data.aggregates.nSE[g, t], data.aggregates.nSS[g, t]
    logw = zeros(Float64, n_states)
    @inbounds for s in 1:n_states
        # The focal was reversed out of the aggregates for its own resample, so
        # its own contribution to the denominator is added back per candidate.
        I, M = s == 3 ? (I_minus + 1, M_minus + 1) :
               s == 4 ? (I_minus, M_minus) : (I_minus, M_minus + 1)
        foi = M == 0 ? 0.0 :
              -expm1(-(model.lambda * model.alpha[g] +
                       model.beta * I / ((M / data.K)^model.q)))
        logw[s] = nSE * log(max(foi, 1e-12)) + nSS * log(max(1.0 - foi, 1e-12))
    end
    logw .-= maximum(logw)          # stabilise before exponentiating
    return exp.(logw)
')

## ---------------------------------------------------------------------------
## 2. The observation process, in both forms
## ---------------------------------------------------------------------------
# The FILTER needs the whole weight vector: it must know a non-capture is
# informative, and that a badger seen alive later cannot be dead now.
#
# The LIKELIHOOD needs one entry -- and a DIFFERENT factor. `etas` is drawn by an
# exact conjugate kernel, so scoring capture in the likelihood too would
# double-count it. Because the weights multiply, dropping a factor from the
# likelihood drops exactly its term.
#
# This is also what keeps `etas` out of the derived depends=: it reaches the
# filter, which is never differentiated, and never reaches the likelihood.

observation <- function(model, data, X, i, t, s) {
  eta <- model$etas[data$season[t]]
  capture_w <- if (data$capture[t, i] == 0) {
    if (s == 4) ifelse(t <= data$last_capture_time[i], 0, 1) else 1 - eta
  } else {
    ifelse(s == 4, 0, eta)
  }
  capture_w * tests_weight(model, data, X, i, t, s)
}

# Tests only, one state. No array, no allocation per cell.
#
# et_one() rather than 1: `w` accumulates a product of test-parameter terms, so
# it must start at the parameter scalar type. Starting at an integer would make
# the Julia accumulator a Union the moment a parameter multiplied into it --
# runtime dispatch in the hottest loop in the package.
tests_weight <- function(model, data, X, i, t, s) {
  if (data$capture[t, i] == 0) return(1)
  if (s == 4) return(1)
  w <- et_one()
  for (j in 1:dim(data$tests)[3]) {
    nneg <- data$n_neg[t, i, j]
    npos <- data$n_pos[t, i, j]
    if (nneg + npos > 0) {
      if (s == 1) {
        phi <- model$phis[j]
        if (npos > 0) w <- w * (1 - phi)^npos
        if (nneg > 0) w <- w * phi^nneg
      } else if (s == 2) {
        tr <- model$thetas[j] * model$rhos[j]
        if (npos > 0) w <- w * tr^npos
        if (nneg > 0) w <- w * (1 - tr)^nneg
      } else {
        th <- model$thetas[j]
        if (npos > 0) w <- w * th^npos
        if (nneg > 0) w <- w * (1 - th)^nneg
      }
    }
  }
  w
}

## ---------------------------------------------------------------------------
## Assemble
## ---------------------------------------------------------------------------

dat <- et_data(
  n_individuals  = raw$n_individuals,
  n_timepoints   = raw$n_timepoints,
  transitions    = trans,
  starting_state = starting_state,
  aggregates     = aggs,
  observation_weight  = observation,     # filter: capture x tests
  likelihood_weight   = tests_weight,    # likelihood: tests only  (2)
  rest_contribution   = coupling,        # O(n_states) coupling    (3)
  coupled_transitions = list(c("S", "E")),                       # (1)
  sampling_period     = raw$sampling_period,
  affected_individuals = et_affected_from_groups(raw$social_group),
  helpers = list(et_helper(tests_weight, "tests_weight")),
  extras = list(
    social_group = raw$social_group, age = raw$age, capture = raw$capture,
    capt_effort = raw$capt_effort, tests = raw$tests,
    n_neg = raw$n_neg, n_pos = raw$n_pos, season = raw$season,
    birth_time = raw$birth_time, nu_times = raw$nu_times,
    last_capture_time = raw$last_capture_time,
    first_capture_time = raw$first_capture_time,
    K = raw$K, k = raw$k))

mod <- et_model(
  data = dat,
  parameters = list(
    tau    = prior(exponential_dist(1 / 100), init = 5.0),
    alpha  = prior(exponential_dist(1), init = rep(0.5, raw$n_groups),
                    n = raw$n_groups),
    lambda = prior(exponential_dist(1), init = 0.5),
    beta   = prior(exponential_dist(1), init = 0.3),
    q      = prior(beta_dist(1, 1), init = 0.2),
    c1     = prior(exponential_dist(1), init = 0.45),
    a2     = prior(exponential_dist(1), init = 0.05),
    b2     = prior(exponential_dist(1), init = 0.3),
    thetas = prior(beta_dist(1, 1), init = rep(0.3, raw$n_tests), n = raw$n_tests),
    rhos   = prior(beta_dist(1, 1), init = rep(0.5, raw$n_tests), n = raw$n_tests),
    phis   = prior(beta_dist(1, 1), init = rep(0.5, raw$n_tests), n = raw$n_tests),
    etas   = prior(beta_dist(1, 1), init = rep(0.3, raw$n_seasons),
                    n = raw$n_seasons),
    nu     = prior(init = rep(0.05, raw$n_nu_times * 2),
                    dim = c(raw$n_nu_times, 2), kind = "latent")),
  entry_time = "first_capture_time")

## ---------------------------------------------------------------------------
## The sampler
## ---------------------------------------------------------------------------
# `n_steps = 15` matches the C++ reference's EXPECTED trajectory length (it
# draws L uniform on 1..30, mean 15.5), not its nominal L = 30.

blocks <- list(
  # AdaptiveHMC learns its metric and step size during warm-up, then holds both
  # fixed -- a constant `n_steps` leapfrog steps per iteration, no tree
  # doubling, and nothing for you to hand-tune.
  #
  # This replaces a hand-tuned et_hmc() metric. A fixed metric has to be right in
  # both ORDER and SCALE, and a scale far too small gives a chain that never
  # moves while the acceptance diagnostics still look healthy - which is exactly
  # what the hand-tuned badger metric did here.
  et_adaptive_hmc(c("tau", "alpha", "lambda", "beta", "q", "c1", "a2", "b2",
                    "thetas", "rhos", "phis"),
                  target_accept = 0.8, n_steps = 15, metric = "dense"),

  et_conjugate_capture_prob("etas",
    caught = "capture", effort = "capt_effort", group = "social_group",
    index  = "season", dead_state = "D", n = raw$n_seasons),

  et_conjugate_initial_state("nu",
    at = "nu_times", states = c("S", "E", "I"), n = raw$n_nu_times,
    eligible = function(X, data, i, t) {
      data$sampling_period[i][1] == t && data$birth_time[i] < t
    }),

  et_iffbs("X"))

## ---------------------------------------------------------------------------
## 5. Residuals
## ---------------------------------------------------------------------------
# Each is a PIT residual: uniform on (0, 1) if the mechanism it targets is
# right, so a departure from uniformity points at THAT mechanism rather than
# saying "the model is bad" in general.
#
# Nothing here is hard-coded to this model. `death = "D"` is a name resolved
# against your declared state space, and the waiting-time residual derives its
# hazard from the transitions -- any declared move gets one for free.

residuals_spec <- list(
  # How long a badger LIVES, against Gompertz-Makeham. Left-truncated because we
  # only ever see badgers that survived to first capture: without conditioning
  # on that, the residual is calibrated against the unconditional lifetime and
  # reports the sampling design as model misfit.
  et_residual_survival(
    origin       = "birth_time",
    condition_on = "first_capture_time",
    censor_at    = "last_capture_time",
    death        = "D",
    name         = "survival"),

  # How long a SUSCEPTIBLE badger waits before infection. The clock starts when
  # observation starts, not at entry to S -- it was already susceptible when we
  # began watching. Death is a competing risk, so it CENSORS rather than
  # dropping: a badger that died before catching it has not falsified anything,
  # and dropping it would bias the residual towards badgers that lived longer.
  et_residual_waiting("S", "E",
    origin = "window_start", censor_at = c("window_end", "D"),
    accumulate = "discrete_product", name = "exposure"),

  # How long an EXPOSED badger takes to become infectious -- the latent period,
  # the residual for `tau`. Its clock starts at entry to E, read from the
  # trajectory.
  et_residual_waiting("E", "I",
    origin = "entry_to_from_state", censor_at = c("window_end", "D"),
    accumulate = "discrete_product", name = "latent_period"))

if (!isTRUE(getOption("etr.example.dry_run"))) {
  et_setup()

  n_sweeps <- as.integer(Sys.getenv("BADGER_SWEEPS", "20"))

  cat("\n=== structural checks (no sampling) ===\n")
  ll <- et_loglik(mod, blocks, x_init = raw$X_init)
  cat(sprintf("log density at init:  epidemic %.4f   observation %.4f\n",
              ll[["epidemic"]], ll[["observation"]]))
  stopifnot(all(is.finite(ll)))

  # Verify the DERIVED depends= before trusting a single gradient. Perturbs each
  # parameter and checks no term it was excluded from actually moves.
  dep <- et_check_depends(mod, blocks, x_init = raw$X_init)
  bad <- dep[dep$verdict == "UNDER-DECLARED", , drop = FALSE]
  if (nrow(bad)) {
    print(bad)
    stop("depends= is under-declared: those gradients would be silently wrong.")
  }

  cat("\n=== sampling ===\n")
  x_path <- file.path(.here, "badger_power_X.jld2")
  fit <- et_sample(mod, blocks, n_sweeps = n_sweeps, seed = 13,
                   x_init = raw$X_init,
                   # Stream the trajectory to disc: it is 161 x 2384 integers
                   # EVERY sweep, so it is kept live for conditioning but out of
                   # the chain object. Thinned, because at 5000 draws the
                   # unthinned archive is ~1.9 GB.
                   save_x = x_path, save_every = 100)
  print(summary(fit))

  cat("\n=== residuals (PIT; Uniform(0,1) if the model is right) ===\n")
  # sync_aggregates because the force of infection reads data$aggregates: a
  # spec-derived hazard needs them to agree with THIS draw's trajectory, not
  # with whatever was current when the fit ended.
  res <- et_residuals(fit, residuals_spec, sync_aggregates = TRUE)
  print(et_residual_summary(res))
}
