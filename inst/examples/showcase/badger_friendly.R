# The badger bTB model in R, written the easy way.
#
# 2384 badgers x 161 quarters. Four states, S -> E -> I, with death reachable
# from any of them. Six imperfect diagnostic tests, capture-recapture
# observation, social groups that badgers move between, and age-dependent
# mortality.
#
# This is the FRIENDLY build: describe the model, let the package do everything
# else. Every performance decision is left at its default. `badger_power.R` is
# the same model with the optimisations turned on; read this one first.
#
# No Julia appears anywhere below. The rate and observation functions are
# ordinary R, transpiled to a single Julia module you can print with
# et_julia_source(mod).
#
# Run:  Rscript inst/examples/showcase/badger_friendly.R
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

## ---------------------------------------------------------------------------
## 1. The states
## ---------------------------------------------------------------------------
# The order fixes the encoding: state k in the trajectory means states[k].

states <- c("S", "E", "I", "D")

## ---------------------------------------------------------------------------
## 2. What to track while the trajectory is resampled
## ---------------------------------------------------------------------------
# The force of infection needs to know how many groupmates are infectious, and
# how many are alive. Declare the arrays and how they update; the package
# allocates them, keeps them in step with the trajectory, and -- because it can
# DERIVE the reverse of each update from its `A[i] <- A[i] + x` shape -- can take
# one individual back out again. That is exactly the leave-one-out count a force
# of infection wants.
#
# The package attaches no meaning to these arrays. They are whatever your rates
# need.

aggs <- et_aggregate(
  states,
  arrays = list(
    n_infectious = et_array("Int", c(raw$n_groups, raw$n_timepoints)),
    n_alive      = et_array("Int", c(raw$n_groups, raw$n_timepoints))),
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
## 3. The rates
## ---------------------------------------------------------------------------
# Ordinary R functions of (model, data, i, t). `model` holds the parameters,
# `data` everything else -- including the aggregates above and any array you
# passed as an extra.

# S -> E. A group-level baseline plus a density-dependent term, with `q`
# interpolating between frequency- and density-dependence.
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

# E -> I. The Erlang(k, tau/k) CDF at one step, so `tau` is the mean latent
# period in quarters and `k` its shape.
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

# Gompertz-Makeham mortality, at the age the badger would reach:
#     P(survive | age) = exp(-c1 + (a2/b2) * (exp(b2*(age-1)) - exp(b2*age)))
survival_fn <- function(model, data, i, t) {
  tt <- t + 1
  age <- if (tt <= data$n_timepoints) data$age[i, tt]
         else data$age[i, data$n_timepoints] + (tt - data$n_timepoints)
  if (age < 0) return(1)
  y1 <- model$b2 * (age - 1)
  y2 <- model$b2 * age
  exp(-model$c1 + (model$a2 / model$b2) * (-exp(y1) * expm1(y2 - y1)))
}

## ---------------------------------------------------------------------------
## 4. The transitions
## ---------------------------------------------------------------------------
# Two moves and a survival declaration. et_survival() scales every non-death
# rate by `survival_fn` and gives EVERY live state a transition to "D" with the
# remaining mass -- including I -> D, which is never written here. The
# self-transitions (S->S, E->E, I->I) take the leftover probability
# automatically.
#
# The arrow strings mirror the Julia line for line.

trans <- et_transitions(states,
  survival = et_survival(survival_fn, death = "D"),
  "S -> E" = infection,
  "E -> I" = progression)

## ---------------------------------------------------------------------------
## 5. The observation process
## ---------------------------------------------------------------------------
# Capture x tests. A badger that was not caught is not tested; one caught later
# cannot have been dead now.
#
# Each READING is scored once, via the repeat-capture multiplicities -- see the
# note in badger_data.R. theta is sensitivity in I, theta*rho in E, phi is
# specificity in S.
#
# `rep(1, ...)` becomes an allocation of the parameter scalar type, so this
# function works whether its parameters are sampled conjugately or sit in an HMC
# block under any autodiff backend. You do not have to think about that.

observation <- function(model, data, X, i, t) {
  w <- rep(1, data$n_states)
  eta <- model$etas[data$season[t]]

  if (data$capture[t, i] == 0) {
    w[1] <- 1 - eta
    w[2] <- 1 - eta
    w[3] <- 1 - eta
    w[4] <- ifelse(t <= data$last_capture_time[i], 0, 1)
    return(w)                      # not caught, so not tested
  }
  w[1] <- eta
  w[2] <- eta
  w[3] <- eta
  w[4] <- 0

  for (j in 1:dim(data$tests)[3]) {
    nneg <- data$n_neg[t, i, j]
    npos <- data$n_pos[t, i, j]
    if (nneg + npos > 0) {
      theta <- model$thetas[j]
      rho   <- model$rhos[j]
      phi   <- model$phis[j]
      tr    <- theta * rho
      if (npos > 0) {
        w[1] <- w[1] * (1 - phi)^npos
        w[2] <- w[2] * tr^npos
        w[3] <- w[3] * theta^npos
      }
      if (nneg > 0) {
        w[1] <- w[1] * phi^nneg
        w[2] <- w[2] * (1 - tr)^nneg
        w[3] <- w[3] * (1 - theta)^nneg
      }
    }
  }
  w
}

# Where a badger starts, for one already alive when watching began.
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
## 6. Assemble
## ---------------------------------------------------------------------------
# Everything the package does not name is yours: pass it in `extras` and read it
# back as data$whatever inside any of the functions above.

dat <- et_data(
  n_individuals  = raw$n_individuals,
  n_timepoints   = raw$n_timepoints,
  transitions    = trans,
  starting_state = starting_state,
  aggregates     = aggs,
  observation_process = observation,
  sampling_period     = raw$sampling_period,
  # Badgers move between groups, so who each one affects varies with time.
  affected_individuals = et_affected_from_groups(raw$social_group),
  extras = list(
    social_group = raw$social_group, age = raw$age, capture = raw$capture,
    capt_effort = raw$capt_effort, tests = raw$tests,
    n_neg = raw$n_neg, n_pos = raw$n_pos, season = raw$season,
    birth_time = raw$birth_time, nu_times = raw$nu_times,
    last_capture_time = raw$last_capture_time,
    first_capture_time = raw$first_capture_time,
    K = raw$K, k = raw$k))

## ---------------------------------------------------------------------------
## 7. Parameters and priors
## ---------------------------------------------------------------------------
# exponential_dist() takes a RATE; Julia's Exponential is parameterised by scale,
# and the package converts. `nu` is owned by a conjugate kernel, so it is
# declared "latent": a placeholder density, a real draw from the kernel.

mod <- et_model(
  data = dat,
  parameters = list(
    tau    = prior(exponential_dist(1 / 100), init = 5.0),
    alpha  = prior(exponential_dist(1), init = rep(0.5, raw$n_groups),
                    n = raw$n_groups),
    lambda = prior(exponential_dist(1), init = 0.5),
    beta   = prior(exponential_dist(1), init = 0.3),
    q      = prior(beta_dist(1, 1), init = 0.2),

    c1 = prior(exponential_dist(1), init = 0.45),   # Makeham constant
    a2 = prior(exponential_dist(1), init = 0.05),   # Gompertz scale
    b2 = prior(exponential_dist(1), init = 0.3),    # Gompertz rate

    thetas = prior(beta_dist(1, 1), init = rep(0.3, raw$n_tests), n = raw$n_tests),
    rhos   = prior(beta_dist(1, 1), init = rep(0.5, raw$n_tests), n = raw$n_tests),
    phis   = prior(beta_dist(1, 1), init = rep(0.5, raw$n_tests), n = raw$n_tests),
    etas   = prior(beta_dist(1, 1), init = rep(0.3, raw$n_seasons),
                    n = raw$n_seasons),
    nu     = prior(init = rep(0.05, raw$n_nu_times * 2),
                    dim = c(raw$n_nu_times, 2), kind = "latent")),
  # Entry gate: a badger is known alive at first capture, but its disease
  # dynamics before then were not watched. The likelihood scores the disease
  # moves pre-entry and divides the survival factor back out.
  entry_time = "first_capture_time")

## ---------------------------------------------------------------------------
## 8. Fit
## ---------------------------------------------------------------------------
# One NUTS block for everything continuous; conjugate kernels where a closed
# form exists (an exact independent draw beats a correlated HMC step); iFFBS for
# the trajectory.
#
# You could omit `blocks` entirely -- every parameter would go into one NUTS
# block and the trajectory would still get iFFBS. It is spelled out here because
# the two conjugate kernels are worth having.

blocks <- list(
  et_nuts(c("tau", "alpha", "lambda", "beta", "q", "c1", "a2", "b2",
            "thetas", "rhos", "phis")),

  et_conjugate_capture_prob("etas",
    caught = "capture", effort = "capt_effort", group = "social_group",
    index  = "season", dead_state = "D", n = raw$n_seasons),

  et_conjugate_initial_state("nu",
    at = "nu_times", states = c("S", "E", "I"), n = raw$n_nu_times,
    eligible = function(X, data, i, t) {
      data$sampling_period[i][1] == t && data$birth_time[i] < t
    }),

  et_iffbs("X"))

if (!isTRUE(getOption("etr.example.dry_run"))) {
  et_setup()

  n_sweeps <- as.integer(Sys.getenv("BADGER_SWEEPS", "20"))

  # A cheap structural check before committing to a fit: a non-finite log
  # density means the model is misspecified before any sampler is involved.
  ll <- et_loglik(mod, blocks, x_init = raw$X_init)
  cat(sprintf("log density at init:  epidemic %.4f   observation %.4f\n",
              ll[["epidemic"]], ll[["observation"]]))
  stopifnot(all(is.finite(ll)))

  fit <- et_sample(mod, blocks, n_sweeps = n_sweeps, seed = 13,
                   x_init = raw$X_init)

  print(summary(fit))
}
