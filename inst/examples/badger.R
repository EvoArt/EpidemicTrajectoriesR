# The badger bTB SEID model, written entirely in R.
#
# 2384 badgers x 161 quarters, four states S -> E -> I with death (D) reachable
# from every live state, six imperfect diagnostic tests, capture-recapture
# observation, time-varying social groups, and Gompertz-Makeham survival.
#
# WHICH VARIANT THIS IS, AND WHY IT MATTERS
# -----------------------------------------
# This is the **Gompertz-Makeham** model with the **repeat-test fix**. That pair
# is what finally reproduced the reference's mean infectious period, recorded in
#   EpidemicTrajectories/latest_badger_ref/HANDOVER_2026-08-02_xi_fix_and_siler.md
# section 4b:
#
#     tau = 14.983  (sd 1.764, ESS 402.8, MCSE 0.088, Rhat 1.003)
#       upstream C++          15.06   (-0.88 MCSE)  indistinguishable
#       semi-Markov reference 14.87   (+1.28 MCSE)  indistinguishable
#
# Two earlier variants were wrong, and both fixes are carried here:
#
#   * BEFORE the repeat-test fix, tau = 16.12. `TestMat.csv` carries MULTIPLE
#     rows for the same (individual, quarter) at 1201 of 12430 cells; 150 of them
#     hold CONTRADICTORY readings. Flattening them to one slot (last row wins)
#     dropped every repeat reading. Both reference implementations score every
#     row. `n_neg` / `n_pos` below restore that.
#   * BEFORE the xi fix, tau = 14.10. The Brock changepoint sampler had gone
#     inert because its multiplicity counts were never moved with it.
#
# NOT Siler. The C++ reference has no early-life (a1, b1) term at all -- `a1` and
# `b1` appear in ZERO .cpp files -- and fitting one lets it absorb signal
# belonging to c1/a2/b2.
#
# LIMITATION, stated plainly: the Brock changepoint `xi` is sampled by a bespoke
# Metropolis kernel (an in-place swap of the test columns plus a windowed MH
# ratio) that the R DSL cannot yet express. This example FIXES it at 101, exactly
# as `EpidemicTrajectories/examples/badger_data.jl` does. The posterior for xi
# visits only 100-102 (98.9% at 100-101), so this is close to the mode -- but it
# is a fixed-xi run, not the sampled-xi run that produced the number above.

library(EpidemicTrajectoriesR)

BADGER_DATA_DIR <- Sys.getenv("BADGER_DATA_DIR", unset = file.path(
  path.expand("~"), ".julia", "dev", "EpidemicTrajectories", "badger_ref", "RData2"))

BROCK_CHANGEPOINT_RAW <- 80   # what the raw TestMat columns encode
BROCK_CHANGEPOINT     <- 101  # what we fix it to

states <- c("S", "E", "I", "D")

## ---------------------------------------------------------------------------
## Data
## ---------------------------------------------------------------------------

# All of this is user code: reading the CSVs and converting the reference's
# conventions into ET's. Two differ and are converted here:
#   * the reference is individual-major (X[i, t]); ET is time-major (X[t, i])
#   * the reference encodes states sparsely (S=0, E=3, I=1, D=9); ET uses the
#     position in the state space, so S=1, E=2, I=3, D=4
load_badger_data <- function(dir = BADGER_DATA_DIR,
                             brock_changepoint = BROCK_CHANGEPOINT,
                             known_sex_only = TRUE) {
  rd <- function(f) as.matrix(utils::read.csv(file.path(dir, f)))

  dims <- utils::read.csv(file.path(dir, "dimensions.csv"))
  m <- dims$m[1]; maxt <- dims$maxt[1]
  n_groups <- dims$G[1]; n_tests <- dims$numTests[1]
  n_seasons <- dims$numSeasons[1]; n_nu_times <- dims$numNuTimes[1]

  Xinit_raw   <- rd("Xinit.csv")
  test_mat    <- rd("TestMat.csv")
  capt_hist   <- rd("CaptHist.csv")
  capt_effort <- rd("CaptEffort.csv")
  birth_time  <- as.vector(rd("birthTimes.csv"))
  start_p     <- as.vector(rd("startSamplingPeriod.csv"))
  end_p       <- as.vector(rd("endSamplingPeriod.csv"))
  nu_times    <- as.vector(rd("nuTimes.csv"))
  sex_raw     <- as.vector(rd("sex.csv"))
  K <- utils::read.csv(file.path(dir, "Kay.csv"))$K[1]
  k <- utils::read.csv(file.path(dir, "k.csv"))$k[1]

  # The reference filters to badgers of known sex; the base model has no sex
  # effects, but keeping the filter makes the dataset comparable with the sex
  # models.
  keep <- if (known_sex_only) which(sex_raw != 0) else seq_len(m)
  old_to_new <- integer(m)
  old_to_new[keep] <- seq_along(keep)

  Xinit_raw  <- Xinit_raw[keep, , drop = FALSE]
  capt_hist  <- capt_hist[keep, , drop = FALSE]
  birth_time <- birth_time[keep]; start_p <- start_p[keep]; end_p <- end_p[keep]

  kept_rows <- old_to_new[test_mat[, 2]] != 0
  test_mat  <- test_mat[kept_rows, , drop = FALSE]
  test_mat[, 2] <- old_to_new[test_mat[, 2]]
  m <- length(keep)

  ref_code <- c("0" = 1L, "3" = 2L, "1" = 3L, "9" = 4L)
  X_init <- matrix(1L, maxt, m)
  for (i in seq_len(m)) {
    codes <- Xinit_raw[i, ]
    X_init[, i] <- ifelse(codes == -10, 1L, ref_code[as.character(codes)])
  }

  # Group membership genuinely varies with t (badgers move); 0 means "absent".
  social_group <- matrix(0L, m, maxt)
  for (i in seq_len(m)) {
    rows <- which(test_mat[, 2] == i)
    if (!length(rows)) next
    times_i  <- test_mat[rows, 1]
    groups_i <- test_mat[rows, 3]
    g <- as.integer(groups_i[which.min(times_i)])
    for (t in max(1, birth_time[i]):maxt) {
      hit <- match(t, times_i)
      if (!is.na(hit)) g <- as.integer(groups_i[hit])
      social_group[i, t] <- g
    }
  }

  age <- matrix(-10L, m, maxt)
  for (i in seq_len(m)) {
    tt <- max(1, birth_time[i]):maxt
    age[i, tt] <- as.integer(tt - birth_time[i])
  }

  # THE REPEAT-TEST FIX. `n_neg` / `n_pos` record how many ROWS read 0 / 1 for
  # each (t, i, j), so a cell with repeat captures in one quarter contributes
  # each reading once -- as both reference implementations do. `tests` keeps the
  # last row's value only for the presence check.
  tests <- array(-1L, dim = c(maxt, m, n_tests))
  n_neg <- array(0L,  dim = c(maxt, m, n_tests))
  n_pos <- array(0L,  dim = c(maxt, m, n_tests))
  for (r in seq_len(nrow(test_mat))) {
    t <- test_mat[r, 1]; i <- test_mat[r, 2]
    if (t < 1 || t > maxt || i < 1 || i > m) next
    for (j in seq_len(n_tests)) {
      v <- test_mat[r, 3 + j]
      if (is.na(v) || v == -10) next
      tests[t, i, j] <- as.integer(v)
      if (v == 1) n_pos[t, i, j] <- n_pos[t, i, j] + 1L
      else        n_neg[t, i, j] <- n_neg[t, i, j] + 1L
    }
  }
  # The Brock changepoint. The raw columns are already correct at xi = 80, and
  # the reference only ever SWAPS the two Brock columns for tests lying between
  # the current and proposed changepoint -- so fixing it at 101 means applying
  # that same swap once, for the tests in [80, 101).
  if (brock_changepoint != BROCK_CHANGEPOINT_RAW) {
    win <- if (brock_changepoint > BROCK_CHANGEPOINT_RAW)
      BROCK_CHANGEPOINT_RAW:(brock_changepoint - 1)
    else brock_changepoint:(BROCK_CHANGEPOINT_RAW - 1)
    win <- win[win >= 1 & win <= maxt]
    swap <- function(A) { tmp <- A[win, , 1]; A[win, , 1] <- A[win, , 2]
                          A[win, , 2] <- tmp; A }
    tests <- swap(tests); n_neg <- swap(n_neg); n_pos <- swap(n_pos)
  }

  capture <- t(capt_hist)
  storage.mode(capture) <- "integer"

  last_capture_time <- vapply(seq_len(m), function(i) {
    w <- which(capture[, i] == 1L); if (length(w)) max(w) else 0L
  }, integer(1))
  first_capture_time <- vapply(seq_len(m), function(i) {
    w <- which(capture[, i] == 1L)
    if (length(w)) min(w) else as.integer(max(birth_time[i], 1))
  }, integer(1))

  cam <- rd("capturesAfterMonit.csv")
  for (r in seq_len(nrow(cam))) {
    oid <- cam[r, 1]
    if (oid < 1 || oid > length(old_to_new) || old_to_new[oid] == 0) next
    nid <- old_to_new[oid]
    last_capture_time[nid] <- max(last_capture_time[nid], cam[r, 2])
  }

  season <- integer(maxt); season[1] <- 1L
  for (t in 2:maxt) season[t] <- if (season[t - 1] < n_seasons) season[t - 1] + 1L else 1L

  list(n_individuals = m, n_timepoints = maxt, n_groups = n_groups,
       n_tests = n_tests, n_seasons = n_seasons, n_nu_times = n_nu_times,
       X_init = X_init, social_group = social_group, age = age,
       capture = capture, capt_effort = matrix(as.integer(capt_effort),
                                               nrow = nrow(capt_effort)),
       tests = tests, n_neg = n_neg, n_pos = n_pos,
       sampling_period = cbind(as.integer(start_p), as.integer(end_p)),
       birth_time = as.integer(birth_time),
       last_capture_time = last_capture_time,
       first_capture_time = first_capture_time,
       season = season, nu_times = as.integer(nu_times),
       K = as.numeric(K), k = as.integer(k))
}

raw <- load_badger_data()
cat(sprintf("Badger data: %d individuals x %d quarters, %d groups, %d tests\n",
            raw$n_individuals, raw$n_timepoints, raw$n_groups, raw$n_tests))

## ---------------------------------------------------------------------------
## The model
## ---------------------------------------------------------------------------

# Per-group alive and infectious counts. During the focal's own filter these are
# leave-one-out, which is exactly what the force of infection needs; ET re-adds
# the focal transiently for the focal's OWN rate, so no `+1` is written here.
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

# S -> E force of infection: a group-level baseline plus a density-dependent
# term, with `q` interpolating between frequency- and density-dependence.
infection <- function(model, data, i, t) {
  g <- data$social_group[i, t]
  if (g == 0) return(0)
  Inf_ <- data$aggregates$n_infectious[g, t]
  M    <- data$aggregates$n_alive[g, t]
  if (M == 0) return(0)
  foi <- model$lambda * model$alpha[g] + model$beta * Inf_ / ((M / data$K)^model$q)
  -expm1(-foi)
}

# E -> I progression: the Erlang(k, tau/k) CDF at one step, so `tau` is the mean
# latent period in quarters and `k` its shape.
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

# GOMPERTZ-MAKEHAM survival for the t -> t+1 step, evaluated at the age at t+1:
#     P(survive | age) = exp(-c1 + (a2/b2) * (exp(b2*(age-1)) - exp(b2*age)))
survival_fn <- function(model, data, i, t) {
  tt <- t + 1
  age <- if (tt <= data$n_timepoints) data$age[i, tt]
         else data$age[i, data$n_timepoints] + (tt - data$n_timepoints)
  if (age < 0) return(1)
  y1 <- model$b2 * (age - 1)
  y2 <- model$b2 * age
  late <- -exp(y1) * expm1(y2 - y1)
  exp(-model$c1 + (model$a2 / model$b2) * late)
}

trans <- et_transitions(states,
  survival = et_survival(survival_fn, death = "D"),
  "S -> E" = infection,
  "E -> I" = progression)

# The initial-state mixing for individuals already alive when watching starts.
# `nu` is owned by a conjugate Dirichlet kernel, so it is a plain Float64 here.
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

# The observation process factorises: capture x tests.
#
#   * the FILTER needs the product (it must know a non-capture is informative,
#     and that a badger seen alive later cannot be dead now);
#   * the LIKELIHOOD scores only the TESTS factor, because `etas` (the capture
#     probabilities) is sampled by an exact conjugate kernel. Scoring capture in
#     both places would double-count it.
#
# et_data() takes both: `observation_weight` for the filter, `likelihood_weight`
# for the likelihood.
capture_weight <- function(model, data, X, i, t, s) {
  eta <- model$etas[data$season[t]]
  if (data$capture[t, i] == 0) {
    if (s == 4) {
      # A badger caught later cannot have been dead now: ban the death state.
      return(ifelse(t <= data$last_capture_time[i], 0, 1))
    }
    return(1 - eta)
  }
  ifelse(s == 4, 0, eta)
}

# Each READING is scored once, via the n_neg / n_pos multiplicities -- the
# repeat-test fix. theta = sensitivity in I, theta*rho = sensitivity in E,
# phi = specificity in S.
tests_weight <- function(model, data, X, i, t, s) {
  if (data$capture[t, i] == 0) return(1)
  if (s == 4) return(1)
  # et_one(), not 1: `w` accumulates a product of test-parameter terms, so it
  # must start at the AD scalar type rather than at an integer.
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

obs_weight <- function(model, data, X, i, t, s) {
  capture_weight(model, data, X, i, t, s) * tests_weight(model, data, X, i, t, s)
}

dat <- et_data(
  n_individuals = raw$n_individuals,
  n_timepoints  = raw$n_timepoints,
  transitions   = trans,
  starting_state = starting_state,
  aggregates    = aggs,
  observation_weight = obs_weight,      # filter: capture x tests
  likelihood_weight  = tests_weight,    # likelihood: tests only (etas conjugate)
  sampling_period    = raw$sampling_period,
  # Badgers move between groups, so who each affects varies with t.
  affected_individuals = et_affected_from_groups(raw$social_group),
  # The only way one badger changes another's fate is by contributing to its
  # force of infection, so only S -> E is coupled. Exact, and worth ~10x.
  coupled_transitions = list(c("S", "E")),
  helpers = list(et_helper(capture_weight, "capture_weight"),
                 et_helper(tests_weight, "tests_weight")),
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
    tau    = et_par(et_exponential(1 / 100), init = 5.0),
    alpha  = et_par(et_exponential(1), init = rep(0.5, raw$n_groups),
                    n = raw$n_groups),
    lambda = et_par(et_exponential(1), init = 0.5),
    beta   = et_par(et_exponential(1), init = 0.3),
    q      = et_par(et_beta(1, 1), init = 0.2),
    c1     = et_par(et_exponential(1), init = 0.45),
    a2     = et_par(et_exponential(1), init = 0.05),
    b2     = et_par(et_exponential(1), init = 0.3),
    thetas = et_par(et_beta(1, 1), init = rep(0.3, raw$n_tests), n = raw$n_tests),
    rhos   = et_par(et_beta(1, 1), init = rep(0.5, raw$n_tests), n = raw$n_tests),
    phis   = et_par(et_beta(1, 1), init = rep(0.5, raw$n_tests), n = raw$n_tests),
    etas   = et_par(et_beta(1, 1), init = rep(0.3, raw$n_seasons),
                    n = raw$n_seasons),
    # Owned by the conjugate Dirichlet kernel: a placeholder density, a real
    # kernel. Declared "latent" so nothing tries to sample it from a prior.
    nu     = et_par(init = rep(0.05, raw$n_nu_times * 2),
                    dim = c(raw$n_nu_times, 2), kind = "latent")),
  # Entry gate: a badger is known alive at first capture, but its disease
  # dynamics before then were not watched. ET scores the disease transitions
  # pre-entry and divides the survival factor back out.
  entry_time = "first_capture_time")

# Per-parameter HMC step sizes, hand-tuned. `n_steps = 15` matches the C++
# reference's EXPECTED trajectory length (it draws L uniform on 1..30, mean
# 15.5), not its nominal L = 30.
hmc_block <- et_hmc(
  c("tau", "alpha", "lambda", "beta", "q", "c1", "a2", "b2",
    "thetas", "rhos", "phis"),
  n_steps = 15,
  step_size = list(tau = 0.002, alpha = 0.2, lambda = 0.01, beta = 0.05,
                   q = 0.05, c1 = 0.02, a2 = 0.001, b2 = 0.001,
                   thetas = 0.005, rhos = 0.005, phis = 0.005))

blocks <- list(
  hmc_block,
  et_conjugate_capture_prob("etas", caught = "capture", effort = "capt_effort",
                            group = "social_group", index = "season",
                            dead_state = "D", n = raw$n_seasons),
  et_conjugate_initial_state("nu", at = "nu_times", states = c("S", "E", "I"),
                             n = raw$n_nu_times,
                             eligible = function(X, data, i, t) {
                               data$sampling_period[i][1] == t &&
                                 data$birth_time[i] < t
                             }),
  et_iffbs("X"))

if (!isTRUE(getOption("etr.example.dry_run"))) {
  et_setup()

  # A converged run is 5000 draws after 1000 adapt -- about 5.5 hours at
  # 3.3 s/sweep. The default here is a short demonstration run.
  n_sweeps <- as.integer(Sys.getenv("BADGER_SWEEPS", "20"))
  n_burn   <- as.integer(Sys.getenv("BADGER_BURN", "0"))

  cat("\n=== Structural checks (no sampling) ===\n")
  ll <- et_loglik(mod, blocks, x_init = raw$X_init)
  cat(sprintf("log density at init:  epidemic %.4f   observation %.4f\n",
              ll[["epidemic"]], ll[["observation"]]))
  stopifnot(all(is.finite(ll)))

  dep <- et_check_depends(mod, blocks, x_init = raw$X_init)
  print(dep[dep$verdict != "ok", , drop = FALSE])

  X1 <- et_iffbs_sweep(mod, raw$X_init, blocks)
  cat(sprintf("one iFFBS sweep changed %d of %d cells\n",
              sum(X1 != raw$X_init), length(X1)))

  cat("\n=== Sampling ===\n")
  fit <- et_sample(mod, blocks = blocks, n_sweeps = n_sweeps, n_burn = n_burn,
                   seed = 13, x_init = raw$X_init)
  print(summary(fit))
  cat("\nReference tau (5000 draws, sampled xi): 14.983; ",
      "upstream C++ 15.06; semi-Markov 14.87\n", sep = "")
}
