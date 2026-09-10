# The cattle E. coli model (Touloupou et al. 2019), written entirely in R.
#
# Two states, S and I, with recurrent S <-> I: an animal is infected at a rate
# driven by how many of its penmates are infected, and recovers back to
# susceptible. Neither state is ever observed directly -- only two imperfect
# diagnostic tests, on a subset of days.
#
# This is the R counterpart of EpidemicTrajectories' `examples/cattle_user_code.jl`
# and reproduces it line for line. Compare the two side by side: that
# correspondence is the design goal.

library(EpidemicTrajectoriesR)

## ---------------------------------------------------------------------------
## The model
## ---------------------------------------------------------------------------

n_pens        <- 10
n_per_pen     <- 8
n_timepoints  <- 80
n_individuals <- n_pens * n_per_pen
group         <- rep(seq_len(n_pens), each = n_per_pen)
observed_days <- seq(1, n_timepoints, by = 6)

# Declaring the state space fixes the numbering, so the package and I agree on
# which state is which: S is 1, I is 2.
states <- c("S", "I")

# The array I want tracked during the latent update, declared together with how
# it updates. The package allocates it and knows nothing about what it means; it
# just runs this update forwards and, when the sampler needs to, in reverse.
aggs <- et_aggregate(
  states,
  arrays = list(n_infected = et_array("Int", c(n_pens, n_timepoints))),
  update = function(model, data, X, state, i, t) {
    n_infected[data$group[i], t] <- n_infected[data$group[i], t] + (state == "I")
  })

# The rates read that array. `n_infected` excludes the individual being resampled
# while it is being resampled, which is exactly the leave-one-out count the force
# of infection needs.
infection <- function(model, data, i, t) {
  I_minus <- data$aggregates$n_infected[data$group[i], t]
  -expm1(-(model$alpha + model$beta * I_minus))
}

recovery <- function(model, data, i, t) 1 / model$m

trans <- et_transitions(states,
  "S -> I" = infection,
  "I -> S" = recovery)

starting_state <- function(model, data, X, i, t) {
  p <- numeric(data$n_states)
  p[1] <- 1 - model$nu
  p[2] <- model$nu
  p
}

# What I observe, and how it relates to the states: two imperfect tests, each
# with its own sensitivity and perfect specificity. A negative entry means the
# animal was not tested that day, so the result says nothing. The package knows
# none of this -- it just multiplies these weights into the filter.
#
# Note `rep(1, ...)`: the transpiler turns that into an allocation of the AD
# element type, which is what keeps this function working whether its parameters
# are sampled conjugately (Float64) or in an HMC block (ForwardDiff.Dual).
cattle_observations <- function(model, data, X, i, t) {
  w <- rep(1, data$n_states)
  y <- data$rams[t, i]
  if (y >= 0) {
    w[1] <- w[1] * ifelse(y == 1, 0, 1)
    w[2] <- w[2] * ifelse(y == 1, model$theta_r, 1 - model$theta_r)
  }
  y <- data$faecal[t, i]
  if (y >= 0) {
    w[1] <- w[1] * ifelse(y == 1, 0, 1)
    w[2] <- w[2] * ifelse(y == 1, model$theta_f, 1 - model$theta_f)
  }
  w
}

## ---------------------------------------------------------------------------
## Simulate data from known parameters
## ---------------------------------------------------------------------------

true_pars <- list(alpha = 0.01, beta = 0.02, m = 6.0, nu = 0.10,
                  theta_r = 0.8, theta_f = 0.5)

# A plain R simulator for the trajectory. Nothing in the package requires this --
# a real analysis reads observed data instead -- but it gives the example a truth
# to recover.
simulate_cattle <- function(seed = 2024) {
  set.seed(seed)
  X <- matrix(1L, n_timepoints, n_individuals)
  X[1, ] <- ifelse(runif(n_individuals) < true_pars$nu, 2L, 1L)
  for (t in 2:n_timepoints) {
    n_inf <- tapply(X[t - 1, ] == 2L, group, sum)
    for (i in seq_len(n_individuals)) {
      g <- group[i]
      if (X[t - 1, i] == 1L) {
        I_minus <- n_inf[[g]] - 0             # penmates infected at t-1
        p <- -expm1(-(true_pars$alpha + true_pars$beta * I_minus))
        X[t, i] <- if (runif(1) < p) 2L else 1L
      } else {
        X[t, i] <- if (runif(1) < 1 / true_pars$m) 1L else 2L
      }
    }
  }
  X
}

simulate_tests <- function(X, theta) {
  Y <- matrix(-1L, n_timepoints, n_individuals)
  for (t in observed_days) {
    infected <- X[t, ] == 2L
    Y[t, ] <- as.integer(runif(n_individuals) < ifelse(infected, theta, 0))
  }
  Y
}

X_true <- simulate_cattle()
rams   <- simulate_tests(X_true, true_pars$theta_r)
faecal <- simulate_tests(X_true, true_pars$theta_f)

cat(sprintf("Simulated herd: %d pens x %d animals x %d days\n",
            n_pens, n_per_pen, n_timepoints))
cat(sprintf("True infection prevalence: %.3f\n", mean(X_true == 2L)))

## ---------------------------------------------------------------------------
## Fit
## ---------------------------------------------------------------------------

dat <- et_data(
  n_individuals = n_individuals,
  n_timepoints  = n_timepoints,
  transitions   = trans,
  starting_state = starting_state,
  aggregates    = aggs,
  group         = group,
  observation_process = cattle_observations,
  # my own data, reachable as data$rams / data$faecal in the functions above
  extras = list(rams = rams, faecal = faecal))

mod <- et_model(
  data = dat,
  parameters = list(
    alpha   = et_par(et_gamma(1, 1), init = 0.05),
    beta    = et_par(et_gamma(1, 1), init = 0.05),
    m_tilde = et_par(et_gamma(2, 4), init = 4.0),
    nu      = et_par(et_beta(1, 1),  init = 0.10),
    theta_r = et_par(et_beta(1, 1),  init = 0.70),
    theta_f = et_par(et_beta(1, 1),  init = 0.40)),
  # m = m_tilde + 1 keeps 1/m < 1, so the recovery probability can never leave
  # (0, 1) and send the sampler into a DomainError.
  derived = list(m = quote(m_tilde + 1)))

if (!isTRUE(getOption("etr.example.dry_run"))) {
  et_setup()

  n_sweeps <- as.integer(Sys.getenv("CATTLE_SWEEPS", "2000"))
  n_burn   <- as.integer(Sys.getenv("CATTLE_BURN", "800"))
  n_adapts <- as.integer(Sys.getenv("CATTLE_ADAPTS", "400"))

  fit <- et_sample(mod,
    blocks = list(
      et_nuts(c("alpha", "beta", "m_tilde")),
      # The initial infection frequency and the two test sensitivities have
      # closed-form conditional posteriors, so they get conjugate kernels rather
      # than NUTS -- an exact independent draw beats a correlated HMC step.
      et_conjugate_initial_state("nu", at = 1, states = c("S", "I")),
      et_conjugate_test_sensitivity("theta_r", y = "rams",   infected_state = "I"),
      et_conjugate_test_sensitivity("theta_f", y = "faecal", infected_state = "I"),
      et_iffbs("X")),
    n_sweeps = n_sweeps, n_burn = n_burn, n_adapts = n_adapts, seed = 7,
    x_init = X_true)

  print(summary(fit))

  cat("\n=== Posterior recovery ===\n")
  draws <- as.data.frame(fit)
  for (nm in c("alpha", "beta", "m", "nu", "theta_r", "theta_f")) {
    post <- draws[[nm]]
    truth <- true_pars[[nm]]
    mn <- mean(post); s <- stats::sd(post)
    ok <- abs(mn - truth) < 3 * (s + s / sqrt(length(post)))
    cat(sprintf("%-8s posterior mean = %-8.4f (sd %.4f)   truth = %-6g %s\n",
                nm, mn, s, truth, if (ok) "OK" else "MISS"))
  }
}
