# Fit the cattle model and cache the result.
#
# Both documents call this rather than embedding a fit: a chain that takes
# minutes to produce should not be tied to the lifetime of a render. The chain
# is cached to `docs/cache/`, and the trajectory draws are streamed to disc
# rather than being held in the chain, which is what `save_x` is for.
#
# Re-run with ETR_REFIT=1 to force a refit.

library(EpidemicTrajectoriesR)

fit_cattle <- function(cache_dir = "docs/cache",
                       n_sweeps = 5000, n_burn = 0, n_adapts = 1000,
                       seed = 7, refit = FALSE) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  chain_file <- file.path(cache_dir, "cattle_chain.rds")
  x_file     <- file.path(cache_dir, "cattle_X.jld2")

  if (file.exists(chain_file) && !refit) {
    message("Using cached cattle chain: ", chain_file)
    return(readRDS(chain_file))
  }

  ## -- the model, exactly as the document builds it -------------------------

  n_pens        <- 10
  n_per_pen     <- 8
  n_timepoints  <- 80
  n_individuals <- n_pens * n_per_pen
  group         <- rep(seq_len(n_pens), each = n_per_pen)
  observed_days <- seq(1, n_timepoints, by = 6)
  states        <- c("S", "I")

  aggs <- et_aggregate(
    states,
    arrays = list(n_infected = et_array("Int", c(n_pens, n_timepoints))),
    update = function(model, data, X, state, i, t) {
      n_infected[data$group[i], t] <- n_infected[data$group[i], t] + (state == "I")
    })

  infection <- function(model, data, i, t) {
    I_minus <- data$aggregates$n_infected[data$group[i], t]
    -expm1(-(model$alpha + model$beta * I_minus))
  }
  recovery <- function(model, data, i, t) 1 / model$m

  trans <- et_transitions(states, "S -> I" = infection, "I -> S" = recovery)

  starting_state <- function(model, data, X, i, t) {
    p <- numeric(data$n_states)
    p[1] <- 1 - model$nu
    p[2] <- model$nu
    p
  }

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

  ## -- simulate a truth to recover -----------------------------------------

  true_pars <- list(alpha = 0.01, beta = 0.02, m = 6.0, nu = 0.10,
                    theta_r = 0.8, theta_f = 0.5)

  set.seed(2024)
  X_true <- matrix(1L, n_timepoints, n_individuals)
  X_true[1, ] <- ifelse(runif(n_individuals) < true_pars$nu, 2L, 1L)
  for (t in 2:n_timepoints) {
    n_inf <- tapply(X_true[t - 1, ] == 2L, group, sum)
    for (i in seq_len(n_individuals)) {
      if (X_true[t - 1, i] == 1L) {
        p <- -expm1(-(true_pars$alpha + true_pars$beta * n_inf[[group[i]]]))
        X_true[t, i] <- if (runif(1) < p) 2L else 1L
      } else {
        X_true[t, i] <- if (runif(1) < 1 / true_pars$m) 1L else 2L
      }
    }
  }

  simulate_tests <- function(X, theta) {
    Y <- matrix(-1L, n_timepoints, n_individuals)
    for (t in observed_days) {
      infected <- X[t, ] == 2L
      Y[t, ] <- as.integer(runif(n_individuals) < ifelse(infected, theta, 0))
    }
    Y
  }
  rams   <- simulate_tests(X_true, true_pars$theta_r)
  faecal <- simulate_tests(X_true, true_pars$theta_f)

  dat <- et_data(
    n_individuals = n_individuals,
    n_timepoints  = n_timepoints,
    transitions   = trans,
    starting_state = starting_state,
    aggregates    = aggs,
    group         = group,
    observation_process = cattle_observations,
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
    derived = list(m = quote(m_tilde + 1)))

  ## -- fit ------------------------------------------------------------------

  et_setup()

  t0 <- Sys.time()
  fit <- et_sample(mod,
    blocks = list(
      # AdaptiveHMC adapts the mass matrix and step size during warm-up and then
      # holds them fixed. A dense metric lets it learn the correlation between
      # alpha and beta, which a diagonal one cannot.
      et_adaptive_hmc(c("alpha", "beta", "m_tilde"),
                      target_accept = 0.8, n_steps = 15, metric = "dense"),
      et_conjugate_initial_state("nu", at = 1, states = c("S", "I")),
      et_conjugate_test_sensitivity("theta_r", y = "rams",   infected_state = "I"),
      et_conjugate_test_sensitivity("theta_f", y = "faecal", infected_state = "I"),
      et_iffbs("X")),
    n_sweeps = n_sweeps, n_burn = n_burn, n_adapts = n_adapts, seed = seed,
    x_init = X_true,
    # The trajectory is 80 x 80 integers per draw. Streamed to disc every 10th
    # sweep it stays inspectable without ever sitting in the chain.
    save_x = x_file, save_every = 10)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  out <- list(fit = fit, truth = true_pars, X_true = X_true,
              elapsed = elapsed, x_file = x_file,
              settings = list(n_sweeps = n_sweeps, n_burn = n_burn,
                              n_adapts = n_adapts, seed = seed))
  saveRDS(out, chain_file)
  message(sprintf("Cattle fit cached to %s (%.1f s)", chain_file, elapsed))
  out
}

if (sys.nframe() == 0L) {
  res <- fit_cattle(refit = nzchar(Sys.getenv("ETR_REFIT")))
  print(summary(res$fit))
  cat("\n=== Posterior recovery ===\n")
  draws <- as.data.frame(res$fit)
  for (nm in c("alpha", "beta", "m", "nu", "theta_r", "theta_f")) {
    post <- draws[[nm]]; truth <- res$truth[[nm]]
    mn <- mean(post); s <- stats::sd(post)
    ok <- abs(mn - truth) < 3 * (s + s / sqrt(length(post)))
    cat(sprintf("%-8s mean = %-8.4f (sd %.4f)  truth = %-6g %s\n",
                nm, mn, s, truth, if (ok) "OK" else "MISS"))
  }
}
