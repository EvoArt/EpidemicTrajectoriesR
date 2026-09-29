# A multistate capture-recapture model fitted with its hidden states summed out.
#
# Animals are seronegative (N), seropositive (P) or dead (D). Each interval a
# live animal survives with probability phi, and a surviving negative converts
# with probability lambda. A live animal is caught with probability p and its
# serostatus read at capture; the first capture is conditioned on.
#
# Nothing couples one animal to another, so each animal's trajectory is an
# independent hidden Markov model given the parameters. likelihood = "marginal"
# sums every trajectory out with the forward algorithm: the fit has no latent
# block, and leave-future-out scoring integrates each animal's state at the
# cutoff and its future path exactly.

library(EpidemicTrajectoriesR)

## ---------------------------------------------------------------------------
## Simulated data
## ---------------------------------------------------------------------------

set.seed(1)
n_animals <- 300; n_times <- 12
truth <- list(phi = 0.8, p = 0.6, lambda = 0.2, nu = 0.25)
first <- sample.int(n_times - 2, n_animals, replace = TRUE)
y <- matrix(0L, n_times, n_animals)       # 0 missed, 1 caught negative, 2 positive
for (i in seq_len(n_animals)) {
  s <- if (runif(1) < truth$nu) 2L else 1L
  for (t in first[i]:n_times) {
    if (t > first[i] && s != 3L) {
      if (runif(1) > truth$phi) s <- 3L
      else if (s == 1L && runif(1) < truth$lambda) s <- 2L
    }
    if (t == first[i] || (s != 3L && runif(1) < truth$p)) y[t, i] <- s
  }
}

## ---------------------------------------------------------------------------
## The model
## ---------------------------------------------------------------------------

states <- c("N", "P", "D")
model <- et_model(
  et_data(
    n_individuals = n_animals, n_timepoints = n_times,
    transitions = et_transitions(
      states,
      "N -> P" = function(model, data, i, t) model$phi * model$lambda,
      "N -> D" = function(model, data, i, t) 1 - model$phi,
      "P -> D" = function(model, data, i, t) 1 - model$phi),
    starting_state = function(model, data, X, i, t) {
      p <- numeric(data$n_states)
      p[1] <- 1 - model$nu
      p[2] <- model$nu
      p
    },
    observation_weight = function(model, data, X, i, t, s) {
      y <- data$y[t, i]
      if (s == 3) return(ifelse(y == 0, 1, 0))
      w <- ifelse(t == data$first[i], 1, ifelse(y > 0, model$p, 1 - model$p))
      if (y == 0) return(w)
      w * ifelse(y == s, 1, 0)
    },
    aggregates = et_aggregate(
      states, arrays = list(n_pos = et_array("Int", c(1L, n_times))),
      update = function(model, data, X, state, i, t) {
        n_pos[1, t] <- n_pos[1, t] + (state == "P")
      }),
    group = seq_len(n_animals),              # one animal per group: no coupling
    sampling_period = cbind(first, n_times),
    extras = list(y = y, first = first)),
  parameters = list(
    phi = prior(beta_dist(2, 2), init = 0.7),
    p = prior(beta_dist(2, 2), init = 0.5),
    lambda = prior(beta_dist(1, 4), init = 0.1),
    nu = prior(beta_dist(1, 4), init = 0.2)),
  likelihood = "marginal")

## ---------------------------------------------------------------------------
## Fit, and score the future exactly
## ---------------------------------------------------------------------------

et_setup()
et_marginal_loglik(model)                       # at the initial values
fit <- et_sample(model, n_sweeps = 1000, n_burn = 500, n_adapts = 500)
sapply(fit$draws, mean)                         # against `truth`

lfo <- et_lfo_cv(model, et_lfo_spec(truncation = et_lfo_truncation(keep = c("y", "first")),
                                    method = "exact_hmm"),
                 L = 6, M = 2, granularity = c("by_individual", "joint"),
                 n_sweeps = 400, n_burn = 200, n_adapts = 200)
et_lfo_windows(lfo)
head(et_lfo_cells(lfo))                         # one row per animal and window
