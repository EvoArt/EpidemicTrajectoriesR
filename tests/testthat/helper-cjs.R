# A small multistate capture-recapture model: seronegative, seropositive, dead.
# Individuals are independent given the parameters, so it can be fitted either
# with iFFBS or with the trajectories summed out, and the two must agree.

cjs_states <- c("N", "P", "D")

# Survival first, then seroconversion among survivors; the first capture is
# conditioned on and always reveals the serostatus.
cjs_sim <- function(n = 60L, n_t = 8L, seed = 1L, phi = 0.8, p = 0.6,
                    lambda = 0.2, nu = 0.25) {
  set.seed(seed)
  first <- sample.int(n_t - 2L, n, replace = TRUE)
  y <- matrix(0L, n_t, n)
  for (i in seq_len(n)) {
    s <- if (runif(1) < nu) 2L else 1L
    for (t in first[i]:n_t) {
      if (t > first[i] && s != 3L) {
        if (runif(1) > phi) s <- 3L
        else if (s == 1L && runif(1) < lambda) s <- 2L
      }
      if (t == first[i] || (s != 3L && runif(1) < p)) y[t, i] <- s
    }
  }
  list(y = y, first = as.integer(first))
}

# A starting trajectory built from the observations alone: the last seen
# status carried forward, dead after the last capture.
cjs_x_init <- function(sim) {
  y <- sim$y
  x <- matrix(1L, nrow(y), ncol(y))
  for (i in seq_len(ncol(y))) {
    s <- y[sim$first[i], i]
    last <- max(which(y[, i] > 0L))
    for (t in sim$first[i]:nrow(y)) {
      if (y[t, i] > 0L) s <- y[t, i]
      x[t, i] <- if (t > last) 3L else s
    }
  }
  x
}

cjs_model <- function(sim, likelihood = "marginal", group = NULL,
                      coupled_transitions = NULL) {
  n_t <- nrow(sim$y); n <- ncol(sim$y)
  et_model(
    et_data(
      n_individuals = n, n_timepoints = n_t,
      transitions = et_transitions(
        cjs_states,
        "N -> P" = function(model, data, i, t) model$phi * model$lambda,
        "N -> D" = function(model, data, i, t) 1 - model$phi,
        "P -> D" = function(model, data, i, t) 1 - model$phi),
      starting_state = function(model, data, X, i, t) {
        p <- numeric(data$n_states)
        p[1] <- 1 - model$nu
        p[2] <- model$nu
        p
      },
      aggregates = et_aggregate(
        cjs_states,
        arrays = list(n_pos = et_array("Int", c(1L, n_t))),
        update = function(model, data, X, state, i, t) {
          n_pos[1, t] <- n_pos[1, t] + (state == "P")
        }),
      group = if (is.null(group)) seq_len(n) else group,
      coupled_transitions = coupled_transitions,
      observation_weight = function(model, data, X, i, t, s) {
        y <- data$y[t, i]
        if (s == 3) return(ifelse(y == 0, 1, 0))
        w <- ifelse(t == data$first[i], 1, ifelse(y > 0, model$p, 1 - model$p))
        if (y == 0) return(w)
        w * ifelse(y == s, 1, 0)
      },
      sampling_period = cbind(sim$first, n_t),
      extras = list(y = sim$y, first = sim$first)),
    parameters = list(
      phi = prior(beta_dist(2, 2), init = 0.7),
      p = prior(beta_dist(2, 2), init = 0.5),
      lambda = prior(beta_dist(1, 4), init = 0.1),
      nu = prior(beta_dist(1, 4), init = 0.2)),
    likelihood = likelihood)
}
