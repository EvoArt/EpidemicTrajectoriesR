# Shared fixtures and the Julia gate.
#
# The suite is in two tiers (DESIGN.md section 7). Tier 1 needs no Julia and is
# the bulk of it; tier 2 needs a session and is skipped where there is none.

tr  <- EpidemicTrajectoriesR:::et_transpile_expr
nc  <- EpidemicTrajectoriesR:::new_ctx
role <- EpidemicTrajectoriesR:::et_transpile_role

# Transpile one expression to Julia source, with a fresh context.
jl <- function(e, ...) tr(e, nc(...))

# Transpile a quoted expression: `jlq(model$alpha + 1)`.
jlq <- function(expr, ...) tr(substitute(expr), nc(...))

julia_available <- function() {
  if (!requireNamespace("JuliaCall", quietly = TRUE)) return(FALSE)
  if (identical(Sys.getenv("ETR_SKIP_JULIA"), "1")) return(FALSE)
  if (EpidemicTrajectoriesR::et_is_ready()) return(TRUE)
  ok <- tryCatch({ EpidemicTrajectoriesR::et_setup(); TRUE },
                 error = function(e) FALSE)
  isTRUE(ok)
}

skip_without_julia <- function() {
  testthat::skip_if_not(julia_available(), "no Julia session")
}

# ---- a small, complete model used across several test files ------------------
#
# Deliberately tiny (4 pens x 3 animals x 12 days) so the tier-2 tests that
# actually compile and run it stay quick.

toy_states <- c("S", "I")
# TOY_T is deliberately NOT equal to TOY_N: with square dimensions a transposed
# trajectory would pass the shape check, and the test for that check would pass
# vacuously.
TOY_PENS <- 4L; TOY_PER_PEN <- 3L; TOY_T <- 10L
TOY_N <- TOY_PENS * TOY_PER_PEN

toy_infection <- function(model, data, i, t) {
  I_minus <- data$aggregates$n_infected[data$group[i], t]
  -expm1(-(model$alpha + model$beta * I_minus))
}
toy_recovery <- function(model, data, i, t) 1 / model$m

toy_start <- function(model, data, X, i, t) {
  p <- numeric(data$n_states)
  p[1] <- 1 - model$nu
  p[2] <- model$nu
  p
}

toy_obs <- function(model, data, X, i, t, s) {
  y <- data$y[t, i]
  if (y < 0) return(1)
  if (s == 1) return(ifelse(y == 1, 0, 1))
  ifelse(y == 1, model$theta, 1 - model$theta)
}

toy_aggregate <- function() {
  et_aggregate(
    toy_states,
    arrays = list(n_infected = et_array("Int", c(TOY_PENS, TOY_T))),
    update = function(model, data, X, state, i, t) {
      n_infected[data$group[i], t] <- n_infected[data$group[i], t] + (state == "I")
    })
}

toy_y <- function(seed = 1) {
  set.seed(seed)
  y <- matrix(-1L, TOY_T, TOY_N)
  for (t in seq(1, TOY_T, by = 3)) y[t, ] <- as.integer(runif(TOY_N) < 0.3)
  y
}

toy_data <- function() {
  et_data(
    n_individuals = TOY_N, n_timepoints = TOY_T,
    transitions = et_transitions(toy_states,
                                 "S -> I" = toy_infection,
                                 "I -> S" = toy_recovery),
    starting_state = toy_start,
    aggregates = toy_aggregate(),
    group = rep(seq_len(TOY_PENS), each = TOY_PER_PEN),
    observation_weight = toy_obs,
    extras = list(y = toy_y()))
}

toy_model <- function() {
  et_model(
    data = toy_data(),
    parameters = list(
      alpha   = prior(gamma_dist(1, 1), init = 0.05),
      beta    = prior(gamma_dist(1, 1), init = 0.05),
      m_tilde = prior(gamma_dist(2, 4), init = 4.0),
      nu      = prior(beta_dist(1, 1), init = 0.1),
      theta   = prior(beta_dist(1, 1), init = 0.6)),
    derived = list(m = quote(m_tilde + 1)))
}

toy_blocks <- function() {
  list(et_nuts(c("alpha", "beta", "m_tilde")),
       et_conjugate_initial_state("nu", at = 1, states = c("S", "I")),
       et_conjugate_test_sensitivity("theta", y = "y", infected_state = "I"),
       et_iffbs("X"))
}
