# Tier 1: the derived `depends=` annotations.
#
# PracticalBayes' `depends=` is the largest measured speed-up in the stack
# (1.96x on 89% of badger runtime) and its one hazard is that UNDER-declaring is
# silent: a real gradient contribution is dropped with no error and no warning.
# This package derives the annotation from the transpiled bodies rather than
# asking for it, so these tests are the guard on that derivation.

deps <- function(model = toy_model(), blocks = toy_blocks()) {
  et_julia_source(model, blocks)$depends
}

test_that("each term declares exactly the parameters its functions read", {
  d <- deps()
  # alpha, beta from the infection rate; m_tilde via the derived `m` in the
  # recovery rate; nu from the starting state. X because both terms read it.
  expect_setequal(d$epidemic, c("alpha", "beta", "m_tilde", "nu", "X"))
  expect_setequal(d$observation, c("theta", "X"))
})

test_that("a derived parameter resolves to the SAMPLED parameter behind it", {
  # The recovery rate reads `model$m`, which is `m_tilde + 1`. Declaring `m`
  # would be useless (no block owns it) and declaring nothing would be wrong.
  expect_true("m_tilde" %in% deps()$epidemic)
  expect_false("m" %in% deps()$epidemic)
})

test_that("X is always declared: both terms read the trajectory", {
  d <- deps()
  expect_true("X" %in% d$epidemic)
  expect_true("X" %in% d$observation)
})

test_that("the terms are disjoint when the model factorises", {
  d <- deps()
  expect_length(intersect(setdiff(d$epidemic, "X"), setdiff(d$observation, "X")), 0L)
})

test_that("reads reached through a HELPER are declared too", {
  # A rate that reaches a parameter only via a helper still depends on it. This
  # transitivity is what makes the derivation trustworthy rather than a
  # first-order approximation.
  foi_helper <- function(model, i_minus) model$alpha + model$beta * i_minus
  rate <- function(model, data, i, t) {
    -expm1(-foi_helper(model, data$aggregates$n_infected[data$group[i], t]))
  }
  d <- et_data(
    n_individuals = TOY_N, n_timepoints = TOY_T,
    transitions = et_transitions(toy_states, "S -> I" = rate,
                                 "I -> S" = toy_recovery),
    starting_state = toy_start, aggregates = toy_aggregate(),
    group = rep(seq_len(TOY_PENS), each = TOY_PER_PEN),
    observation_weight = toy_obs,
    helpers = list(et_helper(foi_helper, "foi_helper")),
    extras = list(y = toy_y()))
  m <- et_model(data = d, parameters = toy_model()$parameters,
                derived = list(m = quote(m_tilde + 1)))
  g <- et_julia_source(m)
  expect_true(all(c("alpha", "beta") %in% g$depends$epidemic))
  # And the prologue for the ROLE promotes over the helper's reads as well, so
  # the function is typed correctly even though it never names them directly.
  expect_match(g$src, "eltype(model.alpha)", fixed = TRUE)
})

test_that("a coupling-only rate does NOT enter depends", {
  # `coupling_trans_mat` is never seen by epidemic_loglik, so a parameter that
  # appears only there contributes nothing to any gradient. Declaring it would
  # be harmless but wrong-headed; the point of the separate spec is that it is
  # structurally unreachable from AD.
  cached <- function(model, data, i, t) model$kappa * 0.1
  d <- et_data(
    n_individuals = TOY_N, n_timepoints = TOY_T,
    transitions = et_transitions(toy_states, "S -> I" = toy_infection,
                                 "I -> S" = toy_recovery),
    coupling_transitions = et_transitions(toy_states, "S -> I" = cached,
                                          "I -> S" = toy_recovery),
    starting_state = toy_start, aggregates = toy_aggregate(),
    group = rep(seq_len(TOY_PENS), each = TOY_PER_PEN),
    observation_weight = toy_obs, extras = list(y = toy_y()))
  m <- et_model(data = d, parameters = c(
    toy_model()$parameters, list(kappa = et_par(et_beta(1, 1), init = 0.5))),
    derived = list(m = quote(m_tilde + 1)))
  g <- et_julia_source(m)
  expect_false("kappa" %in% g$depends$epidemic)
  expect_match(g$src, "coupling_trans_mat=TRANS_COUPLING", fixed = TRUE)
})

test_that("aggregate updates that read parameters are attributed to the epidemic term", {
  a <- et_aggregate(
    toy_states,
    arrays = list(n_infected = et_array("Float64", c(TOY_PENS, TOY_T))),
    update = function(model, data, X, state, i, t) {
      n_infected[data$group[i], t] <-
        n_infected[data$group[i], t] + model$weight * (state == "I")
    })
  expect_equal(a$reads, "weight")
  d <- et_data(
    n_individuals = TOY_N, n_timepoints = TOY_T,
    transitions = et_transitions(toy_states, "S -> I" = toy_infection,
                                 "I -> S" = toy_recovery),
    starting_state = toy_start, aggregates = a,
    group = rep(seq_len(TOY_PENS), each = TOY_PER_PEN),
    observation_weight = toy_obs, extras = list(y = toy_y()))
  m <- et_model(data = d, parameters = c(
    toy_model()$parameters, list(weight = et_par(et_beta(1, 1), init = 0.5))),
    derived = list(m = quote(m_tilde + 1)))
  expect_true("weight" %in% et_julia_source(m)$depends$epidemic)
})

test_that("the emitted annotation matches the computed set", {
  g <- et_julia_source(toy_model(), toy_blocks())
  epi_line <- grep("@addlogprob! loglik_fn", strsplit(g$src, "\n")[[1]], value = TRUE)
  for (p in g$depends$epidemic) expect_match(epi_line, paste0(":", p), fixed = TRUE)
  obs_line <- grep("@addlogprob! obs_loglik_fn", strsplit(g$src, "\n")[[1]], value = TRUE)
  for (p in g$depends$observation) expect_match(obs_line, paste0(":", p), fixed = TRUE)
  # And nothing extra: the observation term must not name an epidemic parameter.
  expect_false(grepl(":alpha", obs_line, fixed = TRUE))
})
