# Tier 1: the model-structure DSL and its validation.

# ---- states ------------------------------------------------------------------

test_that("state spaces are validated", {
  expect_error(et_transitions("S", "S -> S" = function(model, data, i, t) 1),
               "at least two")
  expect_error(et_transitions(c("S", "S"), "S -> S" = function(model, data, i, t) 1),
               "duplicate")
})

# ---- transitions -------------------------------------------------------------

test_that("the arrow-string form parses", {
  f <- function(model, data, i, t) 1
  tr <- et_transitions(c("S", "I"), "S -> I" = f, "I -> S" = f)
  expect_length(tr$transitions, 2L)
  expect_equal(tr$transitions[[1]]$from, "S")
  expect_equal(tr$transitions[[1]]$to, "I")
  expect_true(tr$auto_self)
})

test_that("whitespace around the arrow is tolerated", {
  f <- function(model, data, i, t) 1
  tr <- et_transitions(c("S", "I"), "S->I" = f, "  I  ->  S  " = f)
  expect_equal(vapply(tr$transitions, function(x) x$from, character(1)),
               c("S", "I"))
})

test_that("et_trans() objects are accepted alongside the arrow form", {
  f <- function(model, data, i, t) 1
  tr <- et_transitions(c("S", "I"), et_trans("S", "I", f), "I -> S" = f)
  expect_length(tr$transitions, 2L)
})

test_that("transitions naming an unknown state are refused", {
  f <- function(model, data, i, t) 1
  expect_error(et_transitions(c("S", "I"), "S -> R" = f), "not in the state space")
})

test_that("a malformed arrow string is refused", {
  f <- function(model, data, i, t) 1
  expect_error(et_transitions(c("S", "I"), "S => I" = f), "not a transition")
  expect_error(et_transitions(c("S", "I"), f), "arrow")
})

test_that("duplicate transitions are refused", {
  f <- function(model, data, i, t) 1
  expect_error(et_transitions(c("S", "I"), "S -> I" = f, "S -> I" = f),
               "duplicate")
})

test_that("at least one transition is required", {
  expect_error(et_transitions(c("S", "I")), "at least one")
})

test_that("a death state outside the state space is refused", {
  f <- function(model, data, i, t) 1
  expect_error(
    et_transitions(c("S", "I"), "S -> I" = f,
                   survival = et_survival(f, death = "D")),
    "not in the state space")
})

# ---- aggregates --------------------------------------------------------------

test_that("a simple aggregate emits the reversible += form", {
  a <- toy_aggregate()
  expect_equal(a$lines, "n_infected[data.group[i], t] += (state == :I)")
  expect_equal(a$arrays$n_infected$type, "Int")
  expect_equal(a$arrays$n_infected$dim, c(TOY_PENS, TOY_T))
})

test_that("a guarded aggregate keeps its guard", {
  a <- et_aggregate(
    c("S", "I"),
    arrays = list(n = et_array("Int", c(2, 2))),
    update = function(model, data, X, state, i, t) {
      if (data$g[i, t] > 0) {
        n[data$g[i, t], t] <- n[data$g[i, t], t] + (state == "I")
      }
    })
  expect_equal(a$lines,
               "if data.g[i, t] > 0\n    n[data.g[i, t], t] += (state == :I)\nend")
})

test_that("a multiplicative aggregate emits *=", {
  a <- et_aggregate(
    c("S", "I"),
    arrays = list(w = et_array("Float64", c(2, 2))),
    update = function(model, data, X, state, i, t) {
      w[i, t] <- w[i, t] * 2
    })
  expect_equal(a$lines, "w[i, t] *= 2")
})

test_that("an aggregate that is not reversible is refused, not guessed", {
  # The reverse update is DERIVED from the `A[idx] <- A[idx] + x` shape. Any
  # other shape has no derivable reverse, and guessing one would silently break
  # the aggregates-agree-with-X invariant everything downstream relies on.
  expect_error(
    et_aggregate(c("S", "I"), arrays = list(n = et_array("Int", c(2, 2))),
                 update = function(model, data, X, state, i, t) {
                   n[i, t] <- 5
                 }),
    "each update must read")
  expect_error(
    et_aggregate(c("S", "I"), arrays = list(n = et_array("Int", c(2, 2))),
                 update = function(model, data, X, state, i, t) {
                   n[i, t] <- n[t, i] + 1        # indices differ from the LHS
                 }),
    "each update must read")
  expect_error(
    et_aggregate(c("S", "I"), arrays = list(n = et_array("Int", c(2, 2))),
                 update = function(model, data, X, state, i, t) {
                   x <- 1
                 }),
    "must index a declared array")
  expect_error(
    et_aggregate(c("S", "I"), arrays = list(n = et_array("Int", c(2, 2))),
                 update = function(model, data, X, state, i, t) {
                   other[i, t] <- other[i, t] + 1
                 }),
    "not a declared array")
})

test_that("`else` in an aggregate update is refused", {
  expect_error(
    et_aggregate(c("S", "I"), arrays = list(n = et_array("Int", c(2, 2))),
                 update = function(model, data, X, state, i, t) {
                   if (i > 1) { n[i, t] <- n[i, t] + 1 } else { n[i, t] <- n[i, t] - 1 }
                 }),
    "`else` is not supported")
})

test_that("aggregate arrays must be named et_array()s", {
  expect_error(
    et_aggregate(c("S", "I"), arrays = list(et_array("Int", c(2, 2))),
                 update = function(model, data, X, state, i, t) NULL),
    "NAMED list")
  expect_error(
    et_aggregate(c("S", "I"), arrays = list(n = "not an array"),
                 update = function(model, data, X, state, i, t) NULL),
    "et_array")
})

test_that("et_array validates its dimensions", {
  expect_error(et_array("Int", c(2, -1)), "non-negative")
  expect_error(et_array("Int", c(2, 2.5)), "non-negative")
})

# ---- data --------------------------------------------------------------------

test_that("the state space must agree between transitions and aggregates", {
  f <- function(model, data, i, t) 1
  expect_error(
    et_data(n_individuals = 2, n_timepoints = 2,
            transitions = et_transitions(c("S", "I"), "S -> I" = f),
            starting_state = toy_start,
            aggregates = et_aggregate(c("I", "S"),
              arrays = list(n = et_array("Int", c(1, 1))),
              update = function(model, data, X, state, i, t) {
                n[1, 1] <- n[1, 1] + 1
              })),
    "state space differs")
})

test_that("problem-size mismatches are caught at spec time", {
  d <- function(...) {
    f <- function(model, data, i, t) 1
    et_data(transitions = et_transitions(c("S", "I"), "S -> I" = f),
            starting_state = toy_start, aggregates = toy_aggregate(),
            observation_weight = toy_obs, ...)
  }
  expect_error(d(n_individuals = 4, n_timepoints = 2, group = c(1, 1, 2)),
               "has length 3 but there are 4")
  expect_error(d(n_individuals = 4, n_timepoints = 2,
                 sampling_period = matrix(1L, 3, 2)), "4 x 2 matrix")
  expect_error(d(n_individuals = 0, n_timepoints = 2), "positive whole number")
})

test_that("extras are validated", {
  d <- function(extras) {
    f <- function(model, data, i, t) 1
    et_data(n_individuals = 12, n_timepoints = 12,
            transitions = et_transitions(c("S", "I"), "S -> I" = f),
            starting_state = toy_start, aggregates = toy_aggregate(),
            observation_weight = toy_obs, extras = extras)
  }
  expect_error(d(list(1)), "must be named")
  expect_error(d(list(`bad name` = 1)), "identifier")
  # A name ET already defines on `data` would shadow it silently.
  expect_error(d(list(n_states = 1)), "clash")
})

test_that("a model with no observation process warns rather than passing quietly", {
  f <- function(model, data, i, t) 1
  expect_message(
    et_data(n_individuals = 12, n_timepoints = 12,
            transitions = et_transitions(c("S", "I"), "S -> I" = f),
            starting_state = toy_start, aggregates = toy_aggregate()),
    "no observation process")
})

test_that("likelihood_weight requires a filter weight too", {
  f <- function(model, data, i, t) 1
  expect_error(
    et_data(n_individuals = 12, n_timepoints = 12,
            transitions = et_transitions(c("S", "I"), "S -> I" = f),
            starting_state = toy_start, aggregates = toy_aggregate(),
            likelihood_weight = toy_obs),
    "Supply `observation_weight`")
})

test_that("coupled_transitions is validated against the state space", {
  f <- function(model, data, i, t) 1
  mk <- function(ct) {
    et_data(n_individuals = 12, n_timepoints = 12,
            transitions = et_transitions(c("S", "I"), "S -> I" = f),
            starting_state = toy_start, aggregates = toy_aggregate(),
            observation_weight = toy_obs, coupled_transitions = ct)
  }
  expect_error(mk(list(c("S", "R"))), "unknown state")
  expect_error(mk(list("S")), "c\\(from, to\\)")
  expect_s3_class(mk(list(c("S", "I"))), "et_data")
})
