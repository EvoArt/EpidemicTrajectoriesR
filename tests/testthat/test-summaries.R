# Tier 1: hand-written reversible summaries (et_summary / et_aggregates).
#
# `et_aggregate()` derives the reverse from the update's shape, so it cannot be
# got wrong. `et_summary()` hands that to the user, and the generated Julia has
# to reach the arrays the way ET expects -- through `data.aggregates` -- rather
# than by the bare name the update body uses.

mixed_aggs <- function() {
  et_aggregates(
    c("S", "E", "I", "D"),
    arrays = list(
      n_alive = et_array("Int", c(3, 5)),
      nSE = et_summary("Int", c(3, 5),
        update = function(model, data, X, state, i, t, reverse) {
          g <- data$g[i, t]
          if (g > 0 && t < data$n_timepoints) {
            contrib <- (state == "S") && (X[t + 1, i] == "E")
            if (reverse) {
              nSE[g, t] <- nSE[g, t] - contrib
            } else {
              nSE[g, t] <- nSE[g, t] + contrib
            }
          }
        })),
    update = function(model, data, X, state, i, t) {
      if (data$g[i, t] > 0) {
        n_alive[data$g[i, t], t] <- n_alive[data$g[i, t], t] + (state != "D")
      }
    })
}

test_that("et_aggregates separates derived from hand-written arrays", {
  a <- mixed_aggs()
  expect_equal(a$derived_names, "n_alive")
  expect_equal(a$hand_names, "nSE")
  expect_s3_class(a, "et_aggregate")
})

test_that("aggregate_records captures the pieces, including the guard", {
  a <- mixed_aggs()
  expect_length(a$records, 1L)
  r <- a$records[[1]]
  expect_equal(r$array, "n_alive")
  expect_equal(r$op, "+")
  expect_equal(r$guard, "data.g[i, t] > 0")
  expect_equal(r$contrib, "(state != :D)")
})

test_that("a derived array is emitted as a reversible function, guard intact", {
  a <- mixed_aggs()
  src <- EpidemicTrajectoriesR:::derived_summary_src("n_alive", a$records,
                                                     names(a$arrays))
  expect_match(src, "function et_summary_n_alive(model, data, X, state, i, t, reverse=false)",
               fixed = TRUE)
  expect_match(src, "if data.g[i, t] > 0", fixed = TRUE)
  # Both directions, and the array reached through `data.aggregates` -- the
  # qualification `@aggregate` would otherwise have done.
  expect_match(src, "data.aggregates.n_alive[data.g[i, t], t] += (state != :D)",
               fixed = TRUE)
  expect_match(src, "data.aggregates.n_alive[data.g[i, t], t] -= (state != :D)",
               fixed = TRUE)
})

test_that("a multiplicative update inverts to division, not subtraction", {
  a <- et_aggregates(c("S", "I"),
    arrays = list(w = et_array("Float64", c(2, 2)),
                  h = et_summary("Int", c(2, 2),
                        update = function(model, data, X, state, i, t, reverse) {
                          if (reverse) h[i, t] <- h[i, t] - 1 else h[i, t] <- h[i, t] + 1
                        })),
    update = function(model, data, X, state, i, t) {
      w[i, t] <- w[i, t] * 2
    })
  src <- EpidemicTrajectoriesR:::derived_summary_src("w", a$records, names(a$arrays))
  expect_match(src, "data.aggregates.w[i, t] *= 2", fixed = TRUE)
  expect_match(src, "data.aggregates.w[i, t] /= 2", fixed = TRUE)
})

test_that("a hand-written summary keeps its own body and gets ET's signature", {
  a <- mixed_aggs()
  r <- EpidemicTrajectoriesR:::et_transpile_summary(
    "nSE", a$arrays$nSE, a$states, list(), names(a$arrays))
  expect_equal(r$name, "et_summary_nSE")
  # `reverse` is POSITIONAL with a default: a keyword on a call the compiler
  # cannot resolve forces the kwarg path and allocates per call.
  expect_match(r$src, "function et_summary_nSE(model, data, X, state, i, t, reverse=false)",
               fixed = TRUE)
  expect_match(r$src,
               "data.aggregates.nSE[g, t] = data.aggregates.nSE[g, t] + contrib",
               fixed = TRUE)
  expect_match(r$src,
               "data.aggregates.nSE[g, t] = data.aggregates.nSE[g, t] - contrib",
               fixed = TRUE)
  # It reads X[t+1, i] -- the thing @aggregate cannot express, and the whole
  # reason this path exists. The trajectory holds CODES, so the state name
  # resolves to its index, not to a Symbol: `:E` there would compare an Int to a
  # Symbol and never match, with no error.
  expect_match(r$src, "X[t + 1, i] == 2", fixed = TRUE)
})

test_that("qualification is word-anchored, so a short name cannot corrupt others", {
  q <- EpidemicTrajectoriesR:::qualify_aggregates
  expect_equal(q("n[i, t] += 1", "n"), "data.aggregates.n[i, t] += 1")
  # `nn` merely contains `n` and must be left alone.
  expect_equal(q("nn[i, t] += 1", "n"), "nn[i, t] += 1")
  # An already-qualified name is not double-qualified.
  expect_equal(q("data.aggregates.n[i, t] += 1", "n"),
               "data.aggregates.n[i, t] += 1")
  # A bare mention that is not an index is not a write target.
  expect_equal(q("x = n + 1", "n"), "x = n + 1")
})

test_that("et_summary validates its update's signature", {
  expect_error(
    et_summary("Int", c(2, 2), update = function(model, data, X, state, i, t) NULL),
    "model, data, X, state, i, t, reverse")
  expect_error(et_summary("Int", c(2, 2), update = 42), "must be a function")
  expect_error(et_summary("Int", c(2, -1), update = function(model, data, X, state, i, t, reverse) NULL),
               "non-negative")
})

test_that("et_aggregates refuses inconsistent update declarations", {
  # An et_array() with nowhere to put its update.
  expect_error(
    et_aggregates(c("S", "I"), arrays = list(a = et_array("Int", c(2, 2)))),
    "`update` is required")
  # An update given when every array carries its own.
  expect_error(
    et_aggregates(c("S", "I"),
      arrays = list(h = et_summary("Int", c(2, 2),
                          update = function(model, data, X, state, i, t, reverse) NULL)),
      update = function(model, data, X, state, i, t) NULL),
    "Drop `update`")
  expect_error(
    et_aggregates(c("S", "I"), arrays = list(a = "not an array")),
    "et_array\\(\\) or et_summary\\(\\)")
})

test_that("codegen uses the verbose fallback when any array is hand-written", {
  # The `@aggregate` macro cannot be mixed with plain summaries, so ONE
  # hand-written array switches the whole declaration to ET's fallback.
  toy <- function(aggs) {
    d <- et_data(
      n_individuals = 3, n_timepoints = 5,
      transitions = et_transitions(c("S", "E", "I", "D"),
        "S -> E" = function(model, data, i, t) -expm1(-model$alpha),
        "E -> I" = function(model, data, i, t) model$gamma),
      starting_state = function(model, data, X, i, t) {
        p <- numeric(data$n_states); p[1] <- 1; p
      },
      aggregates = aggs,
      observation_weight = function(model, data, X, i, t, s) 1,
      extras = list(g = matrix(1L, 3, 5)))
    et_model(data = d, parameters = list(
      alpha = prior(beta_dist(1, 1), init = 0.1),
      gamma = prior(beta_dist(1, 1), init = 0.2)))
  }
  src <- et_julia_source(toy(mixed_aggs()))$src
  # A plain NamedTuple of arrays, not the macro.
  expect_match(src, "const AGGS = (; n_alive = zeros(Int, 3, 5)", fixed = TRUE)
  expect_false(grepl("@aggregate", src, fixed = TRUE))
  # Both arrays get a summary function, passed as a Tuple so each keeps its own
  # concrete type through the hot loop.
  expect_match(src, "derived_summaries=(et_summary_n_alive, et_summary_nSE,)",
               fixed = TRUE)
})

test_that("an all-derived declaration still uses the macro", {
  a <- et_aggregates(c("S", "I"),
    arrays = list(n = et_array("Int", c(2, 2))),
    update = function(model, data, X, state, i, t) {
      n[i, t] <- n[i, t] + (state == "I")
    })
  expect_length(a$hand_names, 0L)
  src <- EpidemicTrajectoriesR:::aggregate_src(a)
  expect_match(src, "@aggregate STATES begin", fixed = TRUE)
  expect_match(src, "@array n Int (2, 2)", fixed = TRUE)
})

test_that("a state name compared against the trajectory becomes its CODE", {
  # `state == "I"` is a Symbol comparison; `X[t, i] == "I"` is an integer one.
  # Emitting a Symbol for the second would never match, silently.
  ctx <- EpidemicTrajectoriesR:::new_ctx(state_syms = TRUE,
                                         states = c("S", "E", "I", "D"))
  tr <- EpidemicTrajectoriesR:::et_transpile_expr
  expect_equal(tr(quote(state == "I"), ctx), "state == :I")
  expect_equal(tr(quote(X[t, i] == "I"), ctx), "X[t, i] == 3")
  expect_equal(tr(quote(X[t + 1, i] != "D"), ctx), "X[t + 1, i] != 4")
  expect_equal(tr(quote("E" == X[t, i]), ctx), "2 == X[t, i]")
})

test_that("comparing the trajectory against a non-state is refused", {
  ctx <- EpidemicTrajectoriesR:::new_ctx(state_syms = TRUE,
                                         states = c("S", "I"))
  expect_error(
    EpidemicTrajectoriesR:::et_transpile_expr(quote(X[t, i] == "R"), ctx),
    "'R' is not a state")
})
