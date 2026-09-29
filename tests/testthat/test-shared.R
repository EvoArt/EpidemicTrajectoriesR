# Shared step values: quantities several rates (or every state's observation
# weight) need, computed once per individual and timepoint. Tier 1 checks what is
# generated; tier 2 checks the likelihood is the one written without them.

sim <- cjs_sim()

# The CJS model of helper-cjs.R, with survival and the capture weight shared.
cjs_model_shared <- function(sim) {
  m <- cjs_model(sim)
  d <- m$data
  d$transitions <- et_transitions(
    cjs_states,
    shared = function(model, data, i, t) list(die = 1 - model$phi),
    "N -> P" = function(model, data, i, t, shared) (1 - shared$die) * model$lambda,
    "N -> D" = function(model, data, i, t, shared) shared$die,
    "P -> D" = function(model, data, i, t, shared) shared$die)
  d$observation_shared <- function(model, data, i, t) {
    y <- data$y[t, i]
    list(y = y, w = ifelse(t == data$first[i], 1, ifelse(y > 0, model$p, 1 - model$p)))
  }
  d$observation_weight <- function(model, data, X, i, t, s, shared) {
    if (s == 3) return(ifelse(shared$y == 0, 1, 0))
    if (shared$y == 0) return(shared$w)
    shared$w * ifelse(shared$y == s, 1, 0)
  }
  et_model(d, m$parameters, likelihood = "marginal")
}

test_that("list() transpiles to a NamedTuple, and only with every element named", {
  expect_equal(et_transpile_expr(quote(list(a = 1 - x, b = exp(y)))),
               "(; a = 1 - x, b = exp(y))")
  expect_error(et_transpile_expr(quote(list(1, b = 2))), "name every element")
  expect_error(et_transpile_expr(quote(list())), "name every element")
})

test_that("shared step functions are emitted and passed to every rate and weight", {
  src <- et_julia_source(cjs_model_shared(sim))$src
  expect_match(src, "function et_shared(model, data, i, t)", fixed = TRUE)
  expect_match(src, "(; die = convert(ETR_T, 1 - model.phi))", fixed = TRUE)
  expect_match(src, "@shared et_shared", fixed = TRUE)
  expect_match(src, "function et_rate_N_P(model, data, i, t, shared)", fixed = TRUE)
  expect_match(src, "function et_obs_weight(model, data, X, i, t, s, shared)",
               fixed = TRUE)
  expect_match(src, "observation_shared=et_obs_shared", fixed = TRUE)
  # The shared function returns its list as written, not converted to a scalar.
  expect_false(grepl("convert(ETR_T, (; die", src, fixed = TRUE))
})

test_that("the default simulated LFO score uses the shared-aware observation wrapper", {
  own <- EpidemicTrajectoriesR:::own_cell_logdensity(cjs_model_shared(sim))
  expect_equal(own$src, "log(data.observation_weight(model, data, X, i, t, X[t, i]))")
})

test_that("a rate that reads only `shared` still depends on what it reads", {
  gen <- et_julia_source(cjs_model_shared(sim))
  src <- gen$src
  # N -> D reads no parameter by name; its element type and depends= come from
  # the shared function's reads.
  expect_match(src, "function et_rate_N_D(model, data, i, t, shared)\n    ETR_T = promote_type(eltype(model.phi))",
               fixed = TRUE)
  expect_setequal(gen$depends$marginal, c("phi", "p", "lambda", "nu"))
})

test_that("rates without the fifth argument are refused when shared is given", {
  m <- cjs_model_shared(sim)
  m$data$transitions$transitions[[1]]$rate <- function(model, data, i, t) model$phi
  expect_error(et_julia_source(m), "must be declared as")
})

test_that("reads through helpers of helpers reach depends=", {
  m <- cjs_model(sim)
  m$data$helpers <- list(et_helper(function(model) model$lambda, name = "inner"),
                         et_helper(function(model) inner(model), name = "outer"))
  m$data$transitions$transitions[[1]]$rate <-
    function(model, data, i, t) model$phi * outer(model)
  m$data$transitions$transitions[[2]]$rate <- function(model, data, i, t) 1 - model$phi
  gen <- et_julia_source(m)
  expect_true("lambda" %in% gen$depends$marginal)
  expect_match(gen$src, "ETR_T = promote_type(eltype(model.phi), eltype(model.lambda))",
               fixed = TRUE)
})

test_that("a marginal model defaults to one forward-mode pass over its parameters", {
  m <- cjs_model(sim)
  a <- EpidemicTrajectoriesR:::default_adtype(m, "forwarddiff")
  expect_equal(a$chunksize, 4)
  expect_identical(EpidemicTrajectoriesR:::default_adtype(m, "mooncake"), "mooncake")
  aug <- cjs_model(sim, likelihood = "augmented")
  expect_identical(EpidemicTrajectoriesR:::default_adtype(aug, "forwarddiff"), "forwarddiff")
})

# ---- tier 2 -----------------------------------------------------------------

test_that("the shared model's likelihood is the plain model's", {
  skip_without_julia()
  plain <- et_marginal_loglik(cjs_model(sim), by_individual = TRUE)
  shared <- et_marginal_loglik(cjs_model_shared(sim), by_individual = TRUE)
  expect_equal(shared, plain, tolerance = 1e-12)
  pars <- list(phi = 0.6, p = 0.3, lambda = 0.3, nu = 0.4)
  expect_equal(et_marginal_loglik(cjs_model_shared(sim), pars),
               et_marginal_loglik(cjs_model(sim), pars), tolerance = 1e-12)
})
