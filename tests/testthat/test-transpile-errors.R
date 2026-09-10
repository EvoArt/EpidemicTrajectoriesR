# Tier 1: what the transpiler REFUSES.
#
# These matter more than the positive cases. The failure mode this package is
# built to avoid is emitting something plausible for a construct that does not
# mean the same thing in Julia -- a model that runs and is quietly wrong. Every
# rejection below is therefore a deliberate feature, and each error must NAME the
# offending construct so the user can act on it.

test_that("unsupported loops are refused by name", {
  expect_error(jl(quote(while (a) b)), "while", fixed = TRUE)
  expect_error(jl(quote(repeat b)), "repeat", fixed = TRUE)
})

test_that("an unknown function is refused, with a pointer to et_helper()", {
  expect_error(jlq(sapply(x, f)), "sapply")
  expect_error(jlq(my_own_thing(x)), "et_helper")
})

test_that("a declared helper IS allowed, and is recorded", {
  ctx <- nc(helpers = "siler")
  expect_equal(tr(quote(siler(a, b)), ctx), "siler(a, b)")
  expect_equal(ctx$reads$helpers, "siler")
})

test_that("<<- is refused: a transpiled body owns only its locals", {
  expect_error(jl(quote(x <<- 1)), "<<-", fixed = TRUE)
})

test_that("[[ is refused, pointing at [ and $", {
  expect_error(jlq(x[[1]]), "[[", fixed = TRUE)
})

test_that("empty index slots are refused rather than guessed at", {
  # `A[, j]` is a whole-column slice in R. Emitting `A[, j]` in Julia would be a
  # different object (a view vs a copy) in a context that expects a scalar.
  expect_error(jl(quote(A[, j])), "empty index")
})

test_that("a non-symbol call head is refused", {
  expect_error(jl(quote(f(x)(y))), "call head")
})

test_that("named arguments in a transpiled call are refused", {
  ctx <- nc(helpers = "h")
  expect_error(tr(quote(h(a, b = 2)), ctx), "[Nn]amed arguments")
})

test_that("vector literals are refused: they have no scalar meaning", {
  expect_error(jl(c(1, 2, 3)), "vector literals")
})

test_that("a formals mismatch names the expected signature", {
  bad_rate <- function(pars, dat, i, t) 1
  expect_error(role(bad_rate, "rate", "f"),
               "function(model, data, i, t)", fixed = TRUE)
  bad_obs <- function(model, data, X, i, t) 1
  expect_error(role(bad_obs, "obs_weight", "f"), "model, data, X, i, t, s",
               fixed = TRUE)
})

test_that("an unknown role is refused", {
  expect_error(role(function(model, data, i, t) 1, "not_a_role", "f"),
               "unknown role")
})

test_that("a non-function where a function is expected", {
  expect_error(role(42, "rate", "f"), "must be an R function")
})

test_that("an invalid Julia identifier is refused wherever a name is taken", {
  expect_error(et_helper(function(x) x, name = "my.helper"), "identifier")
  expect_error(et_transitions(c("S", "1I"), "S -> I" = function(model, data, i, t) 1),
               "identifier")
})

test_that("a helper using typed literals without `model` is refused clearly", {
  # `et_one()` names the parameter scalar type, which can only be read off
  # `model`. Silently emitting Float64 would reintroduce the very instability the
  # helper was trying to avoid.
  h <- et_helper(function(a, b) { w <- et_one(); w * a * b }, "bad_helper")
  expect_error(EpidemicTrajectoriesR:::et_transpile_helpers(list(h)),
               "does not take")
})
