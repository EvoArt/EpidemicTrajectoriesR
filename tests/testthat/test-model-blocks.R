# Tier 1: parameters, priors, derived values and Gibbs blocks.

# ---- priors ------------------------------------------------------------------

test_that("priors render as Distributions.jl constructors", {
  expect_equal(EpidemicTrajectoriesR:::dist_to_julia(normal_dist(0, 1)),
               "Normal(0.0, 1.0)")
  expect_equal(EpidemicTrajectoriesR:::dist_to_julia(gamma_dist(2, 4)),
               "Gamma(2.0, 4.0)")
  expect_equal(EpidemicTrajectoriesR:::dist_to_julia(beta_dist(1, 1)),
               "Beta(1.0, 1.0)")
  expect_equal(EpidemicTrajectoriesR:::dist_to_julia(uniform_dist(0, 2)),
               "Uniform(0.0, 2.0)")
  expect_equal(EpidemicTrajectoriesR:::dist_to_julia(dirichlet_dist(c(1, 1, 1))),
               "Dirichlet(Float64[1.0, 1.0, 1.0])")
})

test_that("exponential_dist takes a RATE and emits Julia's SCALE", {
  # Julia's Exponential is parameterised by scale; R users overwhelmingly think
  # in rates. Getting this backwards would be a silent 1/x prior error, so the
  # conversion is explicit and tested.
  expect_equal(EpidemicTrajectoriesR:::dist_to_julia(exponential_dist(1)),
               "Exponential(1.0)")
  expect_equal(EpidemicTrajectoriesR:::dist_to_julia(exponential_dist(1 / 100)),
               "Exponential(100.0)")
  expect_error(exponential_dist(0), "positive")
  expect_error(exponential_dist(-1), "positive")
})

test_that("truncation nests", {
  expect_equal(
    EpidemicTrajectoriesR:::dist_to_julia(truncated_dist(normal_dist(0, 1), lower = 0)),
    "truncated(Normal(0.0, 1.0), 0.0, Inf)")
  expect_equal(
    EpidemicTrajectoriesR:::dist_to_julia(truncated_dist(normal_dist(0, 1), 0, 5)),
    "truncated(Normal(0.0, 1.0), 0.0, 5.0)")
})

test_that("custom_dist is the escape hatch and validates its name", {
  expect_equal(EpidemicTrajectoriesR:::dist_to_julia(custom_dist("Weibull", list(2, 3))),
               "Weibull(2.0, 3.0)")
  expect_error(custom_dist("not a name"), "identifier")
})

# ---- parameters --------------------------------------------------------------

test_that("prior checks init against the declared shape", {
  expect_error(prior(beta_dist(1, 1), init = c(0.1, 0.2), n = 1), "length 2")
  expect_error(prior(beta_dist(1, 1), init = 0.1, n = 3), "length 1")
  expect_error(prior(init = rep(0.1, 5), dim = c(2, 2), kind = "latent"),
               "length 5")
  expect_s3_class(prior(beta_dist(1, 1), init = rep(0.1, 3), n = 3), "prior")
})

test_that("a sampled parameter needs a prior", {
  expect_error(prior(init = 1), "needs a distribution")
})

# ---- model -------------------------------------------------------------------

test_that("the model validates names and derived expressions", {
  d <- toy_data()
  mk <- function(...) et_model(data = d, ...)
  expect_error(mk(parameters = list(alpha = prior(beta_dist(1, 1), init = 0.1)),
                  derived = list(m = quote(nonexistent + 1))),
               "unknown name")
  expect_error(mk(parameters = list(X = prior(beta_dist(1, 1), init = 0.1))),
               "reserved")
  expect_error(mk(parameters = list(a = prior(beta_dist(1, 1), init = 0.1)),
                  derived = list(a = quote(a + 1))),
               "clash")
  expect_error(mk(parameters = list(a = 1)), "prior")
})

test_that("entry_time requires a survival declaration", {
  expect_error(et_model(data = toy_data(),
                        parameters = list(a = prior(beta_dist(1, 1), init = 0.1)),
                        entry_time = rep(1L, TOY_N)),
               "requires the transitions to declare an et_survival")
})

test_that("derived values expand to the SAMPLED parameters behind them", {
  # This is what keeps `depends=` honest when a rate reads `model$m` and `m` is
  # a deterministic function of `m_tilde`.
  m <- toy_model()
  expand <- EpidemicTrajectoriesR:::expand_to_sampled
  expect_equal(expand("m", m), "m_tilde")
  expect_setequal(expand(c("alpha", "m"), m), c("alpha", "m_tilde"))
  # A name that is neither a parameter nor a derived value is a local, and
  # contributes no dependency.
  expect_equal(expand("some_local", m), character())
})

# ---- blocks ------------------------------------------------------------------

test_that("blocks are resolved, with sensible defaults", {
  m <- toy_model()
  rb <- EpidemicTrajectoriesR:::resolve_blocks(m, list())
  kinds <- vapply(rb, function(b) b$kind, character(1))
  # No blocks given: every sampled parameter lands in one NUTS block, and the
  # trajectory gets an iFFBS kernel. That default is what makes a first fit
  # possible without knowing anything about Gibbs blocking.
  expect_setequal(kinds, c("nuts", "iffbs"))
  nuts <- rb[[which(kinds == "nuts")]]
  expect_setequal(nuts$vars, m$par_names)
})

test_that("an omitted iFFBS block is added", {
  m <- toy_model()
  rb <- EpidemicTrajectoriesR:::resolve_blocks(m, list(et_nuts(m$par_names)))
  expect_true(any(vapply(rb, function(b) b$kind == "iffbs", logical(1))))
})

test_that("leftover parameters join a NUTS block", {
  m <- toy_model()
  rb <- EpidemicTrajectoriesR:::resolve_blocks(m, list(et_nuts(c("alpha", "beta"))))
  nuts_vars <- unlist(lapply(rb[vapply(rb, function(b) b$kind == "nuts", logical(1))],
                             function(b) b$vars))
  expect_setequal(nuts_vars, m$par_names)
})

test_that("a parameter in two blocks is refused", {
  m <- toy_model()
  expect_error(
    EpidemicTrajectoriesR:::resolve_blocks(m, list(et_nuts("alpha"), et_nuts("alpha"))),
    "more than one block")
})

test_that("a block naming an unknown parameter is refused", {
  m <- toy_model()
  expect_error(EpidemicTrajectoriesR:::resolve_blocks(m, list(et_nuts("nope"))),
               "unknown parameter")
})

test_that("two iFFBS blocks are refused", {
  m <- toy_model()
  expect_error(
    EpidemicTrajectoriesR:::resolve_blocks(m, list(et_iffbs("X"), et_iffbs("X"))),
    "more than one iFFBS")
})

test_that("a conjugate-owned parameter with no kernel is refused", {
  # Its density is a constant placeholder, so without a kernel nothing would
  # ever inform it -- a silent modelling failure rather than an error.
  d <- toy_data()
  m <- et_model(data = d, parameters = list(
    alpha = prior(gamma_dist(1, 1), init = 0.05),
    nu = prior(init = rep(0.05, 4), dim = c(2, 2), kind = "latent")))
  expect_error(EpidemicTrajectoriesR:::resolve_blocks(m, list()),
               "no conjugate kernel")
})

test_that("a non-block in `blocks` is refused", {
  m <- toy_model()
  expect_error(EpidemicTrajectoriesR:::resolve_blocks(m, list("nuts")),
               "must come from et_nuts")
})

test_that("HMC step sizes expand per parameter, in declaration order", {
  d <- toy_data()
  m <- et_model(data = d, parameters = list(
    a = prior(beta_dist(1, 1), init = 0.1),
    v = prior(beta_dist(1, 1), init = rep(0.1, 3), n = 3)))
  b <- et_hmc(c("a", "v"), n_steps = 5, step_size = list(a = 0.01, v = 0.2))
  expect_equal(EpidemicTrajectoriesR:::expand_step_size(m, b),
               c(0.01, 0.2, 0.2, 0.2))
  # A single number is recycled across the whole block.
  expect_equal(EpidemicTrajectoriesR:::expand_step_size(m, et_hmc(c("a", "v"), 5, 0.1)),
               rep(0.1, 4))
  # A per-element vector is allowed.
  b3 <- et_hmc(c("a", "v"), 5, list(a = 0.01, v = c(0.1, 0.2, 0.3)))
  expect_equal(EpidemicTrajectoriesR:::expand_step_size(m, b3),
               c(0.01, 0.1, 0.2, 0.3))
  # A wrong length is an error, not silent recycling.
  expect_error(EpidemicTrajectoriesR:::expand_step_size(
    m, et_hmc(c("a", "v"), 5, list(a = 0.01, v = c(0.1, 0.2)))), "length 2")
  expect_error(EpidemicTrajectoriesR:::expand_step_size(
    m, et_hmc(c("a", "v"), 5, list(a = 0.01))), "missing entries")
})

test_that("conjugate count bodies must end in the counts", {
  expect_error(
    EpidemicTrajectoriesR:::transpile_count_body(
      function(X, data) { n <- 1 }, "beta"),
    "must be the counts")
  expect_error(
    EpidemicTrajectoriesR:::transpile_count_body(
      function(X, data) c(1, 2, 3), "beta"),
    "c\\(successes, failures\\)")
})

test_that("conjugate count bodies render per family", {
  tcb <- EpidemicTrajectoriesR:::transpile_count_body
  expect_equal(tcb(function(X, data) c(a, b), "beta"), "(a, b)")
  expect_equal(tcb(function(X, data) c(a, b, c), "dirichlet"), "Int[a, b, c]")
  expect_equal(tcb(function(X, data) { n <- 0; c(n, 1) }, "beta"),
               "n = 0\n(n, 1)")
})

# ---- AD backends -------------------------------------------------------------

test_that("AD backends render as ADTypes constructors", {
  a2j <- EpidemicTrajectoriesR:::adtype_to_julia
  expect_equal(a2j("forwarddiff"), "ADTypes.AutoForwardDiff()")
  expect_equal(a2j(et_adtype("forwarddiff", chunksize = 8)),
               "ADTypes.AutoForwardDiff(; chunksize=8)")
  expect_equal(a2j("polyester"),
               "ADTypes.AutoPolyesterForwardDiff(; chunksize=nothing, tag=nothing)")
  expect_equal(a2j("mooncake"), "ADTypes.AutoMooncake()")
  expect_equal(a2j("enzyme"), "ADTypes.AutoEnzyme()")
  expect_equal(a2j("reversediff"), "ADTypes.AutoReverseDiff()")
  expect_error(et_adtype("nope"), "arg")
})

test_that("only the reverse-mode backends need an extra Julia package", {
  ap <- EpidemicTrajectoriesR:::adtype_package
  expect_null(ap("forwarddiff"))
  expect_equal(ap("polyester"), "PolyesterForwardDiff")
  expect_equal(ap("mooncake"), "Mooncake")
})
