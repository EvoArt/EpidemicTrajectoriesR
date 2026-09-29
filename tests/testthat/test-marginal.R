# likelihood = "marginal": the trajectories summed out by the forward
# algorithm. Tier 1 checks what is generated and what is refused; tier 2 fits
# the same data both ways and checks the posteriors agree.

sim <- cjs_sim()

test_that("a collapsed model generates one likelihood term and no latent block", {
  src <- et_julia_source(cjs_model(sim))$src
  expect_match(src, "epidemic_marginal_loglik(DATA; threads)", fixed = TRUE)
  expect_match(src, "@addlogprob! marginal_fn(pars, data)", fixed = TRUE)
  expect_false(grepl("TrajectoryLatent", src, fixed = TRUE))
  expect_false(grepl("iffbs_kernel", src, fixed = TRUE))
  expect_false(grepl("init = (; X=X0", src, fixed = TRUE))
})

test_that("its depends= is the union of every body the marginal reads", {
  gen <- et_julia_source(cjs_model(sim))
  expect_setequal(gen$depends$marginal, c("phi", "p", "lambda", "nu"))
  expect_false("X" %in% gen$depends$marginal)
})

test_that("the augmented model is generated exactly as before", {
  a <- et_julia_source(cjs_model(sim, likelihood = "augmented"))
  expect_match(a$src, "TrajectoryLatent", fixed = TRUE)
  expect_match(a$src, "iffbs_kernel", fixed = TRUE)
  expect_false(grepl("epidemic_marginal_loglik", a$src, fixed = TRUE))
})

test_that("coupling must be ruled out before a model can be collapsed", {
  expect_error(cjs_model(sim, group = rep(1L, ncol(sim$y))), "independent")
  expect_error(cjs_model(sim, group = rep(1L, ncol(sim$y)),
                         coupled_transitions = list(c("N", "P"))), "independent")
  expect_s3_class(cjs_model(sim, group = rep(1L, ncol(sim$y)),
                            coupled_transitions = list()), "et_model")
})

test_that("what a collapsed model cannot carry is refused", {
  m <- cjs_model(sim)
  expect_error(et_julia_source(m, list(et_iffbs())), "no latent trajectory")
  expect_error(et_julia_source(m, list(et_conjugate_initial_state(
    "nu", states = c("N", "P")))), "no latent trajectory")
  d <- m$data
  expect_error(et_model(d, m$parameters, likelihood = "marginal",
                        entry_time = sim$first), "entry_time")
  d2 <- d; d2$likelihood_weight <- d$observation_weight
  expect_error(et_model(d2, m$parameters, likelihood = "marginal"),
               "likelihood_weight")
})

# ---- tier 2 -----------------------------------------------------------------

test_that("the marginal likelihood is finite, per individual and in total", {
  skip_without_julia()
  m <- cjs_model(sim)
  per <- et_marginal_loglik(m, by_individual = TRUE)
  expect_length(per, ncol(sim$y))
  expect_true(all(is.finite(per)))
  expect_equal(sum(per), et_marginal_loglik(m), tolerance = 1e-10)
  expect_equal(unname(et_loglik(m)), et_marginal_loglik(m), tolerance = 1e-10)
  expect_gt(et_marginal_loglik(m, list(p = 0.6)), et_marginal_loglik(m, list(p = 0.05)))
})

test_that("collapsed and augmented fits give the same posterior", {
  skip_without_julia()
  big <- cjs_sim(n = 200L, n_t = 10L, seed = 3L)
  fm <- et_sample(cjs_model(big), n_sweeps = 1000, n_burn = 300, n_adapts = 300,
                  seed = 1, quiet = TRUE)
  fa <- et_sample(cjs_model(big, "augmented"), n_sweeps = 1000, n_burn = 300,
                  n_adapts = 300, seed = 1, x_init = cjs_x_init(big), quiet = TRUE)
  for (nm in c("phi", "p", "lambda", "nu")) {
    a <- fm$draws[[nm]]; b <- fa$draws[[nm]]
    # a generous bound: both chains are autocorrelated
    expect_lt(abs(mean(a) - mean(b)), 4 * sqrt(var(a) / 100 + var(b) / 100),
              label = nm)
  }
})

test_that("an extra named like a Base function cannot shadow generated code", {
  # The CJS fixture has an extra called `first`; the generated helpers must
  # still reach Base.first, Base.time and friends.
  src <- et_julia_source(cjs_model(sim))$src
  expect_match(src, "const first = ", fixed = TRUE)
  runner <- sub(".*# ---- runner", "", src)
  expect_false(grepl("(?<![A-Za-z0-9_.:])(first|time|length|size)[(]", runner,
                     perl = TRUE))
  expect_match(runner, "Base.time()", fixed = TRUE)
})
