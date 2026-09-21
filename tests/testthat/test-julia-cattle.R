# Tier 2: the generated module actually loads, runs, and recovers a truth.
# Skipped where there is no Julia session.

test_that("the generated module loads and gives a finite log density", {
  skip_without_julia()
  m <- toy_model()
  ll <- et_loglik(m, toy_blocks())
  expect_named(ll, c("epidemic", "observation"))
  expect_true(all(is.finite(ll)))
})

test_that("the same model reuses its compiled module", {
  skip_without_julia()
  m <- toy_model()
  a <- et_julia_source(m, toy_blocks())
  et_loglik(m, toy_blocks())
  # Second call must not produce a different module name, or every check would
  # recompile the whole stack.
  b <- et_julia_source(m, toy_blocks())
  expect_equal(a$module, b$module)
})

test_that("one iFFBS sweep moves the trajectory", {
  skip_without_julia()
  m <- toy_model()
  X0 <- matrix(1L, TOY_T, TOY_N)
  X1 <- et_iffbs_sweep(m, X0, toy_blocks())
  expect_equal(dim(X1), dim(X0))
  expect_true(all(X1 %in% c(1L, 2L)))
  # Starting from all-susceptible, the filter should place SOME infection: the
  # test data has positives in it, and a positive is impossible in state S.
  expect_gt(sum(X1 == 2L), 0L)
})

test_that("x_init is validated before it reaches Julia", {
  skip_without_julia()
  m <- toy_model()
  expect_error(et_iffbs_sweep(m, matrix(1L, TOY_N, TOY_T), toy_blocks()),
               "timepoints x individuals")
  bad <- matrix(1L, TOY_T, TOY_N); bad[1, 1] <- 7L
  expect_error(et_iffbs_sweep(m, bad, toy_blocks()), "outside 1..2")
})

test_that("et_check_depends passes on the toy model at a realistic trajectory", {
  skip_without_julia()
  m <- toy_model()
  X0 <- et_iffbs_sweep(m, matrix(1L, TOY_T, TOY_N), toy_blocks())
  res <- et_check_depends(m, toy_blocks(), x_init = X0)
  expect_s3_class(res, "data.frame")
  expect_true(all(res$verdict != "UNDER-DECLARED"))
  # Every parameter is checked against every term.
  expect_equal(nrow(res), length(m$par_names) * 2L)
})

test_that("et_check_depends CATCHES a deliberately under-declared term", {
  # The validator has to be able to fail, or a passing run means nothing. The
  # observation weight is replaced with an et_julia() body that reads `alpha`
  # but declares only `theta`, which is exactly the silent-wrong-gradient case.
  skip_without_julia()
  d <- toy_data()
  d$observation_weight <- et_julia(
    "y = data.y[t, i]
     y < 0 && return one(eltype(model.theta))
     s == 1 && return y == 1 ? zero(eltype(model.theta)) : one(eltype(model.theta))
     w = y == 1 ? model.theta : 1 - model.theta
     return w * exp(-model.alpha)", reads = "theta")     # `alpha` NOT declared
  m <- et_model(data = d, parameters = toy_model()$parameters,
                derived = list(m = quote(m_tilde + 1)))
  X0 <- et_iffbs_sweep(toy_model(), matrix(1L, TOY_T, TOY_N), toy_blocks())
  res <- suppressWarnings(et_check_depends(m, toy_blocks(), x_init = X0))
  bad <- res[res$verdict == "UNDER-DECLARED", ]
  expect_gt(nrow(bad), 0L)
  expect_true("alpha" %in% bad$parameter)
  expect_true(all(bad$term == "observation"))
})

test_that("a short fit runs and returns well-shaped draws", {
  skip_without_julia()
  m <- toy_model()
  fit <- et_sample(m, toy_blocks(), n_sweeps = 25, n_burn = 5, n_adapts = 10,
                   seed = 3, quiet = TRUE)
  expect_s3_class(fit, "et_fit")
  expect_setequal(names(fit$draws), c(m$par_names, "m"))
  for (nm in m$par_names) expect_length(fit$draws[[nm]], 25L)
  # A deterministic is recomputed R-side from the sampled draws, exactly.
  expect_equal(fit$draws$m, fit$draws$m_tilde + 1)
  df <- as.data.frame(fit)
  expect_equal(nrow(df), 25L)
  expect_true(all(is.finite(as.matrix(df))))
  expect_s3_class(summary(fit), "summary.et_fit")
})

test_that("draws obey their support", {
  skip_without_julia()
  fit <- et_sample(toy_model(), toy_blocks(), n_sweeps = 25, n_burn = 5,
                   n_adapts = 10, seed = 4, quiet = TRUE)
  expect_true(all(fit$draws$alpha > 0))
  expect_true(all(fit$draws$beta > 0))
  expect_true(all(fit$draws$m > 1))          # m = m_tilde + 1, m_tilde > 0
  expect_true(all(fit$draws$nu >= 0 & fit$draws$nu <= 1))
  expect_true(all(fit$draws$theta >= 0 & fit$draws$theta <= 1))
})

test_that("the same seed gives the same chain", {
  skip_without_julia()
  m <- toy_model()
  a <- et_sample(m, toy_blocks(), n_sweeps = 10, n_adapts = 5, seed = 11, quiet = TRUE)
  b <- et_sample(m, toy_blocks(), n_sweeps = 10, n_adapts = 5, seed = 11, quiet = TRUE)
  expect_equal(a$draws$alpha, b$draws$alpha)
  expect_equal(a$draws$theta, b$draws$theta)
})

test_that("a fit with no blocks at all still works", {
  # The defaulting is what makes a first fit possible with no knowledge of
  # Gibbs blocking; if it ever breaks, the package's headline claim breaks.
  skip_without_julia()
  d <- toy_data()
  m <- et_model(data = d, parameters = list(
    alpha = prior(gamma_dist(1, 1), init = 0.05),
    beta  = prior(gamma_dist(1, 1), init = 0.05),
    m     = prior(gamma_dist(2, 4), init = 5.0),
    nu    = prior(beta_dist(1, 1), init = 0.1),
    theta = prior(beta_dist(1, 1), init = 0.5)))
  fit <- et_sample(m, n_sweeps = 10, n_adapts = 5, seed = 5, quiet = TRUE)
  expect_length(fit$draws$alpha, 10L)
  kinds <- vapply(fit$blocks, function(b) b$kind, character(1))
  expect_setequal(kinds, c("nuts", "iffbs"))
})

test_that("the cattle example recovers its simulated truth", {
  # The full example, at the settings that make it a real check rather than a
  # smoke test. Env-gated because it takes ~1 minute.
  skip_without_julia()
  skip_if_not(identical(Sys.getenv("ETR_RUN_SLOW"), "1"),
              "set ETR_RUN_SLOW=1 to run the full cattle recovery check")
  ex <- system.file("examples", "cattle.R", package = "EpidemicTrajectoriesR")
  skip_if(!nzchar(ex), "example not installed")
  e <- new.env(parent = globalenv())
  withr_seed <- Sys.setenv(CATTLE_SWEEPS = "1500", CATTLE_BURN = "500",
                           CATTLE_ADAPTS = "300")
  sys.source(ex, envir = e)
  draws <- as.data.frame(e$fit)
  for (nm in c("alpha", "beta", "m", "nu", "theta_r", "theta_f")) {
    truth <- e$true_pars[[nm]]
    post <- draws[[nm]]
    se <- stats::sd(post)
    expect_lt(abs(mean(post) - truth), 3 * (se + se / sqrt(length(post))),
              label = paste("posterior mean of", nm))
  }
})
