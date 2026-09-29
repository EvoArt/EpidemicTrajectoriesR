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

test_that("an exact LFO spec has no simulation settings, and says why", {
  plan <- et_lfo_truncation(keep = c("y", "first"))
  s <- et_lfo_spec(truncation = plan, method = "exact_hmm")
  expect_equal(s$method, "exact_hmm")
  src <- EpidemicTrajectoriesR:::lfo_spec_src(s, "plan")
  expect_match(src, "scorer = ExactHMM()", fixed = TRUE)
  expect_false(grepl("cell_logdensity", src, fixed = TRUE))

  cr <- et_capture_recapture("y", "p")
  expect_error(et_lfo_spec(cr, plan, method = "exact_hmm"), "cell_logdensity")
  expect_error(et_lfo_spec(truncation = plan, method = "exact_hmm",
                           constrain_survival = TRUE), "no forward proposal")
  expect_error(et_lfo_spec(truncation = plan, method = "exact_hmm", n_sim = 5),
               "n_sim")
})

test_that("a simulated score defaults to the model's own observation weight", {
  plan <- et_lfo_truncation(keep = c("y", "first"))
  s <- et_lfo_spec(truncation = plan)
  expect_null(s$cell_logdensity)
  own <- EpidemicTrajectoriesR:::own_cell_logdensity(cjs_model(sim))
  expect_equal(own$src, "log(data.observation_weight(model, data, X, i, t, X[t, i]))")
  expect_match(et_julia_source(cjs_model(sim))$src, "function et_obs_weight(model, data, X, i, t, s)",
               fixed = TRUE)
})

test_that("a collapsed refit keeps parameter draws only and passes the scorer on", {
  src <- EpidemicTrajectoriesR:::lfo_inject_src(cjs_model(sim))
  expect_match(src, "et_the_model(data, et_marginal_for(data; threads=et_marginal_threads(adtype)))",
               fixed = TRUE)
  expect_false(grepl("t.X", src, fixed = TRUE))
  expect_match(src, "scorer=spec.scorer", fixed = TRUE)
  aug <- EpidemicTrajectoriesR:::lfo_inject_src(cjs_model(sim, "augmented"))
  expect_match(aug, "Xs[k] = copy(t.X)", fixed = TRUE)
  expect_match(aug, "scorer=spec.scorer", fixed = TRUE)
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

test_that("exact LFO runs on a collapsed fit, and its cells add up", {
  skip_without_julia()
  m <- cjs_model(sim)
  plan <- et_lfo_truncation(keep = c("y", "first"))
  res <- et_lfo_cv(m, et_lfo_spec(truncation = plan, method = "exact_hmm"),
                   L = 5L, M = 2L, granularity = c("by_individual", "joint"),
                   n_sweeps = 200, n_burn = 100, n_adapts = 100, quiet = TRUE)
  expect_true(is.finite(et_lfo_elpd(res, "by_individual")))
  cells <- et_lfo_cells(res, "by_individual")
  expect_true(is.integer(cells$cell))
  expect_equal(sum(cells$elpd), et_lfo_elpd(res, "by_individual"), tolerance = 1e-8)
  expect_true(all(cells$n_finite == 200L))
  expect_error(et_lfo_cv(m, et_lfo_spec(truncation = plan, method = "exact_hmm"),
                         L = 5L, M = 2L, granularity = "pointwise"), "whole forecast")
})

test_that("a sweep warm-started from another starts on its tuning, from memory or file", {
  skip_without_julia()
  m <- cjs_model(sim)
  spec <- et_lfo_spec(truncation = et_lfo_truncation(keep = c("y", "first")),
                      method = "exact_hmm")
  run <- function(cutoffs = c(5L, 6L), ...)
    et_lfo_cv(m, spec, L = 5L, M = 2L, granularity = "joint", cutoffs = cutoffs,
              quiet = TRUE, ...)
  cold <- run(n_sweeps = 100, n_burn = 100, n_adapts = 100)
  f <- tempfile(fileext = ".jls")
  expect_error(et_lfo_save_sampler(cold, f), "say which with `cutoff`")
  et_lfo_save_sampler(cold, f, cutoff = 6L)
  from_res  <- run(n_sweeps = 20, warm_start = cold, adapt = "full")
  from_file <- run(n_sweeps = 20, warm_start = f, adapt = "full")
  st <- EpidemicTrajectoriesR:::lfo_sampler_dict_src
  same <- function(a, ta, b, tb) JuliaCall::julia_eval(sprintf(
    "let x = %s[%d].tuning, y = %s[%d].tuning; any(!isnothing, x) && x == y end",
    st(a), ta, st(b), tb))
  # n_adapts = 0: each warm fit ends on exactly the tuning it was handed.
  for (t in 5:6) expect_true(same(from_res, t, cold, t))
  for (t in 5:6) expect_true(same(from_file, t, cold, 6L))
  expect_false(same(cold, 5L, cold, 6L))
  expect_error(run(cutoffs = 7L, n_sweeps = 20, warm_start = cold),
               "no fit at cutoff 7")
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

test_that("the observation-guided proposal is an option of the simulated score only", {
  plan <- et_lfo_truncation(keep = c("y", "first"))
  s <- et_lfo_spec(truncation = plan, observation_guided = TRUE)
  expect_match(EpidemicTrajectoriesR:::lfo_spec_src(
    modifyList(s, list(cell_logdensity = et_julia("0.0"))), "plan"),
    "guide = true", fixed = TRUE)
  expect_error(et_lfo_spec(truncation = plan, observation_guided = TRUE,
                           constrain_survival = TRUE, known_present = "y"), "not both")
  expect_error(et_lfo_spec(truncation = plan, method = "exact_hmm",
                           observation_guided = TRUE), "no forward proposal")
})
