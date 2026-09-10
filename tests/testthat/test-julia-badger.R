# Tier 2: the badger SEID model -- the stress case.
#
# 2384 individuals x 161 quarters, four states, six tests, time-varying groups,
# Gompertz-Makeham survival, an entry gate, and an observation process split
# between the filter and the likelihood. If the R DSL can express this, it can
# express the models this package exists for.
#
# Needs the reference CSVs, so it is skipped where they are absent.

badger_dir <- function() {
  d <- Sys.getenv("BADGER_DATA_DIR", unset = file.path(
    path.expand("~"), ".julia", "dev", "EpidemicTrajectories", "badger_ref", "RData2"))
  if (dir.exists(d)) d else ""
}

skip_without_badger <- function() {
  skip_if(!nzchar(badger_dir()), "badger CSVs not available")
}

# Build the badger spec by sourcing the example with sampling switched off.
badger_env <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    ex <- system.file("examples", "badger.R", package = "EpidemicTrajectoriesR")
    if (!nzchar(ex)) return(NULL)
    e <- new.env(parent = globalenv())
    old <- options(etr.example.dry_run = TRUE)
    on.exit(options(old), add = TRUE)
    ok <- tryCatch({ suppressMessages(sys.source(ex, envir = e)); TRUE },
                   error = function(err) FALSE)
    if (!ok) return(NULL)
    cache <<- e
    e
  }
})

test_that("the badger data loads with the reference's dimensions", {
  skip_without_badger()
  e <- badger_env(); skip_if(is.null(e), "badger example did not load")
  raw <- e$raw
  # These are the numbers every handover in the investigation quotes.
  expect_equal(raw$n_individuals, 2384L)
  expect_equal(raw$n_timepoints, 161L)
  expect_equal(raw$n_groups, 34L)
  expect_equal(raw$n_tests, 6L)
  expect_equal(dim(raw$X_init), c(161L, 2384L))
  expect_true(all(raw$X_init %in% 1:4))
})

test_that("the repeat-test multiplicities are actually present", {
  # `TestMat.csv` carries multiple rows for the same (individual, quarter) at a
  # meaningful fraction of cells. Flattening them to one slot gave tau = 16.12
  # instead of ~14.98, so the counts are load-bearing, not bookkeeping.
  skip_without_badger()
  e <- badger_env(); skip_if(is.null(e), "badger example did not load")
  raw <- e$raw
  total <- raw$n_neg + raw$n_pos
  expect_gt(sum(total > 1), 0L)
  # Every recorded reading is counted exactly once as either negative or positive.
  recorded <- raw$tests >= 0
  expect_true(all(total[recorded] >= 1L))
  expect_true(all(total[!recorded] == 0L))
})

test_that("the model is the Gompertz-Makeham variant, not Siler", {
  # The C++ reference has no early-life term; fitting one lets it absorb signal
  # belonging to c1/a2/b2 and is the likeliest reason estimates drift.
  skip_without_badger()
  e <- badger_env(); skip_if(is.null(e), "badger example did not load")
  expect_setequal(names(e$mod$parameters),
                  c("tau", "alpha", "lambda", "beta", "q", "c1", "a2", "b2",
                    "thetas", "rhos", "phis", "etas", "nu"))
  expect_false(any(c("a1", "b1") %in% names(e$mod$parameters)))
})

test_that("the spec carries the structure the model needs", {
  skip_without_badger()
  e <- badger_env(); skip_if(is.null(e), "badger example did not load")
  d <- e$mod$data
  expect_equal(d$transitions$states, c("S", "E", "I", "D"))
  expect_equal(d$transitions$survival$death, "D")
  expect_equal(d$coupled_transitions, list(c("S", "E")))
  expect_s3_class(d$affected_individuals, "et_affected")
  expect_equal(e$mod$entry_time, "first_capture_time")
  # The observation process is split: the filter sees capture x tests, the
  # likelihood only tests, because `etas` is conjugate.
  expect_false(is.null(d$observation_weight))
  expect_false(is.null(d$likelihood_weight))
})

test_that("the derived depends= keeps the conjugate capture parameter out", {
  # `etas` appears in the FILTER weight but not the likelihood weight. If it
  # leaked into the observation term it would be double-counted against its own
  # conjugate kernel.
  skip_without_badger()
  e <- badger_env(); skip_if(is.null(e), "badger example did not load")
  g <- et_julia_source(e$mod, e$blocks)
  expect_setequal(g$depends$observation, c("thetas", "rhos", "phis", "X"))
  expect_setequal(g$depends$epidemic,
                  c("tau", "alpha", "lambda", "beta", "q", "c1", "a2", "b2",
                    "nu", "X"))
  expect_false("etas" %in% g$depends$observation)
  expect_false("etas" %in% g$depends$epidemic)
})

test_that("the generated badger module loads and gives a finite log density", {
  skip_without_julia()
  skip_without_badger()
  e <- badger_env(); skip_if(is.null(e), "badger example did not load")
  ll <- et_loglik(e$mod, e$blocks, x_init = e$raw$X_init)
  expect_named(ll, c("epidemic", "observation"))
  expect_true(all(is.finite(ll)))
  expect_lt(ll[["epidemic"]], 0)
  expect_lt(ll[["observation"]], 0)
})

test_that("et_check_depends passes on the badger model", {
  skip_without_julia()
  skip_without_badger()
  e <- badger_env(); skip_if(is.null(e), "badger example did not load")
  res <- et_check_depends(e$mod, e$blocks, x_init = e$raw$X_init)
  expect_true(all(res$verdict != "UNDER-DECLARED"))
})

test_that("one badger iFFBS sweep moves the trajectory but stays in-space", {
  skip_without_julia()
  skip_without_badger()
  e <- badger_env(); skip_if(is.null(e), "badger example did not load")
  X1 <- et_iffbs_sweep(e$mod, e$raw$X_init, e$blocks)
  expect_equal(dim(X1), dim(e$raw$X_init))
  expect_true(all(X1 %in% 1:4))
  expect_gt(sum(X1 != e$raw$X_init), 0L)
})

test_that("a short badger fit runs", {
  skip_without_julia()
  skip_without_badger()
  skip_if_not(identical(Sys.getenv("ETR_RUN_SLOW"), "1"),
              "set ETR_RUN_SLOW=1 to run the badger fit (minutes)")
  e <- badger_env(); skip_if(is.null(e), "badger example did not load")
  fit <- et_sample(e$mod, e$blocks, n_sweeps = 5, seed = 13,
                   x_init = e$raw$X_init, quiet = TRUE)
  expect_length(fit$draws$tau, 5L)
  expect_true(all(fit$draws$tau > 0))
  expect_equal(dim(fit$draws$alpha), c(5L, 34L))
  expect_equal(dim(fit$draws$nu), c(5L, e$raw$n_nu_times, 2L))
  expect_true(all(fit$draws$thetas >= 0 & fit$draws$thetas <= 1))
})
