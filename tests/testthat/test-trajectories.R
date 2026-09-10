# The trajectory archive: chunk discovery (tier 1), and the `.rds` round trip
# plus chunk-at-a-time residuals (tier 2).

# ---- tier 1: no Julia needed -------------------------------------------------

chunk_files <- EpidemicTrajectoriesR:::et_chunk_files
with_ext <- EpidemicTrajectoriesR:::et_with_ext

make_chunks <- function(dir, stem, starts, ext, size = 100) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  for (s in starts) {
    writeLines("x", file.path(dir, sprintf("%s_iters_%d_to_%d.%s",
                                           stem, s, s + size - 1, ext)))
  }
}

test_that("chunks are ordered by first iteration, not by name", {
  d <- file.path(tempdir(), "etr-chunks-order")
  unlink(d, recursive = TRUE)
  # Spanning 1..1100 so a string sort would give 1, 1001, 101, 1101, 201, 901 --
  # which would score every trajectory against the wrong parameter draw, and do
  # it silently.
  starts <- c(1L, 101L, 201L, 901L, 1001L, 1101L)
  make_chunks(d, "X", rev(starts), "rds")

  f <- chunk_files(file.path(d, "X.rds"))
  got <- as.integer(sub("^X_iters_([0-9]+)_to_.*$", "\\1", basename(f)))
  expect_identical(got, starts)
})

test_that("discovery matches only this stem and extension", {
  d <- file.path(tempdir(), "etr-chunks-decoy")
  unlink(d, recursive = TRUE)
  make_chunks(d, "X", 1L, "rds")
  make_chunks(d, "Y", 1L, "rds")     # different stem
  make_chunks(d, "X", 1L, "jld2")    # different format

  f <- chunk_files(file.path(d, "X.rds"))
  expect_length(f, 1L)
  expect_match(basename(f), "^X_iters_1_to_100[.]rds$")

  expect_length(chunk_files(file.path(d, "X.jld2")), 1L)
})

test_that("a stem containing regex metacharacters is matched literally", {
  d <- file.path(tempdir(), "etr-chunks-meta")
  unlink(d, recursive = TRUE)
  make_chunks(d, "run.v1", 1L, "rds")
  make_chunks(d, "runXv1", 1L, "rds")   # '.' must not match 'X'

  f <- chunk_files(file.path(d, "run.v1.rds"))
  expect_length(f, 1L)
  expect_match(basename(f), "^run[.]v1_", perl = FALSE)
})

test_that("an absent archive is empty rather than an error", {
  expect_length(chunk_files(file.path(tempdir(), "etr-no-such.rds")), 0L)
})

test_that("extension swapping round-trips", {
  expect_equal(with_ext("a/b/X.rds", "jld2"), "a/b/X.jld2")
  expect_equal(with_ext("a/b/X.jld2", "rds"), "a/b/X.rds")
  expect_equal(with_ext("a/b/X.rds", "rds"), "a/b/X.rds")
})

# ---- tier 2: a real fit, archived and scored ---------------------------------

test_that("save_x = '.rds' leaves an R-readable archive and no Julia files", {
  skip_without_julia()
  d <- file.path(tempdir(), "etr-fit-rds")
  unlink(d, recursive = TRUE)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  x <- file.path(d, "X.rds")

  fit <- et_sample(toy_model(), toy_blocks(), n_sweeps = 20, seed = 1,
                   save_x = x, save_every = 5, quiet = TRUE)

  # The user asked for .rds, so that is what they get -- and the Julia
  # intermediates are gone, which is the whole point of the conversion.
  expect_length(chunk_files(x), 4L)
  expect_length(chunk_files(with_ext(x, "jld2")), 0L)

  # Every chunk is a plain R object: readRDS with no package attached.
  one <- readRDS(chunk_files(x)[1])
  expect_type(one, "list")
  expect_length(one, 5L)                     # save_every = 5
  expect_equal(dim(one[[1]]), c(TOY_T, TOY_N))
  expect_true(all(one[[1]] %in% c(1L, 2L)))

  # And the accessor reads them back in sweep order.
  all_x <- et_trajectories(fit)
  expect_length(all_x, 20L)
  expect_equal(dim(all_x[[1]]), c(TOY_T, TOY_N))

  # Selecting draws reads only what it needs.
  some <- et_trajectories(fit, draws = c(1L, 7L))
  expect_length(some, 2L)
  expect_equal(some[[1]], all_x[[1]])
  expect_equal(some[[2]], all_x[[7]])
})

test_that("et_convert_trajectories recovers an interrupted run's chunks", {
  skip_without_julia()
  d <- file.path(tempdir(), "etr-fit-recover")
  unlink(d, recursive = TRUE)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  x <- file.path(d, "X.rds")

  # convert_x = FALSE leaves the archive exactly as a crashed fit would: the
  # Julia chunks written, nothing converted.
  fit <- et_sample(toy_model(), toy_blocks(), n_sweeps = 20, seed = 1,
                   save_x = x, save_every = 5, convert_x = FALSE, quiet = TRUE)
  expect_length(chunk_files(with_ext(x, "jld2")), 4L)
  expect_length(chunk_files(x), 0L)

  out <- et_convert_trajectories(x, quiet = TRUE)
  expect_length(out, 4L)
  expect_length(chunk_files(x), 4L)
  expect_length(chunk_files(with_ext(x, "jld2")), 0L)

  # Re-running is a no-op, not an error: recovery has to be safe to repeat.
  expect_silent(et_convert_trajectories(x, quiet = TRUE))
  expect_length(chunk_files(x), 4L)
})

test_that("residuals agree whichever format the archive is in", {
  skip_without_julia()
  d <- file.path(tempdir(), "etr-fit-resid")
  unlink(d, recursive = TRUE)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

  spec <- list(et_residual_waiting("S", "I", origin = "window_start",
                                   censor_at = "window_end",
                                   accumulate = "discrete_product",
                                   name = "exposure"))

  # Same seed, same sweeps: the only difference is the on-disk format, so the
  # residuals must match to the bit. This is what pins the chunk-at-a-time
  # reader to the behaviour the eager one had.
  xj <- file.path(d, "J.jld2")
  fj <- et_sample(toy_model(), toy_blocks(), n_sweeps = 20, seed = 4,
                  save_x = xj, save_every = 5, quiet = TRUE)
  rj <- et_residuals(fj, spec, sync_aggregates = TRUE, seed = 9)

  xr <- file.path(d, "R.rds")
  fr <- et_sample(toy_model(), toy_blocks(), n_sweeps = 20, seed = 4,
                  save_x = xr, save_every = 5, quiet = TRUE)
  rr <- et_residuals(fr, spec, sync_aggregates = TRUE, seed = 9)

  expect_equal(nrow(rj), nrow(rr))
  expect_equal(rj$value, rr$value)
  expect_equal(max(rj$draw), 20L)
})

test_that("the chunk boundary does not change the residuals", {
  skip_without_julia()
  d <- file.path(tempdir(), "etr-fit-flush")
  unlink(d, recursive = TRUE)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

  spec <- list(et_residual_waiting("S", "I", origin = "window_start",
                                   censor_at = "window_end",
                                   accumulate = "discrete_product",
                                   name = "exposure"))

  # save_every changes how many files the run is split across, and nothing else.
  # If the reader ever mis-paired a draw with its parameters at a boundary, the
  # two would diverge here.
  a <- file.path(d, "A.rds")
  fa <- et_sample(toy_model(), toy_blocks(), n_sweeps = 20, seed = 5,
                  save_x = a, save_every = 5, quiet = TRUE)
  b <- file.path(d, "B.rds")
  fb <- et_sample(toy_model(), toy_blocks(), n_sweeps = 20, seed = 5,
                  save_x = b, save_every = 20, quiet = TRUE)

  expect_length(chunk_files(a), 4L)
  expect_length(chunk_files(b), 1L)
  expect_equal(et_residuals(fa, spec, seed = 9)$value,
               et_residuals(fb, spec, seed = 9)$value)
})

test_that("scoring an archive never materialises the whole run", {
  skip_without_julia()
  d <- file.path(tempdir(), "etr-fit-mem")
  unlink(d, recursive = TRUE)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  x <- file.path(d, "X.rds")

  fit <- et_sample(toy_model(), toy_blocks(), n_sweeps = 40, seed = 6,
                   save_x = x, save_every = 4, quiet = TRUE)
  expect_length(chunk_files(x), 10L)

  # The archive is read one chunk at a time, so every draw is scored without the
  # run ever being materialised. Ten chunks in, forty draws out.
  #
  # This asserts the OBSERVABLE consequence rather than a memory figure: R's
  # gc() cannot see Julia's heap or JuliaCall's marshalling buffers, so a
  # byte-count here measures the wrong runtime (it read ~1.4 GB on a model whose
  # entire archive is a few MB). The laziness itself is pinned where it can be
  # measured honestly -- one chunk resident, correct order, correct pairing --
  # by the iterator's own unit test in the Julia package.
  r <- et_residuals(fit, list(et_residual_waiting(
    "S", "I", origin = "window_start", censor_at = "window_end",
    accumulate = "discrete_product", name = "exposure")), seed = 9)

  expect_equal(max(r$draw), 40L)
  expect_equal(nrow(r), 40L * TOY_N)
  # Every draw produced a scored row rather than being skipped by a reader that
  # ran out of chunks early.
  expect_setequal(unique(r$draw), seq_len(40L))
})

# ---- et_collect(): residuals scored during the fit, no archive ---------------

test_that("et_collect scores residuals without archiving a trajectory", {
  skip_without_julia()
  spec <- list(et_residual_waiting("S", "I", origin = "window_start",
                                   censor_at = "window_end",
                                   accumulate = "discrete_product",
                                   name = "exposure"))
  out <- et_collect(toy_model(), toy_blocks(), spec,
                    n_sweeps = 20, seed = 4, n_keep = 20, quiet = TRUE)

  expect_s3_class(out$fit, "et_fit")
  expect_equal(nrow(as.data.frame(out$fit)), 20L)
  expect_equal(nrow(out$residuals), 20L * TOY_N)
  expect_true(any(!is.na(out$residuals$value)))
  # Nothing was archived -- that is the whole point of this route.
  expect_null(out$fit$x_path)
})

test_that("et_collect draws the same chain as et_sample at the same seed", {
  skip_without_julia()
  spec <- list(et_residual_waiting("S", "I", origin = "window_start",
                                   censor_at = "window_end",
                                   accumulate = "discrete_product",
                                   name = "exposure"))
  out <- et_collect(toy_model(), toy_blocks(), spec,
                    n_sweeps = 20, seed = 4, n_keep = 20, quiet = TRUE)
  fit <- et_sample(toy_model(), toy_blocks(), n_sweeps = 20, seed = 4,
                   quiet = TRUE)
  # Stepping the sampler by hand must not change what it draws. If this ever
  # diverges, the collector is seeing a different chain than the archive route
  # and the two sets of residuals are not comparable.
  expect_equal(as.data.frame(out$fit)$alpha, as.data.frame(fit)$alpha)
  expect_equal(as.data.frame(out$fit)$beta, as.data.frame(fit)$beta)
})

test_that("et_collect rejects a bad residual list", {
  skip_without_julia()
  expect_error(et_collect(toy_model(), toy_blocks(), list(), quiet = TRUE),
               "non-empty")
  expect_error(et_collect(toy_model(), toy_blocks(), list("nope"), quiet = TRUE),
               "et_residual_survival")
})

# ---- et_reproduction_numbers(): R_i ------------------------------------------

rn_components <- function() {
  list(
    background = function(model, data, i, t) -expm1(-model$alpha),
    transmission = function(model, data, i, t) {
      I_minus <- data$aggregates$n_infected[data$group[i], t]
      -expm1(-(model$beta * I_minus))
    })
}

test_that("et_reproduction_numbers returns one R_i per individual per draw", {
  skip_without_julia()
  d <- file.path(tempdir(), "etr-rn")
  unlink(d, recursive = TRUE)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  x <- file.path(d, "X.rds")
  fit <- et_sample(toy_model(), toy_blocks(), n_sweeps = 20, seed = 4,
                   save_x = x, save_every = 5, quiet = TRUE)

  R <- et_reproduction_numbers(fit, components = rn_components(),
                               secondary = "transmission",
                               infection = c("S", "I"), infectious_state = "I")
  expect_equal(nrow(R), 20L * TOY_N)
  expect_setequal(unique(R$draw), seq_len(20L))
  expect_setequal(unique(R$individual), seq_len(TOY_N))
  # A never-infectious individual gets NA (no opportunity), which is different
  # from having had the opportunity and infected nobody (a genuine 0).
  expect_true(any(!is.na(R$value)))
  expect_true(all(R$value[!is.na(R$value)] >= 0))
})

test_that("et_reproduction_numbers can select draws", {
  skip_without_julia()
  d <- file.path(tempdir(), "etr-rn-sel")
  unlink(d, recursive = TRUE)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  x <- file.path(d, "X.rds")
  fit <- et_sample(toy_model(), toy_blocks(), n_sweeps = 20, seed = 4,
                   save_x = x, save_every = 5, quiet = TRUE)

  R <- et_reproduction_numbers(fit, components = rn_components(),
                               secondary = "transmission",
                               infection = c("S", "I"), infectious_state = "I",
                               draws = c(1L, 5L, 9L))
  expect_equal(nrow(R), 3L * TOY_N)
})

test_that("et_reproduction_numbers refuses an undeclared secondary", {
  skip_without_julia()
  d <- file.path(tempdir(), "etr-rn-bad")
  unlink(d, recursive = TRUE)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  x <- file.path(d, "X.rds")
  fit <- et_sample(toy_model(), toy_blocks(), n_sweeps = 5, seed = 4,
                   save_x = x, save_every = 5, quiet = TRUE)

  expect_error(
    et_reproduction_numbers(fit, components = rn_components(),
                            secondary = "not_a_component",
                            infection = c("S", "I"), infectious_state = "I"),
    "must name one of")
  # An unnamed component list cannot say which part is transmission.
  expect_error(
    et_reproduction_numbers(fit, components = list(function(m, d, i, t) 1),
                            secondary = "transmission",
                            infection = c("S", "I"), infectious_state = "I"),
    "NAMED list")
})

test_that("R_i needs an archive", {
  skip_without_julia()
  fit <- et_sample(toy_model(), toy_blocks(), n_sweeps = 5, seed = 1,
                   quiet = TRUE)
  expect_error(
    et_reproduction_numbers(fit, components = rn_components(),
                            secondary = "transmission",
                            infection = c("S", "I"), infectious_state = "I"),
    "did not archive")
})
