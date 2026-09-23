# A leave-future-out refit must see only the data up to its cutoff.
#
# Everything LFO scores rests on that, and no comparison between two of our own
# estimators can check it: a fit that has seen the future still gives finite,
# plausible scores, and gives them to every candidate alike. So the tier-2 tests
# below pin it directly. The tier-1 tests check that the generated code has the
# shape that makes it possible at all.

module_src <- function() et_julia_source(toy_model(), toy_blocks())$src

test_that("every data-dependent artefact is built by a function of the data", {
  src <- module_src()
  for (f in c("et_loglik_for", "et_obsloglik_for", "et_latent_for")) {
    expect_match(src, sprintf("%s(DATA) = ", f), fixed = TRUE)
  }
  expect_match(src, "et_sampler_for(DATA, LATENT!) = Gibbs(", fixed = TRUE)
  # The constants a normal fit uses are those functions applied to DATA.
  expect_match(src, "const LOGLIK = et_loglik_for(DATA)", fixed = TRUE)
  expect_match(src, "const LATENT! = et_latent_for(DATA)", fixed = TRUE)
  expect_match(src, "const SPL = et_sampler_for(DATA, LATENT!)", fixed = TRUE)
})

test_that("the injected fit builds its sampler from the data it is given", {
  src <- EpidemicTrajectoriesR:::lfo_inject_src(toy_model())
  code <- grep("^\\s*#", strsplit(src, "\n")[[1]], value = TRUE, invert = TRUE)
  code <- paste(code, collapse = "\n")
  expect_match(code, "et_the_model(data, et_loglik_for(data), et_obsloglik_for(data))",
               fixed = TRUE)
  expect_match(code, "spl = et_sampler_for(data, et_latent_for(data))", fixed = TRUE)
  # The module constants are built from the full series. The fit must not
  # touch them: SPL's latent block would resample every state against every
  # observation, whatever `data` the model was given.
  expect_false(grepl("\\bSPL\\b", code))
  expect_false(grepl("\\bLOGLIK\\b", code))
  expect_false(grepl("\\bLATENT!", code))
})

# ---- tier 2 ------------------------------------------------------------------

lfo_module <- function() {
  m <- toy_model()
  mod <- EpidemicTrajectoriesR:::et_load_module(et_julia_source(m, toy_blocks()))
  EpidemicTrajectoriesR:::et_lfo_inject(mod, m)
  JuliaCall::julia_command(sprintf(
    "Base.include_string(%s, \"using EpidemicTrajectories: truncation, truncate_data\")",
    mod))
  mod
}

in_module <- function(mod, src) {
  nm <- paste0(mod, "_refit_test_src")
  JuliaCall::julia_assign(nm, src)
  JuliaCall::julia_eval(sprintf("Base.include_string(%s, Main.%s)", mod, nm))
}

test_that("a refit at cutoff t never changes a state after t", {
  skip_without_julia()
  mod <- lfo_module()
  # truncate_data() ends every sampling period at t, and iFFBS writes only
  # inside it, so an honest refit leaves every later state where it started.
  # A sampler built from the full data resamples them.
  ok <- in_module(mod, paste(
    "let t = 5",
    "    X0 = fill(2, DATA.n_timepoints, DATA.n_individuals)",
    "    train = truncate_data(DATA, truncation(keep = (:y,)), t)",
    "    _, Xs = et_lfo_fit(train; n_sweeps = 20, seed = 7, x_init = X0)",
    "    all(X -> X[(t + 1):end, :] == X0[(t + 1):end, :], Xs) &&",
    "        any(X -> X[1:t, :] != X0[1:t, :], Xs)",
    "end", sep = "\n"))
  expect_true(ok)
})

test_that("a refit at cutoff t ignores every observation after t", {
  skip_without_julia()
  mod <- lfo_module()
  # The same fit, once on the real observations and once with everything after
  # the cutoff blanked, must be bit-identical draw for draw. This also reaches
  # the conjugate blocks, which read the observation matrix directly.
  ok <- in_module(mod, paste(
    "let t = 5",
    "    blank = (v, c) -> (w = copy(v); w[(c + 1):end, :] .= -1; w)",
    "    X0 = fill(2, DATA.n_timepoints, DATA.n_individuals)",
    "    a = et_lfo_fit(truncate_data(DATA, truncation(keep = (:y,)), t);",
    "                   n_sweeps = 20, seed = 7, x_init = X0)",
    "    b = et_lfo_fit(truncate_data(DATA, truncation(custom = (:y => blank,)), t);",
    "                   n_sweeps = 20, seed = 7, x_init = X0)",
    "    a[1] == b[1] && a[2] == b[2]",
    "end", sep = "\n"))
  expect_true(ok)
})
