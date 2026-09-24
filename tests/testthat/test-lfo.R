# Leave-future-out: the R surface, and the Julia it generates.
#
# The generated source is the thing most likely to break, and every way it can
# break is a PARSE error -- an assignment in a keyword slot, an unbalanced
# `begin`, a name only defined in the wrong scope. So the whole of this file is
# tier 1: it builds the source with no Julia session, and parse-checks it with
# `Meta.parseall` if a `julia` binary happens to be on PATH. That costs about a
# second and loads no packages.

trunc_src  <- EpidemicTrajectoriesR:::truncation_src
spec_src   <- EpidemicTrajectoriesR:::lfo_spec_src
inject_src <- EpidemicTrajectoriesR:::lfo_inject_src
seen_src   <- EpidemicTrajectoriesR:::cr_seen_within_src
known_src  <- EpidemicTrajectoriesR:::cr_known_present
cr_ld      <- EpidemicTrajectoriesR:::cr_cell_logdensity
unobs_src  <- EpidemicTrajectoriesR:::cr_unobservable_src

# Parse a Julia fragment; returns TRUE/FALSE, or skips where julia is absent.
julia_parses <- function(src) {
  bin <- Sys.which("julia")
  testthat::skip_if(!nzchar(bin), "no julia binary on PATH")
  f <- tempfile(fileext = ".jl"); writeLines(src, f)
  chk <- tempfile(fileext = ".jl")
  writeLines(sprintf(paste0(
    'ex = Meta.parseall(read(raw"%s", String))\n',
    'bad = filter(x -> x isa Expr && x.head === :error, ex.args)\n',
    'exit(isempty(bad) ? 0 : 1)'), f), chk)
  system2(bin, c("--startup-file=no", shQuote(chk)),
          stdout = FALSE, stderr = FALSE) == 0L
}

cr <- et_capture_recapture(caught = "caught", p = "p")
crt <- et_capture_recapture(caught = "caught", p = "p", infected_states = 2L,
                            tests = et_test("tested", "se", "sp"))
cr2 <- et_capture_recapture(caught = "caught", p = "p", infected_states = 2L,
                            tests = list(et_test("t1", "se1", "sp1"),
                                         et_test("t2", "se2", "sp2")))
plan <- et_lfo_truncation(clamp = "last_seen", keep = "sex")

test_that("truncation declarations become a truncation() call", {
  expect_equal(trunc_src(plan),
               "truncation(clamp=(:last_seen,), keep=(:sex,), strict=true)")
  expect_match(trunc_src(et_lfo_truncation(strict = FALSE)), "strict=false")
})

test_that("a capture-recapture score writes its own observation density", {
  s <- et_lfo_spec(cell_logdensity = cr, truncation = plan)
  src <- spec_src(s, plan_expr = "plan")
  # A capture contradicts a dead trajectory, and only that case is -Inf.
  expect_match(src, "-Inf", fixed = TRUE)
  expect_match(src, "log1p(-model.p)", fixed = TRUE)
  expect_true(julia_parses(src))
})

test_that("constrain_survival adds BOTH halves of the survival correction", {
  s <- et_lfo_spec(cell_logdensity = cr, truncation = plan,
                   constrain_survival = TRUE)
  s$known_present <- known_src(cr, 30, 3, 3) # et_lfo_cv() does this once L, M, stride are known
  src <- spec_src(s, plan_expr = "plan")
  # Constraining without the matching weight is a bias that does not shrink
  # with n_sim, so the two must always appear together.
  expect_match(src, "constrain = constrain", fixed = TRUE)
  expect_match(src, "survival_weight = survival_weight", fixed = TRUE)
  # The pair is destructured BEFORE the call: an assignment is not a legal
  # keyword argument, and inlining it silently produces unparseable Julia.
  expect_match(src, "^constrain, survival_weight = survival_constrained")
  expect_true(julia_parses(src))
})

test_that("the naive arm carries neither half", {
  src <- spec_src(et_lfo_spec(cell_logdensity = cr, truncation = plan),
                  plan_expr = "plan")
  expect_false(grepl("survival_constrained", src, fixed = TRUE))
  expect_false(grepl("constrain =", src, fixed = TRUE))
})

test_that("constraining without knowing who was alive is refused", {
  expect_error(
    et_lfo_spec(cell_logdensity = et_julia("0.0"), truncation = plan,
                constrain_survival = TRUE),
    "who was alive")
})

test_that("a hand-written score still works, and needs et_julia", {
  s <- et_lfo_spec(cell_logdensity = et_julia("X[t, i] == 3 ? -Inf : log(0.5)"),
                   truncation = plan, is_informative = et_julia("true"))
  expect_true(julia_parses(spec_src(s, plan_expr = "plan")))
  expect_error(et_lfo_spec(cell_logdensity = "X[t,i]", truncation = plan),
               "et_capture_recapture")
})

test_that("the constraint looks ahead only to the end of the scored window", {
  # A capture after the window's last occasion is OUTSIDE the block being
  # scored. Conditioning on it deletes a death branch carrying real mass -- a
  # lost positive contribution that no reweighting restores. The bound is the
  # end of the window each step belongs to: a fixed lookahead of M - 1 from
  # every step reaches past it for all steps after the first.
  expect_match(known_src(cr, 4, 2, 2)$src, "_et_seen_to_end_caught_L4_S2_M2",
               fixed = TRUE)
  bin <- Sys.which("julia")
  skip_if(!nzchar(bin), "no julia binary on PATH")
  # One individual, T = 9, caught only at occasion 7. L = 4, M = 2, stride = 2:
  # windows are 5:6 and 7:8. Step 6 ends the first window, so the capture at 7
  # must not make it known present; steps 7 and 8's window does contain it.
  f <- tempfile(fileext = ".jl")
  writeLines(c(
    "DATA = (caught = reshape([0, 0, 0, 0, 0, 0, 1, 0, 0], 9, 1),",
    "        n_timepoints = 9, n_individuals = 1)",
    seen_src(cr, 4, 2, 2),
    "got = vec(_et_seen_to_end_caught_L4_S2_M2)[5:8]",
    "exit(got == [false, false, true, false] ? 0 : 1)"), f)
  expect_equal(system2(bin, c("--startup-file=no", shQuote(f)),
                       stdout = FALSE, stderr = FALSE), 0L)
})

test_that("a test makes the observation depend on the INFECTION state", {
  # Capture alone separates alive from dead, so every infection state scores
  # the same and an LFO comparison between transmission models has nothing to
  # work with. The test factor is what puts transmission into the held-out
  # score at all.
  plain <- cr_ld(cr)$src
  with_test <- cr_ld(crt)$src
  expect_false(grepl("model.se", plain, fixed = TRUE))
  expect_match(with_test, "model.se", fixed = TRUE)
  expect_match(with_test, "X[t, i] in (2, )", fixed = TRUE)
  # Scored only where the individual was actually caught.
  expect_match(with_test, "if y == 1", fixed = TRUE)
  expect_true(julia_parses(spec_src(
    et_lfo_spec(crt, plan), plan_expr = "plan")))
})

test_that("a false positive costs probability, it does not kill the draw", {
  # -Inf on a positive test in an uninfected state would discard the whole draw
  # over one test error. With a specificity it is merely improbable.
  expect_match(cr_ld(crt)$src, "log1p(-model.sp)", fixed = TRUE)
  no_spec <- et_capture_recapture("caught", "p", infected_states = 2L,
                                  tests = et_test("tested", "se"))
  expect_match(cr_ld(no_spec)$src, "-Inf", fixed = TRUE)   # perfect test
})

test_that("several tests each contribute their own factor", {
  # One test cannot identify sensitivity, specificity and prevalence at once,
  # so more than one is the usual case rather than the exotic one.
  src <- cr_ld(cr2)$src
  for (nm in c("data.t1", "data.t2", "model.se1", "model.se2",
               "model.sp1", "model.sp2")) {
    expect_match(src, nm, fixed = TRUE)
  }
  # Both are read only where the individual was caught.
  expect_equal(length(gregexpr("if y == 1", src, fixed = TRUE)[[1]]), 2L)
  expect_true(julia_parses(spec_src(et_lfo_spec(cr2, plan),
                                    plan_expr = "plan")))
})

test_that("a test without the parameters to score it is refused", {
  expect_error(et_capture_recapture("caught", "p",
                                    tests = et_test("tested", "se")),
               "infected_states")
})

test_that("the contradicted states are DERIVED, not restated", {
  # The default asks ET for whichever state its transitions make absorbing --
  # the same function the proposal uses -- so the score and the proposal cannot
  # disagree about which state a capture rules out.
  expect_match(unobs_src(cr), "_absorbing_state(DATA)", fixed = TRUE)
  expect_true(julia_parses(unobs_src(cr)))
})

test_that("a capture can be contradicted by states other than the absorbing one", {
  # Death is only the WORST case: it is absorbing, so one bad step condemns the
  # whole draw. A transient state the observation rules out does the same
  # damage for that cell, so the score must be able to name any set of states.
  multi <- et_capture_recapture("caught", "p", unobservable_states = c(3L, 4L))
  expect_match(unobs_src(multi), "Int[3, 4]", fixed = TRUE)
  expect_match(cr_ld(multi)$src, "X[t, i] in _et_unobservable", fixed = TRUE)
  expect_true(julia_parses(unobs_src(multi)))
})

test_that("a model where every state is observable has no contradiction branch", {
  always <- et_capture_recapture("caught", "p", unobservable_states = NULL)
  src <- cr_ld(always)$src
  expect_false(grepl("-Inf", src, fixed = TRUE))
  expect_false(grepl("_et_unobservable", src, fixed = TRUE))
  expect_true(julia_parses(spec_src(et_lfo_spec(always, plan),
                                    plan_expr = "plan")))
})

test_that("the constraint and the score read the same capture matrix", {
  # Derived, not restated: they cannot drift into disagreeing about who was
  # alive, which would break the co-gating the weight depends on.
  expect_match(known_src(cr, 10, 2, 2)$src, cr$caught, fixed = TRUE)
  expect_match(cr_ld(cr)$src, cr$caught, fixed = TRUE)
})

test_that("the injected fit RETAINS the trajectory", {
  src <- inject_src(toy_model())
  # LFO forward-simulates from every draw's X, so `X` must stay in the output.
  # A normal fit drops it (save_states=(X=:buffer,)), which is why this file
  # steps the sampler by hand rather than calling sample().
  expect_match(src, "Xs[k] = copy(t.X)", fixed = TRUE)
  # Stepped by hand rather than through sample(), whose chain drops X.
  expect_match(src, "AbstractMCMC.step", fixed = TRUE)
  code <- grep("^\\s*#", strsplit(src, "\n")[[1]], value = TRUE, invert = TRUE)
  expect_false(any(grepl("save_states", code, fixed = TRUE)))
  # Defined once per module: a second et_lfo_cv() call must not redefine it.
  expect_match(src, "if !isdefined(@__MODULE__, :et_lfo_fit)", fixed = TRUE)
  expect_true(julia_parses(src))
})

test_that("a model with no observation process omits OBSLOGLIK", {
  m <- toy_model()
  m$data$observation_weight <- NULL
  m$data$observation_process <- NULL
  expect_false(grepl("OBSLOGLIK", inject_src(m), fixed = TRUE))
})

diag_src <- EpidemicTrajectoriesR:::lfo_diagnostics_src

test_that("the diagnostics source is valid Julia once its path is filled in", {
  src <- sprintf(diag_src(), '"/tmp/fit_t0030.jls"')
  expect_true(julia_parses(src))
})

test_that("the diagnostics source reaches FlexiChains through PracticalBayes", {
  src <- diag_src()
  # Not a bare `import FlexiChains`: it is PracticalBayes's dependency, not
  # this project's, so importing it directly fails even though it is loaded.
  expect_false(grepl("import FlexiChains", src, fixed = TRUE))
  expect_match(src, "@eval(PracticalBayes, FlexiChains)", fixed = TRUE)
})

test_that("the diagnostics source unwraps FlexiSummary to a scalar", {
  # ess/rhat/mcse each return a 3-D array per parameter; without `only` the
  # data frame silently gets list columns.
  src <- diag_src()
  for (stat in c("e[FC.Parameter(n)]", "m[FC.Parameter(n)]", "r[FC.Parameter(n)]")) {
    expect_match(src, paste0("only(", stat, ")"), fixed = TRUE)
  }
})

test_that("et_lfo_diagnostics() refuses a missing or empty cache", {
  # Needs a session: the function checks for one before it looks at the path.
  skip_without_julia()
  expect_error(et_lfo_diagnostics(tempfile()), "no such cache directory")
  d <- tempfile(); dir.create(d)
  expect_error(et_lfo_diagnostics(d), "no fit_t", fixed = TRUE)
})

test_that("the injected fit takes x_init and defaults to all-susceptible", {
  src <- inject_src(toy_model())
  # The default must stay all-susceptible: existing callers rely on it.
  expect_match(src, "x_init === nothing ? fill(1,", fixed = TRUE)
  # ... but it must be overridable, or a transmission model's rate parameter is
  # unidentified (n_infected stays 0) and two models differing only in that
  # term fit identically.
  expect_match(src, "x_init=nothing", fixed = TRUE)
  expect_match(src, "Matrix{Int}(x_init)", fixed = TRUE)
  # and et_lfo_cv must forward it rather than swallow it
  expect_match(src, "x_init=x_init", fixed = TRUE)
  expect_true(julia_parses(src))
})
