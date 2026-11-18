# Weighted chains with sampled latents (et_latent_data), and the transpiler
# pieces they brought in.

# ---- a small model ------------------------------------------------------------
#
# Two live states would do, but `pre` makes the step into the first timepoint
# do work. 7 individuals x 10 timepoints x 3 years: no two sizes equal.
LT_N <- 7L; LT_T <- 10L; LT_Y <- 3L

lt_year <- function(t) (t + 3) %/% 4          # 4 timepoints a year, the last one short

lt_step <- function(model, data, i, t) {
  if (t == 1) return(NULL)                     # nothing moves into the first timepoint
  phi <- plogis(model$a + data$latents$u[i] + data$latents$S[1, lt_year(t), i] +
                  model$b * data$latents$sex[i])
  list(live = phi, die = 1 - phi)
}

lt_obs <- function(model, data, i, t) {
  y <- data$y[t, i]
  if (y < 0) return(NULL)
  if (y == 1) c(model$p, 0) else c(1 - model$p, 1)
}

lt_fixture <- function(seed = 3) {
  set.seed(seed)
  y <- matrix(sample(c(-1L, 0L, 1L), LT_T * LT_N, replace = TRUE), LT_T, LT_N)
  y[1, ] <- 1L
  S0 <- array(stats::rnorm(LT_Y * LT_N, 0, 0.3), c(1L, LT_Y, LT_N))
  u0 <- stats::rnorm(LT_N, 0, 0.5)
  sex_obs <- c(1, 0, NA, 1, NA, 0, 1)
  sex0 <- c(1, 0, 1, 1, 0, 0, 1)
  list(y = y, S0 = S0, u0 = u0, sex_obs = sex_obs, sex0 = sex0)
}

lt_model <- function(fx = lt_fixture()) {
  tr <- et_transitions(c("alive", "dead"), weights = TRUE, shared = lt_step,
                       "alive -> alive" = ~ shared$live,
                       "alive -> dead"  = ~ shared$die,
                       "dead -> dead"   = 1)
  latents <- et_latents(
    S   = et_field(1, LT_Y, prior = normal_dist(0, 1), index = lt_year, init = fx$S0),
    u   = et_per_individual(normal_dist(0, 1), init = fx$u0),
    sex = et_per_individual(function(model) bernoulli_dist(model$q),
                            observed = fx$sex_obs, init = fx$sex0))
  d <- et_latent_data(
    n_individuals = LT_N, n_timepoints = LT_T, transitions = tr,
    starting_state = function(model, data, i) c(1, 0),
    observation_weights = lt_obs, latents = latents,
    helpers = list(et_helper(lt_year)), extras = list(y = fx$y))
  et_model(d, parameters = list(
    a = prior(normal_dist(0, 2), init = 1.2),
    b = prior(normal_dist(0, 1), init = -0.4),
    p = prior(beta_dist(2, 2), init = 0.6),
    q = prior(beta_dist(1, 1), init = 0.35)))
}

# The same model's per-individual marginal likelihood by the forward algorithm,
# in plain R.
lt_reference <- function(pars, fx) {
  vapply(seq_len(LT_N), function(i) {
    a <- c(1, 0); ll <- 0
    for (t in seq_len(LT_T)) {
      if (t > 1) {
        phi <- stats::plogis(pars$a + fx$u0[i] + fx$S0[1, (t + 3) %/% 4, i] + pars$b * fx$sex0[i])
        a <- c(a[1] * phi, a[1] * (1 - phi) + a[2])
      }
      y <- fx$y[t, i]
      if (y >= 0) a <- a * (if (y == 1) c(pars$p, 0) else c(1 - pars$p, 1))
      s <- sum(a); ll <- ll + log(s); a <- a / s
    }
    ll
  }, numeric(1))
}

# ---- transpiler ---------------------------------------------------------------

test_that("%/% and %% become fld and mod, never div and rem", {
  expect_equal(jlq((t + 7) %/% 8), "fld((t + 7), 8)")
  expect_equal(jlq(k %% 4 + 1), "mod(k, 4) + 1")
})

test_that("plogis/qlogis take one argument only", {
  expect_equal(jlq(plogis(x)), "StatsFuns.logistic(x)")
  expect_equal(jlq(qlogis(model$p)), "StatsFuns.logit(model.p)")
  expect_error(jlq(plogis(x, 1)), "one argument")
})

test_that("dnorm fills in R's defaults and refuses named arguments", {
  expect_equal(jlq(dnorm(x)), "Distributions.pdf(Distributions.Normal(0, 1), x)")
  expect_equal(jlq(dnorm(x, m)), "Distributions.pdf(Distributions.Normal(m, 1), x)")
  expect_equal(jlq(dnorm(x, m, s)), "Distributions.pdf(Distributions.Normal(m, s), x)")
  expect_error(jlq(dnorm(x, m, s, log = TRUE)), "named arguments")
})

test_that("c() is refused outside a weighted-chain role, and must have one entry per state", {
  expect_error(jlq(c(1, 2)), "unsupported function 'c'")
  f <- function(model, data, i) c(model$psi, 1 - model$psi)
  expect_error(role(f, "latent_start", "s", n_states = 3L), "needs 3 values")
  src <- role(f, "latent_start", "s", n_states = 2L)$src
  expect_match(src, "(convert(ETR_T, model.psi), convert(ETR_T, 1 - model.psi))", fixed = TRUE)
})

test_that("the observation role takes `shared` or not, and is always called with it", {
  f4 <- function(model, data, i, t) NULL
  f5 <- function(model, data, i, t, shared) NULL
  expect_match(role(f4, "latent_obs", "o", n_states = 2L)$src, "function o(model, data, i, t, shared)",
               fixed = TRUE)
  expect_match(role(f5, "latent_obs", "o", n_states = 2L)$src, "function o(model, data, i, t, shared)",
               fixed = TRUE)
  expect_error(role(function(model, data, i) NULL, "latent_obs", "o", n_states = 2L),
               "must be declared as")
})

test_that("a shared step's fields come out in one order whatever order they are written in", {
  f1 <- function(model, data, i, t) list(b = 1, a = model$x)
  f2 <- function(model, data, i, t) list(a = model$x, b = 1)
  s1 <- role(f1, "shared", "s")$src; s2 <- role(f2, "shared", "s")$src
  expect_identical(s1, s2)
  expect_match(s1, "(; a = convert(ETR_T, model.x), b = convert(ETR_T, 1))", fixed = TRUE)
})

test_that("et_capture_loglik transpiles to ET's capture_loglik with a centre tuple", {
  expect_equal(
    jlq(et_capture_loglik(data$captures, i, k, x, y, s, r, detectors = data$latents$Xu)),
    "capture_loglik(data.captures, i, k, (x, y), s, r; detectors = data.latents.Xu)")
  expect_equal(jlq(et_capture_loglik(cc, i, k, x, y, s, r, form = "poisson")),
               "capture_loglik(cc, i, k, (x, y), s, r; form = :poisson)")
  expect_error(jlq(et_capture_loglik(cc, i, k, x, s, r)), "positionally")
  expect_error(jlq(et_capture_loglik(cc, i, k, x, y, s, r, form = "binomial")), "hazard")
  expect_error(jlq(et_capture_loglik(cc, i, k, x, y, s, r, sigma = 2)), "positionally")
  expect_error(et_capture_loglik(), "transpiled")
})

test_that("a function without `model` may not read the parameters", {
  pf <- EpidemicTrajectoriesR:::transpile_plain_fn
  expect_match(pf(function(t) (t + 1) %/% 2, "idx", 1L, list(), "index"),
               "function idx(t)\n    fld((t + 1), 2)\nend", fixed = TRUE)
  expect_error(pf(function(t) t * model$k, "idx", 1L, list(), "index"), "cannot read the parameters")
  expect_error(pf(function(a, b) a, "idx", 1L, list(), "index"), "1 argument")
})

test_that("a parameter-dependent prior is one distribution call", {
  pp <- EpidemicTrajectoriesR:::transpile_param_prior
  expect_equal(pp(function(model) bernoulli_dist(model$q), "pr", list(), "prior")$src,
               "pr(model) = Distributions.Bernoulli(model.q)")
  expect_equal(pp(function(model) exponential_dist(model$r), "pr", list(), "prior")$src,
               "pr(model) = Distributions.Exponential(1 / (model.r))")
  expect_error(pp(function(model) model$q, "pr", list(), "prior"), "call to one of")
  expect_error(pp(function(m) bernoulli_dist(m$q), "pr", list(), "prior"), "function\\(model\\)")
})

# ---- transitions ------------------------------------------------------------

test_that("rates may be formulas or constants, and grouped sources expand", {
  tr <- et_transitions(c("A", "B", "D"), weights = TRUE,
                       shared = function(model, data, i, t) list(d = model$m),
                       "(A, B) -> D" = ~ shared$d, "A -> A" = 1, "B -> B" = ~ 1 - shared$d,
                       "D -> D" = 1)
  key <- vapply(tr$transitions, function(x) paste(x$from, x$to), character(1))
  expect_equal(key, c("A D", "B D", "A A", "B B", "D D"))
  expect_equal(names(formals(tr$transitions[[1]]$rate)), c("model", "data", "i", "t", "shared"))
  expect_error(et_transitions(c("A", "B"), "A -> B" = "x"), "formula or a number")
  expect_error(et_transitions(c("A", "B"), "A -> B" = y ~ x), "one-sided")
})

test_that("auto_self = FALSE is refused, and :no_auto_self is never emitted", {
  expect_error(et_transitions(c("S", "I"), "S -> I" = ~ 0.1, auto_self = FALSE),
               "no longer supported")
  tr <- et_transitions(c("S", "I"), "S -> I" = ~ 0.1)
  src <- EpidemicTrajectoriesR:::transitions_src("T", tr, c("S->I" = "f"), NULL)
  expect_false(grepl("auto_self", src))
})

test_that("weights go with et_latent_data only, and not with survival", {
  tr <- et_transitions(c("S", "I"), weights = TRUE, "S -> I" = ~ 0.1, "S -> S" = 1, "I -> I" = 1)
  expect_error(et_data(n_individuals = 2, n_timepoints = 3, transitions = tr,
                       starting_state = toy_start, aggregates = NULL), "belong to")
  expect_error(et_transitions(c("S", "D"), weights = TRUE, "S -> D" = ~ 0.1,
                              survival = et_survival(function(model, data, i, t) 1, "D")),
               "does not go with weights")
})

# ---- declarations -----------------------------------------------------------

test_that("latent declarations are checked against the data", {
  expect_error(et_field(2, 3, prior = normal_dist()), "one of `index` or `touches`")
  expect_error(et_field(2, 3, prior = list(normal_dist()), index = lt_year), "list of 2")
  fx <- lt_fixture()
  m <- lt_model(fx)
  d <- m$data
  bad <- function(...) {
    ls <- d$latents; mods <- list(...)
    for (nm in names(mods)) ls[[nm]] <- mods[[nm]]
    et_latent_data(LT_N, LT_T, d$transitions, d$starting_state, lt_obs,
                   latents = structure(ls, class = "et_latents"), extras = d$extras)
  }
  # a transposed init must not pass: the fixture's sizes all differ
  expect_error(bad(S = et_field(1, LT_Y, prior = normal_dist(), index = lt_year,
                                init = array(0, c(1L, LT_N, LT_Y)))), "must be 1 x 3 x 7")
  expect_error(bad(u = et_per_individual(normal_dist(), init = rep(0, LT_N - 1))), "length 6")
  expect_error(bad(sex = et_per_individual(function(model) bernoulli_dist(model$q))),
               "give `init` or")
  expect_error(et_latent_data(LT_N, LT_T, d$transitions, d$starting_state, lt_obs,
                              latents = d$latents, extras = list(current = 1)), "clash")
})

test_that("spatial captures are checked before they reach Julia", {
  xy <- rbind(c(0, 1, NA), c(0, 1, NA))
  act <- matrix(TRUE, 3, 4)
  caps <- data.frame(i = c(1, 2), occasion = c(1, 3), detector = c(1, 3), n = c(1, 2))
  sc <- et_spatial_captures(xy, act, caps, n_individuals = 5, latent_detectors = 3)
  expect_true(is.nan(sc$detector_xy[1, 3]))
  expect_error(et_spatial_captures(xy, act, caps, n_individuals = 5), "no position")
  expect_error(et_spatial_captures(t(xy), act, caps, 5, 3), "must be 2 x 3")
  act[3, 3] <- FALSE
  expect_error(et_spatial_captures(xy, act, caps, 5, 3), "not active")
  caps$i[1] <- 6
  expect_error(et_spatial_captures(xy, matrix(TRUE, 3, 4), caps, 5, 3), "out of range")
})

test_that("et_model takes latent data, and refuses what it cannot carry", {
  m <- lt_model()
  expect_identical(m$likelihood, "latent")
  expect_error(et_model(m$data, list(a = prior(normal_dist(), init = 0)), entry_time = 1:LT_N),
               "entry_time")
  expect_error(et_model(m$data, list(u = prior(normal_dist(), init = 0))), "clash")
  blocks <- EpidemicTrajectoriesR:::resolve_blocks(m, list())
  expect_equal(vapply(blocks, function(b) b$kind, character(1)), "nuts")
  expect_error(EpidemicTrajectoriesR:::resolve_blocks(m, list(et_iffbs())), "no latent trajectory")
})

# ---- generated source -------------------------------------------------------

test_that("the generated module builds a latent HMM and a latent kernel", {
  m <- lt_model()
  g <- et_julia_source(m, list(et_adaptive_hmc(c("a", "b", "p", "q"))))
  src <- g$src
  expect_match(src, "@transitions STATES :weights begin", fixed = TRUE)
  expect_match(src, "S ~ Field(1, 3; prior=Normal(0.0, 1.0), index=et_index_S, scale=1.0)", fixed = TRUE)
  expect_match(src, "sex ~ PerIndividual(et_prior_sex; observed=_etr_sex_observed, scale=1.0)",
               fixed = TRUE)
  expect_match(src, "const DATA = latent_hmm_data(;", fixed = TRUE)
  expect_match(src, "latents ~ LatentBlock(data)", fixed = TRUE)
  expect_match(src, ":latents => latent_kernel(LATENT!; params=et_params)", fixed = TRUE)
  expect_match(src, "Union{Missing,Int}[k ? v : missing", fixed = TRUE)
  expect_match(src, "const _etr_sex_init = Vector{Int}(vcat(", fixed = TRUE)
  expect_false(grepl("TrajectoryLatent|iffbs", src))
  expect_identical(et_julia_source(m, list(et_adaptive_hmc(c("a", "b", "p", "q"))))$module,
                   g$module)
})

test_that("starting values from the captures are emitted as ET calls", {
  sc <- et_spatial_captures(rbind(c(0, 1, NA), c(0, 1, NA)), matrix(TRUE, 3, 8),
                            data.frame(i = 1, occasion = 2, detector = 3, n = 1),
                            n_individuals = 4, latent_detectors = 3)
  idx <- function(k) (k + 3) %/% 4
  cen <- et_starting_centres(sc, idx, 2, c(-1, 1, -1, 1))
  tr <- et_transitions(c("A", "B"), weights = TRUE, "A -> A" = 1, "B -> B" = 1)
  d <- suppressMessages(et_latent_data(4, 8, tr, function(model, data, i) c(1, 0), NULL,
    latents = et_latents(
      S = et_field(2, 2, prior = normal_dist(), index = idx, init = cen),
      X = et_shared_latent(2, 1, prior = normal_dist(),
                           init = et_starting_detectors(sc, cen, idx, c(-1, 1, -1, 1)))),
    extras = list(cap = sc)))
  src <- et_julia_source(et_model(d, list(a = prior(normal_dist(), init = 0))))$src
  expect_match(src, "S = starting_centres(cap, et_init_index_S, 2; bounds=(-1.0, 1.0, -1.0, 1.0))",
               fixed = TRUE)
  expect_match(src, "X = starting_detectors(cap, starting_centres(cap,", fixed = TRUE)
  # one latent detector is still a vector
  expect_match(src, "const _etr_cap_latent = Vector{Int}(vcat(", fixed = TRUE)
  expect_error(suppressMessages(et_latent_data(4, 8, tr, function(model, data, i) c(1, 0), NULL,
                              latents = et_latents(S = et_field(2, 2, prior = normal_dist(),
                                                                index = idx, init = cen)))),
               "also be passed in `extras`")
})

# ---- tier 2 -----------------------------------------------------------------

skip_without_latent_api <- function() {
  skip_without_julia()
  ok <- isTRUE(et_eval("isdefined(EpidemicTrajectories, :latent_hmm_data)"))
  testthat::skip_if_not(ok, "EpidemicTrajectories has no latent_hmm_data (not the latent-hmm branch)")
}

test_that("fld and mod agree with R's %/% and %% for negative operands", {
  skip_without_julia()
  a <- c(-7, -6, -1, 0, 5, 7, -5.5, 5.5); b <- c(2, -2, 3, -3)
  for (bb in b) {
    jl_q <- et_eval(sprintf("[fld(x, %s) for x in %s]", bb, EpidemicTrajectoriesR:::julia_vector(a, "Float64")))
    jl_r <- et_eval(sprintf("[mod(x, %s) for x in %s]", bb, EpidemicTrajectoriesR:::julia_vector(a, "Float64")))
    expect_equal(as.numeric(jl_q), a %/% bb)
    expect_equal(as.numeric(jl_r), a %% bb)
  }
})

test_that("the latent model's likelihood is the forward algorithm's, plus the parameter-dependent prior", {
  skip_without_latent_api()
  fx <- lt_fixture()
  m <- lt_model(fx)
  pars <- list(a = 0.7, b = 0.3, p = 0.45, q = 0.2)
  ref <- lt_reference(pars, fx)
  expect_equal(et_marginal_loglik(m, pars, by_individual = TRUE), ref, tolerance = 1e-10)
  prior_term <- sum(stats::dbinom(fx$sex0, 1, pars$q, log = TRUE))
  expect_equal(et_marginal_loglik(m, pars), sum(ref) + prior_term, tolerance = 1e-10)
})

test_that("a latent model samples, and its latents come back", {
  skip_without_latent_api()
  m <- lt_model()
  fit <- et_sample(m, list(et_adaptive_hmc(c("a", "b", "p", "q"), n_steps = 4)),
                   n_sweeps = 20, n_adapts = 20, quiet = TRUE)
  expect_length(fit$draws$a, 20)
  expect_true(all(fit$draws$p > 0 & fit$draws$p < 1))
  expect_equal(dim(fit$latents$S), c(1L, LT_Y, LT_N))
  expect_true(all(fit$latents$sex[c(1, 2, 4, 6, 7)] == c(1, 0, 1, 0, 1)))
  expect_error(et_sample(m, x_init = matrix(1L, LT_T, LT_N)), "nothing to start from")
})
