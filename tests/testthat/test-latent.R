# Weighted chains with sampled latents (et_latent_data), and the transpiler
# pieces they brought in.

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
  expect_equal(jlq(dnorm(x)), "(1 > 0 ? exp(-((x - 0) / 1)^2 / 2) / (1 * sqrt(2pi)) : NaN)")
  expect_equal(jlq(dnorm(x, m)), "(1 > 0 ? exp(-((x - m) / 1)^2 / 2) / (1 * sqrt(2pi)) : NaN)")
  expect_equal(jlq(dnorm(x, m, s)), "(s > 0 ? exp(-((x - m) / s)^2 / 2) / (s * sqrt(2pi)) : NaN)")
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

# ---- latents on a sampled-trajectory model ----------------------------------

toy_latent_data <- function(latents) {
  et_data(
    n_individuals = TOY_N, n_timepoints = TOY_T,
    transitions = et_transitions(toy_states,
      "S -> I" = function(model, data, i, t) plogis(model$alpha + data$latents$u[i]),
      "I -> S" = toy_recovery),
    starting_state = toy_start,
    aggregates = et_aggregate(toy_states,
      arrays = list(n_infected = et_array("Int", c(TOY_PENS, TOY_T))),
      update = function(model, data, X, state, i, t) {
        n_infected[data$group[i], t] <- n_infected[data$group[i], t] + (state == "I")
      }),
    group = rep(seq_len(TOY_PENS), each = TOY_PER_PEN),
    observation_weight = toy_obs, extras = list(y = matrix(-1L, TOY_T, TOY_N)),
    latents = latents)
}

test_that("et_data latents are bound, passed to ET, and updated after the iFFBS sweep", {
  d <- toy_latent_data(et_latents(
    u = et_per_individual(normal_dist(0, 1), init = rep(0, TOY_N)),
    S = et_field(1, 3, prior = normal_dist(0, 2), index = function(t) (t + 3) %/% 4)))
  m <- et_model(d, list(alpha = prior(normal_dist(), init = -2), m = prior(gamma_dist(2, 2), init = 4),
                        nu = prior(beta_dist(), init = 0.1), theta = prior(beta_dist(), init = 0.5)))
  src <- et_julia_source(m)$src
  expect_match(src, "const latents = (; u = _etr_latent_u, S = _etr_latent_S)", fixed = TRUE)
  expect_match(src, "    latents=latents,", fixed = TRUE)
  # with latents the observation term is summed without ET's 1e-12 floor
  expect_match(src, "function et_obsloglik_for(DATA)\n    function obs_loglik(model, data, X)",
               fixed = TRUE)
  expect_match(src, "(rng, model, X) -> (sweep(rng, model, X); et_latent_mh!(rng, model, DATA, X))",
               fixed = TRUE)
  # u is read by a rate only, S by nothing: each is scored by what reads it
  expect_match(src, "_et_individual_sweep!(rng, model, DATA, X, latents.u, Normal(0.0, 1.0), :u, false, true)",
               fixed = TRUE)
  expect_match(src, "_et_windows(et_index_S, 3, 10), :S, false, false)", fixed = TRUE)
  # S had no init: drawn from its prior in Julia
  expect_match(src, "const _etr_latent_S = _etr_draw_coords(Normal(0.0, 2.0), (1, 3, 12))", fixed = TRUE)
})

test_that("et_data refuses latents it cannot update given X", {
  expect_error(toy_latent_data(et_latents(X = et_shared_latent(2, 1, prior = normal_dist()))),
               "shared by every individual")
  expect_error(toy_latent_data(et_latents(
    s = et_per_individual(function(model) bernoulli_dist(model$p), init = rep(0, TOY_N)))),
    "depends on the")
  expect_error(toy_latent_data(et_latents(u = et_per_individual(normal_dist(), init = rep(0, 3)))),
               "length 3")
})

test_that("a shared step and survival give the entry gate a four-argument survival", {
  tr <- et_transitions(c("S", "I", "D"),
    shared = function(model, data, i, t) list(p = model$p),
    survival = et_survival(function(model, data, i, t, shared) shared$p, death = "D"),
    "S -> I" = ~ 0.1)
  d <- et_data(n_individuals = 3, n_timepoints = 4, transitions = tr,
               starting_state = function(model, data, X, i, t) {
                 p <- numeric(data$n_states)
                 p[1] <- 1
                 p
               },
               aggregates = et_aggregate(c("S", "I", "D"), arrays = list(a = et_array("Int", 1)),
                 update = function(model, data, X, state, i, t) a[1] <- a[1] + 0),
               observation_weight = function(model, data, X, i, t, s) 1)
  m <- et_model(d, list(p = prior(beta_dist(), init = 0.9)), entry_time = c(1, 2, 3))
  expect_match(et_julia_source(m)$src,
               "et_survival_fn(model, data, i, t) = et_survival_fn(model, data, i, t, et_shared(model, data, i, t))",
               fixed = TRUE)
})

test_that("et_focal() reads ET's focal index, and only as (data, i)", {
  expect_equal(jlq(data$aggregates$n_alive[g, t] + et_focal(data, i)),
               "data.aggregates.n_alive[g, t] + (data._focal[] == i)")
  expect_error(jlq(et_focal(d, i)), "takes \\(data, i\\)")
  expect_error(jlq(et_focal(data)), "takes \\(data, i\\)")
  expect_error(et_focal(), "transpiled rate")
})

test_that("dnorm brackets compound arguments", {
  expect_equal(jlq(dnorm(x + 1, a - b, s)),
               "(s > 0 ? exp(-(((x + 1) - (a - b)) / s)^2 / 2) / (s * sqrt(2pi)) : NaN)")
})

test_that("transpiled dnorm gives R's values, and NaN rather than an error for a bad sd", {
  skip_without_julia()
  cases <- list(c(0.3, 0, 1), c(-1.2, 0.5, 2), c(4, -1, 0.7), c(1, 2 - 3, 1.5))
  for (k in cases) {
    src <- et_transpile_expr(substitute(dnorm(X, M, S), list(X = k[1], M = k[2], S = k[3])))
    expect_equal(et_eval(src), dnorm(k[1], k[2], k[3]), tolerance = 1e-12)
  }
  expect_true(is.nan(et_eval(et_transpile_expr(quote(dnorm(1, 0, NaN))))))
  expect_true(is.nan(et_eval(et_transpile_expr(quote(dnorm(1, 0, -1))))))
})
