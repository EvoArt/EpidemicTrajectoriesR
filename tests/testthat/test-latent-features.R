# What et_latent_data() and et_model() pass on to ET's latent HMM beyond the
# weights: values computed once per evaluation, observation values shared by a
# key, the cache over them, extra log-density terms, and per-element priors.
# The fixtures (lt_model and friends) are in test-latent.R.

# The lt_model, rewritten so the same likelihood goes through per_evaluation
# and a keyed observation_shared. `p` reaches the observations only through
# per_evaluation's `p_obs`, and `a` the weights only through its `a_step`.
lt_feature_model <- function(fx = lt_fixture(), cache = TRUE, key = TRUE,
                             parameters = NULL, ...) {
  step <- function(model, data, i, t) {
    if (t == 1) return(NULL)
    phi <- plogis(model$a_step + data$latents$u[i] + data$latents$S[1, lt_year(t), i] +
                    model$b * data$latents$sex[i])
    list(live = phi, die = 1 - phi)
  }
  tr <- et_transitions(c("alive", "dead"), weights = TRUE, shared = step,
                       "alive -> alive" = ~ shared$live,
                       "alive -> dead"  = ~ shared$die,
                       "dead -> dead"   = 1)
  latents <- et_latents(
    S   = et_field(1, LT_Y, prior = normal_dist(0, 1), index = lt_year, init = fx$S0),
    u   = et_per_individual(normal_dist(0, 1), init = fx$u0),
    sex = et_per_individual(function(model) bernoulli_dist(model$q),
                            observed = fx$sex_obs, init = fx$sex0))
  per_eval <- function(model, data) {
    w <- numeric(2)
    w[1] <- model$p
    w[2] <- 1 - model$p
    list(a_step = model$a + 0, p_obs = w, n_two = 2)
  }
  obs_shared <- function(model, data, i, key) list(caught = model$p_obs[1], missed = model$p_obs[2])
  obs <- function(model, data, i, t, shared) {
    y <- data$y[t, i]
    if (y < 0) return(NULL)
    if (y == 1) c(shared$caught, 0) else c(shared$missed, 1)
  }
  d <- et_latent_data(
    n_individuals = LT_N, n_timepoints = LT_T, transitions = tr,
    starting_state = function(model, data, i) c(1, 0),
    observation_weights = obs, latents = latents,
    helpers = list(et_helper(lt_year)), extras = list(y = fx$y),
    observation_shared = obs_shared,
    shared_key = if (key) lt_year, n_keys = if (key) LT_Y,
    per_evaluation = per_eval, cache = cache)
  et_model(d, parameters = parameters %||% list(
    a = prior(normal_dist(0, 2), init = 1.2),
    b = prior(normal_dist(0, 1), init = -0.4),
    p = prior(beta_dist(2, 2), init = 0.6),
    q = prior(beta_dist(1, 1), init = 0.35)), ...)
}

lt_args <- function(...) {
  fx <- lt_fixture()
  tr <- et_transitions(c("alive", "dead"), weights = TRUE,
                       "alive -> alive" = 1, "dead -> dead" = 1)
  base <- list(n_individuals = LT_N, n_timepoints = LT_T, transitions = tr,
               starting_state = function(model, data, i) c(1, 0),
               observation_weights = function(model, data, i, t, shared) NULL,
               latents = et_latents(u = et_per_individual(normal_dist(0, 1), init = fx$u0)))
  utils::modifyList(base, list(...))
}

# ---- et_latent_data: accepted and refused --------------------------------------

test_that("observation_shared, shared_key and per_evaluation are checked on the way in", {
  f4 <- function(model, data, i, key) 1
  expect_s3_class(do.call(et_latent_data, lt_args(observation_shared = f4)), "et_latent_data")
  expect_error(suppressMessages(do.call(et_latent_data,
                 lt_args(observation_shared = f4, observation_weights = NULL))),
               "give those too")
  expect_error(do.call(et_latent_data, lt_args(observation_shared = 3)), "function\\(model, data")
  expect_error(do.call(et_latent_data, lt_args(shared_key = function(t) t, n_keys = 3)),
               "give that too")
  expect_error(do.call(et_latent_data, lt_args(observation_shared = f4,
                                               shared_key = function(t) t)), "together")
  expect_error(do.call(et_latent_data, lt_args(observation_shared = f4,
                                               shared_key = function(t, u) t, n_keys = 3)),
               "function of the timepoint")
  expect_error(do.call(et_latent_data, lt_args(per_evaluation = list(a = 1))),
               "function\\(model, data\\)")
})

test_that("cache is TRUE, FALSE or the names by hand, and only known latents", {
  expect_null(do.call(et_latent_data, lt_args(cache = FALSE))$cache)
  expect_true(do.call(et_latent_data, lt_args(cache = TRUE))$cache$derive)
  expect_error(do.call(et_latent_data, lt_args(cache = "yes")), "must be TRUE, FALSE")
  expect_error(do.call(et_latent_data, lt_args(cache = list(shrd = "a"))), "must be TRUE, FALSE")
  expect_error(do.call(et_latent_data, lt_args(cache = list(shared = 1))), "character vector")
  expect_error(do.call(et_latent_data, lt_args(cache = list(latents = "S"))),
               "unknown latent")
})

# ---- the generated source -----------------------------------------------------

test_that("the new arguments reach latent_hmm_data", {
  src <- et_julia_source(lt_feature_model())$src
  expect_match(src, "observation_shared=et_observation_shared", fixed = TRUE)
  expect_match(src, "shared_key=et_shared_key,\n    n_keys=3", fixed = TRUE)
  expect_match(src, "per_evaluation=et_per_evaluation", fixed = TRUE)
  expect_match(src, "function et_observation_shared(model, data, i, key)", fixed = TRUE)
  expect_match(src, "function et_per_evaluation(model, data)", fixed = TRUE)
  expect_match(src, "function et_shared_key(t)", fixed = TRUE)
  # the observations' number type covers what the shared values are made of,
  # or a shared value carrying derivatives would be converted to a plain number
  expect_match(src, paste0("function et_observation_weights(model, data, i, t, shared)\n",
                           "    ETR_T = promote_type(Float64, eltype(model.p_obs))"), fixed = TRUE)
  # without a key, ET's default (one key per timepoint) is left alone
  src2 <- et_julia_source(lt_feature_model(key = FALSE))$src
  expect_false(grepl("shared_key=", src2, fixed = TRUE))
  expect_match(src2, "observation_shared=et_observation_shared", fixed = TRUE)
})

test_that("per_evaluation's numbers are settled to the parameter type, its indices left alone", {
  src <- et_julia_source(lt_feature_model())$src
  expect_match(src, paste0("(; a_step = _etr_settle(ETR_T, model.a + 0), ",
                           "n_two = _etr_settle(ETR_T, 2), p_obs = _etr_settle(ETR_T, w))"),
               fixed = TRUE)
  expect_match(src, "_etr_settle(::Type{T}, x::Union{Integer,Bool}) where {T} = x", fixed = TRUE)
  # numeric(n) follows the parameter type, so a dual number can be stored in it
  expect_match(src, "w = zeros(ETR_T, 2)", fixed = TRUE)
  # bodies reading per_evaluation's values floor their number type at Float64
  expect_match(src, "ETR_T = promote_type(Float64, eltype(model.p_obs))", fixed = TRUE)
})

test_that("observation_shared converts a number it returns, and each field of a list", {
  r1 <- role(function(model, data, i, key) {
    if (key == 1) return(0)
    model$x * 2
  }, "latent_obs_shared", "os")$src
  expect_match(r1, "return convert(ETR_T, 0)", fixed = TRUE)
  expect_match(r1, "convert(ETR_T, model.x * 2)", fixed = TRUE)
  r2 <- role(function(model, data, i, key) list(b = 1, a = model$x),
             "latent_obs_shared", "os")$src
  expect_match(r2, "(; a = convert(ETR_T, model.x), b = convert(ETR_T, 1))", fixed = TRUE)
  expect_error(role(function(model, data, i, t) 1, "latent_obs_shared", "os"),
               "must be declared as")
})

test_that("the cache's parameters are traced through per_evaluation, field by field", {
  g <- et_julia_source(lt_feature_model())
  # observation_shared reads p_obs, which per_evaluation computes from p alone;
  # a_step (from a) feeds only the weights, so it must not enter the cache keys.
  expect_identical(g$cache$shared, "p")
  expect_identical(g$cache$weights, character())
  expect_match(g$src, "cache=(; shared = (:p,), weights = ())", fixed = TRUE)
})

test_that("tracing follows locals, element assignments, loops and conditions", {
  tr <- EpidemicTrajectoriesR:::per_eval_field_reads
  f <- function(model, data) {
    g <- numeric(3)
    acc <- 0
    for (k in 1:3) {
      acc <- acc + model$raw[k]
      g[k] <- model$s * acc
    }
    z <- 0
    if (model$flag > 0) z <- model$c
    list(g = g, z = z, k0 = model$k)
  }
  out <- tr(f, list())
  expect_setequal(out$g, c("raw", "s"))
  expect_setequal(out$z, c("flag", "c"))
  expect_identical(out$k0, "k")
  # a read carried round a loop: b is assigned after it is read
  f2 <- function(model, data) {
    a <- 0; b <- 0
    for (k in 1:2) {
      a <- a + b
      b <- model$x
    }
    list(a = a)
  }
  expect_identical(tr(f2, list())$a, "x")
  # helpers count
  h <- et_helper(function(model) model$h, name = "hh")
  expect_identical(tr(function(model, data) list(v = hh(model)), list(h))$v, "h")
  # untraceable: the last value is not a list
  expect_null(tr(function(model, data) model$x, list()))
  expect_null(tr(et_julia("(; a = 1)"), list()))
})

test_that("a hand-written cache is checked against what the bodies read", {
  m_ok <- lt_feature_model(cache = list(shared = c("p", "a"), weights = "b"))
  expect_match(et_julia_source(m_ok)$src, "cache=(; shared = (:p, :a), weights = (:b,))",
               fixed = TRUE)
  expect_error(et_julia_source(lt_feature_model(cache = list(shared = "a"))),
               "shared needs p")
  expect_error(et_julia_source(lt_feature_model(cache = list(shared = "nope"))),
               "unknown parameter")
  m_lat <- lt_feature_model(cache = list(shared = "p", latents = "S"))
  expect_match(et_julia_source(m_lat)$src, "latents = (:S,)", fixed = TRUE)
})

test_that("cache = TRUE refuses what it cannot trace, and per_evaluation may not shadow a parameter", {
  fx <- lt_fixture()
  m <- lt_feature_model(fx)
  d <- m$data
  d$observation_shared <- et_julia("1.0")
  expect_error(et_julia_source(et_model(d, m$parameters)), "et_julia\\(\\) without `reads`")
  d <- m$data
  d$per_evaluation <- function(model, data) list(p = model$p)
  expect_error(et_julia_source(et_model(d, m$parameters)), "replace the parameter")
  d <- m$data
  d$observation_shared <- function(model, data, i, key) model$nothing_here
  expect_error(et_julia_source(et_model(d, m$parameters)), "model\\$nothing_here")
})

# ---- et_model(log_density = ) ------------------------------------------------

test_that("log_density terms are transpiled and added to the model body", {
  ld <- function(model) et_logpdf(normal_dist(0, 2), model$a) - et_logpdf(normal_dist(0, 1), model$b)
  m <- lt_feature_model(log_density = ld)
  src <- et_julia_source(m)$src
  expect_match(src, "function et_log_density_1(model)", fixed = TRUE)
  expect_match(src, paste0("Distributions.logpdf(Distributions.Normal(0, 2), model.a) - ",
                           "Distributions.logpdf(Distributions.Normal(0, 1), model.b)"), fixed = TRUE)
  expect_match(src, "pars = (; a=a, b=b, p=p, q=q)\n    @addlogprob! latent_loglik(pars, data, latents)\n    @addlogprob! et_log_density_1(pars)",
               fixed = TRUE)
  # on a sampled-trajectory model the term carries the depends= its reads give
  toy <- toy_model()
  m2 <- et_model(toy$data, toy$parameters, log_density = list(function(model) -model$alpha^2))
  expect_match(et_julia_source(m2)$src, "@addlogprob! et_log_density_1(pars)  depends=(:alpha,)",
               fixed = TRUE)
  m3 <- et_model(toy$data, toy$parameters,
                 log_density = et_julia("-model.alpha^2"))
  expect_match(et_julia_source(m3)$src, "@addlogprob! et_log_density_1(pars)\n", fixed = TRUE)
})

test_that("log_density is refused in the wrong shape, or reading an unknown name", {
  toy <- toy_model()
  expect_error(et_model(toy$data, toy$parameters, log_density = 3), "function\\(model\\)")
  expect_error(et_model(toy$data, toy$parameters, log_density = list(function(x) 1)),
               "term 1 must be a function\\(model\\)")
  m <- et_model(toy$data, toy$parameters, log_density = function(model) model$zeta)
  expect_error(et_julia_source(m), "model\\$zeta, which is not a parameter")
})

test_that("et_logpdf transpiles each family, and gives the same density in R", {
  expect_equal(jlq(et_logpdf(beta_dist(2, 3), x)),
               "Distributions.logpdf(Distributions.Beta(2, 3), x)")
  expect_equal(jlq(et_logpdf(exponential_dist(2), x)),
               "Distributions.logpdf(Distributions.Exponential(1 / (2)), x)")
  expect_equal(jlq(et_logpdf(normal_dist(log(20), 1), model$a)),
               "Distributions.logpdf(Distributions.Normal(log(20), 1), model.a)")
  expect_error(jlq(et_logpdf(x, 1)), "must be a call to one of")
  expect_error(jlq(et_logpdf(normal_dist(sd = 1), 1)), "positionally")
  expect_error(jlq(et_logpdf(normal_dist(0, 1))), "a distribution and a value")
  expect_equal(et_logpdf(normal_dist(1, 2), 0.3), dnorm(0.3, 1, 2, log = TRUE))
  expect_equal(et_logpdf(exponential_dist(2), 0.3), dexp(0.3, 2, log = TRUE))
  expect_equal(et_logpdf(gamma_dist(2, 3), 0.3), dgamma(0.3, shape = 2, scale = 3, log = TRUE))
})

test_that("NaN is a literal, NA is refused", {
  expect_identical(jlq(dnorm(1, 0, NaN)), "(NaN > 0 ? exp(-((1 - 0) / NaN)^2 / 2) / (NaN * sqrt(2pi)) : NaN)")
  expect_error(EpidemicTrajectoriesR:::julia_num(NA_real_), "cannot emit NA")
})

test_that("as.integer truncates toward zero, and is.finite / is.nan map across", {
  expect_equal(jlq(as.integer((x - x0) / w) + 1), "trunc(Int, (x - x0) / w) + 1")
  expect_equal(jlq(is.finite(v)), "isfinite(v)")
  expect_equal(jlq(is.nan(v)), "isnan(v)")
  expect_error(jlq(as.integer(x, y)), "one argument")
  expect_error(jlq(is.finite(x, y)), "one number")
})

# ---- per-element priors --------------------------------------------------------

test_that("a vector parameter may give each element its own prior", {
  p <- prior(list(beta_dist(2, 3), beta_dist(4, 5), beta_dist(1, 1)), init = c(0.4, 0.4, 0.5))
  expect_identical(p$n, 3L)
  expect_identical(EpidemicTrajectoriesR:::par_dist_src("Se", p),
                   "PracticalBayes.arraydist([Beta(2.0, 3.0), Beta(4.0, 5.0), Beta(1.0, 1.0)])")
  expect_error(prior(list(beta_dist(2, 3), beta_dist(4, 5)), init = c(0.4, 0.4), n = 3),
               "n = 3")
  expect_error(prior(list(beta_dist(2, 3), 5), init = c(0.4, 0.4)), "one distribution per element")
  expect_error(prior(list(beta_dist(2, 3)), init = 0.4), "vector parameter")
  expect_error(prior(list(beta_dist(2, 3), beta_dist(1, 1)), init = 0.4), "has length 1")
  expect_output(print(p), "arraydist")
})

test_that("a prior truncated on one side is emitted with an infinite other bound", {
  expect_identical(EpidemicTrajectoriesR:::dist_to_julia(truncated_dist(normal_dist(0, 0.05), lower = 0)),
                   "truncated(Normal(0.0, 0.05), 0.0, Inf)")
})

# ---- tier 2 ----------------------------------------------------------------------

test_that("per_evaluation, a keyed observation_shared and the cache leave the likelihood unchanged", {
  skip_without_latent_api()
  fx <- lt_fixture()
  pars <- list(a = 0.7, b = 0.3, p = 0.45, q = 0.2)
  ref <- et_marginal_loglik(lt_model(fx), pars)
  for (cache in list(TRUE, FALSE)) for (key in c(TRUE, FALSE)) {
    m <- lt_feature_model(fx, cache = cache, key = key)
    expect_equal(et_marginal_loglik(m, pars), ref, tolerance = 1e-12)
    expect_equal(et_marginal_loglik(m, pars, by_individual = TRUE),
                 lt_reference(pars, fx), tolerance = 1e-12)
  }
})

test_that("a cached model samples, with blocks that do and do not move the cache keys", {
  skip_without_latent_api()
  m <- lt_feature_model(log_density = function(model) -0.5 * model$a^2)
  fit <- et_sample(m, list(et_adaptive_hmc(c("a", "b"), n_steps = 4),
                           et_adaptive_hmc(c("p", "q"), n_steps = 4)),
                   n_sweeps = 20, n_adapts = 20, quiet = TRUE)
  expect_length(fit$draws$p, 20)
  expect_true(length(unique(fit$draws$a)) > 1 && length(unique(fit$draws$p)) > 1)
  expect_true(all(fit$draws$p > 0 & fit$draws$p < 1))
})
