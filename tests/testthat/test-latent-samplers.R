# A field's joint moves and particle Gibbs (et_field(moves =, pgas =)), and
# interweaving parameters with its paths (et_interweave()).

# lt_model() from helper-latent.R, with the field declared as given.
lt_model_field <- function(S, fx = lt_fixture()) {
  tr <- et_transitions(c("alive", "dead"), weights = TRUE, shared = lt_step,
                       "alive -> alive" = ~ shared$live,
                       "alive -> dead"  = ~ shared$die,
                       "dead -> dead"   = 1)
  latents <- et_latents(
    S   = S,
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

lt_anchors <- function() {
  a <- array(0.1, c(1L, LT_Y, LT_N))
  a[1, 3, 2] <- NA
  a
}

test_that("et_field's moves and et_pgas are checked on the way in", {
  expect_error(et_field(1, 3, prior = normal_dist(), index = lt_year, moves = "jump"),
               "must be among")
  expect_error(et_field(1, 3, prior = normal_dist(), touches = function(e) e:e,
                        pgas = et_pgas()), "needs `index`")
  expect_error(et_field(1, 3, prior = normal_dist(), index = lt_year, pgas = list(particles = 4)),
               "et_pgas()")
  expect_error(et_pgas(particles = 1), "at least 2")
  expect_error(et_pgas(scales = c(1, -1)), "positive")
  expect_error(et_pgas(anchor_scale = c(1, 2)), "one number")
  expect_error(et_field(1, 3, prior = normal_dist(), index = lt_year,
                        pgas = et_pgas(anchors = array(0, c(2L, 3L, LT_N)))), "1 x 3")
  p <- et_pgas(anchors = lt_anchors())
  expect_true(is.nan(p$anchors[1, 3, 2]))
})

test_that("moves and pgas reach the Field declaration, the anchors as payload", {
  m <- lt_model_field(et_field(1, LT_Y, prior = normal_dist(0, 1), index = lt_year,
                               init = lt_fixture()$S0, moves = c("pair", "neighbour"),
                               pgas = et_pgas(8, scales = c(0.1, 1), anchors = lt_anchors(),
                                              anchor_scale = 0.5)))
  src <- et_julia_source(m, list(et_adaptive_hmc(c("a", "b", "p", "q"))))$src
  expect_match(src, paste0("S ~ Field(1, 3; prior=Normal(0.0, 1.0), index=et_index_S, scale=1.0, ",
                           "moves=(:pair, :neighbour,), pgas=(; particles=8, scales=(0.1, 1.0,), ",
                           "anchors=_etr_S_anchors, anchor_scale=0.5))"), fixed = TRUE)
  expect_match(src, "const _etr_S_anchors =", fixed = TRUE)
  plain <- et_julia_source(lt_model(), list(et_adaptive_hmc(c("a", "b", "p", "q"))))$src
  expect_match(plain, "index=et_index_S, scale=1.0)", fixed = TRUE)
})

test_that("et_interweave gives the latent kernel its parameters", {
  m <- lt_model_field(et_field(1, LT_Y, prior = normal_dist(0, 1), index = lt_year,
                               init = lt_fixture()$S0))
  step <- function(model, u) {
    if (model$b + u > 3) return(list(b = 0, logjac = -Inf))
    list(b = model$b + u, logjac = 0)
  }
  g <- et_julia_source(m, list(et_interweave("b", field = "S", step = step,
                                             elements = function(i) 1:3),
                               et_adaptive_hmc(c("a", "p", "q"))))
  src <- g$src
  expect_match(src, paste0("(:latents, :b) => interweave_kernel(LATENT!, Interweave((:b,), ",
                           "et_interweave_step; field=:S, elements=et_interweave_elements, ",
                           "radius=Inf, scale=0.05, centred=0.05, n=3); n_adapt=N_ADAPT, params=et_params)"),
               fixed = TRUE)
  expect_match(src, "function et_interweave_step_r(model, u)", fixed = TRUE)
  expect_match(src, "function et_interweave_elements(i)", fixed = TRUE)
  expect_match(src, "et_sampler_for(DATA, LATENT!, N_ADAPT=0) = Gibbs(", fixed = TRUE)
  default <- et_julia_source(m, list(et_interweave("b", field = "S"),
                                     et_adaptive_hmc(c("a", "p", "q"))))$src
  expect_match(default, "et_interweave_step(v, u) = ((; b=v.b .+ u), 0.0)", fixed = TRUE)
  expect_error(et_julia_source(m, list(et_interweave("b", field = "u"),
                                       et_adaptive_hmc(c("a", "p", "q")))),
               "must name an et_field")
  expect_error(et_julia_source(m, list(et_interweave("b", field = "S"),
                                       et_adaptive_hmc(c("a", "b", "p", "q")))),
               "more than one block")
  expect_error(et_interweave("b", field = "S", step = function(model) 1), "function\\(model, u\\)")
  expect_error(et_interweave("b", field = "S", radius = 0), "positive")
  radial <- et_julia_source(m, list(et_interweave("b", field = "S", radius = 2),
                                    et_adaptive_hmc(c("a", "p", "q"))))$src
  expect_match(radial, "elements=nothing, radius=2.0, scale=0.05", fixed = TRUE)
})

test_that("a model samples with joint moves, pgas and interweaving", {
  skip_without_latent_api()
  skip_if_not(isTRUE(et_eval("isdefined(PracticalEpiBayes, :interweave_kernel)")),
              "PracticalEpiBayes has no interweave_kernel")
  m <- lt_model_field(et_field(1, LT_Y, prior = normal_dist(0, 1), index = lt_year,
                               init = lt_fixture()$S0, moves = c("pair", "shift"),
                               pgas = et_pgas(8, scales = 0.5)))
  fit <- et_sample(m, list(et_interweave("b", field = "S", radius = 0.5),
                           et_adaptive_hmc(c("a", "p", "q"), n_steps = 4)),
                   n_sweeps = 30, n_adapts = 20, quiet = TRUE)
  expect_length(fit$draws$b, 30)
  expect_true(all(is.finite(fit$draws$b)))
  expect_true(length(unique(fit$draws$b)) > 1 && length(unique(fit$draws$a)) > 1)
})
