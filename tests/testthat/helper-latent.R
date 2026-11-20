# The small weighted-chain model (et_latent_data) shared by test-latent.R and
# test-latent-features.R, and the gate for ET's latent API.

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

skip_without_latent_api <- function() {
  skip_without_julia()
  ok <- isTRUE(et_eval("isdefined(EpidemicTrajectories, :latent_hmm_data)"))
  testthat::skip_if_not(ok, "EpidemicTrajectories has no latent_hmm_data (not the latent-hmm branch)")
}
