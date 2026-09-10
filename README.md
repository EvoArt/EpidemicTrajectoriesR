# EpidemicTrajectoriesR

Fit discrete-time individual-level epidemic models from R.

You declare the model — its states, the transitions between them, what gets
tracked, how the observations relate to the hidden states — with R objects, and
write the rate and observation functions in ordinary R. Both are translated to
Julia and run against `EpidemicTrajectories.jl`, `PracticalBayes.jl` and
`PracticalEpiBayes.jl`, so the hidden trajectory is resampled by iFFBS once per
Gibbs sweep rather than inside every gradient evaluation.

No Julia is required to use it, and none appears in a model specification.

## Installation

```r
# install.packages("remotes")
remotes::install_github("EvoArt/EpidemicTrajectoriesR")

library(EpidemicTrajectoriesR)
et_setup()     # starts Julia and resolves the pinned project on first use
```

`et_setup()` will install Julia's dependencies into a project private to this
package. The first call takes a few minutes; later ones are immediate.

## A complete model

The two-state S/I cattle model of Touloupou et al. (2019): animals are infected
at a rate driven by how many of their penmates are infected, recover back to
susceptible, and are never observed directly — only two imperfect diagnostic
tests, on a subset of days.

```r
states <- c("S", "I")

# What to track during the latent update, and how it updates. The package
# attaches no meaning to this array; it runs the update forwards and, when the
# sampler needs to, in reverse.
aggs <- et_aggregate(
  states,
  arrays = list(n_infected = et_array("Int", c(n_pens, n_days))),
  update = function(model, data, X, state, i, t) {
    n_infected[data$group[i], t] <- n_infected[data$group[i], t] + (state == "I")
  })

# The rates read that array. During an individual's own update it holds the
# leave-one-out count, which is exactly what a force of infection needs.
infection <- function(model, data, i, t) {
  I_minus <- data$aggregates$n_infected[data$group[i], t]
  -expm1(-(model$alpha + model$beta * I_minus))
}
recovery <- function(model, data, i, t) 1 / model$m

trans <- et_transitions(states, "S -> I" = infection, "I -> S" = recovery)

dat <- et_data(
  n_individuals = n_animals, n_timepoints = n_days,
  transitions = trans, starting_state = starting_state, aggregates = aggs,
  group = pen, observation_process = cattle_observations,
  extras = list(rams = rams, faecal = faecal))   # your own data, as data$rams

mod <- et_model(
  data = dat,
  parameters = list(
    alpha   = et_par(et_gamma(1, 1), init = 0.05),
    beta    = et_par(et_gamma(1, 1), init = 0.05),
    m_tilde = et_par(et_gamma(2, 4), init = 4.0),
    nu      = et_par(et_beta(1, 1),  init = 0.10),
    theta_r = et_par(et_beta(1, 1),  init = 0.70),
    theta_f = et_par(et_beta(1, 1),  init = 0.40)),
  derived = list(m = quote(m_tilde + 1)))

fit <- et_sample(mod,
  blocks = list(
    et_nuts(c("alpha", "beta", "m_tilde")),
    et_conjugate_initial_state("nu", at = 1, states = c("S", "I")),
    et_conjugate_test_sensitivity("theta_r", y = "rams",   infected_state = "I"),
    et_conjugate_test_sensitivity("theta_f", y = "faecal", infected_state = "I")),
  n_sweeps = 2000, n_burn = 800, n_adapts = 400, seed = 7)

summary(fit)
```

`blocks` is optional. Omit it and every parameter goes into one NUTS block, with
the trajectory sampled by iFFBS — a valid fit with no blocking decisions at all.

## Seeing the Julia

The generated module is a first-class artefact, not a hidden intermediate. It is
commented, it runs standalone, and it is what to show a Julia-literate colleague
when something goes wrong.

```r
et_julia_source(mod)      # prints the module
et_source(fit)            # the module behind a fit
et_eval("...")            # run anything in the session
```

## Checking a model before you trust it

```r
et_loglik(mod, blocks)              # log density at the initial values
et_iffbs_sweep(mod, X0, blocks)     # one latent sweep, in isolation
et_check_depends(mod, blocks)       # verify the derived depends= annotations
```

`et_check_depends()` is worth knowing about. Each Gibbs block skips the
likelihood terms its parameters cannot move, which is a large speed-up and a
silent hazard: getting the dependency list wrong drops a real gradient
contribution with no error. This package derives the lists from your function
bodies rather than asking you for them, and `et_check_depends()` perturbs each
parameter to confirm no term it was excluded from actually moves.

## Autodiff backends

```r
et_sample(mod, adtype = "polyester")                 # threaded forward mode
et_setup(ad_backends = "mooncake")                   # resolve a reverse-mode one
et_sample(mod, adtype = "mooncake")
```

Forward mode costs roughly one pass per *chunk* of parameters, so it rewards
small blocks; reverse mode costs one pass whatever the block size. The choice
changes the time taken, never the answer.

## Writing the functions

Rate and observation functions are ordinary R, restricted to what has an exact
Julia equivalent: arithmetic and comparison, `if`/`else`, `for` over a range,
indexing, `$`, `ifelse`, and the usual maths functions. Anything outside that
raises an error naming the construct, rather than emitting something plausible
that means a different thing in Julia.

`et_julia()` takes literal Julia anywhere a function is accepted, and
`et_helper()` declares your own R functions for the bodies to call.

## Examples

- `inst/examples/cattle.R` — the model above, end to end, recovering the truth
  from simulated data.
- `inst/examples/badger.R` — a badger bovine-TB SEID model: 2384 individuals x
  161 quarters, six diagnostic tests, capture-recapture observation,
  time-varying social groups and Gompertz-Makeham survival.

## Licence

MIT.
