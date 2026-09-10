#' EpidemicTrajectoriesR: epidemic models from R, fitted in Julia
#'
#' Specify a discrete-time individual-level epidemic model with R objects, write
#' its rate and observation functions in a restricted subset of R, and fit it
#' with the Julia stack (`EpidemicTrajectories.jl`, `PracticalBayes.jl`,
#' `PracticalEpiBayes.jl`) without writing Julia.
#'
#' @section The shape of a model:
#' \enumerate{
#'   \item [et_aggregate()] -- the arrays tracked during the latent update, with
#'     their reversible updates. The package attaches no meaning to them.
#'   \item [et_transitions()] -- the states, the transitions between them, and the
#'     rate function for each. Optionally an [et_survival()] making every step
#'     conditional on survival.
#'   \item [et_data()] -- problem size, group structure, the observation process,
#'     and your own arrays.
#'   \item [et_model()] -- parameters, priors, and deterministic
#'     reparameterisations.
#'   \item [et_sample()] -- Gibbs blocks and the run.
#' }
#'
#' @section Two things worth knowing up front:
#' The likelihood has **two halves** -- the transitions and the observations --
#' and both are needed. This package always emits both when an observation
#' process is given, and warns when one is not, because omitting the observation
#' term is a silent modelling bug: every observation parameter would then be
#' sampled from its prior with no data reaching it.
#'
#' The `depends=` annotations that let each Gibbs block skip likelihood terms it
#' cannot move are **derived** from your function bodies, not declared by you.
#' [et_check_depends()] verifies them.
#'
#' @section Escape hatches:
#' [et_julia()] accepts literal Julia anywhere a transpiled R function is
#' accepted; [et_kernel_julia()] drops a raw Gibbs kernel into the sampler;
#' [et_eval()] reaches anything in the Julia session; and [et_julia_source()]
#' prints the generated module, which runs standalone.
#'
#' @keywords internal
"_PACKAGE"
