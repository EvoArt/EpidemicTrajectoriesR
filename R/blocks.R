# Gibbs blocks. Each constructor returns a plain description; codegen.R turns
# it into one entry of the Julia `Gibbs(...)`.

new_block <- function(kind, vars, ...) {
  structure(c(list(kind = kind, vars = vars), list(...)), class = "et_block")
}

#' @export
print.et_block <- function(x, ...) {
  cat("<et_block> ", x$kind, ": ", paste(x$vars, collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' A NUTS block.
#' @param vars Character vector of parameter names sampled together.
#' @param target_accept Target acceptance rate.
#' @return An `et_block`.
#' @export
et_nuts <- function(vars, target_accept = 0.8) {
  new_block("nuts", as.character(vars), target_accept = target_accept)
}

#' A fixed-trajectory HMC block.
#'
#' The badger model uses this rather than NUTS, with a hand-tuned diagonal
#' metric. `n_steps = 15` matches the C++ reference's EXPECTED trajectory length
#' (it draws L uniformly on 1..30, mean 15.5), not its nominal `L = 30`.
#'
#' @param vars Character vector of parameter names sampled together.
#' @param n_steps Leapfrog steps per trajectory.
#' @param step_size Either a single number, or a named list giving a step size
#'   per parameter (scalars are recycled across a vector parameter's elements).
#'   Becomes the diagonal metric `DiagEuclideanMetric(step_size^2)`.
#' @return An `et_block`.
#' @export
et_hmc <- function(vars, n_steps = 15, step_size = NULL) {
  new_block("hmc", as.character(vars), n_steps = as.integer(n_steps),
            step_size = step_size)
}

#' A fixed-trajectory HMC block that adapts its metric and step size.
#'
#' The usual choice alongside an iFFBS trajectory kernel, and the one to reach
#' for first.
#'
#' Like [et_hmc()] it costs a constant `n_steps` leapfrog steps per iteration,
#' so the run has a flat predictable cost and no tree doubling. Unlike
#' [et_hmc()] it does not ask you for a step size: it adapts the mass matrix and
#' step size during warm-up exactly as NUTS does, then holds both fixed for
#' sampling.
#'
#' That difference matters more than it sounds. A hand-tuned metric has to be
#' right in both order and scale, and a scale that is far too small produces a
#' chain that never moves while every acceptance diagnostic still looks healthy.
#' Adapting removes that failure mode.
#'
#' @param vars Character vector of parameter names sampled together.
#' @param target_accept Target acceptance rate for step-size dual averaging.
#' @param n_steps Leapfrog steps per iteration, in warm-up and sampling alike.
#' @param metric `"diagonal"` (the default), `"dense"` or `"unit"`. A dense
#'   metric learns correlations between parameters, at the cost of estimating
#'   many more entries during warm-up, worth it when the posterior is strongly
#'   correlated and you can afford a long warm-up.
#' @return An `et_block`.
#' @seealso [et_hmc()] for a fixed metric you supply yourself, [et_nuts()] for
#'   tree-doubling with an adaptive trajectory length.
#' @export
et_adaptive_hmc <- function(vars, target_accept = 0.8, n_steps = 15,
                            metric = c("diagonal", "dense", "unit")) {
  metric <- match.arg(metric)
  new_block("adaptive_hmc", as.character(vars),
            target_accept = target_accept, n_steps = as.integer(n_steps),
            metric = metric)
}

#' The latent-trajectory block, resampled by iFFBS.
#'
#' Always present: if you omit it, one is added for you. The trajectory is stored
#' as one whole-matrix latent, held constant during the gradients and updated
#' once per Gibbs sweep, which is the entire reason the Julia stack exists.
#'
#' @param name The trajectory block's name. Leave as `"X"`.
#' @param mh Use the Metropolis-corrected iFFBS proposal
#'   (`epidemic_latent_sampler(mh = TRUE)`) instead of the plain Gibbs sweep.
#' @return An `et_block`.
#' @export
et_iffbs <- function(name = "X", mh = FALSE) {
  new_block("iffbs", name, mh = isTRUE(mh))
}

#' A conjugate Gibbs block with a user-written count body.
#'
#' The closed form stays open-ended: you write the counting loop in the R subset,
#' with the current trajectory `X` (indexed `X[t, i]`) and `data` in scope, and
#' return the counts the family expects.
#'
#' @param name Parameter name.
#' @param family `"beta"` (return `c(successes, failures)`) or `"dirichlet"`
#'   (return `c(count_1, ..., count_k)`).
#' @param prior `c(a, b)` for Beta, or the concentration vector for Dirichlet.
#' @param count A function `(X, data)`, or `(X, data, k)` when `n > 1`, whose last
#'   expression is the `c(...)` of counts.
#' @param n Number of independent draws collected into a vector.
#' @return An `et_block`.
#' @export
et_conjugate <- function(name, family = c("beta", "dirichlet"), prior, count,
                         n = 1L) {
  family <- match.arg(family)
  if (!is.function(count)) {
    stop("et_conjugate(): `count` must be a function.", call. = FALSE)
  }
  expected <- if (n > 1L) c("X", "data", "k") else c("X", "data")
  check_protocol_formals(count, expected, "conjugate count", name)
  if (family == "beta" && length(prior) != 2L) {
    stop("et_conjugate(): a Beta prior is c(a, b).", call. = FALSE)
  }
  new_block("conjugate", name, family = family, prior = as.numeric(prior),
            count = count, n = as.integer(n))
}

#' Conjugate Beta kernel for a diagnostic test's sensitivity.
#'
#' Wraps `PracticalEpiBayes::test_sensitivity_kernel`. Assumes perfect
#' specificity, susceptibles never test positive, so only infected tested cells
#' contribute.
#'
#' @param name Parameter name.
#' @param y Name of the result matrix in `extras` (`1` positive, `0` negative,
#'   negative = not tested), indexed `[t, i]`.
#' @param infected_state Name of the state counted as infected.
#' @param prior `c(a, b)`.
#' @return An `et_block`.
#' @export
et_conjugate_test_sensitivity <- function(name, y, infected_state,
                                          prior = c(1, 1)) {
  new_block("test_sensitivity", name, y = y, infected_state = infected_state,
            prior = as.numeric(prior))
}

#' Conjugate Beta kernel for a per-index capture probability.
#'
#' Wraps `PracticalEpiBayes::capture_prob_kernel` (the badger `etas`).
#'
#' @param name Parameter name.
#' @param caught,effort,group,index Names of entries in `extras`: capture history
#'   `[t, i]`, trapping effort `[g, t]`, group membership `[i, t]`, and the
#'   per-time index selecting which probability applies.
#' @param dead_state Name of the dead state.
#' @param n Number of probabilities (e.g. seasons).
#' @param prior `c(a, b)`.
#' @return An `et_block`.
#' @export
et_conjugate_capture_prob <- function(name, caught, effort, group, index,
                                      dead_state, n, prior = c(1, 1)) {
  new_block("capture_prob", name, caught = caught, effort = effort,
            group = group, index = index, dead_state = dead_state,
            n = as.integer(n), prior = as.numeric(prior))
}

#' Conjugate Dirichlet kernel for the initial-state mixing of entrants.
#'
#' Wraps `PracticalEpiBayes::initial_state_kernel` (the cattle `nu`, the badger
#' `nu`).
#'
#' @param name Parameter name.
#' @param at A single entry time, the name of an `extras` vector of entry times,
#'   or a numeric vector of them.
#' @param states Character vector of state names forming the Dirichlet
#'   categories, ordered to match `prior`.
#' @param prior Concentration vector.
#' @param eligible Optional predicate `(X, data, i, t)` selecting contributing
#'   individuals; defaults to all.
#' @param n Number of entry cohorts.
#' @return An `et_block`.
#' @export
et_conjugate_initial_state <- function(name, at = 1, states, prior = NULL,
                                       eligible = NULL, n = 1L) {
  if (is.null(prior)) prior <- rep(1, length(states))
  if (length(prior) != length(states)) {
    stop("et_conjugate_initial_state(): `prior` must have one entry per state ",
         "(got ", length(prior), " for ", length(states), " states).",
         call. = FALSE)
  }
  new_block("initial_state", name, at = at, states = states,
            prior = as.numeric(prior), eligible = eligible, n = as.integer(n))
}

#' Drop in a raw Julia Gibbs kernel.
#'
#' The honest answer for a bespoke MH kernel the R DSL cannot express: the
#' badger `xi` changepoint sampler is the worked example. `src` is a Julia
#' expression evaluated inside the generated module, where the transpiled
#' functions, `DATA` and the extras are all in scope.
#'
#' @param vars Character vector of the block's parameter names.
#' @param src Julia source for the kernel expression.
#' @return An `et_block`.
#' @export
et_kernel_julia <- function(vars, src) {
  stopifnot(is.character(src), length(src) == 1L)
  new_block("julia", as.character(vars), src = src)
}

# ---- assembly ---------------------------------------------------------------

# Fill in the blocks the user did not specify, and check the ones they did.
#
# Defaulting is what makes "minimal Julia knowledge" achievable: a user who
# passes no blocks at all still gets a valid sampler -- one NUTS block over every
# sampled parameter, plus the iFFBS trajectory kernel.
resolve_blocks <- function(model, blocks) {
  blocks <- blocks %||% list()
  if (!is.list(blocks)) stop("`blocks` must be a list.", call. = FALSE)
  for (b in blocks) {
    if (!inherits(b, "et_block")) {
      stop("every entry of `blocks` must come from et_nuts(), et_hmc(), ",
           "et_iffbs(), et_conjugate*() or et_kernel_julia().", call. = FALSE)
    }
  }
  is_traj <- vapply(blocks, function(b) b$kind == "iffbs", logical(1))
  traj <- if (any(is_traj)) blocks[[which(is_traj)[1]]] else et_iffbs()
  if (sum(is_traj) > 1L) {
    stop("more than one iFFBS block was given; there is one trajectory.",
         call. = FALSE)
  }
  rest <- blocks[!is_traj]

  claimed <- unlist(lapply(rest, function(b) b$vars), use.names = FALSE)
  dup <- claimed[duplicated(claimed)]
  if (length(dup)) {
    stop("parameter(s) appear in more than one block: ",
         paste(unique(dup), collapse = ", "),
         ". Each parameter belongs to exactly one Gibbs block.", call. = FALSE)
  }
  unknown <- missing_from(claimed, model$par_names)
  if (length(unknown)) {
    stop("block(s) name unknown parameter(s): ", paste(unknown, collapse = ", "),
         ". Declared parameters are: ",
         paste(model$par_names, collapse = ", "), ".", call. = FALSE)
  }
  leftover <- setdiff(model$par_names, claimed)
  # A conjugate-owned ("latent") parameter with no kernel is a real error: its
  # placeholder density is constant, so an HMC block would sample it from
  # nothing at all.
  orphan_latent <- leftover[vapply(leftover,
    function(nm) model$parameters[[nm]]$kind == "latent", logical(1))]
  if (length(orphan_latent)) {
    stop("parameter(s) declared with kind = \"latent\" have no conjugate ",
         "kernel: ", paste(orphan_latent, collapse = ", "),
         ". A latent parameter's density is a constant placeholder, so without ",
         "a kernel it would never be informed by anything.
",
         "  If you are calling et_loglik() / et_check_depends() / ",
         "et_iffbs_sweep(), pass the same `blocks` you will give et_sample().",
         call. = FALSE)
  }
  if (length(leftover)) {
    rest <- c(rest, list(et_nuts(leftover)))
  }
  c(rest, list(traj))
}
