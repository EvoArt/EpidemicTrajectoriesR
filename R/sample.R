#' Fit the model.
#'
#' @param model An [et_model()].
#' @param blocks List of sampler blocks ([et_nuts()], [et_hmc()],
#'   [et_conjugate()] and friends). Any parameter not named in a block is
#'   collected into one NUTS block, and the trajectory always gets an iFFBS
#'   kernel, so omitting `blocks` entirely still gives a valid sampler.
#' @param n_sweeps Kept Gibbs sweeps.
#' @param n_burn Sweeps discarded before collection.
#' @param n_adapts Adaptation steps passed to every HMC/NUTS block.
#' @param seed RNG seed.
#' @param x_init Optional starting trajectory, an `n_timepoints x
#'   n_individuals` integer matrix of 1-based state codes. Defaults to
#'   all-first-state.
#' @param adtype The automatic-differentiation backend: a backend name or an
#'   [et_adtype()]. Forward mode (the default) costs roughly one pass per CHUNK
#'   of parameters, so it rewards small blocks; reverse mode costs one pass
#'   whatever the block size. It changes only the time taken, never the answer.
#' @param save_x Optional path (e.g. `"X.rds"`) to stream the latent trajectory
#'   to during sampling. Required for [et_residuals()]. The trajectory is always
#'   kept live for conditioning; this only decides where the per-sweep copy goes,
#'   and by default it is dropped rather than carried in the chain: it is
#'   `n_timepoints x n_individuals` integers every sweep.
#'
#'   The archive is written as `.rds`, which R reads natively. Sampling streams
#'   through Julia's own format and [et_convert_trajectories()] rewrites it once
#'   sampling finishes; pass a `.jld2` path to skip that and keep the Julia
#'   files. Either way the archive is one file per flush, not one huge file.
#' @param convert_x Convert the archive to `.rds` after sampling. Only consulted
#'   when `save_x` names a `.rds` path. Set `FALSE` to defer it -- the fit still
#'   records where the archive is, so [et_convert_trajectories()] can be run
#'   later.
#' @param save_every Flush the archive every this many sweeps.
#' @param keep_source Optional path to write the generated Julia to.
#' @param quiet Suppress progress messages.
#' @return An object of class `et_fit`.
#' @export
et_sample <- function(model, blocks = list(), n_sweeps = 1000, n_burn = 0,
                      n_adapts = 0, seed = 1, x_init = NULL,
                      adtype = "forwarddiff", save_x = NULL, save_every = 100,
                      convert_x = TRUE, keep_source = NULL, quiet = FALSE) {
  et_require_session()
  gen <- et_julia_source(model, blocks)
  mod <- et_load_module(gen, keep_source)

  x_arg <- "nothing"
  if (!is.null(x_init)) {
    x_init <- as.matrix(x_init)
    storage.mode(x_init) <- "integer"
    check_x_init(x_init, model)
    nm <- paste0(mod, "_xinit")
    JuliaCall::julia_assign(nm, x_init)
    x_arg <- sprintf("Matrix{Int}(Main.%s)", nm)
  }
  if (!quiet) {
    message(sprintf("Running %d sweeps (burn %d, adapt %d) ...",
                    n_sweeps, n_burn, n_adapts))
  }
  # The archive is addressed by template, and the two formats differ only in the
  # extension. Sampling always streams through Julia's own backend; asking for
  # `.rds` means "convert it afterwards", not "have Julia write R files": no
  # Julia package can write `.rds`, and conversion is per flush, so it costs
  # nothing against the sweeps it follows.
  want_rds <- FALSE
  if (!is.null(save_x)) {
    save_x <- normalizePath(save_x, winslash = "/", mustWork = FALSE)
    want_rds <- identical(tolower(tools::file_ext(save_x)), "rds")
    julia_save_x <- if (want_rds) et_with_ext(save_x, "jld2") else save_x
    # `read_states` dispatches through the backend that wrote the files.
    JuliaCall::julia_command("using JLD2")
  } else {
    julia_save_x <- NULL
  }
  call <- sprintf(
    "%s.et_run(; n_sweeps=%s, n_burn=%s, n_adapts=%s, seed=%s, x_init=%s, adtype=%s, save_x=%s, save_every=%s)",
    mod, julia_int(n_sweeps), julia_int(n_burn), julia_int(n_adapts),
    julia_int(seed), x_arg, adtype_to_julia(adtype),
    if (is.null(julia_save_x)) "nothing" else julia_string(julia_path(julia_save_x)),
    julia_int(save_every))
  raw <- JuliaCall::julia_eval(call)

  elapsed <- raw[["_elapsed"]]
  raw[["_elapsed"]] <- NULL
  draws <- restore_shapes(raw, model)
  if (!quiet) {
    message(sprintf("done in %.1f s (%.3f s/sweep)", elapsed,
                    elapsed / max(1, n_sweeps + n_burn)))
  }
  # Convert after sampling, so a crash during the fit still leaves every
  # already-flushed chunk on disc for et_convert_trajectories() to recover.
  if (want_rds && isTRUE(convert_x)) {
    et_convert_trajectories(save_x, quiet = quiet)
  }
  structure(list(draws = draws, model = model, blocks = gen$blocks,
                 depends = gen$depends, source = gen, module = mod,
                 x_path = save_x,
                 elapsed = elapsed, n_sweeps = n_sweeps, n_burn = n_burn,
                 seed = seed), class = "et_fit")
}

check_x_init <- function(x_init, model, caller = "et_sample()") {
  d <- model$data
  if (!identical(dim(x_init), c(d$n_timepoints, d$n_individuals))) {
    stop(caller, ": `x_init` must be ", d$n_timepoints, " x ",
         d$n_individuals, " (timepoints x individuals), got ",
         paste(dim(x_init), collapse = " x "),
         ". ET indexes the trajectory X[t, i].", call. = FALSE)
  }
  ns <- length(d$transitions$states)
  bad <- x_init[x_init < 1 | x_init > ns]
  if (length(bad)) {
    stop(caller, ": `x_init` holds state code(s) outside 1..", ns, ": ",
         paste(utils::head(unique(bad), 5), collapse = ", "),
         ". Codes are 1-based positions in the state space (",
         paste(d$transitions$states, collapse = ", "), ").", call. = FALSE)
  }
  invisible(TRUE)
}

# Julia hands back a vector per scalar parameter and an n_draws x prod(dim)
# matrix otherwise; restore the declared shape so `fit$draws$nu` is an array of
# the right dimensions rather than a flat block.
#
# Only sampled parameters come back from the chain: PracticalBayes' `:=`
# deterministics are not stored under their own key, so the derived values are
# recomputed here from the draws. That is exact: a deterministic is a function
# of the sampled values and nothing else, and it keeps `fit$draws$m` available
# without asking the sampler to carry a redundant column.
restore_shapes <- function(raw, model) {
  out <- list()
  for (nm in names(raw)) {
    v <- raw[[nm]]
    p <- model$parameters[[nm]]
    if (!is.null(p) && !is.null(p$dim) && is.matrix(v)) {
      out[[nm]] <- array(v, dim = c(nrow(v), p$dim))
    } else {
      out[[nm]] <- v
    }
  }
  out <- out[model$par_names]
  for (nm in names(model$derived)) {
    val <- tryCatch(eval(model$derived[[nm]], envir = out, enclos = baseenv()),
                    error = function(e) e)
    if (inherits(val, "error")) {
      warning("could not recompute derived '", nm, "' from the draws (",
              conditionMessage(val), "); it is omitted from fit$draws.",
              call. = FALSE)
      next
    }
    out[[nm]] <- val
  }
  out
}

#' The log density of a model at its initial parameters.
#'
#' A cheap structural check that needs no sampling: it loads the generated
#' module, establishes the aggregates-agree-with-X invariant, and evaluates both
#' likelihood terms. A non-finite value means the model is misspecified before
#' any sampler is involved.
#'
#' @param model An [et_model()].
#' @param blocks The same sampler blocks you will pass to [et_sample()]. They do
#'   not affect the log density, but passing them makes this call generate the
#'   same Julia module the fit will use, so nothing is compiled twice, and a
#'   model with a conjugate-owned parameter needs its kernel declared to build at
#'   all.
#' @param x_init Optional starting trajectory.
#' @return A named numeric vector with an `epidemic` and (if any) `observation`
#'   entry.
#' @export
et_loglik <- function(model, blocks = list(), x_init = NULL) {
  et_require_session()
  gen <- et_julia_source(model, blocks)
  mod <- et_load_module(gen)
  arg <- ""
  if (!is.null(x_init)) {
    x_init <- as.matrix(x_init); storage.mode(x_init) <- "integer"
    check_x_init(x_init, model)
    nm <- paste0(mod, "_llx")
    JuliaCall::julia_assign(nm, x_init)
    arg <- sprintf("Matrix{Int}(Main.%s)", nm)
  }
  res <- JuliaCall::julia_eval(sprintf("%s.et_loglik_at(%s)", mod, arg))
  unlist(res)
}

#' Simulate a trajectory from the model.
#'
#' Draws a state trajectory forward in time using the transitions the model
#' declares, so simulated data and fitted model cannot disagree about the
#' dynamics. Writing the simulator by hand instead means maintaining a second
#' copy of the same rates, and nothing checks that the two still match.
#'
#' Note this simulates the LATENT STATES only. Turning states into observations
#' is the user's own step, because only they know what was observed and how --
#' see the example below.
#'
#' @param model An [et_model()].
#' @param pars Named list of parameter values to simulate at. Defaults to each
#'   parameter's `init`. Derived values are recomputed, not supplied.
#' @param seed RNG seed.
#' @param blocks Sampler blocks, as for [et_loglik()]. They do not affect the
#'   simulation, but passing the same ones reuses the compiled module.
#' @return An `n_timepoints x n_individuals` integer matrix of 1-based state
#'   codes.
#' @export
#' @examples
#' \dontrun{
#' X <- et_simulate(model, pars = list(alpha = 0.005, beta = 0.03))
#'
#' # Observe it: caught with probability p while alive, then tested.
#' caught <- (X != 3L) & matrix(runif(length(X)) < 0.6, nrow(X), ncol(X))
#' }
et_simulate <- function(model, pars = NULL, seed = 1, blocks = list()) {
  et_require_session()
  if (!inherits(model, "et_model")) {
    stop("et_simulate(): `model` must come from et_model().", call. = FALSE)
  }
  gen <- et_julia_source(model, blocks)
  mod <- et_load_module(gen)

  vals <- lapply(model$par_names, function(nm) model$parameters[[nm]]$init)
  names(vals) <- model$par_names
  if (!is.null(pars)) {
    unknown <- setdiff(names(pars), model$par_names)
    if (length(unknown)) {
      stop("et_simulate(): unknown parameter(s): ",
           paste(unknown, collapse = ", "), ". The model declares: ",
           paste(model$par_names, collapse = ", "), ".", call. = FALSE)
    }
    for (nm in names(pars)) {
      if (length(pars[[nm]]) != length(vals[[nm]])) {
        stop("et_simulate(): '", nm, "' has length ", length(pars[[nm]]),
             " but the model declares ", length(vals[[nm]]), ".", call. = FALSE)
      }
      vals[[nm]] <- as.numeric(pars[[nm]])
    }
  }

  et_simulate_inject(mod)
  nt <- paste(sprintf("%s=%s", names(vals),
                      vapply(vals, function(v)
                        if (length(v) > 1) julia_vector(v, "Float64")
                        else julia_float(v), character(1))),
              collapse = ", ")
  out <- JuliaCall::julia_eval(sprintf(
    "%s.et_simulate_at((; %s); seed=%s)", mod, nt, julia_int(seed)))
  matrix(as.integer(out), nrow = model$data$n_timepoints)
}

# The simulator entry point, injected once per module. Not part of
# et_julia_source(), whose text is hashed to fix the module name.
et_simulate_inject <- function(mod) {
  src <- paste0(
    "if !isdefined(@__MODULE__, :et_simulate_at)\n",
    "function et_simulate_at(v; seed=1)\n",
    "    epidemic_simulator(DATA)(StableRNG(seed), et_params(v))\n",
    "end\n",
    "end\n")
  nm <- paste0(mod, "_sim_src")
  JuliaCall::julia_assign(nm, src)
  JuliaCall::julia_command(sprintf("Base.include_string(%s, Main.%s)", mod, nm))
}

#' Run one iFFBS latent sweep in isolation.
#'
#' The latent sampler on its own, at the model's initial parameters. Useful for
#' checking that a trajectory moves at all, and for timing the latent half of a
#' sweep separately from the gradient.
#'
#' @param model An [et_model()].
#' @param x_init Starting trajectory (`n_timepoints x n_individuals`).
#' @param blocks The sampler blocks, as for [et_loglik()].
#' @param seed RNG seed.
#' @return The updated trajectory, an integer matrix.
#' @export
et_iffbs_sweep <- function(model, x_init, blocks = list(), seed = 1) {
  et_require_session()
  gen <- et_julia_source(model, blocks)
  mod <- et_load_module(gen)
  x_init <- as.matrix(x_init); storage.mode(x_init) <- "integer"
  check_x_init(x_init, model)
  nm <- paste0(mod, "_sweepx")
  JuliaCall::julia_assign(nm, x_init)
  out <- JuliaCall::julia_eval(sprintf(
    "%s.et_iffbs_sweep(Matrix{Int}(Main.%s); seed=%s)", mod, nm, julia_int(seed)))
  matrix(as.integer(out), nrow = nrow(x_init))
}
