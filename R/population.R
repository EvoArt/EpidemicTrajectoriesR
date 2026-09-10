# Layer 3: the two diagnostics that do not fit the per-individual archive route.
#
# `et_residuals()` scores an ARCHIVE: sample with `save_x =`, read the
# trajectories back, score them. That is the right default, but it fails at both
# ends of the size range for opposite reasons, and each end has a Julia-side
# answer this file exposes:
#
#   * too big to archive -- `et_collect()` wraps `SummaryCollector`, scoring each
#     draw as it is produced and keeping only the residuals. On the badger model
#     that is ~9.5 MB per summary against ~15 GB of trajectories.
#
#   * not computable one individual at a time -- `et_reproduction_numbers()`
#     wraps `case_reproduction_numbers`. R_i needs every infective present in the
#     group at the moment of an infection, so it takes a whole-population pass
#     and cannot use the `(model, data, X, i, rng)` contract the residuals do.
#
# Both take the SAME declared force-of-infection decomposition idea the rest of
# the package uses: the user names the components as ordinary rate functions and
# says which one is transmission. The package sums and normalises them and never
# needs to know what any of them means.

#' Score residuals during the fit, without ever storing a trajectory.
#'
#' The bounded-memory alternative to `et_sample(save_x =)` + [et_residuals()].
#' Each draw is scored as it is produced and only the residual values are kept,
#' so nothing about the trajectory is retained.
#'
#' Use this when the archive would be too large to keep: the trajectory is
#' `n_timepoints x n_individuals` integers EVERY sweep, and `save_every` is a
#' flush interval rather than thinning, so a long run over a large population
#' archives everything. Storage here is `n_individuals x n_keep` doubles per
#' residual instead -- on the badger model ~9.5 MB per residual against ~15 GB
#' of trajectories.
#'
#' The trade is that the residuals must be decided BEFORE the fit. Archiving
#' lets you add new ones later without refitting; this does not.
#'
#' @param model An [et_model()].
#' @param blocks Sampler blocks, as for [et_sample()].
#' @param residuals List of [et_residual_survival()] / [et_residual_waiting()].
#' @param n_sweeps,n_burn,n_adapts,seed,x_init,adtype As [et_sample()].
#' @param n_keep How many scored draws to store. 500 is far more than a
#'   uniformity test needs; the collector stops storing once it has this many.
#' @param thin Score every `thin`-th sweep. Defaults to spreading `n_keep` draws
#'   over the whole run.
#' @param quiet Suppress progress messages.
#' @return A list with `fit` (an [et_fit()], without an archive) and `residuals`
#'   (the same data frame [et_residuals()] returns).
#' @export
et_collect <- function(model, blocks = list(), residuals,
                       n_sweeps = 1000, n_burn = 0, n_adapts = 0, seed = 1,
                       x_init = NULL, adtype = "forwarddiff",
                       n_keep = 500, thin = NULL, quiet = FALSE) {
  et_require_session()
  if (!is.list(residuals) || !length(residuals)) {
    stop("et_collect(): `residuals` must be a non-empty list.", call. = FALSE)
  }
  for (r in residuals) {
    if (!inherits(r, "et_residual")) {
      stop("et_collect(): every entry must come from et_residual_survival() ",
           "or et_residual_waiting().", call. = FALSE)
    }
  }
  if (is.null(thin)) thin <- max(1L, as.integer(n_sweeps %/% n_keep))

  # The collector has to see every sweep, so the sweep loop is the one place the
  # generated `et_run` cannot own. It is driven from the generated module with
  # the residual specs built in the same scope as everything else they name.
  gen <- et_julia_source(model, blocks)
  mod <- et_load_module(gen, NULL)

  states <- model$data$transitions$states
  specs <- paste(vapply(residuals, residual_to_julia, character(1),
                        states = states),
                 collapse = ",\n    ")
  names_v <- vapply(residuals, function(r) r$name, character(1))

  x_arg <- "nothing"
  if (!is.null(x_init)) {
    x_init <- as.matrix(x_init)
    storage.mode(x_init) <- "integer"
    check_x_init(x_init, model)
    JuliaCall::julia_assign(paste0(mod, "_xinit"), x_init)
    x_arg <- sprintf("Matrix{Int}(Main.%s_xinit)", mod)
  }

  if (!quiet) {
    message(sprintf("Running %d sweeps, scoring every %d (keeping up to %d) ...",
                    n_sweeps, thin, n_keep))
  }
  body <- sprintf(
    "et_collect_run(; specs=(\n    %s,\n  ), names=%s, n_sweeps=%s, n_burn=%s, n_adapts=%s, seed=%s, x_init=%s, adtype=%s, n_keep=%s, thin=%s)",
    specs, julia_symbol_vector(names_v), julia_int(n_sweeps), julia_int(n_burn),
    julia_int(n_adapts), julia_int(seed), x_arg, adtype_to_julia(adtype),
    julia_int(n_keep), julia_int(thin))
  src_sym <- paste0(mod, "_collect_src")
  JuliaCall::julia_assign(src_sym, body)
  raw <- JuliaCall::julia_eval(sprintf("Base.include_string(%s, Main.%s)",
                                       mod, src_sym))

  elapsed <- raw[["_elapsed"]]
  raw[["_elapsed"]] <- NULL
  resid_raw <- raw[["_residuals"]]
  raw[["_residuals"]] <- NULL

  draws <- restore_shapes(raw, model)
  if (!quiet) {
    message(sprintf("done in %.1f s (%.3f s/sweep)", elapsed,
                    elapsed / max(1, n_sweeps + n_burn)))
  }
  fit <- structure(list(draws = draws, model = model, blocks = gen$blocks,
                        depends = gen$depends, source = gen, module = mod,
                        x_path = NULL, elapsed = elapsed, n_sweeps = n_sweeps,
                        n_burn = n_burn, seed = seed), class = "et_fit")
  list(fit = fit, residuals = residual_frame(resid_raw, names_v))
}

# Shared by et_collect() and et_residuals(): the per-name individuals x draws
# matrices Julia returns, as one long data frame.
residual_frame <- function(raw, names_v) {
  out <- do.call(rbind, lapply(names_v, function(nm) {
    m <- raw[[nm]]
    if (is.null(m)) return(NULL)
    m <- as.matrix(m)
    data.frame(residual = nm,
               draw = rep(seq_len(ncol(m)), each = nrow(m)),
               individual = rep(seq_len(nrow(m)), times = ncol(m)),
               value = as.numeric(m), stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

#' Case reproduction numbers, R_i.
#'
#' The expected number of secondary infections caused by each individual, from a
#' fit's archived trajectories.
#'
#' **This is not a residual, and it is not computed one individual at a time.**
#' Attributing an infection needs to know every infective present in the victim's
#' group at that moment, so R_i takes a pass over the whole population before any
#' one individual's value is known. That is why it has its own entry point rather
#' than joining [et_residuals()].
#'
#' ## The declared decomposition
#'
#' The package cannot split a force of infection by itself: `S -> E` is one rate
#' function returning one number, and nothing in it says which part is background
#' and which is transmission. So you declare the parts, exactly as you declare
#' rates and aggregates -- as ordinary R rate functions of `(model, data, i, t)`
#' -- and say which one is transmission via `secondary`. The package sums them,
#' normalises, and never needs to know what any of them means.
#'
#' For each infection event, the transmission component of the victim's force of
#' infection is shared among the infectives in its group in proportion to
#' `weight`. Background infections are attributed to nobody, which is the point
#' of the split. An individual that was never infectious gets `NA` -- it had no
#' opportunity to infect anyone, which differs from having had the opportunity
#' and infected nobody (a genuine `0`).
#'
#' @param fit An [et_fit()] from [et_sample()] run with `save_x =`.
#' @param components Named list of rate functions -- the force-of-infection
#'   decomposition. Each is written in the same R subset as any other rate
#'   function, taking `(model, data, i, t)`.
#' @param secondary Name of the component that is transmission (the part
#'   attributable to other individuals).
#' @param infection The transition marking a new infection, e.g. `c("S", "E")`.
#' @param infectious_state State from which an individual can infect others.
#' @param weight Optional relative infectiousness, a function
#'   `(model, data, X, i, t)`. Defaults to uniform, i.e. equal attribution.
#' @param draws Which archived draws to use. Defaults to all of them.
#' @return A data frame with one row per individual per draw: `draw`,
#'   `individual`, `value` (`NA` where the individual was never infectious).
#' @export
et_reproduction_numbers <- function(fit, components, secondary,
                                    infection = c("S", "E"),
                                    infectious_state = "I",
                                    weight = NULL, draws = NULL) {
  et_require_session()
  if (!inherits(fit, "et_fit")) {
    stop("et_reproduction_numbers(): `fit` must come from et_sample().",
         call. = FALSE)
  }
  if (is.null(fit$x_path)) {
    stop("et_reproduction_numbers(): this fit did not archive its trajectories. ",
         "Re-run et_sample() with `save_x = \"X.rds\"`.", call. = FALSE)
  }
  if (!is.list(components) || !length(components) ||
      is.null(names(components)) || any(!nzchar(names(components)))) {
    stop("et_reproduction_numbers(): `components` must be a NAMED list of rate ",
         "functions.", call. = FALSE)
  }
  if (!secondary %in% names(components)) {
    stop("et_reproduction_numbers(): `secondary` must name one of `components` (",
         paste(names(components), collapse = ", "), "); got '", secondary, "'.",
         call. = FALSE)
  }
  if (length(infection) != 2L) {
    stop("et_reproduction_numbers(): `infection` must be two state names, ",
         "from and to.", call. = FALSE)
  }

  # Each component is written as a RATE function -- same subset, same
  # `(model, data, i, t)` contract as anything in et_transitions() -- so it goes
  # through the same transpiler with the same role. Emitting them as named
  # functions (rather than inline lambdas) keeps a transpile error pointing at
  # the component's own name.
  #
  # ET calls the components with X as well, `(model, data, X, i, t)`, because a
  # component is free to look at the trajectory. The transpiled body does not
  # take X, so each is wrapped to drop it. Keeping the USER-FACING contract equal
  # to every other rate function is worth one wrapper here.
  comp_names <- paste0("etr_rn_comp_", names(components))
  comp_defs <- vapply(seq_along(components), function(k) {
    et_transpile_role(components[[k]], "rate", comp_names[k])$src
  }, character(1))
  comp_src <- sprintf("%s = (mo, da, X, i, t) -> %s(mo, da, i, t)",
                      names(components), comp_names)

  # `weight` is relative infectiousness at (i, t). It is not one of the declared
  # protocols -- it takes X as well as (i, t) -- so it must be an et_julia().
  if (is.null(weight)) {
    wt_src <- "nothing"
  } else if (inherits(weight, "et_julia")) {
    wt_src <- sprintf("(model, data, X, i, t) -> begin\n%s\nend", weight$src)
  } else {
    stop("et_reproduction_numbers(): `weight` must be an et_julia() object -- ",
         "it takes (model, data, X, i, t), which is not one of the transpiler's ",
         "declared protocols. Leave it NULL for equal attribution.",
         call. = FALSE)
  }

  files <- et_chunk_files(et_with_ext(fit$x_path, "rds"))
  from_rds <- length(files) > 0
  if (!from_rds) files <- et_chunk_files(et_with_ext(fit$x_path, "jld2"))
  if (!length(files)) {
    stop("et_reproduction_numbers(): no archived trajectories found for '",
         fit$x_path, "'.", call. = FALSE)
  }

  JuliaCall::julia_command("using JLD2")
  if (from_rds) {
    stage <- file.path(tempdir(), paste0("etr-rstage-", fit$module))
    dir.create(stage, showWarnings = FALSE, recursive = TRUE)
    on.exit(unlink(stage, recursive = TRUE), add = TRUE)
    files <- et_stage_rds_as_jld2(files, stage)
  }

  flat <- flatten_draws(fit)
  par_sym <- paste0(fit$module, "_rn_pars")
  JuliaCall::julia_assign(par_sym, flat$values)
  files_sym <- paste0(fit$module, "_rn_files")
  JuliaCall::julia_assign(files_sym, julia_path(files))

  # The component definitions and the call go across together: they are parsed
  # and evaluated INSIDE the generated module, where `model$...` reads, helpers
  # and the extras they name are already in scope.
  body <- sprintf(paste0(
    "%s\n\n",
    "et_reproduction_numbers(Main.%s; components=(%s,), secondary=%s, ",
    "infection=(%s => %s), infectious_state=%s, weight=%s, ",
    "par_draws=_et_unflatten(Main.%s), draws=%s)"),
    paste(comp_defs, collapse = "\n\n"),
    files_sym, paste(comp_src, collapse = ", "), julia_symbol(secondary),
    julia_symbol(infection[1]), julia_symbol(infection[2]),
    julia_symbol(infectious_state), wt_src, par_sym,
    if (is.null(draws)) "nothing" else julia_vector(as.integer(draws), "Int"))
  src_sym <- paste0(fit$module, "_rn_src")
  JuliaCall::julia_assign(src_sym, body)
  raw <- JuliaCall::julia_eval(sprintf("Base.include_string(%s, Main.%s)",
                                       fit$module, src_sym))

  m <- as.matrix(raw[["R"]])          # individuals x draws
  data.frame(draw = rep(seq_len(ncol(m)), each = nrow(m)),
             individual = rep(seq_len(nrow(m)), times = ncol(m)),
             value = as.numeric(m), stringsAsFactors = FALSE)
}
