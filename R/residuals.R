# Layer 3: residuals -- the diagnostic half.
#
# A fit tells you what the parameters are; a residual tells you whether the model
# that produced them is any good. Each of these is a PIT residual: uniform on
# (0, 1) if the mechanism it targets is right, so a departure from uniformity is
# evidence of misfit in THAT mechanism rather than a diffuse "the model is bad".
#
# The trajectories they score come from the fit. Two routes, and the R surface
# takes the post-hoc one: sample with `save_x =`, which streams the trajectory to
# disc, then compute residuals from the archive afterwards. That decouples the
# diagnostics from the fit -- new residuals can be added without refitting -- and
# keeps `et_sample()` a single call.

#' Residual for how long an individual survives.
#'
#' Wraps ET's `LeftTruncatedSurvivalResidual`. Scores the realised lifetime
#' against the model's survival function.
#'
#' **Left truncation is not optional.** An individual is only ever in the data
#' because it survived long enough to be observed, so the residual has to be
#' conditioned on that. Without it, it is calibrated against the unconditional
#' lifetime and reports the sampling design as model misfit.
#'
#' @param origin Clock origin per individual: the name of an `extras` vector, or
#'   a function of the individual index. Typically birth.
#' @param condition_on A time the individual is KNOWN to have been alive --
#'   typically first capture. The CDF is renormalised by survival to here.
#' @param censor_at When observation ends for an individual that does not die
#'   within the trajectory. Typically last capture.
#' @param death Name of the death state.
#' @param survival Optional: the survival function, when it differs from the one
#'   the transitions declare. Defaults to the transitions' own, which is almost
#'   always right -- they must agree, or the residual scores a different model
#'   than the one you fitted.
#' @param name Result key.
#' @return An object of class `et_residual`.
#' @export
et_residual_survival <- function(origin, condition_on, censor_at,
                                 death = "D", survival = NULL,
                                 name = "survival") {
  check_julia_name(name, "residual name")
  structure(list(kind = "survival", name = name, origin = origin,
                 condition_on = condition_on, censor_at = censor_at,
                 death = death, survival = survival),
            class = "et_residual")
}

#' Residual for how long an individual waits before making a transition.
#'
#' Wraps ET's `WaitingTimeResidual`. Scores the realised waiting time for a
#' `from -> to` move against the hazard the model spec implies.
#'
#' @param from,to State names.
#' @param origin Where the clock starts:
#'   * `"entry_to_from_state"` -- when the individual entered `from`, read from
#'     the trajectory. Right for a latent period: the clock starts at exposure.
#'   * `"window_start"` -- when observation started. Right for an exposure time:
#'     a susceptible individual was already susceptible when watching began, so
#'     its clock cannot start at entry to `S`.
#' @param censor_at Character vector naming what ends observation.
#'   `"window_end"` is the individual's own sampling window; any other entry is
#'   read as a STATE, censoring at the last step before it was entered. The
#'   earliest applicable one wins.
#'
#'   A competing risk belongs here rather than being dropped: an individual that
#'   died before making the move has not falsified the model, it stopped being
#'   observable. Dropping it biases the residual towards individuals who lived
#'   long enough.
#' @param accumulate `"discrete_product"` (a product of per-step probabilities,
#'   matching a discrete-time model) or `"cumulative_hazard"` (the Sellke form).
#'   No default: the two agree only when the per-step probability is
#'   `1 - exp(-hazard)`, and picking the wrong one is a silent miscalibration.
#' @param name Result key. Defaults to `from_to_to`.
#' @return An object of class `et_residual`.
#' @export
et_residual_waiting <- function(from, to,
                                origin = c("entry_to_from_state", "window_start"),
                                censor_at = "window_end",
                                accumulate = c("discrete_product",
                                               "cumulative_hazard"),
                                name = NULL) {
  origin <- match.arg(origin)
  accumulate <- match.arg(accumulate)
  if (is.null(name)) name <- paste0(from, "_to_", to)
  check_julia_name(name, "residual name")
  structure(list(kind = "waiting", name = name, from = from, to = to,
                 origin = origin, censor_at = as.character(censor_at),
                 accumulate = accumulate),
            class = "et_residual")
}

#' @export
print.et_residual <- function(x, ...) {
  cat("<et_residual> ", x$name, "  (", x$kind, ")\n", sep = "")
  invisible(x)
}

# Render one residual as the Julia constructor call.
residual_to_julia <- function(r, states) {
  if (r$kind == "waiting") {
    return(sprintf(
      "WaitingTimeResidual(%s => %s; accumulate=%s, origin=%s, censor_at=(%s), name=%s)",
      julia_symbol(r$from), julia_symbol(r$to), julia_symbol(r$accumulate),
      julia_symbol(r$origin),
      paste0(paste(vapply(r$censor_at, julia_symbol, character(1)),
                   collapse = ", "), ","),
      julia_symbol(r$name)))
  }
  sprintf(paste0("LeftTruncatedSurvivalResidual(survival=%s, origin=%s, ",
                 "condition_on=%s, censor_at=%s, death=%s, name=%s)"),
          survival_src(r$survival),
          index_fn_src(r$origin), index_fn_src(r$condition_on),
          index_fn_src(r$censor_at), julia_symbol(r$death),
          julia_symbol(r$name))
}

# A per-individual time: the name of an extras vector, or a literal vector.
#
# `max(..., 1)` because a clock origin can legitimately predate the observation
# window (a badger born before monitoring began), while the arrays the survival
# function indexes only exist over `1:n_timepoints`. The age is still right --
# `age[i, t]` is `t - birth[i]`, so an individual born earlier is simply older at
# t = 1 -- but the product has to start inside the array.
index_fn_src <- function(x) {
  if (is.character(x) && length(x) == 1L) {
    return(sprintf("i -> max(%s[i], 1)", x))
  }
  if (is.numeric(x)) {
    return(sprintf("i -> max(%s[i], 1)", julia_vector(as.integer(x), "Int")))
  }
  if (inherits(x, "et_julia")) return(x$src)
  stop("a residual time must be the name of an `extras` vector, a numeric ",
       "vector, or an et_julia() expression; got ", class(x)[1], ".",
       call. = FALSE)
}

# The survival function the residual scores against. Defaults to the one the
# transitions declared -- they MUST agree, or the residual scores a different
# model than the one that was fitted.
survival_src <- function(x) {
  if (is.null(x)) return("et_survival_fn")
  if (inherits(x, "et_julia")) return(x$src)
  as.character(x)
}

#' Compute residuals from a fit's archived trajectories.
#'
#' Reads the trajectories [et_sample()] streamed to disc, pairs each with the
#' parameter draw that produced it, and evaluates the residuals.
#'
#' Requires the fit to have been run with `save_x =`. Computing residuals
#' post-hoc means they can be recomputed, or new ones added, without refitting.
#'
#' @param fit An [et_fit()] from [et_sample()] run with `save_x =`.
#' @param residuals List of [et_residual_survival()] / [et_residual_waiting()].
#' @param sync_aggregates Rebuild the tracked arrays for each draw before scoring
#'   it. Needed when any rate reads `data$aggregates` -- a spec-derived hazard
#'   then depends on them agreeing with THIS draw's trajectory, not with whatever
#'   was current when the fit ended. Costs one pass per draw, which is why it is
#'   not unconditional.
#' @param seed RNG seed for the randomised PIT. Its own stream, so it never
#'   perturbs anything else.
#' @return A data frame with one row per residual per scoreable individual per
#'   draw: `residual`, `draw`, `individual`, `value`. `NA` where an individual
#'   had no scoreable event.
#' @export
et_residuals <- function(fit, residuals, sync_aggregates = TRUE, seed = 1) {
  et_require_session()
  if (!inherits(fit, "et_fit")) {
    stop("et_residuals(): `fit` must come from et_sample().", call. = FALSE)
  }
  if (is.null(fit$x_path)) {
    stop("et_residuals(): this fit did not archive its trajectories, so there is ",
         "nothing to score. Re-run et_sample() with `save_x = \"X.jld2\"`.",
         call. = FALSE)
  }
  if (!is.list(residuals) || !length(residuals)) {
    stop("et_residuals(): `residuals` must be a non-empty list.", call. = FALSE)
  }
  for (r in residuals) {
    if (!inherits(r, "et_residual")) {
      stop("et_residuals(): every entry must come from et_residual_survival() ",
           "or et_residual_waiting().", call. = FALSE)
    }
  }

  states <- fit$model$data$transitions$states
  specs <- paste(vapply(residuals, residual_to_julia, character(1),
                        states = states),
                 collapse = ",\n    ")
  names_v <- vapply(residuals, function(r) r$name, character(1))

  # Each archived trajectory must be scored under the parameters that produced
  # it, so the draws go across as ONE flat numeric matrix (draws x scalars) and
  # the generated module rebuilds the NamedTuples from it. A matrix marshals
  # unambiguously; a list of NamedTuples does not.
  flat <- flatten_draws(fit)
  sym <- paste0(fit$module, "_resid_pars")
  JuliaCall::julia_assign(sym, flat$values)

  # Parsed AND evaluated inside the generated module: the residual constructors,
  # the survival function and the `extras` vectors all live there, so building
  # the spec anywhere else would mean qualifying every one of them.
  # `include_string` is what puts both halves in the module's scope -- `Base.eval`
  # with a `quote` would still resolve the names where the quote was written.
  # The per-flush chunks, in sweep order. Discovery happens HERE because R is
  # the side that knows which format the archive is in; Julia is handed an
  # explicit ordered file list and reads them one at a time. Passing the list
  # rather than a template also keeps the ordering decision in one place.
  files <- et_chunk_files(et_with_ext(fit$x_path, "rds"))
  from_rds <- length(files) > 0
  if (!from_rds) files <- et_chunk_files(et_with_ext(fit$x_path, "jld2"))
  if (!length(files)) {
    stop("et_residuals(): no archived trajectories found for '", fit$x_path,
         "'. If the fit was interrupted, the flushes already on disc can be ",
         "recovered with et_convert_trajectories().", call. = FALSE)
  }

  JuliaCall::julia_command("using JLD2")

  if (from_rds) {
    # An `.rds` chunk is an R object, so Julia cannot open it. Rather than
    # SCORE it separately -- a second scoring path that must be kept in step
    # with the first, and silently was not -- the chunks are staged back to
    # temporary `.jld2` and fed to the SAME lazy reader. One scoring
    # implementation, so the two formats cannot disagree by construction.
    #
    # Staging is one chunk at a time and the copies are deleted together
    # afterwards, so peak extra disc is one chunk, and peak memory is unchanged.
    stage <- file.path(tempdir(), paste0("etr-stage-", fit$module))
    dir.create(stage, showWarnings = FALSE, recursive = TRUE)
    on.exit(unlink(stage, recursive = TRUE), add = TRUE)
    files <- et_stage_rds_as_jld2(files, stage)
  }
  files_sym <- paste0(fit$module, "_resid_files")
  JuliaCall::julia_assign(files_sym, julia_path(files))
  body <- sprintf(
    "et_residuals(Main.%s; specs=(\n    %s,\n  ), names=%s, par_draws=_et_unflatten(Main.%s), sync=%s, seed=%s)",
    files_sym, specs, julia_symbol_vector(names_v), sym,
    if (isTRUE(sync_aggregates)) "true" else "false", julia_int(seed))
  # The source goes across as a VALUE, not embedded in another Julia string
  # literal -- nesting one inside the other double-escapes every quote in it.
  src_sym <- paste0(fit$module, "_resid_src")
  JuliaCall::julia_assign(src_sym, body)
  raw <- JuliaCall::julia_eval(sprintf("Base.include_string(%s, Main.%s)",
                                       fit$module, src_sym))

  out <- do.call(rbind, lapply(names_v, function(nm) {
    m <- raw[[nm]]                      # individuals x draws
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

# Flatten a fit's sampled parameter draws to one numeric matrix, draws x scalars,
# in the model's declared parameter order. The generated `_et_unflatten` reads it
# back using the same order and the declared shapes, so the two cannot disagree
# about which column is which.
#
# Only SAMPLED parameters: deterministics are recomputed from them, so carrying
# them would be redundant and could disagree.
flatten_draws <- function(fit) {
  model <- fit$model
  cols <- list()
  for (nm in model$par_names) {
    v <- fit$draws[[nm]]
    if (is.null(v)) {
      stop("et_residuals(): the fit has no draws for '", nm,
           "', so its trajectories cannot be scored.", call. = FALSE)
    }
    cols[[nm]] <- if (is.null(dim(v))) matrix(as.numeric(v), ncol = 1L)
                  else matrix(as.numeric(v), nrow = dim(v)[1])
  }
  list(values = do.call(cbind, cols))
}

#' Summarise residuals.
#'
#' A PIT residual is Uniform(0, 1) under a correct model, so `mean` near 0.5 and
#' a flat distribution are what "no evidence of misfit" looks like. This reports
#' the mean and the proportion in each quartile; a strong skew points at the
#' mechanism that residual targets.
#'
#' @param x A data frame from [et_residuals()].
#' @return A data frame, one row per residual.
#' @export
et_residual_summary <- function(x) {
  stopifnot(is.data.frame(x))
  do.call(rbind, lapply(split(x, x$residual), function(d) {
    v <- d$value[!is.na(d$value)]
    if (!length(v)) {
      return(data.frame(residual = d$residual[1], n = 0L, mean = NA_real_,
                        q1 = NA_real_, q2 = NA_real_, q3 = NA_real_,
                        q4 = NA_real_, stringsAsFactors = FALSE))
    }
    br <- cut(v, breaks = c(-Inf, 0.25, 0.5, 0.75, Inf), labels = FALSE)
    p <- tabulate(br, nbins = 4L) / length(v)
    data.frame(residual = d$residual[1], n = length(v), mean = mean(v),
               q1 = p[1], q2 = p[2], q3 = p[3], q4 = p[4],
               stringsAsFactors = FALSE)
  }))
}

# Stage `.rds` chunks back to temporary `.jld2`, so the ONE lazy Julia reader
# can score them.
#
# The alternative -- a second, R-side scoring loop -- was tried and dropped: it
# silently disagreed with the Julia path (93/240 residual values differed, while
# `sync=FALSE` agreed exactly), because keeping two scoring implementations in
# step is precisely the thing that does not stay true. Staging costs one file
# copy per flush, off the hot path, and cannot drift.
#
# One chunk at a time: read it, write it, drop it. Peak extra disc is one
# chunk's worth, peak memory unchanged.
et_stage_rds_as_jld2 <- function(files, stage_dir) {
  out <- character(length(files))
  for (k in seq_along(files)) {
    chunk <- readRDS(files[k])
    if (!is.list(chunk)) chunk <- list(chunk)
    # Keep the `_iters_x_to_y` name: the Julia side is handed an explicit
    # ordered list, but a faithful name keeps any error message meaningful.
    dest <- file.path(stage_dir, et_with_ext(basename(files[k]), "jld2"))
    if (!length(chunk)) {
      JuliaCall::julia_command(sprintf(
        'JLD2.jldopen(%s, "w") do f; f["states"] = Matrix{Int}[]; end; nothing',
        julia_string(julia_path(dest))))
    } else {
      arr <- array(unlist(lapply(chunk, function(m) {
        storage.mode(m) <- "integer"; m
      })), dim = c(nrow(chunk[[1]]), ncol(chunk[[1]]), length(chunk)))
      JuliaCall::julia_assign("etr_stage_arr", arr)
      JuliaCall::julia_command(sprintf(paste0(
        'JLD2.jldopen(%s, "w") do f; ',
        'f["states"] = Matrix{Int}[Matrix{Int}(Main.etr_stage_arr[:, :, k]) ',
        'for k in axes(Main.etr_stage_arr, 3)]; end; nothing'),
        julia_string(julia_path(dest))))
      rm(arr)
    }
    rm(chunk)
    out[k] <- dest
  }
  JuliaCall::julia_command("Main.etr_stage_arr = nothing; nothing")
  out
}
