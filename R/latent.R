# Weighted chains with sampled latents: the R mirror of ET's latent_hmm_data().
#
# Each individual's hidden chain is summed out exactly, as in a marginal model,
# but its weights and observations may read latent variables (an activity
# centre per year, a birth year, an unmapped detector's position) that are
# sampled, once per Gibbs sweep, outside every gradient. Individuals are
# independent given the parameters and the latents.

# ---- latent declarations ----------------------------------------------------

check_coord_prior <- function(prior, dim, what) {
  if (inherits(prior, "custom_dist")) return(invisible(TRUE))
  if (is.list(prior) && length(prior) == dim &&
      all(vapply(prior, inherits, logical(1), "custom_dist"))) {
    return(invisible(TRUE))
  }
  stop(what, ": `prior` must be one distribution for every coordinate, or a ",
       "list of ", dim, ", one per coordinate.", call. = FALSE)
}

check_dims <- function(dim, n, what) {
  for (v in list(dim, n)) {
    if (!is.numeric(v) || length(v) != 1L || is.na(v) || v != round(v) || v < 1) {
      stop(what, ": `dim` and `n` must be positive whole numbers.", call. = FALSE)
    }
  }
}

#' Latent variables of a weighted chain.
#'
#' The three kinds of latent [et_latent_data()] samples, mirroring ET's
#' `Field`, `PerIndividual` and `Shared`. Collect them with [et_latents()].
#'
#' * `et_field()`: a `dim`-vector for each individual at each of `n` elements,
#'   e.g. an activity centre per year. Held as a `dim x n x n_individuals`
#'   array, read as `data$latents$name[c, e, i]`.
#' * `et_per_individual()`: one value per individual (a random effect, a birth
#'   year, sex), read as `data$latents$name[i]`.
#' * `et_shared_latent()`: `n` vectors common to every individual, e.g. the
#'   positions of unmapped detectors, held as a `dim x n` matrix.
#'
#' @param dim,n Size: `dim` coordinates for each of `n` elements.
#' @param prior A distribution such as [normal_dist()], or for `et_field()` and
#'   `et_shared_latent()` a list of `dim` of them, one per coordinate. For
#'   `et_per_individual()` it may instead be `function(model) <distribution>`
#'   when it depends on the parameters, e.g.
#'   `function(model) bernoulli_dist(model$sex_prop)`; such a prior is part of
#'   the likelihood, for every individual.
#' @param index `function(t)`: the element in use at timepoint `t`, e.g. the
#'   year. The weights and observations at `t` may read elements `index(t)` and
#'   `index(t) - 1`.
#' @param touches Instead of `index`, `function(e)` returning the range of
#'   timepoints (`a:b`) whose weights or observations read element `e`.
#' @param scale Starting random-walk scale, adapted during warm-up.
#' @param init Starting values: an array of the latent's shape, or for a field
#'   [et_starting_centres()] and for a shared latent [et_starting_detectors()].
#'   Left out, the latent is drawn from its prior (a parameter-dependent prior
#'   from the `observed` values).
#' @param observed A vector with `NA` where the value is unknown; known values
#'   are fixed.
#' @param allowed Optional `function(data, i, value)` returning `TRUE` where a
#'   value is allowed; moves that break it are rejected.
#' @return An object of class `et_latent`.
#' @name latent_kinds
NULL

#' @rdname latent_kinds
#' @export
et_field <- function(dim, n, prior, index = NULL, touches = NULL, scale = 1,
                     init = NULL) {
  check_dims(dim, n, "et_field()")
  check_coord_prior(prior, dim, "et_field()")
  if (is.null(index) == is.null(touches)) {
    stop("et_field(): give one of `index` or `touches`.", call. = FALSE)
  }
  if (!is.null(init) && !inherits(init, "et_starting_centres") && !is.numeric(init)) {
    stop("et_field(): `init` must be a numeric array or et_starting_centres().",
         call. = FALSE)
  }
  structure(list(kind = "field", dim = as.integer(dim), n = as.integer(n),
                 prior = prior, index = index, touches = touches,
                 scale = as.numeric(scale), init = init), class = "et_latent")
}

#' @rdname latent_kinds
#' @export
et_per_individual <- function(prior, observed = NULL, allowed = NULL, scale = 1,
                              init = NULL) {
  if (!inherits(prior, "custom_dist") && !is.function(prior)) {
    stop("et_per_individual(): `prior` must be a distribution or ",
         "function(model) <distribution>.", call. = FALSE)
  }
  if (!is.null(allowed) && (!is.function(allowed) || length(formals(allowed)) != 3L)) {
    stop("et_per_individual(): `allowed` must be function(data, i, value).",
         call. = FALSE)
  }
  if (!is.null(observed) && !is.numeric(observed) && !is.logical(observed)) {
    stop("et_per_individual(): `observed` must be a numeric vector with NA ",
         "where the value is unknown.", call. = FALSE)
  }
  structure(list(kind = "per_individual", prior = prior, observed = observed,
                 allowed = allowed, scale = as.numeric(scale), init = init),
            class = "et_latent")
}

#' @rdname latent_kinds
#' @export
et_shared_latent <- function(dim, n, prior, scale = 1, init = NULL) {
  check_dims(dim, n, "et_shared_latent()")
  check_coord_prior(prior, dim, "et_shared_latent()")
  if (!is.null(init) && !inherits(init, "et_starting_detectors") && !is.numeric(init)) {
    stop("et_shared_latent(): `init` must be a numeric matrix or ",
         "et_starting_detectors().", call. = FALSE)
  }
  structure(list(kind = "shared", dim = as.integer(dim), n = as.integer(n),
                 prior = prior, scale = as.numeric(scale), init = init),
            class = "et_latent")
}

#' Collect a model's latent variables.
#'
#' The R counterpart of ET's `@latents` block: one named argument per latent.
#'
#' @param ... Named [latent_kinds] declarations.
#' @return An object of class `et_latents`.
#' @export
#' @examples
#' latents <- et_latents(
#'   S = et_field(2, 10, prior = normal_dist(0, 5), index = function(t) (t + 3) %/% 4),
#'   u = et_per_individual(normal_dist(0, 1), scale = 0.5))
et_latents <- function(...) {
  ls <- list(...)
  if (!length(ls) || is.null(names(ls)) || any(!nzchar(names(ls)))) {
    stop("et_latents(): give each latent as a named argument, e.g. ",
         "S = et_field(...).", call. = FALSE)
  }
  check_julia_name(names(ls), "latent name")
  if (anyDuplicated(names(ls))) {
    stop("et_latents(): duplicate latent name(s).", call. = FALSE)
  }
  for (nm in names(ls)) {
    if (!inherits(ls[[nm]], "et_latent")) {
      stop("et_latents(): '", nm, "' must come from et_field(), ",
           "et_per_individual() or et_shared_latent().", call. = FALSE)
    }
  }
  structure(ls, class = "et_latents")
}

#' @export
print.et_latents <- function(x, ...) {
  cat("<et_latents>\n")
  for (nm in names(x)) {
    k <- x[[nm]]
    cat("  ", nm, " ~ ", switch(k$kind,
      field = sprintf("field %d x %d", k$dim, k$n),
      per_individual = "per individual",
      shared = sprintf("shared %d x %d", k$dim, k$n)), "\n", sep = "")
  }
  invisible(x)
}

# ---- spatial captures -------------------------------------------------------

#' Spatial capture-recapture data.
#'
#' Detectors at known or latent positions, each active on some occasions, and
#' capture counts per individual, occasion and detector. Pass the result in
#' [et_latent_data()]'s `extras`, and score it in an observation function with
#' [et_capture_loglik()].
#'
#' @param detector_xy A `2 x J` matrix of detector positions. Columns of latent
#'   detectors are ignored and may be `NA`.
#' @param active A `J x n_occasions` logical matrix: whether each detector was
#'   trapping on each occasion.
#' @param captures A data frame (or four-column matrix) with columns `i`,
#'   `occasion`, `detector` and `n`: one row for each detector that caught an
#'   individual on an occasion.
#' @param n_individuals The number of individuals, augmented ones included.
#' @param latent_detectors Indices of the detectors whose positions are unknown.
#'   Their positions are a [et_shared_latent()] of `dim = 2`, in this order.
#' @return An object of class `et_spatial_captures`.
#' @export
et_spatial_captures <- function(detector_xy, active, captures, n_individuals,
                                latent_detectors = integer()) {
  detector_xy <- as.matrix(detector_xy)
  active <- as.matrix(active)
  J <- nrow(active)
  if (!identical(dim(detector_xy), c(2L, J))) {
    stop("et_spatial_captures(): `detector_xy` must be 2 x ", J,
         " (one column per row of `active`), got ",
         paste(dim(detector_xy), collapse = " x "), ".", call. = FALSE)
  }
  if (!is.logical(active)) {
    stop("et_spatial_captures(): `active` must be a logical matrix.", call. = FALSE)
  }
  latent_detectors <- as.integer(latent_detectors)
  if (any(is.na(latent_detectors)) || any(latent_detectors < 1L | latent_detectors > J) ||
      anyDuplicated(latent_detectors)) {
    stop("et_spatial_captures(): `latent_detectors` must be distinct indices in ",
         "1..", J, ".", call. = FALSE)
  }
  mapped <- setdiff(seq_len(J), latent_detectors)
  if (any(!is.finite(detector_xy[, mapped]))) {
    stop("et_spatial_captures(): mapped detector(s) ",
         paste(utils::head(mapped[colSums(!is.finite(detector_xy[, mapped, drop = FALSE])) > 0], 5),
               collapse = ", "),
         " have no position; list them in `latent_detectors`.", call. = FALSE)
  }
  detector_xy[!is.finite(detector_xy)] <- NaN
  storage.mode(detector_xy) <- "double"
  captures <- as.data.frame(captures)
  if (!all(c("i", "occasion", "detector", "n") %in% names(captures))) {
    if (ncol(captures) == 4L) names(captures) <- c("i", "occasion", "detector", "n")
    else stop("et_spatial_captures(): `captures` needs columns i, occasion, ",
              "detector and n.", call. = FALSE)
  }
  n_individuals <- as_count(n_individuals, "n_individuals")
  cols <- lapply(captures[c("i", "occasion", "detector", "n")], function(v) {
    if (any(is.na(v)) || any(v != round(v))) {
      stop("et_spatial_captures(): `captures` must hold whole numbers.", call. = FALSE)
    }
    as.integer(v)
  })
  if (any(cols$i < 1L | cols$i > n_individuals) ||
      any(cols$occasion < 1L | cols$occasion > ncol(active)) ||
      any(cols$detector < 1L | cols$detector > J) || any(cols$n < 0L)) {
    stop("et_spatial_captures(): a capture row is out of range (individual, ",
         "occasion or detector), or has a negative count.", call. = FALSE)
  }
  off <- cols$n > 0L & !active[cbind(cols$detector, cols$occasion)]
  if (any(off)) {
    stop("et_spatial_captures(): ", sum(off), " capture(s) at a detector that ",
         "was not active on that occasion.", call. = FALSE)
  }
  structure(list(detector_xy = detector_xy, active = active, captures = cols,
                 n_individuals = n_individuals, latent_detectors = latent_detectors),
            class = "et_spatial_captures")
}

#' @export
print.et_spatial_captures <- function(x, ...) {
  cat("<et_spatial_captures> ", ncol(x$detector_xy), " detectors (",
      length(x$latent_detectors), " latent), ", ncol(x$active), " occasions, ",
      sum(x$captures$n), " captures of ", x$n_individuals, " individuals\n", sep = "")
  invisible(x)
}

#' Starting values from the capture data.
#'
#' `et_starting_centres()` starts an [et_field()] of activity centres at each
#' individual's mean capture position on each element, carried across elements
#' without captures; individuals never caught at a mapped detector start at one
#' uniform position in `bounds`. `et_starting_detectors()` starts the latent
#' detectors at the mean starting centre of the individuals caught there. Both
#' are computed by ET when the model is built.
#'
#' @param captures An [et_spatial_captures()] object, also passed in `extras`.
#' @param index `function(occasion)`: the element in use on an occasion.
#' @param n The number of elements.
#' @param bounds `c(xmin, xmax, ymin, ymax)`.
#' @param centres An [et_starting_centres()] object.
#' @return A marker object, given as a latent's `init`.
#' @export
et_starting_centres <- function(captures, index, n, bounds) {
  check_starting(captures, index, bounds, "et_starting_centres()")
  structure(list(captures = captures, index = index, n = as.integer(n),
                 bounds = as.numeric(bounds)), class = "et_starting_centres")
}

#' @rdname et_starting_centres
#' @export
et_starting_detectors <- function(captures, centres, index, bounds) {
  check_starting(captures, index, bounds, "et_starting_detectors()")
  if (!inherits(centres, "et_starting_centres")) {
    stop("et_starting_detectors(): `centres` must come from et_starting_centres().",
         call. = FALSE)
  }
  structure(list(captures = captures, centres = centres, index = index,
                 bounds = as.numeric(bounds)), class = "et_starting_detectors")
}

check_starting <- function(captures, index, bounds, what) {
  if (!inherits(captures, "et_spatial_captures")) {
    stop(what, ": `captures` must come from et_spatial_captures().", call. = FALSE)
  }
  if (!is.function(index) || length(formals(index)) != 1L) {
    stop(what, ": `index` must be a function of one argument.", call. = FALSE)
  }
  if (!is.numeric(bounds) || length(bounds) != 4L) {
    stop(what, ": `bounds` must be c(xmin, xmax, ymin, ymax).", call. = FALSE)
  }
}

#' The log probability of an individual's captures on one occasion.
#'
#' For use inside an observation function of [et_latent_data()], where it is
#' transpiled to ET's `capture_loglik`. Each active detector `j` has a
#' half-normal kernel `k_j = exp(-d_j^2 / (2 sigma^2))`, `d_j` its distance
#' from the centre `(x, y)`. With `form = "hazard"` the individual is caught at
#' least once with probability `1 - exp(-rate * sum(k))`, each capture at
#' detector `j` with probability `k_j / sum(k)`; with `form = "poisson"` the
#' count at detector `j` is Poisson with mean `rate * k_j`.
#'
#' It has no R implementation: calling it outside a transpiled body is an error.
#'
#' @param captures The [et_spatial_captures()] extra, e.g. `data$captures`.
#' @param i,occasion The individual and the occasion.
#' @param x,y The activity centre.
#' @param sigma The kernel scale.
#' @param rate The capture rate.
#' @param detectors The latent detectors' positions, e.g. `data$latents$Xu`.
#' @param form `"hazard"` or `"poisson"`.
#' @return A log probability.
#' @export
et_capture_loglik <- function(captures, i, occasion, x, y, sigma, rate,
                              detectors = NULL, form = "hazard") {
  stop("et_capture_loglik() is only meaningful inside a transpiled observation ",
       "function.", call. = FALSE)
}


#' Is this the individual being resampled?
#'
#' For a rate function. While iFFBS resamples individual `i`, its own share is
#' taken out of the aggregates, so a rate that reads them (a group's alive count
#' in a force of infection) sees the group without it. By default ET puts it back
#' around every rate call, which re-runs the step's shared values for each rate.
#' With `et_data(focal_self_contribution = FALSE)` it does not, and the rate adds
#' the individual back itself:
#'
#' ```r
#' n <- data$aggregates$n_alive[g, t] + et_focal(data, i)
#' ```
#'
#' It has no R implementation: calling it outside a transpiled body is an error.
#'
#' @param data,i The rate's own arguments.
#' @return `TRUE` while `i` is being resampled.
#' @export
et_focal <- function(data, i) {
  stop("et_focal() is only meaningful inside a transpiled rate function.", call. = FALSE)
}

# ---- the data ---------------------------------------------------------------

#' A weighted chain with sampled latent variables.
#'
#' The R mirror of ET's `latent_hmm_data()`. Each individual's hidden chain has
#' weighted transitions ([et_transitions()] with `weights = TRUE`) and is summed
#' out exactly; the latents its weights and observations read are sampled once
#' per Gibbs sweep. Individuals are independent given the parameters and the
#' latents. Fit it with [et_model()] and [et_sample()]: every parameter is
#' sampled by HMC against the marginal likelihood.
#'
#' Timing: the chain starts in `starting_state` before the first timepoint, and
#' at each timepoint `t` takes the weighted step INTO `t`, then the observation
#' at `t`.
#'
#' @param n_individuals,n_timepoints Problem size.
#' @param transitions An [et_transitions()] spec with `weights = TRUE`.
#' @param starting_state `function(model, data, i)` returning the weight of each
#'   state before the first step, as `c(...)`.
#' @param observation_weights `function(model, data, i, t)` returning the weight
#'   of each state's observations at `t` as `c(...)`, or `NULL` when nothing is
#'   observed at `t`.
#' @param latents An [et_latents()] declaration.
#' @param sampling_period Optional two-column matrix (`first`, `last`) per
#'   individual.
#' @param extras Named list of your own data, reachable as `data$name`. An
#'   [et_spatial_captures()] object may be one of them.
#' @param helpers List of [et_helper()] declarations callable from any body.
#' @param observation_shared Optional `function(model, data, i, key)`: whatever
#'   the observation weights at every timepoint with `shared_key(t) == key`
#'   need, computed once per individual and key and handed to
#'   `observation_weights` as its fifth argument, `shared`. It returns a number,
#'   or a named `list()` read as `shared$name`. Worth it for an expensive term
#'   common to several timepoints, such as a year's capture likelihood.
#' @param shared_key,n_keys `function(t)` giving the key of timepoint `t`, and
#'   the number of keys (keys run from 1 to `n_keys`). Left out, every
#'   timepoint is its own key.
#' @param per_evaluation Optional `function(model, data)` returning a named
#'   `list()` of values computed once per likelihood evaluation rather than at
#'   every step (a random walk's cumulative sum, a table of rates, a
#'   reparameterisation). The values are added to `model`, so every other body
#'   reads them as `model$name`. Vectors are allocated in the parameters' number
#'   type (`numeric(n)`), and every number in the list is converted to it unless
#'   it is a whole number or a logical, which stay as they are so they can index.
#' @param cache `TRUE` to reuse observation values an HMC block cannot change:
#'   a block that moves none of the parameters `observation_shared` reads reuses
#'   the shared values, and one that moves none of the parameters the
#'   observations read reuses the observation weights. Which parameters those are
#'   is worked out from the bodies, as for `depends=`, through `per_evaluation`
#'   and the helpers. A named list `list(shared = , weights = , latents = )`
#'   gives them by hand instead (needed only when a body is an [et_julia()]
#'   without `reads`); it is checked against what the bodies read, and refused
#'   if it leaves any out, since a stale cache is silently wrong. `latents`
#'   names the latents `observation_shared` reads (all of them if left out).
#' @return An object of class `et_latent_data`.
#' @export
et_latent_data <- function(n_individuals, n_timepoints, transitions, starting_state,
                           observation_weights = NULL, latents,
                           sampling_period = NULL, extras = list(), helpers = list(),
                           observation_shared = NULL, shared_key = NULL, n_keys = NULL,
                           per_evaluation = NULL, cache = FALSE) {
  if (!inherits(transitions, "et_transitions") || !isTRUE(transitions$weights)) {
    stop("et_latent_data(): `transitions` must come from et_transitions(..., ",
         "weights = TRUE).", call. = FALSE)
  }
  if (!inherits(latents, "et_latents")) {
    stop("et_latent_data(): `latents` must come from et_latents().", call. = FALSE)
  }
  n_individuals <- as_count(n_individuals, "n_individuals")
  n_timepoints <- as_count(n_timepoints, "n_timepoints")
  if (is.null(observation_weights)) {
    message("et_latent_data(): no observation weights given -- the model will ",
            "have no observations.")
  }
  if (!is.list(extras)) stop("et_latent_data(): `extras` must be a list.", call. = FALSE)
  if (length(extras)) {
    if (is.null(names(extras)) || any(!nzchar(names(extras)))) {
      stop("et_latent_data(): every entry of `extras` must be named.", call. = FALSE)
    }
    check_julia_name(names(extras), "extras name")
    clash <- intersect(names(extras), .latent_data_fields)
    if (length(clash)) {
      stop("et_latent_data(): extras name(s) clash with fields ET already ",
           "defines on `data`: ", paste(clash, collapse = ", "), ".", call. = FALSE)
    }
    for (nm in names(extras)) {
      x <- extras[[nm]]
      if (inherits(x, "et_spatial_captures") && x$n_individuals != n_individuals) {
        stop("et_latent_data(): the captures in `", nm, "` are for ",
             x$n_individuals, " individuals, not ", n_individuals, ".", call. = FALSE)
      }
    }
  }
  for (h in helpers) {
    if (!inherits(h, "et_helper")) {
      stop("et_latent_data(): every entry of `helpers` must come from et_helper().",
           call. = FALSE)
    }
  }
  if (!is.null(sampling_period)) {
    sampling_period <- as.matrix(sampling_period)
    if (ncol(sampling_period) != 2L || nrow(sampling_period) != n_individuals) {
      stop("et_latent_data(): `sampling_period` must be an ", n_individuals,
           " x 2 matrix of (first, last) timepoints.", call. = FALSE)
    }
    storage.mode(sampling_period) <- "integer"
  }
  for (nm in names(latents)) check_latent_values(latents[[nm]], nm, n_individuals, extras)
  is_body <- function(f) is.function(f) || inherits(f, "et_julia")
  if (!is.null(observation_shared)) {
    if (!is_body(observation_shared)) {
      stop("et_latent_data(): `observation_shared` must be a function(model, data, ",
           "i, key) or an et_julia().", call. = FALSE)
    }
    if (is.null(observation_weights)) {
      stop("et_latent_data(): `observation_shared` is handed to ",
           "`observation_weights`, so give those too.", call. = FALSE)
    }
  }
  if (!is.null(shared_key) || !is.null(n_keys)) {
    if (is.null(observation_shared)) {
      stop("et_latent_data(): `shared_key` and `n_keys` key the values of ",
           "`observation_shared`; give that too.", call. = FALSE)
    }
    if (is.null(shared_key) || is.null(n_keys)) {
      stop("et_latent_data(): give `shared_key` and `n_keys` together.", call. = FALSE)
    }
    if (!is.function(shared_key) || length(formals(shared_key)) != 1L) {
      stop("et_latent_data(): `shared_key` must be a function of the timepoint.",
           call. = FALSE)
    }
    n_keys <- as_count(n_keys, "n_keys")
  }
  if (!is.null(per_evaluation) && !is_body(per_evaluation)) {
    stop("et_latent_data(): `per_evaluation` must be a function(model, data) or ",
         "an et_julia().", call. = FALSE)
  }
  cache <- check_cache_arg(cache, latents)
  structure(list(
    n_individuals = n_individuals, n_timepoints = n_timepoints,
    transitions = transitions, starting_state = starting_state,
    observation_weights = observation_weights, latents = latents,
    sampling_period = sampling_period, extras = extras, helpers = helpers,
    observation_shared = observation_shared, shared_key = shared_key,
    n_keys = n_keys, per_evaluation = per_evaluation, cache = cache),
    class = "et_latent_data")
}

# FALSE, TRUE (derive the names), or the names by hand. Parameter names are
# checked against the model in codegen, where the parameters are known.
check_cache_arg <- function(cache, latents) {
  if (isFALSE(cache) || is.null(cache)) return(NULL)
  if (isTRUE(cache)) return(list(derive = TRUE))
  if (!is.list(cache) || is.null(names(cache)) ||
      !all(names(cache) %in% c("shared", "weights", "latents"))) {
    stop("et_latent_data(): `cache` must be TRUE, FALSE, or ",
         "list(shared = , weights = , latents = ).", call. = FALSE)
  }
  for (f in names(cache)) {
    if (!is.null(cache[[f]]) && !is.character(cache[[f]])) {
      stop("et_latent_data(): `cache$", f, "` must be a character vector of names.",
           call. = FALSE)
    }
  }
  bad <- setdiff(cache$latents, names(latents))
  if (length(bad)) {
    stop("et_latent_data(): `cache$latents` names unknown latent(s): ",
         paste(bad, collapse = ", "), ".", call. = FALSE)
  }
  c(list(derive = FALSE), cache)
}

# Names LatentHMMData answers itself, so an extra under one would never be read.
.latent_data_fields <- c("n_individuals", "n_timepoints", "n_states", "state_space",
                         "sampling_period", "trans_mat", "from", "to", "vfrom", "vto",
                         "starting_state", "observation", "observation_shared",
                         "shared_key", "n_keys", "latent_decls", "current",
                         "per_evaluation", "cache_decl", "cache", "extras", "latents")

# Shapes are checked here, in R, because a wrong one reaches Julia as an
# indexing convention rather than an error.
check_latent_values <- function(k, nm, N, extras) {
  what <- sprintf("latent '%s'", nm)
  if (k$kind == "per_individual") {
    if (!is.null(k$observed) && length(k$observed) != N) {
      stop(what, ": `observed` has length ", length(k$observed), ", not ", N, ".",
           call. = FALSE)
    }
    if (is.numeric(k$init) && length(k$init) != N) {
      stop(what, ": `init` has length ", length(k$init), ", not ", N, ".", call. = FALSE)
    }
    if (is.function(k$prior) && is.null(k$init) &&
        (is.null(k$observed) || all(is.na(k$observed)))) {
      stop(what, ": its prior depends on the parameters, so give `init` or ",
           "some `observed` values to start from.", call. = FALSE)
    }
    return(invisible(TRUE))
  }
  shape <- if (k$kind == "field") c(k$dim, k$n, N) else c(k$dim, k$n)
  if (is.numeric(k$init) && !identical(as.integer(dim(k$init)), as.integer(shape))) {
    stop(what, ": `init` must be ", paste(shape, collapse = " x "), ", got ",
         paste(dim(k$init) %||% length(k$init), collapse = " x "), ".", call. = FALSE)
  }
  markers <- if (inherits(k$init, "et_starting_detectors")) list(k$init, k$init$centres)
             else list(k$init)
  for (marker in markers) {
    if (inherits(marker, c("et_starting_centres", "et_starting_detectors")) &&
        is.null(extras_name_of(marker$captures, extras))) {
      stop(what, ": the captures its starting values come from must also be ",
           "passed in `extras`.", call. = FALSE)
    }
  }
  if (inherits(k$init, "et_starting_centres") && k$init$n != k$n) {
    stop(what, ": et_starting_centres() has n = ", k$init$n, " but the field has ",
         k$n, " elements.", call. = FALSE)
  }
  if (inherits(k$init, "et_starting_detectors") &&
      length(k$init$captures$latent_detectors) != k$n) {
    stop(what, ": the captures have ", length(k$init$captures$latent_detectors),
         " latent detectors but the latent has ", k$n, ".", call. = FALSE)
  }
  invisible(TRUE)
}

extras_name_of <- function(x, extras) {
  for (nm in names(extras)) if (identical(extras[[nm]], x)) return(nm)
  NULL
}

#' @export
print.et_latent_data <- function(x, ...) {
  cat("<et_latent_data> ", x$n_individuals, " individuals x ", x$n_timepoints,
      " timepoints, weighted chain\n", sep = "")
  cat("  states: ", paste(x$transitions$states, collapse = ", "), "\n", sep = "")
  cat("  latents: ", paste(names(x$latents), collapse = ", "), "\n", sep = "")
  if (length(x$extras)) {
    cat("  extras: ", paste(names(x$extras), collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

is_latent <- function(model) identical(model$likelihood, "latent")
