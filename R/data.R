# et_data(): the R mirror of ET's epidemic_data(), argument for argument.
# `extras` is ET's `extras...`, reachable as `data$name` in a transpiled body
# and never looked inside.

#' Build the model's fixed structure and tracked state.
#'
#' @param n_individuals,n_timepoints Problem size.
#' @param transitions An [et_transitions()] spec.
#' @param starting_state A function `(model, data, X, i, t)` returning the
#'   probability vector over states at the individual's first timepoint.
#' @param aggregates An [et_aggregate()] declaration.
#' @param group Integer group index per individual. Used only by the default
#'   `affected_individuals` and by whatever your own functions read off it.
#' @param observation_process A function `(model, data, X, i, t)` returning a
#'   per-state weight vector.
#' @param observation_weight A function `(model, data, X, i, t, s)` returning the
#'   scalar weight for state `s`. **This is the fast form**: the likelihood needs
#'   only one entry, so the scalar avoids allocating a weight vector per cell
#'   under AD. Supply either or both; see ET's `epidemic_data` docs for when both
#'   are needed.
#' @param likelihood_weight A function `(model, data, X, i, t, s)` giving the
#'   observation factor scored in the LIKELIHOOD, when that differs from the
#'   filter's `observation_weight`. This is ET's seam for keeping some
#'   observation parameters conjugate: a process that factorises
#'   multiplicatively (`w = capture x tests`) gives the product as
#'   `observation_weight`, so the filter sees everything, and only the
#'   non-conjugate FACTOR here. Because the weights multiply, the log-likelihood
#'   is a sum of the factors' contributions, so dropping one drops exactly its
#'   term. **Keeping a factor in both this and a conjugate block double-counts
#'   it** - the conjugate parameter must not appear in `likelihood_weight`.
#'   Defaults to `observation_weight`.
#' @param sampling_period A two-column matrix (`first`, `last`) per individual.
#'   Defaults to the whole window for everyone.
#' @param affected_individuals Who each individual's state affects. `NULL` uses
#'   groupmates under a fixed `group`; pass [et_affected_from_groups()] for
#'   time-varying membership, or an [et_julia()] expression for anything else.
#' @param coupled_transitions List of `c(from, to)` pairs naming which of a
#'   neighbour's transitions this individual can influence, purely an
#'   optimisation, and often a large one.
#' @param coupling_transitions An optional separate [et_transitions()] used for
#'   the coupling term only (ET's `coupling_trans_mat`). Power users only.
#' @param rest_contribution An [et_julia()] replacing the whole coupling term.
#' @param focal_self_contribution See ET's `epidemic_data`; leave `TRUE` unless
#'   your rate functions already do the focal's own accounting.
#' @param extras Named list of your own arrays, reachable as `data$name`.
#' @param helpers List of [et_helper()] declarations callable from any body.
#' @return An object of class `et_data`.
#' @export
et_data <- function(n_individuals, n_timepoints, transitions, starting_state,
                    aggregates, group = NULL,
                    observation_process = NULL, observation_weight = NULL,
                    likelihood_weight = NULL, sampling_period = NULL, affected_individuals = NULL,
                    coupled_transitions = NULL, coupling_transitions = NULL,
                    rest_contribution = NULL, focal_self_contribution = TRUE,
                    extras = list(), helpers = list()) {
  if (!inherits(transitions, "et_transitions")) {
    stop("et_data(): `transitions` must come from et_transitions().", call. = FALSE)
  }
  if (!inherits(aggregates, "et_aggregate")) {
    stop("et_data(): `aggregates` must come from et_aggregate() or et_aggregates().",
         call. = FALSE)
  }
  if (!identical(transitions$states, aggregates$states)) {
    stop("et_data(): the state space differs between `transitions` (",
         paste(transitions$states, collapse = ", "), ") and `aggregates` (",
         paste(aggregates$states, collapse = ", "), "). They must be identical, ",
         "including order -- the order is what fixes the encoding in the ",
         "trajectory.", call. = FALSE)
  }
  n_individuals <- as_count(n_individuals, "n_individuals")
  n_timepoints <- as_count(n_timepoints, "n_timepoints")

  if (!is.null(likelihood_weight) && is.null(observation_weight) &&
      is.null(observation_process)) {
    stop("et_data(): `likelihood_weight` replaces the likelihood's factor of an ",
         "observation process that the FILTER still needs in full. Supply ",
         "`observation_weight` (or `observation_process`) as well.", call. = FALSE)
  }
  if (is.null(observation_process) && is.null(observation_weight)) {
    message("et_data(): no observation process given -- the model will have no ",
            "observation likelihood. This is ET's honest default, but if you do ",
            "have data it is almost certainly a mistake.")
  }
  if (!is.list(extras)) stop("et_data(): `extras` must be a list.", call. = FALSE)
  if (length(extras)) {
    if (is.null(names(extras)) || any(!nzchar(names(extras)))) {
      stop("et_data(): every entry of `extras` must be named.", call. = FALSE)
    }
    check_julia_name(names(extras), "extras name")
    reserved <- c("aggregates", "state_space", "n_states", "n_individuals",
                  "n_timepoints", "sampling_period", "trans_mat", "group")
    clash <- intersect(names(extras), reserved)
    if (length(clash)) {
      stop("et_data(): extras name(s) clash with fields ET already defines on ",
           "`data`: ", paste(clash, collapse = ", "), ".", call. = FALSE)
    }
  }
  if (!is.null(group)) {
    group <- as.integer(group)
    if (length(group) != n_individuals) {
      stop("et_data(): `group` has length ", length(group), " but there are ",
           n_individuals, " individuals.", call. = FALSE)
    }
  }
  if (!is.null(sampling_period)) {
    sampling_period <- as.matrix(sampling_period)
    if (ncol(sampling_period) != 2L || nrow(sampling_period) != n_individuals) {
      stop("et_data(): `sampling_period` must be an ", n_individuals,
           " x 2 matrix of (first, last) timepoints.", call. = FALSE)
    }
    storage.mode(sampling_period) <- "integer"
  }
  if (!is.null(coupled_transitions)) {
    if (!is.list(coupled_transitions) ||
        !all(vapply(coupled_transitions, function(p) is.character(p) && length(p) == 2L,
                    logical(1)))) {
      stop("et_data(): `coupled_transitions` must be a list of c(from, to) ",
           "character pairs, e.g. list(c(\"S\", \"E\")).", call. = FALSE)
    }
    bad <- missing_from(unlist(coupled_transitions), transitions$states)
    if (length(bad)) {
      stop("et_data(): coupled_transitions names unknown state(s): ",
           paste(bad, collapse = ", "), ".", call. = FALSE)
    }
  }
  for (h in helpers) {
    if (!inherits(h, "et_helper")) {
      stop("et_data(): every entry of `helpers` must come from et_helper().",
           call. = FALSE)
    }
  }

  structure(list(
    n_individuals = n_individuals, n_timepoints = n_timepoints,
    transitions = transitions, starting_state = starting_state,
    aggregates = aggregates, group = group,
    observation_process = observation_process,
    observation_weight = observation_weight,
    likelihood_weight = likelihood_weight,
    sampling_period = sampling_period,
    affected_individuals = affected_individuals,
    coupled_transitions = coupled_transitions,
    coupling_transitions = coupling_transitions,
    rest_contribution = rest_contribution,
    focal_self_contribution = isTRUE(focal_self_contribution),
    extras = extras, helpers = helpers), class = "et_data")
}

#' @export
print.et_data <- function(x, ...) {
  cat("<et_data> ", x$n_individuals, " individuals x ", x$n_timepoints,
      " timepoints\n", sep = "")
  cat("  states: ", paste(x$transitions$states, collapse = ", "), "\n", sep = "")
  cat("  observation: ",
      paste(c(if (!is.null(x$observation_weight)) "weight (scalar)",
              if (!is.null(x$observation_process)) "process (vector)",
              if (is.null(x$observation_weight) && is.null(x$observation_process)) "none"),
            collapse = " + "), "\n", sep = "")
  if (length(x$extras)) {
    cat("  extras: ", paste(names(x$extras), collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

#' Build time-varying `affected_individuals` from a group-membership matrix.
#'
#' Individuals affect their groupmates, with membership allowed to change over
#' time. The focal is excluded from its own list, matching ET's default.
#'
#' @param group_matrix An `n_individuals x n_timepoints` integer matrix of group
#'   indices; `0` means "not present at that time".
#' @return An object of class `et_affected`, accepted by [et_data()].
#' @export
et_affected_from_groups <- function(group_matrix) {
  group_matrix <- as.matrix(group_matrix)
  storage.mode(group_matrix) <- "integer"
  structure(list(group_matrix = group_matrix), class = "et_affected")
}

as_count <- function(x, what) {
  x <- as.numeric(x)
  if (length(x) != 1L || is.na(x) || x != round(x) || x < 1) {
    stop("et_data(): `", what, "` must be a single positive whole number.",
         call. = FALSE)
  }
  as.integer(x)
}
