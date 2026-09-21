# Hand-written reversible summaries, for updates et_aggregate() cannot derive
# ones whose contribution depends on something the `A[i] <- A[i] + x` shape
# does not expose, such as the badger coupling counts, which depend on where
# the individual goes at t+1.
#
# Forward update then reverse must leave the aggregates exactly as they were.
# That invariant is what makes iFFBS correct, and nothing here can check it
# for you.

#' Declare a tracked array with a hand-written reversible update.
#'
#' The escape hatch from [et_aggregate()], for an update whose reverse cannot be
#' derived from the `A[idx] <- A[idx] + contribution` shape.
#'
#' `update` is a function `(model, data, X, state, i, t, reverse)` written in the
#' transpilable R subset. It must apply the contribution when `reverse` is
#' `FALSE` and remove exactly the same contribution when it is `TRUE`. Refer to
#' the array by name; the package allocates it.
#'
#' @param type Julia element type: `"Int"`, `"Float64"` or `"Bool"`.
#' @param dim Integer vector of dimensions.
#' @param update A function `(model, data, X, state, i, t, reverse)`.
#' @return An object of class `et_summary`.
#' @export
#' @examples
#' # How many S individuals in this group went on to E at t -> t+1. The
#' # contribution depends on X[t+1, i], so @aggregate cannot express it.
#' nSE <- et_summary("Int", c(34, 161),
#'   update = function(model, data, X, state, i, t, reverse) {
#'     g <- data$social_group[i, t]
#'     if (g > 0 && t < data$n_timepoints) {
#'       c <- (state == "S") && (X[t + 1, i] == "E")
#'       if (reverse) nSE[g, t] <- nSE[g, t] - c else nSE[g, t] <- nSE[g, t] + c
#'     }
#'   })
et_summary <- function(type = c("Int", "Float64", "Bool"), dim, update) {
  type <- match.arg(type)
  dim <- as.numeric(dim)
  if (!length(dim) || any(is.na(dim)) || any(dim != round(dim)) || any(dim < 0)) {
    stop("et_summary(): `dim` must be non-negative whole numbers.", call. = FALSE)
  }
  if (!is.function(update) && !inherits(update, "et_julia")) {
    stop("et_summary(): `update` must be a function or an et_julia() object.",
         call. = FALSE)
  }
  if (is.function(update)) {
    check_protocol_formals(update,
                           c("model", "data", "X", "state", "i", "t", "reverse"),
                           "summary update", "update")
  }
  structure(list(type = type, dim = as.integer(dim), update = update),
            class = "et_summary")
}

#' @export
print.et_summary <- function(x, ...) {
  cat("<et_summary> ", x$type, " (", paste(x$dim, collapse = ", "),
      ")  hand-written reverse\n", sep = "")
  invisible(x)
}

#' Declare tracked arrays, mixing derived and hand-written updates.
#'
#' Accepts [et_array()] entries (whose reverse the package derives from the
#' update body, exactly as [et_aggregate()] does) alongside [et_summary()]
#' entries (whose reverse you write yourself). Use it when a model needs both.
#'
#' When every array is an [et_array()], prefer [et_aggregate()]: it is the same
#' thing with less ceremony.
#'
#' @param states Character vector of state names, in encoding order.
#' @param arrays Named list of [et_array()] and/or [et_summary()] declarations.
#' @param update A function `(model, data, X, state, i, t)` holding the update
#'   statements for the [et_array()] entries. Omit when there are none.
#' @return An object of class `et_aggregate`, accepted by [et_data()].
#' @export
et_aggregates <- function(states, arrays, update = NULL) {
  states <- check_states(states)
  if (!is.list(arrays) || !length(arrays) || is.null(names(arrays)) ||
      any(!nzchar(names(arrays)))) {
    stop("et_aggregates(): `arrays` must be a non-empty NAMED list.", call. = FALSE)
  }
  check_julia_name(names(arrays), "aggregate array name")

  kinds <- vapply(arrays, function(a) class(a)[1], character(1))
  bad <- names(arrays)[!(kinds %in% c("et_array", "et_summary"))]
  if (length(bad)) {
    stop("et_aggregates(): '", paste(bad, collapse = "', '"),
         "' must be created with et_array() or et_summary().", call. = FALSE)
  }

  derived_names <- names(arrays)[kinds == "et_array"]
  hand_names <- names(arrays)[kinds == "et_summary"]

  lines <- character()
  records <- list()
  reads <- character()
  if (length(derived_names)) {
    if (is.null(update)) {
      stop("et_aggregates(): `update` is required when any array is an ",
           "et_array() -- that is where its update statements go. Arrays needing ",
           "one: ", paste(derived_names, collapse = ", "), ".", call. = FALSE)
    }
    check_protocol_formals(update, c("model", "data", "X", "state", "i", "t"),
                           "aggregate update", "update")
    ctx <- new_ctx(state_syms = TRUE, states = states)
    # Parsed twice, deliberately: `parse_aggregate_body` validates (it is the
    # thing that refuses a non-reversible shape) and renders the `@aggregate`
    # lines; `aggregate_records` returns the same updates in pieces, so the
    # fallback path can emit them as standalone functions.
    lines <- parse_aggregate_body(body(update), derived_names, ctx)
    records <- aggregate_records(body(update), derived_names, ctx)
    reads <- ctx$reads$model
  } else if (!is.null(update)) {
    stop("et_aggregates(): `update` was given but every array is an ",
         "et_summary(), which carries its own update. Drop `update`, or declare ",
         "the arrays it writes with et_array().", call. = FALSE)
  }

  structure(list(states = states, arrays = arrays, lines = lines,
                 records = records, reads = reads,
                 derived_names = derived_names, hand_names = hand_names),
            class = "et_aggregate")
}

# Transpile one hand-written summary into the Julia function ET calls.
#
# ET passes `reverse` positionally with a default, and a keyword on a call the
# compiler cannot resolve forces the kwarg path and allocates a NamedTuple per
# call, which, in a function invoked per individual per timepoint, is a real
# cost. So the emitted signature takes it positionally.
et_transpile_summary <- function(name, spec, states, helpers = list(),
                                 array_names = character()) {
  fname <- paste0("et_summary_", name)
  if (inherits(spec$update, "et_julia")) {
    return(list(src = sprintf(
      "function %s(model, data, X, state, i, t, reverse=false)\n%s\nend",
      fname, indent(spec$update$src)),
      name = fname, reads = spec$update$reads %||% character(),
      opaque = is.null(spec$update$reads)))
  }
  helper_names <- vapply(helpers, function(h) h$name, character(1))
  ctx <- new_ctx(helpers = helper_names, state_syms = TRUE, states = states)
  # A tracked array is written as a bare name (the vocabulary `@aggregate`
  # establishes); ET sees this function directly, so the qualification the macro
  # would have done has to happen here.
  body_src <- qualify_aggregates(transpile_braced_body(body(spec$update), ctx),
                                 array_names)
  list(src = sprintf(
         "function %s(model, data, X, state, i, t, reverse=false)\n%s\n    nothing\nend",
         fname, indent(body_src)),
       name = fname, reads = ctx$reads$model, opaque = FALSE)
}
