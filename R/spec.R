# The model-structure DSL: states, aggregates, transitions, survival. One R
# constructor per ET concept, same names and argument order, so a model reads
# across from the ET documentation.

# ---- aggregates -------------------------------------------------------------

#' Declare a tracked array for [et_aggregate()].
#'
#' @param type Julia element type: `"Int"`, `"Float64"` or `"Bool"`.
#' @param dim Integer vector of dimensions.
#' @return An object of class `et_array`.
#' @export
et_array <- function(type = c("Int", "Float64", "Bool"), dim) {
  type <- match.arg(type)
  dim <- as.numeric(dim)
  if (!length(dim) || any(is.na(dim)) || any(dim != round(dim)) || any(dim < 0)) {
    stop("et_array(): `dim` must be non-negative whole numbers.", call. = FALSE)
  }
  structure(list(type = type, dim = as.integer(dim)), class = "et_array")
}

#' @export
print.et_array <- function(x, ...) {
  cat("<et_array> ", x$type, " (", paste(x$dim, collapse = ", "), ")\n", sep = "")
  invisible(x)
}

#' Declare the arrays tracked during the latent update, and how they update.
#'
#' This is ET's central design rule made available from R: the package attaches
#' no meaning to these arrays, it only runs the update forwards and -- when the
#' latent sampler removes an individual's contribution in reverse. That
#' reversibility is what makes iFFBS both correct and cheap, so the `update` body
#' is restricted to statements from which a reverse can be derived.
#'
#' Each statement must be
#' `A[idx...] <- A[idx...] + contribution` (or `*`), optionally wrapped in
#' `if (cond) { ... }`. `model`, `data`, `X`, `state`, `i` and `t` are in scope;
#' compare states by name, e.g. `state == "I"`.
#'
#' @param states Character vector of state names, in the order that fixes their
#'   encoding.
#' @param arrays Named list of [et_array()] declarations.
#' @param update A function `(model, data, X, state, i, t)` whose body holds the
#'   update statements.
#' @return An object of class `et_aggregate`.
#' @export
#' @examples
#' aggs <- et_aggregate(
#'   c("S", "I"),
#'   arrays = list(n_infected = et_array("Int", c(10, 80))),
#'   update = function(model, data, X, state, i, t) {
#'     n_infected[data$group[i], t] <- n_infected[data$group[i], t] + (state == "I")
#'   })
et_aggregate <- function(states, arrays, update) {
  states <- check_states(states)
  if (!is.list(arrays) || !length(arrays) || is.null(names(arrays)) ||
      any(!nzchar(names(arrays)))) {
    stop("et_aggregate(): `arrays` must be a non-empty NAMED list of et_array().",
         call. = FALSE)
  }
  check_julia_name(names(arrays), "aggregate array name")
  for (nm in names(arrays)) {
    if (!inherits(arrays[[nm]], "et_array")) {
      stop("et_aggregate(): array '", nm, "' must be created with et_array().",
           call. = FALSE)
    }
  }
  if (!is.function(update)) {
    stop("et_aggregate(): `update` must be a function.", call. = FALSE)
  }
  check_protocol_formals(update, c("model", "data", "X", "state", "i", "t"),
                         "aggregate update", "update")

  ctx <- new_ctx(state_syms = TRUE, states = states)
  lines <- parse_aggregate_body(body(update), names(arrays), ctx)
  structure(list(states = states, arrays = arrays, lines = lines,
                 reads = ctx$reads$model), class = "et_aggregate")
}

#' @export
print.et_aggregate <- function(x, ...) {
  cat("<et_aggregate> states: ", paste(x$states, collapse = ", "), "\n", sep = "")
  for (nm in names(x$arrays)) {
    cat("  @array ", nm, " ", x$arrays[[nm]]$type, " (",
        paste(x$arrays[[nm]]$dim, collapse = ", "), ")\n", sep = "")
  }
  cat(paste0("  ", unlist(strsplit(paste(x$lines, collapse = "\n"), "\n")),
             collapse = "\n"), "\n", sep = "")
  invisible(x)
}

# Walk the update body, recognising only reversible shapes. Anything else is an
# error: guessing a reverse would silently break the aggregates-agree-with-X
# invariant that everything downstream relies on.
parse_aggregate_body <- function(node, array_names, ctx, guard = NULL) {
  stmts <- if (is.call(node) && identical(as.character(node[[1]]), "{")) {
    as.list(node)[-1L]
  } else {
    list(node)
  }
  unlist(lapply(stmts, parse_aggregate_stmt, array_names = array_names,
                ctx = ctx, guard = guard))
}

# The same walk, but returning one record per update instead of rendered text:
# `array`, `target`, `op`, `contrib` and the enclosing `guard` (or NULL). This is
# what lets the same declaration be emitted either as an `@aggregate` line or as
# a standalone reversible function, without re-parsing or rewriting strings.
aggregate_records <- function(node, array_names, ctx, guard = NULL) {
  stmts <- if (is.call(node) && identical(as.character(node[[1]]), "{")) {
    as.list(node)[-1L]
  } else {
    list(node)
  }
  out <- list()
  for (e in stmts) {
    op <- if (is.call(e)) as.character(e[[1]]) else ""
    if (op == "if") {
      cond <- et_transpile_expr(e[[2]], ctx)
      inner_guard <- if (is.null(guard)) cond else paste0("(", guard, ") && (", cond, ")")
      out <- c(out, aggregate_records(e[[3]], array_names, ctx, inner_guard))
      next
    }
    lhs <- e[[2]]; rhs <- e[[3]]
    out[[length(out) + 1L]] <- list(
      array   = as.character(lhs[[2]]),
      target  = et_transpile_expr(lhs, ctx),
      op      = as.character(rhs[[1]]),
      contrib = et_transpile_expr(rhs[[3]], ctx),
      guard   = guard)
  }
  out
}

parse_aggregate_stmt <- function(e, array_names, ctx, guard = NULL) {
  if (!is.call(e)) stop(agg_error(e), call. = FALSE)
  op <- as.character(e[[1]])

  if (op == "if") {
    if (length(e) == 4L) {
      stop("et_aggregate(): `else` is not supported in an update body. Write two ",
           "guarded statements instead -- each must be reversible on its own.",
           call. = FALSE)
    }
    cond <- et_transpile_expr(e[[2]], ctx)
    inner <- parse_aggregate_body(e[[3]], array_names, ctx, guard)
    return(paste0("if ", cond, "\n", indent(paste(inner, collapse = "\n")), "\nend"))
  }

  if (!(op %in% c("<-", "="))) stop(agg_error(e), call. = FALSE)
  lhs <- e[[2]]; rhs <- e[[3]]
  if (!is.call(lhs) || !identical(as.character(lhs[[1]]), "[")) {
    stop("et_aggregate(): the left-hand side of an update must index a declared ",
         "array, e.g. `n_infected[g, t] <- ...`. Got: ", deparse(lhs), call. = FALSE)
  }
  arr <- as.character(lhs[[2]])
  if (!(arr %in% array_names)) {
    stop("et_aggregate(): '", arr, "' is not a declared array. Declared: ",
         paste(array_names, collapse = ", "), ".", call. = FALSE)
  }
  if (!is.call(rhs) || !(as.character(rhs[[1]]) %in% c("+", "*")) ||
      length(rhs) != 3L || !identical(rhs[[2]], lhs)) {
    stop("et_aggregate(): each update must read\n",
         "    ", deparse(lhs), " <- ", deparse(lhs), " + <contribution>\n",
         "  (or `*`). The package derives the REVERSE update from that shape, ",
         "which is what lets the latent sampler take an individual's ",
         "contribution back out. Got: ", deparse(e), call. = FALSE)
  }
  sprintf("%s %s= %s", et_transpile_expr(lhs, ctx), as.character(rhs[[1]]),
          et_transpile_expr(rhs[[3]], ctx))
}

agg_error <- function(e) {
  paste0("et_aggregate(): unsupported statement in the update body: ",
         deparse(e)[1], "\n  Allowed: `A[idx] <- A[idx] + expr`, `A[idx] <- ",
         "A[idx] * expr`, and `if (cond) { ... }` around them.")
}

# ---- transitions ------------------------------------------------------------

#' Declare one transition.
#' @param from,to State names.
#' @param rate A rate function `(model, data, i, t)` or an [et_julia()].
#' @return An object of class `et_trans`.
#' @export
et_trans <- function(from, to, rate) {
  structure(list(from = from, to = to, rate = rate), class = "et_trans")
}

#' Declare that every step is conditional on survival.
#'
#' Mirrors ET's `@survival`: each non-death transition's rate is scaled by
#' `fn`, and **every** live state gains a transition to `death` with the
#' remaining mass: including states that never appear as the source of a
#' declared transition.
#'
#' @param fn A survival function `(model, data, i, t)` returning P(survive the
#'   `t` to `t+1` step), or an [et_julia()].
#' @param death The absorbing state's name.
#' @return An object of class `et_survival`.
#' @export
et_survival <- function(fn, death) {
  stopifnot(is.character(death), length(death) == 1L)
  structure(list(fn = fn, death = death), class = "et_survival")
}

#' Declare a model's states, transitions and rates.
#'
#' @param states Character vector of state names. The order fixes the encoding:
#'   state `k` in the trajectory means `states[k]`.
#' @param ... Transitions, either as named arguments whose name is an arrow
#'   string (`"S -> I" = infection`) or as [et_trans()] objects. The arrow form
#'   mirrors the Julia line character for character.
#' @param survival Optional [et_survival()].
#' @param auto_self Fill each state's self-transition with the leftover
#'   probability mass (ET's `:auto_self`, on by default).
#' @return An object of class `et_transitions`.
#' @export
#' @examples
#' infection <- function(model, data, i, t) -expm1(-model$alpha)
#' recovery  <- function(model, data, i, t) 1 / model$m
#' et_transitions(c("S", "I"), "S -> I" = infection, "I -> S" = recovery)
et_transitions <- function(states, ..., survival = NULL, auto_self = TRUE) {
  states <- check_states(states)
  args <- list(...)
  nms <- names(args) %||% rep("", length(args))
  trans <- list()
  for (k in seq_along(args)) {
    a <- args[[k]]
    if (inherits(a, "et_trans")) {
      trans[[length(trans) + 1L]] <- a
      next
    }
    if (!nzchar(nms[k])) {
      stop("et_transitions(): argument ", k, " is neither an et_trans() object ",
           'nor a named arrow argument such as `"S -> I" = infection`.',
           call. = FALSE)
    }
    m <- regmatches(nms[k], regexec("^\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*->\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*$", nms[k]))[[1]]
    if (!length(m)) {
      stop("et_transitions(): '", nms[k], '\' is not a transition. Use the form ',
           '"From -> To" = rate.', call. = FALSE)
    }
    trans[[length(trans) + 1L]] <- et_trans(m[2], m[3], a)
  }
  if (!length(trans)) {
    stop("et_transitions(): declare at least one transition.", call. = FALSE)
  }
  for (tr in trans) {
    bad <- missing_from(c(tr$from, tr$to), states)
    if (length(bad)) {
      stop("et_transitions(): transition ", tr$from, " -> ", tr$to,
           " names state(s) not in the state space: ",
           paste(bad, collapse = ", "), ". States are: ",
           paste(states, collapse = ", "), ".", call. = FALSE)
    }
  }
  if (!is.null(survival)) {
    if (!inherits(survival, "et_survival")) {
      stop("et_transitions(): `survival` must be created with et_survival().",
           call. = FALSE)
    }
    if (!(survival$death %in% states)) {
      stop("et_survival(): death state '", survival$death,
           "' is not in the state space.", call. = FALSE)
    }
  }
  key <- vapply(trans, function(tr) paste0(tr$from, "->", tr$to), character(1))
  if (anyDuplicated(key)) {
    stop("et_transitions(): duplicate transition(s): ",
         paste(unique(key[duplicated(key)]), collapse = ", "), ".", call. = FALSE)
  }
  structure(list(states = states, transitions = trans, survival = survival,
                 auto_self = isTRUE(auto_self)), class = "et_transitions")
}

#' @export
print.et_transitions <- function(x, ...) {
  cat("<et_transitions> states: ", paste(x$states, collapse = ", "), "\n", sep = "")
  if (!is.null(x$survival)) {
    cat("  @survival  death = ", x$survival$death, "\n", sep = "")
  }
  for (tr in x$transitions) cat("  ", tr$from, " -> ", tr$to, "\n", sep = "")
  invisible(x)
}

# ---- shared validation ------------------------------------------------------

check_states <- function(states) {
  if (!is.character(states) || length(states) < 2L) {
    stop("states must be a character vector of at least two state names.",
         call. = FALSE)
  }
  if (anyDuplicated(states)) {
    stop("duplicate state name(s): ",
         paste(unique(states[duplicated(states)]), collapse = ", "), ".",
         call. = FALSE)
  }
  check_julia_name(states, "state name")
  states
}
