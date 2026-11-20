# Wraps a whole R function in the Julia signature ET will call it by, checks
# the formals against that protocol, and emits the element-type prologue for
# the roles returning a per-state vector.

.et_protocols <- list(
  rate           = list(args = c("model", "data", "i", "t"),           vector = FALSE),
  survival       = list(args = c("model", "data", "i", "t"),           vector = FALSE),
  starting_state = list(args = c("model", "data", "X", "i", "t"),      vector = TRUE),
  obs_process    = list(args = c("model", "data", "X", "i", "t"),      vector = TRUE),
  obs_weight     = list(args = c("model", "data", "X", "i", "t", "s"), vector = FALSE),
  # Quantities computed once per (i, t) and handed to every rate, or every
  # state's observation weight, as `shared`. It returns a list(), so its return
  # is left as written rather than converted to the parameter scalar type.
  shared         = list(args = c("model", "data", "i", "t"),           vector = FALSE,
                        keep_return = TRUE),
  rate_shared    = list(args = c("model", "data", "i", "t", "shared"), vector = FALSE),
  survival_shared = list(args = c("model", "data", "i", "t", "shared"), vector = FALSE),
  # A coupling term in R: the neighbours' moves as a weight per state of i.
  rest           = list(args = c("model", "data", "X", "i", "t", "n_states"), vector = TRUE),
  obs_weight_shared = list(args = c("model", "data", "X", "i", "t", "s", "shared"),
                           vector = FALSE),
  # A weighted chain with sampled latents (et_latent_data()): the weights over
  # the states before the first step, and each state's observation weight at
  # `t`, both returned as `c(...)` (a tuple) or NULL for no observation. The
  # observation takes `shared` only when an observation_shared is declared, but
  # is always called with it.
  latent_start   = list(args = c("model", "data", "i"),                vector = TRUE,
                        tuple = TRUE),
  latent_obs     = list(args = c("model", "data", "i", "t", "shared"), vector = TRUE,
                        tuple = TRUE, alt_args = c("model", "data", "i", "t")),
  # What the observations at every timepoint of one key need, once per
  # individual and key: a number (converted to the parameter type, as a scalar
  # role's return is) or a list() (each field converted, as a shared step's).
  latent_obs_shared = list(args = c("model", "data", "i", "key"), vector = FALSE,
                           keep_return = TRUE, wrap_scalar = TRUE),
  # Values computed once per likelihood evaluation and merged into `model`. Its
  # list() keeps whole numbers and logicals as they are, so they can index, and
  # its vectors are allocated in the parameter type.
  per_evaluation = list(args = c("model", "data"), vector = FALSE,
                        keep_return = TRUE, settle_list = TRUE, typed_alloc = TRUE),
  # An extra term of the log density, over the parameters alone.
  log_density    = list(args = "model", vector = FALSE)
)

#' The call protocols a user function can be transpiled into.
#' @return A character vector of role names.
#' @export
et_protocols <- function() names(.et_protocols)

# ---- the escape hatch -------------------------------------------------------

#' Wrap literal Julia source.
#'
#' Accepted anywhere a transpiled R function is. Use it for a construct the R
#' subset does not cover, or to hand-optimise one function without abandoning the
#' R model.
#'
#' The body cannot be scanned for `model$x` reads, so any likelihood term
#' containing an `et_julia()` function is emitted **without** a `depends=`
#' annotation. That is the safe direction (PracticalBayes then evaluates the term
#' unconditionally, exactly as it would with no annotation at all), but it does
#' forgo the speed-up the annotation gives. Pass `reads` to declare
#' the parameters the body touches and keep the annotation.
#'
#' @param src Julia source for the function BODY (not the signature, the
#'   protocol's signature is generated around it).
#' @param reads Optional character vector of parameter names the body reads. When
#'   given, `depends=` derivation continues to work; it is then the caller's
#'   responsibility to be exhaustive.
#' @return An object of class `et_julia`.
#' @export
et_julia <- function(src, reads = NULL) {
  stopifnot(is.character(src), length(src) == 1L)
  structure(list(src = src, reads = reads), class = "et_julia")
}

#' @export
print.et_julia <- function(x, ...) {
  cat("<et_julia>\n"); cat(x$src, "\n")
  if (!is.null(x$reads)) cat("reads:", paste(x$reads, collapse = ", "), "\n")
  invisible(x)
}

#' Declare a helper function callable from other transpiled bodies.
#'
#' Helpers keep their own argument names and have no protocol. They are emitted
#' before the leaf functions that call them, and their `model$x` reads propagate
#' to every caller (which is what makes `depends=` derivation transitive).
#'
#' @param f An R function, or an [et_julia()] object.
#' @param name Julia name to emit. Defaults to the deparsed argument.
#' @param args Argument names, required only when `f` is an [et_julia()].
#' @return An object of class `et_helper`.
#' @export
et_helper <- function(f, name = deparse(substitute(f)), args = NULL) {
  check_julia_name(name, "helper name")
  if (inherits(f, "et_julia")) {
    if (is.null(args)) {
      stop("et_helper(): `args` is required when the helper is et_julia().",
           call. = FALSE)
    }
  } else if (!is.function(f)) {
    stop("et_helper(): `f` must be an R function or an et_julia() object.",
         call. = FALSE)
  } else {
    args <- names(formals(f))
  }
  structure(list(f = f, name = name, args = args), class = "et_helper")
}

#' @export
print.et_helper <- function(x, ...) {
  cat("<et_helper> ", x$name, "(", paste(x$args, collapse = ", "), ")\n", sep = "")
  invisible(x)
}

# ---- function-level transpilation -------------------------------------------

# Transpile one function into a named Julia definition for `role`.
#
# Returns a list:
#   src    - the Julia `function ... end` text
#   reads  - the model parameter names the body reads (transitively via helpers)
#   opaque - TRUE if the body could not be scanned (an et_julia without `reads`)
#
# `extra_reads` are parameters the body depends on without naming them: for a
# rate that reads `shared`, whatever the shared step function reads. They enter
# the element-type prologue and the returned reads alike.
et_transpile_role <- function(f, role, name, helpers = list(),
                              extra_reads = character(), n_states = NULL,
                              float_floor = FALSE) {
  spec <- .et_protocols[[role]]
  if (is.null(spec)) {
    stop("unknown role '", role, "'. One of: ",
         paste(et_protocols(), collapse = ", "), ".", call. = FALSE)
  }
  check_julia_name(name, "function name")
  helper_names <- vapply(helpers, function(h) h$name, character(1))

  if (inherits(f, "et_julia")) {
    src <- sprintf("function %s(%s)\n%s\nend", name,
                   paste(spec$args, collapse = ", "), indent(f$src))
    return(list(src = src, reads = union(f$reads %||% character(), extra_reads),
                opaque = is.null(f$reads), helpers = character()))
  }
  if (!is.function(f)) {
    stop("the ", role, " '", name, "' must be an R function or an et_julia() ",
         "object; got ", class(f)[1], ".", call. = FALSE)
  }
  if (is.null(spec$alt_args) || !identical(names(formals(f)), spec$alt_args)) {
    check_protocol_formals(f, spec$args, role, name)
  }
  if (isTRUE(spec$tuple) && is.null(n_states)) {
    stop("internal: the ", role, " role needs the number of states.", call. = FALSE)
  }

  ctx <- new_ctx(helpers = helper_names,
                 vector_result = spec$vector || isTRUE(spec$typed_alloc),
                 typed_list = isTRUE(spec$keep_return) && !isTRUE(spec$settle_list),
                 settle_list = isTRUE(spec$settle_list),
                 tuple_c = if (isTRUE(spec$tuple)) n_states)
  # A scalar role must return one concrete type. An R body naturally mixes
  # `return(1)` (an Int) with `return(1 - eta)` (whatever the parameters are),
  # which makes the Julia return type a Union, runtime dispatch in the hottest
  # function in the package. Converting every returned value to `ETR_T` fixes it
  # without the user having to think about it, and `ETR_T` is the promotion of
  # exactly the parameters this body reads.
  fbody <- if (spec$vector) body(f)
           else if (isTRUE(spec$wrap_scalar)) wrap_returns_in_T(body(f), skip_list = TRUE)
           else if (isTRUE(spec$keep_return)) body(f)
           else wrap_returns_in_T(body(f))
  body_src <- transpile_braced_body(fbody, ctx)
  reads <- ctx$reads$model
  used_helpers <- ctx$reads$helpers

  # Reads reached through helpers count too, that transitivity is exactly what
  # makes an automatically derived `depends=` trustworthy.
  reads <- union(union(reads, helper_reads(used_helpers, helpers)), extra_reads)
  opaque <- any(vapply(helpers[helper_names %in% helper_closure(used_helpers, helpers)],
                       function(h) inherits(h$f, "et_julia") && is.null(h$f$reads),
                       logical(1)))

  prologue <- eltype_prologue(reads, float_floor)
  full <- if (length(prologue)) paste0(prologue, "\n", body_src) else body_src
  list(src = sprintf("function %s(%s)\n%s\nend", name,
                     paste(spec$args, collapse = ", "), indent(full)),
       reads = reads, opaque = opaque, helpers = used_helpers)
}

# The name of the parameter-scalar-type variable the generator binds. not `T`:
# R spells `TRUE` as `T`, so a body containing a bare `T` would be ambiguous, and
# a user local named `T` would shadow the binding.
ET_TYPE_VAR <- "ETR_T"

# `ETR_T` is the promotion of the element types of every parameter the body
# reads. It is not "the Dual type": what the parameters are made of depends on
# the autodiff backend, and only FORWARD mode makes them Duals --
#
#   ForwardDiff / PolyesterForwardDiff : Dual
#   Mooncake, Enzyme                   : Float64 (the tangent is carried apart
#                                        from the primal, which stays plain)
#   ReverseDiff                        : a tracked scalar type
#
# Reading it off `eltype(model.x)` is the only formulation right for all of them,
# and it is the idiom ET's own code uses. It is also what makes an allocation
# like `rep(1, n)` follow the BLOCKING rather than the model: the same function
# is called with plain Float64 parameters when they are sampled conjugately, and
# with the backend's differentiable type when they sit in an HMC block.
# That distinction is a real bug source.
#
# Everything the generator emits against it, `convert(ETR_T, x)`, `one(ETR_T)`,
# `zeros(ETR_T, n)`: is defined across all three families. A type CONSTRUCTOR
# call would not be.
#
# `float_floor` adds Float64 to the promotion. Parameters are always floating
# point, but a value `per_evaluation` adds to `model` may be a whole number or a
# logical (a grid index, a flag), and a body reading only those would otherwise
# get an integer `ETR_T` and fail to convert its first fraction.
eltype_prologue <- function(reads, float_floor = FALSE) {
  if (!length(reads)) return(paste0(ET_TYPE_VAR, " = Float64"))
  sprintf("%s = promote_type(%s%s)", ET_TYPE_VAR, if (float_floor) "Float64, " else "",
          paste(sprintf("eltype(model.%s)", reads), collapse = ", "))
}

# AST transform: every value the function can return is wrapped in `T(...)`.
#   `return(x)`        -> `return(T(x))`
#   `{ ...; last }`    -> `{ ...; T(last) }`  (unless the last is a return)
#   `if (c) a else b`  -> branches wrapped individually, so the emitted Julia
#                         keeps its readable multi-line `if` form
# A `for` loop or an assignment as the terminal statement returns nothing
# meaningful, so it is left alone.
#
# With `skip_list`, a returned `list(...)` is left alone: its fields are
# converted one by one instead (a role that may return a number or a list).
wrap_returns_in_T <- function(node, skip_list = FALSE) {
  Tcall <- function(x) {
    if (skip_list && is.call(x) && identical(x[[1]], as.symbol("list"))) return(x)
    as.call(list(as.symbol(ET_TYPE_VAR), x))
  }
  if (!is.call(node)) return(Tcall(node))
  head <- as.character(node[[1]])
  if (head == "return" && length(node) == 2L) {
    node[[2]] <- Tcall(node[[2]])
    return(node)
  }
  if (head == "{") {
    n <- length(node)
    if (n >= 2L) {
      # Interior statements may still hold early `return`s.
      for (k in 2:(n - 1)) node[[k]] <- wrap_inner_returns(node[[k]], skip_list)
      node[[n]] <- wrap_returns_in_T(node[[n]], skip_list)
    }
    return(node)
  }
  if (head == "if") {
    node[[3]] <- wrap_returns_in_T(node[[3]], skip_list)
    if (length(node) == 4L) node[[4]] <- wrap_returns_in_T(node[[4]], skip_list)
    return(node)
  }
  if (head %in% c("for", "while", "<-", "=")) return(node)
  Tcall(node)
}

# Wrap `return(x)` anywhere inside a statement, without touching its value.
wrap_inner_returns <- function(node, skip_list = FALSE) {
  if (!is.call(node)) return(node)
  head <- as.character(node[[1]])
  if (head == "return" && length(node) == 2L) {
    v <- node[[2]]
    if (!(skip_list && is.call(v) && identical(v[[1]], as.symbol("list")))) {
      node[[2]] <- as.call(list(as.symbol(ET_TYPE_VAR), v))
    }
    return(node)
  }
  for (k in seq_along(node)[-1]) node[[k]] <- wrap_inner_returns(node[[k]], skip_list)
  node
}

# A function with no protocol beyond its number of arguments, which keep the
# user's names: a latent's `index(t)`, an `allowed(data, i, value)` constraint.
# Its result is returned as computed (an index must stay an Int), and it does
# not take `model`, so a `model$x` read in it is refused rather than left to
# fail at run time.
transpile_plain_fn <- function(f, name, arity, helpers, what) {
  check_julia_name(name, "function name")
  if (inherits(f, "et_julia")) {
    stop(what, ": an et_julia() body is not accepted here; write the function ",
         "in R, or use a helper.", call. = FALSE)
  }
  if (!is.function(f) || length(formals(f)) != arity) {
    stop(what, " must be a function of ", arity, " argument(s).", call. = FALSE)
  }
  args <- names(formals(f))
  check_julia_name(args, "argument name")
  helper_names <- vapply(helpers, function(h) h$name, character(1))
  ctx <- new_ctx(helpers = helper_names)
  src <- transpile_braced_body(body(f), ctx)
  reads <- union(ctx$reads$model, helper_reads(ctx$reads$helpers, helpers))
  if (length(reads) || isTRUE(ctx$reads$uses_type)) {
    stop(what, " cannot read the parameters (it reads ",
         paste0("model$", reads, collapse = ", "), "); it is evaluated without them.",
         call. = FALSE)
  }
  sprintf("function %s(%s)\n%s\nend", name, paste(args, collapse = ", "), indent(src))
}

# A prior that depends on the parameters, `function(model) bernoulli_dist(model$p)`:
# the body is one call to a distribution constructor, emitted as the matching
# Distributions.jl constructor over transpiled arguments.
.et_dist_ctors <- c(normal_dist = "Normal", beta_dist = "Beta", gamma_dist = "Gamma",
                    uniform_dist = "Uniform", bernoulli_dist = "Bernoulli",
                    discrete_uniform_dist = "DiscreteUniform",
                    exponential_dist = "Exponential")

transpile_param_prior <- function(f, name, helpers, what) {
  if (!is.function(f) || !identical(names(formals(f)), "model")) {
    stop(what, ": a prior that depends on the parameters must be ",
         "function(model) <distribution>.", call. = FALSE)
  }
  b <- body(f)
  if (is.call(b) && identical(as.character(b[[1]]), "{")) {
    if (length(b) != 2L) {
      stop(what, ": the prior's body must be a single distribution call.",
           call. = FALSE)
    }
    b <- b[[2]]
  }
  head <- if (is.call(b) && is.symbol(b[[1]])) as.character(b[[1]]) else ""
  if (!head %in% names(.et_dist_ctors)) {
    stop(what, ": the prior's body must be a call to one of ",
         paste0(names(.et_dist_ctors), "()", collapse = ", "), ".", call. = FALSE)
  }
  args <- as.list(b)[-1L]
  if (!is.null(names(args)) && any(nzchar(names(args)))) {
    stop(what, ": pass the distribution's arguments positionally.", call. = FALSE)
  }
  helper_names <- vapply(helpers, function(h) h$name, character(1))
  ctx <- new_ctx(helpers = helper_names)
  src <- vapply(args, et_transpile_expr, character(1), ctx = ctx)
  # exponential_dist() takes a rate; Julia's Exponential takes the scale.
  if (head == "exponential_dist") src <- sprintf("1 / (%s)", if (length(src)) src else "1")
  list(src = sprintf("%s(model) = Distributions.%s(%s)", name, .et_dist_ctors[[head]],
                     paste(src, collapse = ", ")),
       dist = head)
}

helper_reads <- function(used, helpers) {
  if (!length(used) || !length(helpers)) return(character())
  out <- character()
  for (h in helpers) {
    if (!(h$name %in% helper_closure(used, helpers))) next
    if (inherits(h$f, "et_julia")) {
      out <- union(out, h$f$reads %||% character())
      next
    }
    ctx <- new_ctx(helpers = vapply(helpers, function(x) x$name, character(1)))
    transpile_braced_body(body(h$f), ctx)
    out <- union(out, ctx$reads$model)
  }
  out
}

# Every helper reachable from `used`, through helpers calling helpers. A helper
# two calls deep still reads parameters, and missing them would drop a
# parameter from `depends=` -- a Gibbs block holding only that parameter would
# then skip the likelihood and sample its prior.
helper_closure <- function(used, helpers) {
  names_all <- vapply(helpers, function(h) h$name, character(1))
  seen <- intersect(used, names_all)
  frontier <- seen
  while (length(frontier)) {
    nxt <- character()
    for (h in helpers[names_all %in% frontier]) {
      if (inherits(h$f, "et_julia")) next
      ctx <- new_ctx(helpers = names_all)
      transpile_braced_body(body(h$f), ctx)
      nxt <- union(nxt, ctx$reads$helpers)
    }
    frontier <- setdiff(nxt, seen)
    seen <- union(seen, frontier)
  }
  seen
}

# Emit every helper definition, in declaration order.
et_transpile_helpers <- function(helpers) {
  if (!length(helpers)) return(character())
  names_all <- vapply(helpers, function(h) h$name, character(1))
  vapply(helpers, function(h) {
    if (inherits(h$f, "et_julia")) {
      return(sprintf("function %s(%s)\n%s\nend", h$name,
                     paste(h$args, collapse = ", "), indent(h$f$src)))
    }
    ctx <- new_ctx(helpers = names_all)
    body_src <- transpile_braced_body(body(h$f), ctx)
    # A helper gets the `ETR_T` prologue too, but only when its body actually
    # uses the type variable (via et_one() / et_zero() / et_num()). Helpers keep
    # the user's own argument names and need not take `model` at all, so the
    # prologue can only be emitted when they do.
    if (isTRUE(ctx$reads$uses_type)) {
      if (!("model" %in% h$args)) {
        stop("the helper '", h$name, "' uses et_one() / et_zero() / et_num(), ",
             "which name the PARAMETER scalar type - but it does not take ",
             "`model`, so that type cannot be determined. Either add `model` to ",
             "its arguments, or use plain numeric literals and let the caller ",
             "convert the result.", call. = FALSE)
      }
      body_src <- paste0(eltype_prologue(ctx$reads$model), "\n", body_src)
    }
    sprintf("function %s(%s)\n%s\nend", h$name, paste(h$args, collapse = ", "),
            indent(body_src))
  }, character(1))
}

check_protocol_formals <- function(f, expected, role, name) {
  got <- names(formals(f))
  if (!identical(got, expected)) {
    stop("the ", role, " '", name, "' must be declared as\n",
         "    function(", paste(expected, collapse = ", "), ") ...\n",
         "  but its arguments are (", paste(got, collapse = ", "), ").\n",
         "  ET fixes this signature; the names are how the body refers to the ",
         "parameters, the data and the indices.", call. = FALSE)
  }
  invisible(TRUE)
}
