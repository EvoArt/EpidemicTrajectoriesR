# Layer 1/2 boundary: function-level transpilation into ET's call protocols.
#
# `transpile.R` turns an expression into Julia text. This file wraps a whole R
# function in the exact Julia signature ET will call it by, checks the R formals
# against that protocol, and emits the element-type prologue for the roles that
# return a per-state vector.
#
# Unlike BadgeR's twelve fixed engine roles, these are ET's actual call
# conventions -- the ones documented on `epidemic_data`, `@transitions` and
# `epidemic_obs_loglik`. A model's STRUCTURE is declared separately (spec.R), so
# this table does not grow when a user writes a different model.

.et_protocols <- list(
  rate           = list(args = c("model", "data", "i", "t"),           vector = FALSE),
  survival       = list(args = c("model", "data", "i", "t"),           vector = FALSE),
  starting_state = list(args = c("model", "data", "X", "i", "t"),      vector = TRUE),
  obs_process    = list(args = c("model", "data", "X", "i", "t"),      vector = TRUE),
  obs_weight     = list(args = c("model", "data", "X", "i", "t", "s"), vector = FALSE)
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
#' forgo the speed-up described in DESIGN.md section 5. Pass `reads` to declare
#' the parameters the body touches and keep the annotation.
#'
#' @param src Julia source for the function BODY (not the signature -- the
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
et_transpile_role <- function(f, role, name, helpers = list()) {
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
    return(list(src = src, reads = f$reads %||% character(),
                opaque = is.null(f$reads), helpers = character()))
  }
  if (!is.function(f)) {
    stop("the ", role, " '", name, "' must be an R function or an et_julia() ",
         "object; got ", class(f)[1], ".", call. = FALSE)
  }
  check_protocol_formals(f, spec$args, role, name)

  ctx <- new_ctx(helpers = helper_names, vector_result = spec$vector)
  # A SCALAR role must return one concrete type. An R body naturally mixes
  # `return(1)` (an Int) with `return(1 - eta)` (whatever the parameters are),
  # which makes the Julia return type a Union -- runtime dispatch in the hottest
  # function in the package. Converting every returned value to `ETR_T` fixes it
  # without the user having to think about it, and `ETR_T` is the promotion of
  # exactly the parameters this body reads.
  fbody <- if (spec$vector) body(f) else wrap_returns_in_T(body(f))
  body_src <- transpile_braced_body(fbody, ctx)
  reads <- ctx$reads$model
  used_helpers <- ctx$reads$helpers

  # Reads reached through helpers count too -- that transitivity is exactly what
  # makes an automatically derived `depends=` trustworthy.
  reads <- union(reads, helper_reads(used_helpers, helpers))
  opaque <- any(vapply(helpers[helper_names %in% used_helpers],
                       function(h) inherits(h$f, "et_julia") && is.null(h$f$reads),
                       logical(1)))

  prologue <- eltype_prologue(reads)
  full <- if (length(prologue)) paste0(prologue, "\n", body_src) else body_src
  list(src = sprintf("function %s(%s)\n%s\nend", name,
                     paste(spec$args, collapse = ", "), indent(full)),
       reads = reads, opaque = opaque, helpers = used_helpers)
}

# The name of the parameter-scalar-type variable the generator binds. NOT `T`:
# R spells `TRUE` as `T`, so a body containing a bare `T` would be ambiguous, and
# a user local named `T` would shadow the binding.
ET_TYPE_VAR <- "ETR_T"

# `ETR_T` is the promotion of the element types of every parameter the body
# reads. It is NOT "the Dual type": what the parameters are made of depends on
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
# Everything the generator emits against it -- `convert(ETR_T, x)`, `one(ETR_T)`,
# `zeros(ETR_T, n)` -- is defined across all three families. A type CONSTRUCTOR
# call would not be.
eltype_prologue <- function(reads) {
  if (!length(reads)) return(paste0(ET_TYPE_VAR, " = Float64"))
  sprintf("%s = promote_type(%s)", ET_TYPE_VAR,
          paste(sprintf("eltype(model.%s)", reads), collapse = ", "))
}

# AST transform: every value the function can return is wrapped in `T(...)`.
#   `return(x)`        -> `return(T(x))`
#   `{ ...; last }`    -> `{ ...; T(last) }`  (unless the last is a return)
#   `if (c) a else b`  -> branches wrapped individually, so the emitted Julia
#                         keeps its readable multi-line `if` form
# A `for` loop or an assignment as the terminal statement returns nothing
# meaningful, so it is left alone.
wrap_returns_in_T <- function(node) {
  Tcall <- function(x) as.call(list(as.symbol(ET_TYPE_VAR), x))
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
      for (k in 2:(n - 1)) node[[k]] <- wrap_inner_returns(node[[k]])
      node[[n]] <- wrap_returns_in_T(node[[n]])
    }
    return(node)
  }
  if (head == "if") {
    node[[3]] <- wrap_returns_in_T(node[[3]])
    if (length(node) == 4L) node[[4]] <- wrap_returns_in_T(node[[4]])
    return(node)
  }
  if (head %in% c("for", "while", "<-", "=")) return(node)
  Tcall(node)
}

# Wrap `return(x)` anywhere inside a statement, without touching its value.
wrap_inner_returns <- function(node) {
  if (!is.call(node)) return(node)
  head <- as.character(node[[1]])
  if (head == "return" && length(node) == 2L) {
    node[[2]] <- as.call(list(as.symbol(ET_TYPE_VAR), node[[2]]))
    return(node)
  }
  for (k in seq_along(node)[-1]) node[[k]] <- wrap_inner_returns(node[[k]])
  node
}

helper_reads <- function(used, helpers) {
  if (!length(used) || !length(helpers)) return(character())
  out <- character()
  for (h in helpers) {
    if (!(h$name %in% used)) next
    if (inherits(h$f, "et_julia")) {
      out <- union(out, h$f$reads %||% character())
      next
    }
    ctx <- new_ctx(helpers = vapply(helpers, function(x) x$name, character(1)))
    transpile_braced_body(body(h$f), ctx)
    out <- union(out, ctx$reads$model)
    # One level of nesting is followed here; deeper chains are resolved by the
    # fixed point in et_helper_closure().
  }
  out
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
