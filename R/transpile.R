# The R -> Julia expression transpiler. Walks the parsed R AST and emits
# Julia source.
#
# Deliberately narrow: only the constructs needed for the leaf functions ET
# calls (rates, survival, starting states, observation weights) and their
# helpers. Anything else errors, naming the offending node. Emitting
# something plausible would give a model that runs and is wrong.
#
# `model$x` maps to plain `model.x`, since ET's `model` is a NamedTuple of
# values. Each such read is recorded, which is what lets codegen.R work out
# `depends=` instead of asking for it. Indexing passes through unchanged --
# both languages are 1-based and column-major.

# ---- transpiler context -----------------------------------------------------

#' Create a transpiler context.
#'
#' @param helpers Character vector of user helper names callable from a body.
#' @param vector_result `TRUE` inside a role that returns a per-state vector, in
#'   which `rep(1, n)` / `numeric(n)` become `ones(ETR_T, n)` / `zeros(ETR_T, n)`
#'   so the allocation follows the parameter scalar type (see DESIGN.md 3.3).
#' @param state_syms `TRUE` where a string naming a state should be resolved
#'   (aggregate update bodies): against `state` it becomes a Julia Symbol,
#'   against `X[...]` the state's integer code.
#' @param states The state space, in encoding order. Needed to resolve a name
#'   compared against the trajectory to its code.
#' @param reads An environment accumulating the `model$x` names read.
#' @keywords internal
new_ctx <- function(helpers = character(), vector_result = FALSE,
                    state_syms = FALSE, states = character(), reads = NULL) {
  if (is.null(reads)) {
    reads <- new.env(parent = emptyenv())
    reads$model <- character()
    reads$helpers <- character()
    reads$uses_type <- FALSE
  }
  list(helpers = helpers, vector_result = vector_result,
       state_syms = state_syms, states = states, reads = reads)
}

# Maths / builtin calls allowed in a body, mapped R-name -> Julia-name.
.et_builtin_map <- c(
  exp = "exp", log = "log", log1p = "log1p", expm1 = "expm1",
  sqrt = "sqrt", abs = "abs", sum = "sum", prod = "prod",
  min = "min", max = "max", floor = "floor", ceiling = "ceil",
  round = "round", lgamma = "lgamma", gamma = "gamma",
  length = "length", clamp = "clamp",
  # `any`/`all` are only meaningful on the scalar comparisons this subset can
  # build, but they behave identically in both languages there.
  any = "any", all = "all"
)

# Binary operators mapped R -> Julia. R's `&&`/`||` and `&`/`|` mean the same
# things in Julia; `%%` and `%/%` do not have identical negative-operand
# semantics and so are not included.
.et_binop_map <- c(
  "+" = "+", "-" = "-", "*" = "*", "/" = "/", "^" = "^",
  "==" = "==", "!=" = "!=", "<" = "<", ">" = ">",
  "<=" = "<=", ">=" = ">=",
  "&&" = "&&", "||" = "||", "&" = "&", "|" = "|"
)

#' Transpile one R expression (language object) to Julia source.
#'
#' @param e An R language object, symbol or literal, e.g. from `body(f)`.
#' @param ctx A context from [new_ctx()].
#' @return A length-1 character string of Julia source.
#' @export
et_transpile_expr <- function(e, ctx = new_ctx()) {
  if (is.symbol(e)) return(transpile_symbol(e, ctx))
  if (is.numeric(e) || is.logical(e)) return(transpile_literal(e))
  if (is.character(e)) return(julia_string(e))
  if (is.null(e)) return("nothing")
  if (!is.call(e)) {
    stop("transpile: unsupported node of class '", class(e)[1], "'.", call. = FALSE)
  }
  head <- e[[1]]
  if (!is.symbol(head)) {
    stop("transpile: unsupported call head -- only named functions and operators ",
         "are supported (no `f(x)(y)`, no `obj$method()`).", call. = FALSE)
  }
  op <- as.character(head)

  if (op == "{") return(transpile_block(e, ctx))
  if (op == "(") return(paste0("(", et_transpile_expr(e[[2]], ctx), ")"))

  if (op %in% c("<-", "=")) {
    return(paste0(et_transpile_expr(e[[2]], ctx), " = ",
                  et_transpile_expr(e[[3]], ctx)))
  }
  if (op == "<<-") {
    stop("transpile: `<<-` is not supported. A transpiled function may only ",
         "assign to its own locals.", call. = FALSE)
  }
  if (op == "return") {
    if (length(e) == 1L) return("return")
    return(paste0("return ", et_transpile_expr(e[[2]], ctx)))
  }
  if (op == "if") return(transpile_if(e, ctx))
  if (op == "for") return(transpile_for(e, ctx))
  if (op == "next") return("continue")
  if (op == "break") return("break")
  if (op %in% c("while", "repeat")) {
    stop("transpile: `", op, "` loops are not supported (a leaf function that ",
         "needs one is a sign the work belongs in Julia -- use et_julia()).",
         call. = FALSE)
  }

  # unary forms
  if (op %in% c("-", "+") && length(e) == 2L) {
    return(paste0(op, et_transpile_expr(e[[2]], ctx)))
  }
  if (op == "!" && length(e) == 2L) {
    return(paste0("!", et_transpile_expr(e[[2]], ctx)))
  }

  if (op %in% names(.et_binop_map) && length(e) == 3L) {
    return(transpile_binop(e, ctx, op))
  }

  if (op == "ifelse" && length(e) == 4L) {
    return(paste0("(", et_transpile_expr(e[[2]], ctx), " ? ",
                  et_transpile_expr(e[[3]], ctx), " : ",
                  et_transpile_expr(e[[4]], ctx), ")"))
  }

  if (op == "[") return(transpile_index(e, ctx))
  if (op == "[[") {
    stop("transpile: `[[` is not supported; use `[` for arrays and `$` for ",
         "fields of `model`/`data`.", call. = FALSE)
  }
  if (op == "$") return(transpile_dollar(e, ctx))

  if (op == ":" && length(e) == 3L) {
    return(paste0(et_transpile_expr(e[[2]], ctx), ":",
                  et_transpile_expr(e[[3]], ctx)))
  }

  transpile_call(e, ctx, op)
}

# ---- leaves -----------------------------------------------------------------

transpile_symbol <- function(e, ctx) {
  nm <- as.character(e)
  if (nm %in% c("TRUE", "T")) return("true")
  if (nm %in% c("FALSE", "F")) return("false")
  if (nm == "Inf") return("Inf")
  nm
}

transpile_literal <- function(e) {
  if (length(e) != 1L) {
    stop("transpile: vector literals of length ", length(e),
         " are not supported; build the value elementwise.", call. = FALSE)
  }
  if (is.logical(e)) return(if (isTRUE(e)) "true" else "false")
  if (is.integer(e)) return(as.character(e))
  v <- e
  # A whole-number double stays integral so it can still index or bound a loop.
  # Julia promotes to Float64 on contact with a float, so arithmetic is
  # unaffected; a non-integral value keeps full round-trip precision.
  if (is.finite(v) && v == round(v) && abs(v) < 1e15) {
    return(format(v, scientific = FALSE, trim = TRUE, nsmall = 0))
  }
  julia_num(v)
}

# ---- compound forms ---------------------------------------------------------

transpile_block <- function(e, ctx) {
  stmts <- as.list(e)[-1L]
  if (length(stmts) == 0L) return("")
  paste(vapply(stmts, et_transpile_expr, character(1), ctx = ctx), collapse = "\n")
}

transpile_braced_body <- function(node, ctx) {
  if (is.call(node) && identical(as.character(node[[1]]), "{")) {
    transpile_block(node, ctx)
  } else {
    et_transpile_expr(node, ctx)
  }
}

transpile_if <- function(e, ctx) {
  cond <- et_transpile_expr(e[[2]], ctx)
  then_src <- transpile_braced_body(e[[3]], ctx)
  if (length(e) == 3L) {
    return(paste0("if ", cond, "\n", indent(then_src), "\nend"))
  }
  els <- e[[4]]
  if (is.call(els) && identical(as.character(els[[1]]), "if")) {
    else_src <- transpile_if(els, ctx)                 # else-if chaining
    else_src <- sub("^if ", "elseif ", else_src)
    else_src <- sub("\nend$", "", else_src)
    return(paste0("if ", cond, "\n", indent(then_src), "\n", else_src, "\nend"))
  }
  paste0("if ", cond, "\n", indent(then_src), "\nelse\n",
         indent(transpile_braced_body(els, ctx)), "\nend")
}

transpile_for <- function(e, ctx) {
  var <- as.character(e[[2]])
  seqx <- et_transpile_expr(e[[3]], ctx)
  body <- transpile_braced_body(e[[4]], ctx)
  paste0("for ", var, " in ", seqx, "\n", indent(body), "\nend")
}

# `state == "I"` must reach Julia as `state == :I`: inside an @aggregate body ET
# compares the applied state against a Symbol from the declared state space.
# Comparing against a state name. Two different things wear the same syntax, and
# getting them the wrong way round is silent:
#
#   * `state` is the Symbol ET applies or reverses, so `state == "I"` is
#     `state == :I`;
#   * the trajectory holds integer CODES, so `X[t+1, i] == "E"` is
#     `X[t+1, i] == 2`, its position in the state space.
#
# Emitting `:E` for the second would compare an Int to a Symbol: never true, no
# error, and a tracked array that quietly counts nothing.
transpile_binop <- function(e, ctx, op) {
  lhs <- e[[2]]; rhs <- e[[3]]
  if (isTRUE(ctx$state_syms) && op %in% c("==", "!=")) {
    if (identical(lhs, quote(state)) && is.character(rhs)) {
      return(paste0("state ", op, " ", julia_symbol(rhs)))
    }
    if (identical(rhs, quote(state)) && is.character(lhs)) {
      return(paste0(julia_symbol(lhs), " ", op, " state"))
    }
    if (is.character(rhs) && is_trajectory_ref(lhs)) {
      return(paste0(et_transpile_expr(lhs, ctx), " ", op, " ",
                    state_code_or_stop(rhs, ctx)))
    }
    if (is.character(lhs) && is_trajectory_ref(rhs)) {
      return(paste0(state_code_or_stop(lhs, ctx), " ", op, " ",
                    et_transpile_expr(rhs, ctx)))
    }
  }
  paste0(et_transpile_expr(lhs, ctx), " ", .et_binop_map[[op]], " ",
         et_transpile_expr(rhs, ctx))
}

# `X[...]` -- an index into the trajectory, which holds state codes.
is_trajectory_ref <- function(e) {
  is.call(e) && identical(as.character(e[[1]]), "[") &&
    is.symbol(e[[2]]) && identical(as.character(e[[2]]), "X")
}

state_code_or_stop <- function(name, ctx) {
  k <- match(name, ctx$states %||% character())
  if (is.na(k)) {
    stop("transpile: '", name, "' is not a state. States are: ",
         paste(ctx$states, collapse = ", "),
         ". Comparing the trajectory against a name the state space does not ",
         "contain would silently never match.", call. = FALSE)
  }
  as.character(k)
}

transpile_index <- function(e, ctx) {
  obj <- e[[2]]
  idx <- as.list(e)[-(1:2)]
  if (any(vapply(idx, function(a) identical(a, quote(expr = )), logical(1)))) {
    stop("transpile: empty index slots (`x[, j]`) are not supported; write the ",
         "loop out, or use et_julia() for a whole-slice operation.", call. = FALSE)
  }
  # `dim(x)[k]` is the idiomatic R spelling of Julia's `size(x, k)`.
  if (is.call(obj) && identical(as.character(obj[[1]]), "dim") &&
      length(obj) == 2L && length(idx) == 1L) {
    return(sprintf("size(%s, %s)", et_transpile_expr(obj[[2]], ctx),
                   et_transpile_expr(idx[[1]], ctx)))
  }
  idx_src <- vapply(idx, et_transpile_expr, character(1), ctx = ctx)
  paste0(et_transpile_expr(obj, ctx), "[", paste(idx_src, collapse = ", "), "]")
}

# `obj$name` -> `obj.name`. The only special case is `model$x`, which is
# recorded so the generator can derive `depends=` (DESIGN.md section 5).
transpile_dollar <- function(e, ctx) {
  obj <- e[[2]]
  field <- as.character(e[[3]])
  if (identical(obj, quote(model))) {
    ctx$reads$model <- union(ctx$reads$model, field)
    return(paste0("model.", field))
  }
  paste0(et_transpile_expr(obj, ctx), ".", field)
}

transpile_call <- function(e, ctx, op) {
  args <- as.list(e)[-1L]
  if (!is.null(names(args)) && any(nzchar(names(args)))) {
    stop("transpile: named arguments are not supported in a transpiled call ",
         "(`", op, "`); pass them positionally.", call. = FALSE)
  }
  arg_src <- vapply(args, et_transpile_expr, character(1), ctx = ctx)

  # Vector allocation. In a vector-returning role these must follow the PARAMETER
  # scalar type, or the function throws the moment its parameters move into an
  # HMC block under a backend whose parameters are not plain Float64
  # (the observation-eltype gotcha). `ETR_T` is bound by the prologue the
  # function-level transpiler emits.
  if (op %in% c("rep", "numeric", "rep_len")) return(transpile_alloc(op, arg_src, ctx))

  if (op %in% names(.et_builtin_map)) {
    return(paste0(.et_builtin_map[[op]], "(", paste(arg_src, collapse = ", "), ")"))
  }
  if (op == "nrow" && length(arg_src) == 1L) return(paste0("size(", arg_src, ", 1)"))
  if (op == "ncol" && length(arg_src) == 1L) return(paste0("size(", arg_src, ", 2)"))
  if (op == "seq_len" && length(arg_src) == 1L) return(paste0("1:", arg_src))
  if (op == "seq_along" && length(arg_src) == 1L) {
    return(paste0("eachindex(", arg_src, ")"))
  }
  # Typed literals. An accumulator started at a bare `1` is an Int, and becomes a
  # `Union{Int, <parameter type>}` once a parameter multiplies into it, runtime
  # dispatch in the hottest loop in the package. `et_one()` / `et_zero()` /
  # `et_num(x)` start it at the parameter scalar type instead. They are ordinary
  # R functions too, so a body using them still runs unchanged in R.
  if (op %in% c("et_one", "et_zero") && length(arg_src) == 0L) {
    ctx$reads$uses_type <- TRUE
    return(sprintf("%s(%s)", if (op == "et_one") "one" else "zero", ET_TYPE_VAR))
  }
  if (op == "et_num" && length(arg_src) == 1L) {
    ctx$reads$uses_type <- TRUE
    return(sprintf("convert(%s, %s)", ET_TYPE_VAR, arg_src))
  }
  # `ETR_T` is the parameter scalar type bound by the role's prologue; the
  # function-level transpiler inserts these calls, a user body never writes them.
  if (op == ET_TYPE_VAR && length(arg_src) == 1L) {
    ctx$reads$uses_type <- TRUE
    # `convert`, not the type's constructor: it is the idiom every numeric type
    # in the AD ecosystem implements, including the tracked types reverse-mode
    # backends use. A constructor call is not.
    return(sprintf("convert(%s, %s)", ET_TYPE_VAR, arg_src))
  }
  if (op %in% ctx$helpers) {
    ctx$reads$helpers <- union(ctx$reads$helpers, op)
    return(paste0(op, "(", paste(arg_src, collapse = ", "), ")"))
  }
  stop("transpile: call to unsupported function '", op, "'.\n",
       "  Allowed: ", paste(sort(names(.et_builtin_map)), collapse = ", "),
       ", nrow, ncol, seq_len, seq_along, rep, numeric, ifelse,\n",
       "  plus operators and any function declared with et_helper().\n",
       "  If '", op, "' is your own function, declare it as a helper; if it has ",
       "no Julia equivalent, use et_julia().", call. = FALSE)
}

# `rep(x, n)` / `numeric(n)` -> a typed Julia allocation.
transpile_alloc <- function(op, arg_src, ctx) {
  ty <- if (isTRUE(ctx$vector_result)) ET_TYPE_VAR else "Float64"
  if (isTRUE(ctx$vector_result)) ctx$reads$uses_type <- TRUE
  if (op == "numeric") {
    if (length(arg_src) != 1L) {
      stop("transpile: numeric() takes exactly one length argument.", call. = FALSE)
    }
    return(sprintf("zeros(%s, %s)", ty, arg_src[1]))
  }
  if (length(arg_src) != 2L) {
    stop("transpile: ", op, "() is supported only as ", op,
         "(value, times).", call. = FALSE)
  }
  val <- arg_src[1]
  if (val %in% c("1", "1.0")) return(sprintf("ones(%s, %s)", ty, arg_src[2]))
  if (val %in% c("0", "0.0")) return(sprintf("zeros(%s, %s)", ty, arg_src[2]))
  sprintf("fill(%s(%s), %s)", ty, val, arg_src[2])
}
