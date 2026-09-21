# Priors: a naming layer over Distributions.jl. Only the families epidemic
# models actually use are named; custom_dist() covers the rest.

new_dist <- function(name, args, vector_arg = FALSE) {
  structure(list(name = name, args = args, vector_arg = vector_arg),
            class = "custom_dist")
}

#' Prior distributions.
#'
#' These map onto `Distributions.jl` constructors of the same name.
#'
#' @param mean,sd Normal parameters.
#' @param shape,scale Gamma parameters (Julia's `Gamma(shape, scale)`).
#' @param a,b Beta / Uniform parameters.
#' @param rate Exponential rate; emitted as `Exponential(1 / rate)` because
#'   Julia's `Exponential` is parameterised by scale, which is the single most
#'   likely thing for an R user to get backwards.
#' @param alpha Dirichlet concentration vector.
#' @param d A distribution to truncate.
#' @param lower,upper Truncation bounds; `NA` for unbounded on that side.
#' @param name,args Escape hatch: a Julia distribution name and its arguments.
#' @return An object of class `custom_dist`.
#' @name prior_distributions
NULL

#' @rdname prior_distributions
#' @export
normal_dist <- function(mean = 0, sd = 1) new_dist("Normal", list(mean, sd))

#' @rdname prior_distributions
#' @export
gamma_dist <- function(shape, scale) new_dist("Gamma", list(shape, scale))

#' @rdname prior_distributions
#' @export
beta_dist <- function(a = 1, b = 1) new_dist("Beta", list(a, b))

#' @rdname prior_distributions
#' @export
exponential_dist <- function(rate = 1) {
  if (!is.numeric(rate) || length(rate) != 1L || is.na(rate) || rate <= 0) {
    stop("exponential_dist(): `rate` must be a single positive number.", call. = FALSE)
  }
  new_dist("Exponential", list(1 / rate))
}

#' @rdname prior_distributions
#' @export
uniform_dist <- function(a = 0, b = 1) new_dist("Uniform", list(a, b))

#' @rdname prior_distributions
#' @export
dirichlet_dist <- function(alpha) new_dist("Dirichlet", list(as.numeric(alpha)),
                                         vector_arg = TRUE)

#' @rdname prior_distributions
#' @export
truncated_dist <- function(d, lower = NA, upper = NA) {
  if (!inherits(d, "custom_dist")) {
    stop("truncated_dist(): `d` must be a prior from this package.", call. = FALSE)
  }
  structure(list(name = "truncated", args = list(d, lower, upper),
                 vector_arg = FALSE, truncate = TRUE), class = "custom_dist")
}

#' @rdname prior_distributions
#' @export
custom_dist <- function(name, args = list()) {
  check_julia_name(name, "distribution name")
  new_dist(name, as.list(args))
}

#' @export
print.custom_dist <- function(x, ...) {
  cat("<custom_dist> ", dist_to_julia(x), "\n", sep = "")
  invisible(x)
}

# Render a prior as Julia source.
dist_to_julia <- function(d) {
  if (isTRUE(d$truncate)) {
    lo <- if (is.na(d$args[[2]])) "-Inf" else julia_float(d$args[[2]])
    hi <- if (is.na(d$args[[3]])) "Inf" else julia_float(d$args[[3]])
    return(sprintf("truncated(%s, %s, %s)", dist_to_julia(d$args[[1]]), lo, hi))
  }
  if (isTRUE(d$vector_arg)) {
    return(sprintf("%s(%s)", d$name, julia_vector(d$args[[1]], "Float64")))
  }
  args <- vapply(d$args, function(a) {
    if (inherits(a, "custom_dist")) return(dist_to_julia(a))
    if (length(a) > 1L) return(julia_vector(a, "Float64"))
    julia_float(a)
  }, character(1))
  sprintf("%s(%s)", d$name, paste(args, collapse = ", "))
}

# ---- deprecated spellings ---------------------------------------------------
#
# The distributions were `et_beta()` and friends, and a parameter was declared
# with `et_par()`. The prefix was there to avoid collisions, but `beta` and
# `gamma` really are base R functions, so the suffix form keeps the protection
# while reading as a noun. `et_par()` was the worse offender: it said
# "parameter" twice and never said "prior".
#
# Kept working, with a warning, so existing scripts do not break.

.deprecated <- function(old, new) {
  force(old); force(new)
  function(...) {
    warning(old, "() is deprecated; use ", new, "()", call. = FALSE)
    do.call(new, list(...))
  }
}

#' @rdname prior_distributions
#' @export
et_normal <- .deprecated("et_normal", "normal_dist")
#' @rdname prior_distributions
#' @export
et_gamma <- .deprecated("et_gamma", "gamma_dist")
#' @rdname prior_distributions
#' @export
et_beta <- .deprecated("et_beta", "beta_dist")
#' @rdname prior_distributions
#' @export
et_exponential <- .deprecated("et_exponential", "exponential_dist")
#' @rdname prior_distributions
#' @export
et_uniform <- .deprecated("et_uniform", "uniform_dist")
#' @rdname prior_distributions
#' @export
et_dirichlet <- .deprecated("et_dirichlet", "dirichlet_dist")
#' @rdname prior_distributions
#' @export
et_truncated <- .deprecated("et_truncated", "truncated_dist")
#' @rdname prior_distributions
#' @export
et_dist <- .deprecated("et_dist", "custom_dist")
