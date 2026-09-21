# Priors: a naming layer over Distributions.jl. Only the families epidemic
# models actually use are named; et_dist() covers the rest.

new_dist <- function(name, args, vector_arg = FALSE) {
  structure(list(name = name, args = args, vector_arg = vector_arg),
            class = "et_dist")
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
#' @return An object of class `et_dist`.
#' @name et_priors
NULL

#' @rdname et_priors
#' @export
et_normal <- function(mean = 0, sd = 1) new_dist("Normal", list(mean, sd))

#' @rdname et_priors
#' @export
et_gamma <- function(shape, scale) new_dist("Gamma", list(shape, scale))

#' @rdname et_priors
#' @export
et_beta <- function(a = 1, b = 1) new_dist("Beta", list(a, b))

#' @rdname et_priors
#' @export
et_exponential <- function(rate = 1) {
  if (!is.numeric(rate) || length(rate) != 1L || is.na(rate) || rate <= 0) {
    stop("et_exponential(): `rate` must be a single positive number.", call. = FALSE)
  }
  new_dist("Exponential", list(1 / rate))
}

#' @rdname et_priors
#' @export
et_uniform <- function(a = 0, b = 1) new_dist("Uniform", list(a, b))

#' @rdname et_priors
#' @export
et_dirichlet <- function(alpha) new_dist("Dirichlet", list(as.numeric(alpha)),
                                         vector_arg = TRUE)

#' @rdname et_priors
#' @export
et_truncated <- function(d, lower = NA, upper = NA) {
  if (!inherits(d, "et_dist")) {
    stop("et_truncated(): `d` must be a prior from this package.", call. = FALSE)
  }
  structure(list(name = "truncated", args = list(d, lower, upper),
                 vector_arg = FALSE, truncate = TRUE), class = "et_dist")
}

#' @rdname et_priors
#' @export
et_dist <- function(name, args = list()) {
  check_julia_name(name, "distribution name")
  new_dist(name, as.list(args))
}

#' @export
print.et_dist <- function(x, ...) {
  cat("<et_dist> ", dist_to_julia(x), "\n", sep = "")
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
    if (inherits(a, "et_dist")) return(dist_to_julia(a))
    if (length(a) > 1L) return(julia_vector(a, "Float64"))
    julia_float(a)
  }, character(1))
  sprintf("%s(%s)", d$name, paste(args, collapse = ", "))
}
