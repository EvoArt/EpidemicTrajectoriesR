# Layer 1 (substrate): typed literals inside a transpiled body.

#' Typed literals for a transpiled body.
#'
#' An accumulator started at a bare `1` is an integer, and becomes a
#' `Union{Int, Dual}` in Julia the moment a parameter multiplies into it --
#' runtime dispatch inside the autodiff loop, which is the single most expensive
#' thing to get wrong in this stack. `et_one()` and `et_zero()` start it at the
#' automatic-differentiation scalar type instead, and `et_num()` converts a value
#' to it.
#'
#' Use them for a value that will be combined with parameters. Do NOT use them
#' for an index or a loop counter -- those must stay integers.
#'
#' In R these are the identity, so a function using them still runs unchanged
#' outside the transpiler.
#'
#' @param x A value to convert.
#' @return `1`, `0`, or `x`.
#' @name et_typed_literals
#' @examples
#' # w accumulates a product of test-parameter terms, so it starts typed:
#' f <- function(model, data, X, i, t, s) {
#'   w <- et_one()
#'   w * model$theta
#' }
NULL

#' @rdname et_typed_literals
#' @export
et_one <- function() 1

#' @rdname et_typed_literals
#' @export
et_zero <- function() 0

#' @rdname et_typed_literals
#' @export
et_num <- function(x) x
