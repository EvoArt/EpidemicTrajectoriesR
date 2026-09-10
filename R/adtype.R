# Layer 3: choosing the automatic-differentiation backend.
#
# The gradient is ~83% of a badger Gibbs sweep, so the backend is one of the
# few genuinely consequential knobs an R user has. It is also the thing that
# decides what a parameter is MADE of during a gradient, which is why the
# transpiler never assumes Duals -- see ET_TYPE_VAR in transpile_fn.R.

# name -> (ADTypes constructor, Julia package that must be loaded first)
.et_ad_backends <- list(
  forwarddiff = list(ctor = "AutoForwardDiff", pkg = NULL),
  polyester   = list(ctor = "AutoPolyesterForwardDiff", pkg = "PolyesterForwardDiff"),
  reversediff = list(ctor = "AutoReverseDiff", pkg = "ReverseDiff"),
  mooncake    = list(ctor = "AutoMooncake", pkg = "Mooncake"),
  enzyme      = list(ctor = "AutoEnzyme", pkg = "Enzyme")
)

#' Choose the automatic-differentiation backend.
#'
#' @param name One of `"forwarddiff"` (the default), `"polyester"`,
#'   `"reversediff"`, `"mooncake"`, `"enzyme"`.
#' @param chunksize Forward-mode chunk size. `NULL` lets the backend choose.
#'   Cost is roughly linear in the number of CHUNKS, not in the parameter count,
#'   so this interacts strongly with how you block the parameters.
#' @param args Extra named arguments spliced into the `ADTypes` constructor.
#' @return An object of class `et_adtype`.
#' @details
#' **Forward mode** (`forwarddiff`, `polyester`) carries the derivative alongside
#' the value, so during a gradient every parameter is a `ForwardDiff.Dual`. Cost
#' scales with the number of chunk passes over the parameters, which makes it
#' strong for small blocks and weak for large ones. `polyester` is the threaded
#' variant and is what the badger model was tuned with.
#'
#' **Reverse mode** (`mooncake`, `enzyme`, `reversediff`) computes all partials
#' in one reverse pass, so it does not care how many parameters a block has --
#' the usual choice once a block is large. `mooncake` and `enzyme` leave the
#' primal as plain `Float64`; `reversediff` uses a tracked scalar type.
#'
#' Backends other than the two forward ones are not installed by default; add
#' them with `et_setup(ad_backends = "mooncake")`, which resolves them into this
#' package's pinned Julia project.
#'
#' No backend changes the answer -- only the time taken. If two backends disagree
#' on a log density or a gradient, that is a bug, not a tuning question.
#' @export
et_adtype <- function(name = c("forwarddiff", "polyester", "reversediff",
                               "mooncake", "enzyme"),
                      chunksize = NULL, args = list()) {
  name <- match.arg(name)
  structure(list(name = name, chunksize = chunksize, args = args),
            class = "et_adtype")
}

#' @export
print.et_adtype <- function(x, ...) {
  cat("<et_adtype> ", adtype_to_julia(x), "\n", sep = "")
  invisible(x)
}

# Render as an ADTypes constructor call.
adtype_to_julia <- function(x) {
  if (is.character(x)) x <- et_adtype(x)
  if (!inherits(x, "et_adtype")) {
    stop("expected an et_adtype() or a backend name.", call. = FALSE)
  }
  spec <- .et_ad_backends[[x$name]]
  kw <- character()
  if (!is.null(x$chunksize)) {
    kw <- c(kw, sprintf("chunksize=%s", julia_int(x$chunksize)))
  } else if (x$name == "polyester") {
    # AutoPolyesterForwardDiff has no default for these two.
    kw <- c(kw, "chunksize=nothing", "tag=nothing")
  }
  for (nm in names(x$args)) {
    v <- x$args[[nm]]
    kw <- c(kw, sprintf("%s=%s", nm,
                        if (is.character(v)) v else julia_num(v)))
  }
  sprintf("ADTypes.%s(%s)", spec$ctor,
          if (length(kw)) paste0("; ", paste(kw, collapse = ", ")) else "")
}

# The Julia package a backend needs loaded, if any.
adtype_package <- function(x) {
  if (is.character(x)) x <- et_adtype(x)
  .et_ad_backends[[x$name]]$pkg
}

#' Add and load extra automatic-differentiation backends.
#'
#' Resolves the named backends into this package's pinned Julia project and
#' loads them, so [et_adtype()] can select them. Forward mode needs nothing
#' extra; the reverse-mode backends are large, so they are opt-in.
#'
#' @param backends Character vector of backend names, e.g. `"mooncake"`.
#' @return Invisibly `TRUE`.
#' @export
et_add_ad_backend <- function(backends) {
  et_require_session()
  pkgs <- unique(stats::na.omit(vapply(backends, function(b) {
    spec <- .et_ad_backends[[b]]
    if (is.null(spec)) {
      stop("unknown AD backend '", b, "'. One of: ",
           paste(names(.et_ad_backends), collapse = ", "), ".", call. = FALSE)
    }
    spec$pkg %||% NA_character_
  }, character(1))))
  for (p in pkgs) {
    message("resolving Julia package ", p, " (this can take a few minutes) ...")
    JuliaCall::julia_command(sprintf('using Pkg; Pkg.add("%s"; io=devnull)', p))
    JuliaCall::julia_command(sprintf("using %s", p))
  }
  invisible(TRUE)
}
