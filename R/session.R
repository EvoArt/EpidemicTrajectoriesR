# One Julia session per R process. Everything here is idempotent: et_setup()
# can be called repeatedly, and a model can be regenerated and reloaded
# without module-name collisions.

.et_state <- new.env(parent = emptyenv())
.et_state$ready <- FALSE
.et_state$loaded <- character()

#' Start the Julia session and load the modelling packages.
#'
#' @param julia_home Directory containing the `julia` binary. `NULL` lets
#'   JuliaCall find it (respecting `JULIA_HOME`).
#' @param project Path to the Julia project to activate. Defaults to the pinned
#'   environment shipped with this package.
#' @param dev Named character vector of local checkouts to `Pkg.develop`, e.g.
#'   `c(EpidemicTrajectories = "~/.julia/dev/EpidemicTrajectories")`. Use this
#'   when benchmarking or debugging local source changes; without it the pinned
#'   GitHub versions are used, and a local edit is invisible.
#' @param ad_backends Extra automatic-differentiation backends to resolve and
#'   load, e.g. `"mooncake"`. Forward mode works out of the box; the
#'   reverse-mode backends are large, so they are opt-in. See [et_adtype()].
#' @param verbose Print Julia's setup chatter.
#' @param force Re-run even if the session is already up.
#' @return Invisibly `TRUE`.
#' @export
et_setup <- function(julia_home = NULL, project = NULL, dev = NULL,
                     ad_backends = NULL, verbose = FALSE, force = FALSE) {
  if (.et_state$ready && !force && is.null(dev) && is.null(ad_backends)) {
    return(invisible(TRUE))
  }
  if (!requireNamespace("JuliaCall", quietly = TRUE)) {
    stop("EpidemicTrajectoriesR needs the JuliaCall package: ",
         'install.packages("JuliaCall").', call. = FALSE)
  }
  if (is.null(project)) project <- et_julia_project()
  JuliaCall::julia_setup(JULIA_HOME = julia_home, verbose = verbose)
  JuliaCall::julia_command(sprintf('using Pkg; Pkg.activate("%s"; io=devnull)',
                                   julia_path(project)))
  if (!is.null(dev)) {
    for (nm in names(dev)) {
      JuliaCall::julia_command(sprintf('Pkg.develop(path="%s"; io=devnull)',
                                       julia_path(normalizePath(dev[[nm]], mustWork = TRUE))))
    }
  }
  JuliaCall::julia_command("Pkg.instantiate(; io=devnull)")
  for (pkg in c("EpidemicTrajectories", "PracticalBayes", "PracticalEpiBayes",
                "Distributions", "StableRNGs", "AbstractMCMC", "AdvancedHMC",
                "ADTypes")) {
    JuliaCall::julia_command(sprintf("using %s", pkg))
  }
  # PolyesterForwardDiff is a speed option, not a requirement, so a session
  # starts without it. Its DifferentiationInterface extension fails to
  # precompile on some Julia 1.11 resolutions, and loading it unconditionally
  # made that failure fatal to every fit rather than to the one backend nobody
  # had asked for. et_adtype("polyester") reports the absence if it is wanted.
  .et_state$polyester <- isTRUE(tryCatch({
    JuliaCall::julia_command("using PolyesterForwardDiff"); TRUE
  }, error = function(e) FALSE))
  .et_state$ready <- TRUE
  if (!is.null(ad_backends)) et_add_ad_backend(ad_backends)
  invisible(TRUE)
}

#' Is the Julia session up?
#' @return `TRUE` or `FALSE`.
#' @export
et_is_ready <- function() isTRUE(.et_state$ready)

#' The pinned Julia project shipped with this package.
#' @return A path.
#' @export
et_julia_project <- function() {
  p <- system.file("julia", package = "EpidemicTrajectoriesR")
  if (!nzchar(p)) {
    # Running from a source checkout rather than an installed package.
    p <- file.path(getwd(), "inst", "julia")
  }
  p
}

#' Evaluate Julia source in the session.
#'
#' The escape hatch for anything this package does not wrap: ET's residuals
#' layer, `check_iffbs_exact`, a bespoke kernel. A loaded model's module is in
#' scope by name.
#'
#' @param src Julia source.
#' @param need_return Return the value to R (`TRUE`) or run for effect.
#' @return The Julia value, converted by JuliaCall, or invisibly `NULL`.
#' @export
et_eval <- function(src, need_return = TRUE) {
  et_require_session()
  if (need_return) JuliaCall::julia_eval(src) else {
    JuliaCall::julia_command(src); invisible(NULL)
  }
}

et_require_session <- function() {
  if (!et_is_ready()) {
    stop("no Julia session -- call et_setup() first.", call. = FALSE)
  }
  invisible(TRUE)
}

# Julia string literals need forward slashes on Windows.
julia_path <- function(p) gsub("\\\\", "/", p)

# Push the payload into Main and load the generated module. Returns the module
# name. Idempotent per module name, so re-running a fit in one session is cheap.
et_load_module <- function(gen, keep_file = NULL) {
  et_require_session()
  if (gen$module %in% .et_state$loaded) return(gen$module)
  for (nm in names(gen$payload)) {
    JuliaCall::julia_assign(nm, gen$payload[[nm]])
  }
  f <- keep_file %||% tempfile(pattern = gen$module, fileext = ".jl")
  writeLines(gen$src, f, useBytes = TRUE)
  # `include` at Main scope, so the module lands as `Main.<name>`.
  JuliaCall::julia_command(sprintf('include("%s"); nothing', julia_path(f)))
  .et_state$loaded <- c(.et_state$loaded, gen$module)
  gen$module
}
