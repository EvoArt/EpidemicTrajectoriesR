# Layer 3: et_par() / et_model() -- parameters, priors, deterministics, and the
# two likelihood terms.

#' Declare one model parameter.
#'
#' @param prior A prior from [et_priors]. Omit only for `kind = "latent"`.
#' @param init Initial value. Length must match `n` (or `prod(dim)`).
#' @param n Length of a vector parameter; `1` for a scalar. A vector parameter is
#'   emitted as `PracticalBayes.filldist(prior, n)`.
#' @param dim Dimensions of a matrix-valued parameter (`kind = "latent"` only).
#' @param kind `"sampled"` (the default) for a parameter with a prior, or
#'   `"latent"` for one **owned by a conjugate kernel** -- emitted as a
#'   placeholder distribution whose density is constant, exactly as the badger
#'   model's `nu` is.
#' @return An object of class `et_par`.
#' @export
et_par <- function(prior = NULL, init, n = 1L, dim = NULL,
                   kind = c("sampled", "latent")) {
  kind <- match.arg(kind)
  if (kind == "sampled" && !inherits(prior, "et_dist")) {
    stop("et_par(): a sampled parameter needs a `prior` from et_normal(), ",
         "et_gamma(), ... (or et_dist() for anything else).", call. = FALSE)
  }
  if (kind == "latent" && is.null(dim)) dim <- length(init)
  init <- as.numeric(init)
  expected <- if (!is.null(dim)) prod(dim) else n
  if (length(init) != expected) {
    stop("et_par(): `init` has length ", length(init), " but the parameter has ",
         expected, " element(s). A wrong length here becomes a wrong Julia ",
         "indexing convention later, so it is checked now.", call. = FALSE)
  }
  structure(list(prior = prior, init = init, n = as.integer(n),
                 dim = if (is.null(dim)) NULL else as.integer(dim),
                 kind = kind), class = "et_par")
}

#' @export
print.et_par <- function(x, ...) {
  cat("<et_par> ", x$kind,
      if (!is.null(x$prior)) paste0(" ~ ", dist_to_julia(x$prior)) else "",
      "  init=[", paste(signif(x$init, 4), collapse = ", "), "]\n", sep = "")
  invisible(x)
}

#' Assemble the full model: data, parameters, priors and deterministics.
#'
#' @param data An [et_data()] object.
#' @param parameters Named list of [et_par()] declarations.
#' @param derived Named list of `quote()`d expressions over other parameter
#'   names -- PracticalBayes deterministics (`:=`). E.g.
#'   `list(m = quote(m_tilde + 1))`.
#' @param entry_time Optional per-individual entry time: a numeric vector, or the
#'   name of an entry in `data`'s `extras`. Supplying it switches on ET's entry
#'   gate, which scores disease transitions before entry but NOT the survival
#'   factor. Requires the transitions to declare an [et_survival()] -- the same
#'   function is passed through automatically, which is the only way the
#'   subtraction removes exactly what was multiplied in.
#' @return An object of class `et_model`.
#' @export
et_model <- function(data, parameters, derived = list(), entry_time = NULL) {
  if (!inherits(data, "et_data")) {
    stop("et_model(): `data` must come from et_data().", call. = FALSE)
  }
  if (!is.list(parameters) || !length(parameters) || is.null(names(parameters))) {
    stop("et_model(): `parameters` must be a non-empty NAMED list of et_par().",
         call. = FALSE)
  }
  check_julia_name(names(parameters), "parameter name")
  for (nm in names(parameters)) {
    if (!inherits(parameters[[nm]], "et_par")) {
      stop("et_model(): parameter '", nm, "' must be created with et_par().",
           call. = FALSE)
    }
  }
  if ("X" %in% names(parameters)) {
    stop("et_model(): 'X' is reserved for the latent trajectory block.",
         call. = FALSE)
  }
  if (!is.list(derived)) stop("et_model(): `derived` must be a list.", call. = FALSE)
  if (length(derived)) {
    if (is.null(names(derived)) || any(!nzchar(names(derived)))) {
      stop("et_model(): every entry of `derived` must be named.", call. = FALSE)
    }
    check_julia_name(names(derived), "derived name")
    clash <- intersect(names(derived), names(parameters))
    if (length(clash)) {
      stop("et_model(): derived name(s) clash with parameters: ",
           paste(clash, collapse = ", "), ".", call. = FALSE)
    }
    known <- c(names(parameters), names(derived))
    for (nm in names(derived)) {
      unknown <- missing_from(expr_symbols(derived[[nm]]), known)
      unknown <- setdiff(unknown, c(names(.et_builtin_map), "T", "F", "TRUE", "FALSE"))
      if (length(unknown)) {
        stop("et_model(): derived '", nm, "' refers to unknown name(s): ",
             paste(unknown, collapse = ", "),
             ". Only parameters and other derived values are in scope.",
             call. = FALSE)
      }
    }
  }

  if (!is.null(entry_time)) {
    if (is.null(data$transitions$survival)) {
      stop("et_model(): `entry_time` requires the transitions to declare an ",
           "et_survival(). The entry gate removes the survival factor before ",
           "entry, so it has to know which factor that is.", call. = FALSE)
    }
    if (is.character(entry_time)) {
      if (!(entry_time %in% names(data$extras))) {
        stop("et_model(): entry_time '", entry_time, "' is not an entry of ",
             "`extras`. Available: ", paste(names(data$extras), collapse = ", "),
             ".", call. = FALSE)
      }
    } else {
      entry_time <- as.integer(entry_time)
      if (length(entry_time) != data$n_individuals) {
        stop("et_model(): `entry_time` has length ", length(entry_time),
             " but there are ", data$n_individuals, " individuals.", call. = FALSE)
      }
    }
  }

  # The reference set every downstream check uses: which names a parameter
  # expression may legally mention.
  structure(list(data = data, parameters = parameters, derived = derived,
                 entry_time = entry_time,
                 par_names = names(parameters),
                 all_names = c(names(parameters), names(derived))),
            class = "et_model")
}

#' @export
print.et_model <- function(x, ...) {
  cat("<et_model>\n")
  print(x$data)
  cat("  parameters:\n")
  for (nm in x$par_names) {
    p <- x$parameters[[nm]]
    sz <- if (!is.null(p$dim)) paste0("[", paste(p$dim, collapse = "x"), "]")
          else if (p$n > 1L) paste0("[", p$n, "]") else ""
    cat("    ", nm, sz, if (p$kind == "latent") "  (conjugate-owned)"
        else paste0(" ~ ", dist_to_julia(p$prior)), "\n", sep = "")
  }
  if (length(x$derived)) {
    for (nm in names(x$derived)) {
      cat("    ", nm, " := ", deparse(x$derived[[nm]]), "\n", sep = "")
    }
  }
  invisible(x)
}

# Every symbol appearing in an R expression.
expr_symbols <- function(e) {
  if (is.symbol(e)) return(as.character(e))
  if (!is.call(e)) return(character())
  unique(unlist(lapply(as.list(e)[-1L], expr_symbols)))
}

# Expand a set of names to the SAMPLED parameters behind them: a derived value
# resolves to the parameters its expression reads, transitively. This is what
# keeps `depends=` honest when a rate function reads `model$m` and `m` is
# `m_tilde + 1`.
expand_to_sampled <- function(names_in, model) {
  out <- character(); seen <- character()
  frontier <- names_in
  while (length(frontier)) {
    nm <- frontier[1]; frontier <- frontier[-1]
    if (nm %in% seen) next
    seen <- c(seen, nm)
    if (nm %in% model$par_names) { out <- c(out, nm); next }
    if (nm %in% names(model$derived)) {
      frontier <- c(frontier, expr_symbols(model$derived[[nm]]))
    }
    # A name that is neither is not a model variable (a local, a builtin); it
    # contributes no dependency.
  }
  unique(out)
}
