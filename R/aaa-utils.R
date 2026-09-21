# Julia literals from R values. Knows nothing about epidemics or
# PracticalBayes.
#
# Every literal carries its element type explicitly. Letting Julia infer one
# from a promoted mixture gives an abstractly-typed container, which costs
# about 5x downstream and does not show up in inference checks.

# ---- scalars ----------------------------------------------------------------

julia_num <- function(x) {
  if (is.logical(x)) return(if (isTRUE(x)) "true" else "false")
  if (is.na(x)) stop("cannot emit NA as a Julia literal.")
  if (is.integer(x)) return(as.character(x))
  if (!is.finite(x)) {
    return(if (is.nan(x)) "NaN" else if (x > 0) "Inf" else "-Inf")
  }
  # A whole-number double is emitted with a `.0` suffix only where the context
  # wants a Float64; callers that need an index use julia_int().
  #
  # Shortest round-tripping representation: `digits = 17` always round-trips but
  # renders 0.8 as 0.80000000000000004, which makes the generated source much
  # harder to read and to diff. Widen only until the value comes back exactly.
  sci <- abs(x) != 0 && (abs(x) < 1e-4 || abs(x) >= 1e15)
  for (dg in c(15L, 16L, 17L)) {
    s <- format(x, digits = dg, scientific = sci, trim = TRUE)
    if (identical(as.numeric(s), x)) return(s)
  }
  format(x, digits = 17, scientific = sci, trim = TRUE)
}

julia_int <- function(x) {
  x <- as.numeric(x)
  if (any(is.na(x)) || any(x != round(x))) {
    stop("expected whole numbers, got: ", paste(x, collapse = ", "))
  }
  format(x, scientific = FALSE, trim = TRUE, nsmall = 0)
}

# A Float64 literal: always carries a decimal point, so `[1, 2]` cannot come out
# as a Vector{Int} where a Vector{Float64} was meant.
julia_float <- function(x) {
  s <- julia_num(x)
  if (!grepl("[.eE]", s) && !s %in% c("Inf", "-Inf", "NaN", "true", "false")) {
    s <- paste0(s, ".0")
  }
  s
}

julia_string <- function(x) paste0('"', gsub('"', '\\\\"', x, fixed = TRUE), '"')

julia_symbol <- function(x) paste0(":", x)

julia_symbol_tuple <- function(x) {
  if (length(x) == 0) return("()")
  if (length(x) == 1) return(paste0("(", julia_symbol(x), ",)"))
  paste0("(", paste(vapply(x, julia_symbol, character(1)), collapse = ", "), ")")
}

julia_symbol_vector <- function(x) {
  paste0("[", paste(vapply(x, julia_symbol, character(1)), collapse = ", "), "]")
}

# ---- arrays -----------------------------------------------------------------

# Emit a numeric vector as a typed Julia vector literal. `type` is "Int",
# "Float64" or "Bool"; giving it explicitly is what keeps the element type
# concrete regardless of how R happened to store the value.
julia_vector <- function(x, type = NULL) {
  if (is.null(type)) type <- r_julia_type(x)
  fmt <- switch(type, Int = julia_int, Float64 = julia_float, Bool = julia_num,
                stop("unsupported element type '", type, "'."))
  if (length(x) == 0) return(sprintf("%s[]", type))
  sprintf("%s[%s]", type, paste(vapply(x, fmt, character(1)), collapse = ", "))
}

# The Julia element type an R vector should become. R has no integer/double
# distinction at the literal level, so `storage.mode` is the only signal and a
# double that happens to hold whole numbers stays Float64 -- the caller must say
# otherwise if it wants indices.
r_julia_type <- function(x) {
  if (is.logical(x)) return("Bool")
  if (is.integer(x)) return("Int")
  if (is.double(x)) return("Float64")
  stop("cannot map R type '", class(x)[1], "' to a Julia element type.")
}

# Emit an R matrix/array as a Julia array literal. Column-major in both
# languages, so `as.vector()` is already in Julia's storage order and only the
# reshape has to be written out.
julia_array <- function(x, type = NULL) {
  if (is.null(type)) type <- r_julia_type(as.vector(x))
  d <- dim(x)
  if (is.null(d)) return(julia_vector(as.vector(x), type))
  sprintf("reshape(%s, %s)", julia_vector(as.vector(x), type),
          paste(julia_int(d), collapse = ", "))
}

# ---- misc -------------------------------------------------------------------

indent <- function(s, n = 4L) {
  if (!nzchar(s)) return("")
  pad <- strrep(" ", n)
  lines <- strsplit(s, "\n", fixed = TRUE)[[1]]
  paste0(pad, lines, collapse = "\n")
}

# Join non-empty source chunks with blank lines between them.
join_blocks <- function(...) {
  parts <- unlist(list(...), use.names = FALSE)
  parts <- parts[!vapply(parts, function(p) is.null(p) || !nzchar(p), logical(1))]
  paste(parts, collapse = "\n\n")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# A name that is safe as a Julia identifier. Users get told rather than silently
# renamed, because a silently-renamed parameter would not match the draws coming
# back.
check_julia_name <- function(x, what = "name") {
  bad <- x[!grepl("^[A-Za-z_][A-Za-z0-9_]*$", x)]
  if (length(bad)) {
    stop(what, " must be a valid Julia identifier (letters, digits, underscore, ",
         "not starting with a digit). Offending: ",
         paste(sQuote(bad), collapse = ", "), ".", call. = FALSE)
  }
  invisible(x)
}

# Structural equality helper used by validation messages.
missing_from <- function(x, allowed) setdiff(x, allowed)
