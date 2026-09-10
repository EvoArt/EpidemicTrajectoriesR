# Layer 3: the trajectory archive -- and keeping it in a format R can read.
#
# `et_sample(save_x =)` streams the latent trajectory to disc through
# PracticalBayes, whose disk backend is JLD2. That is a Julia format: an R user
# who asked for their trajectories should not be left holding files only Julia
# can open. So the archive is CONVERTED to `.rds` after sampling, in the
# directory they named, and the Julia originals are removed as each one is
# converted.
#
# Conversion is deliberately a separate, exported step rather than something
# buried in `et_sample()`. A fit that crashed part-way through -- as a long
# badger run can -- still leaves every already-flushed chunk complete and
# readable on disc (one file per flush, never appended to), so pointing
# `et_convert_trajectories()` at the template afterwards recovers them. That is
# also why it converts one chunk at a time and only deletes a source file once
# its `.rds` exists: an interrupted conversion loses nothing.
#
# The cost is per FLUSH, not per sweep -- once every `save_every` sweeps against
# thousands of gradient evaluations -- so it is not on any hot path.

# --- chunk discovery ---------------------------------------------------------
#
# PracticalBayes writes one file per flush, named `<base>_iters_<a>_to_<b><ext>`
# (see `chunk_path` in its save_states.jl). Both the converter and the residuals
# reader need to find those files and put them in sweep order, so the logic
# lives here once. Ordering is by FIRST ITERATION as an integer -- sorting the
# names as strings would put `_iters_1000_` before `_iters_101_` and score every
# trajectory against the wrong parameter draw, silently.
et_chunk_files <- function(template) {
  base <- tools::file_path_sans_ext(template)
  ext  <- tools::file_ext(template)
  dir  <- dirname(base)
  if (!nzchar(dir)) dir <- "."
  stem <- basename(base)

  pat <- paste0("^", et_regex_escape(stem), "_iters_(\\d+)_to_(\\d+)\\.",
                et_regex_escape(ext), "$")
  hits <- list.files(dir, pattern = pat)
  if (!length(hits)) return(character(0))

  first <- as.integer(sub(pat, "\\1", hits))
  file.path(dir, hits[order(first)])
}

# Escape a literal string for use inside a regex: a stem or extension may
# contain `.`, which would otherwise match any character. The metacharacters are
# matched as a FIXED set character by character rather than by a regex of their
# own -- writing that pattern is how this went wrong the first time (`{}` inside
# a bracket expression is read as a repetition count by TRE).
et_regex_escape <- function(s) {
  meta <- c(".", "\\", "^", "$", "|", "?", "*", "+", "(", ")", "[", "]", "{", "}")
  chars <- strsplit(s, "", fixed = TRUE)[[1]]
  paste0(ifelse(chars %in% meta, paste0("\\", chars), chars), collapse = "")
}

#' Convert an archived trajectory to R format.
#'
#' [et_sample()] streams the latent trajectory through Julia, whose disk backend
#' writes `.jld2`. This rewrites that archive as `.rds`, which R reads natively,
#' and removes the Julia originals.
#'
#' `et_sample(save_x = "X.rds")` calls this for you once sampling finishes. Call
#' it yourself when a run was interrupted: every flush already on disc is a
#' complete, independently readable file, so a crashed fit's trajectories are
#' still recoverable by pointing this at the template it was writing.
#'
#' Chunks are converted one at a time, and a source file is deleted only once
#' its `.rds` counterpart exists -- so an interrupted conversion loses nothing
#' and can simply be re-run. Only one chunk is ever held in memory.
#'
#' @param template The path given to `et_sample(save_x =)`, with either
#'   extension (e.g. `"X.rds"` or `"X.jld2"`). The per-flush files are found by
#'   their `_iters_x_to_y` names, so the template itself need not exist.
#' @param keep_jld2 Keep the Julia originals instead of deleting them as they
#'   are converted.
#' @param quiet Suppress progress messages.
#' @return Invisibly, the paths of the `.rds` chunks in sweep order.
#' @export
et_convert_trajectories <- function(template, keep_jld2 = FALSE, quiet = FALSE) {
  et_require_session()
  src_template <- et_with_ext(template, "jld2")
  src <- et_chunk_files(src_template)
  if (!length(src)) {
    # Already converted is a no-op, not an error: this is re-runnable by design.
    done <- et_chunk_files(et_with_ext(template, "rds"))
    if (length(done)) {
      if (!quiet) message("Trajectories already in R format (", length(done),
                          " chunks).")
      return(invisible(done))
    }
    stop("et_convert_trajectories(): no archived trajectories found for '",
         template, "'. Looked for files named '",
         basename(tools::file_path_sans_ext(src_template)),
         "_iters_x_to_y.jld2'.", call. = FALSE)
  }

  JuliaCall::julia_command("using JLD2")
  out <- character(length(src))
  for (k in seq_along(src)) {
    dest <- et_with_ext(src[k], "rds")
    # One chunk at a time: a badger flush of 100 sweeps is ~307 MB, and the
    # whole point is never to hold the run at once.
    # Stacked to a 3-D array (time x individual x draw) on the Julia side. A
    # `Vector{Matrix}` comes back as an opaque `JuliaObject` rather than a list
    # of R matrices, so the chunk is marshalled as ONE array -- which crosses as
    # a plain R array -- and split here.
    arr <- JuliaCall::julia_eval(sprintf(paste0(
      "JLD2.jldopen(%s, \"r\") do f; s = f[\"states\"]; ",
      "isempty(s) ? Array{Int}(undef, 0, 0, 0) : ",
      "cat((Array{Int}(x) for x in s)...; dims=3); end"),
      julia_string(julia_path(src[k]))))
    arr <- as.array(arr)
    storage.mode(arr) <- "integer"
    n_draw <- if (length(dim(arr)) == 3L) dim(arr)[3] else 0L
    saveRDS(lapply(seq_len(n_draw), function(j) arr[, , j, drop = FALSE][, , 1]),
            dest)
    rm(arr)
    if (!keep_jld2) unlink(src[k])
    out[k] <- dest
    if (!quiet && (k %% 10 == 0 || k == length(src))) {
      message(sprintf("Converted %d/%d chunks", k, length(src)))
    }
  }
  if (!quiet) {
    message("Trajectories written as .rds: ", dirname(out[1]))
  }
  invisible(out)
}

#' Read archived trajectories back into R.
#'
#' Reads the trajectories [et_sample()] archived, as a list of
#' `n_timepoints x n_individuals` integer matrices of 1-based state codes.
#'
#' Selecting `draws` reads only the chunks that contain them, so inspecting a
#' handful out of thousands does not read the whole archive. Reading them ALL is
#' what the archive exists to avoid -- [et_residuals()] scores them one chunk at
#' a time without ever materialising the run -- so pass `draws` unless the
#' archive is small.
#'
#' @param fit An [et_fit()] from [et_sample()] run with `save_x =`, or a path
#'   template.
#' @param draws Which archived draws to read (indices into the archive, in sweep
#'   order). Defaults to all of them.
#' @return A list of integer matrices.
#' @export
et_trajectories <- function(fit, draws = NULL) {
  template <- if (inherits(fit, "et_fit")) fit$x_path else fit
  if (is.null(template)) {
    stop("et_trajectories(): this fit did not archive its trajectories. ",
         "Re-run et_sample() with `save_x = \"X.rds\"`.", call. = FALSE)
  }
  files <- et_chunk_files(et_with_ext(template, "rds"))
  if (!length(files)) {
    files <- et_chunk_files(et_with_ext(template, "jld2"))
    if (length(files)) {
      stop("et_trajectories(): the archive is still in Julia format. ",
           "Run et_convert_trajectories('", template, "') first.", call. = FALSE)
    }
    stop("et_trajectories(): no archived trajectories found for '", template,
         "'.", call. = FALSE)
  }

  out <- list()
  seen <- 0L
  for (f in files) {
    chunk <- readRDS(f)
    idx <- seen + seq_along(chunk)
    if (is.null(draws)) {
      out <- c(out, chunk)
    } else {
      want <- which(idx %in% draws)
      if (length(want)) out <- c(out, chunk[want])
    }
    seen <- seen + length(chunk)
    # Stop once every requested draw is behind us, so a request for draw 3 does
    # not read a thousand later chunks.
    if (!is.null(draws) && seen >= max(draws)) break
  }
  out
}

# Swap a path's extension. The archive is addressed by template throughout, and
# the two formats differ only here.
et_with_ext <- function(path, ext) {
  paste0(tools::file_path_sans_ext(path), ".", ext)
}
