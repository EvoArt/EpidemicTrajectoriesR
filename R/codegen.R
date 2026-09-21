# Turns an et_model plus its blocks into a self-contained Julia module.
#
# Source text rather than call-by-call API driving, because three of the four
# ET/PB declaration points are macros and JuliaCall cannot reach those. The
# generated module also compiles once and stays concretely typed.
#
# et_julia_source() returns the source; it runs standalone.

#' Generate the Julia source for a model.
#'
#' Works without a Julia session: useful for inspection, for review, and as the
#' on-ramp to writing the model in Julia directly.
#'
#' @param model An [et_model()].
#' @param blocks List of sampler blocks; see [et_nuts()] and friends. Defaults
#'   are filled in exactly as [et_sample()] would.
#' @param module_name Julia module name. Generated if omitted.
#' @return An object of class `et_source`: `$src` is the Julia text, `$payload`
#'   the R arrays that must be assigned into Julia before it is loaded.
#' @export
et_julia_source <- function(model, blocks = list(), module_name = NULL) {
  if (!inherits(model, "et_model")) {
    stop("et_julia_source(): `model` must come from et_model().", call. = FALSE)
  }
  blocks <- resolve_blocks(model, blocks)
  # Generate once against a placeholder, hash the result, and only then fix the
  # module name. Identical models therefore produce identical source and reuse
  # an already-compiled module, instead of recompiling the whole stack on every
  # et_loglik()/et_check_depends()/et_sample() call.
  gen <- generate_module(model, blocks, MODULE_PLACEHOLDER)
  if (is.null(module_name)) module_name <- content_module_name(gen)
  gen <- rename_module(gen, module_name)
  structure(gen, class = "et_source")
}

#' @export
print.et_source <- function(x, ...) {
  cat(x$src)
  invisible(x)
}

#' @export
as.character.et_source <- function(x, ...) x$src

MODULE_PLACEHOLDER <- "ETRGenXXPLACEHOLDERXX"

# A name that is a function of the model, so the same model always lands in the
# same Julia module.
content_module_name <- function(gen) {
  f <- tempfile(fileext = ".bin")
  on.exit(unlink(f), add = TRUE)
  con <- file(f, "wb")
  serialize(list(src = gen$src, payload = gen$payload), con)
  close(con)
  sprintf("ETRGen_%s", substr(unname(tools::md5sum(f)), 1, 12))
}

rename_module <- function(gen, module_name) {
  gen$src <- gsub(MODULE_PLACEHOLDER, module_name, gen$src, fixed = TRUE)
  names(gen$payload) <- gsub(MODULE_PLACEHOLDER, module_name,
                             names(gen$payload), fixed = TRUE)
  gen$module <- module_name
  gen
}

# ---- the generator ----------------------------------------------------------

generate_module <- function(model, blocks, module_name) {
  d <- model$data
  states <- d$transitions$states
  prefix <- paste0(module_name, "_")
  payload <- list()

  # A payload entry is an R array assigned into Julia's `Main` before the module
  # loads; the module re-binds it with an explicit element type, so nothing
  # downstream depends on how JuliaCall happened to marshal it.
  add_payload <- function(name, value) {
    payload[[paste0(prefix, name)]] <<- value
    julia_typed_binding(name, paste0("Main.", prefix, name), value)
  }

  helpers <- d$helpers
  helper_src <- et_transpile_helpers(helpers)

  # --- leaf functions, and the reads that drive `depends=` -------------------
  epi_reads <- character(); epi_opaque <- FALSE
  obs_reads <- character(); obs_opaque <- FALSE
  fn_src <- character()

  rate_names <- character()
  for (tr in d$transitions$transitions) {
    nm <- sprintf("et_rate_%s_%s", tr$from, tr$to)
    r <- et_transpile_role(tr$rate, "rate", nm, helpers)
    fn_src <- c(fn_src, r$src)
    epi_reads <- union(epi_reads, r$reads); epi_opaque <- epi_opaque || r$opaque
    rate_names[paste0(tr$from, "->", tr$to)] <- nm
  }
  surv_name <- NULL
  if (!is.null(d$transitions$survival)) {
    surv_name <- "et_survival_fn"
    r <- et_transpile_role(d$transitions$survival$fn, "survival", surv_name, helpers)
    fn_src <- c(fn_src, r$src)
    epi_reads <- union(epi_reads, r$reads); epi_opaque <- epi_opaque || r$opaque
  }
  r <- et_transpile_role(d$starting_state, "starting_state", "et_starting_state",
                         helpers)
  fn_src <- c(fn_src, r$src)
  epi_reads <- union(epi_reads, r$reads); epi_opaque <- epi_opaque || r$opaque
  # The aggregate update can read parameters too; those reach the likelihood
  # through the arrays the rates then read, so they belong to the epidemic term.
  epi_reads <- union(epi_reads, d$aggregates$reads)

  obs_process_name <- NULL
  if (!is.null(d$observation_process)) {
    obs_process_name <- "et_obs_process"
    r <- et_transpile_role(d$observation_process, "obs_process", obs_process_name,
                           helpers)
    fn_src <- c(fn_src, r$src)
    obs_reads <- union(obs_reads, r$reads); obs_opaque <- obs_opaque || r$opaque
  }
  obs_weight_name <- NULL
  if (!is.null(d$observation_weight)) {
    obs_weight_name <- "et_obs_weight"
    r <- et_transpile_role(d$observation_weight, "obs_weight", obs_weight_name,
                           helpers)
    fn_src <- c(fn_src, r$src)
    obs_reads <- union(obs_reads, r$reads); obs_opaque <- obs_opaque || r$opaque
  }

  # The likelihood's own observation factor, when it differs from the filter's.
  # Its reads, not the filter weight's, are what drive the observation term's
  # `depends=`: the filter is never differentiated, so a parameter that appears
  # only there contributes nothing to any gradient.
  lik_weight_name <- NULL
  if (!is.null(d$likelihood_weight)) {
    lik_weight_name <- "et_obs_lik_weight"
    r <- et_transpile_role(d$likelihood_weight, "obs_weight", lik_weight_name,
                           helpers)
    fn_src <- c(fn_src, r$src)
    obs_reads <- r$reads          # replaces, not unions
    obs_opaque <- r$opaque
  }

  # A user-supplied coupling term replaces the whole default. It is an
  # et_julia() body, so it needs wrapping in ET's `rest_contribution` signature
  # Splicing the body straight into the keyword would put statements inside a
  # call's parentheses.
  #
  # `affected_override` is part of that signature and is usually ignored; the
  # default is what lets ET call it with six arguments or seven.
  if (inherits(d$rest_contribution, "et_julia")) {
    fn_src <- c(fn_src, sprintf(
      "function et_rest_contribution(model, data, X, i, t, n_states, affected_override=nothing)\n%s\nend",
      indent(d$rest_contribution$src)))
  }

  # A separate coupling spec: the rates ET evaluates for the coupling term only,
  # which is never differentiated. Its reads are not added to `depends`, because
  # `epidemic_loglik` never sees it.
  coupling_names <- NULL
  if (!is.null(d$coupling_transitions)) {
    coupling_names <- character()
    for (tr in d$coupling_transitions$transitions) {
      nm <- sprintf("et_couprate_%s_%s", tr$from, tr$to)
      rr <- et_transpile_role(tr$rate, "rate", nm, helpers)
      fn_src <- c(fn_src, rr$src)
      coupling_names[paste0(tr$from, "->", tr$to)] <- nm
    }
    if (!is.null(d$coupling_transitions$survival)) {
      rr <- et_transpile_role(d$coupling_transitions$survival$fn, "survival",
                              "et_coupsurvival_fn", helpers)
      fn_src <- c(fn_src, rr$src)
    }
  }

  # --- payload bindings ------------------------------------------------------
  bindings <- character()
  for (nm in names(d$extras)) bindings <- c(bindings, add_payload(nm, d$extras[[nm]]))
  if (!is.null(d$group)) bindings <- c(bindings, add_payload("group", d$group))
  if (!is.null(d$sampling_period)) {
    bindings <- c(bindings, add_payload("sampling_period_raw", d$sampling_period))
    bindings <- c(bindings, sprintf(
      "const sampling_period = [(sampling_period_raw[i, 1], sampling_period_raw[i, 2]) for i in 1:%s]",
      julia_int(d$n_individuals)))
  }
  if (inherits(d$affected_individuals, "et_affected")) {
    bindings <- c(bindings,
                  add_payload("affected_groups", d$affected_individuals$group_matrix),
                  affected_builder_src(d$n_timepoints))
  }
  if (is.character(model$entry_time)) {
    entry_expr <- sprintf("Vector{Int}(%s)", model$entry_time)
  } else if (!is.null(model$entry_time)) {
    bindings <- c(bindings, add_payload("entry_time", model$entry_time))
    entry_expr <- "entry_time"
  } else {
    entry_expr <- NULL
  }

  # --- spec blocks -----------------------------------------------------------
  states_src <- sprintf("const STATES = %s", julia_symbol_vector(states))
  trans_src <- transitions_src("TRANS", d$transitions, rate_names, surv_name)
  coupling_src <- if (!is.null(coupling_names)) {
    transitions_src("TRANS_COUPLING", d$coupling_transitions, coupling_names,
                    if (!is.null(d$coupling_transitions$survival))
                      "et_coupsurvival_fn" else NULL)
  } else NULL
  # When any tracked array carries a hand-written reverse, the `@aggregate` macro
  # cannot be used at all, so every array is emitted through ET's verbose
  # fallback: a plain NamedTuple plus a tuple of reversible summary functions.
  agg <- d$aggregates
  summary_names <- character()
  summary_src <- character()
  if (length(agg$hand_names %||% character())) {
    for (nm in agg$derived_names) {
      summary_src <- c(summary_src, derived_summary_src(nm, agg$records,
                                                        names(agg$arrays)))
      summary_names <- c(summary_names, paste0("et_summary_", nm))
    }
    for (nm in agg$hand_names) {
      r <- et_transpile_summary(nm, agg$arrays[[nm]], states, helpers,
                                names(agg$arrays))
      summary_src <- c(summary_src, r$src)
      summary_names <- c(summary_names, r$name)
      epi_reads <- union(epi_reads, r$reads)
      epi_opaque <- epi_opaque || r$opaque
    }
  }
  aggs_src <- aggregate_src(agg)
  data_src <- data_src_for(d, obs_process_name, obs_weight_name,
                           !is.null(coupling_names), summary_names)

  loglik_src <- if (!is.null(entry_expr)) {
    sprintf("const LOGLIK = epidemic_loglik(DATA; entry_time=%s, survival=%s)",
            entry_expr, surv_name)
  } else {
    "const LOGLIK = epidemic_loglik(DATA)"
  }
  has_obs <- !is.null(obs_process_name) || !is.null(obs_weight_name)
  obsll_src <- if (!has_obs) NULL
    else if (!is.null(lik_weight_name)) sprintf(
      "const OBSLOGLIK = epidemic_obs_loglik(DATA; observation_weight=%s)",
      lik_weight_name)
    else "const OBSLOGLIK = epidemic_obs_loglik(DATA)"
  latent_src <- sprintf("const LATENT! = epidemic_latent_sampler(DATA%s)",
                        if (isTRUE(traj_block(blocks)$mh)) "; mh=true" else "")

  # --- the PracticalBayes model ---------------------------------------------
  dep_epi <- if (epi_opaque) NULL
             else c(expand_to_sampled(epi_reads, model), "X")
  dep_obs <- if (obs_opaque) NULL
             else c(expand_to_sampled(obs_reads, model), "X")
  model_src <- pb_model_src(model, has_obs, dep_epi, dep_obs)
  placeholder_src <- placeholder_defs(model)

  # --- sampler and runner ----------------------------------------------------
  gibbs_src <- gibbs_src_for(model, blocks)
  run_src <- run_src_for(model, has_obs)
  population_src <- population_src_for(model, has_obs)

  # no timestamp. The module name is a hash of this text, so a header that
  # varied with the clock would make every generation a fresh module and force a
  # full recompile of the stack on each call.
  header <- sprintf(paste0(
    "# Generated by EpidemicTrajectoriesR. DO NOT EDIT BY HAND.\n",
    "#\n",
    "# This file is self-contained: it runs as an ordinary Julia script once the\n",
    "# payload arrays named `%s*` exist in `Main`. Editing the R model and\n",
    "# regenerating is the supported workflow; editing here is the on-ramp to\n",
    "# writing the model in Julia directly.\n"), prefix)

  src <- paste0(
    header, "\n",
    sprintf("module %s\n", module_name),
    join_blocks(
      julia_usings(),
      section("payload from R", bindings),
      section("state space", states_src),
      section("user functions", c(helper_src, fn_src)),
      section("transitions", c(trans_src, coupling_src)),
      section("aggregates", c(summary_src, aggs_src)),
      section("data", data_src),
      section("generated artefacts", c(loglik_src, obsll_src, latent_src)),
      section("placeholders for conjugate-owned parameters", placeholder_src),
      section("the model", model_src),
      section("sampler", gibbs_src),
      section("runner", run_src),
      section("population", population_src),
      section("depends= validator", check_src_for(model, has_obs)),
      section("parameter draws from R", unflatten_src(model)),
      extract_helper_src()
    ),
    sprintf("\n\nend # module %s\n", module_name))

  list(src = src, payload = payload, module = module_name, blocks = blocks,
       depends = list(epidemic = dep_epi, observation = dep_obs),
       model = model)
}

section <- function(title, body) {
  body <- body[!vapply(body, function(b) is.null(b) || !nzchar(b), logical(1))]
  if (!length(body)) return(NULL)
  paste0("# ---- ", title, " ",
         strrep("-", max(3, 74 - nchar(title))), "\n",
         paste(body, collapse = "\n\n"))
}

julia_usings <- function() paste(
  "using EpidemicTrajectories",
  "using PracticalBayes",
  "using PracticalEpiBayes",
  "using Distributions",
  "using Random: AbstractRNG",
  # PracticalBayes re-exports HMC/NUTS but not the integrator or the metric, and
  # a hand-tuned diagonal metric is the badger model's whole HMC configuration.
  "using AdvancedHMC: Leapfrog, DiagEuclideanMetric",
  "using StableRNGs: StableRNG",
  "import AbstractMCMC",
  # The AD backend is a runtime argument to et_run, not a compile-time constant,
  # so switching backend does not force the module to be regenerated.
  "import ADTypes",
  sep = "\n")

# ---- payload typing ---------------------------------------------------------

# Bind an R-marshalled value under an explicit Julia type. Concrete types are
# the first performance rule in the Julia stack, and JuliaCall's own conversion is not
# something to depend on: an R integer may arrive as Int32, a 1-column matrix as
# a vector. Converting here makes the module's types a property of the generated
# source rather than of the marshalling.
julia_typed_binding <- function(name, expr, value) {
  el <- if (is.integer(value)) "Int" else if (is.logical(value)) "Bool" else "Float64"
  d <- dim(value)
  if (is.null(d)) {
    if (length(value) == 1L) {
      return(sprintf("const %s = %s(%s)", name, el, expr))
    }
    return(sprintf("const %s = Vector{%s}(vec(%s))", name, el, expr))
  }
  if (length(d) == 2L) {
    return(sprintf("const %s = Matrix{%s}(%s)", name, el, expr))
  }
  sprintf("const %s = Array{%s,%d}(%s)", name, el, length(d), expr)
}

affected_builder_src <- function(n_timepoints) {
  sprintf(paste0(
    "# Who each individual affects at each time, from time-varying group\n",
    "# membership. The focal is excluded from its own list, matching ET's default.\n",
    "const affected_individuals = let gm = affected_groups, nt = %s\n",
    "    m = size(gm, 1)\n",
    "    ng = maximum(gm)\n",
    "    by_group = [Int[] for _ in 1:ng, _ in 1:nt]\n",
    "    for t in 1:nt, i in 1:m\n",
    "        g = gm[i, t]\n",
    "        g > 0 && push!(by_group[g, t], i)\n",
    "    end\n",
    "    out = Matrix{Vector{Int}}(undef, nt, m)\n",
    "    for t in 1:nt, i in 1:m\n",
    "        g = gm[i, t]\n",
    "        out[t, i] = g == 0 ? Int[] : [j for j in by_group[g, t] if j != i]\n",
    "    end\n",
    "    out\n",
    "end"), julia_int(n_timepoints))
}

# ---- spec blocks ------------------------------------------------------------

transitions_src <- function(const_name, spec, rate_names, surv_name) {
  lines <- character()
  if (!is.null(spec$survival)) {
    lines <- c(lines, sprintf("@survival %s death=%s", surv_name,
                              julia_symbol(spec$survival$death)))
  }
  for (tr in spec$transitions) {
    lines <- c(lines, sprintf("%s -> %s = %s", tr$from, tr$to,
                              rate_names[[paste0(tr$from, "->", tr$to)]]))
  }
  sprintf("const %s = @transitions STATES%s begin\n%s\nend", const_name,
          if (spec$auto_self) "" else " :no_auto_self",
          indent(paste(lines, collapse = "\n")))
}

# Emit the aggregates.
#
# Two shapes, decided by whether any array carries a hand-written reverse:
#
#   * all derived  -> the `@aggregate` macro, which allocates the arrays and
#     writes both directions of each update from its `+=` shape.
#   * any hand-written -> ET's verbose fallback: a plain NamedTuple of arrays
#     plus a tuple of `(model, data, X, s, i, t, reverse)` functions. The macro
#     cannot be mixed with it, so when one array needs a hand-written reverse
#     every array is emitted this way and the derived ones get generated
#     functions of the same shape.
aggregate_src <- function(aggs, summary_fns = list()) {
  if (!length(aggs$hand_names %||% character())) {
    decls <- vapply(names(aggs$arrays), function(nm) {
      a <- aggs$arrays[[nm]]
      sprintf("@array %s %s (%s)", nm, a$type,
              paste(julia_int(a$dim), collapse = ", "))
    }, character(1))
    return(sprintf("const AGGS = @aggregate STATES begin\n%s\nend",
                   indent(paste(c(decls, aggs$lines), collapse = "\n"))))
  }

  alloc <- vapply(names(aggs$arrays), function(nm) {
    a <- aggs$arrays[[nm]]
    sprintf("%s = zeros(%s, %s)", nm, a$type,
            paste(julia_int(a$dim), collapse = ", "))
  }, character(1))
  sprintf("const AGGS = (; %s)", paste(alloc, collapse = ", "))
}

# The generated counterpart of the `@aggregate` lines for one array, used on the
# fallback path: a reversible summary function of the shape ET calls.
#
# Built from the parsed records rather than by rewriting the rendered text, so a
# guard stays a guard and the inverse operator is chosen structurally (`+` -> `-`,
# `*` -> `/`) rather than by string substitution.
derived_summary_src <- function(name, records, array_names) {
  recs <- Filter(function(r) r$array == name, records)
  body <- vapply(recs, function(r) {
    inv <- if (r$op == "+") "-" else "/"
    tgt <- qualify_aggregates(r$target, array_names)
    # `if reverse ... else ... end` per update, so several updates to the same
    # array each keep their own guard.
    core <- sprintf("if reverse\n    %s %s= %s\nelse\n    %s %s= %s\nend",
                    tgt, inv, r$contrib, tgt, r$op, r$contrib)
    if (is.null(r$guard)) core
    else sprintf("if %s\n%s\nend", r$guard, indent(core))
  }, character(1))
  sprintf(paste0(
    "function et_summary_%s(model, data, X, state, i, t, reverse=false)\n",
    "%s\n    nothing\nend"),
    name, indent(paste(body, collapse = "\n")))
}

# A tracked array is written as a bare name in an update body. That is the
# vocabulary `@aggregate` establishes, and the macro rewrites it to
# `data.aggregates.<name>` before ET ever sees it. A hand-written summary reaches
# ET directly, so the same rewrite has to happen here.
#
# Word-boundary anchored so an array called `n` does not corrupt an unrelated
# identifier that merely contains it.
qualify_aggregates <- function(src, array_names) {
  for (nm in array_names) {
    src <- gsub(sprintf("(?<![A-Za-z0-9_.])%s(?=\\s*\\[)", nm),
                sprintf("data.aggregates.%s", nm), src, perl = TRUE)
  }
  src
}

data_src_for <- function(d, obs_process_name, obs_weight_name, has_coupling,
                         summary_names = character()) {
  kw <- c(
    sprintf("n_individuals=%s", julia_int(d$n_individuals)),
    sprintf("n_timepoints=%s", julia_int(d$n_timepoints)),
    "trans_mat=TRANS",
    if (has_coupling) "coupling_trans_mat=TRANS_COUPLING",
    "starting_state=et_starting_state",
    "aggregates=AGGS",
    # ET's verbose fallback: a plain NamedTuple of arrays needs its updates
    # supplied explicitly. A Tuple, not a Vector, so each summary keeps its own
    # concrete type through the hot loop.
    if (length(summary_names)) sprintf("derived_summaries=(%s,)",
                                       paste(summary_names, collapse = ", ")),
    if (!is.null(d$group)) "group=group",
    if (!is.null(obs_process_name)) sprintf("observation_process=%s", obs_process_name),
    if (!is.null(obs_weight_name)) sprintf("observation_weight=%s", obs_weight_name),
    if (!is.null(d$sampling_period)) "sampling_period=sampling_period",
    if (inherits(d$affected_individuals, "et_affected"))
      "affected_individuals=affected_individuals"
    else if (inherits(d$affected_individuals, "et_julia"))
      sprintf("affected_individuals=%s", d$affected_individuals$src),
    if (!is.null(d$coupled_transitions)) sprintf(
      "coupled_transitions=[%s]",
      paste(vapply(d$coupled_transitions, function(p)
        sprintf("(%s, %s)", julia_symbol(p[1]), julia_symbol(p[2])), character(1)),
        collapse = ", ")),
    if (inherits(d$rest_contribution, "et_julia")) "rest_contribution=et_rest_contribution",
    if (!d$focal_self_contribution) "focal_self_contribution=false",
    vapply(names(d$extras), function(nm) sprintf("%s=%s", nm, nm), character(1))
  )
  kw <- kw[!vapply(kw, is.null, logical(1))]
  sprintf("const DATA = epidemic_data(;\n%s,\n)",
          indent(paste(kw, collapse = ",\n")))
}

# ---- the PracticalBayes model ----------------------------------------------

placeholder_defs <- function(model) {
  out <- character()
  for (nm in model$par_names) {
    p <- model$parameters[[nm]]
    if (p$kind != "latent") next
    if (length(p$dim) != 2L) {
      stop("et_par(): kind = \"latent\" is currently supported only for ",
           "matrix-valued parameters (dim of length 2). Parameter '", nm,
           "' has dim of length ", length(p$dim),
           ". A scalar or vector parameter sampled by a conjugate kernel should ",
           "keep its real prior -- the kernel overrides it, exactly as the cattle ",
           "model's nu does.", call. = FALSE)
    }
    ty <- paste0("_ETPlaceholder_", nm)
    out <- c(out, sprintf(paste0(
      "# `%s` is owned by a conjugate kernel. Its density is a constant, so the\n",
      "# placeholder only has to supply a VALID value: PracticalBayes evaluates\n",
      "# the model body with `rand` during Gibbs-coverage validation, before the\n",
      "# kernel or `init` takes over.\n",
      "struct %s <: Distributions.DiscreteMatrixDistribution end\n",
      "Base.size(::%s) = (%s, %s)\n",
      "Distributions.logpdf(::%s, ::AbstractMatrix) = 0.0\n",
      "Distributions.rand(::AbstractRNG, ::%s) = %s"),
      nm, ty, ty, julia_int(p$dim[1]), julia_int(p$dim[2]), ty, ty,
      julia_array(matrix(p$init, p$dim[1], p$dim[2]), "Float64")))
  }
  out
}

par_dist_src <- function(nm, p) {
  if (p$kind == "latent") return(sprintf("_ETPlaceholder_%s()", nm))
  if (p$n > 1L) {
    return(sprintf("PracticalBayes.filldist(%s, %s)", dist_to_julia(p$prior),
                   julia_int(p$n)))
  }
  dist_to_julia(p$prior)
}

pb_model_src <- function(model, has_obs, dep_epi, dep_obs) {
  lines <- character()
  for (nm in model$par_names) {
    lines <- c(lines, sprintf("%s ~ %s", nm, par_dist_src(nm, model$parameters[[nm]])))
  }
  for (nm in names(model$derived)) {
    lines <- c(lines, sprintf("%s := %s", nm,
                              et_transpile_expr(model$derived[[nm]], new_ctx())))
  }
  lines <- c(lines, sprintf("X ~ TrajectoryLatent(%s, %s)",
                            julia_int(model$data$n_timepoints),
                            julia_int(model$data$n_individuals)))
  pars <- c(model$par_names, names(model$derived))
  lines <- c(lines, sprintf("pars = (; %s)",
                            paste(sprintf("%s=%s", pars, pars), collapse = ", ")))

  dep <- function(v) if (is.null(v)) "" else paste0("  depends=", julia_symbol_tuple(v))
  lines <- c(lines, sprintf("@addlogprob! loglik_fn(pars, data, X)%s", dep(dep_epi)))
  if (has_obs) {
    lines <- c(lines, sprintf("@addlogprob! obs_loglik_fn(pars, data, X)%s",
                              dep(dep_obs)))
  }
  note <- paste0(
    "# `depends=` is DERIVED, not declared: every transpiled body records which\n",
    "# `model$x` it reads, transitively through helpers, and the union is what\n",
    "# appears below. Under-declaring would silently drop a gradient\n",
    "# contribution, so it is never left to hand maintenance. Check it with\n",
    "# et_check_depends().\n")
  paste0(note,
         sprintf("@model function et_the_model(data, loglik_fn%s)\n%s\nend",
                 if (has_obs) ", obs_loglik_fn" else "",
                 indent(paste(lines, collapse = "\n"))))
}

# ---- sampler ----------------------------------------------------------------

traj_block <- function(blocks) {
  blocks[[which(vapply(blocks, function(b) b$kind == "iffbs", logical(1)))[1]]]
}

gibbs_src_for <- function(model, blocks) {
  entries <- vapply(blocks, function(b) block_entry(model, b), character(1))
  paste0(
    params_closure_src(model), "\n\n",
    sprintf("const SPL = Gibbs(\n%s,\n)", indent(paste(entries, collapse = ",\n"))))
}

# The one piece that varies between models: how the Gibbs state maps to the
# parameter NamedTuple the rate functions expect (PracticalEpiBayes takes it as
# a closure precisely so the reparameterisations live outside the kernel).
params_closure_src <- function(model) {
  subs <- stats::setNames(paste0("v.", model$par_names), model$par_names)
  parts <- c(
    sprintf("%s=v.%s", model$par_names, model$par_names),
    vapply(names(model$derived), function(nm)
      sprintf("%s=%s", nm,
              et_transpile_expr(subst_symbols(model$derived[[nm]], subs), new_ctx())),
      character(1)))
  sprintf("et_params(v) = (; %s)", paste(parts, collapse = ", "))
}

block_vars_tuple <- function(vars) {
  if (length(vars) == 1L) return(julia_symbol(vars))
  julia_symbol_tuple(vars)
}

block_entry <- function(model, b) {
  key <- block_vars_tuple(b$vars)
  switch(b$kind,
    nuts = sprintf("%s => NUTS(%s)", key, julia_float(b$target_accept)),
    hmc = sprintf("%s => %s", key, hmc_kernel_src(model, b)),
    adaptive_hmc = sprintf(
      "%s => AdaptiveHMC(%s; n_leapfrog=%s, metric=%s)",
      key, julia_float(b$target_accept), julia_int(b$n_steps),
      julia_symbol(b$metric)),
    iffbs = sprintf("%s => iffbs_kernel(LATENT!; params=et_params)", key),
    conjugate = sprintf("%s => %s", key, conjugate_src(model, b)),
    test_sensitivity = sprintf(
      "%s => test_sensitivity_kernel(%s; Y=%s, infected_state=%s, prior=(%s, %s))",
      key, julia_symbol(b$vars), b$y, state_index(model, b$infected_state),
      julia_float(b$prior[1]), julia_float(b$prior[2])),
    capture_prob = sprintf(paste0(
      "%s => capture_prob_kernel(%s; caught=%s, effort=%s, group=%s, index=%s,\n",
      "        dead_state=%s, n=%s, prior=(%s, %s))"),
      key, julia_symbol(b$vars), b$caught, b$effort, b$group, b$index,
      state_index(model, b$dead_state), julia_int(b$n),
      julia_float(b$prior[1]), julia_float(b$prior[2])),
    initial_state = sprintf(
      "%s => initial_state_kernel(%s; at=%s, eligible=%s, states=(%s), prior=%s, n=%s)",
      key, julia_symbol(b$vars), initial_at_src(b),
      eligible_src(b$eligible),
      paste(vapply(b$states, function(s) state_index(model, s), character(1)),
            collapse = ", "),
      julia_vector(b$prior, "Float64"), julia_int(b$n)),
    julia = sprintf("%s => %s", key, b$src),
    stop("unknown block kind '", b$kind, "'.", call. = FALSE))
}

hmc_kernel_src <- function(model, b) {
  if (is.null(b$step_size)) {
    return(sprintf("HMC(%s)", julia_int(b$n_steps)))
  }
  eps <- expand_step_size(model, b)
  sprintf(paste0("HMC(%s; integrator=Leapfrog(1.0),\n",
                 "        metric=DiagEuclideanMetric(%s .^ 2))"),
          julia_int(b$n_steps), julia_vector(eps, "Float64"))
}

# The metric is a flat vector over the block's parameters, in declaration order,
# so a named list has to be expanded to each parameter's length.
expand_step_size <- function(model, b) {
  ss <- b$step_size
  if (is.numeric(ss) && length(ss) == 1L) {
    total <- sum(vapply(b$vars, function(v) par_length(model, v), numeric(1)))
    return(rep(ss, total))
  }
  if (is.numeric(ss)) {
    total <- sum(vapply(b$vars, function(v) par_length(model, v), numeric(1)))
    if (length(ss) != total) {
      stop("et_hmc(): `step_size` has length ", length(ss),
           " but the block has ", total, " scalar parameters. Give one value, ",
           "a full-length vector, or a named list.", call. = FALSE)
    }
    return(ss)
  }
  if (!is.list(ss) || is.null(names(ss))) {
    stop("et_hmc(): `step_size` must be a number, a numeric vector, or a NAMED ",
         "list with one entry per parameter in the block.", call. = FALSE)
  }
  missing <- setdiff(b$vars, names(ss))
  if (length(missing)) {
    stop("et_hmc(): `step_size` is missing entries for: ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  unlist(lapply(b$vars, function(v) {
    n <- par_length(model, v); val <- as.numeric(ss[[v]])
    if (length(val) == 1L) rep(val, n)
    else if (length(val) == n) val
    else stop("et_hmc(): step_size$", v, " has length ", length(val),
              " but '", v, "' has ", n, " element(s).", call. = FALSE)
  }), use.names = FALSE)
}

par_length <- function(model, nm) {
  p <- model$parameters[[nm]]
  if (is.null(p)) stop("unknown parameter '", nm, "'.", call. = FALSE)
  if (!is.null(p$dim)) prod(p$dim) else p$n
}

state_index <- function(model, state_name) {
  states <- model$data$transitions$states
  k <- match(state_name, states)
  if (is.na(k)) {
    stop("'", state_name, "' is not a state. States are: ",
         paste(states, collapse = ", "), ".", call. = FALSE)
  }
  julia_int(k)
}

initial_at_src <- function(b) {
  if (is.character(b$at)) return(b$at)
  if (length(b$at) == 1L) return(julia_int(b$at))
  julia_vector(as.integer(b$at), "Int")
}

eligible_src <- function(eligible) {
  if (is.null(eligible)) return("(X, data, i, t) -> true")
  if (inherits(eligible, "et_julia")) return(eligible$src)
  check_protocol_formals(eligible, c("X", "data", "i", "t"),
                         "eligible predicate", "eligible")
  ctx <- new_ctx()
  sprintf("(X, data, i, t) -> begin\n%s\nend",
          indent(transpile_braced_body(body(eligible), ctx)))
}

conjugate_src <- function(model, b) {
  fam <- if (b$family == "beta") "Beta" else "Dirichlet"
  prior <- if (b$family == "beta") {
    sprintf("Beta(%s, %s)", julia_float(b$prior[1]), julia_float(b$prior[2]))
  } else {
    sprintf("Dirichlet(%s)", julia_vector(b$prior, "Float64"))
  }
  lhs <- if (b$n > 1L) sprintf("%s[1:%s]", b$vars, julia_int(b$n)) else b$vars
  sprintf("@conjugate %s ~ %s begin\n%s\nend", lhs, prior,
          indent(transpile_count_body(b$count, b$family)))
}

# The count body's terminal `c(...)` becomes the shape the family expects:
# a tuple `(successes, failures)` for Beta, a vector of category counts for
# Dirichlet. Requiring the explicit `c(...)` is what lets the error be specific.
transpile_count_body <- function(f, family) {
  b <- body(f)
  stmts <- if (is.call(b) && identical(as.character(b[[1]]), "{")) as.list(b)[-1L]
           else list(b)
  last <- stmts[[length(stmts)]]
  if (!is.call(last) || !identical(as.character(last[[1]]), "c")) {
    stop("et_conjugate(): the count body's last expression must be the counts, ",
         "written as c(...) -- c(successes, failures) for Beta, or ",
         "c(count_1, ..., count_k) for Dirichlet. Got: ", deparse(last)[1],
         call. = FALSE)
  }
  ctx <- new_ctx()
  head_src <- if (length(stmts) > 1L) {
    paste(vapply(stmts[-length(stmts)], et_transpile_expr, character(1), ctx = ctx),
          collapse = "\n")
  } else ""
  parts <- vapply(as.list(last)[-1L], et_transpile_expr, character(1), ctx = ctx)
  tail_src <- if (family == "beta") {
    if (length(parts) != 2L) {
      stop("et_conjugate(): a Beta count must be c(successes, failures); got ",
           length(parts), " values.", call. = FALSE)
    }
    sprintf("(%s, %s)", parts[1], parts[2])
  } else {
    sprintf("Int[%s]", paste(parts, collapse = ", "))
  }
  if (nzchar(head_src)) paste0(head_src, "\n", tail_src) else tail_src
}

# ---- runner -----------------------------------------------------------------

run_src_for <- function(model, has_obs) {
  init_parts <- vapply(model$par_names, function(nm) {
    p <- model$parameters[[nm]]
    val <- if (!is.null(p$dim)) {
      julia_array(matrix(p$init, p$dim[1], p$dim[2]), "Float64")
    } else if (p$n > 1L) {
      julia_vector(p$init, "Float64")
    } else {
      julia_float(p$init)
    }
    sprintf("%s=%s", nm, val)
  }, character(1))

  sprintf(paste0(
    "const INIT_PARS = (; %s)\n\n",
    "# One Gibbs sweep updates the trajectory, so the aggregates must already\n",
    "# agree with the STARTING trajectory before the first likelihood call. That\n",
    "# invariant is established once here and preserved by the sampler.\n",
    "function et_prepare!(X0)\n",
    "    reset_aggregates!(DATA)\n",
    "    apply_derived_summaries!(et_params(INIT_PARS), DATA, X0)\n",
    "    X0\n",
    "end\n\n",
    "function et_run(; n_sweeps, n_burn=0, n_adapts=0, seed=1, x_init=nothing,\n",
    "                  adtype=ADTypes.AutoForwardDiff(), save_x=nothing,\n",
    "                  save_every::Int=100)\n",
    "    X0 = x_init === nothing ? fill(1, %s, %s) : Matrix{Int}(x_init)\n",
    "    et_prepare!(X0)\n",
    "    m = et_the_model(DATA, LOGLIK%s)\n",
    "    init = (; X=X0, INIT_PARS...)\n",
    "    # The trajectory is n_timepoints x n_individuals Ints EVERY sweep. It is\n",
    "    # always kept live for conditioning; `save_states` only decides where the\n",
    "    # per-sweep COPY goes -- `:buffer` drops it, a (path, every) pair streams\n",
    "    # it to disc for post-hoc residuals.\n",
    "    save = save_x === nothing ? (X = :buffer,) : (X = (save_x, save_every),)\n",
    "    t0 = time()\n",
    "    chn = AbstractMCMC.sample(StableRNG(seed), m, SPL, n_sweeps;\n",
    "                              init=init, n_adapts=n_adapts, adtype=adtype,\n",
    "                              save_states=save, discard_initial=n_burn)\n",
    "    out = _et_extract(chn, %s)\n",
    "    out[\"_elapsed\"] = time() - t0\n",
    "    out\n",
    "end\n\n",
    "# Residuals from the archived trajectories. `specs` are ET residual objects\n",
    "# built R-side; `names` keys the result.\n",
    "#\n",
    "# `sync` rebuilds the tracked arrays for each draw before scoring it. A\n",
    "# spec-derived hazard reading `data.aggregates` needs them to agree with THIS\n",
    "# draw's trajectory, not with whatever was current when the fit ended --\n",
    "# `trajectory_summaries` deliberately does not do it unconditionally, since\n",
    "# for rates that ignore the aggregates it is a wasted pass per draw.\n",
    "# `par_draws` is one NamedTuple per archived sweep, in the same order: a\n",
    "# residual has to be scored under the parameters that PRODUCED its\n",
    "# trajectory, not under a single fixed set.\n",
    "# `files` are the per-flush chunks in sweep order, discovered R-side (which\n",
    "# also knows the format). They are read ONE AT A TIME: `read_states` would\n",
    "# stitch the whole run into memory, which on the badger model is ~15 GB and\n",
    "# segfaults -- its own docstring says to prefer one flush at a time.\n",
    "#\n",
    "# The result is a LAZY iterator, so only the current chunk is ever resident.\n",
    "# `trajectory_summaries` consumes an iterator directly, but pre-sizes its\n",
    "# output from `length`, which is why this declares HasLength.\n",
    "struct _ChunkedDraws\n",
    "    files::Vector{String}\n",
    "    par_draws::Vector{Any}\n",
    "    n::Int\n",
    "end\n",
    "Base.length(c::_ChunkedDraws) = c.n\n",
    "Base.IteratorSize(::Type{_ChunkedDraws}) = Base.HasLength()\n",
    "Base.IteratorEltype(::Type{_ChunkedDraws}) = Base.EltypeUnknown()\n",
    "# State is (next file, current chunk, index within it, draws emitted). The\n",
    "# `while` skips any empty chunk rather than assuming one file means one draw.\n",
    "function Base.iterate(c::_ChunkedDraws, st=(1, Vector{Matrix{Int}}(), 0, 0))\n",
    "    fi, buf, bi, done = st\n",
    "    while bi >= length(buf)\n",
    "        fi > length(c.files) && return nothing\n",
    "        buf = _et_load_chunk(c.files[fi])\n",
    "        fi += 1; bi = 0\n",
    "    end\n",
    "    bi += 1; done += 1\n",
    "    ((et_params(c.par_draws[done]), buf[bi]), (fi, buf, bi, done))\n",
    "end\n\n",
    "# One flush file -> its trajectories as Matrix{Int}. `Int8` is a STORAGE\n",
    "# format (see archive_draw); widening stops here, per draw, not for the run.\n",
    "function _et_load_chunk(file)\n",
    "    # JLD2 is an OPTIONAL backend (a PracticalBayes weakdep), so it is\n",
    "    # resolved at CALL time from Main -- where the wrapper loads it before\n",
    "    # sampling -- not with a `using` this module could not satisfy when no\n",
    "    # archive is in play.\n",
    "    jld2 = getfield(Main, :JLD2)\n",
    "    states = jld2.jldopen(file, \"r\") do f\n",
    "        f[\"states\"]\n",
    "    end\n",
    "    Matrix{Int}[Matrix{Int}(x) for x in states]\n",
    "end\n\n",
    "# Scoring one chunk R has already loaded (the `.rds` path): same specs, same\n",
    "# parameters, one chunk's worth of draws. Returns the per-draw residual\n",
    "# matrices, which R concatenates -- residuals are m x draws Float64, i.e.\n",
    "# tiny beside the trajectories, so THOSE are safe to accumulate.\n",
    "function et_residuals_chunk(Xs; specs, names, par_draws, sync=true, seed=1)\n",
    "    draws = [(et_params(par_draws[d]), Matrix{Int}(Xs[d])) for d in eachindex(Xs)]\n",
    "    it = sync ? aggregate_synced_draws(DATA, draws) : draws\n",
    "    R = trajectory_summaries(specs, DATA, it; rng=StableRNG(seed))\n",
    "    out = Dict{String,Any}()\n",
    "    for nm in names\n",
    "        out[String(nm)] = R[nm]\n",
    "    end\n",
    "    out\n",
    "end\n\n",
    "function et_residuals(files; specs, names, par_draws, sync=true, seed=1)\n",
    "    # A ONE-element character vector arrives from R as a bare `String`, not a\n",
    "    # `Vector{String}` -- and `collect(String, \"a/b.jld2\")` then iterates its\n",
    "    # CHARACTERS. Wrap the scalar case before collecting.\n",
    "    fs = files isa AbstractString ? [String(files)] : collect(String, files)\n",
    "    isempty(fs) && error(\"et_residuals: no trajectories were archived\")\n",
    "    draws = _ChunkedDraws(fs, collect(Any, par_draws), length(par_draws))\n",
    "    it = sync ? aggregate_synced_draws(DATA, draws) : draws\n",
    "    R = trajectory_summaries(specs, DATA, it; rng=StableRNG(seed))\n",
    "    out = Dict{String,Any}()\n",
    "    for nm in names\n",
    "        out[String(nm)] = R[nm]\n",
    "    end\n",
    "    out\n",
    "end\n\n",
    "# The log density at the initial parameters, for a cheap sanity check that\n",
    "# does not require running the sampler.\n",
    "function et_loglik_at(X0=fill(1, %s, %s))\n",
    "    et_prepare!(X0)\n",
    "    p = et_params(INIT_PARS)\n",
    "    (epidemic=LOGLIK(p, DATA, X0)%s)\n",
    "end\n\n",
    "# One latent sweep in isolation, at the initial parameters.\n",
    "function et_iffbs_sweep(X0; seed=1)\n",
    "    et_prepare!(X0)\n",
    "    iffbs!(et_params(INIT_PARS), DATA, X0, StableRNG(seed))\n",
    "    X0\n",
    "end"),
    paste(init_parts, collapse = ", "),
    julia_int(model$data$n_timepoints), julia_int(model$data$n_individuals),
    if (has_obs) ", OBSLOGLIK" else "",
    julia_symbol_vector(model$par_names),
    julia_int(model$data$n_timepoints), julia_int(model$data$n_individuals),
    if (has_obs) ", observation=OBSLOGLIK(p, DATA, X0)" else "")
}

# The `depends=` validator PERF_REPORT_2026-08-25.md section 6.1 proposed: move
# one parameter at a time and see which likelihood terms move with it. A term
# that moves when a parameter not in its declared set is perturbed has an
# under-declared `depends=`, which would silently zero that parameter's gradient
# contribution.
#
# The perturbation shrinks toward zero rather than growing, so a probability
# parameter stays inside (0, 1) and the check does not fail for the wrong reason.
check_src_for <- function(model, has_obs) {
  sprintf(paste0(
    "function et_perturb(nm::Symbol, factor)\n",
    "    v = getproperty(INIT_PARS, nm)\n",
    "    merge(INIT_PARS, NamedTuple{(nm,)}((v .* factor,)))\n",
    "end\n\n",
    "function et_check_depends(; factor=0.7, x_init=nothing)\n",
    "    X0 = x_init === nothing ? fill(1, %s, %s) : Matrix{Int}(x_init)\n",
    "    et_prepare!(X0)\n",
    "    p0 = et_params(INIT_PARS)\n",
    "    base_epi = LOGLIK(p0, DATA, X0)\n",
    "    base_obs = %s\n",
    "    out = Dict{String,Any}()\n",
    "    for nm in %s\n",
    "        p = et_params(et_perturb(nm, factor))\n",
    "        d_epi = LOGLIK(p, DATA, X0) - base_epi\n",
    "        d_obs = %s - base_obs\n",
    "        out[String(nm)] = Float64[d_epi, d_obs]\n",
    "    end\n",
    "    out\n",
    "end"),
    julia_int(model$data$n_timepoints), julia_int(model$data$n_individuals),
    if (has_obs) "OBSLOGLIK(p0, DATA, X0)" else "0.0",
    julia_symbol_vector(model$par_names),
    if (has_obs) "OBSLOGLIK(p, DATA, X0)" else "0.0")
}

# Chain -> plain arrays. A scalar parameter becomes a length-`n_draws` vector; a
# vector or matrix parameter becomes an `n_draws x prod(dim)` matrix, flattened
# column-major so R can restore the shape.
# Rebuild the sampled-parameter NamedTuples from the flat draws matrix R sends.
#
# The column layout is the model's declared parameter order and shapes, written
# out here rather than inferred, so the reader and `flatten_draws` cannot drift:
# a mismatch would silently score every trajectory under the wrong parameters.
unflatten_src <- function(model) {
  offset <- 0L
  fields <- character()
  for (nm in model$par_names) {
    p <- model$parameters[[nm]]
    n <- if (!is.null(p$dim)) prod(p$dim) else p$n
    rng <- sprintf("%s:%s", julia_int(offset + 1L), julia_int(offset + n))
    fields <- c(fields, if (!is.null(p$dim)) {
      sprintf("%s = reshape(view(M, d, %s), %s)", nm, rng,
              paste(julia_int(p$dim), collapse = ", "))
    } else if (p$n > 1L) {
      sprintf("%s = collect(view(M, d, %s))", nm, rng)
    } else {
      sprintf("%s = M[d, %s]", nm, julia_int(offset + 1L))
    })
    offset <- offset + n
  }
  sprintf(paste0(
    "# Columns: %s (total %s).\n",
    "function _et_unflatten(raw)\n",
    "    M = Matrix{Float64}(raw)\n",
    "    size(M, 2) == %s || error(\"_et_unflatten: expected %s columns, got \" *\n",
    "                              string(size(M, 2)))\n",
    "    return [(; %s) for d in axes(M, 1)]\n",
    "end"),
    paste(model$par_names, collapse = ", "), julia_int(offset),
    julia_int(offset), julia_int(offset), paste(fields, collapse = ", "))
}

extract_helper_src <- function() {
  paste0("# ---- extraction ", strrep("-", 60), "\n",
"function _et_flatten(raw)
    v = collect(raw)
    isempty(v) && return Float64[]
    if first(v) isa Real
        return Float64[Float64(x) for x in v]
    end
    d = length(first(v))
    out = Matrix{Float64}(undef, length(v), d)
    for (k, x) in enumerate(v)
        out[k, :] = Float64.(vec(collect(x)))
    end
    out
end

function _et_extract(chn, names)
    out = Dict{String,Any}()
    for nm in names
        out[String(nm)] = _et_flatten(vec(chn[nm]))
    end
    out
end")
}

# Replace symbols in an R expression according to a named character map.
subst_symbols <- function(e, map) {
  if (is.symbol(e)) {
    nm <- as.character(e)
    if (nm %in% names(map)) return(as.symbol(map[[nm]]))
    return(e)
  }
  if (!is.call(e)) return(e)
  for (k in seq_along(e)[-1]) e[[k]] <- subst_symbols(e[[k]], map)
  e
}


# The population-level entry points: residuals collected DURING a fit (never
# storing a trajectory), and R_i across archived draws. Split into their own
# generator because R's sprintf refuses a format string over 8192 chars, and
# folding these into run_src_for() crossed it.
population_src_for <- function(model, has_obs) {
  sprintf(paste0(
    "# Residuals scored DURING the fit, never storing a trajectory. This is the\n",
    "# one place the sweep loop cannot live inside `sample`: the collector has to\n",
    "# see each draw as it is produced, so the loop is written out and the chain\n",
    "# is rebuilt from what it kept.\n",
    "#\n",
    "# `save_states = (X = :buffer,)` keeps X live for conditioning while dropping\n",
    "# the per-sweep copy -- the whole point being that nothing about X is stored.\n",
    "function et_collect_run(; specs, names, n_sweeps, n_burn=0, n_adapts=0,\n",
    "                          seed=1, x_init=nothing,\n",
    "                          adtype=ADTypes.AutoForwardDiff(),\n",
    "                          n_keep=500, thin=1)\n",
    "    X0 = x_init === nothing ? fill(1, %s, %s) : Matrix{Int}(x_init)\n",
    "    et_prepare!(X0)\n",
    "    m = et_the_model(DATA, LOGLIK%s)\n",
    "    init = (; X=X0, INIT_PARS...)\n",
    "    coll = SummaryCollector(specs, DATA, n_keep; rng=StableRNG(seed), thin=thin)\n",
    "    t0 = time()\n",
    "    # The sweep loop is stepped BY HAND rather than run through `sample`.\n",
    "    # PracticalBayes has its own loop and forwards `kwargs` only to `step`; it\n",
    "    # never invokes AbstractMCMC's `callback`. A callback passed to `sample`\n",
    "    # would therefore be ACCEPTED AND NEVER CALLED -- the collector would come\n",
    "    # back empty with no error at all. Stepping manually is what actually sees\n",
    "    # every draw.\n",
    "    rng = StableRNG(seed)\n",
    "    t, state = AbstractMCMC.step(rng, m, SPL; init=init, adtype=adtype,\n",
    "                                 n_adapts=n_adapts)\n",
    "    for _ in 1:n_burn\n",
    "        t, state = AbstractMCMC.step(rng, m, SPL, state; n_adapts=n_adapts)\n",
    "    end\n",
    "    kept = Vector{Any}()\n",
    "    collect_summaries!(coll, et_params(t), t.X)\n",
    "    push!(kept, t)\n",
    "    for _ in 2:n_sweeps\n",
    "        t, state = AbstractMCMC.step(rng, m, SPL, state; n_adapts=n_adapts)\n",
    "        collect_summaries!(coll, et_params(t), t.X)\n",
    "        push!(kept, t)\n",
    "    end\n",
    "    out = Dict{String,Any}()\n",
    "    for nm in %s\n",
    "        out[String(nm)] = _et_flatten([getproperty(k, nm) for k in kept])\n",
    "    end\n",
    "    out[\"_elapsed\"] = time() - t0\n",
    "    R = finish(coll)\n",
    "    res = Dict{String,Any}()\n",
    "    for nm in names\n",
    "        res[String(nm)] = R[nm]\n",
    "    end\n",
    "    out[\"_residuals\"] = res\n",
    "    out\n",
    "end\n\n",
    "# R_i across archived draws. A population-level quantity, so it goes through\n",
    "# `summarize_population` rather than the per-individual residual driver -- see\n",
    "# case_reproduction_numbers for why it cannot use that contract.\n",
    "#\n",
    "# The trajectories are read one flush file at a time, exactly as et_residuals\n",
    "# does, so a large archive is never materialised.\n",
    "function et_reproduction_numbers(files; components, secondary, infection,\n",
    "                                 infectious_state, weight, par_draws,\n",
    "                                 draws=nothing)\n",
    "    fs = files isa AbstractString ? [String(files)] : collect(String, files)\n",
    "    isempty(fs) && error(\"et_reproduction_numbers: no trajectories were archived\")\n",
    "    wt = weight === nothing ? (mo, d, x, i, t) -> 1.0 : weight\n",
    "    keep = draws === nothing ? nothing : Set(Int.(draws))\n",
    "    cols = Vector{Vector{Union{Float64,Missing}}}()\n",
    "    d = 0\n",
    "    for f in fs\n",
    "        for X in _et_load_chunk(f)\n",
    "            d += 1\n",
    "            (keep === nothing || d in keep) || continue\n",
    "            mo = et_params(par_draws[d])\n",
    "            reset_aggregates!(DATA)\n",
    "            apply_derived_summaries!(mo, DATA, X)\n",
    "            push!(cols, case_reproduction_numbers(mo, DATA, X;\n",
    "                      components=components, secondary=secondary,\n",
    "                      infection=infection, infectious_state=infectious_state,\n",
    "                      weight=wt))\n",
    "        end\n",
    "    end\n",
    "    isempty(cols) && error(\"et_reproduction_numbers: no draws selected\")\n",
    "    Dict{String,Any}(\"R\" => reduce(hcat, cols))\n",
    "end\n\n",
    ""),
    julia_int(model$data$n_timepoints), julia_int(model$data$n_individuals),
    if (has_obs) ", OBSLOGLIK" else "",
    julia_symbol_vector(model$par_names))
}
