# The Julia module for an et_latent_data() model: ET's latent_hmm_data() with a
# PracticalBayes model whose one likelihood term is the marginal over every
# individual's weighted chain, given the latent block, and a Gibbs sampler that
# updates the latent block once per sweep before the HMC blocks.
#
# Written against ET's settled latent API (latent_hmm_data, @latents and the
# latent sampler) and PracticalEpiBayes' LatentBlock / latent_loglik /
# latent_kernel, rather than the @latent_model and latent_gibbs sugar over
# them, so a rename there does not reach the generated text.

# The latent block's name in the PracticalBayes model.
LATENT_BLOCK <- "latents"

generate_latent_module <- function(model, blocks, module_name) {
  d <- model$data
  states <- d$transitions$states
  K <- length(states)
  prefix <- paste0(module_name, "_")
  payload <- list()
  add_payload <- function(name, value, el = NULL, vector = FALSE) {
    payload[[paste0(prefix, name)]] <<- value
    julia_typed_binding(name, paste0("Main.", prefix, name), value, el, vector)
  }

  helpers <- d$helpers
  fn_src <- character()
  # With per_evaluation, a body may read a whole number or a logical it added
  # to `model`, so its number type is floored at Float64.
  pe <- !is.null(d$per_evaluation)
  role <- function(f, r, nm, extra = character(), n = NULL) {
    et_transpile_role(f, r, nm, helpers, extra, n_states = n, float_floor = pe)
  }

  # --- once per evaluation ----------------------------------------------------
  if (pe) {
    pe_names <- names(per_eval_field_reads(d$per_evaluation, helpers))
    clash <- intersect(pe_names, model$all_names)
    if (length(clash)) {
      stop("et_latent_data(): per_evaluation returns ", paste(clash, collapse = ", "),
           ", which would replace the parameter(s) of that name.", call. = FALSE)
    }
    fn_src <- c(fn_src, SETTLE_SRC,
                role(d$per_evaluation, "per_evaluation", "et_per_evaluation")$src)
  }

  # --- weights, starting state, observations --------------------------------
  shared <- if (is.null(d$transitions$shared)) shared_role(NULL, "et_shared", helpers)
            else {
              r <- role(d$transitions$shared, "shared", "et_shared")
              list(name = "et_shared", src = r$src, reads = r$reads, opaque = r$opaque)
            }
  fn_src <- c(fn_src, shared$src)
  rate_role <- if (is.null(shared$name)) "rate" else "rate_shared"
  rate_names <- character()
  for (tr in d$transitions$transitions) {
    nm <- sprintf("et_weight_%s_%s", tr$from, tr$to)
    fn_src <- c(fn_src, role(tr$rate, rate_role, nm, shared$reads)$src)
    rate_names[paste0(tr$from, "->", tr$to)] <- nm
  }
  fn_src <- c(fn_src, role(d$starting_state, "latent_start", "et_starting_state", n = K)$src)
  obs_shared <- list(reads = character(), opaque = FALSE)
  if (!is.null(d$observation_shared)) {
    obs_shared <- role(d$observation_shared, "latent_obs_shared", "et_observation_shared")
    fn_src <- c(fn_src, obs_shared$src)
  }
  obs <- list(reads = character(), opaque = FALSE)
  if (!is.null(d$observation_weights)) {
    # The observations read `shared`, so their number type must cover whatever
    # observation_shared reads, as a rate's covers its shared step's: otherwise a
    # shared value carrying derivatives is converted to a plain number.
    obs <- role(d$observation_weights, "latent_obs", "et_observation_weights",
                extra = obs_shared$reads, n = K)
    obs$opaque <- obs$opaque || obs_shared$opaque
    fn_src <- c(fn_src, obs$src)
  }
  if (!is.null(d$shared_key)) {
    fn_src <- c(fn_src, transpile_plain_fn(d$shared_key, "et_shared_key", 1L, helpers,
                                           "et_latent_data(): `shared_key`"))
  }
  cache <- cache_src(model, obs_shared, obs)

  # --- extras ---------------------------------------------------------------
  bindings <- character()
  for (nm in names(d$extras)) {
    x <- d$extras[[nm]]
    bindings <- c(bindings, if (inherits(x, "et_spatial_captures"))
      captures_binding(nm, x, add_payload) else add_payload(nm, x))
  }
  if (!is.null(d$sampling_period)) {
    bindings <- c(bindings, add_payload("_etr_sampling_period", d$sampling_period),
                  sprintf(paste0("const _etr_sampling_period_pairs = [(_etr_sampling_period[i, 1], ",
                                 "_etr_sampling_period[i, 2]) for i in 1:%s]"),
                          julia_int(d$n_individuals)))
  }

  # --- latents --------------------------------------------------------------
  decl <- character(); init <- character()
  for (nm in names(d$latents)) {
    r <- latent_src(nm, d$latents[[nm]], d, helpers, add_payload)
    fn_src <- c(fn_src, r$fns)
    bindings <- c(bindings, r$bindings)
    decl <- c(decl, r$decl)
    if (!is.null(r$init)) init <- c(init, sprintf("%s = %s", nm, r$init))
  }
  latents_src <- sprintf("const LATENTS = @latents begin\n%s\nend",
                         indent(paste(decl, collapse = "\n")))
  init_src <- sprintf("const LATENT_INIT = (; %s)", paste(init, collapse = ",\n    "))

  states_src <- sprintf("const STATES = %s", julia_symbol_vector(states))
  trans_src <- transitions_src("TRANS", d$transitions, rate_names, NULL, shared$name)

  kw <- c(sprintf("n_individuals=%s", julia_int(d$n_individuals)),
          sprintf("n_timepoints=%s", julia_int(d$n_timepoints)),
          "trans_mat=TRANS", "starting_state=et_starting_state",
          if (!is.null(d$observation_weights)) "observation_weights=et_observation_weights",
          if (!is.null(d$observation_shared)) "observation_shared=et_observation_shared",
          if (!is.null(d$shared_key)) c("shared_key=et_shared_key",
                                        sprintf("n_keys=%s", julia_int(d$n_keys))),
          if (pe) "per_evaluation=et_per_evaluation",
          if (!is.null(cache)) sprintf("cache=%s", cache$src),
          "latents=LATENTS", "init=LATENT_INIT",
          if (!is.null(d$sampling_period)) "sampling_period=_etr_sampling_period_pairs",
          vapply(names(d$extras), function(nm) sprintf("%s=%s", nm, nm), character(1)))
  data_src <- sprintf("const DATA = latent_hmm_data(;\n%s,\n)",
                      indent(paste(kw, collapse = ",\n")))

  header <- sprintf(paste0(
    "# Generated by EpidemicTrajectoriesR. DO NOT EDIT BY HAND.\n",
    "#\n",
    "# This file is self-contained: it runs as an ordinary Julia script once the\n",
    "# payload arrays named `%s*` exist in `Main`.\n"), prefix)

  src <- paste0(
    header, "\n",
    sprintf("module %s\n", module_name),
    join_blocks(
      julia_usings(),
      section("payload from R", qualify_base(bindings)),
      section("state space", states_src),
      section("user functions", c(et_transpile_helpers(helpers), fn_src,
                                  log_density_src(model))),
      section("weighted transitions", trans_src),
      section("latent variables", c(latents_src, init_src)),
      section("data", data_src),
      section("the model", pb_latent_model_src(model)),
      section("sampler", latent_gibbs_src(model, blocks, helpers)),
      section("runner", qualify_base(run_src_latent(model))),
      section("depends= validator", qualify_base(check_src_marginal(model))),
      section("parameter draws from R", qualify_base(unflatten_src(model))),
      qualify_base(extract_helper_src())
    ),
    sprintf("\n\nend # module %s\n", module_name))

  list(src = src, payload = payload, module = module_name, blocks = blocks,
       depends = c(list(latent = NULL), log_density_depends(model)), model = model,
       cache = if (!is.null(cache)) cache[c("shared", "weights")])
}

# ---- what each per_evaluation value reads -----------------------------------
#
# The cache is keyed on the parameters a body reads, and a body may read a value
# `per_evaluation` added to `model`. Charging such a read with everything
# `per_evaluation` reads would make the capture term depend on every movement
# parameter, and the cache would never be used, so each returned field is
# traced back to the parameters it is computed from: through locals, element
# assignments, loops, the conditions it is computed under, and helpers.
# Over-counting is the safe direction throughout (a stale cache is silently
# wrong), so a local's reads only ever grow, and the pass repeats until they
# stop growing, which carries reads around a loop.
#
# Returns a named list, field -> parameter names, or NULL when the body cannot
# be traced (an et_julia(), or a return that is not a list()).
per_eval_field_reads <- function(f, helpers) {
  if (is.null(f) || inherits(f, "et_julia")) return(NULL)
  helper_names <- vapply(helpers, function(h) h$name, character(1))
  deps <- new.env(parent = emptyenv())
  fields <- list()
  traceable <- TRUE
  changed <- FALSE

  add <- function(nm, r) {
    old <- if (exists(nm, deps, inherits = FALSE)) get(nm, deps) else character()
    new <- union(old, r)
    if (length(new) > length(old)) changed <<- TRUE
    assign(nm, new, deps)
  }
  reads_of <- function(e) {
    if (is.symbol(e)) {
      nm <- as.character(e)
      return(if (exists(nm, deps, inherits = FALSE)) get(nm, deps) else character())
    }
    if (!is.call(e)) return(character())
    if (identical(e[[1]], as.symbol("$")) && identical(e[[2]], quote(model))) {
      return(as.character(e[[3]]))
    }
    out <- character()
    if (is.symbol(e[[1]]) && as.character(e[[1]]) %in% helper_names) {
      out <- helper_reads(as.character(e[[1]]), helpers)
    }
    for (a in as.list(e)[-1L]) out <- union(out, reads_of(a))
    out
  }
  record <- function(lst, ctrl) {
    nms <- names(lst)[-1L]
    for (k in seq_along(nms)) {
      fields[[nms[k]]] <<- union(fields[[nms[k]]], union(reads_of(lst[[k + 1L]]), ctrl))
    }
  }
  is_list_call <- function(e) is.call(e) && identical(e[[1]], as.symbol("list"))
  walk <- function(e, ctrl, last) {
    if (!is.call(e)) {
      if (last) traceable <<- FALSE
      return(invisible())
    }
    head <- as.character(e[[1]])
    if (head == "{") {
      n <- length(e)
      for (k in seq_len(n)[-1L]) walk(e[[k]], ctrl, last && k == n)
    } else if (head %in% c("<-", "=")) {
      lhs <- e[[2]]
      r <- union(reads_of(e[[3]]), ctrl)
      if (is.call(lhs)) {           # x[k] <- v: the index matters too
        r <- union(r, reads_of(lhs))
        lhs <- lhs[[2]]
      }
      add(as.character(lhs), r)
    } else if (head == "if") {
      c2 <- union(ctrl, reads_of(e[[2]]))
      walk(e[[3]], c2, last)
      if (length(e) == 4L) walk(e[[4]], c2, last) else if (last) traceable <<- FALSE
    } else if (head == "for") {
      r <- union(reads_of(e[[3]]), ctrl)
      add(as.character(e[[2]]), r)
      walk(e[[4]], r, FALSE)
    } else if (head == "return") {
      if (length(e) == 2L && is_list_call(e[[2]])) record(e[[2]], ctrl) else traceable <<- FALSE
    } else if (head == "list") {
      if (last) record(e, ctrl)
    } else if (last) {
      traceable <<- FALSE
    }
    invisible()
  }
  repeat {
    changed <- FALSE
    walk(body(f), character(), TRUE)
    if (!changed) break
  }
  if (!traceable || !length(fields)) return(NULL)
  fields
}

# The sampled parameters behind a body's reads, a value per_evaluation adds
# resolved to what that value is computed from.
cache_reads <- function(reads, model, pe_fields, what) {
  out <- expand_to_sampled(intersect(reads, model$all_names), model)
  for (nm in setdiff(reads, model$all_names)) {
    if (is.null(pe_fields[[nm]])) {
      stop("et_latent_data(): ", what, " reads model$", nm, ", which is not a ",
           "parameter, a derived value, or a value per_evaluation returns.",
           call. = FALSE)
    }
    out <- union(out, expand_to_sampled(pe_fields[[nm]], model))
  }
  out
}

# The `cache = (; ...)` keyword, or NULL for none.
cache_src <- function(model, obs_shared, obs) {
  d <- model$data
  cc <- d$cache
  if (is.null(cc)) return(NULL)
  pe_fields <- per_eval_field_reads(d$per_evaluation, d$helpers)
  traced <- !obs_shared$opaque && !obs$opaque &&
    (is.null(d$per_evaluation) || !is.null(pe_fields) ||
       !length(setdiff(union(obs_shared$reads, obs$reads), model$all_names)))
  if (traced) {
    shared <- cache_reads(obs_shared$reads, model, pe_fields, "observation_shared")
    weights <- setdiff(cache_reads(obs$reads, model, pe_fields, "observation_weights"),
                       shared)
  }
  if (isTRUE(cc$derive)) {
    if (!traced) {
      stop("et_latent_data(): cache = TRUE cannot work out what the observations ",
           "read, because ", if (obs_shared$opaque || obs$opaque)
             "a body is an et_julia() without `reads`"
           else "per_evaluation does not end in a list() of its values",
           ". Give `cache = list(shared = , weights = )` by hand.", call. = FALSE)
    }
  } else {
    given <- list(shared = cc$shared %||% character(), weights = cc$weights %||% character())
    unknown <- setdiff(unlist(given), model$all_names)
    if (length(unknown)) {
      stop("et_latent_data(): `cache` names unknown parameter(s): ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    if (traced) {
      miss_s <- setdiff(shared, given$shared)
      miss_w <- setdiff(weights, union(given$shared, given$weights))
      if (length(miss_s) || length(miss_w)) {
        stop("et_latent_data(): `cache` leaves out parameters the observations ",
             "read, so their values would go stale: ",
             if (length(miss_s)) paste0("shared needs ", paste(miss_s, collapse = ", "), "; "),
             if (length(miss_w)) paste0("weights needs ", paste(miss_w, collapse = ", "), "; "),
             "or use cache = TRUE.", call. = FALSE)
      }
    }
    shared <- expand_to_sampled(given$shared, model)
    weights <- setdiff(expand_to_sampled(given$weights, model), shared)
  }
  parts <- c(sprintf("shared = %s", julia_symbol_tuple(shared)),
             sprintf("weights = %s", julia_symbol_tuple(weights)),
             if (!is.null(cc$latents)) sprintf("latents = %s", julia_symbol_tuple(cc$latents)))
  list(src = sprintf("(; %s)", paste(parts, collapse = ", ")),
       shared = shared, weights = weights)
}

# Values per_evaluation adds to `model`: a number in the parameter type, a whole
# number or a logical as it is (it may index or branch), arrays and tuples
# element by element.
SETTLE_SRC <- paste(
  "_etr_settle(::Type{T}, x::Union{Integer,Bool}) where {T} = x",
  "_etr_settle(::Type{T}, x::Real) where {T} = convert(T, x)",
  "_etr_settle(::Type{T}, x::Tuple) where {T} = map(v -> _etr_settle(T, v), x)",
  "_etr_settle(::Type{T}, x::AbstractArray{<:Union{Integer,Bool}}) where {T} = x",
  "_etr_settle(::Type{T}, x::AbstractArray{<:Real}) where {T} = convert(Array{T}, x)",
  "_etr_settle(::Type{T}, x) where {T} = x", sep = "\n")

# A spatial-captures extra: its arrays as payload, then ET's constructor.
captures_binding <- function(nm, x, add_payload) {
  p <- paste0("_etr_", nm, "_")
  latent <- if (length(x$latent_detectors)) {
    add_payload(paste0(p, "latent"), x$latent_detectors, "Int", vector = TRUE)
  } else sprintf("const %slatent = Int[]", p)
  c(add_payload(paste0(p, "xy"), x$detector_xy, "Float64"),
    add_payload(paste0(p, "active"), x$active, "Bool"),
    vapply(c("i", "occasion", "detector", "n"), function(f)
      add_payload(paste0(p, f), x$captures[[f]], "Int", vector = TRUE), character(1)),
    latent,
    sprintf(paste0(
      "const %s = spatial_captures(; detector_xy=%sxy, active=%sactive,\n",
      "    n_individuals=%s, latent_detectors=%slatent,\n",
      "    captures=zip(%si, %soccasion, %sdetector, %sn))"),
      nm, p, p, julia_int(x$n_individuals), p, p, p, p, p))
}

coord_prior_src <- function(prior) {
  if (inherits(prior, "custom_dist")) return(dist_to_julia(prior))
  paste0("(", paste(vapply(prior, dist_to_julia, character(1)), collapse = ", "), ")")
}

# One latent: its declaration line, the functions it names, its payload, and
# the Julia expression for its starting value (or NULL to draw from the prior).
latent_src <- function(nm, k, d, helpers, add_payload) {
  fns <- character(); bindings <- character(); init <- NULL
  what <- sprintf("latent '%s'", nm)
  if (k$kind == "field") {
    if (!is.null(k$index)) {
      fns <- transpile_plain_fn(k$index, paste0("et_index_", nm), 1L, helpers,
                                paste0(what, ": `index`"))
      where <- sprintf("index=et_index_%s", nm)
    } else {
      fns <- transpile_plain_fn(k$touches, paste0("et_touches_", nm), 1L, helpers,
                                paste0(what, ": `touches`"))
      where <- sprintf("touches=et_touches_%s", nm)
    }
    decl <- sprintf("%s ~ Field(%s, %s; prior=%s, %s, scale=%s%s)", nm, julia_int(k$dim),
                    julia_int(k$n), coord_prior_src(k$prior), where, julia_float(k$scale),
                    field_extras_src(nm, k, d, add_payload, function(b) bindings <<- c(bindings, b)))
  } else if (k$kind == "shared") {
    decl <- sprintf("%s ~ Shared(%s, %s; prior=%s, scale=%s)", nm, julia_int(k$dim),
                    julia_int(k$n), coord_prior_src(k$prior), julia_float(k$scale))
  } else {
    if (is.function(k$prior)) {
      r <- transpile_param_prior(k$prior, paste0("et_prior_", nm), helpers,
                                 paste0(what, ": `prior`"))
      fns <- c(fns, r$src)
      prior <- paste0("et_prior_", nm)
      discrete <- .et_dist_ctors[[r$dist]] %in% .discrete_families
    } else {
      prior <- dist_to_julia(k$prior)
      discrete <- is_discrete_dist(k$prior)
    }
    el <- if (discrete) "Int" else "Float64"
    opts <- character()
    if (!is.null(k$observed)) {
      known <- !is.na(k$observed)
      vals <- k$observed; vals[!known] <- 0
      if (discrete && any(vals != round(vals))) {
        stop(what, ": `observed` holds non-whole values for a discrete prior.",
             call. = FALSE)
      }
      p <- paste0("_etr_", nm, "_observed_")
      bindings <- c(bindings,
        add_payload(paste0(p, "values"), as.numeric(vals), el, vector = TRUE),
        add_payload(paste0(p, "known"), known, "Bool", vector = TRUE),
        sprintf("const _etr_%s_observed = Union{Missing,%s}[k ? v : missing for (v, k) in zip(%svalues, %sknown)]",
                nm, el, p, p))
      opts <- c(opts, sprintf("observed=_etr_%s_observed", nm))
    }
    if (!is.null(k$allowed)) {
      fns <- c(fns, transpile_plain_fn(k$allowed, paste0("et_allowed_", nm), 3L, helpers,
                                       paste0(what, ": `allowed`")))
      opts <- c(opts, sprintf("allowed=et_allowed_%s", nm))
    }
    opts <- c(opts, sprintf("scale=%s", julia_float(k$scale)))
    decl <- sprintf("%s ~ PerIndividual(%s; %s)", nm, prior, paste(opts, collapse = ", "))
    if (is.numeric(k$init)) {
      if (discrete && any(k$init != round(k$init))) {
        stop(what, ": `init` holds non-whole values for a discrete prior.", call. = FALSE)
      }
      bindings <- c(bindings, add_payload(paste0("_etr_", nm, "_init"), as.numeric(k$init),
                                          el, vector = TRUE))
      init <- sprintf("_etr_%s_init", nm)
    }
    return(list(decl = decl, fns = fns, bindings = bindings, init = init))
  }

  if (is.numeric(k$init)) {
    bindings <- c(bindings, add_payload(paste0("_etr_", nm, "_init"), k$init, "Float64"))
    init <- sprintf("_etr_%s_init", nm)
  } else if (inherits(k$init, "et_starting_centres")) {
    r <- starting_centres_src(k$init, d, helpers, paste0("et_init_index_", nm), what)
    fns <- c(fns, r$fns); init <- r$src
  } else if (inherits(k$init, "et_starting_detectors")) {
    s <- k$init
    r <- starting_centres_src(s$centres, d, helpers, paste0("et_init_centres_index_", nm), what)
    fns <- c(fns, r$fns,
             transpile_plain_fn(s$index, paste0("et_init_index_", nm), 1L, helpers,
                                paste0(what, ": the starting detectors' `index`")))
    init <- sprintf("starting_detectors(%s, %s, et_init_index_%s; bounds=%s)",
                    extras_name_of(s$captures, d$extras), r$src, nm, bounds_src(s$bounds))
  }
  list(decl = decl, fns = fns, bindings = bindings, init = init)
}

starting_centres_src <- function(s, d, helpers, fn_name, what) {
  list(fns = transpile_plain_fn(s$index, fn_name, 1L, helpers,
                                paste0(what, ": the starting centres' `index`")),
       src = sprintf("starting_centres(%s, %s, %s; bounds=%s)",
                     extras_name_of(s$captures, d$extras), fn_name, julia_int(s$n),
                     bounds_src(s$bounds)))
}

bounds_src <- function(b) {
  paste0("(", paste(vapply(b, julia_float, character(1)), collapse = ", "), ")")
}

# The priors, the deterministics, the latent block, and one likelihood term:
# the marginal over every weighted chain given the block. The block's own
# priors are scored by the sampler, and the parameter-dependent ones by the
# likelihood, so the block has no density here.
pb_latent_model_src <- function(model) {
  lines <- character()
  for (nm in model$par_names) {
    lines <- c(lines, sprintf("%s ~ %s", nm, par_dist_src(nm, model$parameters[[nm]])))
  }
  for (nm in names(model$derived)) {
    lines <- c(lines, sprintf("%s := %s", nm,
                              et_transpile_expr(model$derived[[nm]], new_ctx())))
  }
  pars <- c(model$par_names, names(model$derived))
  lines <- c(lines,
             sprintf("%s ~ LatentBlock(data)", LATENT_BLOCK),
             sprintf("pars = (; %s)", paste(sprintf("%s=%s", pars, pars), collapse = ", ")),
             sprintf("@addlogprob! latent_loglik(pars, data, %s)", LATENT_BLOCK),
             # Like the likelihood term, these carry no depends=: every block
             # evaluates them.
             log_density_lines(model, with_depends = FALSE))
  sprintf("@model function et_the_model(data)\n%s\nend",
          indent(paste(lines, collapse = "\n")))
}

# A field's joint moves and particle Gibbs, as Field keywords (with a leading
# comma), the anchors passed as payload.
field_extras_src <- function(nm, k, d, add_payload, add_binding) {
  out <- character()
  if (length(k$moves)) {
    out <- c(out, sprintf("moves=(%s,)", paste(julia_symbol(k$moves), collapse = ", ")))
  }
  p <- k$pgas
  if (!is.null(p)) {
    parts <- sprintf("particles=%s", julia_int(p$particles))
    if (!is.null(p$scales)) {
      parts <- c(parts, sprintf("scales=(%s,)", paste(vapply(p$scales, julia_float, character(1)),
                                                     collapse = ", ")))
    }
    if (!is.null(p$anchors)) {
      if (dim(p$anchors)[3] != d$n_individuals) {
        stop(sprintf("latent '%s': `pgas` anchors must be %d x %d x %d.", nm, k$dim, k$n,
                     d$n_individuals), call. = FALSE)
      }
      add_binding(add_payload(paste0("_etr_", nm, "_anchors"), p$anchors, "Float64"))
      parts <- c(parts, sprintf("anchors=_etr_%s_anchors", nm))
    }
    if (!is.null(p$anchor_scale)) {
      parts <- c(parts, sprintf("anchor_scale=%s", julia_float(p$anchor_scale)))
    }
    out <- c(out, sprintf("pgas=(; %s)", paste(parts, collapse = ", ")))
  }
  if (length(out)) paste0(", ", paste(out, collapse = ", ")) else ""
}

# The latent kernel first, so each sweep's HMC blocks see the latents it drew.
# With an et_interweave() block, the latent kernel also owns that block's
# parameters and moves them with the field's paths.
latent_gibbs_src <- function(model, blocks, helpers = list()) {
  is_iw <- vapply(blocks, function(b) b$kind == "interweave", logical(1))
  if (sum(is_iw) > 1L) {
    stop("more than one et_interweave() block was given.", call. = FALSE)
  }
  iw_src <- NULL
  if (any(is_iw)) {
    b <- blocks[[which(is_iw)]]
    iw_src <- interweave_src(model, b, helpers)
    first <- sprintf("%s => interweave_kernel(LATENT!, %s; n_adapt=N_ADAPT, params=et_params)",
                     julia_symbol_tuple(c(LATENT_BLOCK, b$vars)), iw_src$ctor)
  } else {
    first <- sprintf("%s => latent_kernel(LATENT!; params=et_params)", julia_symbol(LATENT_BLOCK))
  }
  entries <- c(first, vapply(blocks[!is_iw], function(b) block_entry(model, b), character(1)))
  paste0(
    params_closure_src(model), "\n\n",
    if (!is.null(iw_src)) paste0(iw_src$fns, "\n\n"),
    "# `LATENT!` is the latent sampler, built per run: its proposal scales adapt\n",
    "# over the run's warm-up of N_ADAPT sweeps.\n",
    sprintf("et_sampler_for(DATA, LATENT!, N_ADAPT=0) = Gibbs(\n%s,\n)",
            indent(paste(entries, collapse = ",\n"))))
}

# An et_interweave() block: the step (transpiled, or each parameter plus u),
# the elements to stretch, and the Interweave constructor.
interweave_src <- function(model, b, helpers) {
  field <- model$data$latents[[b$field]]
  if (is.null(field) || field$kind != "field") {
    stop("et_interweave(): `field` must name an et_field() latent; got '", b$field, "'.",
         call. = FALSE)
  }
  pars <- julia_symbol_tuple(b$vars)
  if (is.null(b$step)) {
    fns <- sprintf(paste0(
      "# Interweaving: each of %s moved by u, with no Jacobian.\n",
      "et_interweave_step(v, u) = ((; %s), 0.0)"), paste(b$vars, collapse = ", "),
      paste(sprintf("%s=v.%s .+ u", b$vars, b$vars), collapse = ", "))
  } else {
    r <- et_transpile_role(b$step, "interweave_step", "et_interweave_step_r", helpers)
    fns <- paste0(r$src, "\n\n", sprintf(paste0(
      "# Interweaving: the step written in R, given the parameters and u; a\n",
      "# logjac that is not finite rejects the move.\n",
      "function et_interweave_step(v, u)\n",
      "    r = et_interweave_step_r(et_params(v), u)\n",
      "    isfinite(r.logjac) || return nothing\n",
      "    (NamedTuple{%s}(map(p -> getproperty(r, p), %s)), Float64(r.logjac))\n",
      "end"), pars, pars))
  }
  el <- "nothing"
  if (!is.null(b$elements)) {
    fns <- paste0(fns, "\n\n", transpile_plain_fn(b$elements, "et_interweave_elements", 1L, helpers,
                                                  "et_interweave(): `elements`"))
    el <- "et_interweave_elements"
  }
  list(fns = fns, ctor = sprintf(
    "Interweave(%s, et_interweave_step; field=%s, elements=%s, scale=%s, centred=%s, n=%s)",
    pars, julia_symbol(b$field), el, julia_float(b$scale), julia_float(b$centred), julia_int(b$n)))
}

run_src_latent <- function(model) {
  init_parts <- vapply(model$par_names, function(nm) {
    p <- model$parameters[[nm]]
    val <- if (p$n > 1L) julia_vector(p$init, "Float64") else julia_float(p$init)
    sprintf("%s=%s", nm, val)
  }, character(1))
  sprintf(paste0(
    "const INIT_PARS = (; %s)\n\n",
    "et_prepare!(X0) = X0\n\n",
    "function et_run(; n_sweeps, n_burn=0, n_adapts=0, seed=1, x_init=nothing,\n",
    "                  adtype=ADTypes.AutoForwardDiff(), save_x=nothing,\n",
    "                  save_every::Int=100)\n",
    "    x_init === nothing && save_x === nothing ||\n",
    "        error(\"a weighted-chain model has no trajectory to start or save\")\n",
    "    m = et_the_model(DATA)\n",
    "    spl = et_sampler_for(DATA, epidemic_latent_sampler(DATA; n_adapt=n_adapts, seed=seed), n_adapts)\n",
    "    t0 = time()\n",
    "    # The latent block is kept live for conditioning but not stored per sweep.\n",
    "    chn = AbstractMCMC.sample(StableRNG(seed), m, spl, n_sweeps;\n",
    "                              init=INIT_PARS, n_adapts=n_adapts, adtype=adtype,\n",
    "                              save_states=(%s = :buffer,), discard_initial=n_burn)\n",
    "    out = _et_extract(chn, %s)\n",
    "    out[\"_elapsed\"] = time() - t0\n",
    "    out[\"_latents\"] = et_latent_values()\n",
    "    out\n",
    "end\n\n",
    "# The latents the data object holds: the starting values, or the last ones\n",
    "# the likelihood saw.\n",
    "et_latent_values() = Dict{String,Any}(String(k) => v for (k, v) in pairs(current_state(DATA).values))\n\n",
    "# The marginal log likelihood at a parameter set and the current latents, in\n",
    "# total (with the parameter-dependent latent priors) or per individual\n",
    "# (without them).\n",
    "et_marginal_at(v) = epidemic_marginal_loglik(DATA)(et_params(v), DATA, current_state(DATA))\n",
    "et_hmm_logliks_at(v) = Float64.(hmm_logliks(et_params(v), DATA, current_state(DATA)))\n",
    "et_loglik_at(X0=nothing) = (latent=et_marginal_at(INIT_PARS),)"),
    paste(init_parts, collapse = ", "), LATENT_BLOCK,
    julia_symbol_vector(model$par_names))
}

# ---- latents on a sampled-trajectory model (et_data(latents = )) ------------
#
# Given the trajectory X, individuals' latents are independent of each other,
# so each one is updated by random-walk Metropolis from its own terms: its
# observation weights at the timepoints that read the latent, and its
# transition probabilities into and out of them. The update runs after the
# iFFBS sweep in the same Gibbs step, so the parameters' blocks then see both.

augmented_latents_src <- function(d, helpers, add_payload, obs_call, reads) {
  if (!length(d$latents)) return(NULL)
  N <- d$n_individuals; T <- d$n_timepoints
  bindings <- character(); fns <- character(); arrays <- character(); kinds <- character()
  for (nm in names(d$latents)) {
    k <- d$latents[[nm]]
    p <- paste0("_etr_latent_", nm)
    if (is.numeric(k$init)) {
      bindings <- c(bindings, if (k$kind == "field") add_payload(p, k$init, "Float64")
                              else add_payload(p, as.numeric(k$init), "Float64", vector = TRUE))
    } else if (k$kind == "field") {
      bindings <- c(bindings, sprintf(
        "const %s = _etr_draw_coords(%s, (%s, %s, %s))", p, coord_prior_src(k$prior),
        julia_int(k$dim), julia_int(k$n), julia_int(N)))
    } else {
      bindings <- c(bindings, sprintf("const %s = rand(StableRNG(1), %s, %s)", p,
                                      dist_to_julia(k$prior), julia_int(N)))
    }
    arrays <- c(arrays, sprintf("%s = %s", nm, p))
    if (k$kind == "field") {
      fns <- c(fns, transpile_plain_fn(k$index, paste0("et_index_", nm), 1L, helpers,
                                       sprintf("latent '%s': `index`", nm)))
      kinds <- c(kinds, sprintf(
        "_et_field_sweep!(rng, model, DATA, X, latents.%s, %s, _et_windows(et_index_%s, %s, %s), %s, %s, %s)",
        nm, coord_prior_src(k$prior), nm, julia_int(k$n), julia_int(T), julia_symbol(nm),
        julia_num(reads[[nm]][["obs"]]), julia_num(reads[[nm]][["trans"]])))
    } else {
      kinds <- c(kinds, sprintf(
        "_et_individual_sweep!(rng, model, DATA, X, latents.%s, %s, %s, %s, %s)",
        nm, dist_to_julia(k$prior), julia_symbol(nm),
        julia_num(reads[[nm]][["obs"]]), julia_num(reads[[nm]][["trans"]])))
    }
  }
  scales <- vapply(names(d$latents), function(nm)
    sprintf("%s => %s", julia_symbol(nm), julia_float(d$latents[[nm]]$scale)), character(1))
  bindings <- c(bindings, sprintf("const latents = (; %s)", paste(arrays, collapse = ", ")))
  mh <- paste0(
"# Random-walk Metropolis for the latents, given the trajectory X. Scales adapt
# towards an acceptance rate of 0.44 while `_ET_ADAPT[]` sweeps remain.
const _ET_SCALE = Dict{Symbol,Float64}(", paste(scales, collapse = ", "), ")
const _ET_ADAPT = Ref(0)

_etr_draw_coords(prior, shape) = _etr_draw_coords((prior,), shape)
function _etr_draw_coords(prior::Tuple, shape)
    A = zeros(shape); rng = StableRNG(1)
    for c in 1:shape[1], e in 1:shape[2], i in 1:shape[3]
        A[c, e, i] = rand(rng, prior[min(c, length(prior))])
    end
    A
end
_et_coord_logpdf(prior, x, c) = logpdf(prior, x)
_et_coord_logpdf(prior::Tuple, x, c) = logpdf(prior[c], x)

# The timepoints whose observation may read element e: from the first t with
# index(t) == e to the first with index(t) == e + 1 (a move reads last year's).
function _et_windows(index, n, T)
    first_t = [something(findfirst(t -> index(t) == e, 1:T), 0) for e in 1:n]
    [first_t[e] == 0 ? (1:0) : first_t[e]:(e < n && first_t[e + 1] > 0 ? first_t[e + 1] : T)
     for e in 1:n]
end

# Individual i's log terms over timepoints lo..hi: observation weights there
# (when the latent is read by an observation), and the transitions out of
# lo - 1 .. hi (when read by a transition), clipped to its sampling period.
function _et_terms(model, DATA, X, i, lo, hi, obs, trans)
    a, b = DATA.sampling_period[i]
    ll = 0.0
    if obs
        for t in max(lo, a):min(hi, b)
            ll += log(", obs_call, ")
        end
    end
    if trans
        for t in max(lo - 1, a):min(hi, b - 1)
            ll += log(transition_prob(TRANS, model, DATA, X, i, t, X[t, i], X[t + 1, i]))
        end
    end
    ll
end

_et_accept(rng, l0, l1) = l0 == -Inf ? l1 > -Inf : log(rand(rng)) < l1 - l0

function _et_adapt!(name, accepted, proposed)
    _ET_ADAPT[] > 0 && proposed > 0 || return
    _ET_SCALE[name] *= exp(accepted / proposed - 0.44)
end

# An element no term reads, given X, has its prior as its conditional: a
# proposal that leaves its terms exactly unchanged shows that, and the element
# is then drawn from its prior. Only the others count towards adapting the
# step, since the prior-only ones accept almost anything and would push it far
# past what the constrained ones take.
_et_coord_rand(rng, prior, c) = rand(rng, prior)
_et_coord_rand(rng, prior::Tuple, c) = rand(rng, prior[c])

function _et_field_sweep!(rng, model, DATA, X, A, prior, windows, name, obs, trans)
    s = _ET_SCALE[name]; acc = 0; prop = 0
    old = zeros(size(A, 1))
    for i in axes(A, 3), e in axes(A, 2)
        w = windows[e]
        isempty(w) && continue
        lp(i, e) = sum(_et_coord_logpdf(prior, A[c, e, i], c) for c in axes(A, 1))
        t0 = _et_terms(model, DATA, X, i, first(w), last(w), obs, trans)
        l0 = t0 + lp(i, e)
        for c in axes(A, 1)
            old[c] = A[c, e, i]; A[c, e, i] += s * randn(rng)
        end
        t1 = _et_terms(model, DATA, X, i, first(w), last(w), obs, trans)
        if t1 == t0
            for c in axes(A, 1); A[c, e, i] = _et_coord_rand(rng, prior, c); end
            continue
        end
        prop += 1
        if _et_accept(rng, l0, t1 + lp(i, e))
            acc += 1
        else
            for c in axes(A, 1); A[c, e, i] = old[c]; end
        end
    end
    _et_adapt!(name, acc, prop)
end

function _et_individual_sweep!(rng, model, DATA, X, v, prior, name, obs, trans)
    s = _ET_SCALE[name]; acc = 0; prop = 0
    T = DATA.n_timepoints
    for i in eachindex(v)
        old = v[i]
        t0 = _et_terms(model, DATA, X, i, 1, T, obs, trans)
        v[i] = old + s * randn(rng)
        t1 = _et_terms(model, DATA, X, i, 1, T, obs, trans)
        if t1 == t0
            v[i] = rand(rng, prior)
            continue
        end
        prop += 1
        _et_accept(rng, t0 + logpdf(prior, old), t1 + logpdf(prior, v[i])) ? (acc += 1) : (v[i] = old)
    end
    _et_adapt!(name, acc, prop)
end

function et_latent_mh!(rng, model, DATA, X)
", paste0("    ", kinds, collapse = "\n"), "
    _ET_ADAPT[] = max(_ET_ADAPT[] - 1, 0)
    nothing
end")
  list(bindings = bindings, fns = fns, mh = mh)
}
