# Tier 1: the transpiler's supported subset. No Julia needed.

test_that("literals and symbols", {
  expect_equal(jlq(x), "x")
  expect_equal(jlq(1), "1")
  expect_equal(jlq(1L), "1")
  expect_equal(jlq(2.5), "2.5")
  expect_equal(jlq(TRUE), "true")
  expect_equal(jlq(FALSE), "false")
  # A whole-number double stays integral so it can still index or bound a loop.
  expect_equal(jlq(3.0), "3")
  expect_equal(jl("abc"), '"abc"')
  expect_equal(jl(NULL), "nothing")
})

test_that("non-integral doubles round-trip exactly", {
  for (v in c(0.1, 1/3, 1e-9, 6.022e23, pi)) {
    src <- jl(v)
    expect_equal(as.numeric(src), v, tolerance = 0,
                 info = paste("value", format(v, digits = 17)))
  }
})

test_that("arithmetic and comparison operators map across", {
  expect_equal(jlq(a + b), "a + b")
  expect_equal(jlq(a - b), "a - b")
  expect_equal(jlq(a * b), "a * b")
  expect_equal(jlq(a / b), "a / b")
  expect_equal(jlq(a^b), "a ^ b")
  expect_equal(jlq(a == b), "a == b")
  expect_equal(jlq(a != b), "a != b")
  expect_equal(jlq(a < b), "a < b")
  expect_equal(jlq(a <= b), "a <= b")
  expect_equal(jlq(a > b), "a > b")
  expect_equal(jlq(a >= b), "a >= b")
  expect_equal(jlq(a && b), "a && b")
  expect_equal(jlq(a || b), "a || b")
  expect_equal(jlq(a & b), "a & b")
  expect_equal(jlq(a | b), "a | b")
})

test_that("unary forms", {
  expect_equal(jlq(-x), "-x")
  expect_equal(jlq(+x), "+x")
  expect_equal(jlq(!x), "!x")
  expect_equal(jlq(-(a + b)), "-(a + b)")
})

test_that("grouping is preserved, so precedence cannot silently change", {
  expect_equal(jlq((a + b) * c), "(a + b) * c")
  expect_equal(jlq(a + b * c), "a + b * c")
  expect_equal(jlq(-expm1(-(x + y))), "-expm1(-(x + y))")
})

test_that("maths builtins map to their Julia names", {
  expect_equal(jlq(exp(x)), "exp(x)")
  expect_equal(jlq(log(x)), "log(x)")
  expect_equal(jlq(log1p(x)), "log1p(x)")
  expect_equal(jlq(expm1(x)), "expm1(x)")
  expect_equal(jlq(sqrt(x)), "sqrt(x)")
  expect_equal(jlq(abs(x)), "abs(x)")
  expect_equal(jlq(lgamma(x)), "lgamma(x)")
  expect_equal(jlq(length(x)), "length(x)")
  expect_equal(jlq(min(a, b)), "min(a, b)")
  expect_equal(jlq(max(a, b)), "max(a, b)")
  # R's `ceiling` is Julia's `ceil`, which is the kind of rename that has to be
  # in the table rather than assumed.
  expect_equal(jlq(ceiling(x)), "ceil(x)")
  expect_equal(jlq(floor(x)), "floor(x)")
})

test_that("shape queries become Julia's size/length", {
  expect_equal(jlq(nrow(A)), "size(A, 1)")
  expect_equal(jlq(ncol(A)), "size(A, 2)")
  expect_equal(jlq(dim(A)[3]), "size(A, 3)")
  expect_equal(jlq(seq_len(n)), "1:n")
  expect_equal(jlq(seq_along(v)), "eachindex(v)")
})

test_that("indexing passes through unchanged at every arity", {
  # Both languages are 1-based and column-major; this is the fact the whole
  # approach rests on, so it is asserted rather than assumed.
  expect_equal(jlq(x[i]), "x[i]")
  expect_equal(jlq(A[i, j]), "A[i, j]")
  expect_equal(jlq(B[i, j, k]), "B[i, j, k]")
  expect_equal(jlq(x[i + 1]), "x[i + 1]")
  expect_equal(jlq(A[data$g[i], t]), "A[data.g[i], t]")
})

test_that("field access, including chains and indexed fields", {
  expect_equal(jlq(data$x), "data.x")
  expect_equal(jlq(data$x[t, i]), "data.x[t, i]")
  expect_equal(jlq(data$aggregates$n_infected), "data.aggregates.n_infected")
  expect_equal(jlq(data$aggregates$n_infected[g, t]),
               "data.aggregates.n_infected[g, t]")
  expect_equal(jlq(data$sampling_period[i][1]), "data.sampling_period[i][1]")
})

test_that("model$x becomes model.x and is RECORDED", {
  ctx <- nc()
  expect_equal(tr(quote(model$alpha), ctx), "model.alpha")
  expect_equal(tr(quote(model$beta[g]), ctx), "model.beta[g]")
  expect_equal(tr(quote(model$nu[i, 1]), ctx), "model.nu[i, 1]")
  expect_setequal(ctx$reads$model, c("alpha", "beta", "nu"))
})

test_that("assignment", {
  expect_equal(jlq(x <- 1), "x = 1")
  # `=` as an assignment operator, built explicitly: `jlq(x = 1)` would be
  # parsed as a named argument to jlq() rather than as an expression.
  expect_equal(jl(call("=", as.symbol("x"), 1)), "x = 1")
  expect_equal(jlq(w[1] <- w[1] * 2), "w[1] = w[1] * 2")
})

test_that("return", {
  expect_equal(jlq(return(x)), "return x")
  expect_equal(jlq(return(a + b)), "return a + b")
})

test_that("if / else / else-if chains", {
  expect_equal(jl(quote(if (a > 0) b)), "if a > 0\n    b\nend")
  expect_equal(jl(quote(if (a > 0) b else c)),
               "if a > 0\n    b\nelse\n    c\nend")
  expect_equal(jl(quote(if (a) 1 else if (b) 2 else 3)),
               "if a\n    1\nelseif b\n    2\nelse\n    3\nend")
})

test_that("if with a braced body", {
  expect_equal(jl(quote(if (a) { x <- 1; y <- 2 })),
               "if a\n    x = 1\n    y = 2\nend")
})

test_that("for loops, next and break", {
  expect_equal(jl(quote(for (i in 1:n) s <- s + i)),
               "for i in 1:n\n    s = s + i\nend")
  expect_equal(jl(quote(for (i in 1:n) { if (i > 2) next })),
               "for i in 1:n\n    if i > 2\n        continue\n    end\nend")
  expect_equal(jl(quote(for (i in 1:n) { if (i > 2) break })),
               "for i in 1:n\n    if i > 2\n        break\n    end\nend")
})

test_that("ranges", {
  expect_equal(jlq(1:n), "1:n")
  expect_equal(jlq(1:(k - 1)), "1:(k - 1)")
})

test_that("ifelse becomes the Julia ternary", {
  expect_equal(jlq(ifelse(a > 0, b, c)), "(a > 0 ? b : c)")
})

test_that("blocks are newline-joined statements", {
  expect_equal(jl(quote({ a <- 1; b <- 2; a + b })), "a = 1\nb = 2\na + b")
})

test_that("vector allocation is typed by the ROLE, not the value", {
  # Outside a vector-returning role there is no parameter type in scope, so a
  # plain Float64 allocation is correct.
  expect_equal(jlq(rep(1, n)), "ones(Float64, n)")
  expect_equal(jlq(rep(0, n)), "zeros(Float64, n)")
  expect_equal(jlq(numeric(n)), "zeros(Float64, n)")
  # Inside one, it must follow the parameter scalar type.
  expect_equal(jlq(rep(1, n), vector_result = TRUE), "ones(ETR_T, n)")
  expect_equal(jlq(rep(0, n), vector_result = TRUE), "zeros(ETR_T, n)")
  expect_equal(jlq(numeric(n), vector_result = TRUE), "zeros(ETR_T, n)")
  expect_equal(jlq(rep(0.5, n), vector_result = TRUE), "fill(ETR_T(0.5), n)")
})

test_that("typed literals emit convert/one/zero, not a constructor call", {
  # `convert` and `one`/`zero` are defined for every numeric type in the AD
  # ecosystem, including the tracked types reverse-mode backends use. A type
  # CONSTRUCTOR call is not.
  expect_equal(jlq(et_one()), "one(ETR_T)")
  expect_equal(jlq(et_zero()), "zero(ETR_T)")
  expect_equal(jlq(et_num(x)), "convert(ETR_T, x)")
})

test_that("state comparisons become Symbols only in an aggregate context", {
  expect_equal(jlq(state == "I", state_syms = TRUE), "state == :I")
  expect_equal(jlq("I" == state, state_syms = TRUE), ":I == state")
  expect_equal(jlq(state != "D", state_syms = TRUE), "state != :D")
  # Elsewhere a string is just a string.
  expect_equal(jlq(state == "I"), 'state == "I"')
})

test_that("the R-side typed literals are the identity, so bodies still run in R", {
  expect_identical(et_one(), 1)
  expect_identical(et_zero(), 0)
  expect_identical(et_num(3.5), 3.5)
})
