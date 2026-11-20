# Every mgcv basis type usable as a t2() marginal (soap films are not: they
# carry two penalties, t2() allows one per marginal) builds and fits.

# ---- Local data generators ---------------------------------------------------

make_basis_scalar <- function(n = 80L, seed = 1000L) {
  set.seed(seed)
  x <- sort(stats::runif(n))
  y <- sin(2 * pi * x) + stats::rnorm(n, sd = 0.3)
  data.frame(y = y, x = x)
}

make_basis_bivariate <- function(n = 100L, seed = 1001L) {
  set.seed(seed)
  x <- stats::runif(n)
  z <- stats::runif(n)
  y <- sin(2 * pi * x) * cos(pi * z) + stats::rnorm(n, sd = 0.3)
  data.frame(y = y, x = x, z = z)
}

make_cyclic_data <- function(n = 80L, seed = 1002L) {
  set.seed(seed)
  # Predictor on [0, 2*pi]; the periodic boundary condition matches the range
  x <- sort(stats::runif(n, 0, 2 * pi))
  y <- cos(x) + stats::rnorm(n, sd = 0.3)
  data.frame(y = y, x = x)
}

# Shared minimal control arguments for fast initialisation tests
lc_basis <- gb_learner_ctrl(step_size = 0.1)
sc_basis <- gb_stop_ctrl("none", max_iter = 8L)
ec_basis <- gb_edf_ctrl("none")
ic_basis <- gb_uq_ctrl("none")

fit_basis_mod <- function(formula, data, knots = NULL, ...) {
  lc <- if (is.null(knots)) {
    lc_basis
  } else {
    gb_learner_ctrl(step_size = 0.1, knots = knots)
  }
  gbgam(
    formula,
    data            = data,
    family          = stats::gaussian(),
    learner_control = lc,
    stop_control    = sc_basis,
    edf_control     = ec_basis,
    uq_control      = ic_basis,
    verbose         = FALSE,
    ...
  )
}


# ---- Standard spline bases: tp, ts, bs, ps, gp -------------------------------
# These are all standard penalised-regression spline types that differ only in
# how the basis and penalty are constructed.  All are valid marginal bases for
# a univariate t2() or a bivariate tensor product.

test_that("univariate t2() initialises for standard spline bases: tp, ts, bs, ps, gp", {
  skip_on_cran()
  dat <- make_basis_scalar()
  for (bs in c("tp", "ts", "bs", "ps", "gp")) {
    mod <- fit_basis_mod(y ~ t2(x, k = 8L, bs = bs), data = dat)
    expect_model_contract(mod, nrow(dat))
  }
})

test_that("bivariate t2() tensor product initialises for standard spline bases: tp, ts, bs, ps", {
  skip_on_cran()
  dat <- make_basis_bivariate()
  for (bs in c("tp", "ts", "bs", "ps")) {
    mod <- fit_basis_mod(
      y ~ t2(x, z, k = c(5L, 5L), bs = bs),
      data = dat
    )
    expect_model_contract(mod, nrow(dat))
  }
})


# ---- Duchon splines (ds) -----------------------------------------------------
# Duchon splines generalise thin-plate splines by allowing fractional smoothness
# orders.  They require no special data setup beyond a scalar/vector predictor.

test_that("univariate t2(bs='ds') Duchon splines initialise", {
  skip_on_cran()
  dat <- make_basis_scalar()
  mod <- fit_basis_mod(y ~ t2(x, k = 8L, bs = "ds"), data = dat)
  expect_model_contract(mod, nrow(dat))
})

test_that("bivariate t2(bs=c('ds','ds')) Duchon tensor product initialises", {
  skip_on_cran()
  dat <- make_basis_bivariate()
  mod <- fit_basis_mod(
    y ~ t2(x, z, k = c(5L, 5L), bs = c("ds", "ds")),
    data = dat
  )
  expect_model_contract(mod, nrow(dat))
})


# ---- Cubic spline with shrinkage (cs) ----------------------------------------
# Like cr but augmented with an extra penalty that shrinks the entire smooth to
# zero, enabling exact zero selection.

test_that("univariate t2(bs='cs') cubic-regression-with-shrinkage initialises", {
  skip_on_cran()
  dat <- make_basis_scalar()
  mod <- fit_basis_mod(y ~ t2(x, k = 8L, bs = "cs"), data = dat)
  expect_model_contract(mod, nrow(dat))
})


# ---- Cyclic cubic splines (cc) -----------------------------------------------
# Enforces periodic boundary conditions.  The knots argument must specify the
# period endpoints so mgcv can match the first and last basis functions.

test_that("univariate t2(bs='cc') cyclic cubic splines initialise on periodic data", {
  skip_on_cran()
  dat <- make_cyclic_data()
  mod <- fit_basis_mod(
    y ~ t2(x, k = 8L, bs = "cc"),
    data  = dat,
    knots = list(x = c(0, 2 * pi))
  )
  expect_model_contract(mod, nrow(dat))
})

test_that("bivariate t2() with cc cyclic marginal and cr scalar marginal initialises", {
  skip_on_cran()
  dat      <- make_cyclic_data()
  dat$z    <- stats::rnorm(nrow(dat))
  mod <- fit_basis_mod(
    y ~ t2(x, z, k = c(6L, 5L), bs = c("cc", "cr")),
    data  = dat,
    knots = list(x = c(0, 2 * pi))
  )
  expect_model_contract(mod, nrow(dat))
})


# ---- Cyclic P-splines (cp) ---------------------------------------------------
# P-spline analogue of cc: periodic B-spline basis with cyclic difference
# penalty.

test_that("univariate t2(bs='cp') cyclic P-splines initialise on periodic data", {
  skip_on_cran()
  dat <- make_cyclic_data()
  mod <- fit_basis_mod(
    y ~ t2(x, k = 8L, bs = "cp"),
    data  = dat,
    knots = list(x = c(0, 2 * pi))
  )
  expect_model_contract(mod, nrow(dat))
})


# ---- Spherical splines (sos) -------------------------------------------------
# Isotropic smooths on the sphere, parameterised by latitude and longitude.
# Using d=2 tells t2() that both lat and lon together form a single bivariate
# marginal component (the spherical basis is inherently 2D).

test_that("t2(lat, lon, bs='sos', d=2) spherical splines initialise for lat/lon data", {
  skip_on_cran()
  set.seed(1003L)
  n   <- 150L
  lat <- stats::runif(n, -80, 80)
  lon <- stats::runif(n, -180, 180)
  y   <- sin(lat * pi / 180) + stats::rnorm(n, sd = 0.3)
  dat <- data.frame(y = y, lat = lat, lon = lon)

  mod <- fit_basis_mod(
    y ~ t2(lat, lon, d = 2L, bs = "sos"),
    data = dat
  )
  expect_model_contract(mod, n)
})

test_that("space-time t2() with sos spherical and cr temporal marginals initialises", {
  # d=c(2, 1) groups lat/lon as one bivariate spatial component and t as the
  # independent temporal marginal.
  skip_on_cran()
  set.seed(1004L)
  n   <- 200L
  lat <- stats::runif(n, -80, 80)
  lon <- stats::runif(n, -180, 180)
  t   <- stats::runif(n)
  y   <- sin(lat * pi / 180) * (0.5 + t) + stats::rnorm(n, sd = 0.3)
  dat <- data.frame(y = y, lat = lat, lon = lon, t = t)

  mod <- fit_basis_mod(
    y ~ t2(lat, lon, t, d = c(2L, 1L), bs = c("sos", "cr"), k = c(6L, 5L)),
    data = dat
  )
  expect_model_contract(mod, n)
})


# ---- Random effects (re) -----------------------------------------------------
# "re" treats factor levels as i.i.d. Gaussian random effects.
# t2(group, bs="re") gives random intercepts; t2(x, group, bs=c("cr","re"))
# gives group-level random slopes for the smooth of x.

test_that("t2(group, bs='re') random intercepts initialise", {
  skip_on_cran()
  set.seed(1005L)
  n     <- 120L
  group <- factor(rep(letters[1:5], each = 24L))
  x     <- stats::runif(n)
  alpha <- c(a = 0, b = 0.5, c = -0.3, d = 0.8, e = -0.2)
  y     <- 0.5 * x + alpha[as.character(group)] + stats::rnorm(n, sd = 0.2)
  dat   <- data.frame(y = y, x = x, group = group)

  mod <- fit_basis_mod(
    y ~ t2(x, k = 5L, bs = "cr") + t2(group, bs = "re"),
    data = dat
  )
  expect_model_contract(mod, n)
})

test_that("t2(x, group, bs=c('cr','re')) random slopes initialise", {
  skip_on_cran()
  set.seed(1006L)
  n     <- 150L
  group <- factor(rep(letters[1:5], each = 30L))
  x     <- stats::runif(n)
  slope <- c(a = 0.5, b = 1.0, c = -0.5, d = 1.5, e = -1.0)
  y     <- slope[as.character(group)] * x + stats::rnorm(n, sd = 0.3)
  dat   <- data.frame(y = y, x = x, group = group)

  mod <- fit_basis_mod(
    y ~ t2(x, group, bs = c("cr", "re")),
    data = dat
  )
  expect_model_contract(mod, n)
})


# ---- Markov random field (mrf) -----------------------------------------------
# MRF smooths are defined over a graph of discrete areas.
# Here we use a simple 5-node chain graph: 1--2--3--4--5.

make_chain_penalty <- function(K = 5L, labels = letters[seq_len(K)]) {
  S <- matrix(0, K, K, dimnames = list(labels, labels))
  for (i in seq_len(K - 1L)) {
    S[i, i + 1L] <- -1
    S[i + 1L, i] <- -1
  }
  diag(S) <- -rowSums(S)
  S
}

test_that("t2(area, bs='mrf') Markov random field splines initialise", {
  skip_on_cran()
  K <- 5L
  S <- make_chain_penalty(K)

  set.seed(1007L)
  n    <- 200L
  area <- factor(sample(letters[seq_len(K)], n, replace = TRUE),
                 levels = letters[seq_len(K)])
  alpha <- c(a = 0, b = 0.5, c = -0.3, d = 0.8, e = -0.2)
  y    <- alpha[as.character(area)] + stats::rnorm(n, sd = 0.3)
  dat  <- data.frame(y = y, area = area)

  mod <- fit_basis_mod(
    y ~ t2(area, bs = "mrf", xt = list(penalty = S)),
    data = dat
  )
  expect_model_contract(mod, n)
})

test_that("space-time t2() with mrf spatial and cr temporal marginals initialises", {
  # xt for multi-marginal t2() is a list of per-marginal xt objects;
  # NULL (or an empty list) is passed for marginals that need no xt.
  skip_on_cran()
  K <- 5L
  S <- make_chain_penalty(K)

  set.seed(1008L)
  n    <- 200L
  area <- factor(sample(letters[seq_len(K)], n, replace = TRUE),
                 levels = letters[seq_len(K)])
  t    <- stats::runif(n)
  alpha <- c(a = 0, b = 0.5, c = -0.3, d = 0.8, e = -0.2)
  y    <- alpha[as.character(area)] * (0.5 + t) + stats::rnorm(n, sd = 0.3)
  dat  <- data.frame(y = y, area = area, t = t)

  mod <- fit_basis_mod(
    y ~ t2(area, t, k = c(K, 5L), bs = c("mrf", "cr"),
           xt = list(list(penalty = S), NULL)),
    data = dat
  )
  expect_model_contract(mod, n)
})

simulate_by_data <- function(n = 300L, seed = 11L) {
  set.seed(seed)
  dat <- data.frame(
    x   = stats::runif(n),
    z   = stats::runif(n, 0.5, 1.5),
    fac = factor(sample(c("a", "b", "c"), n, replace = TRUE))
  )
  dat$f_fac <- 2 * sin(2 * pi * dat$x) * (dat$fac == "a") + as.numeric(dat$fac)
  dat$f_num <- dat$z * sin(2 * pi * dat$x)
  dat$y     <- dat$f_fac + stats::rnorm(n, sd = 0.3)
  dat$y2    <- dat$f_num + stats::rnorm(n, sd = 0.3)
  dat
}


# ---- By-variables ------------------------------------------------------------

test_that("t2(x, by = factor) keeps a smooth per level and predicts correctly", {
  dat <- simulate_by_data()
  mod <- gbgam(
    y ~ fac + t2(x, by = fac, k = 8, bs = "cr"), data = dat,
    learner_control = gb_learner_ctrl(step_size = 0.1), verbose = FALSE
  )
  labels <- vapply(mod$learners, \(l) l$var_label, "")
  for (lev in levels(dat$fac)) {
    expect_true(any(grepl(paste0("t2(x):fac", lev), labels, fixed = TRUE)),
                info = lev)
  }
  # Penalised + linear nullspace learner per level (the constant column
  # within a level is aliased with the fac main effect and is dropped)
  expect_equal(sum(grepl("t2(x):fac", labels, fixed = TRUE)), 6L)
  expect_length(mod$ind0var, nlevels(dat$fac))

  expect_equal(as.vector(predict(mod, newdata = dat)),
               as.vector(mod$linear.predictors), tolerance = 1e-10)
  expect_gt(stats::cor(mod$fitted.values, dat$f_fac), 0.99)
})

test_that("t2(x, by = numeric) fits and predicts correctly", {
  dat <- simulate_by_data()
  mod <- gbgam(
    y2 ~ t2(x, by = z, k = 8, bs = "cr"), data = dat,
    learner_control = gb_learner_ctrl(step_size = 0.1), verbose = FALSE
  )
  labels <- vapply(mod$learners, \(l) l$var_label, "")
  expect_true(any(grepl("t2(x):z", labels, fixed = TRUE)))
  expect_true(any(vapply(mod$learners, inherits, NA, "demmler_reinsch")))
  expect_equal(as.vector(predict(mod, newdata = dat)),
               as.vector(mod$linear.predictors), tolerance = 1e-10)
  expect_gt(stats::cor(mod$fitted.values, dat$f_num), 0.99)
})
