# Control constructors gb_*_ctrl(): fields, defaults and argument checks.


# ---- Constructors and defaults -----------------------------------------------

test_that("gb_learner_ctrl returns a well-formed list with expected fields", {
  lc <- gb_learner_ctrl()
  expect_true(is.list(lc))
  expect_true(all(c(
    "step_size", "initialise_intercept", "learners_to_fit_after_each_iter",
    "selection_loss_type", "use_step_size_for_selection",
    "use_learner_step_size_for_update", "knots", "tau",
    "preselected_learner_indices", "preselected_rhos",
    "excluded_learner_indices"
  ) %in% names(lc)))
})

test_that("gb_stop_ctrl returns a well-formed list with expected fields", {
  sc <- gb_stop_ctrl("none", max_iter = 10)
  expect_true(is.list(sc))
  expect_true(all(c(
    "stopping_method", "iter_no_improv", "iter_no_improv_rel", "true_f",
    "AIC_k",
    "max_iter", "max_runtime_mins"
  ) %in% names(sc)))
})

test_that("gb_edf_ctrl returns a well-formed list with expected fields", {
  ec <- gb_edf_ctrl()
  expect_true(is.list(ec))
  expect_true(all(c("edf_compute_method", "z", "cweights") %in% names(ec)))
})

test_that("every control has its class", {
  expect_s3_class(gb_learner_ctrl(), "gb_learner_ctrl")
  expect_s3_class(gb_stop_ctrl(), "gb_stop_ctrl")
  expect_s3_class(gb_edf_ctrl(), "gb_edf_ctrl")
  expect_s3_class(gb_param_boot_ctrl(), "gb_param_boot_ctrl")
  expect_s3_class(gb_bayesian_ctrl(), "gb_bayesian_ctrl")
  expect_s3_class(gb_uq_ctrl(), "gb_uq_ctrl")
})

test_that("gb_uq_ctrl holds the method and its control", {
  ic_none <- gb_uq_ctrl()
  expect_s3_class(ic_none, "gb_uq_ctrl")
  expect_named(ic_none, c("method", "control"))
  expect_equal(ic_none$method, "none")
  expect_null(ic_none$control)

  ic_boot <- gb_uq_ctrl("bootstrap")
  expect_named(ic_boot, c("method", "control"))
  expect_identical(ic_boot$control, gb_param_boot_ctrl())
  expect_named(ic_boot$control, c("n_samples", "mode", "alpha", "correction",
                                  "phi", "max_iter"))

  ic_bayes <- gb_uq_ctrl("bayes")
  expect_identical(ic_bayes$control, gb_bayesian_ctrl())
  expect_named(ic_bayes$control, c("n_samples", "n_burnin", "linear",
                                   "undersmooth", "edf_floor", "phi"))
  # One default per method, the same as the method's function
  expect_identical(ic_boot$control$n_samples,
                   formals(bootstrap_gbgam)$n_samples)
  expect_identical(ic_bayes$control$n_samples,
                   formals(bayesian_backfitting)$n_samples)
})

test_that("gb_uq_ctrl takes the settings as control or in ...", {
  pb <- gb_param_boot_ctrl(n_samples = 30L, mode = "free_selection",
                           correction = "smooth")
  expect_identical(gb_uq_ctrl("bootstrap", control = pb),
                   gb_uq_ctrl("bootstrap", n_samples = 30L,
                              mode = "free_selection",
                              correction = "smooth"))
  # The method is taken from the control
  expect_identical(gb_uq_ctrl(control = pb)$method, "bootstrap")
  expect_identical(gb_uq_ctrl(control = gb_bayesian_ctrl())$method,
                   "bayes")
  expect_error(gb_uq_ctrl(control = list(n_samples = 3)),
               "gb_param_boot_ctrl\\(\\) or gb_bayesian_ctrl")
  # Both, the wrong class, foreign or unnamed settings, and "none" with any
  expect_error(gb_uq_ctrl("bootstrap", control = pb, n_samples = 5L),
               "not both")
  expect_error(gb_uq_ctrl("bayes", control = pb),
               "gb_bayesian_ctrl\\(\\) for method = \"bayes\"")
  expect_error(gb_uq_ctrl("bootstrap", n_burnin = 5L),
               "n_burnin is not a setting of gb_param_boot_ctrl")
  expect_error(gb_uq_ctrl("bayes", mode = "fixed"),
               "mode is not a setting of gb_bayesian_ctrl")
  expect_error(gb_uq_ctrl("bootstrap", NULL, 30L), "must be named")
  expect_error(gb_uq_ctrl("none", n_samples = 10L), "no settings")
  expect_error(gb_uq_ctrl("none", control = pb), "no settings")
  # A stored control is rebuilt by do.call()
  for (ic in list(gb_uq_ctrl(), gb_uq_ctrl(control = pb),
                  gb_uq_ctrl("bayes", linear = "kept")))
    expect_identical(do.call("gb_uq_ctrl", ic), ic)
})

test_that("gb_param_boot_ctrl() checks the mode and max_iter", {
  expect_identical(gb_param_boot_ctrl()$mode, "fixed")
  for (md in c("fixed", "free_smoothing", "free_selection", "free_stopping"))
    expect_identical(gb_param_boot_ctrl(mode = md)$mode, md)
  expect_error(gb_param_boot_ctrl(mode = "free"), "should be one of")
  expect_error(gb_param_boot_ctrl(mode = character()), "length 1")
  expect_error(gb_param_boot_ctrl(mode = c("fixed", "free_smoothing")),
               "length 1")
  # The default vector of choices, as forwarded by bootstrap_gbgam()
  expect_identical(
    gb_param_boot_ctrl(mode = eval(formals(bootstrap_gbgam)$mode))$mode,
    "fixed"
  )
  expect_identical(gb_param_boot_ctrl(mode = "free_stopping",
                                      max_iter = 50)$max_iter, 50)
  expect_error(gb_param_boot_ctrl(mode = "free_selection", max_iter = 50),
               "only used with mode = \"free_stopping\"")
  expect_error(gb_param_boot_ctrl(mode = "free_stopping", max_iter = 0),
               "max_iter")
})

test_that("gb_edf_ctrl default logic is correct", {
  expect_equal(gb_edf_ctrl("hutchpp")$z, 80)
  expect_null(gb_edf_ctrl("exact")$z)
  expect_null(gb_edf_ctrl("none")$z)
})


# ---- Argument checks ---------------------------------------------------------

test_that("gb_uq_ctrl validates inputs", {
  expect_identical(gb_uq_ctrl()$method, "none")
  ic_b <- gb_uq_ctrl("bootstrap", n_samples = 50L)
  expect_identical(ic_b$method, "bootstrap")
  expect_identical(ic_b$control$n_samples, 50L)
  expect_identical(ic_b$control$alpha, 0.05)
  ic_bb <- gb_uq_ctrl("bayes", n_samples = 100L, n_burnin = 50L)
  expect_identical(ic_bb$method, "bayes")
  expect_identical(ic_bb$control[c("n_samples", "n_burnin")],
                   list(n_samples = 100L, n_burnin = 50L))
  # alpha = NULL is kept (no widening), not replaced by the default
  ic_0 <- gb_uq_ctrl("bootstrap", alpha = NULL)
  expect_true("alpha" %in% names(ic_0$control))
  expect_null(ic_0$control$alpha)
  expect_error(gb_uq_ctrl("bootstrap", n_samples = 0L),
               regexp = "n_samples")
  expect_error(gb_uq_ctrl("bayes", n_burnin = -1L),
               regexp = "n_burnin")
  expect_error(gb_uq_ctrl("bootstrap", alpha = 2), regexp = "alpha")
  # alpha must lie strictly inside (0, 1)
  expect_error(gb_uq_ctrl("bootstrap", alpha = 0), regexp = "alpha")
  expect_error(gb_uq_ctrl("bootstrap", alpha = 1), regexp = "alpha")
  expect_error(gb_uq_ctrl("bootstrap", alpha = NA_real_),
               regexp = "alpha")
  # cross-method argument errors
  expect_error(gb_uq_ctrl("bootstrap", n_burnin = 5L), regexp = "n_burnin")
  expect_error(gb_uq_ctrl("bayes", alpha = 0.05),      regexp = "alpha")
  expect_error(gb_uq_ctrl("bayes", mode = "fixed"),    regexp = "mode")
  expect_error(gb_uq_ctrl("none", n_samples = 10L),    regexp = "none")
})

test_that("gb_stop_ctrl enforces logical stopping conditions", {
  expect_error(gb_stop_ctrl("none"),        "No stopping condition provided")
  expect_error(gb_stop_ctrl("true_f", max_iter = 10), "true_f is required")

  expect_no_error(gb_stop_ctrl("aic",  max_iter = 10))
  expect_no_error(gb_stop_ctrl("none", max_iter = 10))
})

test_that("missing and non-positive control values are rejected", {
  expect_error(gb_stop_ctrl(max_iter = NA), "max_iter")
  expect_error(gb_stop_ctrl(iter_no_improv = NA), "iter_no_improv")
  expect_error(gb_learner_ctrl(step_size = NA_real_), "step_size")
  expect_error(gb_learner_ctrl(step_size = Inf), "step_size")
  expect_error(gb_learner_ctrl(step_size = c(0.1, 0.2)), "step_size")
  expect_error(gb_learner_ctrl(step_size = 0), "strictly positive")
  expect_error(gb_learner_ctrl(step_size = -0.1), "step_size")
  # The step_size > 1 warning is given once, by gbgam() (tested below)
  expect_no_warning(gb_learner_ctrl(step_size = 1.5))
  expect_no_error(gb_learner_ctrl(step_size = 1))
})

test_that("gb_learner_ctrl validates excluded_learner_indices", {
  expect_no_error(gb_learner_ctrl(excluded_learner_indices = NULL))
  expect_no_error(gb_learner_ctrl(excluded_learner_indices = c(2L, 3L)))
  expect_error(gb_learner_ctrl(excluded_learner_indices = 0L),  regexp = "excluded_learner_indices")
  expect_error(gb_learner_ctrl(excluded_learner_indices = -1L), regexp = "excluded_learner_indices")
  expect_error(gb_learner_ctrl(excluded_learner_indices = "a"), regexp = "excluded_learner_indices")
})

test_that("gb_learner_ctrl validates preselected path arguments", {
  expect_error(
    gb_learner_ctrl(preselected_rhos = rep(0, 5)),
    "supply preselected_learner_indices"
  )

  expect_error(
    gb_learner_ctrl(preselected_learner_indices = rep(1L, 5),
                    preselected_rhos = rep(0, 4)),
    "must have the same length"
  )

  expect_no_error(
    gb_learner_ctrl(preselected_learner_indices = rep(1L, 5),
                    preselected_rhos = rep(0, 5))
  )
})

test_that("iter_no_improv must be at least 1", {
  expect_error(gb_stop_ctrl(iter_no_improv = 0), "iter_no_improv")
  expect_no_error(gb_stop_ctrl(iter_no_improv = 1))
})

test_that("the patience defaults to max(50, ceiling(0.1 * m*))", {
  sc <- gb_stop_ctrl()
  expect_identical(sc$iter_no_improv, 50)
  expect_identical(sc$iter_no_improv_rel, 0.1)
  expect_identical(gb_stop_ctrl("cv")$iter_no_improv_rel, 0.1)
  expect_null(gb_stop_ctrl("none", max_iter = 5)$iter_no_improv_rel)
  expect_identical(gb_stop_ctrl(iter_no_improv_rel = 0)$iter_no_improv_rel,
                   0)
})

test_that("iter_no_improv_rel must be a finite non-negative number", {
  for (bad in list(-0.1, Inf, NA_real_, c(0.1, 0.2), "a")) {
    expect_error(gb_stop_ctrl(iter_no_improv_rel = bad),
                 "iter_no_improv_rel")
  }
  expect_no_error(gb_stop_ctrl(iter_no_improv_rel = 2.5))
})

test_that("patience() is max(iter_no_improv, ceiling(rel * m_best))", {
  mod <- list(iter_no_improv = 50, iter_no_improv_rel = 0.1)
  expect_identical(patience(mod, 0), 50)
  expect_identical(patience(mod, 500), 50)
  expect_identical(patience(mod, 501), 51)
  expect_identical(patience(mod, 1e4), 1000)
  # Rounding error in the product does not add an iteration
  expect_identical(patience(list(iter_no_improv = 1,
                                 iter_no_improv_rel = 0.07), 100), 7)
  expect_identical(patience(list(iter_no_improv = 5,
                                 iter_no_improv_rel = 0), 1e4), 5)
  expect_error(patience(list(iter_no_improv = 5), 10))
})

test_that("unused arguments of the controls are errors", {
  expect_error(gb_stop_ctrl(true_f = 1:3), "true_f is not used")
  expect_error(gb_stop_ctrl("cv", AIC_k = 3), "AIC_k is not used")
  expect_error(gb_stop_ctrl("none", max_iter = 5, iter_no_improv = 5),
               "iter_no_improv is not used")
  expect_error(gb_stop_ctrl("none", max_iter = 5, iter_no_improv_rel = 0),
               "iter_no_improv_rel is not used")
  expect_error(gb_edf_ctrl("exact", z = 20), "z is only used")
  expect_error(gb_edf_ctrl("none", cweights = c(1, 1, 1)),
               "cweights is only used")
  expect_identical(gb_edf_ctrl("hutchpp")$cweights, c(1, 2, 1))
  expect_null(gb_edf_ctrl("exact")$cweights)
  # Stored control lists are re-validated without errors
  for (ctl in list(gb_stop_ctrl("none", max_iter = 5), gb_stop_ctrl(),
                   gb_stop_ctrl("true_f", true_f = 1:3),
                   gb_stop_ctrl("cv")))
    expect_identical(do.call("gb_stop_ctrl", ctl), ctl)
  for (ctl in list(gb_edf_ctrl(), gb_edf_ctrl("exact"), gb_edf_ctrl("none")))
    expect_identical(do.call("gb_edf_ctrl", ctl), ctl)
})

test_that("step_size > 1 warns once, without the call", {
  d <- core_data()
  ws <- list()
  withCallingHandlers(
    gbgam(ff_core, data = d, verbose = FALSE,
          learner_control = gb_learner_ctrl(step_size = 1.2),
          stop_control = gb_stop_ctrl("none", max_iter = 3),
          edf_control = gb_edf_ctrl("none")),
    warning = function(w) {
      ws[[length(ws) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })
  expect_length(ws, 1L)
  expect_match(conditionMessage(ws[[1L]]), "step_size > 1 is not recommended")
  expect_null(conditionCall(ws[[1L]]))
})

test_that("a missing z is rejected with an error naming z", {
  expect_error(gb_edf_ctrl(z = NA), "z \\(Hutch\\+\\+ samples\\)")
  expect_error(gb_edf_ctrl(z = NA_integer_), "Contains missing values")
})

test_that("cweights must be positive and give non-empty Hutch++ blocks", {
  expect_error(gb_edf_ctrl(cweights = c(0, 1, 1)), "strictly positive")
  expect_error(gb_edf_ctrl(cweights = c(1, 1, 0)), "strictly positive")
  expect_error(gb_edf_ctrl("hutchpp", z = 3, cweights = c(1, 10, 1)),
               "at least one probe vector")
  expect_error(gb_edf_ctrl("hutchpp", z = 3, cweights = c(1e-9, 1, 1)),
               "at least one probe vector")
  m  <- matrix(seq_len(12), 2, 6)
  sp <- split_mat_cols(m, c(0, 3, 3))
  expect_equal(ncol(sp[[1]]), 0L)
  expect_equal(sp[[2]], m[, 1:3])
  expect_equal(sp[[3]], m[, 4:6])
})

test_that("gb_uq_ctrl() requires n_samples >= 2 for both methods", {
  expect_error(gb_uq_ctrl("bootstrap", n_samples = 1L), "n_samples")
  expect_error(gb_uq_ctrl("bayes", n_samples = 1L), "n_samples")
  expect_identical(
    gb_uq_ctrl("bootstrap", n_samples = 2L)$control$n_samples, 2L)
  expect_identical(
    gb_uq_ctrl("bayes", n_samples = 2L)$control$n_samples, 2L)
})

test_that("gb_uq_ctrl() validates the Bayesian options", {
  ctl <- gb_uq_ctrl("bayes", linear = "kept", undersmooth = TRUE,
                    edf_floor = 2)
  expect_identical(ctl$control$linear, "kept")
  expect_true(ctl$control$undersmooth)
  expect_identical(ctl$control$edf_floor, 2)
  expect_identical(gb_uq_ctrl("bayes")$control$edf_floor, 1)
  expect_error(gb_uq_ctrl("bayes", edf_floor = 0), "edf_floor")
  expect_error(gb_uq_ctrl("bayes", edf_floor = NA), "edf_floor")
  expect_error(gb_uq_ctrl("bootstrap", edf_floor = 1),
               "not a setting of gb_param_boot_ctrl")
  expect_identical(do.call("gb_uq_ctrl", ctl), ctl)
  expect_error(gb_uq_ctrl("bootstrap", linear = "kept"),
               "not a setting of gb_param_boot_ctrl")
  expect_error(gb_uq_ctrl("bayes", correction = "smooth"),
               "not a setting of gb_bayesian_ctrl")
})


# ---- Bootstrap modes ---------------------------------------------------------

test_that("the stored bootstrap control holds the mode", {
  ctl <- gb_uq_ctrl("bootstrap", mode = "free_stopping",
                    correction = "learner", phi = 0.2)
  expect_identical(ctl$control$mode, "free_stopping")
  expect_identical(ctl$control$correction, "learner")
  expect_identical(do.call("gb_uq_ctrl", ctl), ctl)
  expect_identical(names(boot_modes),
                   c("fixed", "free_smoothing", "free_selection",
                     "free_stopping"))
})
