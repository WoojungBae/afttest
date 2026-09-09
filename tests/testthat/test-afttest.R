test_that("afttest linApprox=TRUE runs correctly", {
  datgen <- function(n = 100) {
    z1 <- rbinom(n, 1, 0.5)
    z2 <- rnorm(n)
    e <- rnorm(n)
    tt <- exp(2 + z1 + z2 + 0.5*z2^{2}+ e)
    cen <- runif(n, 0, 100)
    data.frame(Time = pmin(tt, cen), status = 1 * (tt < cen),
               z1 = z1, z2 = z2, id = 1:n)
  }
  set.seed(1)
  simdata = datgen(300)
  
  # linApprox = TRUE
  result = afttest(object = Surv(Time, status) ~ z1 + z2, data = simdata,
                   npath = 100, testType = "covForm", estMethod = "rr",
                   eqType = "ns", covTested = "z2", npathsave = 50,
                   linApprox = TRUE, seed = 1)
  expect_equal(result$p_value, 0.00, tolerance=1e-1)
  expect_equal(result$p_std_value, 0.00, tolerance=1e-1)
  
  result = afttest(object = Surv(Time, status) ~ z1 + z2, data = simdata,
                   npath = 100, testType = "covForm", estMethod = "rr",
                   eqType = "is", covTested = "z2", npathsave = 50,
                   linApprox = TRUE, seed = 1)
  expect_equal(result$p_value, 0.00, tolerance=1e-1)
  expect_equal(result$p_std_value, 0.00, tolerance=1e-1)
  
  result = afttest(object = Surv(Time, status) ~ z1 + z2, data = simdata,
                   npath = 100, testType = "covForm", estMethod = "ls",
                   eqType = "ls", covTested = "z2", npathsave = 50,
                   linApprox = TRUE, seed = 1)
  expect_equal(result$p_value, 0.01, tolerance=1e-1)
  expect_equal(result$p_std_value, 0.00, tolerance=1e-1)
})


test_that("afttest linApprox=FALSE runs correctly", {
  # This block is slow, so we SKIP it on CRAN.
  testthat::skip_on_cran()
  
  datgen <- function(n = 100) {
    z1 <- rbinom(n, 1, 0.5)
    z2 <- rnorm(n)
    e <- rnorm(n)
    tt <- exp(2 + z1 + z2 + 0.5*z2^{2}+ e)
    cen <- runif(n, 0, 100)
    data.frame(Time = pmin(tt, cen), status = 1 * (tt < cen),
               z1 = z1, z2 = z2, id = 1:n)
  }
  set.seed(1)
  simdata = datgen(300)
  
  # linApprox = FALSE
  result = afttest(object = Surv(Time, status) ~ z1 + z2, data = simdata,
                   npath = 100, testType = "covForm", estMethod = "rr",
                   eqType = "ns", covTested = "z2", npathsave = 50,
                   linApprox = FALSE, seed = 1)
  expect_equal(result$p_value, 0.00, tolerance=1e-1)
  expect_equal(result$p_std_value, 0.00, tolerance=1e-1)
  
  result = afttest(object = Surv(Time, status) ~ z1 + z2, data = simdata,
                   npath = 100, testType = "covForm", estMethod = "rr",
                   eqType = "is", covTested = "z2", npathsave = 50,
                   linApprox = FALSE, seed = 1)
  expect_equal(result$p_value, 0.00, tolerance=1e-1)
  expect_equal(result$p_std_value, 0.00, tolerance=1e-1)
  
  result = afttest(object = Surv(Time, status) ~ z1 + z2, data = simdata,
                   npath = 100, testType = "covForm", estMethod = "ls",
                   eqType = "ls", covTested = "z2", npathsave = 50,
                   linApprox = FALSE, seed = 1)
  expect_equal(result$p_value, 0.01, tolerance=1e-1)
  expect_equal(result$p_std_value, 0.00, tolerance=1e-1)
})

test_that("formula input supports factor covariates and saved-path limits", {
  set.seed(11)
  dat <- data.frame(
    time = rexp(60) + 0.1,
    status = rbinom(60, 1, 0.7),
    group = factor(rep(c("a", "b", "c"), 20)),
    x = rnorm(60)
  )

  fit <- NULL
  expect_warning(
    fit <- afttest(
      Surv(time, status) ~ group + x,
      data = dat,
      npath = 50,
      npathsave = 100,
      seed = 1
    ),
    "npathsave exceeds npath"
  )

  expect_s3_class(fit, "afttest")
  expect_identical(fit$npathsave, 50L)
  expect_length(fit$apprx_npath, 50)
  expect_gt(fit$p_value, 0)
  expect_lt(fit$p_value, 1)
  expect_gt(fit$p_std_value, 0)
  expect_lt(fit$p_std_value, 1)
  unstd_max <- vapply(fit$apprx_npath, function(path) max(abs(path)), 0.0)
  std_max <- vapply(fit$apprx_std_npath, function(path) max(abs(path)), 0.0)
  expect_equal(
    fit$p_value,
    (sum(unstd_max >= max(abs(fit$obs_npath))) + 0.5) / 51
  )
  expect_equal(
    fit$p_std_value,
    (sum(std_max >= max(abs(fit$obs_std_npath))) + 0.5) / 51
  )
  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(fit, npath = 50))
})

test_that("invalid resampling arguments fail before model fitting", {
  dat <- data.frame(time = 1:60, status = rep(1, 60), x = 1:60)

  expect_error(
    afttest(Surv(time, status) ~ x, data = dat, npath = 49),
    "npath must be a single integer"
  )
  expect_error(
    afttest(Surv(time, status) ~ x, data = dat, npathsave = -1),
    "npathsave must be a single nonnegative integer"
  )
  expect_error(
    afttest(Surv(time, status) ~ x, data = dat, linApprox = NA),
    "linApprox must be a single logical value"
  )
})

test_that("a fixed seed reproduces the diagnostic calculation", {
  set.seed(21)
  dat <- data.frame(
    time = rexp(60) + 0.1,
    status = rbinom(60, 1, 0.7),
    x = rnorm(60)
  )

  fit1 <- afttest(
    Surv(time, status) ~ x,
    data = dat,
    npath = 50,
    npathsave = 0,
    seed = 99
  )
  fit2 <- afttest(
    Surv(time, status) ~ x,
    data = dat,
    npath = 50,
    npathsave = 0,
    seed = 99
  )

  expect_identical(fit1$p_value, fit2$p_value)
  expect_identical(fit1$p_std_value, fit2$p_std_value)
})
