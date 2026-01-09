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
                   linApprox = TRUE)
  expect_equal(result$p_value, 0.00, tolerance=1e-1)
  expect_equal(result$p_std_value, 0.00, tolerance=1e-1)

  result = afttest(object = Surv(Time, status) ~ z1 + z2, data = simdata,
                   npath = 100, testType = "covForm", estMethod = "rr",
                   eqType = "is", covTested = "z2", npathsave = 50,
                   linApprox = TRUE)
  expect_equal(result$p_value, 0.00, tolerance=1e-1)
  expect_equal(result$p_std_value, 0.00, tolerance=1e-1)

  result = afttest(object = Surv(Time, status) ~ z1 + z2, data = simdata,
                   npath = 100, testType = "covForm", estMethod = "ls",
                   eqType = "ls", covTested = "z2", npathsave = 50,
                   linApprox = TRUE)
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
                   linApprox = FALSE)
  expect_equal(result$p_value, 0.00, tolerance=1e-1)
  expect_equal(result$p_std_value, 0.00, tolerance=1e-1)

  result = afttest(object = Surv(Time, status) ~ z1 + z2, data = simdata,
                   npath = 100, testType = "covForm", estMethod = "rr",
                   eqType = "is", covTested = "z2", npathsave = 50,
                   linApprox = FALSE)
  expect_equal(result$p_value, 0.00, tolerance=1e-1)
  expect_equal(result$p_std_value, 0.00, tolerance=1e-1)

  result = afttest(object = Surv(Time, status) ~ z1 + z2, data = simdata,
                   npath = 100, testType = "covForm", estMethod = "ls",
                   eqType = "ls", covTested = "z2", npathsave = 50,
                   linApprox = FALSE)
  expect_equal(result$p_value, 0.01, tolerance=1e-1)
  expect_equal(result$p_std_value, 0.00, tolerance=1e-1)
})