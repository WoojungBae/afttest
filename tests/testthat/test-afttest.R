test_that("test afttest", {
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
  simdata = datgen(n = 100)
  
  result = afttest(object = Surv(Time, status) ~ z1 + z2, data = simdata, 
                   npath = 100, testType = "covform", estMethod = "rr", 
                   eqType = "ns", cov.tested = "z2", npathsave = 50)
  
  expect_equal(result$p_value, 0.03, tolerance=5e-2)
  expect_equal(result$p_std_value, 0.05, tolerance=5e-2)
  
  # plot(result, std = TRUE)
  # plot(result, std = FALSE)
})