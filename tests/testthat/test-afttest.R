test_that("test afttest", {
  ## Simulate data from an AFT model
  datgen <- function(n = 100) {
    z1 <- rbinom(n, 1, 0.5)
    z2 <- rnorm(n)
    e <- rnorm(n)
    tt <- exp(2 + z1 + z2 + e)
    cen <- runif(n, 0, 100)
    data.frame(Time = pmin(tt, cen), status = 1 * (tt < cen),
               z1 = z1, z2 = z2, id = 1:n)
  }
  set.seed(0)
  simdata = datgen(n = 20)
  
  X = simdata$Time
  D = simdata$status
  z1 = simdata$z1
  z2 = simdata$z2
  
  path = 200; testType = "link"; eqType = "mns"; optimType = "DFSANE"; form = "z2"
  
  result = afttest(formula = Surv(X, D) ~ z1 + z2, path = path, 
                   testType = testType, eqType = eqType, 
                   optimType = optimType, form = form, pathsave = 100)
  
  expect_equal(result$p_value, 0.585, tolerance=5e-2)
  expect_equal(result$p_std_value, 0.435, tolerance=5e-2)
  
  expect_equal(result$path, path)
  expect_equal(result$testType, testType)
  expect_equal(result$eqType, eqType)
  expect_equal(result$optimType, optimType)
})
