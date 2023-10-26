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
simdata <- datgen(n = 20)

X <- simdata$Time
D <- simdata$status
z1 <- simdata$z1
z2 <- simdata$z2

result <- afttest(Surv(X, D) ~ z1 + z2, testType="link", eqType="mns")
afttestplot(result)