
#############################################
# Doubly Robust Estimation under NMAR
# Reproducibility Code (Single File)
#############################################

# Required libraries
library(KernSmooth)

#############################################
# 1. Data Generation under NMAR
#############################################
generate_data <- function(n, beta = c(1, -1), m_fun, alpha) {
  X1 <- runif(n, -1, 1)
  X2 <- runif(n, -1, 1)
  V  <- runif(n, -1, 1)
  X  <- cbind(X1, X2)

  eps <- rnorm(n)
  Y   <- as.numeric(X %*% beta + m_fun(V) + eps)

  lin_pred <- alpha[1] + alpha[2] * Y +
              alpha[3] * X1 + alpha[4] * X2 +
              alpha[5] * V

  pi <- plogis(lin_pred)
  R  <- rbinom(n, 1, pi)

  list(Y = Y, X = X, V = V, R = R, pi = pi)
}

#############################################
# 2. Selection Model Estimation
#############################################
estimate_pi <- function(Y, X, V, R) {
  dat <- data.frame(R = R, Y = Y,
                    X1 = X[,1], X2 = X[,2], V = V)
  fit <- glm(R ~ Y + X1 + X2 + V, family = binomial, data = dat)
  fitted(fit)
}

#############################################
# 3. Kernel Estimation of m(V)
#############################################
estimate_m <- function(V, Y, R, h) {
  idx <- which(R == 1)
  fit <- locpoly(V[idx], Y[idx], bandwidth = h, degree = 1)
  approx(fit$x, fit$y, xout = V, rule = 2)$y
}

#############################################
# 4. Estimators
#############################################
estimate_beta_CC <- function(Y, X, R) {
  idx <- which(R == 1)
  solve(t(X[idx,]) %*% X[idx,], t(X[idx,]) %*% Y[idx])
}

estimate_beta_IPW <- function(Y, X, R, pi_hat) {
  W <- R / pi_hat
  solve(t(X * W) %*% X, t(X * W) %*% Y)
}

estimate_beta_IMP <- function(X, m_hat) {
  solve(t(X) %*% X, t(X) %*% m_hat)
}

estimate_beta_AIPW <- function(Y, X, V, R, pi_hat, m_hat) {
  W <- 1 / pi_hat
  Y_star <- R * Y + (1 - R) * m_hat
  solve(t(X * W) %*% X, t(X * W) %*% Y_star)
}

#############################################
# 5. Sandwich Variance Estimator
#############################################
sandwich_var <- function(Y, X, V, R, pi_hat, m_hat, beta_hat) {
  eps_hat <- Y - X %*% beta_hat - m_hat
  W <- 1 / pi_hat
  psi <- (R * eps_hat + (1 - R) * eps_hat) * W * X
  A <- colMeans(W * X * X)
  B <- cov(psi)
  solve(A) %*% B %*% solve(A) / length(Y)
}

#############################################
# 6. Monte Carlo Simulation
#############################################
run_simulation <- function(n = 500, B = 1000) {
  beta0 <- c(1, -1)
  m_fun <- function(v) sin(pi * v)
  alpha <- c(0, -1, 0.5, -0.5, 0.5)

  res_AIPW <- matrix(NA, B, 2)

  for (b in 1:B) {
    dat <- generate_data(n, beta0, m_fun, alpha)
    pi_hat <- estimate_pi(dat$Y, dat$X, dat$V, dat$R)
    h <- bw.ucv(dat$V[dat$R == 1])
    m_hat <- estimate_m(dat$V, dat$Y, dat$R, h)
    beta_hat <- estimate_beta_AIPW(dat$Y, dat$X, dat$V,
                                   dat$R, pi_hat, m_hat)
    res_AIPW[b,] <- beta_hat
  }

  list(
    mean = colMeans(res_AIPW),
    bias = colMeans(res_AIPW) - beta0,
    sd   = apply(res_AIPW, 2, sd),
    rmse = sqrt(colMeans((res_AIPW - beta0)^2))
  )
}

#############################################
# 7. Example Run (Uncomment to execute)
#############################################
# results <- run_simulation(n = 500, B = 1000)
# print(results)
