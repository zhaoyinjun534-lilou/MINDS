#' Generate Simulated Mixed-Modality Dataset
#'
#' Generates a simulated dataset with binary and continuous outcomes
#' based on the MINDS data generation mechanism.
#'
#' @param seed.no Integer random seed for data generation.
#'
#' @return A named list with:
#' \itemize{
#'   \item `y_1`: Binary outcome matrix (subjects × binary items).
#'   \item `y_2`: Continuous outcome matrix (subjects × continuous items).
#' }
#'
#' @examples
#' \dontrun{
#' sim_data <- generate_data_mixed(103)
#' str(sim_data)
#' }
#'
#' @export
generate_data_mixed <- function(seed.no) {
  # Parameters for data generation (same as data-raw/make_data_mixed.R)
  Nd1 <- 10
  Nd2 <- 12
  Nb <- 1000
  Nc <- 4
  Nt <- 3
  
  sigma2_b.true <- rep(0.1, Nt)
  theta.true <- c(0.3, 0.3, 0.2, 0.2)
  sigma2_y_2.true <- rep(1, Nd2)
  
  # Generate true parameter values (internal version)
  set.seed(4)
  a_1.true <- rnorm(Nd1, mean = 0, sd = 0.5) * (-1)
  a_2.true <- rnorm(Nd2, mean = 0, sd = 0.5) * (-1)
  
  x.true <- replicate(Nt, runif(Nc, 0, 2))
  x.true <- apply(x.true, 2, scale) * 2
  
  V.raw.item <- replicate(Nt, runif(Nd1, 0, 2))
  U.raw.item <- replicate(Nt, runif(Nd2, 0, 2))
  
  V.true <- t(V.raw.item)
  U.true <- t(U.raw.item)
  
  # Generate latent cluster assignments and random effects
  set.seed(4)
  Z <- rmultinom(Nb, 1, theta.true)
  b <- MASS::mvrnorm(Nb, rep(0, Nt), Sigma = diag(sigma2_b.true))
  
  Z.true <- Z
  b.true <- b
  
  # Generate binary and continuous outcomes
  set.seed(seed.no)
  temp <- (t(Z.true) %*% x.true + b.true) %*% V.true
  temp1 <- t(apply(temp, 1, function(t) {t + a_1.true}))
  p <- exp(temp1) / (1 + exp(temp1))
  
  y_1 <- apply(p, c(1, 2), function(a) rbinom(1, 1, a))
  
  temp <- (t(Z.true) %*% x.true + b.true) %*% U.true
  temp1 <- t(apply(temp, 1, function(t) {t + a_2.true}))
  
  y_2 <- t(sapply(seq(Nb), 
                   function(i) MASS::mvrnorm(1, temp1[i, ], diag(sigma2_y_2.true))))
  
  list(y_1 = y_1, y_2 = y_2)
}
