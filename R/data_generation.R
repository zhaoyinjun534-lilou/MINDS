# Internal helper: generate true structural parameters (deterministic given seed 4).
.par_fun <- function(Nc, Nd1, Nd2, Nt) {
  set.seed(4)

  # Item-specific fixed effects (intercepts)
  a_1.true <- rnorm(Nd1, mean = 0, sd = 0.5) * (-1)
  a_2.true <- rnorm(Nd2, mean = 0, sd = 0.5) * (-1)

  # Cluster centers in latent space
  x.true <- replicate(Nt, runif(Nc, 0, 2))
  x.true <- apply(x.true, 2, scale) * 2

  # Modality-specific linear maps from latent space to item space
  V.raw.item <- replicate(Nt, runif(Nd1, 0, 2))
  U.raw.item <- replicate(Nt, runif(Nd2, 0, 2))
  V.raw <- t(V.raw.item)
  U.raw <- t(U.raw.item)

  # 1. QR decomposition with positive diagonal on R
  L.raw  <- rbind(t(V.raw), t(U.raw))
  qr.L   <- qr(L.raw)
  Q      <- qr.Q(qr.L)[, seq_len(Nt), drop = FALSE]
  R      <- qr.R(qr.L)[seq_len(Nt), seq_len(Nt), drop = FALSE]
  sign.d <- sign(diag(R))
  Q      <- sweep(Q, 2, sign.d, "*")
  R      <- diag(sign.d) %*% R
  x.true <- x.true %*% t(R)
  V.true <- t(Q[seq_len(Nd1),       , drop = FALSE])
  U.true <- t(Q[Nd1 + seq_len(Nd2), , drop = FALSE])

  # Fix dimension identifiability by ordering latent dimensions using L2 norms.
  x_l2      <- apply(x.true, 2, function(col_j) sqrt(sum(col_j^2)))
  dim_order <- order(x_l2, decreasing = TRUE)
  x.true    <- x.true[, dim_order, drop = FALSE]
  V.true    <- V.true[dim_order, , drop = FALSE]
  U.true    <- U.true[dim_order, , drop = FALSE]
  R         <- R[dim_order, dim_order, drop = FALSE]

  # 2. Location constraint: sum_k omega_k X_k = 0
  x.shift  <- apply(x.true, 2, mean)
  x.true   <- sweep(x.true, 2, x.shift, FUN = "-")
  a_1.true <- a_1.true + as.numeric(x.shift %*% V.true)
  a_2.true <- a_2.true + as.numeric(x.shift %*% U.true)

  list(x.true = x.true, V.true = V.true, U.true = U.true,
       a_1.true = a_1.true, a_2.true = a_2.true, R = R)
}

# Internal helper: generate latent class assignments and random effects.
.par_fun2 <- function(Nb, Nt, theta, sigma2_b, R) {
  set.seed(4)
  Z.true <- rmultinom(Nb, 1, theta) # Nc x Nb
  b.true <- MASS::mvrnorm(Nb, rep(0, Nt), Sigma = diag(sigma2_b))
  b.true <- b.true %*% t(R)
  list(Z.true = Z.true, b.true = b.true)
}

# Internal helper: generate observed binary and continuous outcomes.
.y_fun <- function(seed.no, Nb, Nd2, Z.true, x.true, b.true,
                   V.true, U.true, a_1.true, a_2.true, sigma2_y_2.true) {
  set.seed(seed.no)

  temp  <- (t(Z.true) %*% x.true + b.true) %*% V.true  # Nb x Nd1
  temp1 <- t(apply(temp, 1, function(t) t + a_1.true))
  p     <- exp(temp1) / (1 + exp(temp1))
  y_1   <- apply(p, c(1, 2), function(a) rbinom(1, 1, a))

  temp  <- (t(Z.true) %*% x.true + b.true) %*% U.true  # Nb x Nd2
  temp1 <- t(apply(temp, 1, function(t) t + a_2.true))
  y_2   <- t(sapply(seq_len(Nb),
                    function(i) MASS::mvrnorm(1, temp1[i, ], diag(sigma2_y_2.true))))

  list(y_1 = y_1, y_2 = y_2)
}

#' Generate Simulated Mixed-Modality Dataset
#'
#' Generates a simulated dataset for the MINDS model, including binary outcomes
#' (`y_1`) and continuous outcomes (`y_2`).
#'
#' @param seed.no Integer random seed for data generation. Default `103`.
#' @param Nb Number of subjects. Default `1000`.
#' @param Nc Number of latent clusters. Default `4`.
#' @param Nd1 Number of binary items. Default `10`.
#' @param Nd2 Number of continuous items. Default `12`.
#' @param Nt Number of latent constructs. Default `3`.
#' @param sigma2_b_element Scalar variance of the random effects. Default `0.1`.
#' @param theta Numeric vector of length `Nc` giving cluster mixture weights.
#'   If `NULL` (default), uses `c(0.3, 0.3, 0.2, 0.2)` when `Nc == 4`,
#'   otherwise equal weights `rep(1/Nc, Nc)`.
#'
#' @return A named list with:
#' \itemize{
#'   \item `y_1`: Binary outcome matrix (subjects x binary items).
#'   \item `y_2`: Continuous outcome matrix (subjects x continuous items).
#' }
#'
#' @examples
#' sim_data <- generate_data_mixed(seed.no = 103)
#' str(sim_data)
#'
#' @export
generate_data_mixed <- function(seed.no          = 103,
                                Nb               = 1000,
                                Nc               = 4,
                                Nd1              = 10,
                                Nd2              = 12,
                                Nt               = 3,
                                sigma2_b_element = 0.1,
                                theta            = NULL) {
  if (is.null(theta)) {
    theta <- if (Nc == 4) c(0.3, 0.3, 0.2, 0.2) else rep(1 / Nc, Nc)
  }
  sigma2_b.true   <- rep(sigma2_b_element, Nt)
  sigma2_y_2.true <- rep(1, Nd2)

  raw.par  <- .par_fun(Nc, Nd1, Nd2, Nt)
  raw.par2 <- .par_fun2(Nb, Nt, theta, sigma2_b.true, raw.par$R)

  y.sim <- .y_fun(seed.no, Nb, Nd2,
                  raw.par2$Z.true, raw.par$x.true, raw.par2$b.true,
                  raw.par$V.true,  raw.par$U.true,
                  raw.par$a_1.true, raw.par$a_2.true,
                  sigma2_y_2.true)

  list(y_1 = y.sim$y_1, y_2 = y.sim$y_2)
}