#' Draw True MINDS Parameters For A Simulation Setting
#'
#' Draws the true cluster centres, item loadings and item intercepts used by
#' [generate_data_mixed()]. Item-level parameters are drawn once at a fixed
#' "item bank" size (`Nd1_max` / `Nd2_max`) and then truncated, so settings that
#' differ only in the number of items share nested loadings.
#'
#' @param param.seed Integer seed for the parameter draws. Five consecutive
#'   seeds (`param.seed` .. `param.seed + 4`) are used, one per parameter block.
#' @param Nd1 Integer number of binary items to keep.
#' @param Nd2 Integer number of continuous items to keep.
#' @param Nc Integer number of latent clusters.
#' @param Nt Integer number of latent constructs.
#' @param Nd1_max,Nd2_max Integer size of the binary / continuous item bank the
#'   parameters are drawn at before truncation. Set these at least as large as
#'   the biggest item count in a sweep to keep loadings nested across settings.
#'   Default `max(Nd1, 100)` / `max(Nd2, 100)`.
#'
#' @return A named list with `a_1.true`, `a_2.true`, `x.true`, `V.true` and
#'   `U.true`.
#'
#' @examples
#' pars <- generate_true_params(param.seed = 4, Nd1 = 10, Nd2 = 12, Nc = 5, Nt = 3)
#' str(pars)
#'
#' @export
generate_true_params <- function(param.seed, Nd1, Nd2, Nc, Nt,
                                 Nd1_max = max(Nd1, 100L),
                                 Nd2_max = max(Nd2, 100L)) {
  # NOTE ON NESTING: item-level parameters (intercepts a_*, loadings V/U) are
  # drawn once at the maximum item count (Nd1_max / Nd2_max) and then the first
  # Nd1 / Nd2 items are kept. This makes settings with different item counts
  # NESTED -- the loadings for, say, 20 items are exactly the first 20 rows of
  # the loadings for any larger setting -- so varying the number of items no
  # longer changes which random loadings are drawn. To guarantee nesting across
  # an entire sweep, set Nd1_max / Nd2_max >= the largest item count you use
  # (the default of 100 covers any sweep with <= 100 items per modality).
  set.seed(param.seed);     x.true      <- replicate(Nt, runif(Nc,  0, 2))
  x.true <- apply(x.true, 2, scale) * 2
  set.seed(param.seed + 1); a_1.full    <- rnorm(Nd1_max, mean = 0, sd = 0.5)
  set.seed(param.seed + 2); a_2.full    <- rnorm(Nd2_max, mean = 0, sd = 0.5)
  set.seed(param.seed + 3); V.full.item <- replicate(Nt, runif(Nd1_max, 0, 2))
  set.seed(param.seed + 4); U.full.item <- replicate(Nt, runif(Nd2_max, 0, 2))

  # Subset to the requested number of items (nested prefix).
  a_1.true   <- a_1.full[seq_len(Nd1)]
  a_2.true   <- a_2.full[seq_len(Nd2)]
  V.raw.item <- V.full.item[seq_len(Nd1), , drop = FALSE]
  U.raw.item <- U.full.item[seq_len(Nd2), , drop = FALSE]

  V.true <- t(V.raw.item)
  U.true <- t(U.raw.item)



  list(a_1.true = a_1.true, a_2.true = a_2.true,
       x.true = x.true, V.true = V.true, U.true = U.true)
}

#' Generate Simulated Mixed-Modality Dataset
#'
#' Generates a simulated dataset with binary and continuous outcomes based on
#' the MINDS data generation mechanism, together with the true parameter values
#' that produced it (so a fit can be scored against the truth).
#'
#' @param seed.no Integer random seed for the subject-level draws: cluster
#'   labels `Z`, random effects `b`, and the outcomes themselves.
#' @param param.seed Integer seed for the true parameters, passed to
#'   [generate_true_params()]. Holding `param.seed` fixed while varying
#'   `seed.no` gives independent replicates of the same setting.
#' @param Nd1 Integer number of binary items. Default 10.
#' @param Nd2 Integer number of continuous items. Default 12.
#' @param Nb Integer number of subjects. Default 1000.
#' @param Nc Integer number of latent clusters. Default 5.
#' @param Nt Integer number of latent constructs. Default 3.
#' @param Nd1_max,Nd2_max Item-bank sizes for nested loadings, passed to
#'   [generate_true_params()].
#' @param sigma2_b Variance of the subject random effects (per latent
#'   construct); sets `sigma2_b.true <- rep(sigma2_b, Nt)`. Default 0.1.
#' @param theta Numeric vector of true cluster weights; must have length `Nc`.
#'   Default `c(0.3, 0.15, 0.15, 0.2, 0.2)`, the five-cluster simulation setting.
#'
#' @return A named list with:
#' \itemize{
#'   \item `y_1`: Binary outcome matrix (subjects x binary items).
#'   \item `y_2`: Continuous outcome matrix (subjects x continuous items).
#'   \item `x.true`: True cluster centre matrix (clusters x latent constructs).
#'   \item `V.true`: True loading matrix for binary items (constructs x items).
#'   \item `U.true`: True loading matrix for continuous items (constructs x items).
#'   \item `a_1.true`: True intercepts for binary items.
#'   \item `a_2.true`: True intercepts for continuous items.
#'   \item `Z.true`: True cluster assignment matrix (clusters x subjects).
#'   \item `b.true`: True subject random effects (subjects x latent constructs).
#'   \item `theta.true`: True cluster weight vector.
#'   \item `sigma2_y_2.true`: True residual variances for continuous items.
#'   \item `sigma2_b.true`: True latent variances (per construct).
#'   \item `Nc`, `Nt`, `Nb`: The setting's dimensions, echoed back.
#' }
#'
#' @examples
#' \dontrun{
#' sim_data <- generate_data_mixed(103)
#' str(sim_data)
#' }
#'
#' @seealso [generate_true_params()], [MINDS_algorithm()]
#'
#' @export
generate_data_mixed <- function(seed.no, param.seed = 4,
                                Nd1 = 10, Nd2 = 12, Nb = 1000, Nc = 5, Nt = 3,
                                Nd1_max = max(Nd1, 100L),
                                Nd2_max = max(Nd2, 100L),
                                sigma2_b = 0.1,
                                theta = c(0.3, 0.15, 0.15, 0.2, 0.2)) {

  # The default weights are the five-cluster setting used throughout the
  # simulations. Nc is a user argument, so a mismatch would otherwise surface
  # far downstream as a dimension error between Z (length(theta) x Nb) and
  # x.true (Nc x Nt).
  if (length(theta) != Nc) {
    stop("`theta` must have length Nc (got ", length(theta), " for Nc = ", Nc,
         "). Supply cluster weights matching Nc.")
  }

  sigma2_b.true <- rep(sigma2_b, Nt)
  theta.true <- theta
  sigma2_y_2.true <- rep(1, Nd2)

  params <- generate_true_params(param.seed, Nd1, Nd2, Nc, Nt,
                                 Nd1_max = Nd1_max, Nd2_max = Nd2_max)
  a_1.true <- params$a_1.true
  a_2.true <- params$a_2.true
  x.true   <- params$x.true
  V.true   <- params$V.true
  U.true   <- params$U.true

  # Generate latent cluster assignments and random effects.
  # Seeded by seed.no so replicates have independent Z and b (seed.no = 4
  # reproduces the previous fixed draw, so existing single-fit scripts are
  # unchanged).
  set.seed(seed.no)
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

  list(y_1 = y_1, y_2 = y_2,
       x.true = x.true, V.true = V.true, U.true = U.true,
       a_1.true = a_1.true, a_2.true = a_2.true, Z.true = Z.true,
       b.true = b.true, theta.true = theta.true,
       sigma2_y_2.true = sigma2_y_2.true, sigma2_b.true = sigma2_b.true,
       Nc = Nc, Nt = Nt, Nb = Nb)
}
