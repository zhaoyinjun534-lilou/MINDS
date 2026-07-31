#' Fit MINDS Algorithm For Mixed Binary-Continuous Data
#'
#' Runs the Bayesian MINDS model using binary outcomes `y_1` and continuous
#' outcomes `y_2`, with user-specified prior hyperparameters.
#'
#' @param y_1 Binary outcome matrix with subjects in rows and items in columns.
#'   Must be a matrix (not a data.frame). `NA` entries are allowed: they are
#'   masked out of every conditional update via `M_1` and excluded from the
#'   likelihood, i.e. an observed-data (MAR) analysis with no imputation in the
#'   fit itself. Items with no observed variation are dropped from the KAMILA
#'   warm start only.
#' @param y_2 Continuous outcome matrix with subjects in rows and items in
#'   columns. Same matrix / `NA` handling as `y_1`, via `M_2`.
#' @param Nc Integer number of latent clusters.
#' @param Nt Integer number of latent constructs.
#' @param iter.max Integer number of MCMC iterations.
#' @param mu_x Prior mean for cluster centers `x`.
#' @param sigma2_x Prior variance for cluster centers `x`.
#' @param IG_b.shape Inverse-gamma shape for random-effect variance `sigma2_b`.
#' @param IG_b.scale Inverse-gamma scale for random-effect variance `sigma2_b`.
#' @param p.theta.prior Dirichlet prior parameters for cluster weights.
#' @param IG_y_2.shape Inverse-gamma shape for `sigma2_y_2`.
#' @param IG_y_2.scale Inverse-gamma scale for `sigma2_y_2`.
#' @param sigma2_a_1 Prior variance for binary intercepts.
#' @param sigma2_a_2 Prior variance for continuous intercepts.
#' @param sigma2_v Prior variance for binary loadings.
#' @param mu_v Prior mean for binary loadings.
#' @param sigma2_u Prior variance for continuous loadings.
#' @param mu_u Prior mean for continuous loadings.
#' @param true.values Optional list of true parameter values used for initialization and equation solving.
#'   Expected list components are `x`, `V`, `U`, `a_1`, `a_2`, `Z`, and `theta`.
#' @param ref_cluster Integer in 1:Nc; cluster fixed at the latent-space origin
#'   (location constraint, passed to `canonicalize_minds()`). For a truth-vs-
#'   estimate comparison the reported estimates are re-centred so this cluster,
#'   in the truth-aligned labelling, sits at the origin.
#' @param var_b Latent variance imposed on each axis by the identifiability
#'   constraint (passed to `canonicalize_minds()`). Default `1`; set to the
#'   data-generating scale (e.g. `0.1`) to report loadings on that metric. Use
#'   the same value for the fit and for any truth-vs-estimate comparison.
#' @param init_seed Seed used for initialization and imputation.
#' @param plot_trace Logical; if `TRUE`, draw the likelihood trace.
#' @param init_params Optional list of starting values overriding the defaults.
#'   Recognised components are `x`, `V`, `U`, `a_1`, `a_2`, `sigma2_y_2`,
#'   `sigma2_b`, `b`, `Z` and `theta`. Supplying `Z` skips the `mice` + KAMILA
#'   warm start entirely. Combined with the `update_*` flags below this gives a
#'   frozen-parameter decode: parameters whose flag is `FALSE` stay fixed at the
#'   values given here for the whole run.
#' @param update_mu Logical; if `FALSE`, hold the cluster centres `x` fixed.
#' @param update_gamma_U Logical; if `FALSE`, hold the loadings `V`, `U` and the
#'   intercepts `a_1`, `a_2` fixed.
#' @param update_sigma Logical; if `FALSE`, hold `sigma2_b` and `sigma2_y_2`
#'   fixed.
#' @param update_b Logical; if `FALSE`, hold the subject random effects `b`
#'   fixed.
#'
#' @return A named list with the posterior summaries and model-selection
#'   criteria:
#' \itemize{
#'   \item `membership`: integer vector of estimated cluster labels (length Nb).
#'   \item `cluster center`: Nc x Nt matrix of canonical cluster centres.
#'   \item `loading to binary modality`, `loading to continuous modality`:
#'     Nt x Nd1 and Nt x Nd2 canonical loading matrices.
#'   \item `binary modality intercept`, `continuous modality intercept`: item
#'     intercept vectors.
#'   \item `membership weight`: estimated cluster weights.
#'   \item `b`: Nb x Nt subject random effects.
#'   \item `sigma2_y_2`, `sigma2_b`: residual and latent variances.
#'   \item `axis order`, `axis sign`: modal latent-axis permutation and signs
#'     across the retained draws.
#'   \item `credible intervals`: 95% intervals for the centres, loadings and
#'     intercepts.
#'   \item `likelihood trace plot`: the per-iteration log-likelihood trace.
#'   \item `dic_complete`: the primary model-selection criterion, equal to
#'     `dic4_complete`.
#'   \item `dic4_complete`, `dic5_complete`, `dic6_complete` and the matching
#'     `pD4_complete`, `pD5_complete`, `pD6_complete`, `Dbar_complete`: the
#'     complete-data DIC family of Celeux et al. (2006). Smaller is better.
#'     DIC4 is the recommended criterion; DIC5 and DIC6 are diagnostics.
#' }
#'
#' @references
#' Celeux, G., Forbes, F., Robert, C. P. and Titterington, D. M. (2006).
#' Deviance Information Criteria for Missing Data Models.
#' \emph{Bayesian Analysis}, 1(4), 651-706.
#'
#' @seealso [generate_data_mixed()], [canonicalize_minds()]
#'
#' @export
MINDS_algorithm <- function(
    y_1,
    y_2,
    Nc,
    Nt,
    iter.max = 800,
    mu_x = 0,
    sigma2_x = 1,
    IG_b.shape = 11,
    IG_b.scale = 1,
    p.theta.prior = rep(1, Nc),
    IG_y_2.shape = 2,
    IG_y_2.scale = 1,
    sigma2_a_1 = 10,
    sigma2_a_2 = 10,
    sigma2_v = 10,
    mu_v = 0,
    sigma2_u = 10,
    mu_u = 0,
    true.values = NULL,
    var_b = 1,
    ref_cluster = 1L,
    init_seed = 123,
    plot_trace = TRUE,
    init_params = NULL,
    update_mu = TRUE,
    update_gamma_U = TRUE,
    update_sigma = TRUE,
    update_b = TRUE) {
  Nd1 <- ncol(y_1)
  Nd2 <- ncol(y_2)
  Nb <- nrow(y_1)

  if (nrow(y_2) != Nb) {
    stop("y_1 and y_2 must have the same number of rows.")
  }

  if (!is.null(true.values)) {
    x.true       <- true.values$x
    V.true       <- true.values$V
    U.true       <- true.values$U
    a_1.true     <- true.values$a_1
    a_2.true     <- true.values$a_2
    Z.true       <- true.values$Z
    theta.true   <- true.values$theta
    b.true       <- true.values$b
    sigma2_b.true <- true.values$sigma2_b

    canon.true <- canonicalize_minds(
      x = x.true, b = b.true, V = V.true, U = U.true,
      a_1 = a_1.true, a_2 = a_2.true,
      sigma2_b = sigma2_b.true,
      var_b = var_b, ref_cluster = ref_cluster
    )
    x.true.c   <- canon.true$x
    V.true.c   <- canon.true$V
    U.true.c   <- canon.true$U
    a_1.true.c <- canon.true$a_1
    a_2.true.c <- canon.true$a_2
    b.true.c   <- canon.true$b
  }

  x.ini <- MASS::mvrnorm(Nc, rep(mu_x, Nt), Sigma = diag(rep(sigma2_x, Nt)))
  sigma2_b.ini <- EDISON::rinvgamma(Nt, IG_b.shape, IG_b.scale)
  sigma2_y_2.ini <- EDISON::rinvgamma(Nd2, IG_y_2.shape, IG_y_2.scale)
  b.ini <- MASS::mvrnorm(Nb, rep(0, Nt), Sigma = diag(sigma2_b.ini))

  # An item with no observed values gives mean(NA, na.rm = TRUE) = NaN. That
  # NaN reaches a_1.ini and then the first BayesLogit::rpg() call, which fails
  # on a non-finite z. Fall back to the neutral value for such items (0.5 for a
  # rate, 0 for a location); the Gibbs updates sample them from their priors.
  y_1_mean <- apply(y_1, 2, function(a) mean(a, na.rm = TRUE))
  y_1_mean[!is.finite(y_1_mean)] <- 0.5
  y_1_mean <- pmin(pmax(y_1_mean, 1e-6), 1 - 1e-6)
  # Log-odds, log(p / (1 - p)) -- the natural scale of the binary intercepts,
  # since the model's linear predictor is (x_k + b_i) V + a_1 on the logit
  # scale.
  logit.fun <- function(a) log(a) - log(1 - a)
  a_1.ini <- logit.fun(y_1_mean) - logit.fun(y_1_mean[1])

  y_2_mean <- apply(y_2, 2, function(a) mean(a, na.rm = TRUE))
  y_2_mean[!is.finite(y_2_mean)] <- 0
  a_2.ini <- y_2_mean - y_2_mean[1]

  set.seed(init_seed)
  V.ini <- t(replicate(Nt, runif(Nd1, 0, 1)))
  U.ini <- t(replicate(Nt, runif(Nd2, 0, 1)))
  
 

  # When init_params supplies Z, skip the mice + KAMILA warm start and
  # initialise membership/weights from it (frozen-parameter decode).
  if (is.null(init_params) || is.null(init_params$Z)) {
  # Items with no observed values -- or no observed variation -- carry no
  # information for the warm start, and mice cannot impute an all-NA column.
  # The resulting NA/NaN column reaches KAMILA and errors on the continuous
  # side (scale() -> NaN) or segfaults on the binary side. Drop such items from
  # the warm start ONLY; the Gibbs updates below handle them correctly via the
  # M_1 / M_2 masks, so no data is discarded from the fit itself.
  informative <- function(Y) vapply(seq_len(ncol(Y)), function(j) {
    v <- Y[, j][!is.na(Y[, j])]
    length(v) > 1L && stats::var(v) > 0
  }, logical(1))
  keep_1 <- informative(y_1)
  keep_2 <- informative(y_2)
  if (!any(keep_1) || !any(keep_2)) {
    stop("Warm start needs at least one binary and one continuous item with ",
         "observed variation (found ", sum(keep_1), " and ", sum(keep_2),
         "). Supply init_params$Z to skip the mice + KAMILA warm start.")
  }
  if (any(!keep_1) || any(!keep_2)) {
    message("  (warm start: dropping ", sum(!keep_1), " binary and ",
            sum(!keep_2), " continuous item(s) with no observed variation)")
  }
  y_1.w <- y_1[, keep_1, drop = FALSE]
  y_2.w <- y_2[, keep_2, drop = FALSE]
  Nd1.w <- ncol(y_1.w)

  y <- data.frame(y_1.w == 1, y_2.w)
  imp <- mice::mice(y, m = 5, maxit = 50, method = "pmm", seed = init_seed, printFlag = FALSE)
  y.complete <- mice::complete(imp, 1)
  
  
  
  # Your inputs (edit as needed)
  Y1 <- as.matrix(y.complete[, seq_len(Nd1.w), drop = FALSE])
  Y2 <- as.matrix(y.complete[, -seq_len(Nd1.w), drop = FALSE])
  # Initial clustering for Z via KAMILA (semiparametric mixed continuous +
  # categorical clustering), replacing the iClusterBayes warm start. Continuous
  # items are the Gaussian modality (standardised); binary items are treated as
  # 2-level factors. numClust = Nc (KAMILA clusters directly, no latent factors).
  con_var    <- as.data.frame(scale(Y2))
  cat_factor <- as.data.frame(lapply(seq_len(Nd1.w), function(j) {
    factor(as.integer(Y1[, j]), levels = c(0, 1))
  }))
  names(cat_factor) <- paste0("b", seq_len(Nd1.w))

  km_init   <- kamila::kamila(conVar = con_var, catFactor = cat_factor,
                              numClust = Nc, numInit = 10)
  Z.cat.est <- as.integer(km_init$finalMemb)

  # f5 <- psych::fa(y.cor, nfactors = min(5, ncol(y.cor)))
  # my.score <- psych::factor.scores(y.complete, f5, method = "tenBerge")$scores
  # 
  # cl <- kmeans(my.score, Nc)
  # Z.cat.est <- cl$cluster
  Z.est <- Zcat_to_Zest(Z.cat.est, Nc, Nb)
  Z.ini <- Z.est
  theta.count <- tabulate(Z.cat.est, nbins = Nc)
  theta.ini <- theta.count / sum(theta.count)
  } else {
    Z.ini <- init_params$Z
    theta.ini <- if (!is.null(init_params$theta)) init_params$theta else
      rowSums(Z.ini) / sum(Z.ini)
  }
  alpha <- t(Z.ini)

  # Override initial values with any supplied in init_params. Parameters whose
  # update_* flag is FALSE stay fixed at these values for the whole run
  # (frozen-parameter decoding, analogous to ecm_algorithm's init_params).
  if (!is.null(init_params)) {
    if (!is.null(init_params$x))          x.ini          <- init_params$x
    if (!is.null(init_params$V))          V.ini          <- init_params$V
    if (!is.null(init_params$U))          U.ini          <- init_params$U
    if (!is.null(init_params$a_1))        a_1.ini        <- init_params$a_1
    if (!is.null(init_params$a_2))        a_2.ini        <- init_params$a_2
    if (!is.null(init_params$sigma2_y_2)) sigma2_y_2.ini <- init_params$sigma2_y_2
    if (!is.null(init_params$sigma2_b))   sigma2_b.ini   <- init_params$sigma2_b
    if (!is.null(init_params$b))          b.ini          <- init_params$b
  }


 ## MINDS 1st loop----
  x <- x.ini
  Z <- Z.ini
  b <- b.ini
  sigma2_b <- sigma2_b.ini
  sigma2_y_2 <- sigma2_y_2.ini
  theta <- theta.ini
  a_1 <- a_1.ini
  a_2 <- a_2.ini
  V <- V.ini
  U <- U.ini

  x.all <- NULL
  Z.all <- NULL
  alpha.all <- NULL
  theta.all <- NULL
  a_1.all <- NULL
  a_2.all <- NULL
  V.all <- NULL
  U.all <- NULL
  b.all <- NULL
  sigma2_y_2.all <- NULL
  sigma2_b.all <- NULL
  llk.all <- NULL

  i.iter <- 0

  if (sum(is.na(y_1)) == 0) {
    M_1 <- matrix(1, nrow = Nb, ncol = Nd1)
  } else {
    M_1 <- 1 * !is.na(y_1)
  }

  if (sum(is.na(y_2)) == 0) {
    M_2 <- matrix(1, nrow = Nb, ncol = Nd2)
  } else {
    M_2 <- 1 * !is.na(y_2)
  }

  while (i.iter <= iter.max) {
    i.iter <- i.iter + 1

    temp1 <- (t(Z) %*% x + b) %*% V
    temp <- t(apply(temp1, 1, function(ti) ti + a_1))
    w <- apply(temp, c(1, 2), function(s) BayesLogit::rpg(num = 1, h = 1, z = s))

    if (update_gamma_U) {
      V_a <- 1 / (apply(w * M_1, 2, sum) + 1 / sigma2_a_1)
      temp <- (t(Z) %*% x + b) %*% V
      r <-  (y_1 - 1 / 2) / w - temp
      Mu_a <- apply(w * M_1 * r, 2, function(a) sum(a, na.rm = TRUE)) * V_a
      a_1 <- sapply(seq_len(Nd1), function(i) rnorm(1, Mu_a[i], sqrt(V_a[i])))

      V_a <- 1 / (apply(M_2 * rep(1 / sigma2_y_2, each = nrow(M_2)), 2, sum) + 1 / sigma2_a_2)
      temp <- (t(Z) %*% x + b) %*% U
      r <- y_2 - temp
      Mu_a <- apply(r * M_2, 2, function(a) sum(a, na.rm = TRUE)) / sigma2_y_2 * V_a
      a_2 <- sapply(seq_len(Nd2), function(i) rnorm(1, Mu_a[i], sqrt(V_a[i])))

      for (ti in seq_len(Nt)) {
        index_x <- rep(0, Nt)
        index_x[ti] <- 1

        temp_A <- (t(Z) %*% x + b) %*% index_x
        temp <- apply(w * M_1, 2, function(s) sum(s * temp_A^2))
        B_V_ti <- 1 / (temp + 1 / sigma2_v)

        x.sub <- x
        x.sub[, ti] <- 0
        b.sub <- b
        b.sub[, ti] <- 0

        temp1 <- (t(Z) %*% x.sub + b.sub) %*% V
        temp2 <- t(apply(temp1, 1, function(tt) tt + a_1))
        phi <- (y_1 - 1 / 2) / w - temp2

        temp3 <- apply(w * M_1, 2, function(s) s * temp_A)
        D_V_ti <- apply(temp3 * phi, 2, function(a) sum(a, na.rm = TRUE)) + mu_v / sigma2_v

        Mu_V_ti <- B_V_ti * D_V_ti
        V_ti <- MASS::mvrnorm(1, Mu_V_ti, diag(B_V_ti))
        V[ti, ] <- V_ti
      }
    }  # end update_gamma_U: a_1, a_2, V

    if (update_mu) {
      for (ti in seq_len(Nt)) {
        x.sub <- x
        x.sub[, ti] <- 0
        B_ti <- diag(sigma2_x, Nc)

        temp <- (t(Z) %*% x.sub + b) %*% V
        temp1 <- t(apply(temp, 1, function(tt) tt + a_1))
        phi_1 <- (y_1 - 1 / 2) / w - temp1
        phi_1 <- ifelse(is.na(phi_1), 0, phi_1)

        colnames(w) <- NULL
        colnames(phi_1) <- NULL
        dt <- abind::abind("z" = Z, "w" = t(w), "phi" = t(phi_1), along = 1)

        aa <- lapply(seq_len(Nb), function(i) {
          Z_i <- dt[startsWith(rownames(dt), "z"), i]
          w_i <- dt[startsWith(rownames(dt), "w"), i]
          phi_i <- dt[startsWith(rownames(dt), "phi"), i]

          tempi <- Z_i %*% t(V[ti, ])
          list(
            v_temp = tempi %*% diag(w_i * M_1[i, ]) %*% t(tempi),
            mu_temp = tempi %*% diag(w_i * M_1[i, ]) %*% phi_i
          )
        })

        v_sum_1 <- Reduce("+", lapply(seq_len(ncol(dt)), function(i) aa[[i]]$v_temp))
        mu_sum_1 <- Reduce("+", lapply(seq_len(ncol(dt)), function(i) aa[[i]]$mu_temp))

        temp <- (t(Z) %*% x.sub + b) %*% U
        temp1 <- t(apply(temp, 1, function(tt) tt + a_2))
        Phi_1 <- y_2 - temp1
        Phi_1 <- ifelse(is.na(Phi_1), 0, Phi_1)

        dt <- abind::abind("z" = Z, "Phi" = t(Phi_1), along = 1)
        aa <- lapply(seq_len(Nb), function(i) {
          Z_i <- dt[startsWith(rownames(dt), "z"), i]
          Phi_i <- dt[startsWith(rownames(dt), "Phi"), i]

          tempi <- Z_i %*% t(U[ti, ])
          list(
            v_temp = tempi %*% diag(1 / sigma2_y_2 * M_2[i, ]) %*% t(tempi),
            mu_temp = tempi %*% diag(1 / sigma2_y_2 * M_2[i, ]) %*% Phi_i
          )
        })

        v_sum_2 <- Reduce("+", lapply(seq_len(ncol(dt)), function(i) aa[[i]]$v_temp))
        mu_sum_2 <- Reduce("+", lapply(seq_len(ncol(dt)), function(i) aa[[i]]$mu_temp))

        V_x_ti_inverse <- solve(B_ti) + v_sum_1 + v_sum_2
        V_x_ti <- solve(V_x_ti_inverse)
        Mu_x_ti <- V_x_ti %*% (mu_sum_1 + mu_sum_2 + solve(B_ti) %*% rep(mu_x, Nc))

        x_ti <- MASS::mvrnorm(1, Mu_x_ti, V_x_ti)
        x[, ti] <- x_ti
      }
    }  # end update_mu: x

    if (update_b) {
      for (ti in seq_len(Nt)) {
        b.sub <- b
        b.sub[, ti] <- 0

        temp <- (t(Z) %*% x + b.sub) %*% V
        temp1 <- t(apply(temp, 1, function(tt) tt + a_1))
        C_ti <- (y_1 - 1 / 2) / w - temp1

        temp <- (t(Z) %*% x + b.sub) %*% U
        temp1 <- t(apply(temp, 1, function(tt) tt + a_2))
        D_ti <- y_2 - temp1

        V_b_ti_diag <- 1 / (apply(w * M_1, 1, function(w_i) sum(w_i * V[ti, ]^2)) +
                              apply(M_2, 1, function(M2_i) sum(M2_i * 1 / sigma2_y_2 * U[ti, ]^2)) +
                              1 / sigma2_b[ti])

        Mu_b_ti <- (apply(w * M_1 * C_ti, 1, function(ti2) sum(ti2 * V[ti, ], na.rm = TRUE)) +
                      apply(D_ti * M_2, 1, function(ti2) {
                        sum(ti2 * U[ti, ] * 1 / sigma2_y_2, na.rm = TRUE)
                      })) * V_b_ti_diag

        b_ti <- sapply(seq_len(Nb), function(bi) rnorm(1, Mu_b_ti[bi], sqrt(V_b_ti_diag[bi])))
        b[, ti] <- b_ti
      }
    }  # end update_b: b

    if (update_gamma_U) {
      for (ti in seq_len(Nt)) {
        index_x <- rep(0, Nt)
        index_x[ti] <- 1

        temp_A <- (t(Z) %*% x + b) %*% index_x
        temp_B <- apply(M_2, 2, function(row) row %*% temp_A^2)
        B_U_ti <- 1 / (temp_B / sigma2_y_2 + 1 / sigma2_u)

        x.sub <- x
        x.sub[, ti] <- 0
        b.sub <- b
        b.sub[, ti] <- 0

        temp1 <- (t(Z) %*% x.sub + b.sub) %*% U
        temp2 <- t(apply(temp1, 1, function(tt) tt + a_2))
        phi <- y_2 - temp2

        temp3 <- apply(phi * M_2, 2, function(s) sum(s * temp_A, na.rm = TRUE))
        D_U_ti <- temp3 / sigma2_y_2 + mu_u / sigma2_u

        Mu_U_ti <- B_U_ti * D_U_ti
        U_ti <- MASS::mvrnorm(1, Mu_U_ti, diag(B_U_ti))
        U[ti, ] <- U_ti
      }
    }  # end update_gamma_U: U


    if (i.iter >= iter.max * 0.2) {
      alpha <- NULL
      for (k in seq_len(Nc)) {
        temp1 <- t(apply(b, 1, function(b_i) x[k, ] + b_i)) %*% V
        temp2 <- t(apply(temp1, 1, function(tt) tt + a_1)) - (y_1 - 1 / 2) / w
        temp2 <- ifelse(is.na(temp2), 0, temp2)
        p1_k <- exp(-w * M_1 / 2 * temp2^2)

        temp1 <- t(apply(b, 1, function(b_i) x[k, ] + b_i)) %*% U
        temp2 <- t(apply(temp1, 1, function(tt) tt + a_2))

        sigma2_y_2_matrix <- matrix(rep(sigma2_y_2, Nb), nrow = Nb, byrow = TRUE)
        aa <- temp2 - y_2
        aa <- ifelse(is.na(aa), 0, aa)
        p2_k <- (1 / sqrt(2 * pi * sigma2_y_2_matrix))^M_2 * exp(-(aa)^2 / 2 / sigma2_y_2_matrix)

        alpha.k <- apply(p1_k, 1, prod) * apply(p2_k, 1, prod) * theta[k]
        alpha <- cbind(alpha, alpha.k)
      }

      alpha <- alpha / apply(alpha, 1, sum)
      alpha[is.nan(alpha) | is.na(alpha)] <- 1 / Nc
      Z <- apply(alpha, 1, function(p) rmultinom(1, 1, p))
    }

      theta <- LaplacesDemon::rdirichlet(1, apply(Z, 1, sum) + p.theta.prior)
      if (update_sigma) {
        sigma2_b <- apply(b, 2, function(b_ti) EDISON::rinvgamma(1, Nb / 2 + IG_b.shape, sum(b_ti^2) / 2 + IG_b.scale))

        temp <- (t(Z) %*% x + b) %*% U
        temp1 <- t(apply(temp, 1, function(tt) tt + a_2))
        C <- y_2 - temp1
        aa <- apply(C^2 * M_2, 2, function(a) sum(a, na.rm = TRUE)) / 2
        # Shape must count the OBSERVED values of item j, not Nb: aa[j] only
        # accumulates residuals where M_2 == 1, so using Nb over-counts the
        # sample size and shrinks sigma2_y_2 by roughly (1 - missing rate).
        n_obs_2 <- colSums(M_2)
        sigma2_y_2 <- unlist(lapply(seq_len(Nd2), function(j) {
          EDISON::rinvgamma(1, n_obs_2[j] / 2 + IG_y_2.shape, aa[j] + IG_y_2.scale)
        }))
      }

    temp <- (t(Z) %*% x + b) %*% V
    temp1 <- t(apply(temp, 1, function(tt) tt + a_1))
    p <- exp(temp1) / (1 + exp(temp1))
    # Clamp away from {0,1} so a saturated p never makes llk = -Inf; keeps the
    # stored per-draw llk (llk.all) finite for the chain-based D_bar estimate.
    p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
    llk_binary <- sum(apply(M_1 * y_1 * log(p) + log(1 - p) * M_1 * (1 - y_1), 1,
                            function(a) sum(a, na.rm = TRUE)), na.rm = TRUE)

    temp <- (t(Z) %*% x + b) %*% U
    temp1 <- t(apply(temp, 1, function(tt) tt + a_2))
    logd <- sapply(seq_len(Nd2), function(j) {
      M_2[, j] * dnorm(y_2[, j], mean = temp1[, j], sd = sqrt(sigma2_y_2[j]), log = TRUE)
    })
    llk_continuous <- sum(apply(logd, 2, function(a) sum(a, na.rm = TRUE)), na.rm = TRUE)
    llk <- llk_binary + llk_continuous

    x.all <- abind::abind(x.all, x, along = 3)
    Z.all <- abind::abind(Z.all, Z, along = 3)
    V.all <- abind::abind(V.all, V, along = 3)
    U.all <- abind::abind(U.all, U, along = 3)
    b.all <- abind::abind(b.all, b, along = 3)

    theta.all <- rbind(theta.all, theta)
    sigma2_y_2.all <- cbind(sigma2_y_2.all, sigma2_y_2)
    a_2.all <- cbind(a_2.all, a_2)
    a_1.all <- cbind(a_1.all, a_1)
    alpha.all <- abind::abind(alpha.all, alpha, along = 3)
    llk.all <- cbind(llk.all, llk)
    sigma2_b.all <- cbind(sigma2_b.all, sigma2_b)
  }

  iter.all <- dim(x.all)[3]
  n.burn <- min(round(iter.all * 0.8, 0), iter.all - 1)
  post.idx <- seq.int(n.burn + 1, iter.all)

  # --- Identify every retained draw, then summarise --------------------------
  # canonicalize_minds() maps each posterior draw to a common identified
  # representative by removing the model's only continuous indeterminacies
  # (per-axis scale, per-axis sign, axis permutation, and latent-space location).
  # It is a closed-form function of the draw -- no optimisation and no use of the
  # true values -- so both the point estimates (median) and the 95% credible
  # intervals below are well defined and independent of the initialisation. 
  nS <- length(post.idx)
  x.canon.all   <- array(NA_real_, c(Nc,  Nt,  nS))
  V.canon.all   <- array(NA_real_, c(Nt,  Nd1, nS))
  U.canon.all   <- array(NA_real_, c(Nt,  Nd2, nS))
  b.canon.all   <- array(NA_real_, c(Nb,  Nt,  nS))
  a_1.canon.all <- matrix(NA_real_, Nd1, nS)
  a_2.canon.all <- matrix(NA_real_, Nd2, nS)
  ord.all       <- matrix(NA_integer_, Nt, nS)   # axis permutation per draw
  sgn.all       <- matrix(NA_integer_, Nt, nS)   # axis sign per draw

  for (s in seq_len(nS)) {
    it <- post.idx[s]
    cs <- canonicalize_minds(
      x = x.all[, , it], b = b.all[, , it],
      V = V.all[, , it], U = U.all[, , it],
      a_1 = a_1.all[, it], a_2 = a_2.all[, it],
      sigma2_b = sigma2_b.all[, it],
      var_b = var_b, ref_cluster = ref_cluster
    )
    x.canon.all[, , s] <- cs$x
    V.canon.all[, , s] <- cs$V
    U.canon.all[, , s] <- cs$U
    b.canon.all[, , s] <- cs$b
    a_1.canon.all[, s] <- cs$a_1
    a_2.canon.all[, s] <- cs$a_2
    ord.all[, s]       <- cs$order
    sgn.all[, s]       <- cs$sign
  }

  # representative (modal) axis order/sign across draws; the mode is a faithful
  # summary of the energy-based canonicalisation applied to each draw.
  mode_int   <- function(v) as.integer(names(which.max(table(v))))
  axis.order <- apply(ord.all, 1, mode_int)
  axis.sign  <- apply(sgn.all, 1, mode_int)

  med2 <- function(A) apply(A, c(1, 2), median)
  ci   <- function(A, m) list(lower = apply(A, m, stats::quantile, 0.025),
                              upper = apply(A, m, stats::quantile, 0.975))

  x.est   <- med2(x.canon.all)
  V.est   <- med2(V.canon.all)
  U.est   <- med2(U.canon.all)
  b.est   <- med2(b.canon.all)
  a_1.est <- apply(a_1.canon.all, 1, median)
  a_2.est <- apply(a_2.canon.all, 1, median)
  b       <- b.est

  x.ci   <- ci(x.canon.all,   c(1, 2))
  V.ci   <- ci(V.canon.all,   c(1, 2))
  U.ci   <- ci(U.canon.all,   c(1, 2))
  a_1.ci <- ci(a_1.canon.all, 1)
  a_2.ci <- ci(a_2.canon.all, 1)

  sigma2_y_2.est <- apply(sigma2_y_2.all[, post.idx, drop = FALSE], 1, median)
  sigma2_b.est   <- rep(var_b, Nt) # latent variance fixed to var_b by canonicalisation
  sigma2_b       <- sigma2_b.est   # keep DIC consistent with the canonical scale

  alpha.est <- apply(Z.all[, , post.idx, drop = FALSE], c(1, 2), median)
  theta.est <- round(apply(theta.all[post.idx, , drop = FALSE], 2, median), 3)
  Z.cat.est <- apply(alpha.est, 2, which.max)
  Z.est     <- Zcat_to_Zest(Z.cat.est, Nc, Nb)  # needed by DIC; reordered below if truth given

# The following block is only used in simulations where true parameter values are
# known. It aligns estimated cluster labels to the truth to enable error reporting.
if (!is.null(true.values)) {
  # Resolve cluster label-switching against the truth (for error reporting only).
  error.result <- error.fun(
    Z.cat.est = Z.cat.est,
    Nc = Nc,
    Nb = Nb,
    Z.true = Z.true,
    theta.true = theta.true,
    cindex_fun = cindex
  )

  ord        <- error.result$ind.reorder
  Z.est      <- error.result$Z.est.reorder
  Z.cat.est  <- error.result$Z.cat.reorder   # keep membership labels consistent with reordered x
  x.est      <- x.est[ord, , drop = FALSE]
  x.ci$lower <- x.ci$lower[ord, , drop = FALSE]
  x.ci$upper <- x.ci$upper[ord, , drop = FALSE]
  theta.est  <- as.numeric(theta.est)[ord]

  # Re-centre so the reference cluster sits at the origin in the truth-aligned
  # labelling (the per-draw reference used raw Gibbs labels). Shift is absorbed
  # into the intercepts so the predictor x V + a / x U + a is unchanged.
  sh         <- x.est[ref_cluster, ]
  x.est      <- sweep(x.est, 2, sh, "-")
  x.ci$lower <- sweep(x.ci$lower, 2, sh, "-")
  x.ci$upper <- sweep(x.ci$upper, 2, sh, "-")
  a_1.est    <- a_1.est + as.numeric(sh %*% V.est)
  a_2.est    <- a_2.est + as.numeric(sh %*% U.est)

  # Align estimated latent axes to the true axes (permutation + sign) so close-norm
  # axes that the canonical ordering may swap are compared like-for-like.
  al <- align_axes_to_truth(V.est, U.est, x.est, V.true.c, U.true.c)
  V.est <- al$V; U.est <- al$U; x.est <- al$x
  b.est <- sweep(b.est[, al$perm, drop = FALSE], 2, al$sign, "*")  # keep b consistent for llk
}




  llk.trace <- as.numeric(llk.all[1, ])
  if (plot_trace) {
    plot(llk.trace, ylim = c(min(llk.trace), max(llk.trace)), col = "#1B9E77", cex = 0.4)
  }
  llk.plot <- llk.trace

  # The complete-data DIC likelihoods (llk.fun for the conditional data
  # likelihood; logmix.fun / EZ_complete_llk.fun for the responsibility-based
  # DIC6 term) read the data objects below as free variables, so rebind them
  # into a local env that sees this call's data.
  dic_env <- new.env(parent = environment(llk.fun))
  dic_env$llk.fun             <- llk.fun
  dic_env$logmix.fun          <- logmix.fun
  dic_env$EZ_complete_llk.fun <- EZ_complete_llk.fun
  environment(dic_env$llk.fun)             <- dic_env
  environment(dic_env$logmix.fun)          <- dic_env
  environment(dic_env$EZ_complete_llk.fun) <- dic_env

  assign("y_1", y_1, envir = dic_env)
  assign("y_2", y_2, envir = dic_env)
  assign("Nb",  Nb,  envir = dic_env)
  assign("Nd1", Nd1, envir = dic_env)
  assign("Nd2", Nd2, envir = dic_env)
  assign("Nc",  Nc,  envir = dic_env)
  assign("M_1", M_1, envir = dic_env)
  assign("M_2", M_2, envir = dic_env)

  # Posterior point-estimate bundle used by the complete-data DIC plug-ins.
  par.hat <- list(
    "Z" = Z.est, "x" = x.est, "b" = b.est, "V" = V.est, "U" = U.est,
    "a_1" = a_1.est, "a_2" = a_2.est,
    "sigma2_y_2" = sigma2_y_2.est, "sigma2_b" = sigma2_b.est
  )

  # --- Complete-data DIC (Celeux et al., 2006: DIC4/DIC5/DIC6) ---------------
  # Treats the latent subtype labels Z as MISSING DATA and builds the deviance
  # from the complete-data joint likelihood f(y,Z|theta)=f(y|Z,theta)*prod_i
  # w_{Z_i}. DIC4 plugs in the posterior-mean estimate and integrates the
  # labels over the retained posterior draws (paper-recommended criterion);
  # DIC5 plugs in the joint MAP (best retained draw; diagnostic); DIC6 plugs
  # in the posterior-mean estimate and integrates the labels under their
  # conditional distribution given it (responsibilities). See cal_DIC_v7.R.
  cdic <- tryCatch(
    complete_DIC.fun(
      post.idx       = post.idx,
      Z.all          = Z.all,          theta.all      = theta.all,
      x.all          = x.all,          b.all          = b.all,
      V.all          = V.all,          U.all          = U.all,
      a_1.all        = a_1.all,        a_2.all        = a_2.all,
      sigma2_y_2.all = sigma2_y_2.all,
      par.hat        = par.hat,        theta.hat      = theta.est,
      llk_cond       = dic_env$llk.fun,
      llk_EZ         = dic_env$EZ_complete_llk.fun
    ),
    error = function(e) {
      message("  (complete DIC skipped: ", conditionMessage(e), ")")
      list(dic4 = NA_real_, dic5 = NA_real_, dic6 = NA_real_,
           pD4 = NA_real_, pD5 = NA_real_, pD6 = NA_real_,
           Dbar_complete = NA_real_)
    }
  )

  return(list(
    "membership" = Z.cat.est,
    "cluster center" = x.est,
    "loading to binary modality" = V.est,
    "loading to continuous modality" = U.est,
    "binary modality intercept" = a_1.est,
    "continuous modality intercept" = a_2.est,
    "membership weight" = theta.est,
    "b" = b,
    "sigma2_y_2" = sigma2_y_2.est,
    "sigma2_b" = sigma2_b.est,
    "axis order" = axis.order,
    "axis sign" = axis.sign,
    "credible intervals" = list(
      "cluster center" = x.ci,
      "loading to binary modality" = V.ci,
      "loading to continuous modality" = U.ci,
      "binary modality intercept" = a_1.ci,
      "continuous modality intercept" = a_2.ci
    ),
    "likelihood trace plot" = llk.plot,
    # Complete-data DIC (Celeux et al., 2006). dic_complete == DIC4 (the
    # paper's recommended criterion for missing-data models; posterior-mean
    # plug-in, labels integrated over the posterior draws). DIC5 (joint-MAP
    # plug-in) and DIC6 (responsibility-based label integration) are reported
    # as diagnostics. pD_* are the corresponding effective dimensions.
    "dic_complete" = cdic$dic4,
    "dic4_complete" = cdic$dic4, "pD4_complete" = cdic$pD4,
    "dic5_complete" = cdic$dic5, "pD5_complete" = cdic$pD5,
    "dic6_complete" = cdic$dic6, "pD6_complete" = cdic$pD6,
    "Dbar_complete" = cdic$Dbar_complete
  ))
}
