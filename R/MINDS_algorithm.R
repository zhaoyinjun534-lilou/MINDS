enforce_identifiability <- function(x.est, V.est, U.est, a_1.est, a_2.est,
                                    theta.est, Nd1, Nd2, Nt) {
  L.raw <- rbind(t(V.est), t(U.est))
  qr.L <- qr(L.raw)
  Q <- qr.Q(qr.L)[, seq_len(Nt), drop = FALSE]
  R <- qr.R(qr.L)[seq_len(Nt), seq_len(Nt), drop = FALSE]

  sign.d <- sign(diag(R))
  sign.d[sign.d == 0] <- 1
  Q <- sweep(Q, 2, sign.d, "*")
  R <- diag(sign.d) %*% R

  x.est <- x.est %*% t(R)
  V.est <- t(Q[seq_len(Nd1), , drop = FALSE])
  U.est <- t(Q[Nd1 + seq_len(Nd2), , drop = FALSE])

  theta.weight <- as.numeric(theta.est)
  theta.weight <- theta.weight / sum(theta.weight)
  x.shift <- as.numeric(theta.weight %*% x.est)
  x.est <- sweep(x.est, 2, x.shift, FUN = "-")
  a_1.est <- a_1.est - as.numeric(x.shift %*% V.est)
  a_2.est <- a_2.est - as.numeric(x.shift %*% U.est)

  list(
    x.est = x.est,
    V.est = V.est,
    U.est = U.est,
    a_1.est = a_1.est,
    a_2.est = a_2.est
  )
}

#' Fit MINDS Algorithm For Mixed Binary-Continuous Data
#'
#' Runs the Bayesian MINDS model using binary outcomes `y_1` and continuous
#' outcomes `y_2`, with user-specified prior hyperparameters.
#'
#' @param y_1 Binary outcome matrix with subjects in rows and items in columns.
#' @param y_2 Continuous outcome matrix with subjects in rows and items in columns.
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
#' @param init_seed Seed used for initialization and imputation.
#' @param plot_trace Logical; if `TRUE`, draw the likelihood trace.
#'
#' @return A named list with elements:
#' \itemize{
#'   \item `membership`: Estimated cluster assignment for each subject.
#'   \item `cluster center`: Estimated cluster centers.
#'   \item `loading to binary modality`: Estimated binary loading.
#'   \item `loading to continuous modality`: Estimated continuous loading.
#'   \item `binary modality intercept`: Estimated binary intercepts.
#'   \item `continuous modality intercept`: Estimated continuous intercepts.
#'   \item `membership weight`: Estimated cluster weights.
#'   \item `likelihood trace plot`: Numeric likelihood trace across iterations.
#'   \item `ic`: Information criterion from `dic.fun(par.est)$ic`.
#' }
#'
#' @examples
#' \dontrun{
#' data(data_mixed)
#' out <- MINDS_algorithm(
#'   y_1 = data_mixed$y_1,
#'   y_2 = data_mixed$y_2,
#'   Nc = 4,
#'   Nt = 3,
#'   iter.max = 50,
#'   plot_trace = FALSE
#' )
#' names(out)
#' }
#'
#' @export
MINDS_algorithm <- function(
    y_1,
    y_2,
    Nc,
    Nt,
    iter.max = 500,
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
    init_seed = 123,
    plot_trace = TRUE) {
  Nd1 <- ncol(y_1)
  Nd2 <- ncol(y_2)
  Nb <- nrow(y_1)

  if (nrow(y_2) != Nb) {
    stop("y_1 and y_2 must have the same number of rows.")
  }

  x.ini <- MASS::mvrnorm(Nc, rep(mu_x, Nt), Sigma = diag(rep(sigma2_x, Nt)))
  sigma2_b.ini <- EDISON::rinvgamma(Nt, IG_b.shape, IG_b.scale)
  sigma2_y_2.ini <- EDISON::rinvgamma(Nd2, IG_y_2.shape, IG_y_2.scale)
  b.ini <- MASS::mvrnorm(Nb, rep(0, Nt), Sigma = diag(sigma2_b.ini))

  y_1_mean <- apply(y_1, 2, function(a) mean(a, na.rm = TRUE))
  odds.fun <- function(a) log(a) / log(1 - a)
  y_1_mean <- pmin(pmax(y_1_mean, 1e-6), 1 - 1e-6)
  a_1.ini <- odds.fun(y_1_mean[1]) - odds.fun(y_1_mean)

  y_2_mean <- apply(y_2, 2, function(a) mean(a, na.rm = TRUE))
  a_2.ini <- y_2_mean[1] - y_2_mean

  set.seed(init_seed)
  V.ini <- t(replicate(Nt, runif(Nd1, 0, 1)))
  U.ini <- t(replicate(Nt, runif(Nd2, 0, 1)))

  y <- data.frame(apply(y_1, 2, as.logical), y_2)
  imp <- mice::mice(y, m = 5, maxit = 50, method = "pmm", seed = init_seed, printFlag = FALSE)
  y.complete <- mice::complete(imp, 1)
  y.cor <- cor(y.complete)

  f5 <- psych::fa(y.cor, nfactors = min(5, ncol(y.cor)))
  my.score <- psych::factor.scores(y.complete, f5, method = "tenBerge")$scores

  cl <- kmeans(my.score, Nc)
  Z.cat.est <- cl$cluster
  Z.est <- matrix(0, Nc, Nb)
  for (i in seq_len(Nb)) {
    Z.est[Z.cat.est[i], i] <- 1
  }
  Z.ini <- Z.est

  theta.count <- tabulate(Z.cat.est, nbins = Nc)
  theta.ini <- theta.count / sum(theta.count)

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
    temp <- t(apply(temp1, 1, function(ti) ti - a_1))
    w <- apply(temp, c(1, 2), function(s) BayesLogit::rpg(num = 1, h = 1, z = s))

    V_a <- 1 / (apply(w * M_1, 2, sum) + 1 / sigma2_a_1)
    temp <- (t(Z) %*% x + b) %*% V
    r <- temp - (y_1 - 1 / 2) / w
    Mu_a <- apply(w * M_1 * r, 2, function(a) sum(a, na.rm = TRUE)) * V_a
    a_1 <- sapply(seq_len(Nd1), function(i) rnorm(1, Mu_a[i], sqrt(V_a[i])))

    V_a <- 1 / (apply(M_2 * rep(1 / sigma2_y_2, each = nrow(M_2)), 2, sum) + 1 / sigma2_a_2)
    temp <- (t(Z) %*% x + b) %*% U
    r <- temp - y_2
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
      temp2 <- t(apply(temp1, 1, function(tt) tt - a_1))
      phi <- (y_1 - 1 / 2) / w - temp2

      temp3 <- apply(w * M_1, 2, function(s) s * temp_A)
      D_V_ti <- apply(temp3 * phi, 2, function(a) sum(a, na.rm = TRUE)) + mu_v / sigma2_v

      Mu_V_ti <- B_V_ti * D_V_ti
      V_ti <- MASS::mvrnorm(1, Mu_V_ti, diag(B_V_ti))
      V[ti, ] <- V_ti
    }

    for (ti in seq_len(Nt)) {
      x.sub <- x
      x.sub[, ti] <- 0
      B_ti <- diag(sigma2_x, Nc)

      temp <- (t(Z) %*% x.sub + b) %*% V
      temp1 <- t(apply(temp, 1, function(tt) tt - a_1))
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
      temp1 <- t(apply(temp, 1, function(tt) tt - a_2))
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

    for (ti in seq_len(Nt)) {
      b.sub <- b
      b.sub[, ti] <- 0

      temp <- (t(Z) %*% x + b.sub) %*% V
      temp1 <- t(apply(temp, 1, function(tt) tt - a_1))
      C_ti <- (y_1 - 1 / 2) / w - temp1

      temp <- (t(Z) %*% x + b.sub) %*% U
      temp1 <- t(apply(temp, 1, function(tt) tt - a_2))
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
      temp2 <- t(apply(temp1, 1, function(tt) tt - a_2))
      phi <- y_2 - temp2

      temp3 <- apply(phi * M_2, 2, function(s) sum(s * temp_A, na.rm = TRUE))
      D_U_ti <- temp3 / sigma2_y_2 + mu_u / sigma2_u

      Mu_U_ti <- B_U_ti * D_U_ti
      U_ti <- MASS::mvrnorm(1, Mu_U_ti, diag(B_U_ti))
      U[ti, ] <- U_ti
    }

    if (i.iter < iter.max * 0.2) {
      Z <- Z.ini
      alpha <- t(Z.ini)
    }

    if (i.iter >= iter.max * 0.2) {
      alpha <- NULL
      for (k in seq_len(Nc)) {
        temp1 <- t(apply(b, 1, function(b_i) x[k, ] + b_i)) %*% V
        temp2 <- t(apply(temp1, 1, function(tt) tt - a_1)) - (y_1 - 1 / 2) / w
        temp2 <- ifelse(is.na(temp2), 0, temp2)
        p1_k <- exp(-w * M_1 / 2 * temp2^2)

        temp1 <- t(apply(b, 1, function(b_i) x[k, ] + b_i)) %*% U
        temp2 <- t(apply(temp1, 1, function(tt) tt - a_2))

        sigma2_y_2_matrix <- matrix(rep(sigma2_y_2, Nb), nrow = Nb, byrow = TRUE)
        aa <- temp2 - y_2
        aa <- ifelse(is.na(aa), 0, aa)
        p2_k <- (1 / sqrt(2 * pi * sigma2_y_2_matrix))^M_2 * exp(-(aa)^2 / 2 / sigma2_y_2_matrix)

        alpha.k <- apply(p1_k, 1, prod) * apply(p2_k, 1, prod) * theta[k]
        alpha <- cbind(alpha, alpha.k)
      }

      alpha <- alpha / apply(alpha, 1, sum)
      Z <- apply(alpha, 1, function(p) rmultinom(1, 1, p))
    }

    theta <- LaplacesDemon::rdirichlet(1, apply(Z, 1, sum) + p.theta.prior)
    sigma2_b <- apply(b, 2, function(b_ti) EDISON::rinvgamma(1, Nb / 2 + IG_b.shape, sum(b_ti^2) / 2 + IG_b.scale))

    temp <- (t(Z) %*% x + b) %*% U
    temp1 <- t(apply(temp, 1, function(tt) tt - a_2))
    C <- y_2 - temp1
    aa <- apply(C^2 * M_2, 2, function(a) sum(a, na.rm = TRUE)) / 2
    sigma2_y_2 <- unlist(lapply(seq_len(Nd2), function(j) {
      EDISON::rinvgamma(1, Nb / 2 + IG_y_2.shape, aa[j] + IG_y_2.scale)
    }))

    temp <- (t(Z) %*% x + b) %*% V
    temp1 <- t(apply(temp, 1, function(tt) tt - a_1))
    p <- exp(temp1) / (1 + exp(temp1))
    llk_binary <- sum(apply(M_1 * y_1 * log(p) + log(1 - p) * M_1 * (1 - y_1), 1,
                            function(a) sum(a, na.rm = TRUE)), na.rm = TRUE)

    temp <- (t(Z) %*% x + b) %*% U
    temp1 <- t(apply(temp, 1, function(tt) tt - a_2))
    logd <- sapply(seq_len(Nd2), function(j) {
      M_2[, j] * dnorm(y_2[, j], mean = temp1[, j], sd = sqrt(sigma2_y_2[j]), log = TRUE)
    })
    llk_continuous <- sum(apply(logd, 2, function(a) sum(a, na.rm = TRUE)), na.rm = TRUE)
    llk <- llk_binary + llk_continuous

    x.all <- abind::abind(x.all, x, along = 3)
    Z.all <- abind::abind(Z.all, Z, along = 3)
    V.all <- abind::abind(V.all, V, along = 3)
    U.all <- abind::abind(U.all, U, along = 3)

    theta.all <- rbind(theta.all, theta)
    sigma2_y_2.all <- cbind(sigma2_y_2.all, sigma2_y_2)
    a_2.all <- cbind(a_2.all, a_2)
    a_1.all <- cbind(a_1.all, a_1)
    alpha.all <- abind::abind(alpha.all, alpha, along = 3)
    llk.all <- cbind(llk.all, llk)
    sigma2_b.all <- rbind(sigma2_b.all, sigma2_b)
  }

  iter.all <- dim(x.all)[3]
  n.burn <- min(round(iter.all * 0.8, 0), iter.all - 1)
  post.idx <- seq.int(n.burn + 1, iter.all)

  x.est <- apply(x.all[, , post.idx, drop = FALSE], c(1, 2), median)
  a_2.est <- apply(a_2.all[, post.idx, drop = FALSE], 1, median)
  a_1.est <- apply(a_1.all[, post.idx, drop = FALSE], 1, median)
  V.est <- apply(V.all[, , post.idx, drop = FALSE], c(1, 2), median)
  U.est <- apply(U.all[, , post.idx, drop = FALSE], c(1, 2), median)

  alpha.est <- apply(Z.all[, , post.idx, drop = FALSE], c(1, 2), median)
  theta.est <- round(apply(theta.all[post.idx, , drop = FALSE], 2, median), 3)
  Z.cat.est <- apply(alpha.est, 2, which.max)

  iden.par <- enforce_identifiability(
    x.est = x.est,
    V.est = V.est,
    U.est = U.est,
    a_1.est = a_1.est,
    a_2.est = a_2.est,
    theta.est = theta.est,
    Nd1 = Nd1,
    Nd2 = Nd2,
    Nt = Nt
  )

  x.est <- iden.par$x.est
  V.est <- iden.par$V.est
  U.est <- iden.par$U.est
  a_1.est <- iden.par$a_1.est
  a_2.est <- iden.par$a_2.est

  llk.trace <- as.numeric(llk.all[1, ])
  if (plot_trace) {
    plot(llk.trace, ylim = c(min(llk.trace), max(llk.trace)), col = "#1B9E77", cex = 0.4)
  }
  llk.plot <- llk.trace

  par.est <- list(
    "Z" = Z.est,
    "x" = x.est,
    "b" = b,
    "V" = V.est,
    "U" = U.est,
    "a_1" = a_1.est,
    "a_2" = a_2.est,
    "sigma2_y_2" = sigma2_y_2,
    "sigma2_b" = sigma2_b
  )

  if (!exists("dic.fun", mode = "function")) {
    stop("dic.fun is not available. Ensure cal_DIC_v7.R is loaded in the package namespace.")
  }

  # Use mutable local env for DIC helpers because package namespace is locked.
  dic_env <- new.env(parent = environment(dic.fun))
  dic_env$llk.fun <- llk.fun
  dic_env$par.simu.fun <- par.simu.fun
  dic_env$dic.fun <- dic.fun
  environment(dic_env$llk.fun) <- dic_env
  environment(dic_env$par.simu.fun) <- dic_env
  environment(dic_env$dic.fun) <- dic_env

  assign("y_1", y_1, envir = dic_env)
  assign("y_2", y_2, envir = dic_env)
  assign("Nb", Nb, envir = dic_env)
  assign("Nd1", Nd1, envir = dic_env)
  assign("Nd2", Nd2, envir = dic_env)
  assign("Nt", Nt, envir = dic_env)
  assign("Nc", Nc, envir = dic_env)
  assign("M_1", M_1, envir = dic_env)
  assign("M_2", M_2, envir = dic_env)
  assign("sigma2_b", sigma2_b, envir = dic_env)
  assign("theta", theta, envir = dic_env)
  assign("sigma2_a_1", sigma2_a_1, envir = dic_env)
  assign("sigma2_a_2", sigma2_a_2, envir = dic_env)
  assign("sigma2_v", sigma2_v, envir = dic_env)
  assign("mu_v", mu_v, envir = dic_env)
  assign("mu_x", mu_x, envir = dic_env)
  assign("sigma2_x", sigma2_x, envir = dic_env)
  assign("sigma2_u", sigma2_u, envir = dic_env)
  assign("mu_u", mu_u, envir = dic_env)
  assign("IG_b.shape", IG_b.shape, envir = dic_env)
  assign("IG_b.scale", IG_b.scale, envir = dic_env)
  assign("p.theta.prior", p.theta.prior, envir = dic_env)
  assign("IG_y_2.shape", IG_y_2.shape, envir = dic_env)
  assign("IG_y_2.scale", IG_y_2.scale, envir = dic_env)

  ic.value <- dic_env$dic.fun(par.est)$ic

  return(list(
    "membership" = Z.cat.est,
    "cluster center" = x.est,
    "loading to binary modality" = V.est,
    "loading to continuous modality" = U.est,
    "binary modality intercept" = a_1.est,
    "continuous modality intercept" = a_2.est,
    "memberhip weight" = theta.est,
    "likelihood trace plot" = llk.plot,
    "ic" = ic.value
  ))
}
