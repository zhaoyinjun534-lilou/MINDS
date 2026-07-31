#' Permutation Utilities and Cluster Assignment Helpers
#'
#' This file contains utility functions for the MINDS algorithm,
#' including permutation generation and cluster assignment conversion.

#' Generate All Permutations Recursively
#'
#' @param n Integer; number of elements to permute.
#'
#' @return A matrix where each row is a permutation of 1:n.
permute_indices <- function(n) {
  if (n == 1) return(matrix(1, nrow = 1))
  prev <- permute_indices(n - 1)
  out <- NULL
  for (i in seq_len(nrow(prev))) {
    v <- prev[i, ]
    for (pos in seq_len(n)) {
      out <- rbind(out, append(v, n, after = pos - 1))
    }
  }
  out
}

#' List All Permutations
#'
#' @param n Integer; number of elements.
#'
#' @return A list of permutations, each as an integer vector.
permn <- function(n) {
  perms <- permute_indices(n)
  lapply(seq_len(nrow(perms)), function(i) as.integer(perms[i, ]))
}

#' Convert Cluster Labels to Binary Matrix
#'
#' Converts categorical cluster assignments to a binary membership matrix.
#'
#' @param Z.cat.est Integer vector of cluster assignments.
#' @param Nc Integer; number of clusters.
#' @param Nb Integer; number of subjects.
#'
#' @return An Nc x Nb binary matrix where `Z.est[k, i] = 1` if subject i is in cluster k.
Zcat_to_Zest <- function(Z.cat.est, Nc, Nb) {
  Z.est <- matrix(0, Nc, Nb)
  for (i in seq_len(Nb)) {
    Z.est[Z.cat.est[i], i] <- 1
  }
  Z.est
}

#' Permutation-based Error Evaluation
#'
#' Evaluates clustering errors under all possible label permutations.
#'
#' @param Z.cat.est Integer vector of estimated cluster assignments.
#' @param Nc Integer; number of clusters.
#' @param Nb Integer; number of subjects.
#' @param Z.true Binary matrix of true cluster assignments.
#' @param theta.true Numeric vector of true cluster weights.
#' @param cindex_fun Function to compute clustering index (e.g., Jaccard).
#' @param Z.cat.true Optional; integer vector of true cluster assignments.
#'
#' @return A list of length Nc! with error metrics for each permutation.
permu.fun <- function(Z.cat.est, Nc, Nb, Z.true, theta.true, cindex_fun, Z.cat.true = NULL) {
  if (is.null(Z.cat.true)) {
    Z.cat.true <- apply(Z.true, 2, which.max)
  }
  perm_list <- permn(Nc)
  lapply(seq_along(perm_list), function(ind) {
    index.est <- perm_list[[ind]]
    
    Z.est <- Zcat_to_Zest(Z.cat.est, Nc, Nb)
    
    Z.est.reorder <- apply(Z.est, 2, function(z) z[index.est])
    Z.cat.reorder <- apply(Z.est.reorder, 2, which.max)
    
    B_error <- mean(apply(Z.est.reorder != Z.true, 2, function(a) a %*% theta.true))
    C_error <- mean(Z.cat.reorder != Z.cat.true)
    
    names(Z.cat.reorder) <- seq_len(Nb)
    names(Z.cat.true) <- seq_len(Nb)
    ci <- cindex_fun(
      clV1 = Z.cat.true,
      clV2 = Z.cat.reorder,
      self = FALSE,
      minSZ = 1,
      method = "jaccard"
    )
    JCC_error <- 1 - ci$Jaccard_Index
    
    list(C_error = C_error, B_error = B_error, JCC_error = JCC_error, index.est = index.est)
  })
}

#' Find Best Permutation of Cluster Labels
#'
#' Selects the label permutation that minimizes clustering error.
#'
#' @param Z.cat.est Integer vector of estimated cluster assignments.
#' @param Nc Integer; number of clusters.
#' @param Nb Integer; number of subjects.
#' @param Z.true Binary matrix of true cluster assignments.
#' @param theta.true Numeric vector of true cluster weights.
#' @param cindex_fun Function to compute clustering index.
#' @param Z.cat.true Optional; integer vector of true cluster assignments.
#'
#' @return A list containing reordered estimates and error metrics.
error.fun <- function(Z.cat.est, Nc, Nb, Z.true, theta.true, cindex_fun, Z.cat.true = NULL) {
  if (is.null(Z.cat.true)) {
    Z.cat.true <- apply(Z.true, 2, which.max)
  }
  Z.est <- Zcat_to_Zest(Z.cat.est, Nc, Nb)
  
  aa <- permu.fun(
    Z.cat.est = Z.cat.est,
    Nc = Nc,
    Nb = Nb,
    Z.true = Z.true,
    theta.true = theta.true,
    cindex_fun = cindex_fun,
    Z.cat.true = Z.cat.true
  )
  ind <- which.min(sapply(aa, function(a) a$C_error))
  
  Z.est.reorder <- apply(Z.est, 2, function(z) z[aa[[ind]]$index.est])
  Z.cat.reorder <- apply(Z.est.reorder, 2, which.max)
  
  list(
    B_error = aa[[ind]]$B_error,
    C_error = aa[[ind]]$C_error,
    JCC_error = aa[[ind]]$JCC_error,
    Z.est.reorder = Z.est.reorder,
    Z.cat.reorder = Z.cat.reorder,
    ind.reorder = aa[[ind]]$index.est
  )
}

#' Generate True Parameter Values
#'
#' Generates true parameter values for simulation: cluster centers (x),
#' factor loadings (V, U), and intercepts (a_1, a_2).
#'
#' @param seed.no Integer random seed.
#' @param Nd1 Integer; number of binary items.
#' @param Nd2 Integer; number of continuous items.
#' @param Nt Integer; number of latent constructs.
#' @param Nc Integer; number of latent clusters.
#'
#' @return A list containing:
#' \itemize{
#'   \item `x`: Cluster centers (Nc x Nt).
#'   \item `V`: Binary loadings (Nt x Nd1).
#'   \item `U`: Continuous loadings (Nt x Nd2).
#'   \item `a_1`: Binary intercepts (Nd1).
#'   \item `a_2`: Continuous intercepts (Nd2).
#' }
par.fun <- function(seed.no, Nd1, Nd2, Nt, Nc){
  
  set.seed(seed.no)
  
  # Item-specific fixed effects (intercepts)
  a_1 <- rnorm(Nd1, mean = 0, sd = 0.5)*(-1)
  a_2 <- rnorm(Nd2, mean = 0, sd = 0.5)*(-1)

  x <- replicate(Nt,
                 runif(Nc, 0, 2)) # for Nt latent constructs
  x <- apply(x, 2, scale)*2
  

  # Modality-specific linear maps from latent space to item space
  V.raw.item <-  (replicate(Nt, runif(Nd1, 0, 2)))
  U.raw.item <-  (replicate(Nt, runif(Nd2, 0, 2)))
  
  V <- t(V.raw.item)
  U <- t(U.raw.item)
  
  return(list(
    'x'   = x,
    'V'   = V,
    'U'   = U,
    'a_1' = a_1,
    'a_2' = a_2
  ))
}

#' Map MINDS parameters to their identified (canonical) form
#'
#' The MINDS likelihood depends on the parameters only through the invariant
#' linear predictors \eqn{m_i = (Z_i^T x + b_i) L + a}, where \eqn{L = [V | U]}.
#' Fixing the latent variance to a common value `var_b` makes the latent
#' covariance \eqn{var_b \cdot I}, which is invariant to any orthogonal rotation;
#' X, V and U are then identified only up to a rotation (plus per-axis sign and an
#' overall location shift), even though \eqn{X L + a} -- and hence the intercepts
#' a -- are well determined. This function removes all of these indeterminacies
#' deterministically, so no optimisation and no knowledge of the true values is
#' required. Cluster label-switching is handled separately (e.g. by `error.fun`).
#'
#' Canonical conventions applied, in order:
#'   1. scale    : each latent axis rescaled so that var(b) = var_b on that axis;
#'   2. rotation : rotate to the loadings' principal axes (eigenbasis of L L^T),
#'                 removing the orthogonal-rotation freedom and ordering axes by
#'                 descending loading energy;
#'   3. sign     : the largest-magnitude loading on each axis is made positive;
#'   4. location : the reference cluster centre is set to 0.
#'
#' Every step is orthogonal or diagonal in the latent space, so the independence
#' (diagonal covariance) of `b` is preserved.
#'
#' @param x Nc x Nt matrix of cluster centres.
#' @param b Nb x Nt matrix of subject random effects.
#' @param V Nt x Nd1 binary loading matrix.
#' @param U Nt x Nd2 continuous loading matrix.
#' @param a_1,a_2 length-Nd1 / length-Nd2 intercept vectors.
#' @param sigma2_b length-Nt latent variances (used for the scale step).
#' @param var_b target latent variance imposed on each axis by the scale step
#'   (the identifiability constraint on the latent metric). Default `1`
#'   (unit-variance constructs); set e.g. `0.1` to match the data-generating scale.
#' @param ref_cluster integer in 1:Nc; the cluster whose centre is fixed at the
#'   origin (location constraint). Other centres are read as contrasts against it.
#' @return list with canonical `x`, `b`, `V`, `U`, `a_1`, `a_2`, and `sigma2_b`
#'   (all equal to `var_b` after canonicalisation), plus the `order` and `sign`.
#' @export
canonicalize_minds <- function(x, b, V, U, a_1, a_2, sigma2_b,
                               var_b = 1, ref_cluster = 1L) {
  Nt  <- ncol(x); Nc <- nrow(x); Nd1 <- ncol(V); Nd2 <- ncol(U)
  if (ref_cluster < 1L || ref_cluster > Nc) stop("`ref_cluster` must be in 1:Nc.")
  L   <- cbind(V, U)
  a   <- c(a_1, a_2)
  s2b <- as.numeric(sigma2_b)
  ## (1) scale: fix the latent variance on each axis to `var_b`
  ##     (diagonal rescaling -> keeps b independent). var_b = 1 gives
  ##     unit-variance latent constructs; e.g. var_b = 0.1 fixes the constraint
  ##     at the latent scale used to generate the data. The factor is
  ##     s_t = sqrt(sigma2_b_t / var_b), applied as x,b /= s and L_row *= s, so
  ##     Var(b_t) = sigma2_b_t / s_t^2 = var_b.
  s   <- sqrt(pmax(s2b, .Machine$double.eps) / var_b)
  x   <- sweep(x, 2, s, "/")
  b   <- sweep(b, 2, s, "/")
  L   <- s * L                       # row t of L multiplied by s[t]
  s2b <- rep(var_b, Nt)

  ## (2) rotation: rotate to the loadings' principal axes. After equal-variance
  ##     whitening the latent covariance is var_b * I, which is invariant to any
  ##     orthogonal rotation; X, V, U are then identified only up to a rotation
  ##     (they "mingle"), even though their product X L + a -- and hence a -- is
  ##     well determined. Rotating to the eigenbasis of L L^T removes this
  ##     freedom and orders axes by descending loading energy (the eigenvalues).
  ##     Requires distinct axis energies; ties leave a residual rotation within
  ##     the tied subspace.
  ev <- eigen(L %*% t(L), symmetric = TRUE)   # eigenvalues in decreasing order
  W  <- ev$vectors                            # eta -> eta W ;  L -> t(W) L
  x  <- x %*% W
  b  <- b %*% W
  L  <- t(W) %*% L
  ord <- seq_len(Nt)                          # axes already ordered by energy

  ## (3) sign: make the largest-magnitude loading on each axis positive.
  sgn <- rep(1, Nt)
  for (t in seq_len(Nt)) {
    j <- which.max(abs(L[t, ]))
    if (L[t, j] < 0) {
      L[t, ] <- -L[t, ]; x[, t] <- -x[, t]; b[, t] <- -b[, t]; sgn[t] <- -1
    }
  }

  ## (4) location: place the reference cluster `ref_cluster` at the origin
  ##     (its center becomes 0; the shift is absorbed into the intercepts a so
  ##     the predictor x L + a is unchanged). Other centers read as contrasts
  ##     against the reference cluster.
  sh <- x[ref_cluster, ]
  x <- sweep(x, 2, sh, "-")
  a <- a + as.numeric(sh %*% L)

  list(
    x = x, b = b,
    V = L[, seq_len(Nd1), drop = FALSE],
    U = L[, Nd1 + seq_len(Nd2), drop = FALSE],
    a_1 = a[seq_len(Nd1)], a_2 = a[Nd1 + seq_len(Nd2)],
    sigma2_b = s2b, order = ord, sign = sgn
  )
}

# All permutations of 1:n (used to align latent axes; n = Nt is small).
.perms <- function(n) {
  if (n == 1L) return(matrix(1L, 1, 1))
  sub <- .perms(n - 1L)
  out <- NULL
  for (i in seq_len(n)) {
    rest <- sub
    rest[rest >= i] <- rest[rest >= i] + 1L
    out <- rbind(out, cbind(i, rest))
  }
  out
}

#' Align estimated latent axes to the true axes (simulation evaluation only)
#'
#' The canonical form orders axes by loading norm and signs them by the largest
#' loading. When two true axes have similar loading norms, estimation noise can
#' swap their order (or flip a sign) relative to the truth, so a like-for-like
#' RMSE needs the estimated axes matched to the true ones. This finds the axis
#' permutation and per-axis signs minimising the loading distance to the truth
#' (the axis analogue of resolving cluster label-switching) and applies them to
#' `V`, `U` and `x`. Intercepts `a` are per item, not per axis, so are unchanged.
#' Uses the true values, so it is for simulation comparison, not real-data fits.
#'
#' @param V.est,U.est,x.est canonical estimates (V,U: Nt x Nd; x: Nc x Nt).
#' @param V.true,U.true canonical true loadings (Nt x Nd).
#' @return list with aligned `V`, `U`, `x`, and the chosen `perm` and `sign`.
#' @export
align_axes_to_truth <- function(V.est, U.est, x.est, V.true, U.true) {
  Nt    <- nrow(V.true)
  Lest  <- cbind(V.est, U.est)
  Ltrue <- cbind(V.true, U.true)
  perms <- .perms(Nt)
  signs <- as.matrix(expand.grid(rep(list(c(1, -1)), Nt)))
  best <- Inf; best.p <- seq_len(Nt); best.s <- rep(1, Nt)
  for (pi in seq_len(nrow(perms))) {
    p <- perms[pi, ]
    for (si in seq_len(nrow(signs))) {
      s <- signs[si, ]
      d <- sum((s * Lest[p, , drop = FALSE] - Ltrue)^2)  # row t scaled by s[t]
      if (d < best) { best <- d; best.p <- p; best.s <- s }
    }
  }
  list(
    V    = best.s * V.est[best.p, , drop = FALSE],
    U    = best.s * U.est[best.p, , drop = FALSE],
    x    = sweep(x.est[, best.p, drop = FALSE], 2, best.s, "*"),
    perm = best.p, sign = best.s
  )
}

