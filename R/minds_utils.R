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
#' @return An Nc × Nb binary matrix where Z.est[k, i] = 1 if subject i is in cluster k.
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
#'   \item `x`: Cluster centers (Nc × Nt).
#'   \item `V`: Binary loadings (Nt × Nd1).
#'   \item `U`: Continuous loadings (Nt × Nd2).
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
