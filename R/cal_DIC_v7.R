
# #Definition of DIC
# $DIC=-4E_{{\theta|y}}[logp(y|\theta)]+ 2logp(y|\bar\theta)$, where $\bar\theta$ is the posterior mean, i.e. Bayes estimator, and $E_{{\theta|y}}[logp(y|\theta)]$ is calculated using simulations $\theta ^s$ as $\frac{1}{S}\sum logp(y|\theta^s))$. 
# 
# Take the average of deviance values over all posterior samples:
#   $\overline{D} = \frac{1}{S} \sum_{s=1}^{S} D(\theta_s)$, where S is the number of MCMC samples.
# 
# #All parameters
# 
# Z, x, b, V, U, a_1, a_2, $\sigma_j$, $\sigma_b$, $\theta$

#log-likelihood function----

llk.fun <- function(par.est){ #log-likelihood function
  Z <- par.est$Z
  x <- par.est$x
  b <- par.est$b
  V <- par.est$V
  U <- par.est$U
  a_1 <- par.est$a_1
  a_2 <- par.est$a_2
  sigma2_y_2 <- par.est$sigma2_y_2
  sigma2_b <- par.est$sigma2_b
  
  temp <- (t(Z) %*% x + b) %*% V; #(Z^T*x+b)*V, Nb*Nd1
  temp1 <- t(apply(temp, 1, function (t) {t + a_1}));#(Z^T*x+b)V+a
  p <- exp(temp1)/(1 + exp(temp1))
  # Clamp away from {0,1}: a saturated p makes log(p) or log(1-p) = -Inf, which
  # na.rm does NOT drop, poisoning the whole sum to -Inf (and hence ic = +Inf).
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)

  aa1 <- sum(apply(M_1*y_1*log(p)+log(1-p)*M_1*(1-y_1), 1,
                   function(a) sum(a, na.rm = TRUE)), 
             na.rm = TRUE)#binary
  
  
  temp <- (t(Z) %*% x + b) %*% U; #(Z^T*x+b)*U, Nb*Nd2
  temp1 <- t(apply(temp, 1, function (t) {t + a_2}));#(Z^T*x+b)U+a
  
  logd <- sapply(1:Nd2, function(j) {
    M_2[, j] * dnorm(y_2[, j], mean = temp1[, j], sd = sqrt(sigma2_y_2[j]), log = TRUE)
  })
  aa2 <- sum(apply(logd,2,function(a) sum(a, na.rm = TRUE)), na.rm = TRUE)
  aa <- aa1 + aa2
  return(list('llk_total'=aa, 'llk_binary'=aa1, 'llk_continuous'=aa2))
}

# Observed-data (label-marginalised) log-likelihood -------------------------
#   log p(y_i | .) = log sum_k theta_k * prod_j p(y_ij | Z_i = k, .)
# Marginalising the latent cluster labels is the appropriate likelihood for
# DIC/BPIC model selection in mixture models (Celeux et al., 2006). The earlier
# `llk.fun` conditions on a plugged-in Z (complete-data likelihood), which is
# label-dependent and biases the choice of Nc. This function does NOT condition
# on Z, so no membership matrix is required. It reads y_1, y_2, M_1, M_2, Nb,
# Nc, Nd1, Nd2 from its enclosing environment (same as `llk.fun`).
#
# @param par  list with x (Nc x Nt), b (Nb x Nt), V (Nt x Nd1), U (Nt x Nd2),
#             a_1 (len Nd1), a_2 (len Nd2), sigma2_y_2 (len Nd2).
# @param theta length-Nc mixture weights.

# Per-subject, per-cluster complete-data log-density matrix (Nb x Nc),
#   logmix[i, k] = log(theta_k) + log p(y_i | Z_i = k, par),
# shared by the marginalised likelihood (llk_obs.fun) and the responsibility-
# based DIC6 plug-in term (EZ_complete_llk.fun).
logmix.fun <- function(par, theta) {
  x   <- par$x;   b <- par$b;   V <- par$V;   U <- par$U
  a_1 <- par$a_1; a_2 <- par$a_2
  s2  <- par$sigma2_y_2

  bV  <- b %*% V                                   # Nb x Nd1
  bU  <- b %*% U                                   # Nb x Nd2
  sd2 <- matrix(sqrt(s2), Nb, Nd2, byrow = TRUE)   # Nb x Nd2

  logmix <- matrix(NA_real_, Nb, Nc)               # log(theta_k) + loglik_ik
  for (k in seq_len(Nc)) {
    # binary part: logistic Bernoulli, linear predictor (x_k + b_i) V + a_1
    xkV   <- matrix(as.numeric(x[k, ] %*% V), Nb, Nd1, byrow = TRUE)
    lin_b <- sweep(bV + xkV, 2, a_1, "+")
    p     <- 1 / (1 + exp(-lin_b))
    p     <- pmin(pmax(p, 1e-12), 1 - 1e-12)
    lb    <- y_1 * log(p) + (1 - y_1) * log(1 - p)
    lb    <- ifelse(M_1 == 1 & is.finite(lb), lb, 0)   # mask missing / non-finite
    ll_bin <- rowSums(lb)

    # continuous part: Gaussian, linear predictor (x_k + b_i) U + a_2
    xkU   <- matrix(as.numeric(x[k, ] %*% U), Nb, Nd2, byrow = TRUE)
    lin_c <- sweep(bU + xkU, 2, a_2, "+")
    ld    <- matrix(dnorm(as.vector(y_2), mean = as.vector(lin_c),
                          sd = as.vector(sd2), log = TRUE), Nb, Nd2)
    ld    <- ifelse(M_2 == 1 & is.finite(ld), ld, 0)
    ll_con <- rowSums(ld)

    logmix[, k] <- log(theta[k]) + ll_bin + ll_con
  }
  logmix
}

llk_obs.fun <- function(par, theta) {
  logmix <- logmix.fun(par, theta)
  # log-sum-exp over clusters, per subject, then sum over subjects
  m    <- apply(logmix, 1, max)
  ll_i <- m + log(rowSums(exp(logmix - m)))
  sum(ll_i)
}

# Responsibility-based label expectation -- the DIC6 plug-in term (Celeux et
# al., 2006, p.657; mixture approximation p.664):
#   E_Z[ log f(y, Z | theta) | y, theta ]
#     = sum_i sum_k t_ik * [ log(w_k) + log p(y_i | Z_i = k, theta) ],
# with responsibilities under the PLUG-IN estimate,
#   t_ik = w_k p(y_i|k) / sum_k' w_k' p(y_i|k').
# NB: this label integration conditions on theta-hat. Averaging over the
# posterior label draws Z^s ~ p(Z|y) instead belongs to DIC4, not DIC6.
EZ_complete_llk.fun <- function(par, theta) {
  logmix <- logmix.fun(par, theta)
  m    <- apply(logmix, 1, max)
  t_ik <- exp(logmix - m)
  t_ik <- t_ik / rowSums(t_ik)
  # A cluster with t_ik = 0 has logmix = -Inf there; 0 * -Inf = NaN, but such
  # cells contribute exactly 0 to the expectation.
  sum(ifelse(t_ik > 0, t_ik * logmix, 0))
}

# ==========================================================================
# Complete-data DIC (Celeux, Forbes, Robert & Titterington, 2006,
#   "Deviance Information Criteria for Missing Data Models",
#    Bayesian Analysis 1(4):651-706)
# --------------------------------------------------------------------------
# The COMPLETE DICs (paper's DIC4/DIC5/DIC6) treat the latent subtype labels Z
# as MISSING DATA and build the deviance from the COMPLETE-DATA likelihood
#
#     f(y, Z | theta) = f(y | Z, theta) * prod_i w_{Z_i}
#
# i.e. the conditional data likelihood f(y|Z,theta) (== `llk.fun` above) TIMES
# the mixture-weight / label term prod_i w_{Z_i}. This is distinct from:
#   * the observed-data (label-marginalised) DIC/BPIC, which integrates the
#     labels OUT (Celeux DIC1-DIC3; llk_obs.fun -- retained as a utility but
#     no longer computed per fit, removed 2026-07); and
#   * the purely CONDITIONAL DIC7/DIC8, which drop the label term.
#
# All three complete DICs share the same first term,
#     -4 * E_{theta,Z}[ log f(y,Z|theta) | y ],
# and differ in the plug-in deviance D(theta-tilde) -- in particular in HOW
# the labels are integrated there:
#   DIC4 : plug in the complete-data posterior mean E_theta[theta | y, Z] and
#          integrate the labels over their POSTERIOR p(Z|y) (the retained
#          draws Z^s). The exact form re-estimates E_theta[theta|y,Z] per
#          label draw (closed form only in the conjugate case, paper p.663);
#          MINDS is a non-conjugate logistic + Gaussian factor model, so we
#          approximate it by the overall posterior-mean estimate `par.hat` /
#          `theta.hat`. DIC4 is the paper's recommended criterion for
#          missing-data models (Conclusion, p.671) -> PRIMARY here.
#   DIC5 : plug in the joint MAP pair (z-hat, theta-hat), approximated by the
#          retained MCMC draw with the highest complete-data joint
#          log-likelihood (diffuse-prior argmax, paper p.656). By this
#          construction pD5 = 2*(max - mean) >= 0 always, BUT the paper finds
#          DIC5 unreliable for mixtures: it treats the n labels as extra
#          parameters, so pD5 scales with Nb and over-penalises (pp.657, 661;
#          Table 3 shows negative/inflated pD5). Diagnostic only.
#   DIC6 : fix theta at the posterior-mean estimate and integrate the labels
#          under their CONDITIONAL distribution given that estimate,
#          E_Z[ log f(y,Z|theta-hat) | y, theta-hat ], computed exactly via
#          the responsibilities t_ik(y, theta-hat) (paper pp.657, 664; see
#          EZ_complete_llk.fun). NB this label integration differs from
#          DIC4's E_Z[ . | y ]: DIC4 and DIC6 are distinct criteria.
#
# Convention:
#     DIC = -4 * E[log f(y,Z|theta)] + 2 * log f(y,Z|theta-tilde)
#         = Dbar + pD,  with Dbar = -2*E[log f], pD = 2*(logf.plug - E[log f]).
# Smaller DIC = better.
# ==========================================================================

# Complete-data joint log-likelihood  log f(y, Z | theta) for one draw.
#   par       : list with x, b, V, U, a_1, a_2, sigma2_y_2 (sigma2_b ignored).
#   Zmat      : Nc x Nb one-hot membership matrix (columns = subjects).
#   w         : length-Nc mixture weights (theta).
#   llk_cond  : the conditional data log-likelihood function (llk.fun), bound to
#               the same data environment (reads y_1, y_2, M_1, M_2, Nb, Nc,
#               Nd1, Nd2 as free variables).
complete_data_llk.fun <- function(par, Zmat, w, llk_cond) {
  par$Z   <- Zmat
  ll_cond <- llk_cond(par)$llk_total               # log f(y | Z, theta)
  nk      <- rowSums(Zmat)                          # subjects per subtype (Nc)
  w       <- pmax(as.numeric(w), 1e-12)             # guard empty components
  ll_cond + sum(nk * log(w))                        # + sum_i log w_{Z_i}
}

# Complete-data DIC4/DIC5/DIC6 from the retained MCMC draws.
#   post.idx        : indices of retained (post-burn-in) draws.
#   Z.all           : Nc x Nb x iter array of sampled one-hot labels.
#   theta.all       : iter x Nc matrix of sampled mixture weights.
#   x/b/V/U/a_1/a_2/sigma2_y_2 .all : stored per-draw parameter arrays.
#   par.hat         : posterior-mean/median point estimate (list; Z overwritten).
#   theta.hat       : posterior point estimate of the mixture weights.
#   llk_cond        : conditional log-likelihood function (llk.fun).
#   llk_EZ          : responsibility-based label-expectation function
#                     (EZ_complete_llk.fun) for the DIC6 plug-in term.
complete_DIC.fun <- function(post.idx, Z.all, theta.all,
                             x.all, b.all, V.all, U.all,
                             a_1.all, a_2.all, sigma2_y_2.all,
                             par.hat, theta.hat, llk_cond, llk_EZ) {

  # per-draw complete-data joint log-likelihood  log f(y, Z^s | theta^s)
  CL <- vapply(post.idx, function(it) {
    par <- list(x = x.all[, , it], b = b.all[, , it],
                V = V.all[, , it], U = U.all[, , it],
                a_1 = a_1.all[, it], a_2 = a_2.all[, it],
                sigma2_y_2 = sigma2_y_2.all[, it])
    complete_data_llk.fun(par, Z.all[, , it], theta.all[it, ], llk_cond)
  }, numeric(1))
  CL  <- CL[is.finite(CL)]
  ECL <- mean(CL)                       # E_{theta,Z}[ log f(y,Z|theta) | y ]

  # DIC5 -- plug in the joint MAP (best retained draw; diffuse-prior approx)
  CL.map <- max(CL)
  dic5   <- -4 * ECL + 2 * CL.map
  pD5    <- -2 * ECL + 2 * CL.map

  # DIC4 (posterior-mean approximation) -- fix theta at the point estimate and
  # integrate the labels over the retained POSTERIOR draws Z^s ~ p(Z|y),
  #   E_Z[ log f(y,Z|theta-hat) | y ] ~= mean_s log f(y, Z^s | theta-hat).
  # (Exact DIC4 would re-estimate E_theta[theta|y,Z^s] per draw; non-conjugate
  # here, so the overall posterior-mean plug-in is used instead.)
  CL.pm <- vapply(post.idx, function(it) {
    complete_data_llk.fun(par.hat, Z.all[, , it], theta.hat, llk_cond)
  }, numeric(1))
  CL.pm <- mean(CL.pm[is.finite(CL.pm)])
  dic4  <- -4 * ECL + 2 * CL.pm
  pD4   <- -2 * ECL + 2 * CL.pm

  # DIC6 -- fix theta at the point estimate and integrate the labels under
  # their CONDITIONAL distribution given that estimate (responsibilities),
  #   E_Z[ log f(y,Z|theta-hat) | y, theta-hat ]  (exact, no MC error).
  CL.ez <- llk_EZ(par.hat, theta.hat)
  dic6  <- -4 * ECL + 2 * CL.ez
  pD6   <- -2 * ECL + 2 * CL.ez

  list(E_complete_loglik = ECL,
       Dbar_complete     = -2 * ECL,
       dic4 = dic4, pD4 = pD4,
       dic5 = dic5, pD5 = pD5,
       dic6 = dic6, pD6 = pD6)
}
