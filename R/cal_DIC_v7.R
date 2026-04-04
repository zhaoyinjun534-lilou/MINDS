
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
  temp1 <- t(apply(temp, 1, function (t) {t - a_1}));#(Z^T*x+b)V-a
  p <- exp(temp1)/(1 + exp(temp1)) 

  aa1 <- sum(apply(M_1*y_1*log(p)+log(1-p)*M_1*(1-y_1), 1, 
                   function(a) sum(a, na.rm = TRUE)), 
             na.rm = TRUE)#binary
  
  
  temp <- (t(Z) %*% x + b) %*% U; #(Z^T*x+b)*U, Nb*Nd2
  temp1 <- t(apply(temp, 1, function (t) {t - a_2}));#(Z^T*x+b)U-a
  
  logd <- sapply(1:Nd2, function(j) {
    M_2[, j] * dnorm(y_2[, j], mean = temp1[, j], sd = sqrt(sigma2_y_2[j]), log = TRUE)
  })
  aa2 <- sum(apply(logd,2,function(a) sum(a, na.rm = TRUE)), na.rm = TRUE)
  aa <- aa1 + aa2
  return(list('llk_total'=aa, 'llk_binary'=aa1, 'llk_continuous'=aa2))
}

par.simu.fun <- function(par.est){
  n.rep <- 20
  
  Z <- par.est$Z
  x <- par.est$x
  b <- par.est$b
  V <- par.est$V
  U <- par.est$U
  a_1 <- par.est$a_1
  a_2 <- par.est$a_2
  sigma2_y_2 <- par.est$sigma2_y_2
  
  # sample polya-gamma distributed w
  
  temp1 <- (t(Z) %*% x + b) %*% V; #(Z^T*x+b)*V
  temp <- t(apply(temp1, 1, function (t) {t-a_1}));#(Z^T*x+b)*V-a_1
  w <- apply(temp, c(1,2), function(s) (rpg(num=1, h=1, z=s))) #random sample from polya-gamma distribution
  colnames(w) <- NULL
  
  ####obtain conditional posterior for a
  
  #a_1----
  V_a <- 1/(apply(w*M_1, 2, sum) + 1/sigma2_a_1)
  
  temp <- (t(Z) %*% x + b) %*% V; #(Z^T*x+b)*V, Nb*Nd
  r <- temp - (y_1-1/2)/w

  Mu_a <- apply(w*M_1*r, 2, function(a) sum(a, na.rm = T)) *V_a
  
  
  a_1.rep <- replicate(n.rep, sapply(seq(Nd1), function(i) rnorm(1, Mu_a[i], sqrt(V_a[i]))))
  
  #a_2----
  V_a <- 1/(apply(M_2 * rep(1/sigma2_y_2, each = nrow(M_2)), 2, sum) + 1/sigma2_a_2)
  
  temp <- (t(Z) %*% x + b) %*% U; #(Z^T*x+b)*V, Nb*Nd
  r <- temp - y_2

  Mu_a <- apply(r*M_2, 2, function(a) sum(a, na.rm = T))/sigma2_y_2 *V_a
  
  
  a_2.rep <- replicate(n.rep, sapply(seq(Nd2), function(i) rnorm(1, Mu_a[i], sqrt(V_a[i]))))
  
  #V_ti----
  V.rep.temp <- NULL
  #obtain conditional posterior for V_ti, loading of the ti th latent construct
  for (ti in 1:Nt) {
    
    index_x <- rep(0, Nt)
    index_x[ti] <- 1
    
    temp_A <- (t(Z) %*% x + b) %*% index_x
    temp <- apply(w*M_1, 2, function(s) sum(s * temp_A^2)) #allow missingness
    
    B_V_ti <- 1/(temp + 1/sigma2_v)
    
    x.sub <- x
    x.sub[, ti] <- 0
    b.sub <- b
    b.sub[,ti] <- 0
    
    temp1 <- (t(Z) %*% x.sub + b.sub) %*% V
    temp2 <- t(apply(temp1, 1, function (t) {t - a_1}));#(Z^T*x+b)V-a
    phi <- (y_1-1/2)/w - temp2
    
    
    temp3 <- apply(w*M_1, 2, function(s) s * temp_A)
    
    D_V_ti <- apply(temp3 *phi, 2, function(a) sum(a, na.rm = T)) + mu_v/sigma2_v
    
    
    Mu_V_ti <- B_V_ti * D_V_ti
    
    V_ti.rep <- replicate(n.rep, mvrnorm(1, Mu_V_ti, diag(B_V_ti)))
    
    V.rep.temp <- append(V.rep.temp, list(V_ti.rep))
  }
  

  temp1 <- do.call(rbind, V.rep.temp)  
  V.rep <- lapply(seq(1:n.rep), function(i) {matrix(temp1[,i], byrow=T, nrow=Nt)})
  
  #X---- 
  x.rep.temp <- NULL
  ####obtain conditional posterior for X_ti, cluster center for ti th latent construct
  for (ti in 1:Nt) {
    
    index_x <- rep(0, Nt)
    index_x[ti] <- 1
    
    x.sub <- x
    x.sub[, ti] <- 0
    
    B_ti <- diag(sigma2_x, Nc)
    
    #binary part
    
    temp <- (t(Z) %*% x.sub + b) %*% V; #c(temp1)# convert 'temp1' from a matrix into a vector
    temp1 <- t(apply(temp, 1, function (t) {t-a_1}));#temp #a needs to be a vector; (Z^T*x+b)V-a
    
    phi_1 <- (y_1-1/2)/w - temp1
    phi_1 <- ifelse(is.na(phi_1), 0, phi_1)
    
    colnames(w) <- NULL
    colnames(phi_1) <- NULL
    
    dt <- abind('z'=Z, 'w'=t(w), 'phi'=t(phi_1), along = 1) #combine z, w and phi

    aa <- lapply(seq(Nb), function(i) {
      
      Z_i <- dt[startsWith(rownames(dt),'z'), i]
      w_i <- dt[startsWith(rownames(dt),'w'), i]
      phi_i <- dt[startsWith(rownames(dt),'phi'), i]
      
      temp <- Z_i %*% t(V[ti, ])
      
      return (list('v_temp'=temp %*%  diag(w_i*  M_1[i,]) %*% t(temp),#allows missingness
                   'mu_temp' = temp %*%  diag(w_i*  M_1[i,]) %*% phi_i
      ))}
    )
    
    v_sum_1 <- Reduce("+", lapply(seq(ncol(dt)),
                                  function(i) aa[[i]]$v_temp))
    mu_sum_1 <- Reduce("+", lapply(seq(ncol(dt)),
                                   function(i) aa[[i]]$mu_temp))
    
    # #continuous part
    temp <- (t(Z) %*% x.sub + b) %*% U; #c(temp1)# convert 'temp1' from a matrix into a vector
    temp1 <- t(apply(temp, 1, function (t) {t-a_2}));#temp #a needs to be a vector; (Z^T*x+b)V-a
    
    Phi_1 <- y_2- temp1
    Phi_1 <- ifelse(is.na(Phi_1), 0, Phi_1)
    
    
    dt <- abind('z'=Z, 'Phi'=t(Phi_1), along = 1) #combine z and phi
    
    aa <- lapply(seq(Nb), function(i) {
      
      Z_i <- dt[startsWith(rownames(dt),'z'), i]
      Phi_i <- dt[startsWith(rownames(dt),'Phi'), i]
      
      temp <- Z_i %*% t(U[ti, ])
      
      return (list('v_temp'=temp %*%  diag(1/sigma2_y_2*M_2[i,]) %*% t(temp),
                   'mu_temp' = temp %*%  diag(1/sigma2_y_2*M_2[i,]) %*% Phi_i
      ))}
    )
    
 
    v_sum_2 <- Reduce("+", lapply(seq(ncol(dt)),
                                  function(i) aa[[i]]$v_temp))
    mu_sum_2 <- Reduce("+", lapply(seq(ncol(dt)),
                                   function(i) aa[[i]]$mu_temp))
    
    V_x_ti_inverse <- solve(B_ti) + v_sum_1 + v_sum_2
    V_x_ti <- solve(V_x_ti_inverse);
    Mu_x_ti <- V_x_ti %*% (mu_sum_1 + mu_sum_2 + solve(B_ti) %*% rep(mu_x, Nc));
    

    x_ti.rep <- replicate(n.rep, mvrnorm(1, Mu_x_ti, V_x_ti))
    
    x.rep.temp <- append(x.rep.temp, list(x_ti.rep))
  }
  
  temp1 <- do.call(rbind, x.rep.temp)  
  x.rep <- lapply(seq(1:n.rep), function(i) {matrix(temp1[,i], byrow=F, nrow=Nc)})
  
  #b_ti ----
  b.rep.temp <- NULL
  for (ti in 1:Nt) {
    
    b.sub <- b
    b.sub[,ti] <- 0
    
    #binary part
    temp <- (t(Z) %*% x + b.sub) %*% V
    temp1 <- t(apply(temp, 1, function (t) {t-a_1}))# (Z^T*x+b)V-a
    C_ti <- (y_1-1/2)/w - temp1
    
    #continuous part
    temp <- (t(Z) %*% x + b.sub) %*% U
    temp1 <- t(apply(temp, 1, function (t) {t-a_2}))# (Z^T*x+b)V-a
    D_ti <- y_2 - temp1
    
    V_b_ti_diag <- 1/(apply(w*M_1,1, function(w_i) {sum(w_i * V[ti,]^2)}) +
                        apply(M_2,1, function(M2_i) {sum(M2_i * 1/sigma2_y_2*U[ti,]^2)}) +
                        1/sigma2_b[ti]) #diagonal elements of cov(b)
    
    Mu_b_ti <- (apply(w*M_1*C_ti, 1, function(t) {sum(t * V[ti,], na.rm = T)}) +
                  apply(D_ti*M_2, 1, function(t) 
                  {sum(t * U[ti,] * 1/sigma2_y_2, na.rm = T)}))* (V_b_ti_diag)
    
    
    b_ti.rep <- replicate(n.rep,
                          sapply(seq(Nb), 
                                 function(bi) rnorm(1, Mu_b_ti[bi], sqrt(V_b_ti_diag[bi]))))
    
    b.rep.temp <- append(b.rep.temp, list(b_ti.rep))
  }
  
  temp1 <- do.call(rbind, b.rep.temp)  
  b.rep <- lapply(seq(1:n.rep), function(i) {matrix(temp1[,i], byrow=F, nrow=Nb)})
  
  #U----
  #obtain conditional posterior for U_ti, loading of ti th latent construct
  U.rep.temp<- NULL
  for (ti in 1:Nt) {
    
    index_x <- rep(0, Nt)
    index_x[ti] <- 1
    
    temp_A <- (t(Z) %*% x + b) %*% index_x
    
    temp_B <- apply(M_2, 2, function(row) row %*% temp_A^2)
    
    B_U_ti <- 1/(temp_B /sigma2_y_2 + 1/sigma2_u)
    
    
    x.sub <- x
    x.sub[, ti] <- 0
    b.sub <- b
    b.sub[,ti] <- 0
    
    temp1 <- (t(Z) %*% x.sub + b.sub) %*% U
    temp2 <- t(apply(temp1, 1, function (t) {t - a_2}));#(Z^T*x+b)V-a
    phi <- y_2 - temp2
    
    temp3 <- apply(phi*M_2, 2, function(s) sum(s * temp_A, na.rm = T))
    
    D_U_ti <- temp3/sigma2_y_2 + mu_u/sigma2_u
    
    Mu_U_ti <- B_U_ti * D_U_ti
    U_ti.rep <- replicate(n.rep, mvrnorm(1, Mu_U_ti, diag(B_U_ti)))
    
    U.rep.temp <- append(U.rep.temp, list(U_ti.rep))
  }
  
  
  temp1 <- do.call(rbind, U.rep.temp)  
  U.rep <- lapply(seq(1:n.rep), function(i) {matrix(temp1[,i], byrow=T, nrow=Nt)})

  
  # Z----
  alpha <- NULL
  for (k in 1:Nc){
    
    # binary items
    temp1 <- t(apply(b, 1, function(b_i) {x[k, ] + b_i}) ) %*% V
    temp2 <- t(apply(temp1, 1, function (t) {t-a_1})) - (y_1-1/2)/w
    temp2 <- ifelse(is.na(temp2), 0, temp2)
    
    p1_k <- exp(-w*M_1/2 * temp2^2)
    
    # continuous items
    temp1 <- t(apply(b, 1, function(b_i) {x[k, ] + b_i}) ) %*% U
    temp2 <- t(apply(temp1, 1, function (t) {t-a_2}))
    sigma2_y_2_matrix <- matrix(rep(sigma2_y_2, Nb), nrow = Nb, byrow = T)
    

    aa <- temp2 - y_2
    aa <- ifelse(is.na(aa), 0, aa)
    
    
    p2_k <- (1/sqrt(2*pi*sigma2_y_2_matrix))^M_2* exp(-(aa)^2/2/sigma2_y_2_matrix)
    
    # integrated
    alpha.k <- apply(p1_k, 1, prod) *apply(p2_k, 1, prod)* theta[k]
    alpha <- cbind(alpha, alpha.k)    
  }
  
  alpha <- alpha/apply(alpha, 1, sum); #alpha #standardize
  
  Z.rep <- replicate(n.rep, apply(alpha, 1, function(pi) rmultinom(1, 1, pi)))
  
  
  #obtain conditional posterior for theta
  theta.rep <- replicate(n.rep, rdirichlet(1, apply(Z, 1, sum) + p.theta.prior))
  
  #obtain conditional posterior for sigma2_y_2
  
  temp <- (t(Z) %*% x + b) %*% U; #(Z^T*x+b)*U, Nb*Nd2
  temp1 <- t(apply(temp, 1, function (t) {t - a_2}));#(Z^T*x+b)U-a
  C <- y_2 - temp1

  aa <- apply(C^2*M_2, 2, function(a) sum(a, na.rm = T))/2
  
  
  sigma2_y_2.rep <- replicate(n.rep, unlist(lapply(seq(Nd2), function(j) {
    rinvgamma(1, Nb/2 + IG_y_2.shape, aa[j] + IG_y_2.scale)})))
  
  par.rep <- list('Z'=Z.rep, 'x'=x.rep, 'b'=b.rep, 'V'=V.rep, 'U'=U.rep, 'a_1'=a_1.rep, 'a_2'=a_2.rep, 'sigma2_y_2'=sigma2_y_2.rep)
  return(par.rep)
}


dic.fun <- function(par.est){
  n.rep <- 20
  par.rep <- par.simu.fun(par.est)
  
  llk.mean <- mean(sapply(1:n.rep, function(i){
    par.s <- list('Z'= par.rep$Z[,,i], 'x'= par.rep$x[[i]], 'b'= par.rep$b[[i]], 
                  'V'= par.rep$V[[i]], 
                  'U'=par.rep$U[[i]], 'a_1'= par.rep$a_1[,i], 'a_2'= par.rep$a_2[,i],
                  'sigma2_y_2'= par.rep$sigma2_y_2[,i]) 
    llk.fun(par.s)$llk_total
  })    )
  #add penalty to parameters
  ic <- -6*llk.mean + 4*llk.fun(par.est)$llk_total
  dic <- -4*llk.mean + 2*llk.fun(par.est)$llk_total
  return(list("dic"=dic, "ic"=ic)) #smaller, better
}







