Nd1 <- 10
Nd2 <- 12
iter.max <- 2000
seed.no <- 103
Nb <- 1000
ind.simu <- 6003
n.rep <-2
sigma2_b.true.element <- 0.1
Nc <- 4
Nt <- 3


#true values of the parameters
sigma2_b.true <- rep(sigma2_b.true.element, Nt) #var of true b
theta.true <- rep(1 / Nc, Nc)

if(Nc==4){
  theta.true <- c(0.3, 0.3, 0.2, 0.2) #Nc
}



sigma2_y_2.true <- rep(1, Nd2)


par.fun <- function(){
  
  set.seed(4)
  
  # Item-specific fixed effects (intercepts)
  a_1.true <- rnorm(Nd1, mean = 0, sd = 0.5)
  a_2.true <- rnorm(Nd2, mean = 0, sd = 0.5)

  # Cluster centers in latent space
  # x.true <- matrix(rnorm(Nc * Nt, mean = 0, sd = 0.2), nrow = Nc, ncol = Nt)
  # x.true <- x.true + (seq_len(Nc) - 1) * 1.2
  
  x.true <- replicate(Nt,
                      runif(Nc, 0, 2)) # for Nt latent constructs
  x.true <- apply(x.true, 2, scale)*2
  
 

  # Modality-specific linear maps from latent space to item space
  V.raw.item <-  (replicate(Nt, runif(Nd1, 0, 2)))
  U.raw.item <-  (replicate(Nt, runif(Nd2, 0, 2)))
  
  
  # V.raw.item <- matrix(rnorm(Nd1 * Nt, mean = 0, sd = 0.5), nrow = Nd1, ncol = Nt)
  # U.raw.item <- matrix(rnorm(Nd2 * Nt, mean = 0, sd = 0.5), nrow = Nd2, ncol = Nt)
  V.raw <- t(V.raw.item)
  U.raw <- t(U.raw.item)
  
  
  # 1. QR decomposition with positive diagonal on R
  L.raw <- rbind(t(V.raw), t(U.raw))          # (p1+p2) x d
  qr.L  <- qr(L.raw)
  Q     <- qr.Q(qr.L)[, seq_len(Nt), drop = FALSE]
  R     <- qr.R(qr.L)[seq_len(Nt), seq_len(Nt), drop = FALSE]
  # Ensure strictly positive diagonal of R
  sign.d <- sign(diag(R))
  Q      <- sweep(Q, 2, sign.d, "*")
  R      <- diag(sign.d) %*% R
  # Reparameterize: X <- X R^T, L <- Q
  x.true <- x.true %*% t(R)
  V.true <- t(Q[seq_len(Nd1),       , drop = FALSE])
  U.true <- t(Q[Nd1 + seq_len(Nd2), , drop = FALSE])

  # Fix dimension identifiability by ordering latent dimensions using
  # L2 norms of x.true across clusters.
  x_l2 <- apply(x.true, 2, function(col_j) sqrt(sum(col_j^2)))
  dim_order <- order(x_l2, decreasing = TRUE)
  x.true <- x.true[, dim_order, drop = FALSE]
  V.true <- V.true[dim_order, , drop = FALSE]
  U.true <- U.true[dim_order, , drop = FALSE]
  R <- R[dim_order, dim_order, drop = FALSE]
  
  # 2. Location constraint: sum_k omega_k X_k = 0
  x.shift <- apply(x.true, 2, mean)
  x.true  <- sweep(x.true, 2, x.shift, FUN = "-")
  a_1.true <- a_1.true - as.numeric(x.shift %*% V.true)
  a_2.true <- a_2.true - as.numeric(x.shift %*% U.true)
  
  
  
  return(list(
    'x.true'     = x.true,
    'V.true'     = V.true,
    'U.true'     = U.true,
    'a_1.true'   = a_1.true,
    'a_2.true'   = a_2.true,
    'R' = R))
}

raw.par    <- par.fun()
x.true     <- raw.par$x.true
V.true     <- raw.par$V.true
U.true     <- raw.par$U.true
a_1.true   <- raw.par$a_1.true
a_2.true   <- raw.par$a_2.true
R <- raw.par$R

par.fun2<- function(){
  set.seed(4)
  Z.true <- rmultinom(Nb, 1, theta.true); #Nc*Nb
  b.true <- MASS::mvrnorm(Nb, rep(0, Nt), Sigma = diag(sigma2_b.true))
  b.true <- b.true %*% t(R)
  return(list("Z.true"=Z.true, "b.true"=b.true))
}
Z.true <- par.fun2()$Z.true
b.true <- par.fun2()$b.true

Z.cat.true <- apply(Z.true, 2, function (a) which(a==1))




#obtain data, this needs seed.no
y.fun <- function(seed.no){
  set.seed(seed.no)
  temp <- (t(Z.true) %*% x.true + b.true) %*% V.true; #(Z^T*x+b)*V, Nb*Nd1
  temp1 <- t(apply(temp, 1, function (t) {t - a_1.true}));#(Z^T*x+b)V-a
  p <- exp(temp1)/(1 + exp(temp1))
  
  y_1 <- apply(p, c(1,2), function(a) (rbinom(1, 1, a)))
  
  temp <- (t(Z.true) %*% x.true + b.true) %*% U.true; #(Z^T*x+b)*U, Nb*Nd2
  temp1 <- t(apply(temp, 1, function (t) {t - a_2.true}));#(Z^T*x+b)U-a
  
  y_2 <- t(sapply(seq(Nb), 
                  function(i) MASS::mvrnorm(1, temp1[i,], diag(sigma2_y_2.true))))
  return(list('y_1'= y_1, 'y_2' = y_2))
}


# Generate the observed mixed outcomes used by downstream code.
#' Generate Simulated Mixed-Modality Dataset
#'
#' Generates the simulated dataset used by MINDS, including binary outcomes
#' (`y_1`) and continuous outcomes (`y_2`), based on the globally defined
#' simulation parameters and latent variables in this script.
#'
#' @param seed.no Integer random seed for data generation.
#'
#' @return A named list with:
#' \itemize{
#'   \item `y_1`: Binary outcome matrix.
#'   \item `y_2`: Continuous outcome matrix.
#' }
#'
#' @examples
#' \dontrun{
#' source("MINDS_Rpackage/data_generation.R")
#' sim_data <- generate_data_mixed(103)
#' str(sim_data)
#' }
#'
#' @export
generate_data_mixed <- function(seed.no){
  y.sim <- y.fun(seed.no)
  data_mixed <- list('y_1' = y.sim$y_1, 'y_2' = y.sim$y_2)
  return(data_mixed)
}



data_mixed <- generate_data_mixed(seed.no)
y_1 <- data_mixed$y_1
y_2 <- data_mixed$y_2