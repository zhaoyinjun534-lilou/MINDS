\name{MINDS_algorithm}
\alias{MINDS_algorithm}
\title{Fit MINDS Algorithm For Mixed Binary-Continuous Data}
\usage{
MINDS_algorithm(
  y_1,
  y_2,
  Nc,
  Nt,
  iter.max = 1000,
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
  init_seed = 123,
  plot_trace = TRUE,
  n.rep = 20
)
}
\arguments{
  \item{y_1}{Binary outcome matrix with subjects in rows and items in columns.}
  \item{y_2}{Continuous outcome matrix with subjects in rows and items in columns.}
  \item{Nc}{Integer number of latent clusters.}
  \item{Nt}{Integer number of latent constructs.}
  \item{iter.max}{Integer number of MCMC iterations. Default \code{1000}.}
  \item{mu_x}{Prior mean for cluster centers \code{x}. Default \code{0}.}
  \item{sigma2_x}{Prior variance for cluster centers \code{x}. Default \code{1}.}
  \item{IG_b.shape}{Inverse-gamma shape for random-effect variance \code{sigma2_b}. Default \code{11}.}
  \item{IG_b.scale}{Inverse-gamma scale for random-effect variance \code{sigma2_b}. Default \code{1}.}
  \item{p.theta.prior}{Dirichlet prior parameters for cluster weights. Default \code{rep(1, Nc)}.}
  \item{IG_y_2.shape}{Inverse-gamma shape for \code{sigma2_y_2}. Default \code{2}.}
  \item{IG_y_2.scale}{Inverse-gamma scale for \code{sigma2_y_2}. Default \code{1}.}
  \item{sigma2_a_1}{Prior variance for binary intercepts. Default \code{10}.}
  \item{sigma2_a_2}{Prior variance for continuous intercepts. Default \code{10}.}
  \item{sigma2_v}{Prior variance for binary loadings. Default \code{10}.}
  \item{mu_v}{Prior mean for binary loadings. Default \code{0}.}
  \item{sigma2_u}{Prior variance for continuous loadings. Default \code{10}.}
  \item{mu_u}{Prior mean for continuous loadings. Default \code{0}.}
  \item{true.values}{Optional list of true parameter values for initialization. 
    Expected list components are \code{x.true}, \code{V.true}, \code{U.true}, 
    \code{a_1.true}, and \code{a_2.true}. Default \code{NULL}.}
  \item{init_seed}{Seed used for initialization and imputation. Default \code{123}.}
  \item{plot_trace}{Logical; if \code{TRUE}, draw the likelihood trace. Default \code{TRUE}.}
  \item{n.rep}{Integer number of posterior draws used by DIC computation. Default \code{20}.}
}
\value{
A named list with MCMC posterior samples and summary statistics.
}
\description{
Runs the Bayesian MINDS model using binary outcomes \code{y_1} and continuous
outcomes \code{y_2}, with user-specified prior hyperparameters.
}
\details{
MINDS (Mixed INtegrative Data Subtyping) jointly models binary and continuous
modalities in a unified Bayesian hierarchical framework using Pólya-Gamma 
augmentation for efficient posterior computation.
}
\examples{\dontrun{
data(data_mixed)
out <- MINDS_algorithm(
  y_1 = data_mixed$y_1,
  y_2 = data_mixed$y_2,
  Nc = 4,
  Nt = 3,
  iter.max = 50,
  plot_trace = FALSE
)
names(out)
}
}
