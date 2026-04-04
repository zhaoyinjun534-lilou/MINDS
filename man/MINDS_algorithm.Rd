\name{MINDS_algorithm}
\alias{MINDS_algorithm}
\title{Fit MINDS Algorithm For Mixed Binary-Continuous Data}
\usage{
MINDS_algorithm(
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
  plot_trace = TRUE
)
}
\arguments{
  \item{y_1}{Binary outcome matrix.}
  \item{y_2}{Continuous outcome matrix.}
  \item{Nc}{Number of latent clusters.}
  \item{Nt}{Number of latent constructs.}
  \item{iter.max}{Number of MCMC iterations.}
  \item{mu_x}{Prior mean for cluster centers.}
  \item{sigma2_x}{Prior variance for cluster centers.}
  \item{IG_b.shape}{Inverse-gamma shape for sigma2_b.}
  \item{IG_b.scale}{Inverse-gamma scale for sigma2_b.}
  \item{p.theta.prior}{Dirichlet prior for cluster weights.}
  \item{IG_y_2.shape}{Inverse-gamma shape for sigma2_y_2.}
  \item{IG_y_2.scale}{Inverse-gamma scale for sigma2_y_2.}
  \item{sigma2_a_1}{Prior variance for binary intercepts.}
  \item{sigma2_a_2}{Prior variance for continuous intercepts.}
  \item{sigma2_v}{Prior variance for binary loadings.}
  \item{mu_v}{Prior mean for binary loadings.}
  \item{sigma2_u}{Prior variance for continuous loadings.}
  \item{mu_u}{Prior mean for continuous loadings.}
  \item{init_seed}{Seed for initialization and imputation.}
  \item{plot_trace}{If TRUE, plots likelihood trace.}
}
\value{
A named list containing membership, cluster center, modality loadings,
intercepts, membership weights, likelihood trace, and ic.
}
\description{
Runs the Bayesian MINDS model for mixed binary-continuous outcomes,
designed for integrative subtyping with multimodal data.
}
\details{
MINDS (Mixed INtegrative Data Subtyping) jointly models binary and continuous
modalities in a unified Bayesian hierarchical framework. The method performs
simultaneous latent structure learning and clustering, and uses
P\'olya-Gamma augmentation for efficient posterior computation.

The function returns subtype membership, latent cluster centers, modality-specific
loadings/intercepts, a likelihood trace summary, and an information criterion
(`ic`) derived from the DIC helper.
}
\references{
Zhao Y, Wang Y, Liu Y. An Integrative Approach for Subtyping Mental Disorders
Using Multimodal Data. (Manuscript/PDF source used for package manual update.)
}
\examples{
data(data_mixed)
out <- MINDS_algorithm(
  y_1 = data_mixed$y_1,
  y_2 = data_mixed$y_2,
  Nc = 4,
  Nt = 3,
  iter.max = 2,
  plot_trace = FALSE
)
names(out)
}
