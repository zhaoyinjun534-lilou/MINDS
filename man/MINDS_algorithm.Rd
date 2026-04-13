\name{MINDS_algorithm}
\alias{MINDS_algorithm}
\title{Fit MINDS Algorithm For Mixed Binary-Continuous Data}
\usage{
MINDS_algorithm(
  y_1,
  y_2,
  Nc,
  Nt,
  iter.max = 2000,
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
  plot_trace = TRUE,
  use_true_init = FALSE,
  use_true_Z_init = FALSE,
  fix_x = FALSE,
  fix_Z = FALSE,
  fix_others_at_true = FALSE
)
}
\arguments{
  \item{y_1}{Binary outcome matrix with subjects in rows and items in columns.}
  \item{y_2}{Continuous outcome matrix with subjects in rows and items in columns.}
  \item{Nc}{Integer number of latent clusters.}
  \item{Nt}{Integer number of latent constructs.}
  \item{iter.max}{Integer number of MCMC iterations.}
  \item{mu_x}{Prior mean for cluster centers \code{x}.}
  \item{sigma2_x}{Prior variance for cluster centers \code{x}.}
  \item{IG_b.shape}{Inverse-gamma shape for random-effect variance \code{sigma2_b}.}
  \item{IG_b.scale}{Inverse-gamma scale for random-effect variance \code{sigma2_b}.}
  \item{p.theta.prior}{Dirichlet prior parameters for cluster weights.}
  \item{IG_y_2.shape}{Inverse-gamma shape for \code{sigma2_y_2}.}
  \item{IG_y_2.scale}{Inverse-gamma scale for \code{sigma2_y_2}.}
  \item{sigma2_a_1}{Prior variance for binary intercepts.}
  \item{sigma2_a_2}{Prior variance for continuous intercepts.}
  \item{sigma2_v}{Prior variance for binary loadings.}
  \item{mu_v}{Prior mean for binary loadings.}
  \item{sigma2_u}{Prior variance for continuous loadings.}
  \item{mu_u}{Prior mean for continuous loadings.}
  \item{init_seed}{Seed used for initialization and imputation.}
  \item{plot_trace}{Logical; if \code{TRUE}, draw the likelihood trace.}
  \item{use_true_init}{Logical; if \code{TRUE}, initialize parameters from available \code{*.true} objects.}
  \item{use_true_Z_init}{Logical; if \code{TRUE}, initialize \code{Z} from \code{Z.true} while leaving other initials unchanged.}
  \item{fix_x}{Logical; if \code{TRUE}, keep cluster centers \code{x} fixed to initialized values.}
  \item{fix_Z}{Logical; if \code{TRUE}, keep cluster assignments \code{Z} fixed to initialized values.}
  \item{fix_others_at_true}{Logical; if \code{TRUE}, keep all non-\code{x} parameters fixed at available true values and update only \code{x}.}
}
\value{
A named list with elements:
\itemize{
\item \code{membership}: Estimated cluster assignment for each subject.
\item \code{cluster center}: Estimated cluster centers.
\item \code{x trace}: MCMC trace array for \code{x} across iterations.
\item \code{loading to binary modality}: Estimated binary loadings.
\item \code{loading to continuous modality}: Estimated continuous loadings.
\item \code{binary modality intercept}: Estimated binary intercepts.
\item \code{continuous modality intercept}: Estimated continuous intercepts.
\item \code{membership weight}: Estimated cluster weights.
\item \code{likelihood trace plot}: Numeric likelihood trace.
\item \code{ic}: Information criterion from DIC helper.
}
}
\description{
Runs the Bayesian MINDS model using binary outcomes \code{y_1} and continuous
outcomes \code{y_2}, with user-specified prior hyperparameters.
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
  iter.max = 50,
  plot_trace = FALSE
)
names(out)
}
