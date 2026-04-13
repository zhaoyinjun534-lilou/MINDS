\name{generate_data_mixed}
\alias{generate_data_mixed}
\title{Generate Simulated Mixed-Modality Dataset}
\usage{
generate_data_mixed(
  seed.no = 103,
  Nb = 1000,
  Nc = 4,
  Nd1 = 10,
  Nd2 = 12,
  Nt = 3,
  sigma2_b_element = 0.1,
  theta = NULL
)
}
\arguments{
  \item{seed.no}{Integer random seed for data generation. Default \code{103}.}
  \item{Nb}{Number of subjects. Default \code{1000}.}
  \item{Nc}{Number of latent clusters. Default \code{4}.}
  \item{Nd1}{Number of binary items. Default \code{10}.}
  \item{Nd2}{Number of continuous items. Default \code{12}.}
  \item{Nt}{Number of latent constructs. Default \code{3}.}
  \item{sigma2_b_element}{Scalar variance of the random effects. Default \code{0.1}.}
  \item{theta}{Numeric vector of length \code{Nc} giving cluster mixture weights.
    If \code{NULL} (default), uses \code{c(0.3, 0.3, 0.2, 0.2)} when \code{Nc == 4},
    otherwise equal weights \code{rep(1/Nc, Nc)}.}
}
\value{
A named list with:
\itemize{
\item \code{y_1}: Binary outcome matrix (subjects x binary items).
\item \code{y_2}: Continuous outcome matrix (subjects x continuous items).
}
}
\description{
Generates a simulated dataset for the MINDS model, including binary outcomes
(\code{y_1}) and continuous outcomes (\code{y_2}).
}
\examples{
sim_data <- generate_data_mixed(seed.no = 103)
str(sim_data)
}
