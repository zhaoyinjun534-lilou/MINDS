\name{generate_data_mixed}
\alias{generate_data_mixed}
\title{Generate Simulated Mixed-Modality Dataset}
\usage{
generate_data_mixed(seed.no)
}
\arguments{
  \item{seed.no}{Integer random seed for data generation.}
}
\value{
A named list with:
\itemize{
\item \code{y_1}: Binary outcome matrix.
\item \code{y_2}: Continuous outcome matrix.
}
}
\description{
Generates the simulated dataset used by MINDS, including binary outcomes
(\code{y_1}) and continuous outcomes (\code{y_2}), based on the simulation
parameters and latent variables defined in the package's data generation code.
}
\examples{
sim_data <- generate_data_mixed(103)
str(sim_data)
}
