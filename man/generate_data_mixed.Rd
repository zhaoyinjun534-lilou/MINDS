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
