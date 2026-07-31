#' Example Simulated Mixed-Modality Dataset
#'
#' Simulated data generated from the MINDS data generation mechanism, produced
#' by `generate_data_mixed(seed.no = 103)` at the package defaults (1000
#' subjects, 10 binary items, 12 continuous items, 5 clusters, 3 latent
#' constructs). The true parameter values are included so a fit can be scored
#' against the truth. See `data-raw/make_data_mixed.R`.
#'
#' @format A list with the outcome matrices, the true parameters, and the
#' dimensions of the setting:
#' \describe{
#'   \item{y_1}{Binary matrix, subjects in rows and binary items in columns.}
#'   \item{y_2}{Continuous matrix, subjects in rows and continuous items in columns.}
#'   \item{x.true}{True cluster centres (clusters x latent constructs).}
#'   \item{V.true, U.true}{True loadings for the binary / continuous items
#'     (latent constructs x items).}
#'   \item{a_1.true, a_2.true}{True binary / continuous item intercepts.}
#'   \item{Z.true}{True cluster assignments (clusters x subjects, one-hot).}
#'   \item{b.true}{True subject random effects (subjects x latent constructs).}
#'   \item{theta.true}{True cluster weights.}
#'   \item{sigma2_y_2.true}{True residual variances for the continuous items.}
#'   \item{sigma2_b.true}{True latent variances, one per construct.}
#'   \item{Nc, Nt, Nb}{Numbers of clusters, latent constructs, and subjects.}
#' }
#' @seealso [generate_data_mixed()]
#' @usage data(data_mixed)
"data_mixed"
