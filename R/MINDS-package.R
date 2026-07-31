#' MINDS: Bayesian Integrative Clustering For Mixed-Type Multimodal Data
#'
#' Fits the MINDS model, a Bayesian hierarchical mixture that clusters subjects
#' jointly on a binary and a continuous modality via a shared low-dimensional
#' latent construct space, using Polya-Gamma augmentation for the binary part.
#' Model selection over the number of clusters `Nc` and constructs `Nt` uses the
#' complete-data DIC family of Celeux et al. (2006).
#'
#' @keywords internal
#'
#' @importFrom BayesLogit rpg
#' @importFrom EDISON rinvgamma
#' @importFrom LaplacesDemon rdirichlet
#' @importFrom MASS mvrnorm
#' @importFrom abind abind
#' @importFrom kamila kamila
#' @importFrom mice mice complete
#' @importFrom stats dnorm median quantile rbinom rmultinom rnorm runif var
#' @importFrom graphics plot
#' @importFrom utils globalVariables
"_PACKAGE"
