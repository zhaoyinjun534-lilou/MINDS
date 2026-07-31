# MINDS

MINDS is an R package for Bayesian integrative clustering and subtyping using mixed-type multimodal data.

The method is motivated by psychiatric disorder subtyping where heterogeneous sources (for example, symptom indicators, cognitive measures, and imaging-derived features) are analyzed jointly. The MINDS framework performs:
- joint modeling of binary and continuous modalities,
- latent structure estimation,
- subtype clustering,
- and uncertainty-aware inference via a Bayesian hierarchical model with Pólya-Gamma augmentation.

## Installation

From the project root:

```r
# install.packages("devtools")
devtools::install("MINDS")
```

Or from a shell:

```bash
R CMD INSTALL MINDS
```

Dependencies are all on CRAN: `BayesLogit`, `EDISON`, `LaplacesDemon`, `MASS`, `abind`, `kamila`, `mice`.

## Included data

The package includes a simulated dataset:

```r
library(MINDS)
data(data_mixed)
str(data_mixed)
```

`data_mixed` is `generate_data_mixed(seed.no = 103)` at its defaults — 1000 subjects, 10 binary items, 12 continuous items, 5 clusters, 3 latent constructs. Besides the outcome matrices `y_1` and `y_2` it carries the true parameters (`x.true`, `V.true`, `U.true`, `a_1.true`, `a_2.true`, `Z.true`, `b.true`, `theta.true`, `sigma2_y_2.true`, `sigma2_b.true`) so a fit can be scored against the truth.

## Main function

```r
out <- MINDS_algorithm(
  y_1 = data_mixed$y_1,
  y_2 = data_mixed$y_2,
  Nc = 5,
  Nt = 3,
  iter.max = 500,
  plot_trace = FALSE
)

names(out)
out$dic_complete
```

### Key inputs

- `y_1`: binary modality matrix (subjects x binary items). Must be a matrix, not a data frame. `NA` is allowed and handled as missing-at-random: missing cells are masked out of every conditional update and excluded from the likelihood, with no imputation in the fit itself.
- `y_2`: continuous modality matrix (subjects x continuous items). Same handling.
- `Nc`: number of latent clusters (subtypes)
- `Nt`: latent construct dimension
- `var_b`, `ref_cluster`: the identifiability constraints (latent variance per axis, and the cluster pinned at the latent-space origin). Use the same values when comparing a fit against known truth.

Initial cluster memberships come from a `mice` imputation followed by a KAMILA warm start on the mixed continuous + categorical data. Supplying `init_params$Z` skips that warm start; combined with the `update_*` flags it gives a frozen-parameter decode.

Default prior hyperparameters are set to weakly/non-informative values for most components, with informative shrinkage on selected variance terms.

## Identifiability

The MINDS likelihood depends on the parameters only through the linear predictors `(Z_i^T x + b_i) L + a` with `L = [V | U]`, so `x`, `V` and `U` are identified only up to a per-axis scale, an orthogonal rotation, a per-axis sign, and a location shift. `canonicalize_minds()` removes all four in closed form — no optimisation, no use of the true values — and is applied to every retained posterior draw before summarising. Point estimates and credible intervals are therefore well defined and independent of the initialisation.

Cluster label switching is a separate matter and is resolved only when true labels are supplied (`true.values`), for error reporting in simulations.

## Returned output

`MINDS_algorithm()` returns a named list containing:
- `membership`, `membership weight`
- `cluster center`
- `loading to binary modality`, `loading to continuous modality`
- `binary modality intercept`, `continuous modality intercept`
- `b`, `sigma2_y_2`, `sigma2_b`
- `axis order`, `axis sign`
- `credible intervals` (95% intervals for the centres, loadings and intercepts)
- `likelihood trace plot`
- `dic_complete`, `dic4_complete`, `dic5_complete`, `dic6_complete`, `pD4_complete`, `pD5_complete`, `pD6_complete`, `Dbar_complete`

## Model selection

Selecting `Nc` and `Nt` uses the **complete-data DIC** family of Celeux, Forbes, Robert & Titterington (2006), *Deviance Information Criteria for Missing Data Models*, Bayesian Analysis 1(4):651-706. These treat the latent subtype labels `Z` as missing data and build the deviance from the complete-data joint likelihood `f(y, Z | theta)`.

- `dic4_complete` is the **primary** criterion — the paper's recommendation for missing-data models. It plugs in the posterior-mean estimate and integrates the labels over the retained posterior draws. `dic_complete` is an alias for it.
- `dic5_complete` plugs in the joint MAP draw. Diagnostic only: the paper shows it over-penalises for mixtures, since `pD5` scales with the number of subjects.
- `dic6_complete` fixes `theta` at the posterior mean and integrates the labels under their conditional distribution given that estimate (responsibilities, computed exactly).

Smaller is better in all cases. The observed-data (label-marginalised) likelihood is retained as the utility `llk_obs.fun()` but is no longer computed per fit.
