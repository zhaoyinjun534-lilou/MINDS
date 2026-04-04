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

## Included data

The package includes a simulated dataset:

```r
library(MINDS)
data(data_mixed)
str(data_mixed)
```

`data_mixed` is a list with:
- `y_1`: binary outcome matrix
- `y_2`: continuous outcome matrix

This dataset is intended as a reproducible toy example for the modeling workflow in the manuscript-style simulation setting.

## Main function

```r
out <- MINDS_algorithm(
  y_1 = data_mixed$y_1,
  y_2 = data_mixed$y_2,
  Nc = 4,
  Nt = 3,
  iter.max = 50,
  plot_trace = FALSE
)

names(out)
out$ic
```

### Key inputs

- `y_1`: binary modality matrix (subjects x binary items)
- `y_2`: continuous modality matrix (subjects x continuous items)
- `Nc`: number of latent clusters (subtypes)
- `Nt`: latent construct dimension

Default prior hyperparameters are set to weakly/non-informative values for most components, with informative shrinkage on selected variance terms.

## Returned output

`MINDS_algorithm()` returns a named list containing:
- `membership`
- `cluster center`
- `loading to binary modality`
- `loading to continuous modality`
- `binary modality intercept`
- `continuous modality intercept`
- `memberhip weight`
- `likelihood trace plot`
- `ic`

The `ic` field is computed from the DIC helper (`dic.fun`) and can be used for model comparison across candidate specifications.
