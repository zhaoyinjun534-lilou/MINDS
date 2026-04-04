# MINDS

MINDS is an R package for Bayesian integrative clustering of mixed binary and continuous outcomes.

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
