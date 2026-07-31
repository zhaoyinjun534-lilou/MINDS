# Rebuild data/data_mixed.rda. Run from the package root:
#   Rscript data-raw/make_data_mixed.R
# The dataset is generate_data_mixed() at its defaults (Nb = 1000, Nd1 = 10,
# Nd2 = 12, Nc = 5, Nt = 3) and now carries the true parameters alongside the
# outcome matrices, so an example fit can be scored against the truth.
pkgload::load_all(".", quiet = TRUE)
data_mixed <- generate_data_mixed(seed.no = 103)
dir.create("data", recursive = TRUE, showWarnings = FALSE)
save(data_mixed, file = "data/data_mixed.rda", compress = "xz")
