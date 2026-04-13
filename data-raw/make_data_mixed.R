pkgload::load_all("MINDS")
data_mixed <- generate_data_mixed(seed.no = 103)
dir.create("MINDS/data", recursive = TRUE, showWarnings = FALSE)
save(data_mixed, file = "MINDS/data/data_mixed.rda")
