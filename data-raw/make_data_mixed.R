source("MINDS_Rpackage/data_generation.R")
dir.create("MINDS/data", recursive = TRUE, showWarnings = FALSE)
save(data_mixed, file = "MINDS/data/data_mixed.rda")
