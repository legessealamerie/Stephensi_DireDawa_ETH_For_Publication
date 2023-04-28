library(dcifer)

filtered_dat <- readr::read_rds("processed/filtered_dat.rds")

mydsmp <- dcifer::formatDat(filtered_dat, "sample_id", "locus", "allele")
coi <- dcifer::getCOI(mydsmp)
afreq <- dcifer::calcAfreq(mydsmp, coi)
dres <- dcifer::ibdDat(mydsmp, coi, afreq,
    pval = TRUE, confint = TRUE, rnull = 0,
    alpha = 0.05, nr = 1e3
)

readr::write_rds(dres, "processed/dcifer_results.rds")
