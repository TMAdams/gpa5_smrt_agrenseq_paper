#!/usr/bin/env Rscript

# Load libraries

library("parallel")
library("ape") # requirement for phangorn
library("phangorn") # load alignment and build phylogeny
library("optparse")

# Parse CLI arguments

opt_list <- list(
    make_option("--cores", type = "double",
    help = "Number of CPU cores to use")
)

opt <- parse_args(OptionParser(option_list = opt_list))
cores <- opt$cores

# Load in the object with the un-bootstrapped tree

fit_clean_opt <- readRDS("/path/to/write/fit_clean_opt_pars.rds")

# Specify random seed

set.seed(123)

bs_clean <- bootstrap.pml(fit_clean_opt, bs = 1000, optNni = TRUE,
                          multicore = TRUE, mc.cores = cores)

# Save output

saveRDS(bs_clean, "/path/to/write/bs_clean_pars1000.rds")
