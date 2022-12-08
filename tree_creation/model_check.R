#!/usr/bin/env Rscript

# Import libraries

library("ape") # requirement for phangorn
library("phangorn") # load alignment and build phylogeny
library("optparse")

# Set up graphical parameters

opar <- par(no.readonly = TRUE)
newpar <- opar
newpar$mar <- c(3, 3, 2, 1)
newpar$mgp <- c(2, .7, 0)
newpar$tcl <- NA
newpar$las <- 1
newpar$pch <- 19
par(newpar)

# Parse CLI arguments

opt_list <- list(
    make_option("--inp", type = "character",
    help = "Path to input alignment file")
)

opt <- parse_args(OptionParser(option_list = opt_list))
inp <- opt$inp

# Read in alignment

msa_clean <- read.phyDat(
  inp,
  format = "fasta",
  type = "AA"
)

# Test different substitution models for the tree

mt_clean_pars <- modelTest(
  object = msa_clean,
  model = c("JTT", "LG", "WAG"),
  multicore = FALSE
)

# Save this data for later use

saveRDS(mt_clean_pars, file = "/path/to/write/mt_clean_pars.rds")
