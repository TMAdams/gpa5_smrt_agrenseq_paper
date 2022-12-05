#!/usr/bin/env Rscript

# Load libraries

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
    help = "Path to input alignment file"),
    make_option("--inv", type = "logical",
    help = "Logical value for inv parameter for tree model, TRUE or FALSE"),
    make_option("--gamma", type = "logical",
    help = "Logical value for gamma parameter for tree model, TRUE or FALSE"),
    make_option("--outp", type = "character",
    help = "Location to write output tree")
)

opt <- parse_args(OptionParser(option_list = opt_list))
inp <- opt$inp
inv <- opt$inv
gamma <- opt$gamma
outp <- opt$outp

# trimmed and sequences removed

msa_clean <- read.phyDat(
  inp,
  format = "fasta",
  type = "AA"
)

mt_clean_pars <- readRDS("results/intermediate/mt_clean_pars.rds")

# Create a pml object to build a maximum likelihood tree

fit_clean_pars <- as.pml(mt_clean_pars, "BIC")

# Optimise the model for the tree

fit_clean_opt_pars <- optim.pml(
  object = fit_clean_pars,
  optInv = inv,
  optGamma = gamma,
  rearrangement = "NNI"
)

# Write out files

saveRDS(fit_clean_opt_pars, "results/intermediate/fit_clean_opt_pars.rds")

write.tree(fit_clean_opt_pars$tree, file = outp)
