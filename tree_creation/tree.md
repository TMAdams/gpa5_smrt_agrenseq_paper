# Make a tree

## Example alignment construction, all Innovator NLRs plus a reference set of NLRs and NRCs

### Run interproscan over the NLR fasta to predict functional domains

```bash
cd ~/scratch
input_fasta=Innovator_redo/NLR_Annotator/Innovator_NLR_Annotator.fa
outdir=Innovator_redo/interproscan
mkdir -p $outdir
sbatch git_repos/JHI_Code/Pan_NLRome/run_interproscan_1_file.sh $input_fasta $outdir
```

### Extract NB domains from the XML output of interproscan

Be aware that interproscan may change their XML schema. This assumes interproscan version 5.x-x.x is being used with XML schema version 4.5.
This was the most recent at the time, but if your version number is HIGHER than 5.54-87.0 you'll need to check if the namespace heading is different and if the schema version is changed.
This only looks for hits found by hmmr3, if you need ones called by a different analysis you will need to change the script.

```bash
xml=Innovator_redo/interproscan/Innovator_NLR_Annotator.fa.xml
feature=IPR002182
pep_full=Innovator_redo/NLR_Annotator/NLR_Annotator_full.pep.fasta
bed_out=Innovator_redo/NLR_Annotator/NLR_Annotator_NBs.bed
domain_pep=Innovator_redo/NLR_Annotator/NLR_Annotator_NBs.pep.fasta
sbatch git_repos/JHI_Code/Gpa5_RenSeq_paper_prep/Extract_Domain_AA_from_Nuc.sh $xml $feature $pep_full $bed_out $domain_pep
```

### Perform the alignment

```bash
input=Innovator_redo/NLR_Annotator/NLR_Annotator_NBs.pep.fasta
profile=reference_NLRs/combined_plus_outgroups_aligned.pep.fa
output=Innovator_redo/Innovator_aligned_reference_outgroups_clustal_omega.fa
iterations=10 # Feel free to change this, I haven't titrated it
logs=Innovator_redo/Clustal_Omega.log
sbatch git_repos/JHI_Code/H1_Analyses/Clustal_Omega_Profile_Align.sh $profile $input $output $iterations $logs
```
## Begin tree building process

### Check alignment

```bash
screen -
srsh
conda activate r_phylogeny
mkdir -p results/intermediate

R
```

```R
# Load libraries

library("ape")
library("phangorn")

# Initialise graphical parameters

opar <- par(no.readonly = TRUE)
newpar <- opar
newpar$mar <- c(3, 3, 2, 1)
newpar$mgp <- c(2, .7, 0)
newpar$tcl <- NA
newpar$las <- 1
newpar$pch <- 19
par(newpar)

# Load alignment

msa_raw <- read.phyDat(
  "/mnt/shared/scratch/tadams/Innovator_redo/Innovator_aligned_reference_outgroups_clustal_omega.fa",
  format = "fasta",
  type = "AA"
)

# Get initial stats

glance(msa_raw)
unique(msa_raw)

# Check missing data

freqs <- data.frame(
  Frequency = baseFreq(msa_raw, all = TRUE, drop.unused.levels = TRUE)
)
round(freqs, digits = 3)
# I have 55.7% missing data

# Check for duplicate sequence IDs

labels_raw <- names(msa_raw)
anyDuplicated(labels_raw)
# All are unique!

# Visualise the alignment

image(msa_raw, xlab = "Site", ylab = "Sequence", xaxt = "n", show.labels = FALSE,
      xlim = c(0,450))
axis(side = 1, at = seq(from = 0, to = 450, by = 50))
# Yikes that's a lot of gaps

# Quantify missing sites

msa_char <- as.AAbin(msa_raw) |> as.character() # for easy crunching
missing_cols <- apply(msa_char, 2, \(x) sum(x == "-") / nrow(msa_char))
range(missing_cols)
# Between 17.06 and 99.96% missing data

# Calculate quantiles

round(quantile(missing_cols, probs = seq(0, 1, .1)), digits = 3)
sum(missing_cols < 0.85) / 415 # Number to divide by is total number of sites, seen from glance()
# 61% of sites have less than 85% missing data

# Plot a histogram of sites with missing data

hist(
  missing_cols,
  breaks = seq(0, 1, 0.01),
  main = "Proportion missing across sites"
)
abline(v = mean(missing_cols), lwd = 3, lty = 1, col = "black")
abline(v = median(missing_cols), lwd = 3, lty = 2, col = "tomato")
# 85% looks okay

# Remove sites with 85% or more missing data

msa_select <- subset(msa_raw, select = missing_cols < 0.85, site.pattern = FALSE)

# Check for genes missing a lot of data

msa_char <- as.AAbin(msa_select) |> as.character() # re-assign here
missing_rows <- apply(msa_char, 1, \(x) sum(x == "-") / ncol(msa_char))
range(missing_rows)
# between 3 and 87% missing

# Calculate quantiles

round(quantile(missing_rows, probs = seq(0, 1, .1)), digits = 3)
sum(missing_rows < 0.85) / 2837 # Number to divide by is total number of samples, seen from glance()
# 99.7% of genes have less than 85% missing sites.

# Plot a histogram of genes with missing data

hist(
  missing_rows,
  breaks = seq(0, 1, 0.01),
  main = "Proportion missing across sequences"
)
abline(v = mean(missing_rows), lwd = 3, lty = 1, col = "black")
abline(v = median(missing_rows), lwd = 3, lty = 2, col = "tomato")
# 85% looks fine

# Remove sequences with high amounts of missing data

msa_select <- subset(
  msa_select, subset = missing_rows < 0.85, site.pattern = FALSE
)
labels_missing <- setdiff(names(msa_raw), names(msa_select))
labels_removed <- data.frame(label = labels_missing, reason = "missing data")

# Remove identical sequences

msa_unique <- unique(msa_select)
labels_unique <- names(msa_unique)
labels_dup <- setdiff(labels_missing, labels_unique)

labels_removed <- rbind(
  labels_removed,
  data.frame(label = labels_dup, reason = "duplicated")
)

# Write names of removed sequences to file

write.csv(
labels_removed,
"results/labels_removed.csv",
row.names = FALSE
)

# Final checks, see how improved it is
# Basic stats

glance(msa_unique)
unique(msa_unique)

# Check missing data

freqs <- data.frame(
  Frequency = baseFreq(msa_unique, all = TRUE, drop.unused.levels = TRUE)
)
round(freqs, digits = 3)
# Drop to 29%, not bad

# Visualise alignment

image(msa_unique, xlab = "Site", ylab = "Sequence", xaxt = "n", show.labels = FALSE,
      xlim = c(0,300))
axis(side = 1, at = seq(from = 0, to = 300, by = 50))
# That looks better

# Write out new alignment

write.phyDat(msa_unique, file = "results/msa_filtered.fasta",
             format = "fasta")
quit()
```

```bash
exit
exit
```

### Check the various models

```bash
path_to_R_script=/mnt/shared/scratch/tadams/git_repos/JHI_Code/Gpa5_RenSeq_paper_prep/tree_creation/model_check.R
input_file=results/msa_filtered.fasta
sbatch /mnt/shared/scratch/tadams/git_repos/JHI_Code/Gpa5_RenSeq_paper_prep/tree_creation/model_check.sh $path_to_R_script $input_file
```

### Load model testing results and pick which to use

```bash
screen -
srsh
conda activate r_phylogeny
R
```

```R
# Load libraries

library("ape")
library("phangorn")

# Load in previous data and view the model scores

mt_clean_pars <- readRDS("results/intermediate/mt_clean_pars.rds")
mt_clean_pars[order(mt_clean_pars$BIC, decreasing = FALSE), ]

# Select the model with the lowest BIC value, in this case JTT+G

quit()
```

```bash
exit
exit
```

### Create an initial tree

```bash
path_to_R_script=/mnt/shared/scratch/tadams/git_repos/JHI_Code/Gpa5_RenSeq_paper_prep/tree_creation/build_tree.R
input=results/msa_filtered.fasta
inv=FALSE
gamma=TRUE
output=results/tree_clean_noBS_pars.tre
sbatch /mnt/shared/scratch/tadams/git_repos/JHI_Code/Gpa5_RenSeq_paper_prep/tree_creation/build_tree.sh $path_to_R_script $input $inv $gamma $output
```

### OPTIONAL - make a re-rooted DRAFT tree as a quick first-look

```bash
screen -
srsh
conda activate r_phylogeny
R
```

```R
# Load required libraries

library("ape")
library("phangorn")

# Load DRAFT tree

tree_clean_pars <- read.tree("results/tree_clean_noBS_pars.tre")

# Identify outgroup sequences

outgroup1 <- grep(pattern = "Human", x = tree_clean_pars$tip.label, value = TRUE)
outgroup2 <- grep(pattern = "CElegans", x = tree_clean_pars$tip.label, value = TRUE)
outgroup <- c(outgroup1, outgroup2)

# Re-root tree

tree_clean_pars_root <- root(tree_clean_pars, outgroup = outgroup, resolve.root = TRUE)

# Save re-rooted tree

write.tree(tree_clean_pars_root, file = "results/tree_rooted_clean_noBS_pars.tre")

quit()
```

```bash
exit
exit
```

### Bootstrap the tree - for the exemplar tree this took 15 hours

```bash
path_to_R_script=/mnt/shared/scratch/tadams/git_repos/JHI_Code/Gpa5_RenSeq_paper_prep/tree_creation/bootstrap_tree.R
sbatch /mnt/shared/scratch/tadams/git_repos/JHI_Code/Gpa5_RenSeq_paper_prep/tree_creation/bootstrap_tree.sh $path_to_R_script
```

### Load bootstraps onto tree and re-root it

```bash
screen -
srsh
conda activate r_phylogeny
R
```

```R
# Load required libraries

library("ape")
library("phangorn")

# Load previous results

fit_clean_opt <- readRDS("results/intermediate/fit_clean_opt_pars.rds")
bs_clean_pars <- readRDS("results/bs_clean_pars1000.rds")

# Adds bootstraps to tree

tree_clean_pars_bs <- plotBS(fit_clean_opt$tree, bs_clean_pars)

# Identify outgroup sequences

outgroup1 <- grep(pattern = "Human", x = tree_clean_pars_bs$tip.label, value = TRUE)
outgroup2 <- grep(pattern = "CElegans", x = tree_clean_pars_bs$tip.label, value = TRUE)
outgroup <- c(outgroup1, outgroup2)

# Re-root tree

tree_clean_pars_bs_root <- root(tree_clean_pars_bs, outgroup = outgroup, resolve.root = TRUE)

# Save re-rooted tree

write.tree(tree_clean_pars_bs_root, file = "results/tree_rooted_clean_BS1000_pars.tre")

# Remove bootstrap values to allow partitioning

tree_clean_pars_bs_root$node.label <- NULL

# Save tree without bootstrap values

write.tree(tree_clean_pars_bs_root, file = "results/tree_rooted_clean_BS1000_noBS_pars.tre")

quit()
```

```bash
exit
exit
```

### Partition the tree

```bash
path_to_jar=/mnt/shared/scratch/tadams/apps/PhyloPart_v2.1/PhyloPart_v2.1.jar
tree=results/tree_rooted_clean_BS1000_noBS_pars.tre
threshold=0.05
output=results/tree_part_0_05.txt
sbatch /mnt/shared/scratch/tadams/git_repos/H1_Analyses/run_PhyloPart.sh $tree $threshold $output $path_to_jar
```
