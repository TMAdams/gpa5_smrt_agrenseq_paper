# Make a tree

All code was run on the UK Crop Diversity HPC <https://www.cropdiversity.ac.uk/>, we recommend running the commands on a similar high-performance system to reduce runtime due to computationally demanding analyses.

## Example alignment construction, all Innovator NLRs plus a reference set of NLRs and NRCs

### Run interproscan over the NLR fasta to predict functional domains

```bash
cd ~/scratch
input_fasta=/path/to/input/fasta
outdir=/path/to/output/directory
mkdir -p $outdir
interproscan.sh -goterms -iprlookup -i $input_fasta -d $outdir -dp -pa -t n -T /path/to/tmp/dir
```

### Extract NB domains from the XML output of interproscan

Be aware that interproscan may change their XML schema. This assumes interproscan version 5.x-x.x is being used with XML schema version 4.5.
This was the most recent at the time, but if your version number is HIGHER than 5.54-87.0 you'll need to check if the namespace heading is different and if the schema version is changed.
This only looks for hits found by hmmr3, if you need ones called by a different analysis you will need to change the script.

```bash
xml=/path/to/interproscan/xml
feature=IPR002182
pep_full=/path/to/write/full/protein/sequences
bed_out=/path/to/write/BED/of/NB/domains
domain_pep=/path/to/write/protein/sequences/of/NB/domains

python /path/to/parse_ipr_xml.py --xml $xml --feature $feature --pep_out $pep_full --bed_out $bed_out

bedtools getfasta -fo $domain_pep -fi $pep_full -bed $bed_out
```

### Perform the alignment

#### Align high-confidence reference sequences

```bash
input=/path/to/high/confidence/sequences
output=/path/to/write/initial/alignment
iterations=10
logs=/path/to/write/log
clustalo -i $input -o $output --iterations $iterations --threads /number/of/threads --log $logs
```

#### Add lower confidence sequences to the alignment

```bash
input=/path/to/NB/fasta
profile=/path/to/alignment/of/references
output=/path/to/write/aligned/output
iterations=10
logs=/path/to/write/log
clustalo --profile1 $profile -i $input -o $output --iterations $iterations --threads /number/of/threads --log $logs
```
## Begin tree building process

### Check alignment

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
  "/path/to/alignment/fasta",
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
# We have 55.7% missing data initially

# Check for duplicate sequence IDs

labels_raw <- names(msa_raw)
anyDuplicated(labels_raw)
# All are unique!

# Visualise the alignment

image(msa_raw, xlab = "Site", ylab = "Sequence", xaxt = "n", show.labels = FALSE,
      xlim = c(0,450))
axis(side = 1, at = seq(from = 0, to = 450, by = 50))
# There are many gaps in the alignment

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
# 85% is a reasonable threshold

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
# 85% is a reasonable threshold

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
"/path/to/write/csv/of/removed/samples",
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
# Drops to 29%

# Visualise alignment

image(msa_unique, xlab = "Site", ylab = "Sequence", xaxt = "n", show.labels = FALSE,
      xlim = c(0,300))
axis(side = 1, at = seq(from = 0, to = 300, by = 50))
# Alignment now appears much more informative

# Write out new alignment

write.phyDat(msa_unique, file = "/path/to/write/filtered/alignment",
             format = "fasta")
```

### Check the various models

Executes bootstrap_tree.R

```bash
input_file=/path/to/cleaned/alignment
/path/to/model_check.R --inp $input_file
```

### Load model testing results and pick which to use

```R
# Load libraries

library("ape")
library("phangorn")

# Load in previous data and view the model scores

mt_clean_pars <- readRDS("/path/to/cleaned/data/R/object")
mt_clean_pars[order(mt_clean_pars$BIC, decreasing = FALSE), ]

# Select the model with the lowest BIC value, in this case JTT+G
```

### Create an initial tree

Executes build_tree.R

```bash
input=/path/to/cleaned/fasta
inv=FALSE
gamma=TRUE
output=/path/to/write/initial/tree
/path/to/build_tree.R --inp $input --inv $inv --gamma $gamma --outp $output
```

### OPTIONAL - make a re-rooted DRAFT tree as a quick first-look

```R
# Load required libraries

library("ape")
library("phangorn")

# Load DRAFT tree

tree_clean_pars <- read.tree("/path/to/initial/tree")

# Identify outgroup sequences

outgroup1 <- grep(pattern = "Human", x = tree_clean_pars$tip.label, value = TRUE)
outgroup2 <- grep(pattern = "CElegans", x = tree_clean_pars$tip.label, value = TRUE)
outgroup <- c(outgroup1, outgroup2)

# Re-root tree

tree_clean_pars_root <- root(tree_clean_pars, outgroup = outgroup, resolve.root = TRUE)

# Save re-rooted tree

write.tree(tree_clean_pars_root, file = "/path/to/rooted/tree")
```

### Bootstrap the tree - for this tree it took 15 hours

Executes bootstrap_tree.R

```bash
/path/to/bootstrap_tree.R --cores /number/of/threads/to/use
```

### Load bootstraps onto tree and re-root it

```R
# Load required libraries

library("ape")
library("phangorn")

# Load previous results

fit_clean_opt <- readRDS("/path/to/tree/R/data/object")
bs_clean_pars <- readRDS("/path/to/bootstrapped/data/object")

# Adds bootstraps to tree

tree_clean_pars_bs <- plotBS(fit_clean_opt$tree, bs_clean_pars)

# Identify outgroup sequences

outgroup1 <- grep(pattern = "Human", x = tree_clean_pars_bs$tip.label, value = TRUE)
outgroup2 <- grep(pattern = "CElegans", x = tree_clean_pars_bs$tip.label, value = TRUE)
outgroup <- c(outgroup1, outgroup2)

# Re-root tree

tree_clean_pars_bs_root <- root(tree_clean_pars_bs, outgroup = outgroup, resolve.root = TRUE)

# Save re-rooted tree

write.tree(tree_clean_pars_bs_root, file = "/path/to/write/bootstrapped/tree")

# Remove bootstrap values to allow partitioning

tree_clean_pars_nobs_vals_root <- tree_clean_pars_bs_root
tree_clean_pars_nobs_vals_root$node.label <- NULL

# Save tree without bootstrap values

write.tree(tree_clean_pars_nobs_vals_root, file = "/path/to/write/tree/without/bootstrap/values/in/newick")
```

### Partition the tree

```bash
path_to_jar=/path/to/PhyloPart_v2.1.jar
tree=/path/to/bootstrapped/tree/without/bootstraps/in/newick
threshold=0.05
output=/path/to/write/partitioning/results
java -jar $path_to_jar $tree $threshold -o"$output"
```
