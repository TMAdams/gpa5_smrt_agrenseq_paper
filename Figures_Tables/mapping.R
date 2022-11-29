#!/usr/bin/env Rscript

#3 files are needed: 1. .gff file for backgroud bar; 2. Chromose size file with header (format: chromosome	length(bp)); 3. .bed file for candidate contigs

#load library
library(ggplot2)
require(grid)
library(ggrepel)



input_drenseq <- args[2]
input_phenotype <- args[3]
input_group_list <- args[4]
output_name <- args[5]
args <- commandArgs(TRUE)
#set working directory 
pwd <- args[1]
setwd(pwd)
#loading the backgroud bar with gff file
reference_genome_gff <-args[2]
#loading chromose size file with format: chromosome	length(bp) 
chromose_size <- args[3]
#loading .bed file for candidate contigs
bed_file <- args[4]
#setting of output file and title of figure
gene_name <- args[5]
output_file<-"Gpa5_F1_mapping.png"
fig_title<- "Gpa5 candidates locations"

#store data
genSeqFile <- read.delim(reference_genome_gff, header=FALSE, stringsAsFactors=FALSE)
DM_chromoSizes <- read.delim(chromose_size)
candidates_R_genes <- read.delim(bed_file,header=FALSE, stringsAsFactors=FALSE)
names(candidates_R_genes) <- c("chr","start","end")

#x-axis setting
roundUp <- function(x) 10^ceiling(log10(x))


#set limits of x-axis according to chromosome size
dummy <- data.frame(DM_chromoSizes$chromosome,DM_chromoSizes$length.bp./1000000)
names(dummy) <- c("chr","start")
dummy <- dummy[ which(dummy$chr!="ch00"), ]
dummy$plot <- 0

#import NLR-Annotator estimated gene locations
genSeqGenes_plot <- data.frame(genSeqFile$V1)
names(genSeqGenes_plot) <- c("chr")
#change the data to number
genSeqGenes_plot$start <- as.integer(genSeqFile$V4)

#binSize given in kb
binSize <- 1000

#########################making data table#########################
#chromosome list
all_ch <- c("chr01","chr02","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12")
all_ch_bins <- data.frame()
count = 1
ch_num = 1
binLoop = binSize*1000
for (ch in all_ch) {
  print(ch)
  loop <- 0
  ch_size <- DM_chromoSizes[ch_num,2] 
  while (loop <= ch_size) {
    all_ch_bins[count,1] <- ch
    all_ch_bins[count,2] <- (loop + (binLoop/2))/1000000
    all_ch_bins[count,3] <- length(genSeqGenes_plot[ which(genSeqGenes_plot$chr==ch&genSeqGenes_plot$start>loop&genSeqGenes_plot$start<loop+binLoop), ]$chr)
    all_ch_bins[count,4] <- length(candidates_R_genes[ which(candidates_R_genes$chr==ch&candidates_R_genes$start>loop&candidates_R_genes$start<loop+binLoop), ]$chr)
    loop <- loop + binLoop
    count = count + 1
  }
  ch_num = ch_num + 1
} 
names(all_ch_bins) <- c("chr","start","backgroundbars","candidates")
all_ch_bins[all_ch_bins == 0] <- NA
all_ch_bins$backgroundbars[is.na(all_ch_bins$backgroundbars)] <- 0
all_ch_bins$chr <- factor(all_ch_bins$chr, levels = c("chr01","chr02","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12"), labels = c("chr01","chr02","chr03","chr04","chr05","chr06","chr07","chr08","chr09","chr10","chr11","chr12"))


#########################plot graphs#########################
  plot_graphs <- TRUE
  if (plot_graphs) {
    ggplot(all_ch_bins,aes(x=start,y=backgroundbars,group=chr)) + 
      facet_grid(~chr,scales="free_x",space="free_x") + 
      geom_ribbon(aes(ymin=0,ymax=backgroundbars,fill=chr), alpha=0.35) +
      geom_point(aes(x=start,y=candidates,group=chr,colour=chr)) + 
      #geom_text(aes(x=start, y=candidates, label=paste(candidates, 'candidates'), vjust = -1, hjust = 'left'), size = 3) +
      #ggrepel::geom_text_repel(aes(x=start, y = candidates, label=paste(candidates, 'candidates')), all_ch_bins, size = 3, box.padding = 0.5, point.padding = 0.8, min.segment.length = 1, segment.color = "black", show.legend = F) +
      geom_text_repel(aes(x=start, y = candidates, label=paste(candidates, 'candidates')),
                      size = 3,
                      direction="y") +
      theme(panel.grid.minor = element_line(colour="white"), 
            panel.grid.major = element_line(colour="white"), 
            panel.background = element_rect(fill="white"), 
            panel.spacing = unit(c(0.2),"cm"), 
            legend.position="none", strip.background = element_blank(), 
            strip.text.x = element_text(colour = "red", size = 10, hjust = 0.5), 
            axis.line.x = element_line(colour = "black"), 
            axis.line.y = element_line(colour = "black"), 
            axis.text.x = element_text(colour="grey20",size=8),
            plot.title = element_text(hjust = 0.5)) + 
      scale_x_continuous(name="Physical location (Mb)",breaks=seq(0,roundUp(max(all_ch_bins$start)),roundUp(max(all_ch_bins$start))/10)) + 
      scale_y_continuous("Number of candidates per 1 Mb",expand=c(0.01,0)) +
      #geom_label(aes(label=candidates), all_ch_bins, alpha=0, nudge_y=3) +
      ggtitle(fig_title)
    ggsave(filename = paste(output_file,sep=""), units="cm", width=36, height=12, dpi = 600)
  }
