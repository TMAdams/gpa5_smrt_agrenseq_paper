#!/usr/bin/env Rscript

library(ggplot2)
library(ggrepel)
library(cli)
library(dplyr)

args <- commandArgs(TRUE)
agrenseq_raw_table <- args[1]
agrenseq_table <- read.table(agrenseq_raw_table, sep="\t")
figure_title <- args[2]
png_name <- args[3]

names(agrenseq_table) <- c("contig_name", "contig", "association_score", "number_of_kmers")
#getting the candidate list
h_score <- max(agrenseq_table$association_score)
filter_score <- ceiling(h_score/2)
candidate_list <- agrenseq_table$contig[agrenseq_table$association_score >= filter_score]
candidate_unique_list <- unique(candidate_list)
candidate_name_list <- agrenseq_table$contig_name[agrenseq_table$association_score >= filter_score]
candidate_unique_name_list <- unique(candidate_name_list)
agrenseq_table <- agrenseq_table %>% mutate(candidate = ifelse(agrenseq_table$contig %in% candidate_list, 1, -1))

x_max <- agrenseq_table$contig[length(agrenseq_table$contig)]
r_gene <- as.numeric(levels(factor(agrenseq_table$contig[agrenseq_table$contig_name == 'tig00000301'])))
y_r_gene <- max(as.numeric(levels(factor(agrenseq_table$association_score[agrenseq_table$contig_name == 'tig00000301']))))



f <- ggplot(agrenseq_table, aes(x = contig, y = association_score, size = number_of_kmers, label = contig_name)) + 
  geom_point(aes(color = candidate >= 0)) +
  guides(colour = "none") +
  labs(size = "Number of k-mers") + 
  scale_x_continuous(name="Assembled NLR-containing contigs based on PacBio HiFi reads from cultivar Innovator")+
  scale_y_continuous(name="Association for k-mers per contig") + 
  ggtitle(figure_title) +
  scale_size_area(max_size = 18) +
  scale_color_manual(values=c("black","red")) + 
  coord_fixed(ratio = 9) +
  theme_bw() + 
  theme_classic()+
  theme(plot.title = element_text(size = 40, face = "bold", hjust = 0.5)) +
  theme(text = element_text(size=30)) +
  geom_hline(aes(yintercept=filter_score), colour="tomato", linetype="dashed") +
  geom_hline(aes(yintercept=h_score), colour="red", linetype="dashed") 
  #annotate(geom = "segment", x = r_gene + 10, y = y_r_gene + 12, xend = r_gene, yend = y_r_gene + 2,arrow = arrow(length = unit(5, "mm"))) +
  #annotate(geom = "text", x = r_gene + 10.5, y = y_r_gene + 14, label = "tig00000301/Gpa5", size = 10, hjust = 'left')

png(filename = png_name,width = 2000,height = 1000, type="cairo")
f
dev.off()
