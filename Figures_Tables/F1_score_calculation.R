#!/usr/bin/env Rscript

#########################calculate F1 score##################################
matrix_value_calculate <- function(input_table, phe_table, start_col, group_number){
  negative_start <- number_positive + 1
  TP <- sum(input_table[1:number_positive, start_col:group_number] > 97)
  FN <- sum(input_table[1:number_positive, start_col:group_number] < 97)
  TN <- sum(input_table[negative_start:nrow(input_table), start_col:group_number] < 97)
  FP <- sum(input_table[negative_start:nrow(input_table), start_col:group_number] > 97)
  PPV <- TP/(TP + FP)
  TPR <- TP/(TP + FN)
  F1 <- 2*TP/(2*TP + FP + FN)
  value_matrix <- c(F1, PPV, TPR)
  return(value_matrix)
}

#####################F1 for table##########################
creat_list <- function(looplist, header, cal_type, confusion_matrix){
  F1_list <- c()
  PPV_list <- c()
  TPR_list <- c()
  candidate_number <- c()
  calculate_number <- 0
  for (l in looplist) {
    calculate_number <- calculate_number + l
    if(cal_type == "group"){s_col <- 1} else {s_col <- calculate_number}
    candidate_number <- c(candidate_number, calculate_number)
    total_list <- matrix_value_calculate(confusion_matrix, phenotypice_table, s_col, calculate_number)
    F1_list <- c(F1_list, total_list[1])
    PPV_list <- c(PPV_list, total_list[2])
    TPR_list <- c(TPR_list, total_list[3])
  }
  total_table <- data.frame(header, F1_list, candidate_number, PPV_list, TPR_list)
  return(total_table)
}

###################order table#################
table_order <- function(phenotypice_table, drenseq_score){
  phenotypice_table$phenotype <- as.numeric(phenotypice_table$phenotype)
  dren_phe <- merge(phenotypice_table, drenseq_score)
  order_table <- arrange(dren_phe, desc(phenotype))
  row.names(order_table) <- order_table$gene
  order_table <- order_table[,-1]
  t_order_table <- t(order_table)
  t_order_table <- dplyr::as_tibble(t_order_table, rownames = "gene")
  all_table <- merge(group_list, t_order_table)
  order_all_table <- arrange(all_table,desc(association_score))
  final_table <- t(order_all_table)
  colnames(final_table) <- final_table[1,]
  final_table <- data.frame(final_table)
  final_table <- final_table %>% select(phenotype, everything())
  final_table[2,1] <- ""
  return(list(group_con=final_table, individual=order_all_table))
}

###############confusion matrix#########################
confusion_matrix_creat <- function(group_list, final_table){
  groups_name <- levels(as.factor(group_list[,2]))
  ordered_group_list <- rev(sort(as.numeric(groups_name)))
  confusion_matrix <- as.data.frame(lapply(data.frame(final_table[4:nrow(final_table)-1, 2:ncol(final_table)], stringsAsFactors = F), as.numeric))
  return(list(loop_list=ordered_group_list, con_matrix=confusion_matrix))
}

#################################group F1####################################
group_F1_creat <- function(group_list, ordered_group_list, confusion_matrix){
  calculate_list <- c()
  for (n in ordered_group_list) {
    calculate_list <- c(calculate_list, sum(group_list[,2][-1] == n))
}
  group_table <- creat_list(calculate_list, ordered_group_list, "group", confusion_matrix)
  names(group_table)[1] <- "Filter group"
  return(group_table)
}

#################################group F1 figure################################
group_F1_figure <- function(group_table, ordered_group_list){
  candidates <- group_table$candidate_number
  f <- ggplot(group_table, aes(x=ordered_group_list)) +
    geom_line(aes(y=F1_list, color="F1 scores"), size=1) +
    geom_point(aes(y=F1_list), size=1, color="black") +
    geom_line(aes(y=rescale(candidates, c(0.5, 1)), color="Candidates"), size=1) +
    geom_point(aes(y=rescale(candidates, c(0.5, 1))), size=1, color="black") +
    scale_y_continuous(name="F1 scores", breaks = pretty_breaks(10), sec.axis = sec_axis(~rescale(.,c(min(candidates),max(candidates))), name="Number of candidates"))+
    scale_x_continuous(name="Association scores", breaks = seq(min(ordered_group_list), max(ordered_group_list)), 1) +
    theme_bw() +
    theme_classic()+
    theme(
      axis.title.y = element_text(size=13),
      axis.title.y.right = element_text(size=13),
      legend.title=element_blank(),
      legend.position = "right", 
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    ) +
    scale_colour_manual(name="Confusion value", values = c("F1 scores" = rgb(0.2, 0.6, 0.9, 1), "Candidates" = "#69b3a2")) +
    ggtitle(figure_title) 
return(f)
}

#############################individual F1######################################
individual_F1_creat <- function(confusion_matrix, order_all_table){
  individual_list <- rep(1, times=ncol(confusion_matrix))
  candidate_list <- colnames(confusion_matrix)
  individual_table <- creat_list(individual_list, candidate_list, "individual", confusion_matrix)
  individual_table$candidate_number <- order_all_table$association_score[2:length(order_all_table$association_score)]
  names(individual_table)[1] <- "Candidates"
  names(individual_table)[3] <- "Association_score"
  individual_table <- individual_table %>% select(Candidates, Association_score,everything())
  individual_table <- arrange(individual_table, desc(F1_list))
  F1_table <- data.frame(individual_table[,1:3])
  names(F1_table)[1] <- "gene" 
  F1_table <- F1_table[,-2]
  F1_table[nrow(F1_table) + 1, ] = c("phenotype", "")
  all_table_r <- order_all_table
  all_table_m <- merge(all_table_r, F1_table)
  all_table_m <- arrange(all_table_m, desc(F1_list), desc(association_score))
  all_table_m <- all_table_m %>% select(gene, association_score,F1_list, everything())
  t_all_table_m <- data.frame(t(all_table_m))
  colnames(t_all_table_m) <- t_all_table_m[1,]
  t_all_table_m <- t_all_table_m %>% select(phenotype, everything())
  t_all_table_m[2,1] <- ""
  #t_all_table_m[4,1] <- ""
  return(list(f1=individual_table,all=t_all_table_m))
}


##################################Start###################################
library("dplyr")
library("ggplot2") 
library("tibble")
library("scales")
library("magrittr")
library("lubridate")
library("reshape2")
args <- commandArgs(TRUE)
pwd <- args[1]
input_drenseq <- args[2]
input_phenotype <- args[3]
input_group_list <- args[4]
output_name <- args[5]

drenseq_score <- read.table(input_drenseq, sep = "\t", stringsAsFactors = F, header =T)
phenotypice_table <- read.table(input_phenotype, sep = "\t", header = T)
group_list <- read.table(input_group_list, sep = "\t", header = T)
names(group_list) <- c("gene","association_score")
max_association_score <- max(group_list$association_score)
phenotype_merge_row <- c("phenotype", max_association_score)
group_list <- rbind(phenotype_merge_row,group_list)
names(phenotypice_table) <- c("gene","phenotype")

output_table_name <- paste(output_name,"_drenseq_F1_table.csv",sep="")
confusion_table_name <- paste(output_name,"_confuse_table.csv",sep="")
individual_confusion_table_name <- paste(output_name,"_individual_confusion_table.csv",sep="")
png_name <- paste(output_name,"_F1.png",sep="")
figure_title <- paste("F1 scores and number of candidates of ", output_name,sep="")
number_positive <- sum(phenotypice_table$phenotype >= 0)


first_table <- table_order(phenotypice_table, drenseq_score)
second_table <- confusion_matrix_creat(group_list, first_table$group_con)
figure_table <- group_F1_creat(group_list, second_table$loop_list, second_table$con_matrix)
f <- group_F1_figure(figure_table,second_table$loop_list)
ggsave(filename = png_name, f, units="cm", width=30, height=12, dpi = 600)
#png(filename = png_name,width = 800,height = 300, type="cairo")
#f
#dev.off()
F1_in_table <- individual_F1_creat(second_table$con_matrix, first_table$individual)
write.table(F1_in_table$all, output_table_name, sep = ",", col.names = F, row.names = T)
write.table(F1_in_table$f1, individual_confusion_table_name, sep = ",", col.names = T, row.names = F)
write.table(figure_table, confusion_table_name, sep = ",", col.names = T, row.names = F)

