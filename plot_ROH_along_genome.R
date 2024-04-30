#!/usr/bin/env Rscript

library(argparse)
library(ggplot2)


args = commandArgs(trailingOnly=TRUE)


roh <- read.table(args[1], header = F, col.names=c('Chrom', 'Start', 'End', 'Size', 'Heterozygosity'))
chrom <- read.table(args[2], header = F, sep=",", col.names = c('GenBank', 'Chrom', 'Length'))

autosomes <- subset(chrom, chrom$Chrom != 'X' & chrom$Chrom != 'Y')

# Consider only ROHs longer than 10 Kb
roh_100kb <- subset(roh, roh$Size >= 10)
roh_len <- sum(roh_100kb$Size)
roh_len * 10000
dim(roh_100kb)

# Plot distribution of ROHs along the genome
plot <- ggplot()+geom_segment(aes(x=0, y=as.numeric(autosomes$Chrom), xend=autosomes$Length, yend=as.numeric(autosomes$Chrom)), color = '#f0f0f0', size = 5)+geom_segment(aes(x=roh_100kb$Start, y=roh_100kb$Chrom, xend=roh_100kb$End, yend=roh_100kb$Chrom), color = '#3182bd', size = 5)+theme_classic()+theme(axis.text.y = element_text(color="black"), axis.text.x =  element_text(color='black'))+ylab("Chromosome")+xlab("Position along the genome (bp)")+scale_y_continuous(breaks=roh_100kb$Chrom, labels=roh_100kb$Chrom)
figure <- sub('\\.txt$', '.pdf', args[1])
ggsave(figure, width = 25, height = 15, units = "cm")

