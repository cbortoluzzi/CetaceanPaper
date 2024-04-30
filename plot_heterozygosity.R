#!/usr/bin/env Rscript

library(argparse)
library(ggplot2)


args = commandArgs(trailingOnly=TRUE)


df <- read.table(args[1], header = F, col.names=c('Chrom', 'Start', 'End', 'Nsites', 'NSNP', 'Heterozygosity'))

# Set color based on chromosome number
df$color <- NA
df$color[df$Chrom == 1 | df$Chrom == 3 | df$Chrom == 5 | df$Chrom == 7 | df$Chrom == 9 | df$Chrom == 11 | df$Chrom == 13 | df$Chrom == 15 | df$Chrom == 17 | df$Chrom == 19 | df$Chrom == 21 | df$Chrom == 23] <- "#0868ac"
df$color[df$Chrom == 2 | df$Chrom == 4 | df$Chrom == 6 | df$Chrom == 8 | df$Chrom == 10 | df$Chrom == 12 | df$Chrom == 14 | df$Chrom == 16 | df$Chrom == 18 | df$Chrom == 20 | df$Chrom == 22] <- "#7a0177"


# Plot genome-wide heterozygosity
het <- subset(df, df$Nsites >= 600000)
avg <- mean(het$Heterozygosity)
std <- sd(het$Heterozygosity)
avg
std
p <- ggplot(het, aes(x=Start, y=Heterozygosity/1000000))+geom_bar(stat = 'identity', fill = het$color)+facet_grid(~Chrom, scales='free_x', space='free_x', switch = 'x')+theme_classic()+scale_x_continuous(expand=c(0,0))+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+ylab("Heterozygosity per bp")+xlab("")+theme(axis.text.y = element_text(color="black"), strip.text = element_text(size = 8), strip.background = element_blank(), panel.spacing = unit(0, "lines"))
figure <- sub('\\.txt$', '.pdf', args[1])
ggsave(figure, width = 25, height = 8, units = "cm")


# Plot distribution of heterozygosity
d <- ggplot(het, aes(x=Heterozygosity/1000000))+geom_histogram(alpha=.8, fill = '#43a2ca', color = '#0868ac')+theme_classic()+ylab("Count")+xlab("Heterozygosity per bp")+theme(axis.text.y = element_text(color="black"), axis.text.x =  element_text(color="black"))
figure <- sub('\\.txt$', '.density.pdf', args[1])
ggsave(figure, width = 10, height = 10, units = "cm")

