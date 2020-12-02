
### this scripts uses as an input a similarity matrix of protein sequences corresponding to the top 100 transcripts at baseline in SB093
### the matrix was generated via multiple Muscle alignments 
### this code generates the plot displayed in Supplementary Figure 1C



### load libraries ####
library(RColorBrewer)
library(DESeq2)
#library(plyr)
library(dplyr)
library(tibble)
library(pheatmap)
library(Rsubread)
library(clusterProfiler)
library(stringr)
library(ggplot2)
library(grDevices)
library(genefilter)
library(xlsx)
library(tidyr)
library(openxlsx)
library(githubinstall)
library(yingtools2)
library(phyloseq)
library(reshape2)
library(ggrepel)
library(gplots)
library(qgraph)
library(ggfortify)
library(cluster)
library(cowplot)
library(ggpubr)

setwd("/MyFolder/")

df <- read.delim('4mix_top100_Matrix.txt')

colnames(df) <- df$X
df <- df[-1,]
colnames(df)[1] <- 'Organism'
df <-  melt(df, id = c('Organism')) %>% 
  mutate(Group_1 = str_extract(pattern='BP|BS|CB|PD', Organism)) %>%
           mutate(Group_2 = str_extract(pattern='BP|BS|CB|PD', variable)) 

df <- df %>% group_by(Organism, Group_2) %>%
  filter(value==max(value)) %>%
  ungroup()

pdf('Baseline_identity.pdf', useDingbats = F)
df %>%
  ggplot(aes(x=factor(Group_1, levels=c('CB', 'BP', 'BS', 'PD')), y=value, fill=Group_1)) +
  geom_boxplot(notch=TRUE, outlier.size = 0.1) +
  #geom_jitter(size=0.1) +
  scale_fill_manual(values=c(    'BS'='palegreen4', 
                                  'BP'='pink',
                                  'CB'='firebrick2',
                                  'PD'='darkgray')) +

  theme_bw(base_size = 16) +
  labs(y='Gene product identity', x='4mix member', title='Identity of top 100 baseline-expressed proteins', colour='4mix member') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_blank(),
        aspect.ratio = 2.5) +
  facet_grid(.~factor(Group_2, levels=c('CB', 'BP', 'BS', 'PD')))
dev.off()



