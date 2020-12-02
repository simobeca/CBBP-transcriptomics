library(RColorBrewer)
library(DESeq2)
#library(plyr)

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
#githubinstall("yingtools2")
library(yingtools2)
#source("https://bioconductor.org/biocLite.R")
#biocLite("phyloseq")
library(phyloseq)
library(reshape2)
library(ggrepel)
library(gplots)
library(qgraph)
library(ggfortify)
library(cluster)
library(cowplot)
library(ggpubr)
library(dplyr)
library(gtools)

setwd("/Users/becattis/Desktop/16s _for_transcriptome_paper/")
#load("/Users/becattis/Desktop/SB093_rerun/SB093_rerun.RData")
load("16s_for_transcriptome.RData")

##### FROM SB093 original code Transcriptional activity levels #####


BP_reads <- data.frame('Reads' = colSums(BP_FtCounts$stat[,-c(1,2)])) %>%
  rownames_to_column(var='Sample') %>%
  mutate(Organism = 'BP')
BS_reads <- data.frame('Reads' = colSums(BS_FtCounts$stat[,-c(1,2)])) %>%
  rownames_to_column(var='Sample') %>%
  mutate(Organism = 'BS')
CB_reads <- data.frame('Reads' = colSums(CB_FtCounts$stat[,-c(1,2)])) %>%
  rownames_to_column(var='Sample') %>%
  mutate(Organism = 'CB') 
PD_reads <- data.frame('Reads' = colSums(PD_FtCounts$stat[,-c(1,2)])) %>%
  rownames_to_column(var='Sample') %>%
  mutate(Organism = 'PD')   

All_reads <- dplyr::bind_rows(BP_reads, BS_reads, CB_reads, PD_reads) %>%
  dplyr::mutate(Sample = gsub("(_IGO.*)", "", Sample)) %>%
  dplyr::mutate(Sample = gsub('Sample_m', '', Sample)) %>%
  dplyr::mutate(Sample = gsub('Sample_sb93_m', '', Sample)) %>%
  dplyr::filter(!grepl(pattern='REPL|RNAlat', All_reads$Sample)) 

All_reads_pct <- All_reads %>% 
  group_by(Sample) %>%
  mutate(pct_reads=Reads/sum(Reads)) %>%
  ungroup() %>%
  dplyr::rename(sample=Sample) %>%
  mutate(sample=paste0('m',sample, sep=''))
  

SB093_flag_pct <- SB093_flag %>%
  filter(!Species=='Other') %>%
  group_by(sample) %>%
  mutate(pctseqs=pctseqs/sum(pctseqs)) %>%
  ungroup()

SB093_aCD3_pct <- SB093_aCD3 %>%
  filter(!Species=='Other') %>%
  group_by(sample) %>%
  mutate(pctseqs=pctseqs/sum(pctseqs)) %>%
  ungroup()

SB093_16s_all_pct <- bind_rows(SB093_flag_pct, SB093_aCD3_pct) %>%
  distinct() %>% # remove 0h which is now there twice
  mutate(Species = case_when(
    Species == 'Bacteroides sartorii' ~ 'BS',
    Species == 'Parabacteroides distasonis' ~ 'PD',
    Species == 'Clostridium bolteae' ~ 'CB',
    Species == 'Blautia producta' ~ 'BP'
  )) %>%
  dplyr::rename(Organism=Species)

### according to PATRIC, PD and BS have 7 copies of 16s, BP and CB have 5, 
### take those into account by dividing BS and PD times 7 and BP CB times 5
### then sum/normalize again
SB093_16s_all_pct_norm <- SB093_16s_all_pct %>% 
  mutate(pctseqs = case_when(
    Organism == 'BS' ~ pctseqs/7,
    Organism == 'PD' ~ pctseqs/7,
    Organism == 'CB'  ~ pctseqs/5,
    Organism == 'BP' ~ pctseqs/5)) %>%
  group_by(sample) %>%
  mutate(pctseqs=pctseqs/sum(pctseqs)) %>%
            ungroup()

SB093_ratio_rna_dna <- SB093_16s_all_pct_norm %>% 
  left_join(All_reads_pct, by=c('sample', 'Organism')) %>%
  mutate(RNA_DNA_ratio =pct_reads/pctseqs )

pdf('SB093_RNA_DNA_ratio_flag.pdf')
SB093_ratio_rna_dna %>%
  mutate(Time_point = str_extract(pattern='0h|6h|24h', sample)) %>%
  filter(grepl('0h|flag', sample)) %>%
  ggplot(aes(fill=factor(Organism), 
             y=RNA_DNA_ratio, 
             x=factor(Time_point, levels=c('0h', '6h', '24h'))))  + 
  geom_point(shape=21, color='black', size=4) +
  theme_bw(base_size = 16) + 
  scale_fill_manual(values=c(
    'BS'='palegreen4', 
    'BP'='pink',
    'CB'='firebrick2',
    'PD'='seashell3'))  +
  theme(strip.background = element_blank(),
        aspect.ratio = 0.75) +
  labs(fill='Species',
       x ='Sample',
       title = 'SB093 RNA flag reads') +
  facet_wrap(~factor(Organism, levels=c('CB', 'BP', 'BS', 'PD')), ncol=1 , scales = "fixed", strip.position = 'right') +
  expand_limits(y=0)
dev.off()

pdf('SB093_RNA_DNA_ratio_aCD3.pdf')
SB093_ratio_rna_dna %>%
  mutate(Time_point = str_extract(pattern='0h|6h|24h', sample)) %>%
  filter(grepl('0h|aCD3', sample)) %>%
  ggplot(aes(fill=factor(Organism), 
             y=RNA_DNA_ratio, 
             x=factor(Time_point, levels=c('0h', '6h', '24h'))))  + 
  geom_point(shape=21, color='black', size=4) +
  theme_bw(base_size = 16) + 
  scale_fill_manual(values=c(
    'BS'='palegreen4', 
    'BP'='pink',
    'CB'='firebrick2',
    'PD'='seashell3'))  +
  theme(strip.background = element_blank(),
        aspect.ratio = 0.7) +
  labs(fill='Species',
       x ='Sample',
       title = 'SB093 RNA flag reads') +
  facet_wrap(~factor(Organism, levels=c('CB', 'BP', 'BS', 'PD')), ncol=1 , scales = "fixed", strip.position = 'right') +
  expand_limits(y=0)
dev.off()















pdf('SB093_RNA_reads_overtime.pdf')
All_reads %>%
  ggplot(aes(fill=factor(Organism, levels = c(   'CB',
                                                 'BP',
                                                 'BS',
                                                 'PD')), 
             y=Reads, 
             x=factor(Sample, levels=mixedsort(Sample)))) + 
  geom_bar(position="fill", stat="identity", color='black', size=0.1) +
  theme_bw(base_size = 16) + 
  scale_fill_manual(values=c(
    'BS'='palegreen4', 
    'BP'='pink',
    'CB'='firebrick2',
    'PD'='seashell3'))  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_blank()) +
  labs(fill='Species',
       x ='Sample',
       title = 'SB093 RNA reads') 
dev.off()

pdf('SB093_RNA_flagellin_reads_overtime.pdf')
All_reads %>%
  mutate(Group=str_extract(pattern='0h|6h_flag|24h_flag|6h_aCD3|24h_aCD3', Sample)) %>%
  filter(Group %in% c('0h', '6h_flag', '24h_flag')) %>%
  ggplot(aes(fill=factor(Organism, levels = c(   'BP',
                                                 'BS',
                                                 'CB',
                                                 'PD')), 
             y=Reads, 
             x=factor(Sample, levels=mixedsort(Sample)))) + 
  geom_bar(position="fill", stat="identity", color='black', size=0.1) +
  theme_bw(base_size = 16) + 
  scale_fill_manual(values=c(
    'BS'='palegreen4', 
    'BP'='pink',
    'CB'='firebrick2',
    'PD'='seashell3'))  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_blank()) +
  labs(fill='Species',
       x ='Sample',
       title = 'SB093 RNA flagellin reads') +
  facet_wrap(~factor(Group, levels=c('0h', '6h_flag', '24h_flag')), nrow =1 ,scales = "free_x")
dev.off()

pdf('SB093_RNA_aCD3_reads_overtime.pdf')
All_reads %>%
  mutate(Group=str_extract(pattern='0h|6h_flag|24h_flag|6h_aCD3|24h_aCD3', Sample)) %>%
  filter(Group %in% c('0h', '6h_aCD3', '24h_aCD3')) %>%
  ggplot(aes(fill=factor(Organism, levels = c(   'BP',
                                                 'BS',
                                                 'CB',
                                                 'PD')), 
             y=Reads, 
             x=factor(Sample, levels=mixedsort(Sample)))) + 
  geom_bar(position="fill", stat="identity", color='black', size=0.1) +
  theme_bw(base_size = 16) + 
  scale_fill_manual(values=c(
    'BS'='palegreen4', 
    'BP'='pink',
    'CB'='firebrick2',
    'PD'='seashell3'))  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_blank()) +
  labs(fill='Species',
       x ='Sample',
       title = 'SB093 RNA aCD3 reads') +
  facet_wrap(~factor(Group, levels=c('0h', '6h_aCD3', '24h_aCD3')), nrow =1 ,scales = "free_x")
dev.off()


All_reads_pct_means <- All_reads %>%
  group_by(Sample) %>%
  dplyr::mutate(pctreads = Reads/sum(Reads)) %>%
  ungroup() %>%
  mutate(Time_point = str_extract(pattern='0h|6h|24h', Sample)) %>%
  filter(grepl('0h|aCD3', Sample)) %>%
  group_by(Organism, Time_point) %>%
  dplyr::summarise(meanreads = mean(pctreads)) %>%
  ungroup() 

pdf('SB093_RNAreads_aCD3.pdf')
All_reads %>%
  group_by(Sample) %>%
  dplyr::mutate(pctreads = Reads/sum(Reads)) %>%
  ungroup() %>%
  mutate(Time_point = str_extract(pattern='0h|6h|24h', Sample)) %>%
  filter(grepl('0h|aCD3', Sample)) %>%
  ggplot(aes(fill=factor(Organism, levels = c('BP',
                                              'BS',
                                              'CB',
                                              'PD')), 
             y=pctreads, 
             x=factor(Time_point, levels=mixedsort(Time_point)))) + 
  geom_point(shape=21, color='black', size=4) +
  theme_bw(base_size = 16) + 
  scale_fill_manual(values=c(
    'BS'='palegreen4', 
    'BP'='pink',
    'CB'='firebrick2',
    'PD'='seashell3'))  +
  theme(axis.text.x = element_blank(),
        strip.background = element_blank()) +
  labs(fill='Species',
       x ='Sample',
       title = 'SB093 RNA aCD3 reads') +
  facet_wrap(~Organism, ncol=1 , scales = "free") +
  expand_limits(y=0)
dev.off()

pdf('SB093_RNAreads_flag.pdf')
All_reads %>%
  group_by(Sample) %>%
  dplyr::mutate(pctreads = Reads/sum(Reads)) %>%
  ungroup() %>%
  mutate(Time_point = str_extract(pattern='0h|6h|24h', Sample)) %>%
  filter(grepl('0h|flag', Sample)) %>%
  ggplot(aes(fill=factor(Organism, levels = c('BP',
                                              'BS',
                                              'CB',
                                              'PD')), 
             y=pctreads, 
             x=factor(Time_point, levels=mixedsort(Time_point)))) + 
  geom_point(shape=21, color='black', size=4) +
  theme_bw(base_size = 16) + 
  scale_fill_manual(values=c(
    'BS'='palegreen4', 
    'BP'='pink',
    'CB'='firebrick2',
    'PD'='seashell3'))  +
  theme(axis.text.x = element_blank(),
        strip.background = element_blank()) +
  labs(fill='Species',
       x ='Sample',
       title = 'SB093 RNA flag reads') +
  facet_wrap(~Organism, ncol=1 , scales = "free") +
  expand_limits(y=0)
dev.off()