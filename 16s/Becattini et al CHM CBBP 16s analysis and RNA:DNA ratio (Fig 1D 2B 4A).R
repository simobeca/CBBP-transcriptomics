# Becattini et al. 'Rapid transcriptional and metabolic adaptation of intestinal microbes to host immune activation'
# Codes used to obtain Figures 1d, 2b and 4a = 16s rRNA gene analysis
# Use the 'DADA2_SB093_SB172_table.txt' obtained through DADA2 analyses of the sequencing data

library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gtools)
library(grid)
library(gridExtra)
library(stringr)

setwd('/Users/becattis/Desktop/CHM revision/Scripts to be uploaded/CBBP-transcriptomics/16s/')

df <- read.delim('DADA2_SB093_SB172_table.txt')


#### SB093 SPF mice flag/aCD3 ####

SB093_flag <- df %>% 
  filter(grepl('SB093', sample)) %>%
  mutate(sample=gsub('\\..pool781', '', sample)) %>%
  mutate(sample=gsub('SB093.', '', sample)) %>%
  select(sample, Species, pctseqs) %>%
  filter(grepl('0h|flag', sample)) %>%
  dplyr::mutate(Group = str_extract(sample, pattern="0h|6h|24h")) %>%
  dplyr::mutate(Species=ifelse(Species %in% c(
    'Blautia marasmi',
    'Bacteroides sartorii',
    '[Clostridium] bolteae',
    'Parabacteroides distasonis'), as.character(Species), "Other")) %>%
  mutate(Species=gsub('\\[Clostridium\\]', 'Clostridium', Species)) %>%
  mutate(Species=gsub('marasmi', 'producta', Species)) %>%
  group_by(sample,Species, Group) %>%
  dplyr::summarize(pctseqs=sum(pctseqs)) %>%
  ungroup() 

Figure_1d <- SB093_flag %>%
  ggplot(aes(fill=factor(Species, levels = c(    'Clostridium bolteae',
                                                 'Blautia producta',
                                                 'Bacteroides sartorii',
                                                 'Parabacteroides distasonis',
                                                 'Other')), 
             y=pctseqs, 
             x=factor(sample))) + 
  geom_bar(position="stack", stat="identity", color='black', size=0.1) +
  theme_bw(base_size = 16) + 
  scale_fill_manual(values=c(
    'Bacteroides sartorii'='palegreen4', 
    'Blautia producta'='pink',
    'Clostridium bolteae'='firebrick2',
    'Parabacteroides distasonis'='seashell3',
    'Other'='black'))  +
  theme(axis.text.x = element_blank(),
        strip.background = element_blank()) +
  labs(fill='Species',
       x ='Sample',
       title = 'SB093 Flagellin cecum post-treatment') +
  facet_wrap(~factor(Group, levels=c('0h', '6h', '24h')),nrow = 1,scales = "free_x")


Figure_1d



SB093_aCD3 <- df %>% 
  filter(grepl('SB093', sample)) %>%
  mutate(sample=gsub('\\..pool781', '', sample)) %>%
  mutate(sample=gsub('SB093.', '', sample)) %>%
  select(sample, Species, pctseqs) %>%
  filter(grepl('0h|aCD3', sample)) %>%
  dplyr::mutate(Group = str_extract(sample, pattern="0h|6h|24h")) %>%
  dplyr::mutate(Species=ifelse(Species %in% c(
    'Blautia marasmi',
    'Bacteroides sartorii',
    '[Clostridium] bolteae',
    'Parabacteroides distasonis'), as.character(Species), "Other")) %>%
  mutate(Species=gsub('\\[Clostridium\\]', 'Clostridium', Species)) %>%
  mutate(Species=gsub('marasmi', 'producta', Species)) %>%
  group_by(sample,Species, Group) %>%
  dplyr::summarize(pctseqs=sum(pctseqs)) %>%
  ungroup() 

Figure_2b <- SB093_aCD3 %>%
  ggplot(aes(fill=factor(Species, levels = c(    'Clostridium bolteae',
                                                 'Blautia producta',
                                                 'Bacteroides sartorii',
                                                 'Parabacteroides distasonis',
                                                 'Other')), 
             y=pctseqs, 
             x=factor(sample))) + 
  geom_bar(position="stack", stat="identity", color='black', size=0.1) +
  theme_bw(base_size = 16) + 
  scale_fill_manual(values=c(
    'Bacteroides sartorii'='palegreen4', 
    'Blautia producta'='pink',
    'Clostridium bolteae'='firebrick2',
    'Parabacteroides distasonis'='seashell3',
    'Other'='black'))  +
  theme(axis.text.x = element_blank(),
        strip.background = element_blank()) +
  labs(fill='Species',
       x ='Sample',
       title = 'SB093 aCD3 cecum post-treatment') +
  facet_wrap(~factor(Group, levels=c('0h', '6h', '24h')),nrow = 1,scales = "free_x")

Figure_2b

#### SB172 GF-mice ####

SB172 <- df %>% 
  filter(grepl('3mix.pre|3mix.post|4mix.pre|4mix.post|SB172', sample)) %>%
  select(sample, Species, pctseqs) %>%
  dplyr::mutate(Group = case_when(
    grepl('3mix|3MIX', sample) ~ '3-mix',
    grepl('4mix|4MIX', sample) ~ '4-mix'
  )) %>%
  dplyr::mutate(Species=ifelse(Species %in% c(
    'Blautia marasmi',
    'Bacteroides sartorii',
    '[Clostridium] bolteae',
    'Parabacteroides distasonis'), as.character(Species), "Other")) %>%
  mutate(Species=gsub('\\[Clostridium\\]', 'Clostridium', Species)) %>%
  mutate(Species=gsub('marasmi', 'producta', Species)) %>%
  group_by(sample,Species, Group) %>%
  dplyr::summarize(pctseqs=sum(pctseqs)) %>%
  ungroup() 

SB172_pre_gg <- SB172 %>%
  filter(grepl('3mix.pre|4mix.pre|SB172', sample)) %>%
  mutate(sample=gsub('\\..pool957', '', sample)) %>%
  mutate(sample=gsub('3mix.pre.', '', sample)) %>%
  mutate(sample=gsub('4mix.pre.', '', sample)) %>%
  mutate(sample=gsub('SB172.INPUT3MIX', '3mix_input', sample)) %>%
  mutate(sample=gsub('SB172.INPUT4MIX', '4mix_input', sample)) %>%
  ggplot(aes(fill=factor(Species, levels = c(    'Clostridium bolteae',
                                                 'Blautia producta',
                                                 'Bacteroides sartorii',
                                                 'Parabacteroides distasonis',
                                                 'Other')), 
             y=pctseqs, 
             x=factor(sample, levels=c('3mix_input','1','2','3', '4', '5', '6',
                                       '4mix_input', '10', '11', '12', '7', '8', '9')))) + 
  geom_bar(position="stack", stat="identity", color='black', size=0.1) +
  theme_bw(base_size = 16) + 
  scale_fill_manual(values=c(
    'Bacteroides sartorii'='palegreen4', 
    'Blautia producta'='pink',
    'Clostridium bolteae'='firebrick2',
    'Parabacteroides distasonis'='seashell3',
    'Other'='black'))  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_blank()) +
  labs(fill='Species',
       x ='Sample',
       title = 'SB172 GF pellets pre-treatment') +
  facet_wrap(~Group,nrow = 1,scales = "free_x")


SB172_pre_gg


Figure_4a <- SB172 %>%
  filter(grepl('3mix.post|4mix.post|SB172', sample)) %>%
  mutate(sample=gsub('\\..pool957', '', sample)) %>%
  mutate(sample=gsub('3mix.post.', '', sample)) %>%
  mutate(sample=gsub('4mix.post.', '', sample)) %>%
  mutate(sample=gsub('SB172.INPUT3MIX', '3mix_input', sample)) %>%
  mutate(sample=gsub('SB172.INPUT4MIX', '4mix_input', sample)) %>%
  ggplot(aes(fill=factor(Species, levels = c(    'Clostridium bolteae',
                                                 'Blautia producta',
                                                 'Bacteroides sartorii',
                                                 'Parabacteroides distasonis',
                                                 'Other')), 
             y=pctseqs, 
             x=factor(sample, levels=c('3mix_input','1','2','3', '4', '5', '6',
                                       '4mix_input', '10', '11', '12', '7', '8', '9')))) + 
  geom_bar(position="stack", stat="identity", color='black', size=0.1) +
  theme_bw(base_size = 16) + 
  scale_fill_manual(values=c(
    'Bacteroides sartorii'='palegreen4', 
    'Blautia producta'='pink',
    'Clostridium bolteae'='firebrick2',
    'Parabacteroides distasonis'='seashell3',
    'Other'='black'))  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.background = element_blank()) +
  labs(fill='Species',
       x ='Sample',
       title = 'SB172 GF pellets post-treatment') +
  facet_wrap(~Group,nrow = 1,scales = "free_x")

Figure_4a


########### RATIO RNA/DNA

### read numbers were extracted for each bacterium from FeatureCounts as explained in the main script (SB093)
### each corresponding df was build as indicated below for BP, then the 4 dfs were attached to make an All_reads file, here provided as a .txt for convenience

#BP_reads <- data.frame('Reads' = colSums(BP_FtCounts$stat[,-c(1,2)])) %>%
#  rownames_to_column(var='Sample') %>%
#  mutate(Organism = 'BP')

All_reads <- read.delim2('All_reads.txt', dec='.')

All_reads_pct <- All_reads %>% 
  dplyr::group_by(Sample) %>%
  dplyr::mutate(pct_reads=Reads/sum(Reads)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(sample=Sample) %>%
  dplyr::mutate(sample=paste0('m',sample, sep=''))


SB093_flag_pct <- SB093_flag %>%
  dplyr::filter(!Species=='Other') %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(pctseqs=pctseqs/sum(pctseqs)) %>%
  dplyr::ungroup()

SB093_aCD3_pct <- SB093_aCD3 %>%
  dplyr::filter(!Species=='Other') %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(pctseqs=pctseqs/sum(pctseqs)) %>%
  dplyr::ungroup()

SB093_16s_all_pct <- bind_rows(SB093_flag_pct, SB093_aCD3_pct) %>%
  dplyr::distinct() %>% # remove 0h which is now there twice
  dplyr::mutate(Species = case_when(
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
  dplyr:: mutate(pctseqs = case_when(
    Organism == 'BS' ~ pctseqs/7,
    Organism == 'PD' ~ pctseqs/7,
    Organism == 'CB'  ~ pctseqs/5,
    Organism == 'BP' ~ pctseqs/5)) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(pctseqs=pctseqs/sum(pctseqs)) %>%
  dplyr::ungroup()

SB093_ratio_rna_dna <- SB093_16s_all_pct_norm %>% 
  dplyr::left_join(All_reads_pct, by=c('sample', 'Organism')) %>%
  dplyr::mutate(RNA_DNA_ratio =pct_reads/pctseqs )


Supplementary_Figure_1D <- SB093_ratio_rna_dna %>%
  dplyr::mutate(Time_point = str_extract(pattern='0h|6h|24h', sample)) %>%
  dplyr::filter(grepl('0h|flag', sample)) %>%
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


Supplementary_Figure_2E <- SB093_ratio_rna_dna %>%
  dplyr::mutate(Time_point = str_extract(pattern='0h|6h|24h', sample)) %>%
  dplyr::filter(grepl('0h|aCD3', sample)) %>%
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
       title = 'SB093 RNA aCD3 reads') +
  facet_wrap(~factor(Organism, levels=c('CB', 'BP', 'BS', 'PD')), ncol=1 , scales = "fixed", strip.position = 'right') +
  expand_limits(y=0)





