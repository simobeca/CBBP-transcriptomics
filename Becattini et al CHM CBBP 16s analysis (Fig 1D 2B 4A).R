# Becattini et al CHM
# Codes used to obtain Figures 1d, 2b and 4a = 16s rRNA gene analysis
# Use the '16s _for_transcriptome_paper' obtained through DADA2 analyses of the sequencing data

library(RColorBrewer)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gtools)
library(grid)
library(gridExtra)
library(stringr)


setwd('/Users/becattis/Desktop/16s _for_transcriptome_paper/')
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



save.image('16s_for_transcriptome.RData')