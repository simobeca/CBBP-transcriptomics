library(RColorBrewer)
library(plyr)
library(stringr)
library(ggplot2)
library(grDevices)
library(tidyr)
library(yingtools2)
library(reshape2)
library(ggrepel)
library(gplots)
library(gtools)
library(dplyr)
library(ggsci)

setwd('/Users/becattis/Desktop/CHM revision/Scripts to be uploaded/CBBP-transcriptomics/Operons/')

# read in the Operon list obtained through Rockhopper
operons_RH <- read.delim2('operons_RH_full.txt') %>% 
  select(-Start, -Stop, -Strand) %>%
  mutate(Operon_composition = case_when(
    !is.na(Operon_composition) ~ Operon_composition,
    is.na(Operon_composition) ~ ID
    )) %>%
  rename(Gene=ID)
 

operon_id <- operons_RH %>%
  #arrange(Gene) %>%
  group_by(Organism) %>%
  distinct(Operon_composition) %>%
  mutate(operon_id = paste0(Organism, '_', row_number())) %>%
  ungroup()

operons_RH <- left_join(operons_RH, operon_id, by = c('Organism', "Operon_composition"))


##### main figure graphs ####

number_hypoth_in_operons <- operons_RH %>%
  filter(Number.of.Genes >= 2) %>%
  group_by(Organism, Operon_composition) %>% 
  summarize(Perc_hypot = sum(grepl('hypothetical protein', Name))/sum(!grepl('hypothetical protein', Name))) %>%
  mutate(Perc_hypot = case_when(
    Perc_hypot == 0 ~ 'All annotated',
    Perc_hypot == 'Inf' ~ 'All unknown',
    Perc_hypot > 0 ~ 'Some unknown'
  )) 
number_hypoth_in_operons <- number_hypoth_in_operons %>%
  group_by(Organism, Perc_hypot) %>%
  summarise(Counts = n()) %>%
  ungroup()

Figure_5D <- number_hypoth_in_operons %>%
  ggplot(aes(x=factor(Organism, levels = c('CB', 'BP', 'BS', 'PD')), 
             y=Counts, 
             fill=factor(Perc_hypot, levels = c('All unknown', 'Some unknown', 'All annotated')))) +
  geom_bar(stat="identity", width=0.75, position = 'fill') +
  scale_fill_manual(values=c('lavenderblush4', 'lavenderblush3','lavenderblush2' )) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        aspect.ratio = 1,
        legend.position = 'right',
        legend.title = element_blank()) 

perc_genes_in_operons <- operons_RH %>% 
  group_by(Organism) %>% 
  summarize(sum(Number.of.Genes>1)/sum(Number.of.Genes>0)) %>%
  ungroup()
colnames(perc_genes_in_operons) <- c('Organism', 'Percentage_operons') 

n_operons <- operons_RH %>%
  dplyr::filter(Number.of.Genes >= 2) %>%
  dplyr::group_by(Organism) %>%
  dplyr::summarize(N_Operons = length(unique(operon_id))) %>% 
  ungroup()


n_operons %>% ggplot(aes(x=factor(Organism, levels=c('CB', 'BP', 'BS', 'PD')), y=N_Operons, fill=Organism)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=c(    'BS'='palegreen4', 
                                 'BP'='pink',
                                 'CB'='firebrick2',
                                 'PD'='darkgray')) +
  theme_bw(base_size = 16) +
  theme(strip.background = element_blank(),
        aspect.ratio = 2,
        legend.position = "none") 

Figure_5A <- perc_genes_in_operons %>% 
  ggplot(aes(x=factor(Organism, levels=c('CB', 'BP', 'BS', 'PD')), y=Percentage_operons, fill=Organism)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=c(    'BS'='palegreen4', 
                                 'BP'='pink',
                                 'CB'='firebrick2',
                                 'PD'='darkgray')) +
  theme_bw(base_size = 16) +
  theme(strip.background = element_blank(),
        aspect.ratio = 2,
        legend.position = "none") +
  geom_text(aes(label=n_operons$N_Operons), position=position_dodge(width=0.9), vjust=-0.25, size=7) +
  ylim(0,1)

length_operons_boxplot <- operons_RH %>%
  filter(Number.of.Genes >= 2) %>%
  group_by(Organism, Operon_composition) %>%
  summarize(Number_of_Genes=mean(Number.of.Genes)) %>%
  ungroup() 

Figure_5B <- length_operons_boxplot %>% 
  ggplot(aes(x=factor(Organism, levels=c('CB', 'BP', 'BS', 'PD')), y=Number_of_Genes, fill=Organism)) +
  geom_boxplot(notch=T, outlier.shape = NA) +
  stat_summary(fun=mean, colour="black", geom="point", 
               shape=16, size=5, show_legend = T) +
  ylim(1, 8) + 
  scale_fill_manual(values=c(    'BS'='palegreen4', 
                                 'BP'='pink',
                                 'CB'='firebrick2',
                                 'PD'='darkgray')) +
  theme_bw(base_size = 16) +
  theme(strip.background = element_blank(),
        aspect.ratio = 2,
        legend.position = "none") 


distr <- length_operons_boxplot %>%
  group_by(Organism, Number_of_Genes) %>%
  summarise(Counts = n()) %>%
  ungroup() 

Figure_5C <- distr %>%
  ggplot(aes(x=Number_of_Genes, y= Counts, fill=Organism)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values=c(    'BS'='palegreen4', 
                                 'BP'='pink',
                                 'CB'='firebrick2',
                                 'PD'='darkgray')) +
  theme_bw(base_size = 14) +
  theme(strip.background = element_blank(),
        aspect.ratio = 0.3,
        legend.position = "none") +
  scale_y_log10() +
  # ylim(1, 1000) +
  geom_text(aes(label=Counts), position=position_dodge(width=0.9), vjust=-0.25, size=7) +
  facet_wrap(.~ factor(Organism, levels=c('CB', 'BP', 'BS', 'PD')), ncol=1, strip.position = 'right') 


####### Supplementary Figure 5  #####

# read in results df from aCD3-treated mice in SB093 and SB172 (SPF and GF, respectively). Here provided as a table 'all_data
# select genes significant in both experiments 
# retain only those with abs(log2FC)>1 and concordant sign

alldata <- read.delim2('all_data.txt', dec = '.') %>%
  left_join(operons_RH, by=c('Organism','Gene')) %>%
  dplyr::group_by(Gene,Organism) %>%
  dplyr::filter(sum(padj<0.1)==2) %>% 
  dplyr::filter((sum(log2FoldChange>1) ==2)|(sum(log2FoldChange<(-1)) ==2)) %>% 
  dplyr::ungroup()

# calculate mean and sd for gene regulation
alldata2 <- alldata %>%
  dplyr::mutate(Gene_Name=gsub('utilization system for glycans and polysaccharides', '', Gene_Name)) %>% #just to make the name shorter for graphical purposes
  dplyr::mutate(operon2=coalesce(operon_id,Gene)) %>%
  dplyr::group_by(Gene,Gene_Name,operon2,Organism) %>%
  dplyr::summarize(mean.log2foldchange=mean(log2FoldChange),
                   sd.log2foldchange=sd(log2FoldChange)/sqrt(n())) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(sd.log2foldchange)) %>% #removes duplicates due to Gene_name
  dplyr::group_by(operon2) %>%
  dplyr::mutate(max.mean=max(mean.log2foldchange)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(operon2) %>% ### the next three lines get rid of operons with genes going in opposite direction
  dplyr::filter((sum(mean.log2foldchange>0 ) == n()) | (sum(mean.log2foldchange<0 ) == n())) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(max.mean,operon2,mean.log2foldchange) %>%
  dplyr::mutate(Gene2=factor(Gene,levels=unique(Gene)),
                operon2=factor(operon2,levels=unique(operon2))) 

# generate key to assign colors
key <- alldata2 %>%
  dplyr::group_by(Organism) %>%
  dplyr::select(Organism, operon2) %>%
  unique() %>%
  dplyr::mutate(num=row_number(operon2)) %>%
  ungroup()

# assign colors
alldata2 <- alldata2 %>%
  dplyr::left_join(key, by='operon2') %>%
  dplyr::mutate(color_gg=ifelse(num %% 2,'black','red')) %>%
  dplyr::arrange(max.mean,operon2,mean.log2foldchange) 

# plot
pdf("Supplementary_Figure_5.pdf", height=130, width=27)
ggplot(alldata2) +
  geom_col(aes(x=Gene2,y=mean.log2foldchange, fill=color_gg)) +
  scale_fill_manual(values=c(red='firebrick2', black='gray60')) +
  geom_errorbar((aes(x=Gene2,ymin=mean.log2foldchange-sd.log2foldchange,ymax=mean.log2foldchange+sd.log2foldchange)), width=0, size=0.3, color='gray20') +
  coord_flip() +
  geom_text(aes(x=Gene2, y=0, label=Gene_Name, hjust=ifelse(mean.log2foldchange>0,1,0), colour=color_gg), size=7) +
  scale_colour_manual(values=c(red='firebrick2', black='#2D2926FF')) +
  theme_bw(base_size = 24, base_line_size = 0) +
  theme(strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(legend.position="none") +
  facet_grid(Organism.x ~. ,space="free",scales="free") 
dev.off()

#write.xlsx(alldata2, file='top operon list.xlsx')
