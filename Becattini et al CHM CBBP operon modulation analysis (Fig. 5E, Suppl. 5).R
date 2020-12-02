library(RColorBrewer)
library(DESeq2)
library(plyr)
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
library(gtools)
library(gsubfn)
library(VennDiagram)
library(grid)
library(dplyr)

setwd("/MyFolder/")

# read in the Operon list obtained through Rockhopper
operons_RH <- read.delim2('/MyFolder/operons_RH_full.txt') %>% 
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

# read in results df from SaCD3-treated mice in B093 and SB172 (SPF and GF, respectively)
SB093_BP <- read.table("SB093_BP_aCD36h_vs_0h.txt") %>% mutate(experiment="EXP_1") %>% mutate(Organism="BP") %>% mutate (Gene_number = gsub('ID=fig\\|33035.10.', '', Gene))
SB172_BP <- read.table("SB172_BP_4mix_6h_DEG_allsamples.txt") %>% mutate(experiment="EXP_2") %>% mutate(Organism="BP") %>% mutate(Gene_Name = Gene_name) %>% select(-Gene_name) %>% mutate(Gene = gsub(";Name=.*", "", Gene))
SB093_BS <- read.table("SB093_BS_aCD36h_vs_0h.txt") %>% mutate(experiment="EXP_1") %>% mutate(Organism="BS") %>% mutate (Gene_number = gsub('ID=fig\\|671267.8.', '', Gene))
SB172_BS <- read.table("SB172_BS_4mix_6h_DEG_allsamples.txt") %>% mutate(experiment="EXP_2") %>% mutate(Organism="BS") %>% mutate(Gene_Name = Gene_name) %>% select(-Gene_name) %>% mutate(Gene = gsub(";Name=.*", "", Gene)) %>% mutate (Gene_number = gsub('ID=fig\\|671267.8.', '', Gene))
SB093_CB <- read.table("SB093_CB_aCD36h_vs_0h.txt") %>% mutate(experiment="EXP_1") %>% mutate(Organism="CB") %>% mutate (Gene_number = gsub('ID=fig\\|208479.6.', '', Gene)) 
SB172_CB <- read.table("SB172_CB_4mix_6h_DEG_allsamples.txt") %>% mutate(experiment="EXP_2") %>% mutate(Organism="CB") %>% mutate(Gene_Name = Gene_name) %>% select(-Gene_name) %>% mutate(Gene = gsub(";Name=.*", "", Gene)) %>% mutate (Gene_number = gsub('ID=fig\\|208479.6.', '', Gene)) 
SB093_PD <- read.table("SB093_PD_aCD36h_vs_0h.txt") %>% mutate(experiment="EXP_1") %>% mutate(Organism="PD") %>% mutate (Gene_number = gsub('ID=fig\\|823.92.', '', Gene))
SB172_PD <- read.table("SB172_PD_4mix_6h_DEG_allsamples.txt") %>% mutate(experiment="EXP_2") %>% mutate(Organism="PD") %>% mutate(Gene_Name = Gene_name) %>% select(-Gene_name) %>% mutate(Gene = gsub(";Name=.*", "", Gene)) %>% mutate (Gene_number = gsub('ID=fig\\|823.92.', '', Gene))

# generate cumulative df
  alldata <- bind_rows(SB093_BP, SB172_BP, SB093_BS, SB172_BS, SB093_CB, SB172_CB, SB093_PD, SB172_PD) %>%
 # dplyr::mutate(Gene_name = gsub("(.*);Name=", "", Gene)) %>% mutate(Gene_name = gsub("(\\(E.*)", "", Gene_name))# %>%
 # mutate(Gene = gsub("(.*);Name=", "", Gene))
 # mutate (Gene_number = gsub('ID=fig|208479.6.', '', Gene_Name) %>%
  mutate(Gene_Name = ifelse(grepl("hypothetical protein|Hypothetical protein", Gene_Name), 'hypothetical protein', Gene_Name)) %>%
  mutate(Gene_Name = ifelse(grepl("hypothetical protein", Gene_Name), paste0(Gene_Name, " ", Gene_number), Gene_Name)) %>%
  mutate(Gene_Name = gsub('%2C', ',' , Gene_Name)) %>% 
  mutate(Gene = gsub('ID=', '', Gene)) %>%
  left_join(operons_RH, by=c('Organism','Gene')) %>%
  mutate(pvalue = coalesce(pvalue, 1)) #replaceing NAs in pvalue column with arbitary '1' to avoid problems downstream
  
  # select genes significant in both experiments 
  # retain only those with abs(log2FC)>1 and concordant sign
alldata <-   alldata %>% 
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
pdf("SB093_SB172_RHoperons_FDR0.1.pdf", height=130, width=27)
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
