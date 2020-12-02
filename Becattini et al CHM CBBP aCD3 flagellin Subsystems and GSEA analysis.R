
### the following codes were used to perform the analyses and generated panels shown in
### to run, this nead loading of the main SB093 RData and access to the 'PATRIC' folder containing updated gff files and annotations
### these codes generate panels shown in Figure 1G, 1F, 2C, 2D, Supplementary Figure 1B, 1F



library(RColorBrewer)
library(DESeq2)
#library(plyr)
library(countToFPKM)
library(tibble)
library(pheatmap)
library(Rsubread)
library(stringr)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(tidyverse)
library(dplyr)
library(viridis)
library(clusterProfiler)
library(reshape2)
library(pheatmap)
library(Pigengene)


setwd('/Users/becattis/Desktop/Kegg_graphs/')
load('/Users/becattis/Desktop/SB093_rerun/SB093_rerun.RData')
#save.image('SB093_baseline_subsystems.RData')
#load('SB093_baseline_subsystems.RData')

##### all subsystem files ####
### note: the following steps were necessary as in order to obtain subsystem annotation, which is a recent feature in PATRIC,
### I had to re-run the annotation, and got some additional gene annotation, which changed the numbering of peg
### therefore I merged the new obtained .gff for each CBBP member, and corresponding subsystem annotation, with the old .gff,
### based on gene position in the genome (start/end/strand) to make sure each gene was being correctly annotated

BP_subs_gff <- read.delim('./PATRIC_Subsystems/Blautia producta SC subsystems.gff', header = F, stringsAsFactors = F)[-1, -c(2,6)] %>%
  `colnames<-` (c("Chr", "Type", "Start", "End", "Strand","Frame", "GeneID")) %>%
  dplyr::mutate(PATRIC.ID = gsub('ID=', '', GeneID)) %>%
  dplyr::mutate(PATRIC.ID = gsub("(;.*)" , "", PATRIC.ID)) %>% 
  select(Start, End, Strand, PATRIC.ID) 
BP_subsystems <- read.delim('./PATRIC_Subsystems/BP_PATRIC_subsystems.txt', header=T, stringsAsFactors=F)
BP_subs_gff2 <- BP_subs_gff %>% left_join(BP_subsystems, by='PATRIC.ID') %>%
  select(-PATRIC.ID) %>% right_join(BP, by=c('Start', 'End', 'Strand')) %>%
  dplyr::mutate(Gene = gsub("(;.*)" , "", GeneID)) %>%
  dplyr::mutate(Gene= gsub('ID=', '', Gene)) %>%
  select(-Chr, -Type, -Start, -End, -Strand, -Frame, -GeneID) %>% mutate(Organims='BP')

#to see what percentage of these genes has a subsystem associated with it. Note: it's basically the same result if annotating with KEGG or COG!
BP_perc_subs <- (BP_subs_gff2 %>% filter(!is.na(Subsystem.Name)) %>% distinct(Gene) %>% nrow())/(BP_subs_gff2 %>% filter(grepl('peg', Gene)) %>%distinct(Gene) %>% nrow())

BS_subs_gff <- read.delim('./PATRIC_Subsystems/Bacteroides sartorii SC subsystems.gff', header = F, stringsAsFactors = F)[-1, -c(2,6)] %>%
  `colnames<-` (c("Chr", "Type", "Start", "End", "Strand","Frame", "GeneID")) %>%
  dplyr::mutate(PATRIC.ID = gsub('ID=', '', GeneID)) %>%
  dplyr::mutate(PATRIC.ID = gsub("(;.*)" , "", PATRIC.ID)) %>% 
  select(Start, End, Strand, PATRIC.ID) 
BS_subsystems <- read.delim('./PATRIC_Subsystems/BS_PATRIC_subsystems.txt', header=T, stringsAsFactors=F)
BS_subs_gff2 <- BS_subs_gff %>% left_join(BS_subsystems, by='PATRIC.ID') %>%
  select(-PATRIC.ID) %>% right_join(BS, by=c('Start', 'End', 'Strand')) %>%
  dplyr::mutate(Gene = gsub("(;.*)" , "", GeneID)) %>%
  dplyr::mutate(Gene= gsub('ID=', '', Gene)) %>%
  select(-Chr, -Type, -Start, -End, -Strand, -Frame, -GeneID) %>% mutate(Organims='BS')

BS_perc_subs <- (BS_subs_gff2 %>% filter(!is.na(Subsystem.Name)) %>% distinct(Gene) %>% nrow())/(BS_subs_gff2 %>% filter(grepl('peg', Gene)) %>% distinct(Gene) %>% nrow())

CB_subs_gff <- read.delim('./PATRIC_Subsystems/Clostridium bolteae SC susbsystem.gff', header = F, stringsAsFactors = F)[-1, -c(2,6)] %>%
  `colnames<-` (c("Chr", "Type", "Start", "End", "Strand","Frame", "GeneID")) %>%
  dplyr::mutate(PATRIC.ID = gsub('ID=', '', GeneID)) %>%
  dplyr::mutate(PATRIC.ID = gsub("(;.*)" , "", PATRIC.ID)) %>% 
  select(Start, End, Strand, PATRIC.ID) 
CB_subsystems <- read.delim('./PATRIC_Subsystems/CB_PATRIC_subsystems.txt', header=T, stringsAsFactors=F)
CB_subs_gff2 <- CB_subs_gff %>% left_join(CB_subsystems, by='PATRIC.ID') %>%
  select(-PATRIC.ID) %>% right_join(CB, by=c('Start', 'End', 'Strand')) %>%
  dplyr::mutate(Gene = gsub("(;.*)" , "", GeneID)) %>%
  dplyr::mutate(Gene= gsub('ID=', '', Gene)) %>%
  select(-Chr, -Type, -Start, -End, -Strand, -Frame, -GeneID) %>% mutate(Organims='CB')

CB_perc_subs <- (CB_subs_gff2 %>% filter(!is.na(Subsystem.Name)) %>% distinct(Gene) %>% nrow())/(CB_subs_gff2 %>% filter(grepl('peg', Gene)) %>% distinct(Gene) %>% nrow())

PD_subs_gff <- read.delim('./PATRIC_Subsystems/Parabacteroides distasonis SC susbystem.gff', header = F, stringsAsFactors = F)[-1, -c(2,6)] %>%
  `colnames<-` (c("Chr", "Type", "Start", "End", "Strand","Frame", "GeneID")) %>%
  dplyr::mutate(PATRIC.ID = gsub('ID=', '', GeneID)) %>%
  dplyr::mutate(PATRIC.ID = gsub("(;.*)" , "", PATRIC.ID)) %>% 
  select(Start, End, Strand, PATRIC.ID) 
PD_subsystems <- read.delim('./PATRIC_Subsystems/PD_PATRIC_subsystems.txt', header=T, stringsAsFactors=F)
PD_subs_gff2 <- PD_subs_gff %>% left_join(PD_subsystems, by='PATRIC.ID') %>%
  select(-PATRIC.ID) %>% right_join(PD, by=c('Start', 'End', 'Strand')) %>%
  dplyr::mutate(Gene = gsub("(;.*)" , "", GeneID)) %>%
  dplyr::mutate(Gene= gsub('ID=', '', Gene)) %>%
  select(-Chr, -Type, -Start, -End, -Strand, -Frame, -GeneID) %>% mutate(Organims='PD')

PD_perc_subs <- (PD_subs_gff2 %>% filter(!is.na(Subsystem.Name)) %>% distinct(Gene) %>% nrow())/(PD_subs_gff2 %>% filter(grepl('peg', Gene)) %>% distinct(Gene) %>% nrow())


### bind all df
all_subs_gff <- bind_rows(BP_subs_gff2, BS_subs_gff2, CB_subs_gff2, PD_subs_gff2) %>%
  mutate(Subclass = ifelse(Subclass == '', Subsystem.Name, Subclass))


#### RPKM ####

BP_Counts <- BP_Counts %>% rownames_to_column('GeneID') %>% mutate(Organism = 'BP') %>% 
  dplyr::rename_all(funs(stringr::str_replace_all( ., "BP_", "" )))

BS_Counts <- BS_Counts %>% rownames_to_column('GeneID') %>% mutate(Organism = 'BS') %>% 
  dplyr::rename_all(funs(stringr::str_replace_all( ., "BS_", "" )))

CB_Counts <- CB_Counts %>% rownames_to_column('GeneID') %>% mutate(Organism = 'CB') %>% 
  dplyr::rename_all(funs(stringr::str_replace_all( ., "CB_", "" )))

PD_Counts <- PD_Counts %>% rownames_to_column('GeneID') %>% mutate(Organism = 'PD') %>% 
  dplyr::rename_all(funs(stringr::str_replace_all( ., "PD_", "" )))

All_Counts <- bind_rows(BP_Counts, BS_Counts, CB_Counts, PD_Counts) 

All_gff <- bind_rows(BP, BS, CB, PD) %>% mutate(Length = End - Start + 1) %>% select(Length, GeneID) 

All_Counts <- dplyr::left_join(All_Counts, All_gff, by='GeneID') #%>%
  #select(-contains(c('RNAlat', 'REPL'))) 


### calculate fpkm by dividing counts for total reads(sum of columns), then for gene length, and then multiply by 10^6 and 10^3 (gene length should be expressed in Kb)
All_counts_fpkm <- All_Counts %>%
  dplyr::group_by(Organism) %>%
  dplyr::mutate_at(vars(2:16), funs(((./sum(.))/Length)*10^9)) %>%
  ungroup() %>%
  dplyr::mutate(Gene = gsub("(;.*)" , '', GeneID)) %>%
  dplyr::mutate(Gene= gsub('ID=', '', Gene))

### melt for convenience
All_counts_fpkm_slim <- All_counts_fpkm %>%
  select(-Length) %>%
  melt(variable.name = 'Sample') %>%
  mutate(Time_point = str_extract(pattern='0h|6h|24h', Sample)) %>%
  mutate(Treatment = str_extract(pattern='flag|aCD3', Sample)) %>%
  mutate(Treatment = replace_na(Treatment, 'none'))

All_counts_fpkm_slim <- All_counts_fpkm_slim %>% left_join(all_subs_gff, by='Gene')

All_counts_fpkm_slim_norm <- All_counts_fpkm_slim %>% 
  filter(!is.na(Class)) %>%
  mutate(rpkm = value) %>%
  group_by(Sample, Organism) %>%
  mutate(rpkm = rpkm/sum(rpkm)) %>%
  ungroup()


### generate color palette
nb.cols <- 27
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)

### graph for RNA reads distribution across CLasses
pdf('SB093_baseline_fpkm.pdf', height = 10, width = 7)
All_counts_fpkm_slim %>%
  filter(!is.na(Class)) %>%
  filter(Time_point == '0h') %>%
  ggplot(aes(x=factor(Organism, levels = c('CB', 'BP', 'BS', 'PD')), y=as.numeric(value))) +
  geom_bar(aes(fill=Class),
           stat="identity",
           position = 'fill') +
  theme_bw(base_size = 14) +
  theme(axis.text.x=element_text(angle=0),
       # legend.position="right",
        strip.text.y=element_text(angle=0),
        strip.text.x=element_text(angle=90),
        aspect.ratio = 3) +
  guides(fill=guide_legend(ncol=1)) +
  scale_fill_manual(values=mycolors) +
  ylab("normalized rpkm")
dev.off()

### graph for CBBP genomes breakdown into CLasses
### this was obtained by giving all genes an arbitrary FPKM of 1
pdf('SB093_baseline_DNA.pdf',  height = 10, width = 7)
All_counts_fpkm_slim %>%
  filter(!is.na(Class)) %>%
  filter(Time_point == '0h') %>%
  mutate(value=1) %>%
  ggplot(aes(x=factor(Organism, levels = c('CB', 'BP', 'BS', 'PD')), y=as.numeric(value))) +
  geom_bar(aes(fill=Class),
           stat="identity",
           position = 'fill') +
  theme_bw(base_size = 14) +
  theme(axis.text.x=element_text(angle=0),
        # legend.position="right",
        strip.text.y=element_text(angle=0),
        strip.text.x=element_text(angle=90),
        aspect.ratio = 3) +
  guides(fill=guide_legend(ncol=1)) +
  scale_fill_manual(values=mycolors) +
  ylab("normalized rpkm")
dev.off()

#### these are zoom ins into specifical Classes reported in Supplementary Figure 1
pdf('SB093_baseline_zoomin_Membrane_transport.pdf',  height = 5, width = 12)
All_counts_fpkm_slim %>%
  filter(Class=='Membrane Transport') %>%
  filter(Time_point == '0h') %>%
  ggplot(aes(x=factor(Organism, levels = c('PD', 'BS', 'BP', 'CB')), y=as.numeric(value))) +
  geom_bar(aes(fill=Subclass),
           stat="identity") +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_text(angle=0),
        # legend.position="right",
        strip.text.y=element_text(angle=0),
        strip.text.x=element_text(angle=90),
        aspect.ratio = 1/4) +
  guides(fill=guide_legend(ncol=1)) +
  scale_fill_manual(values=mycolors) +
  ylab("normalized rpkm") +
  coord_flip()
dev.off()

pdf('SB093_baseline_zoomin_Membrane_transport_DNA.pdf',  height = 10, width = 7)
All_counts_fpkm_slim %>%
  filter(Class=='Membrane Transport') %>%
  filter(Time_point == '0h') %>%
  mutate(value=1) %>%
  ggplot(aes(x=factor(Organism, levels = c('CB', 'BP', 'BS', 'PD')), y=as.numeric(value))) +
  geom_bar(aes(fill=Subclass),
           stat="identity",
           position = 'fill') +
  theme_bw(base_size = 14) +
  theme(axis.text.x=element_text(angle=0),
        # legend.position="right",
        strip.text.y=element_text(angle=0),
        strip.text.x=element_text(angle=90),
        aspect.ratio = 4) +
  guides(fill=guide_legend(ncol=1)) +
  scale_fill_manual(values=mycolors) +
  ylab("normalized rpkm")
dev.off()

pdf('SB093_baseline_zoomin_Stress_response.pdf',  height = 5, width = 12)
All_counts_fpkm_slim %>%
  filter(Class=='Stress Response, Defense and Virulence') %>%
  filter(Time_point == '0h') %>%
  ggplot(aes(x=factor(Organism, levels = c('PD', 'BS', 'BP', 'CB')), y=as.numeric(value))) +
  geom_bar(aes(fill=Subclass),
           stat="identity") +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_text(angle=0),
        # legend.position="right",
        strip.text.y=element_text(angle=0),
        strip.text.x=element_text(angle=90),
        aspect.ratio = 1/4) +
  guides(fill=guide_legend(ncol=1)) +
  scale_fill_manual(values=mycolors) +
  ylab("normalized rpkm") +
  coord_flip()
dev.off()

pdf('SB093_baseline_zoomin_Stress_response_DNA.pdf',  height = 10, width = 7)
All_counts_fpkm_slim %>%
  filter(Class=='Stress Response, Defense and Virulence') %>%
  filter(Time_point == '0h') %>%
  mutate(value=1) %>%
  ggplot(aes(x=factor(Organism, levels = c('CB', 'BP', 'BS', 'PD')), y=as.numeric(value))) +
  geom_bar(aes(fill=Subclass),
           stat="identity",
           position = 'fill') +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_text(angle=0),
        # legend.position="right",
        strip.text.y=element_text(angle=0),
        strip.text.x=element_text(angle=90),
        aspect.ratio = 4) +
  guides(fill=guide_legend(ncol=1)) +
  scale_fill_manual(values=mycolors) +
  ylab("normalized rpkm")
dev.off()

pdf('SB093_baseline_zoomin_Carbohydrates_fpkm.pdf',  height = 5, width = 12)
All_counts_fpkm_slim %>%
  filter(Class=='Carbohydrates') %>%
  filter(Time_point == '0h') %>%
  ggplot(aes(x=factor(Organism, levels = c('PD', 'BS', 'BP', 'CB')), y=as.numeric(value))) +
  geom_bar(aes(fill=Subclass),
           stat="identity") +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_text(angle=0),
        # legend.position="right",
        strip.text.y=element_text(angle=0),
        strip.text.x=element_text(angle=90),
        aspect.ratio = 1/4) +
  guides(fill=guide_legend(ncol=1)) +
  scale_fill_manual(values=mycolors) +
  ylab("normalized rpkm") +
  coord_flip()
dev.off()

pdf('SB093_baseline_zoomin_Carbohydrates_DNA.pdf',  height = 10, width = 7)
All_counts_fpkm_slim %>%
  filter(Class=='Carbohydrates') %>%
  filter(Time_point == '0h') %>%
  mutate(value = 1) %>%
  ggplot(aes(x=factor(Organism, levels = c('CB', 'BP', 'BS', 'PD')), y=as.numeric(value))) +
  geom_bar(aes(fill=Subclass),
           stat="identity",
           position = 'fill') +
  theme_bw(base_size = 14) +
  theme(axis.text.x=element_text(angle=0),
        # legend.position="right",
        strip.text.y=element_text(angle=0),
        strip.text.x=element_text(angle=90),
        aspect.ratio = 4) +
  guides(fill=guide_legend(ncol=1)) +
  scale_fill_manual(values=mycolors) +
  ylab("normalized rpkm")
dev.off()

pdf('SB093_baseline_zoomin_AAs_fpkm.pdf',  height = 5, width = 12)
All_counts_fpkm_slim %>%
  filter(Class=='Amino Acids and Derivatives') %>%
  filter(Time_point == '0h') %>%
  ggplot(aes(x=factor(Organism, levels = c('PD', 'BS', 'BP', 'CB')), y=as.numeric(value))) +
  geom_bar(aes(fill=Subclass),
           stat="identity") +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_text(angle=0),
        # legend.position="right",
        strip.text.y=element_text(angle=0),
        strip.text.x=element_text(angle=90),
        aspect.ratio = 1/4) +
  guides(fill=guide_legend(ncol=1)) +
  scale_fill_manual(values=mycolors) +
  ylab("normalized rpkm") +
  coord_flip()
dev.off()

pdf('SB093_baseline_zoomin_AAs_DNA.pdf',  height = 10, width = 7)
All_counts_fpkm_slim %>%
  filter(Class=='Amino Acids and Derivatives') %>%
  filter(Time_point == '0h') %>%
  mutate(value=1) %>%
  ggplot(aes(x=factor(Organism, levels = c('CB', 'BP', 'BS', 'PD')), y=as.numeric(value))) +
  geom_bar(aes(fill=Subclass),
           stat="identity",
           position = 'fill') +
  theme_bw(base_size = 14) +
  theme(axis.text.x=element_text(angle=0),
        # legend.position="right",
        strip.text.y=element_text(angle=0),
        strip.text.x=element_text(angle=90),
        aspect.ratio = 4) +
  guides(fill=guide_legend(ncol=1)) +
  scale_fill_manual(values=mycolors) +
  ylab("normalized rpkm")
dev.off()

pdf('SB093_baseline_zoomin_Fatty_acids_fpkm.pdf',  height = 5, width = 12)
All_counts_fpkm_slim %>%
  filter(Class=='Fatty Acids, Lipids, and Isoprenoids') %>%
  ggplot(aes(x=factor(Organism, levels = c('PD', 'BS', 'BP', 'CB')), y=as.numeric(value))) +
  geom_bar(aes(fill=Subclass),
           stat="identity") +
  theme_bw(base_size = 16) +
  theme(axis.text.x=element_text(angle=0),
        # legend.position="right",
        strip.text.y=element_text(angle=0),
        strip.text.x=element_text(angle=90),
        aspect.ratio = 1/4) +
  guides(fill=guide_legend(ncol=1)) +
  scale_fill_manual(values=mycolors) +
  ylab("normalized rpkm") +
  coord_flip()
dev.off()

pdf('SB093_baseline_zoomin_Fatty_acids_DNA.pdf',  height = 10, width = 7)
All_counts_fpkm_slim %>%
  filter(Class=='Fatty Acids, Lipids, and Isoprenoids') %>%
  filter(Time_point == '0h') %>%
  mutate(value=1) %>%
  ggplot(aes(x=factor(Organism, levels = c('CB', 'BP', 'BS', 'PD')), y=as.numeric(value))) +
  geom_bar(aes(fill=Subclass),
           stat="identity",
           position = 'fill') +
  theme_bw(base_size = 14) +
  theme(axis.text.x=element_text(angle=0),
        # legend.position="right",
        strip.text.y=element_text(angle=0),
        strip.text.x=element_text(angle=90),
        aspect.ratio = 4) +
  guides(fill=guide_legend(ncol=1)) +
  scale_fill_manual(values=mycolors) +
  ylab("normalized rpkm")
dev.off()



############ GSEA  ##########

### ClusterProfiler was used to perform GSEA using custom annotations (RAST SUbsystems)
### first a term_to_gene and a term_to_name files are generated
### then a rankmetric is chosen: here the log2FC was used but, as reported, the sign-adjusted inverse padj could also be used

#### flagellin

BS_term_to_gene <- all_subs_gff %>% filter(Organims == 'BS') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass, Gene) %>% rename(Term = Subclass, Gene = Gene) 
BS_term_to_name <- all_subs_gff %>% filter(Organims == 'BS') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass) %>% mutate(Term = Subclass) %>% rename(Name = Subclass)

BS_flag6h_vs_0h$rankmetric<-BS_flag6h_vs_0h$log2FoldChange
BS_flag6h_vs_0h<-BS_flag6h_vs_0h[order(BS_flag6h_vs_0h$rankmetric, decreasing=TRUE), ]
BS_flag6h_vs_0h$ID <- BS_flag6h_vs_0h$Gene
BS_flag6h_vs_0h$ID <- gsub("ID=","", BS_flag6h_vs_0h$ID)
BS_flag6h_vs_0h$ID <- gsub("(;.*)" , "", BS_flag6h_vs_0h$ID) #remove anything after;
genelist <- structure(c(BS_flag6h_vs_0h$rankmetric),.Names=BS_flag6h_vs_0h$ID)
BS_GSEA_flag_6h = GSEA(genelist, TERM2GENE=BS_term_to_gene, TERM2NAME=BS_term_to_name, minGSSize=5, pvalueCutoff = 0.05)
BS_GSEA_flag_6h <- as.data.frame(BS_GSEA_flag_6h) %>% mutate(Organims = 'BS') %>% mutate(Treatment = 'flag')

BP_term_to_gene <- all_subs_gff %>% filter(Organims == 'BP') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass, Gene) %>% rename(Term = Subclass, Gene = Gene) 
BP_term_to_name <- all_subs_gff %>% filter(Organims == 'BP') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass) %>% mutate(Term = Subclass) %>% rename(Name = Subclass)

BP_flag6h_vs_0h$rankmetric<-BP_flag6h_vs_0h$log2FoldChange
BP_flag6h_vs_0h<-BP_flag6h_vs_0h[order(BP_flag6h_vs_0h$rankmetric, decreasing=TRUE), ]
BP_flag6h_vs_0h$ID <- BP_flag6h_vs_0h$Gene
BP_flag6h_vs_0h$ID <- gsub("ID=","", BP_flag6h_vs_0h$ID)
BP_flag6h_vs_0h$ID <- gsub("(;.*)" , "", BP_flag6h_vs_0h$ID) #remove anything after;
genelist <- structure(c(BP_flag6h_vs_0h$rankmetric),.Names=BP_flag6h_vs_0h$ID)
BP_GSEA_flag_6h = GSEA(genelist, TERM2GENE=BP_term_to_gene, TERM2NAME=BP_term_to_name, minGSSize=5, pvalueCutoff = 0.05)
BP_GSEA_flag_6h <- as.data.frame(BP_GSEA_flag_6h) %>% mutate(Organims = 'BP') %>% mutate(Treatment = 'flag')

CB_term_to_gene <- all_subs_gff %>% filter(Organims == 'CB') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass, Gene) %>% rename(Term = Subclass, Gene = Gene) 
CB_term_to_name <- all_subs_gff %>% filter(Organims == 'CB') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass) %>% mutate(Term = Subclass) %>% rename(Name = Subclass)

CB_flag6h_vs_0h$rankmetric<-CB_flag6h_vs_0h$log2FoldChange
CB_flag6h_vs_0h<-CB_flag6h_vs_0h[order(CB_flag6h_vs_0h$rankmetric, decreasing=TRUE), ]
CB_flag6h_vs_0h$ID <- CB_flag6h_vs_0h$Gene
CB_flag6h_vs_0h$ID <- gsub("ID=","", CB_flag6h_vs_0h$ID)
CB_flag6h_vs_0h$ID <- gsub("(;.*)" , "", CB_flag6h_vs_0h$ID) #remove anything after;
genelist <- structure(c(CB_flag6h_vs_0h$rankmetric),.Names=CB_flag6h_vs_0h$ID)
CB_GSEA_flag_6h = GSEA(genelist, TERM2GENE=CB_term_to_gene, TERM2NAME=CB_term_to_name, minGSSize=5, pvalueCutoff = 0.05)
CB_GSEA_flag_6h <- as.data.frame(CB_GSEA_flag_6h) %>% mutate(Organims = 'CB') %>% mutate(Treatment = 'flag')

PD_term_to_gene <- all_subs_gff %>% filter(Organims == 'PD') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass, Gene) %>% rename(Term = Subclass, Gene = Gene) 
PD_term_to_name <- all_subs_gff %>% filter(Organims == 'PD') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass) %>% mutate(Term = Subclass) %>% rename(Name = Subclass)

PD_flag6h_vs_0h$rankmetric<-PD_flag6h_vs_0h$log2FoldChange
PD_flag6h_vs_0h<-PD_flag6h_vs_0h[order(PD_flag6h_vs_0h$rankmetric, decreasing=TRUE), ]
PD_flag6h_vs_0h$ID <- PD_flag6h_vs_0h$Gene
PD_flag6h_vs_0h$ID <- gsub("ID=","", PD_flag6h_vs_0h$ID)
PD_flag6h_vs_0h$ID <- gsub("(;.*)" , "", PD_flag6h_vs_0h$ID) #remove anything after;
genelist <- structure(c(PD_flag6h_vs_0h$rankmetric),.Names=PD_flag6h_vs_0h$ID)
PD_GSEA_flag_6h = GSEA(genelist, TERM2GENE=PD_term_to_gene, TERM2NAME=PD_term_to_name, minGSSize=5, pvalueCutoff = 0.05)
PD_GSEA_flag_6h <- as.data.frame(PD_GSEA_flag_6h) %>% mutate(Organims = 'PD') %>% mutate(Treatment = 'flag')

#### aCD3

BS_aCD36h_vs_0h$rankmetric<-BS_aCD36h_vs_0h$log2FoldChange
BS_aCD36h_vs_0h<-BS_aCD36h_vs_0h[order(BS_aCD36h_vs_0h$rankmetric, decreasing=TRUE), ]
BS_aCD36h_vs_0h$ID <- BS_aCD36h_vs_0h$Gene
BS_aCD36h_vs_0h$ID <- gsub("ID=","", BS_aCD36h_vs_0h$ID)
BS_aCD36h_vs_0h$ID <- gsub("(;.*)" , "", BS_aCD36h_vs_0h$ID) #remove anything after;
genelist <- structure(c(BS_aCD36h_vs_0h$rankmetric),.Names=BS_aCD36h_vs_0h$ID)
BS_GSEA_aCD3_6h = GSEA(genelist, TERM2GENE=BS_term_to_gene, TERM2NAME=BS_term_to_name, minGSSize=5, pvalueCutoff = 0.05)
BS_GSEA_aCD3_6h <- as.data.frame(BS_GSEA_aCD3_6h) %>% mutate(Organims = 'BS') %>% mutate(Treatment = 'aCD3')

BP_aCD36h_vs_0h$rankmetric<-BP_aCD36h_vs_0h$log2FoldChange
BP_aCD36h_vs_0h<-BP_aCD36h_vs_0h[order(BP_aCD36h_vs_0h$rankmetric, decreasing=TRUE), ]
BP_aCD36h_vs_0h$ID <- BP_aCD36h_vs_0h$Gene
BP_aCD36h_vs_0h$ID <- gsub("ID=","", BP_aCD36h_vs_0h$ID)
BP_aCD36h_vs_0h$ID <- gsub("(;.*)" , "", BP_aCD36h_vs_0h$ID) #remove anything after;
genelist <- structure(c(BP_aCD36h_vs_0h$rankmetric),.Names=BP_aCD36h_vs_0h$ID)
BP_GSEA_aCD3_6h = GSEA(genelist, TERM2GENE=BP_term_to_gene, TERM2NAME=BP_term_to_name, minGSSize=5, pvalueCutoff = 0.05)
BP_GSEA_aCD3_6h <- as.data.frame(BP_GSEA_aCD3_6h) %>% mutate(Organims = 'BP') %>% mutate(Treatment = 'aCD3')

CB_aCD36h_vs_0h$rankmetric<-CB_aCD36h_vs_0h$log2FoldChange
CB_aCD36h_vs_0h<-CB_aCD36h_vs_0h[order(CB_aCD36h_vs_0h$rankmetric, decreasing=TRUE), ]
CB_aCD36h_vs_0h$ID <- CB_aCD36h_vs_0h$Gene
CB_aCD36h_vs_0h$ID <- gsub("ID=","", CB_aCD36h_vs_0h$ID)
CB_aCD36h_vs_0h$ID <- gsub("(;.*)" , "", CB_aCD36h_vs_0h$ID) #remove anything after;
genelist <- structure(c(CB_aCD36h_vs_0h$rankmetric),.Names=CB_aCD36h_vs_0h$ID)
CB_GSEA_aCD3_6h = GSEA(genelist, TERM2GENE=CB_term_to_gene, TERM2NAME=CB_term_to_name, minGSSize=5, pvalueCutoff = 0.05)
CB_GSEA_aCD3_6h <- as.data.frame(CB_GSEA_aCD3_6h) %>% mutate(Organims = 'CB') %>% mutate(Treatment = 'aCD3')

PD_aCD36h_vs_0h$rankmetric<-PD_aCD36h_vs_0h$log2FoldChange
PD_aCD36h_vs_0h<-PD_aCD36h_vs_0h[order(PD_aCD36h_vs_0h$rankmetric, decreasing=TRUE), ]
PD_aCD36h_vs_0h$ID <- PD_aCD36h_vs_0h$Gene
PD_aCD36h_vs_0h$ID <- gsub("ID=","", PD_aCD36h_vs_0h$ID)
PD_aCD36h_vs_0h$ID <- gsub("(;.*)" , "", PD_aCD36h_vs_0h$ID) #remove anything after;
genelist <- structure(c(PD_aCD36h_vs_0h$rankmetric),.Names=PD_aCD36h_vs_0h$ID)
PD_GSEA_aCD3_6h = GSEA(genelist, TERM2GENE=PD_term_to_gene, TERM2NAME=PD_term_to_name, minGSSize=5, pvalueCutoff = 0.05)
PD_GSEA_aCD3_6h <- as.data.frame(PD_GSEA_aCD3_6h) %>% mutate(Organims = 'PD') %>% mutate(Treatment = 'aCD3')

### generate a df containing all GSEA results
GSEA_all <- bind_rows(BS_GSEA_flag_6h, BS_GSEA_aCD3_6h, BP_GSEA_aCD3_6h, BP_GSEA_flag_6h, CB_GSEA_aCD3_6h, CB_GSEA_flag_6h, PD_GSEA_aCD3_6h, PD_GSEA_flag_6h)

### plot
pdf('SB093_GSEA_flag_all.pdf', useDingbats = F, height=9, width=14)
GSEA_all %>%
  filter(Treatment == 'flag') %>%
  ggplot(aes(x=reorder(ID, NES), y=NES, size = setSize)) +
  ylim(-3,3) +
  geom_point(stat='identity', aes(color=p.adjust)) +
  scale_color_gradient(low='red', high='grey', trans = "log10") +
  theme_bw(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        aspect.ratio = 1/4) +
  facet_wrap(.~factor(Organims, levels=c('CB', 'BP', 'BS', 'PD')), ncol = 1, strip.position="right")
    dev.off()

pdf('SB093_GSEA_aCD3_all.pdf', useDingbats = F, height=9, width=14)
    GSEA_all %>%
      filter(Treatment == 'aCD3') %>%
      ggplot(aes(x=reorder(ID, NES), y=NES, size = setSize)) +
      ylim(-3,3) +
      geom_point(stat='identity', aes(color=p.adjust)) +
      scale_color_gradient(low='red', high='grey', trans = "log10") +
      theme_bw(base_size = 16) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.background = element_blank(),
            aspect.ratio = 1/4) +
      facet_wrap(.~factor(Organims, levels=c('CB', 'BP', 'BS', 'PD')), ncol = 1, strip.position="right")
    dev.off()    

    
############ HEATMAPS ####
    ###### PATRIC HM subsystems do the following  #####
    
    SEED_Classes <- all_subs_gff %>% distinct(Class) %>% filter(!is.na(Class)) 
    
    nb.cols <- 26
    Class <- colorRampPalette(brewer.pal(12, 'Set3'))(nb.cols) #select the n of colors you need, make them into a named vector
    names(Class) = sort(SEED_Classes$Class) #names are the cubsystem classes, sort them so the vector will be arranged by descending alohabetical order
    
    
    my_colour = list(
      Time_point = c('0h'="#868686FF", '6h'="#EFC000FF", '24h'='#7AA6DCFF'),
      Treatment = c(flag='bisque2', aCD3='lightskyblue3', None='plum4'),
      Class = (Class)) # make a composite list with all of your annotations
    
    
    BS_CD3_select_genes_manual <- BS_all_aCD3 %>% filter((padj.x<0.05 & abs(log2FoldChange.x) >1 )|(padj.y<0.05 & abs(log2FoldChange.y) >1)|(padj<0.05 & abs(log2FoldChange)>1)) %>% pull(Gene) %>% unique()
    BS_CD3_mat_manual <- assay(BS_rld)[BS_CD3_select_genes_manual,] 
    BS_CD3_mat_manual <- BS_CD3_mat_manual - rowMeans(BS_CD3_mat_manual)
    BS_CD3_mat_manual <- as.data.frame(BS_CD3_mat_manual) %>% select(-contains('flag')) %>%
      rownames_to_column('Gene') 
    BS_CD3_mat_manual$Gene <- gsub("(;.*)" , "", BS_CD3_mat_manual$Gene)
    BS_CD3_mat_manual$Gene <- gsub("ID=" , "", BS_CD3_mat_manual$Gene)
    
    BS_genelist <- data.frame('Gene'=BS_CD3_mat_manual$Gene)
    
    mat_manual_BS_col_factors <- data.frame(Sample=colnames(BS_CD3_mat_manual)) %>%
      mutate(Time_point=str_extract(Sample, pattern="0h|6h|24h"))  %>% 
      column_to_rownames('Sample')
    
    mat_manual_BS_row_factors = BS_subs_gff2 %>% filter(!is.na(Class)) %>%
      inner_join(BS_genelist, by='Gene') %>%
      select(Class, Gene) %>%
      group_by(Gene) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      group_by(Class)%>%
      filter(n()>1) %>% #get rid of classes with only one gene, else pheatmap.type won't work, needs at least 2 to do internal clustering
      ungroup() %>%
      column_to_rownames('Gene')
    
    BS_CD3_mat_manual <-BS_CD3_mat_manual %>% filter(Gene %in% rownames(mat_manual_BS_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    BS_CD3_mat_manual_ordered <- BS_CD3_mat_manual[rownames(mat_manual_BS_row_factors), ]
    
    gaps <- mat_manual_BS_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    pdf('Subsystems_only_BS_aCD3_heatmap.pdf')
    pheatmap.type(BS_CD3_mat_manual_ordered,  annRow = mat_manual_BS_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color = bluered(50), 
                  annotation_col = mat_manual_BS_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour, 
                  cellwidth = 16, cellheight = 2,
                  main = "B. sartorii aCD3 top DEG" )
    dev.off()
    
    ###### PATRIC HM subsystems including unclassified do the following  #####
    
    SEED_Classes <- all_subs_gff %>% distinct(Class) %>% mutate(Class = replace_na(Class,'Unclassified'))
    
    nb.cols <- 27
    Class <- colorRampPalette(brewer.pal(12, 'Set3'))(nb.cols) #select the n of colors you need, make them into a named vector
    names(Class) = sort(SEED_Classes$Class) #names are the cubsystem classes, sort them so the vector will be arranged by descending alohabetical order
    
    
    my_colour = list(
      Time_point = c('0h'="#868686FF", '6h'="#EFC000FF", '24h'='#7AA6DCFF'),
      Treatment = c(flag='bisque2', aCD3='lightskyblue3', None='plum4'),
      Class = (Class)) # make a composite list with all of your annotations
    
    # to draw HM, select genes that are significant in at least 1 comparison between time points
    
    BS_CD3_select_genes_manual <- BS_all_aCD3 %>% filter((padj.x<0.05 & abs(log2FoldChange.x) >1 )|(padj.y<0.05 & abs(log2FoldChange.y) >1)|(padj<0.05 & abs(log2FoldChange)>1)) %>% pull(Gene) %>% unique()
    BS_CD3_mat_manual <- assay(BS_rld)[BS_CD3_select_genes_manual, c(1,2,3,7,8,9,13,14,15)] # select aCD3 only, use log normalized counts
    BS_CD3_mat_manual <- BS_CD3_mat_manual - rowMeans(BS_CD3_mat_manual) # scale values by subtracting mean
    BS_CD3_mat_manual <- as.data.frame(BS_CD3_mat_manual) %>% select(-contains('flag')) %>%
      rownames_to_column('Gene') 
    BS_CD3_mat_manual$Gene <- gsub("(;.*)" , "", BS_CD3_mat_manual$Gene)
    BS_CD3_mat_manual$Gene <- gsub("ID=" , "", BS_CD3_mat_manual$Gene)
    
    BS_genelist <- data.frame('Gene'=BS_CD3_mat_manual$Gene) # make list of genes of interest
    
    mat_manual_BS_col_factors <- data.frame(Sample=colnames(BS_CD3_mat_manual)) %>% # create factors for HM
      mutate(Time_point=str_extract(Sample, pattern="0h|6h|24h"))  %>% 
      column_to_rownames('Sample')
    
    mat_manual_BS_row_factors = BS_subs_gff2 %>% 
      mutate(Class = replace_na(Class,'Unclassified')) %>%
      inner_join(BS_genelist, by='Gene') %>%
      select(Class, Gene) %>%
      group_by(Gene) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      group_by(Class)%>%
      filter(n()>1) %>% #get rid of classes with only one gene, else pheatmap.type won't work, needs at least 2 to do internal clustering
      ungroup() %>%
      column_to_rownames('Gene')
    
    #filter genes based on list obtained abouve
    BS_CD3_mat_manual <-BS_CD3_mat_manual %>% filter(Gene %in% rownames(mat_manual_BS_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    BS_CD3_mat_manual_ordered <- BS_CD3_mat_manual[rownames(mat_manual_BS_row_factors), ]
    
    gaps <- mat_manual_BS_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many genes for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    ## draw HM including custom annotations
    pdf('SEED_Class_BS_aCD3_heatmap.pdf')
    pheatmap.type(BS_CD3_mat_manual_ordered,  annRow = mat_manual_BS_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color  = colorRampPalette(colors = c("blue", "white", "red"))(17), breaks=seq(-4, 4, by = 0.5), 
                  annotation_col = mat_manual_BS_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "B. sartorii aCD3 top DEG" )
    dev.off()
    
    
    ### same was done for the other CBBP members 
    
    BP_CD3_select_genes_manual <- BP_all_aCD3 %>% filter((padj.x<0.05 & abs(log2FoldChange.x) >1 )|(padj.y<0.05 & abs(log2FoldChange.y) >1)|(padj<0.05 & abs(log2FoldChange)>1)) %>% pull(Gene) %>% unique()
    BP_CD3_mat_manual <- assay(BP_rld)[BP_CD3_select_genes_manual, c(1,2,3,7,8,9,13,14,15)] 
    BP_CD3_mat_manual <- BP_CD3_mat_manual - rowMeans(BP_CD3_mat_manual)
    BP_CD3_mat_manual <- as.data.frame(BP_CD3_mat_manual) %>% select(-contains('flag')) %>%
      rownames_to_column('Gene') 
    BP_CD3_mat_manual$Gene <- gsub("(;.*)" , "", BP_CD3_mat_manual$Gene)
    BP_CD3_mat_manual$Gene <- gsub("ID=" , "", BP_CD3_mat_manual$Gene)
    
    BP_genelist <- data.frame('Gene'=BP_CD3_mat_manual$Gene)
    
    mat_manual_BP_col_factors <- data.frame(Sample=colnames(BP_CD3_mat_manual)) %>%
      mutate(Time_point=str_extract(Sample, pattern="0h|6h|24h"))  %>% 
      column_to_rownames('Sample')
    
    mat_manual_BP_row_factors = BP_subs_gff2 %>% 
      mutate(Class = replace_na(Class,'Unclassified')) %>%
      inner_join(BP_genelist, by='Gene') %>%
      select(Class, Gene) %>%
      group_by(Gene) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      group_by(Class)%>%
      filter(n()>1) %>% #get rid of classes with only one gene, else pheatmap.type won't work, needs at least 2 to do internal clustering
      ungroup() %>%
      column_to_rownames('Gene')
    
    BP_CD3_mat_manual <-BP_CD3_mat_manual %>% filter(Gene %in% rownames(mat_manual_BP_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    BP_CD3_mat_manual_ordered <- BP_CD3_mat_manual[rownames(mat_manual_BP_row_factors), ]
    
    gaps <- mat_manual_BP_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    pdf('SEED_Class_BP_aCD3_heatmap.pdf')
    pheatmap.type(BP_CD3_mat_manual_ordered,  annRow = mat_manual_BP_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color = colorRampPalette(colors = c("blue", "white", "red"))(17), breaks=seq(-4, 4, by = 0.5), 
                  annotation_col = mat_manual_BP_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "B. producta aCD3 top DEG" )
    dev.off()
    
    
    CB_CD3_select_genes_manual <- CB_all_aCD3 %>% filter((padj.x<0.05 & abs(log2FoldChange.x) >1 )|(padj.y<0.05 & abs(log2FoldChange.y) >1)|(padj<0.05 & abs(log2FoldChange)>1)) %>% pull(Gene) %>% unique()
    CB_CD3_mat_manual <- assay(CB_rld)[CB_CD3_select_genes_manual, c(1,2,3,7,8,9,13,14,15)] 
    CB_CD3_mat_manual <- CB_CD3_mat_manual - rowMeans(CB_CD3_mat_manual)
    CB_CD3_mat_manual <- as.data.frame(CB_CD3_mat_manual) %>% select(-contains('flag')) %>%
      rownames_to_column('Gene') 
    CB_CD3_mat_manual$Gene <- gsub("(;.*)" , "", CB_CD3_mat_manual$Gene)
    CB_CD3_mat_manual$Gene <- gsub("ID=" , "", CB_CD3_mat_manual$Gene)
    
    CB_genelist <- data.frame('Gene'=CB_CD3_mat_manual$Gene)
    
    mat_manual_CB_col_factors <- data.frame(Sample=colnames(CB_CD3_mat_manual)) %>%
      mutate(Time_point=str_extract(Sample, pattern="0h|6h|24h"))  %>% 
      column_to_rownames('Sample')
    
    mat_manual_CB_row_factors = CB_subs_gff2 %>% 
      mutate(Class = replace_na(Class,'Unclassified')) %>%
      inner_join(CB_genelist, by='Gene') %>%
      select(Class, Gene) %>%
      group_by(Gene) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      group_by(Class)%>%
      filter(n()>1) %>% #get rid of classes with only one gene, else pheatmap.type won't work, needs at least 2 to do internal clustering
      ungroup() %>%
      column_to_rownames('Gene')
    
    CB_CD3_mat_manual <-CB_CD3_mat_manual %>% filter(Gene %in% rownames(mat_manual_CB_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    CB_CD3_mat_manual_ordered <- CB_CD3_mat_manual[rownames(mat_manual_CB_row_factors), ]
    
    gaps <- mat_manual_CB_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    pdf('SEED_Class_CB_aCD3_heatmap.pdf')
    pheatmap.type(CB_CD3_mat_manual_ordered,  annRow = mat_manual_CB_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color = colorRampPalette(colors = c("blue", "white", "red"))(21), breaks=seq(-5, 5, by = 0.5),
                  annotation_col = mat_manual_CB_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "C. bolteae aCD3 top DEG" )
    dev.off()
    
    
    PD_CD3_select_genes_manual <- PD_all_aCD3 %>% filter((padj.x<0.05 & abs(log2FoldChange.x) >1 )|(padj.y<0.05 & abs(log2FoldChange.y) >1)|(padj<0.05 & abs(log2FoldChange)>1)) %>% pull(Gene) %>% unique()
    PD_CD3_mat_manual <- assay(PD_rld)[PD_CD3_select_genes_manual, c(1,2,3,7,8,9,13,14,15)] 
    PD_CD3_mat_manual <- PD_CD3_mat_manual - rowMeans(PD_CD3_mat_manual)
    PD_CD3_mat_manual <- as.data.frame(PD_CD3_mat_manual) %>% select(-contains('flag')) %>%
      rownames_to_column('Gene') 
    PD_CD3_mat_manual$Gene <- gsub("(;.*)" , "", PD_CD3_mat_manual$Gene)
    PD_CD3_mat_manual$Gene <- gsub("ID=" , "", PD_CD3_mat_manual$Gene)
    
    PD_genelist <- data.frame('Gene'=PD_CD3_mat_manual$Gene)
    
    mat_manual_PD_col_factors <- data.frame(Sample=colnames(PD_CD3_mat_manual)) %>%
      mutate(Time_point=str_extract(Sample, pattern="0h|6h|24h"))  %>% 
      column_to_rownames('Sample')
    
    mat_manual_PD_row_factors = PD_subs_gff2 %>% 
      mutate(Class = replace_na(Class,'Unclassified')) %>%
      inner_join(PD_genelist, by='Gene') %>%
      select(Class, Gene) %>%
      group_by(Gene) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      group_by(Class)%>%
      filter(n()>1) %>% #get rid of classes with only one gene, else pheatmap.type won't work, needs at least 2 to do internal clustering
      ungroup() %>%
      column_to_rownames('Gene')
    
    PD_CD3_mat_manual <-PD_CD3_mat_manual %>% filter(Gene %in% rownames(mat_manual_PD_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    PD_CD3_mat_manual_ordered <- PD_CD3_mat_manual[rownames(mat_manual_PD_row_factors), ]
    
    gaps <- mat_manual_PD_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    pdf('SEED_Class_PD_aCD3_heatmap.pdf')
    pheatmap.type(PD_CD3_mat_manual_ordered,  annRow = mat_manual_PD_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color = colorRampPalette(colors = c("blue", "white", "red"))(17), breaks=seq(-4, 4, by = 0.5), 
                  annotation_col = mat_manual_PD_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "P. distasonis aCD3 top DEG" )
    dev.off()
    
    
    ##### for flagellin treatment 2 separate HM were made to separate unclassified and classified genes 
    
    BP_flag_select_genes_manual <- BP_all_flag %>% filter((padj.x<0.05 & abs(log2FoldChange.x) >1 )|(padj.y<0.05 & abs(log2FoldChange.y) >1)|(padj<0.05 & abs(log2FoldChange)>1)) %>% pull(Gene) %>% unique()
    BP_flag_mat_manual <- assay(BP_rld)[BP_flag_select_genes_manual, c(1,2,3,4,5,6,10,11,12)]
    BP_flag_mat_manual <- BP_flag_mat_manual - rowMeans(BP_flag_mat_manual)
    BP_flag_mat_manual <- as.data.frame(BP_flag_mat_manual) %>% select(-contains('CD3')) %>%
      rownames_to_column('Gene') 
    BP_flag_mat_manual$Gene <- gsub("(;.*)" , "", BP_flag_mat_manual$Gene)
    BP_flag_mat_manual$Gene <- gsub("ID=" , "", BP_flag_mat_manual$Gene)
    
    BP_genelist <- data.frame('Gene'=BP_flag_mat_manual$Gene)
    
    mat_manual_BP_col_factors <- data.frame(Sample=colnames(BP_flag_mat_manual)) %>%
      mutate(Time_point=str_extract(Sample, pattern="0h|6h|24h"))  %>% 
      column_to_rownames('Sample')
    
    mat_manual_BP_row_factors = BP_subs_gff2 %>% 
      mutate(Class = replace_na(Class,'Unclassified')) %>%
      filter(!Class == 'Unclassified') %>% #select Classified only
      inner_join(BP_genelist, by='Gene') %>%
      select(Class, Gene) %>%
      group_by(Gene) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      group_by(Class)%>%
      filter(n()>1) %>% #get rid of classes with only one gene, else pheatmap.type won't work, needs at least 2 to do internal clustering
      ungroup() %>%
      column_to_rownames('Gene')
    
    BP_flag_mat_manual1 <-BP_flag_mat_manual %>% filter(Gene %in% rownames(mat_manual_BP_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    BP_flag_mat_manual_ordered <- BP_flag_mat_manual1[rownames(mat_manual_BP_row_factors), ]
    
    gaps <- mat_manual_BP_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    pdf('SEED_Class_BP_flag_heatmap_classified.pdf')
    pheatmap.type(BP_flag_mat_manual_ordered,  annRow = mat_manual_BP_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color = colorRampPalette(colors = c("blue", "white", "red"))(21), breaks=seq(-5, 5, by = 0.5), 
                  annotation_col = mat_manual_BP_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "B. producta flag top DEG" )
    dev.off()
    
    
    ### this below if for the unclassified
    mat_manual_BP_row_factors = BP_subs_gff2 %>% 
      mutate(Class = replace_na(Class,'Unclassified')) %>%
      filter(Class == 'Unclassified') %>% #select unclassified only
      inner_join(BP_genelist, by='Gene') %>%
      select(Class, Gene) %>%
      group_by(Gene) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      group_by(Class)%>%
      filter(n()>1) %>% #get rid of classes with only one gene, else pheatmap.type won't work, needs at least 2 to do internal clustering
      ungroup() %>%
      column_to_rownames('Gene')
    
    BP_flag_mat_manual1 <-BP_flag_mat_manual %>% filter(Gene %in% rownames(mat_manual_BP_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    BP_flag_mat_manual_ordered <- BP_flag_mat_manual1[rownames(mat_manual_BP_row_factors), ]
    
    gaps <- mat_manual_BP_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    pdf('SEED_Class_BP_flag_heatmap_UNclassified.pdf')
    pheatmap.type(BP_flag_mat_manual_ordered,  annRow = mat_manual_BP_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color = colorRampPalette(colors = c("blue", "white", "red"))(17), breaks=seq(-4, 4, by = 0.5), 
                  annotation_col = mat_manual_BP_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "B. producta flag top DEG" )
    dev.off()
    
    
    
    CB_flag_select_genes_manual <- CB_all_flag %>% filter((padj.x<0.05 & abs(log2FoldChange.x) >1 )|(padj.y<0.05 & abs(log2FoldChange.y) >1)|(padj<0.05 & abs(log2FoldChange)>1)) %>% pull(Gene) %>% unique()
    CB_flag_mat_manual <- assay(CB_rld)[CB_flag_select_genes_manual, c(1,2,3,4,5,6,10,11,12) ]
    CB_flag_mat_manual <- CB_flag_mat_manual - rowMeans(CB_flag_mat_manual)
    CB_flag_mat_manual <- as.data.frame(CB_flag_mat_manual) %>% select(-contains('CD3')) %>%
      rownames_to_column('Gene') 
    CB_flag_mat_manual$Gene <- gsub("(;.*)" , "", CB_flag_mat_manual$Gene)
    CB_flag_mat_manual$Gene <- gsub("ID=" , "", CB_flag_mat_manual$Gene)
    
    CB_genelist <- data.frame('Gene'=CB_flag_mat_manual$Gene)
    
    mat_manual_CB_col_factors <- data.frame(Sample=colnames(CB_flag_mat_manual)) %>%
      mutate(Time_point=str_extract(Sample, pattern="0h|6h|24h"))  %>% 
      column_to_rownames('Sample')
    
    mat_manual_CB_row_factors = CB_subs_gff2 %>% 
      mutate(Class = replace_na(Class,'Unclassified')) %>%
      filter(!Class == 'Unclassified') %>% #select Classified only
      inner_join(CB_genelist, by='Gene') %>%
      select(Class, Gene) %>%
      group_by(Gene) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      group_by(Class)%>%
      filter(n()>1) %>% #get rid of classes with only one gene, else pheatmap.type won't work, needs at least 2 to do internal clustering
      ungroup() %>%
      column_to_rownames('Gene')
    
    CB_flag_mat_manual1 <-CB_flag_mat_manual %>% filter(Gene %in% rownames(mat_manual_CB_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    CB_flag_mat_manual_ordered <- CB_flag_mat_manual1[rownames(mat_manual_CB_row_factors), ]
    
    gaps <- mat_manual_CB_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    pdf('SEED_Class_CB_flag_heatmap_classified.pdf')
    pheatmap.type(CB_flag_mat_manual_ordered,  annRow = mat_manual_CB_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color = colorRampPalette(colors = c("blue", "white", "red"))(21), breaks=seq(-5, 5, by = 0.5), 
                  annotation_col = mat_manual_CB_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "C. bolteae flag top DEG" )
    dev.off()
    
    
    ### this below if for the unclassified
    mat_manual_CB_row_factors = CB_subs_gff2 %>% 
      mutate(Class = replace_na(Class,'Unclassified')) %>%
      filter(Class == 'Unclassified') %>% #select unclassified only
      inner_join(CB_genelist, by='Gene') %>%
      select(Class, Gene) %>%
      group_by(Gene) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      group_by(Class)%>%
      filter(n()>1) %>% #get rid of classes with only one gene, else pheatmap.type won't work, needs at least 2 to do internal clustering
      ungroup() %>%
      column_to_rownames('Gene')
    
    CB_flag_mat_manual1 <-CB_flag_mat_manual %>% filter(Gene %in% rownames(mat_manual_CB_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    CB_flag_mat_manual_ordered <- CB_flag_mat_manual1[rownames(mat_manual_CB_row_factors), ]
    
    gaps <- mat_manual_CB_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    pdf('SEED_Class_CB_flag_heatmap_UNclassified.pdf')
    pheatmap.type(CB_flag_mat_manual_ordered,  annRow = mat_manual_CB_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color = colorRampPalette(colors = c("blue", "white", "red"))(21), breaks=seq(-5, 5, by = 0.5), 
                  annotation_col = mat_manual_CB_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "C. bolteae flag top DEG" )
    dev.off()
    
    
    PD_flag_select_genes_manual <- PD_all_flag %>% filter((padj.x<0.05 & abs(log2FoldChange.x) >1 )|(padj.y<0.05 & abs(log2FoldChange.y) >1)|(padj<0.05 & abs(log2FoldChange)>1)) %>% pull(Gene) %>% unique()
    PD_flag_mat_manual <- assay(PD_rld)[PD_flag_select_genes_manual, c(1,2,3,4,5,6,10,11,12)]
    PD_flag_mat_manual <- PD_flag_mat_manual - rowMeans(PD_flag_mat_manual)
    PD_flag_mat_manual <- as.data.frame(PD_flag_mat_manual) %>% select(-contains('CD3')) %>%
      rownames_to_column('Gene') 
    PD_flag_mat_manual$Gene <- gsub("(;.*)" , "", PD_flag_mat_manual$Gene)
    PD_flag_mat_manual$Gene <- gsub("ID=" , "", PD_flag_mat_manual$Gene)
    
    PD_genelist <- data.frame('Gene'=PD_flag_mat_manual$Gene)
    
    mat_manual_PD_col_factors <- data.frame(Sample=colnames(PD_flag_mat_manual)) %>%
      mutate(Time_point=str_extract(Sample, pattern="0h|6h|24h"))  %>% 
      column_to_rownames('Sample')
    
    mat_manual_PD_row_factors = PD_subs_gff2 %>% 
      mutate(Class = replace_na(Class,'Unclassified')) %>%
      filter(!Class == 'Unclassified') %>% #select Classified only
      inner_join(PD_genelist, by='Gene') %>%
      select(Class, Gene) %>%
      group_by(Gene) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      group_by(Class)%>%
      filter(n()>1) %>% #get rid of classes with only one gene, else pheatmap.type won't work, needs at least 2 to do internal clustering
      ungroup() %>%
      column_to_rownames('Gene')
    
    PD_flag_mat_manual1 <-PD_flag_mat_manual %>% filter(Gene %in% rownames(mat_manual_PD_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    PD_flag_mat_manual_ordered <- PD_flag_mat_manual1[rownames(mat_manual_PD_row_factors), ]
    
    gaps <- mat_manual_PD_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    pdf('SEED_Class_PD_flag_heatmap_classified.pdf')
    pheatmap.type(PD_flag_mat_manual_ordered,  annRow = mat_manual_PD_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color = colorRampPalette(colors = c("blue", "white", "red"))(13), breaks=seq(-3, 3, by = 0.5), 
                  annotation_col = mat_manual_PD_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "P. distasonis flag top DEG" )
    dev.off()
    ### this below if for the unclassified
    mat_manual_PD_row_factors = PD_subs_gff2 %>% 
      mutate(Class = replace_na(Class,'Unclassified')) %>%
      filter(Class == 'Unclassified') %>% #select unclassified only
      inner_join(PD_genelist, by='Gene') %>%
      select(Class, Gene) %>%
      group_by(Gene) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      group_by(Class)%>%
      filter(n()>1) %>% #get rid of classes with only one gene, else pheatmap.type won't work, needs at least 2 to do internal clustering
      ungroup() %>%
      column_to_rownames('Gene')
    
    PD_flag_mat_manual1 <-PD_flag_mat_manual %>% filter(Gene %in% rownames(mat_manual_PD_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    PD_flag_mat_manual_ordered <- PD_flag_mat_manual1[rownames(mat_manual_PD_row_factors), ]
    
    gaps <- mat_manual_PD_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    pdf('SEED_Class_PD_flag_heatmap_UNclassified.pdf')
    pheatmap.type(PD_flag_mat_manual_ordered,  annRow = mat_manual_PD_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color = colorRampPalette(colors = c("blue", "white", "red"))(17), breaks=seq(-4, 4, by = 0.5), 
                  annotation_col = mat_manual_PD_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "P. distasonis flag top DEG" )
    dev.off()
    
    BS_flag_select_genes_manual <- BS_all_flag %>% filter((padj.x<0.05 & abs(log2FoldChange.x) >1 )|(padj.y<0.05 & abs(log2FoldChange.y) >1)|(padj<0.05 & abs(log2FoldChange)>1)) %>% pull(Gene) %>% unique()
    BS_flag_mat_manual <- assay(BS_rld)[BS_flag_select_genes_manual, c(1,2,3,4,5,6,10,11,12)]
    BS_flag_mat_manual <- BS_flag_mat_manual - rowMeans(BS_flag_mat_manual)
    BS_flag_mat_manual <- as.data.frame(BS_flag_mat_manual) %>% select(-contains('CD3')) %>%
      rownames_to_column('Gene') 
    BS_flag_mat_manual$Gene <- gsub("(;.*)" , "", BS_flag_mat_manual$Gene)
    BS_flag_mat_manual$Gene <- gsub("ID=" , "", BS_flag_mat_manual$Gene)
    
    BS_genelist <- data.frame('Gene'=BS_flag_mat_manual$Gene)
    
    mat_manual_BS_col_factors <- data.frame(Sample=colnames(BS_flag_mat_manual)) %>%
      mutate(Time_point=str_extract(Sample, pattern="0h|6h|24h"))  %>% 
      column_to_rownames('Sample')
    
    mat_manual_BS_row_factors = BS_subs_gff2 %>% 
      mutate(Class = replace_na(Class,'Unclassified')) %>%
      filter(!Class == 'Unclassified') %>% #select Classified only
      inner_join(BS_genelist, by='Gene') %>%
      select(Class, Gene) %>%
      group_by(Gene) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      group_by(Class)%>%
      filter(n()>1) %>% #get rid of classes with only one gene, else pheatmap.type won't work, needs at least 2 to do internal clustering
      ungroup() %>%
      column_to_rownames('Gene')
    
    BS_flag_mat_manual1 <-BS_flag_mat_manual %>% filter(Gene %in% rownames(mat_manual_BS_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    BS_flag_mat_manual_ordered <- BS_flag_mat_manual1[rownames(mat_manual_BS_row_factors), ]
    
    gaps <- mat_manual_BS_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    pdf('SEED_Class_BS_flag_heatmap_classified.pdf')
    pheatmap.type(BS_flag_mat_manual_ordered,  annRow = mat_manual_BS_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color = colorRampPalette(colors = c("blue", "white", "red"))(21), breaks=seq(-5, 5, by = 0.5), 
                  annotation_col = mat_manual_BS_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "B. sartorii flag top DEG" )
    dev.off()
    
    unique(mat_manual_BS_row_factors$Class)
    my_colour[]
    ### this below if for the unclassified
    mat_manual_BS_row_factors = BS_subs_gff2 %>% 
      mutate(Class = replace_na(Class,'Unclassified')) %>%
      filter(Class == 'Unclassified') %>% #select unclassified only
      inner_join(BS_genelist, by='Gene') %>%
      select(Class, Gene) %>%
      group_by(Gene) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      group_by(Class)%>%
      filter(n()>1) %>% #get rid of classes with only one gene, else pheatmap.type won't work, needs at least 2 to do internal clustering
      ungroup() %>%
      column_to_rownames('Gene')
    
    BS_flag_mat_manual1 <-BS_flag_mat_manual %>% filter(Gene %in% rownames(mat_manual_BS_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    BS_flag_mat_manual_ordered <- BS_flag_mat_manual1[rownames(mat_manual_BS_row_factors), ]
    
    gaps <- mat_manual_BS_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    pdf('SEED_Class_BS_flag_heatmap_UNclassified.pdf')
    pheatmap.type(BS_flag_mat_manual_ordered,  annRow = mat_manual_BS_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color = colorRampPalette(colors = c("blue", "white", "red"))(21), breaks=seq(-5, 5, by = 0.5), 
                  annotation_col = mat_manual_BS_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "B. sartorii flag top DEG" )
    dev.off()
    
    