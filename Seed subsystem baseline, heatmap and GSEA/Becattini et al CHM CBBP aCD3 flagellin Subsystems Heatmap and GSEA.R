### Becattini et al., CHM 2021
### these codes generate panels shown in Figure 1G, 1F, 2C, 2D, Supplementary Figure 1B, 1F

library(RColorBrewer)
library(DESeq2)
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


setwd('/Users/becattis/Desktop/CHM revision/Scripts to be uploaded/CBBP-transcriptomics/Seed subsystem baseline, heatmap and GSEA/')


##### all subsystem files ####

# read in gff file containing all SEED subsystem annotations for CBBP bacteria, based on individual gff files obtained through PATRIC

all_subs_gff <- read.delim2('all_subs_gff.txt')

#### calculate RPKM starting from counts####

# read in counts as calculated by FeatureCounts in the RNAseq script, to which gene length was added (calculated as end-start+1 in the gff)

All_Counts <- read.delim2('All_counts.txt')

# calculate fpkm by dividing counts for total reads(sum of columns), then for gene length, and then multiply by 10^6 and 10^3 (gene length should be expressed in Kb)
All_counts_fpkm <- All_Counts %>%
  dplyr::group_by(Organism) %>%
  dplyr::mutate_at(vars(2:16), funs(((./sum(.))/Length)*10^9)) %>%
  ungroup() %>%
  dplyr::mutate(Gene = gsub("(;.*)" , '', GeneID)) %>%
  dplyr::mutate(Gene= gsub('ID=', '', Gene))

# melt for convenience
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

# generate color palette
nb.cols <- 27
mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)

# graph for RNA reads distribution across CLasses

Figure_1B_right <- All_counts_fpkm_slim %>%
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

# graph for CBBP genomes breakdown into CLasses
# this was obtained by giving all genes an arbitrary FPKM of 1

Figure_1B_left <- All_counts_fpkm_slim %>%
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


#these are zoom ins into specifical Classes reported in Supplementary Figure 1

Membrane_transport <- All_counts_fpkm_slim %>%
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

Stress_responses <- All_counts_fpkm_slim %>%
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

Carbohydrates <- All_counts_fpkm_slim %>%
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

AA <- All_counts_fpkm_slim %>%
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

Fatty_acids <- All_counts_fpkm_slim %>%
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



############ GSEA  ##########

### ClusterProfiler was used to perform GSEA using custom annotations (SEED Subsystems)
### first a term_to_gene and a term_to_name files are generated
### then a rankmetric is chosen: here the log2FC was used but, as reported, the sign-adjusted inverse padj could also be used, obtaining comparable results

# read in df with all results from SB093 (a-CD3/flagellin in ampicillin-treated, CBBP-reconstituted mice)

SB093_complete_df <- read.delim2('SB093_complete_df.txt') %>%
  mutate(log2FoldChange = as.numeric(log2FoldChange), padj = as.numeric(padj))

#### flagellin

#BS

BS_term_to_gene <- all_subs_gff %>% filter(Organims == 'BS') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass, Gene) %>% rename(Term = Subclass) %>% mutate(Gene = gsub('fig\\|671267.8.','', Gene))
BS_term_to_name <- all_subs_gff %>% filter(Organims == 'BS') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass) %>% mutate(Term = Subclass) %>% rename(Name = Subclass)

BS_flag_6h <- SB093_complete_df %>% 
  filter(Organism == 'BS' & Time_point == '6h' & Treatment == 'flag') %>%
  mutate(Gene = gsub('ID=fig671267.8.','', Gene_number)) %>%
  mutate(log2FoldChange = as.numeric(log2FoldChange)) %>%
  arrange(desc(log2FoldChange))
  
genelist <- structure(c(BS_flag_6h$log2FoldChange),.Names=BS_flag_6h$Gene)

BS_GSEA_flag_6h = GSEA(genelist, TERM2GENE=BS_term_to_gene, TERM2NAME=BS_term_to_name, minGSSize=5, pvalueCutoff = 0.05)
BS_GSEA_flag_6h <- as.data.frame(BS_GSEA_flag_6h) %>% mutate(Organims = 'BS') %>% mutate(Treatment = 'flag')

#BP

BP_term_to_gene <- all_subs_gff %>% filter(Organims == 'BP') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass, Gene) %>% rename(Term = Subclass) %>% mutate(Gene = gsub('fig\\|33035.10.','', Gene))
BP_term_to_name <- all_subs_gff %>% filter(Organims == 'BP') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass) %>% mutate(Term = Subclass) %>% rename(Name = Subclass)

BP_flag_6h <- SB093_complete_df %>% 
  filter(Organism == 'BP' & Time_point == '6h' & Treatment == 'flag') %>%
  mutate(Gene = gsub('fig33035.10.','', Gene_number)) %>%
  mutate(log2FoldChange = as.numeric(log2FoldChange)) %>%
  arrange(desc(log2FoldChange))

genelist <- structure(c(BP_flag_6h$log2FoldChange),.Names=BP_flag_6h$Gene)

BP_GSEA_flag_6h = GSEA(genelist, TERM2GENE=BP_term_to_gene, TERM2NAME=BP_term_to_name, minGSSize=5, pvalueCutoff = 0.05)
BP_GSEA_flag_6h <- as.data.frame(BP_GSEA_flag_6h) %>% mutate(Organims = 'BP') %>% mutate(Treatment = 'flag')

# CB

CB_term_to_gene <- all_subs_gff %>% filter(Organims == 'CB') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass, Gene) %>% rename(Term = Subclass) %>% mutate(Gene = gsub('fig\\|208479.6.','', Gene))
CB_term_to_name <- all_subs_gff %>% filter(Organims == 'CB') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass) %>% mutate(Term = Subclass) %>% rename(Name = Subclass)

CB_flag_6h <- SB093_complete_df %>% 
  filter(Organism == 'CB' & Time_point == '6h' & Treatment == 'flag') %>%
  mutate(Gene = gsub('ID=fig208479.6.','', Gene_number)) %>%
  mutate(log2FoldChange = as.numeric(log2FoldChange)) %>%
  arrange(desc(log2FoldChange))

genelist <- structure(c(CB_flag_6h$log2FoldChange),.Names=CB_flag_6h$Gene)

CB_GSEA_flag_6h = GSEA(genelist, TERM2GENE=CB_term_to_gene, TERM2NAME=CB_term_to_name, minGSSize=5, pvalueCutoff = 0.05)
CB_GSEA_flag_6h <- as.data.frame(CB_GSEA_flag_6h) %>% mutate(Organims = 'CB') %>% mutate(Treatment = 'flag')


#PD

PD_term_to_gene <- all_subs_gff %>% filter(Organims == 'PD') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass, Gene) %>% rename(Term = Subclass) %>% mutate(Gene = gsub('fig\\|823.92.','', Gene))
PD_term_to_name <- all_subs_gff %>% filter(Organims == 'PD') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass) %>% mutate(Term = Subclass) %>% rename(Name = Subclass)

PD_flag_6h <- SB093_complete_df %>% 
  filter(Organism == 'PD' & Time_point == '6h' & Treatment == 'flag') %>%
  mutate(Gene = gsub('ID=fig823.92.','', Gene_number)) %>%
  mutate(log2FoldChange = as.numeric(log2FoldChange)) %>%
  arrange(desc(log2FoldChange))

genelist <- structure(c(PD_flag_6h$log2FoldChange),.Names=PD_flag_6h$Gene)

PD_GSEA_flag_6h = GSEA(genelist, TERM2GENE=PD_term_to_gene, TERM2NAME=PD_term_to_name, minGSSize=5, pvalueCutoff = 0.05)
PD_GSEA_flag_6h <- as.data.frame(PD_GSEA_flag_6h) %>% mutate(Organims = 'PD') %>% mutate(Treatment = 'flag')


#### aCD3

#BS

BS_term_to_gene <- all_subs_gff %>% filter(Organims == 'BS') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass, Gene) %>% rename(Term = Subclass) %>% mutate(Gene = gsub('fig\\|671267.8.','', Gene))
BS_term_to_name <- all_subs_gff %>% filter(Organims == 'BS') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass) %>% mutate(Term = Subclass) %>% rename(Name = Subclass)

BS_aCD3_6h <- SB093_complete_df %>% 
  filter(Organism == 'BS' & Time_point == '6h' & Treatment == 'aCD3') %>%
  mutate(Gene = gsub('ID=fig671267.8.','', Gene_number)) %>%
  mutate(log2FoldChange = as.numeric(log2FoldChange)) %>%
  arrange(desc(log2FoldChange))

genelist <- structure(c(BS_aCD3_6h$log2FoldChange),.Names=BS_aCD3_6h$Gene)

BS_GSEA_aCD3_6h = GSEA(genelist, TERM2GENE=BS_term_to_gene, TERM2NAME=BS_term_to_name, minGSSize=5, pvalueCutoff = 0.05)
BS_GSEA_aCD3_6h <- as.data.frame(BS_GSEA_aCD3_6h) %>% mutate(Organims = 'BS') %>% mutate(Treatment = 'aCD3')

#BP

BP_term_to_gene <- all_subs_gff %>% filter(Organims == 'BP') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass, Gene) %>% rename(Term = Subclass) %>% mutate(Gene = gsub('fig\\|33035.10.','', Gene))
BP_term_to_name <- all_subs_gff %>% filter(Organims == 'BP') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass) %>% mutate(Term = Subclass) %>% rename(Name = Subclass)

BP_aCD3_6h <- SB093_complete_df %>% 
  filter(Organism == 'BP' & Time_point == '6h' & Treatment == 'aCD3') %>%
  mutate(Gene = gsub('fig33035.10.','', Gene_number)) %>%
  mutate(log2FoldChange = as.numeric(log2FoldChange)) %>%
  arrange(desc(log2FoldChange))

genelist <- structure(c(BP_aCD3_6h$log2FoldChange),.Names=BP_aCD3_6h$Gene)

BP_GSEA_aCD3_6h = GSEA(genelist, TERM2GENE=BP_term_to_gene, TERM2NAME=BP_term_to_name, minGSSize=5, pvalueCutoff = 0.05)
BP_GSEA_aCD3_6h <- as.data.frame(BP_GSEA_aCD3_6h) %>% mutate(Organims = 'BP') %>% mutate(Treatment = 'aCD3')

# CB

CB_term_to_gene <- all_subs_gff %>% filter(Organims == 'CB') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass, Gene) %>% rename(Term = Subclass) %>% mutate(Gene = gsub('fig\\|208479.6.','', Gene))
CB_term_to_name <- all_subs_gff %>% filter(Organims == 'CB') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass) %>% mutate(Term = Subclass) %>% rename(Name = Subclass)

CB_aCD3_6h <- SB093_complete_df %>% 
  filter(Organism == 'CB' & Time_point == '6h' & Treatment == 'aCD3') %>%
  mutate(Gene = gsub('ID=fig208479.6.','', Gene_number)) %>%
  mutate(log2FoldChange = as.numeric(log2FoldChange)) %>%
  arrange(desc(log2FoldChange))

genelist <- structure(c(CB_aCD3_6h$log2FoldChange),.Names=CB_aCD3_6h$Gene)

CB_GSEA_aCD3_6h = GSEA(genelist, TERM2GENE=CB_term_to_gene, TERM2NAME=CB_term_to_name, minGSSize=5, pvalueCutoff = 0.05)
CB_GSEA_aCD3_6h <- as.data.frame(CB_GSEA_aCD3_6h) %>% mutate(Organims = 'CB') %>% mutate(Treatment = 'aCD3')


#PD

PD_term_to_gene <- all_subs_gff %>% filter(Organims == 'PD') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass, Gene) %>% rename(Term = Subclass) %>% mutate(Gene = gsub('fig\\|823.92.','', Gene))
PD_term_to_name <- all_subs_gff %>% filter(Organims == 'PD') %>% 
  filter(!is.na(Subclass)) %>% select(Subclass) %>% mutate(Term = Subclass) %>% rename(Name = Subclass)

PD_aCD3_6h <- SB093_complete_df %>% 
  filter(Organism == 'PD' & Time_point == '6h' & Treatment == 'aCD3') %>%
  mutate(Gene = gsub('ID=fig823.92.','', Gene_number)) %>%
  mutate(log2FoldChange = as.numeric(log2FoldChange)) %>%
  arrange(desc(log2FoldChange))

genelist <- structure(c(PD_aCD3_6h$log2FoldChange),.Names=PD_aCD3_6h$Gene)

PD_GSEA_aCD3_6h = GSEA(genelist, TERM2GENE=PD_term_to_gene, TERM2NAME=PD_term_to_name, minGSSize=5, pvalueCutoff = 0.05)
PD_GSEA_aCD3_6h <- as.data.frame(PD_GSEA_aCD3_6h) %>% mutate(Organims = 'PD') %>% mutate(Treatment = 'aCD3')


### generate a df containing all GSEA results
GSEA_all <- bind_rows(BS_GSEA_flag_6h, BS_GSEA_aCD3_6h, BP_GSEA_aCD3_6h, BP_GSEA_flag_6h, CB_GSEA_aCD3_6h, CB_GSEA_flag_6h, PD_GSEA_aCD3_6h, PD_GSEA_flag_6h)

### plot

Figure_1G <- GSEA_all %>%
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


Supplementary_Figure_2F <-  GSEA_all %>%
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
   
    
############ HEATMAPS ####

    ###### PATRIC HM subsystems 
    
    # read in normalized counts obtained with the rlog function in Deseq2, this will be used to build the heatmap
    Norm_rld_counts <- read.delim2('Norm_rld_counts.txt', dec = '.')
    
    SEED_Classes <- all_subs_gff %>% distinct(Class) %>% mutate(Class = replace_na(Class,'Unclassified'))
    
    nb.cols <- 27
    Class <- colorRampPalette(brewer.pal(12, 'Set3'))(nb.cols) #select the n of colors you need, make them into a named vector
    names(Class) = sort(SEED_Classes$Class) #names are the subsystem classes, sort them so the vector will be arranged by descending alphabetical order
    
    
    my_colour = list(
      Time_point = c('0h'="#868686FF", '6h'="#EFC000FF", '24h'='#7AA6DCFF'),
      Treatment = c(flag='bisque2', aCD3='lightskyblue3', None='plum4'),
      Class = (Class)) # make a composite list with all of your annotations
    
    
    # to draw HM, select genes that are significant at 6h or 24h, then make a matrix with centered value and annotate with subsystems
    
    #### a-CD3 
    
    # BS
    
    BS_CD3_select_genes_manual <- SB093_complete_df %>% 
      filter(Organism=='BS' & Treatment == 'aCD3' & Time_point =='6h') %>%
      left_join(SB093_complete_df %>% 
                  filter(Organism=='BS' & Treatment == 'aCD3' & Time_point =='24h'), by='Gene') %>% 
      filter((padj.x<0.05 & abs(log2FoldChange.x) >1 )|(padj.y<0.05 & abs(log2FoldChange.y) >1)) %>% pull(Gene) %>% unique()
    
    BS_CD3_mat_manual <- Norm_rld_counts %>% 
      filter(Organism == 'BS') %>% select(-contains('flag'), -Organism) %>%
      filter(Gene %in% BS_CD3_select_genes_manual) %>%
      mutate(Gene = gsub("(;.*)" , "", Gene)) %>%
      mutate(Gene = gsub("ID=" , "", Gene)) %>%
      column_to_rownames('Gene') 
   
    BS_CD3_mat_manual <- BS_CD3_mat_manual - rowMeans(BS_CD3_mat_manual) # scale values by subtracting mean for each row=gene
    BS_CD3_mat_manual <- as.data.frame(BS_CD3_mat_manual) %>% select(-contains('flag')) %>%
      rownames_to_column('Gene') 
    BS_CD3_mat_manual$Gene <- gsub("(;.*)" , "", BS_CD3_mat_manual$Gene)
    BS_CD3_mat_manual$Gene <- gsub("ID=" , "", BS_CD3_mat_manual$Gene)
    
    BS_genelist <- data.frame('Gene'=BS_CD3_mat_manual$Gene) # make list of genes of interest
    
    mat_manual_BS_col_factors <- data.frame(Sample=colnames(BS_CD3_mat_manual)) %>% # create factors for HM
      mutate(Time_point=str_extract(Sample, pattern="0h|6h|24h"))  %>% 
      column_to_rownames('Sample')
    
    mat_manual_BS_row_factors = all_subs_gff %>%
      filter(Organism == 'BS') %>%
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
    
    #filter genes based on list obtained above
    BS_CD3_mat_manual <-BS_CD3_mat_manual %>% filter(Gene %in% rownames(mat_manual_BS_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    BS_CD3_mat_manual_ordered <- BS_CD3_mat_manual[rownames(mat_manual_BS_row_factors), ]
    
    gaps <- mat_manual_BS_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many genes for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    ## draw HM including custom annotations
    # Figure 2D
    pheatmap.type(BS_CD3_mat_manual_ordered,  annRow = mat_manual_BS_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color  = colorRampPalette(colors = c("blue", "white", "red"))(17), breaks=seq(-4, 4, by = 0.5), 
                  annotation_col = mat_manual_BS_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "B. sartorii aCD3 top DEG" )
  
    # BP
    
    BP_CD3_select_genes_manual <- SB093_complete_df %>% 
      filter(Organism=='BP' & Treatment == 'aCD3' & Time_point =='6h') %>%
      left_join(SB093_complete_df %>% 
                  filter(Organism=='BP' & Treatment == 'aCD3' & Time_point =='24h'), by='Gene') %>% 
      filter((padj.x<0.05 & abs(log2FoldChange.x) >1 )|(padj.y<0.05 & abs(log2FoldChange.y) >1)) %>% pull(Gene) %>% unique()
    
    BP_CD3_mat_manual <- Norm_rld_counts %>% 
      filter(Organism == 'BP') %>% select(-contains('flag'), -Organism) %>%
      filter(Gene %in% BP_CD3_select_genes_manual) %>%
      mutate(Gene = gsub("(;.*)" , "", Gene)) %>%
      mutate(Gene = gsub("ID=" , "", Gene)) %>%
      column_to_rownames('Gene') 
    
    BP_CD3_mat_manual <- BP_CD3_mat_manual - rowMeans(BP_CD3_mat_manual) # scale values by subtracting mean for each row=gene
    BP_CD3_mat_manual <- as.data.frame(BP_CD3_mat_manual) %>% select(-contains('flag')) %>%
      rownames_to_column('Gene') 
    BP_CD3_mat_manual$Gene <- gsub("(;.*)" , "", BP_CD3_mat_manual$Gene)
    BP_CD3_mat_manual$Gene <- gsub("ID=" , "", BP_CD3_mat_manual$Gene)
    
    BP_genelist <- data.frame('Gene'=BP_CD3_mat_manual$Gene) # make list of genes of interest
    
    mat_manual_BP_col_factors <- data.frame(Sample=colnames(BP_CD3_mat_manual)) %>% # create factors for HM
      mutate(Time_point=str_extract(Sample, pattern="0h|6h|24h"))  %>% 
      column_to_rownames('Sample')
    
    mat_manual_BP_row_factors = all_subs_gff %>%
      filter(Organism == 'BP') %>%
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
    
    #filter genes based on list obtained above
    BP_CD3_mat_manual <-BP_CD3_mat_manual %>% filter(Gene %in% rownames(mat_manual_BP_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    BP_CD3_mat_manual_ordered <- BP_CD3_mat_manual[rownames(mat_manual_BP_row_factors), ]
    
    gaps <- mat_manual_BP_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many genes for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    ## draw HM including custom annotations
    # Figure 2D
    pheatmap.type(BP_CD3_mat_manual_ordered,  annRow = mat_manual_BP_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color  = colorRampPalette(colors = c("blue", "white", "red"))(17), breaks=seq(-4, 4, by = 0.5), 
                  annotation_col = mat_manual_BP_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "B. producta aCD3 top DEG" )
    
    # CB
    
    CB_CD3_select_genes_manual <- SB093_complete_df %>% 
      filter(Organism=='CB' & Treatment == 'aCD3' & Time_point =='6h') %>%
      left_join(SB093_complete_df %>% 
                  filter(Organism=='CB' & Treatment == 'aCD3' & Time_point =='24h'), by='Gene') %>% 
      filter((padj.x<0.05 & abs(log2FoldChange.x) >1 )|(padj.y<0.05 & abs(log2FoldChange.y) >1)) %>% pull(Gene) %>% unique()
    
    CB_CD3_mat_manual <- Norm_rld_counts %>% 
      filter(Organism == 'CB') %>% select(-contains('flag'), -Organism) %>%
      filter(Gene %in% CB_CD3_select_genes_manual) %>%
      mutate(Gene = gsub("(;.*)" , "", Gene)) %>%
      mutate(Gene = gsub("ID=" , "", Gene)) %>%
      column_to_rownames('Gene') 
    
    CB_CD3_mat_manual <- CB_CD3_mat_manual - rowMeans(CB_CD3_mat_manual) # scale values by subtracting mean for each row=gene
    CB_CD3_mat_manual <- as.data.frame(CB_CD3_mat_manual) %>% select(-contains('flag')) %>%
      rownames_to_column('Gene') 
    CB_CD3_mat_manual$Gene <- gsub("(;.*)" , "", CB_CD3_mat_manual$Gene)
    CB_CD3_mat_manual$Gene <- gsub("ID=" , "", CB_CD3_mat_manual$Gene)
    
    CB_genelist <- data.frame('Gene'=CB_CD3_mat_manual$Gene) # make list of genes of interest
    
    mat_manual_CB_col_factors <- data.frame(Sample=colnames(CB_CD3_mat_manual)) %>% # create factors for HM
      mutate(Time_point=str_extract(Sample, pattern="0h|6h|24h"))  %>% 
      column_to_rownames('Sample')
    
    mat_manual_CB_row_factors = all_subs_gff %>%
      filter(Organism == 'CB') %>%
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
    
    #filter genes based on list obtained above
    CB_CD3_mat_manual <-CB_CD3_mat_manual %>% filter(Gene %in% rownames(mat_manual_CB_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    CB_CD3_mat_manual_ordered <- CB_CD3_mat_manual[rownames(mat_manual_CB_row_factors), ]
    
    gaps <- mat_manual_CB_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many genes for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    ## draw HM including custom annotations
    # Figure 2D
    pheatmap.type(CB_CD3_mat_manual_ordered,  annRow = mat_manual_CB_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color  = colorRampPalette(colors = c("blue", "white", "red"))(17), breaks=seq(-4, 4, by = 0.5), 
                  annotation_col = mat_manual_CB_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "C. bolteae aCD3 top DEG" )
   
    # PD
    
    PD_CD3_select_genes_manual <- SB093_complete_df %>% 
      filter(Organism=='PD' & Treatment == 'aCD3' & Time_point =='6h') %>%
      left_join(SB093_complete_df %>% 
                  filter(Organism=='PD' & Treatment == 'aCD3' & Time_point =='24h'), by='Gene') %>% 
      filter((padj.x<0.05 & abs(log2FoldChange.x) >1 )|(padj.y<0.05 & abs(log2FoldChange.y) >1)) %>% pull(Gene) %>% unique()
    
    PD_CD3_mat_manual <- Norm_rld_counts %>% 
      filter(Organism == 'PD') %>% select(-contains('flag'), -Organism) %>%
      filter(Gene %in% PD_CD3_select_genes_manual) %>%
      mutate(Gene = gsub("(;.*)" , "", Gene)) %>%
      mutate(Gene = gsub("ID=" , "", Gene)) %>%
      column_to_rownames('Gene') 
    
    PD_CD3_mat_manual <- PD_CD3_mat_manual - rowMeans(PD_CD3_mat_manual) # scale values by subtracting mean for each row=gene
    PD_CD3_mat_manual <- as.data.frame(PD_CD3_mat_manual) %>% select(-contains('flag')) %>%
      rownames_to_column('Gene') 
    PD_CD3_mat_manual$Gene <- gsub("(;.*)" , "", PD_CD3_mat_manual$Gene)
    PD_CD3_mat_manual$Gene <- gsub("ID=" , "", PD_CD3_mat_manual$Gene)
    
    PD_genelist <- data.frame('Gene'=PD_CD3_mat_manual$Gene) # make list of genes of interest
    
    mat_manual_PD_col_factors <- data.frame(Sample=colnames(PD_CD3_mat_manual)) %>% # create factors for HM
      mutate(Time_point=str_extract(Sample, pattern="0h|6h|24h"))  %>% 
      column_to_rownames('Sample')
    
    mat_manual_PD_row_factors = all_subs_gff %>%
      filter(Organism == 'PD') %>%
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
    
    #filter genes based on list obtained above
    PD_CD3_mat_manual <-PD_CD3_mat_manual %>% filter(Gene %in% rownames(mat_manual_PD_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    PD_CD3_mat_manual_ordered <- PD_CD3_mat_manual[rownames(mat_manual_PD_row_factors), ]
    
    gaps <- mat_manual_PD_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many genes for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    ## draw HM including custom annotations
    # Figure 2D
    pheatmap.type(PD_CD3_mat_manual_ordered,  annRow = mat_manual_PD_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color  = colorRampPalette(colors = c("blue", "white", "red"))(17), breaks=seq(-4, 4, by = 0.5), 
                  annotation_col = mat_manual_PD_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "P. distasonis aCD3 top DEG" )
    
    
    #### flagellin
    
    # BS
    
    BS_flag_select_genes_manual <- SB093_complete_df %>% 
      filter(Organism=='BS' & Treatment == 'flag' & Time_point =='6h') %>%
      left_join(SB093_complete_df %>% 
                  filter(Organism=='BS' & Treatment == 'flag' & Time_point =='24h'), by='Gene') %>% 
      filter((padj.x<0.05 & abs(log2FoldChange.x) >1 )|(padj.y<0.05 & abs(log2FoldChange.y) >1)) %>% pull(Gene) %>% unique()
    
    BS_flag_mat_manual <- Norm_rld_counts %>% 
      filter(Organism == 'BS') %>% select(-contains('CD3'), -Organism) %>%
      filter(Gene %in% BS_flag_select_genes_manual) %>%
      mutate(Gene = gsub("(;.*)" , "", Gene)) %>%
      mutate(Gene = gsub("ID=" , "", Gene)) %>%
      column_to_rownames('Gene') 
    
    BS_flag_mat_manual <- BS_flag_mat_manual - rowMeans(BS_flag_mat_manual) # scale values by subtracting mean for each row=gene
    BS_flag_mat_manual <- as.data.frame(BS_flag_mat_manual) %>% select(-contains('CD3')) %>%
      rownames_to_column('Gene') 
    BS_flag_mat_manual$Gene <- gsub("(;.*)" , "", BS_flag_mat_manual$Gene)
    BS_flag_mat_manual$Gene <- gsub("ID=" , "", BS_flag_mat_manual$Gene)
    
    BS_genelist <- data.frame('Gene'=BS_flag_mat_manual$Gene) # make list of genes of interest
    
    mat_manual_BS_col_factors <- data.frame(Sample=colnames(BS_flag_mat_manual)) %>% # create factors for HM
      mutate(Time_point=str_extract(Sample, pattern="0h|6h|24h"))  %>% 
      column_to_rownames('Sample')
    
    mat_manual_BS_row_factors = all_subs_gff %>%
      filter(Organism == 'BS') %>%
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
    
    #filter genes based on list obtained above
    BS_flag_mat_manual <-BS_flag_mat_manual %>% filter(Gene %in% rownames(mat_manual_BS_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    BS_flag_mat_manual_ordered <- BS_flag_mat_manual[rownames(mat_manual_BS_row_factors), ]
    
    gaps <- mat_manual_BS_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many genes for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    ## draw HM including custom annotations
    # Figure 2D
    pheatmap.type(BS_flag_mat_manual_ordered,  annRow = mat_manual_BS_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color  = colorRampPalette(colors = c("blue", "white", "red"))(17), breaks=seq(-4, 4, by = 0.5), 
                  annotation_col = mat_manual_BS_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "B. sartorii flag top DEG" )
    
    # BP
    
    BP_flag_select_genes_manual <- SB093_complete_df %>% 
      filter(Organism=='BP' & Treatment == 'flag' & Time_point =='6h') %>%
      left_join(SB093_complete_df %>% 
                  filter(Organism=='BP' & Treatment == 'flag' & Time_point =='24h'), by='Gene') %>% 
      filter((padj.x<0.05 & abs(log2FoldChange.x) >1 )|(padj.y<0.05 & abs(log2FoldChange.y) >1)) %>% pull(Gene) %>% unique()
    
    BP_flag_mat_manual <- Norm_rld_counts %>% 
      filter(Organism == 'BP') %>% select(-contains('CD3'), -Organism) %>%
      filter(Gene %in% BP_flag_select_genes_manual) %>%
      mutate(Gene = gsub("(;.*)" , "", Gene)) %>%
      mutate(Gene = gsub("ID=" , "", Gene)) %>%
      column_to_rownames('Gene') 
    
    BP_flag_mat_manual <- BP_flag_mat_manual - rowMeans(BP_flag_mat_manual) # scale values by subtracting mean for each row=gene
    BP_flag_mat_manual <- as.data.frame(BP_flag_mat_manual) %>% select(-contains('CD3')) %>%
      rownames_to_column('Gene') 
    BP_flag_mat_manual$Gene <- gsub("(;.*)" , "", BP_flag_mat_manual$Gene)
    BP_flag_mat_manual$Gene <- gsub("ID=" , "", BP_flag_mat_manual$Gene)
    
    BP_genelist <- data.frame('Gene'=BP_flag_mat_manual$Gene) # make list of genes of interest
    
    mat_manual_BP_col_factors <- data.frame(Sample=colnames(BP_flag_mat_manual)) %>% # create factors for HM
      mutate(Time_point=str_extract(Sample, pattern="0h|6h|24h"))  %>% 
      column_to_rownames('Sample')
    
    mat_manual_BP_row_factors = all_subs_gff %>%
      filter(Organism == 'BP') %>%
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
    
    #filter genes based on list obtained above
    BP_flag_mat_manual <-BP_flag_mat_manual %>% filter(Gene %in% rownames(mat_manual_BP_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    BP_flag_mat_manual_ordered <- BP_flag_mat_manual[rownames(mat_manual_BP_row_factors), ]
    
    gaps <- mat_manual_BP_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many genes for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    ## draw HM including custom annotations
    # Figure 2D
    pheatmap.type(BP_flag_mat_manual_ordered,  annRow = mat_manual_BP_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color  = colorRampPalette(colors = c("blue", "white", "red"))(17), breaks=seq(-4, 4, by = 0.5), 
                  annotation_col = mat_manual_BP_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "B. producta flag top DEG" )
    
    # CB
    
    CB_flag_select_genes_manual <- SB093_complete_df %>% 
      filter(Organism=='CB' & Treatment == 'flag' & Time_point =='6h') %>%
      left_join(SB093_complete_df %>% 
                  filter(Organism=='CB' & Treatment == 'flag' & Time_point =='24h'), by='Gene') %>% 
      filter((padj.x<0.05 & abs(log2FoldChange.x) >1 )|(padj.y<0.05 & abs(log2FoldChange.y) >1)) %>% pull(Gene) %>% unique()
    
    CB_flag_mat_manual <- Norm_rld_counts %>% 
      filter(Organism == 'CB') %>% select(-contains('CD3'), -Organism) %>%
      filter(Gene %in% CB_flag_select_genes_manual) %>%
      mutate(Gene = gsub("(;.*)" , "", Gene)) %>%
      mutate(Gene = gsub("ID=" , "", Gene)) %>%
      column_to_rownames('Gene') 
    
    CB_flag_mat_manual <- CB_flag_mat_manual - rowMeans(CB_flag_mat_manual) # scale values by subtracting mean for each row=gene
    CB_flag_mat_manual <- as.data.frame(CB_flag_mat_manual) %>% select(-contains('CD3')) %>%
      rownames_to_column('Gene') 
    CB_flag_mat_manual$Gene <- gsub("(;.*)" , "", CB_flag_mat_manual$Gene)
    CB_flag_mat_manual$Gene <- gsub("ID=" , "", CB_flag_mat_manual$Gene)
    
    CB_genelist <- data.frame('Gene'=CB_flag_mat_manual$Gene) # make list of genes of interest
    
    mat_manual_CB_col_factors <- data.frame(Sample=colnames(CB_flag_mat_manual)) %>% # create factors for HM
      mutate(Time_point=str_extract(Sample, pattern="0h|6h|24h"))  %>% 
      column_to_rownames('Sample')
    
    mat_manual_CB_row_factors = all_subs_gff %>%
      filter(Organism == 'CB') %>%
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
    
    #filter genes based on list obtained above
    CB_flag_mat_manual <-CB_flag_mat_manual %>% filter(Gene %in% rownames(mat_manual_CB_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    CB_flag_mat_manual_ordered <- CB_flag_mat_manual[rownames(mat_manual_CB_row_factors), ]
    
    gaps <- mat_manual_CB_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many genes for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    ## draw HM including custom annotations
    # Figure 2D
    pheatmap.type(CB_flag_mat_manual_ordered,  annRow = mat_manual_CB_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color  = colorRampPalette(colors = c("blue", "white", "red"))(17), breaks=seq(-4, 4, by = 0.5), 
                  annotation_col = mat_manual_CB_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "C. bolteae flag top DEG" )
    
    # PD
    
    PD_flag_select_genes_manual <- SB093_complete_df %>% 
      filter(Organism=='PD' & Treatment == 'flag' & Time_point =='6h') %>%
      left_join(SB093_complete_df %>% 
                  filter(Organism=='PD' & Treatment == 'flag' & Time_point =='24h'), by='Gene') %>% 
      filter((padj.x<0.05 & abs(log2FoldChange.x) >1 )|(padj.y<0.05 & abs(log2FoldChange.y) >1)) %>% pull(Gene) %>% unique()
    
    PD_flag_mat_manual <- Norm_rld_counts %>% 
      filter(Organism == 'PD') %>% select(-contains('CD3'), -Organism) %>%
      filter(Gene %in% PD_flag_select_genes_manual) %>%
      mutate(Gene = gsub("(;.*)" , "", Gene)) %>%
      mutate(Gene = gsub("ID=" , "", Gene)) %>%
      column_to_rownames('Gene') 
    
    PD_flag_mat_manual <- PD_flag_mat_manual - rowMeans(PD_flag_mat_manual) # scale values by subtracting mean for each row=gene
    PD_flag_mat_manual <- as.data.frame(PD_flag_mat_manual) %>% select(-contains('CD3')) %>%
      rownames_to_column('Gene') 
    PD_flag_mat_manual$Gene <- gsub("(;.*)" , "", PD_flag_mat_manual$Gene)
    PD_flag_mat_manual$Gene <- gsub("ID=" , "", PD_flag_mat_manual$Gene)
    
    PD_genelist <- data.frame('Gene'=PD_flag_mat_manual$Gene) # make list of genes of interest
    
    mat_manual_PD_col_factors <- data.frame(Sample=colnames(PD_flag_mat_manual)) %>% # create factors for HM
      mutate(Time_point=str_extract(Sample, pattern="0h|6h|24h"))  %>% 
      column_to_rownames('Sample')
    
    mat_manual_PD_row_factors = all_subs_gff %>%
      filter(Organism == 'PD') %>%
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
    
    #filter genes based on list obtained above
    PD_flag_mat_manual <-PD_flag_mat_manual %>% filter(Gene %in% rownames(mat_manual_PD_row_factors)) %>% column_to_rownames('Gene') %>% as.matrix()
    PD_flag_mat_manual_ordered <- PD_flag_mat_manual[rownames(mat_manual_PD_row_factors), ]
    
    gaps <- mat_manual_PD_row_factors %>% group_by(Class) %>% summarise(gaps=n()) %>% arrange(Class) %>% pull(gaps) # count how many genes for each type of class there are in order to know where to insert the gaps
    gaps <- cumsum(gaps) # do the cumulative sum of gaps so that you get the actual number of row where you need to insert the gap
    
    ## draw HM including custom annotations
    # Figure 2D
    pheatmap.type(PD_flag_mat_manual_ordered,  annRow = mat_manual_PD_row_factors, 
                  cluster_cols =F, show_rownames = F, show_colnames = F, border_color = NA, color  = colorRampPalette(colors = c("blue", "white", "red"))(17), breaks=seq(-4, 4, by = 0.5), 
                  annotation_col = mat_manual_PD_col_factors, gaps_row = gaps, gaps_col = c(3,6), annotation_colors = my_colour,
                  main = "P. distasonis flag top DEG" )
    
    