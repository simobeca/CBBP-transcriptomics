### Beacttini et al, CHM
### codes for analyses and plots shown in Figures:
### Figure 1E, Figure 2C
### Figure 3A, Supplementary Figure 3A
### Figure 3B
### Figures 3C, Supplementary Figure 3B
### Figures 1D, 2B
### Supplementary Figures 1E, 2E


##### load needed libraries ####
library(RColorBrewer)
library(DESeq2)
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
library(dplyr)
library(ggsci)

#save.image("SB093_rerun.RData")
#load("./SB093_rerun.RData")

### to run this need .gff files in working directory, 
###as well as a 'SAM_BAM_ALIGNMENTS' that contains all bowtie2 alignments files (.bam) obtained with the enclosed sh script

# prepare .gff files 

BP <- read.delim("BP.gff", header=FALSE)[-1, -c(2,6)]
BS <- read.delim("Bacteroides sartorii S.C..gff", header=FALSE)[-1, -c(2,6)]
CB <- read.delim("CB.gff", header=FALSE)[-1, -c(2,6)]
PD <- read.delim("PD.gff", header=FALSE)[-1, -c(2,6)]

colnames(BP) <- c("Chr", "Type", "Start", "End", "Strand","Frame", "GeneID")
colnames(BS) <- c("Chr", "Type", "Start", "End", "Strand","Frame", "GeneID")
colnames(CB) <- c("Chr", "Type", "Start", "End", "Strand","Frame", "GeneID")
colnames(PD) <- c("Chr", "Type", "Start", "End", "Strand","Frame", "GeneID")

BS_list <- list.files("/MyFolder/SAM_BAM_ALIGNMENTS", pattern = "BS.only_aligned.bam")
BP_list <- list.files("/MyFolder/SAM_BAM_ALIGNMENTS", pattern = "BP.only_aligned.bam")
CB_list <- list.files("/MyFolder/SAM_BAM_ALIGNMENTS", pattern = "CB.only_aligned.bam")
PD_list <- list.files("/MyFolder/SAM_BAM_ALIGNMENTS", pattern = "PD.only_aligned.bam")


#### the following codes generate raw data for each CBBP members 
#### following are codes to generate cumulative figures



##### BP Count ####

# count aligning reads using Feature counts

BP_FtCounts <- featureCounts(BP_list, annot.ext = BP, isGTFAnnotationFile=FALSE, largestOverlap =T, countChimericFragments=TRUE)
BP_Counts <- as.data.frame(BP_FtCounts$counts) %>% select(-1)

colnames(BP_Counts) <- str_replace(colnames(BP_Counts), '_sb93', '')
colnames(BP_Counts) <- str_extract(colnames(BP_Counts), regex(pattern="Sample_m[0-9]{1,2}_[0-9]{1,2}h_flag|Sample_m[0-9]{1,2}_[0-9]{1,2}h_aCD3|Sample_m[0-9]{1}_0h"))
colnames(BP_Counts) <- str_replace(colnames(BP_Counts), 'Sample_', 'BP_') 

# reorder columns
BP_Counts <- BP_Counts[,c("BP_m1_0h", "BP_m2_0h", "BP_m3_0h", "BP_m4_6h_flag", "BP_m5_6h_flag", "BP_m6_6h_flag", "BP_m7_6h_aCD3", "BP_m8_6h_aCD3", "BP_m9_6h_aCD3", "BP_m10_24h_flag", "BP_m11_24h_flag" , "BP_m12_24h_flag" , "BP_m13_24h_aCD3", "BP_m14_24h_aCD3", "BP_m15_24h_aCD3")]

# remove tRNA/rRNA reads
BP_Counts$feature <- row.names(BP_Counts)
BP_Counts <- BP_Counts %>% filter(!grepl(pattern="tRNA-[A-Z, a-z]{3}-[A-Z, a-z]{3}" , BP_Counts$feature))
BP_Counts <- BP_Counts %>% filter(!grepl(pattern="LSU rRNA|SSU rRNA|5S rRNA" , BP_Counts$feature))
rownames(BP_Counts) <- BP_Counts$feature #need to do this as filter drops rownames!
BP_Counts$feature <- NULL


##### Run BP DeSeq2 ####
BP_colData <- as.data.frame(matrix(c("0h", "0h", "0h", "6h_flag","6h_flag", "6h_flag", "6h_aCD3", "6h_aCD3", "6h_aCD3" ,"24h_flag", "24h_flag", "24h_flag","24h_aCD3", "24h_aCD3", "24h_aCD3"),nrow=15,ncol=1))
row.names(BP_colData) <-c("BP_m1_0h", "BP_m2_0h", "BP_m3_0h", "BP_m4_6h_flag", "BP_m5_6h_flag", "BP_m6_6h_flag", "BP_m7_6h_aCD3", "BP_m8_6h_aCD3", "BP_m9_6h_aCD3", "BP_m10_24h_flag", "BP_m11_24h_flag" , "BP_m12_24h_flag" , "BP_m13_24h_aCD3", "BP_m14_24h_aCD3", "BP_m15_24h_aCD3")
colnames(BP_colData)<-c('BP_condition')
BP_countData <-  as.matrix(select(BP_Counts, BP_m1_0h, BP_m2_0h, BP_m3_0h, BP_m4_6h_flag, BP_m5_6h_flag, BP_m6_6h_flag, BP_m7_6h_aCD3, BP_m8_6h_aCD3, BP_m9_6h_aCD3, BP_m10_24h_flag, BP_m11_24h_flag , BP_m12_24h_flag , BP_m13_24h_aCD3, BP_m14_24h_aCD3, BP_m15_24h_aCD3), rownames.force=T )

## run Deseq excluding low count entries
BP_dds <- DESeqDataSetFromMatrix(countData=BP_countData, colData=BP_colData, design= ~ BP_condition)
BP_dds <- BP_dds[rowSums(counts(BP_dds)) > 10, ]
BP_dds <- DESeq(BP_dds, fitType = 'parametric')

BP_table_counts_normalized <- counts(BP_dds, normalized=TRUE)

##### BP PCA Figure 3A and Supplementary Figure 3A ####

BP_rld <- rlog(BP_dds, blind=FALSE)
BP_norm_counts <- as.data.frame(assay(BP_rld))

# visualize PCA using all counts
BP_PCA <- plotPCA(BP_rld, intgroup=c("BP_condition")) + theme_bw(base_size = 18) + geom_point(size=4) + theme(aspect.ratio = 1) + ggtitle('BP') + geom_text_repel(aes(label = name)) 

pdf("BP_PCA.pdf", useDingbats = F)
BP_PCA
dev.off()

## make a custom PCAwith only significantly different genes using autoplot

# select only genes with a significant difference in any comparison
BP_LRT_dds <- DESeq(BP_dds, test="LRT", reduced=~1) 
BP_res <- as.data.frame(results(BP_LRT_dds))
BP_res$gene <- rownames(BP_res)
# make a vector with names of those significant genes
BP_res_sign <- filter(BP_res, padj < .05 & abs(log2FoldChange) > 1) 
BP_res_sign <- BP_res_sign$gene 
BP_norm_counts_LRT <- BP_norm_counts[which(rownames(BP_norm_counts) %in% BP_res_sign), ] 

t_BP_norm_counts_LRT <- as.data.frame(t(BP_norm_counts_LRT)) %>% 
  mutate(Time_point=str_extract(rownames(BP_LRT_prcomp$x), pattern="0h|6h|24h"))  %>% 
  mutate(Treatment=str_extract(rownames(BP_LRT_prcomp$x), pattern="aCD3|flag")) %>% 
  replace_na(list(Treatment="None")) 

pdf("BP_PCA_LRT_autoplot.pdf", useDingbats=FALSE)
BP_autoplot <- autoplot(prcomp(t_BP_norm_counts_LRT[-rev(seq_len(ncol(t_BP_norm_counts_LRT)))[1:2]]), data = t_BP_norm_counts_LRT) +
  geom_point(size=5,aes(shape=Treatment,color=Time_point)) +
  theme_bw(base_size=18) + 
  scale_color_manual(values=c("24h"="steelblue2","6h"="salmon","0h"="gray33")) +
  ggtitle("B. producta PCA") +
  theme(text = element_text(size=16)) +
  theme(aspect.ratio = 2/3) +
  stat_chull(aes(color = Time_point, fill = Time_point), 
             alpha = 0.1, geom = "polygon")
dev.off()


##### generate df containing DES for all BP comparisons ####
BP_flag6h_vs_0h <- as.data.frame(results(BP_dds, contrast=c("BP_condition","6h_flag","0h"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="BP_flag6h_vs_0h") 
BP_aCD36h_vs_0h <- as.data.frame(results(BP_dds, contrast=c("BP_condition","6h_aCD3","0h"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="BP_aCD36h_vs_0h") 
BP_aCD324h_vs_0h <- as.data.frame(results(BP_dds, contrast=c("BP_condition","24h_aCD3", "0h"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="BP_aCD324h_vs_0h") 
BP_flag24h_vs_0h  <- as.data.frame(results(BP_dds, contrast=c("BP_condition","24h_flag","0h"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="BP_flag24h_vs_0h") 
BP_aCD324h_vs_aCD36h <- as.data.frame(results(BP_dds, contrast=c("BP_condition","24h_aCD3", "6h_aCD3"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="BP_aCD324hh_vs_aCD36h") 
BP_flag24h_vs_flag6h <- as.data.frame(results(BP_dds, contrast=c("BP_condition","24h_flag", "6h_flag"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="BP_flag24hh_vs_flag6h") 
BP_aCD36h_vs_flag6h <- as.data.frame(results(BP_dds, contrast=c("BP_condition","6h_aCD3","6h_flag"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="BP_aCD36h_vs_flag6h")
BP_aCD324h_vs_flag24h <- as.data.frame(results(BP_dds, contrast=c("BP_condition","24h_aCD3", "24h_flag"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="BP_aCD324h_vs_flag24h") 

#### make a unique df with all the BP comparisons
BP_complete_df <- rbind(BP_flag6h_vs_0h, BP_aCD36h_vs_0h, BP_aCD324h_vs_0h, BP_flag24h_vs_0h, BP_aCD36h_vs_flag6h, BP_aCD324h_vs_flag24h) %>% 
  mutate(Time_point = str_extract(Comparison, pattern ="6h|24h")) %>% 
  mutate(Treatment = case_when(
    Comparison %in% c("BP_flag6h_vs_0h", "BP_flag24h_vs_0h") ~ "flag",
    Comparison %in% c("BP_aCD36h_vs_0h", "BP_aCD324h_vs_0h") ~ "aCD3",
    Comparison %in% c("BP_aCD36h_vs_flag6h" ,"BP_aCD324h_vs_flag24h") ~ "aCD3vsflag"))  %>% 
  mutate(Gene_name = gsub("(.*);Name=", "", Gene)) %>% mutate(Gene_name = gsub("(\\(E.*)", "", Gene_name)) %>%
  mutate(Gene_number = gsub("(;.*)" , "", Gene)) %>% mutate(Gene_number = gsub("ID=fig\\|33035.10.", "", Gene_number)) %>% mutate(Gene_number = gsub("\\|", "", Gene_number))


### STATS: 1) up/down aCD3, flag, shared, unique
### graphs: 1) total heatmap (aCD3 + flag), volcano plots, meanplot, comparison (dotplot)

BP_DEG_list <- list(list_6h = list(
  aCD3_up = BP_complete_df %>% filter (padj < 0.05 & log2FoldChange > 1) %>% filter(Time_point == '6h') %>% filter(Treatment == 'aCD3') %>% pull(Gene),
  aCD3_down = BP_complete_df %>% filter (padj < 0.05 & log2FoldChange < (-1)) %>% filter(Time_point == '6h') %>% filter(Treatment == 'aCD3') %>% pull(Gene),
  flag_up = BP_complete_df %>% filter (padj < 0.05 & log2FoldChange > 1) %>% filter(Time_point == '6h') %>% filter(Treatment == 'flag') %>% pull(Gene),
  flag_down = BP_complete_df %>% filter (padj < 0.05 & log2FoldChange < (-1)) %>% filter(Time_point == '6h') %>% filter(Treatment == 'flag') %>% pull(Gene),
  aCD3vsflag_up = BP_complete_df %>% filter (padj < 0.05 & log2FoldChange > 1) %>% filter(Time_point == '6h') %>% filter(Treatment == 'aCD3vsflag') %>% pull(Gene),
  aCD3vsflag_down = BP_complete_df %>% filter (padj < 0.05 & log2FoldChange < (-1)) %>% filter(Time_point == '6h') %>% filter(Treatment == 'aCD3vsflag') %>% pull(Gene)),
  list_24h = list(
    aCD3_up = BP_complete_df %>% filter (padj < 0.05 & log2FoldChange > 1) %>% filter(Time_point == '24h') %>% filter(Treatment == 'aCD3') %>% pull(Gene),
    aCD3_down = BP_complete_df %>% filter (padj < 0.05 & log2FoldChange < (-1)) %>% filter(Time_point == '24h') %>% filter(Treatment == 'aCD3') %>% pull(Gene),
    flag_up = BP_complete_df %>% filter (padj < 0.05 & log2FoldChange > 1) %>% filter(Time_point == '24h') %>% filter(Treatment == 'flag') %>% pull(Gene),
    flag_down = BP_complete_df %>% filter (padj < 0.05 & log2FoldChange < (-1)) %>% filter(Time_point == '24h') %>% filter(Treatment == 'flag') %>% pull(Gene),
    aCD3vsflag_up = BP_complete_df %>% filter (padj < 0.05 & log2FoldChange > 1) %>% filter(Time_point == '24h') %>% filter(Treatment == 'aCD3vsflag') %>% pull(Gene),
    aCD3vsflag_down = BP_complete_df %>% filter (padj < 0.05 & log2FoldChange < (-1)) %>% filter(Time_point == '24h') %>% filter(Treatment == 'aCD3vsflag') %>% pull(Gene)))


BP_DEG_stats <- data.frame("6h" = c(
  length(BP_DEG_list$list_6h$aCD3_up),
  length(BP_DEG_list$list_6h$aCD3_down),
  length(BP_DEG_list$list_6h$flag_up),
  length(BP_DEG_list$list_6h$flag_down),
  length(BP_DEG_list$list_6h$aCD3vsflag_up),
  length(BP_DEG_list$list_6h$aCD3vsflag_down),
  length(intersect(BP_DEG_list$list_6h$aCD3_up, BP_DEG_list$list_6h$flag_up)),
  length(intersect(BP_DEG_list$list_6h$aCD3_down, BP_DEG_list$list_6h$flag_down))), 
  "24" = c(
    length(BP_DEG_list$list_24h$aCD3_up),
    length(BP_DEG_list$list_24h$aCD3_down),
    length(BP_DEG_list$list_24h$flag_up),
    length(BP_DEG_list$list_24h$flag_down),
    length(BP_DEG_list$list_24h$aCD3vsflag_up),
    length(BP_DEG_list$list_24h$aCD3vsflag_down),
    length(intersect(BP_DEG_list$list_24h$aCD3_up, BP_DEG_list$list_24h$flag_up)),
    length(intersect(BP_DEG_list$list_24h$aCD3_down, BP_DEG_list$list_24h$flag_down))),
  "Treatment" = rep(c('aCD3', 'flag', 'aCD3vsflag', 'both'),each = 2),
  "Up_Down" = rep(c("UP", "DOWN"), 4))

##### generate DF joining relevant time points or conditions

BP_all_6h <- BP_aCD36h_vs_0h %>% full_join(BP_flag6h_vs_0h, by='Gene') %>% full_join(BP_aCD36h_vs_flag6h, by='Gene') %>%
  mutate(Gene_name = gsub("(.*);Name=", "", Gene)) %>% mutate(Gene_name = gsub("(\\(E.*)", "", Gene_name)) %>%
  mutate(Gene_number = gsub("(;.*)" , "", Gene)) %>% mutate(Gene_number = gsub("ID=fig\\|33035.10.", "", Gene_number)) %>% mutate(Gene_number = gsub("\\|", "", Gene_number))

BP_all_24h <- BP_aCD324h_vs_0h %>% full_join(BP_flag24h_vs_0h, by='Gene') %>% full_join(BP_aCD324h_vs_flag24h, by='Gene') %>%
  mutate(Gene_name = gsub("(.*);Name=", "", Gene)) %>% mutate(Gene_name = gsub("(\\(E.*)", "", Gene_name)) %>%
  mutate(Gene_number = gsub("(;.*)" , "", Gene)) %>% mutate(Gene_number = gsub("ID=fig\\|33035.10.", "", Gene_number)) %>% mutate(Gene_number = gsub("\\|", "", Gene_number))

BP_all_aCD3 <- BP_aCD36h_vs_0h %>% full_join(BP_aCD324h_vs_0h, by='Gene') %>% full_join(BP_aCD324h_vs_aCD36h, by='Gene') %>%
  mutate(Gene_name = gsub("(.*);Name=", "", Gene)) %>% mutate(Gene_name = gsub("(\\(E.*)", "", Gene_name)) %>%
  mutate(Gene_number = gsub("(;.*)" , "", Gene)) %>% mutate(Gene_number = gsub("ID=fig\\|33035.10.", "", Gene_number)) %>% mutate(Gene_number = gsub("\\|", "", Gene_number))

BP_all_flag <- BP_flag6h_vs_0h %>% full_join(BP_flag24h_vs_0h, by='Gene') %>% full_join(BP_flag24h_vs_flag6h, by='Gene') %>%
  mutate(Gene_name = gsub("(.*);Name=", "", Gene)) %>% mutate(Gene_name = gsub("(\\(E.*)", "", Gene_name)) %>%
  mutate(Gene_number = gsub("(;.*)" , "", Gene)) %>% mutate(Gene_number = gsub("ID=fig\\|33035.10.", "", Gene_number)) %>% mutate(Gene_number = gsub("\\|", "", Gene_number))


BP_aCD3_vs_flag_6h_forgg <- BP_all_6h %>% mutate(Rel_Sign = case_when(
  padj < 0.05 ~ "true_dif",
  padj.x < 0.05 & padj.y < 0.05 ~ "both",
  padj.x < 0.05 & (padj.y > 0.05 | padj.y==NA) ~ "aCD3",
  (padj.x > 0.05 | padj.x==NA) & padj.y < 0.05 ~ "flag",
  (padj.x > 0.05 | padj.x==NA) & (padj.y > 0.05 | padj.y==NA) ~ "none")) %>%
  mutate(FC_Sign = case_when(
    log2FoldChange.x > 0 & log2FoldChange.y < 0 ~ "disc, aCD>0", 
    log2FoldChange.x < 0 & log2FoldChange.y > 0 ~ "disc, flag>0",
    log2FoldChange.x > 0 & log2FoldChange.y > 0  ~ "both >0",
    log2FoldChange.x < 0 & log2FoldChange.y < 0 ~ "both <0"))

only_flag_high <- BP_aCD3_vs_flag_6h_forgg %>% filter(Rel_Sign == 'flag') %>% filter(log2FoldChange.y > 1.5) %>% select(Gene, log2FoldChange.y, padj.y, Gene_name, Gene_number)
write.xlsx(only_flag_high, 'BP_up_in_glagellin_only.xlsx')



##### BS Count ####

# count aligning reads using Feature counts

BS_FtCounts <- featureCounts(BS_list, annot.ext = BS, isGTFAnnotationFile=FALSE, largestOverlap =T, countChimericFragments=TRUE)
BS_Counts <- as.data.frame(BS_FtCounts$counts) %>% select(-1)

colnames(BS_Counts) <- str_replace(colnames(BS_Counts), '_sb93', '')
colnames(BS_Counts) <- str_extract(colnames(BS_Counts), regex(pattern="Sample_m[0-9]{1,2}_[0-9]{1,2}h_flag|Sample_m[0-9]{1,2}_[0-9]{1,2}h_aCD3|Sample_m[0-9]{1}_0h"))
colnames(BS_Counts) <- str_replace(colnames(BS_Counts), 'Sample_', 'BS_') 

# reorder columns
BS_Counts <- BS_Counts[,c("BS_m1_0h", "BS_m2_0h", "BS_m3_0h", "BS_m4_6h_flag", "BS_m5_6h_flag", "BS_m6_6h_flag", "BS_m7_6h_aCD3", "BS_m8_6h_aCD3", "BS_m9_6h_aCD3", "BS_m10_24h_flag", "BS_m11_24h_flag" , "BS_m12_24h_flag" , "BS_m13_24h_aCD3", "BS_m14_24h_aCD3", "BS_m15_24h_aCD3")]

# remove tRNA/rRNA reads
BS_Counts$feature <- row.names(BS_Counts)
BS_Counts <- BS_Counts %>% filter(!grepl(pattern="tRNA-[A-Z, a-z]{3}-[A-Z, a-z]{3}" , BS_Counts$feature))
BS_Counts <- BS_Counts %>% filter(!grepl(pattern="LSU rRNA|SSU rRNA|5S rRNA" , BS_Counts$feature))
rownames(BS_Counts) <- BS_Counts$feature #need to do this as filter drops rownames!
BS_Counts$feature <- NULL


##### Run BS DeSeq2 ####
BS_colData <- as.data.frame(matrix(c("0h", "0h", "0h", "6h_flag","6h_flag", "6h_flag", "6h_aCD3", "6h_aCD3", "6h_aCD3" ,"24h_flag", "24h_flag", "24h_flag","24h_aCD3", "24h_aCD3", "24h_aCD3"),nrow=15,ncol=1))
row.names(BS_colData) <-c("BS_m1_0h", "BS_m2_0h", "BS_m3_0h", "BS_m4_6h_flag", "BS_m5_6h_flag", "BS_m6_6h_flag", "BS_m7_6h_aCD3", "BS_m8_6h_aCD3", "BS_m9_6h_aCD3", "BS_m10_24h_flag", "BS_m11_24h_flag" , "BS_m12_24h_flag" , "BS_m13_24h_aCD3", "BS_m14_24h_aCD3", "BS_m15_24h_aCD3")
colnames(BS_colData)<-c('BS_condition')
BS_countData <-  as.matrix(select(BS_Counts, BS_m1_0h, BS_m2_0h, BS_m3_0h, BS_m4_6h_flag, BS_m5_6h_flag, BS_m6_6h_flag, BS_m7_6h_aCD3, BS_m8_6h_aCD3, BS_m9_6h_aCD3, BS_m10_24h_flag, BS_m11_24h_flag , BS_m12_24h_flag , BS_m13_24h_aCD3, BS_m14_24h_aCD3, BS_m15_24h_aCD3), rownames.force=T )

## run Deseq excluding low count entries
BS_dds <- DESeqDataSetFromMatrix(countData=BS_countData, colData=BS_colData, design= ~ BS_condition)
BS_dds <- BS_dds[rowSums(counts(BS_dds)) > 10, ]
BS_dds <- DESeq(BS_dds, fitType = 'parametric')

BS_table_counts_normalized <- counts(BS_dds, normalized=TRUE)

##### BS PCA Figure 3A and Supplementary Figure 3A ####

BS_rld <- rlog(BS_dds, blind=FALSE)
BS_norm_counts <- as.data.frame(assay(BS_rld))

# visualize PCA using all counts
BS_PCA <- plotPCA(BS_rld, intgroup=c("BS_condition")) + theme_bw(base_size = 18) + geom_point(size=4) + theme(aspect.ratio = 1) + ggtitle('BS') + geom_text_repel(aes(label = name)) 

pdf("BS_PCA.pdf", useDingbats = F)
BS_PCA
dev.off()

##make a custom PCAwith only significantly different genes using autoplot
# select only genes with a significant difference in any comparison
BS_LRT_dds <- DESeq(BS_dds, test="LRT", reduced=~1) 
BS_res <- as.data.frame(results(BS_LRT_dds))
BS_res$gene <- rownames(BS_res)
# make a vector with names of those significant genes
BS_res_sign <- filter(BS_res, padj < .05 & abs(log2FoldChange) > 1) 
BS_res_sign <- BS_res_sign$gene 
BS_norm_counts_LRT <- BS_norm_counts[which(rownames(BS_norm_counts) %in% BS_res_sign), ] 

t_BS_norm_counts_LRT <- as.data.frame(t(BS_norm_counts_LRT)) %>% 
  mutate(Time_point=str_extract(rownames(BS_LRT_prcomp$x), pattern="0h|6h|24h"))  %>% 
  mutate(Treatment=str_extract(rownames(BS_LRT_prcomp$x), pattern="aCD3|flag")) %>% 
  replace_na(list(Treatment="None")) 

pdf("BS_PCA_LRT_autoplot.pdf", useDingbats=FALSE)
BS_autoplot <- autoplot(prcomp(t_BS_norm_counts_LRT[-rev(seq_len(ncol(t_BS_norm_counts_LRT)))[1:2]]), data = t_BS_norm_counts_LRT) +
  geom_point(size=5,aes(shape=Treatment,color=Time_point)) +
  theme_bw(base_size=18) + 
  scale_color_manual(values=c("24h"="steelblue2","6h"="salmon","0h"="gray33")) +
  ggtitle("B. producta PCA") +
  theme(text = element_text(size=16)) +
  theme(aspect.ratio = 2/3) +
  stat_chull(aes(color = Time_point, fill = Time_point), 
             alpha = 0.1, geom = "polygon")
dev.off()


##### generate df containing DES for all BS comparisons ####
BS_flag6h_vs_0h <- as.data.frame(results(BS_dds, contrast=c("BS_condition","6h_flag","0h"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="BS_flag6h_vs_0h") 
BS_aCD36h_vs_0h <- as.data.frame(results(BS_dds, contrast=c("BS_condition","6h_aCD3","0h"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="BS_aCD36h_vs_0h") 
BS_aCD324h_vs_0h <- as.data.frame(results(BS_dds, contrast=c("BS_condition","24h_aCD3", "0h"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="BS_aCD324h_vs_0h") 
BS_flag24h_vs_0h  <- as.data.frame(results(BS_dds, contrast=c("BS_condition","24h_flag","0h"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="BS_flag24h_vs_0h") 
BS_aCD324h_vs_aCD36h <- as.data.frame(results(BS_dds, contrast=c("BS_condition","24h_aCD3", "6h_aCD3"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="BS_aCD324hh_vs_aCD36h") 
BS_flag24h_vs_flag6h <- as.data.frame(results(BS_dds, contrast=c("BS_condition","24h_flag", "6h_flag"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="BS_flag24hh_vs_flag6h") 
BS_aCD36h_vs_flag6h <- as.data.frame(results(BS_dds, contrast=c("BS_condition","6h_aCD3","6h_flag"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="BS_aCD36h_vs_flag6h")
BS_aCD324h_vs_flag24h <- as.data.frame(results(BS_dds, contrast=c("BS_condition","24h_aCD3", "24h_flag"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="BS_aCD324h_vs_flag24h") 

#### make a unique df with all the BS comparisons
BS_complete_df <- rbind(BS_flag6h_vs_0h, BS_aCD36h_vs_0h, BS_aCD324h_vs_0h, BS_flag24h_vs_0h, BS_aCD36h_vs_flag6h, BS_aCD324h_vs_flag24h) %>% 
  mutate(Time_point = str_extract(Comparison, pattern ="6h|24h")) %>% 
  mutate(Treatment = case_when(
    Comparison %in% c("BS_flag6h_vs_0h", "BS_flag24h_vs_0h") ~ "flag",
    Comparison %in% c("BS_aCD36h_vs_0h", "BS_aCD324h_vs_0h") ~ "aCD3",
    Comparison %in% c("BS_aCD36h_vs_flag6h" ,"BS_aCD324h_vs_flag24h") ~ "aCD3vsflag"))  %>% 
  mutate(Gene_name = gsub("(.*);Name=", "", Gene)) %>% mutate(Gene_name = gsub("(\\(E.*)", "", Gene_name)) %>%
  mutate(Gene_number = gsub("(;.*)" , "", Gene)) %>% mutate(Gene_number = gsub("ID=fig\\|671267.8.", "", Gene_number)) %>% mutate(Gene_number = gsub("\\|", "", Gene_number))




##### CB Count ####

# count aligning reads using Feature counts

CB_FtCounts <- featureCounts(CB_list, annot.ext = CB, isGTFAnnotationFile=FALSE, largestOverlap =T, countChimericFragments=TRUE)
CB_Counts <- as.data.frame(CB_FtCounts$counts) %>% select(-1)

colnames(CB_Counts) <- str_replace(colnames(CB_Counts), '_sb93', '')
colnames(CB_Counts) <- str_extract(colnames(CB_Counts), regex(pattern="Sample_m[0-9]{1,2}_[0-9]{1,2}h_flag|Sample_m[0-9]{1,2}_[0-9]{1,2}h_aCD3|Sample_m[0-9]{1}_0h"))
colnames(CB_Counts) <- str_replace(colnames(CB_Counts), 'Sample_', 'CB_') 

# reorder columns
CB_Counts <- CB_Counts[,c("CB_m1_0h", "CB_m2_0h", "CB_m3_0h", "CB_m4_6h_flag", "CB_m5_6h_flag", "CB_m6_6h_flag", "CB_m7_6h_aCD3", "CB_m8_6h_aCD3", "CB_m9_6h_aCD3", "CB_m10_24h_flag", "CB_m11_24h_flag" , "CB_m12_24h_flag" , "CB_m13_24h_aCD3", "CB_m14_24h_aCD3", "CB_m15_24h_aCD3")]

# remove tRNA/rRNA reads
CB_Counts$feature <- row.names(CB_Counts)
CB_Counts <- CB_Counts %>% filter(!grepl(pattern="tRNA-[A-Z, a-z]{3}-[A-Z, a-z]{3}" , CB_Counts$feature))
CB_Counts <- CB_Counts %>% filter(!grepl(pattern="LSU rRNA|SSU rRNA|5S rRNA" , CB_Counts$feature))
rownames(CB_Counts) <- CB_Counts$feature #need to do this as filter drops rownames!
CB_Counts$feature <- NULL


##### Run CB DeSeq2 ####
CB_colData <- as.data.frame(matrix(c("0h", "0h", "0h", "6h_flag","6h_flag", "6h_flag", "6h_aCD3", "6h_aCD3", "6h_aCD3" ,"24h_flag", "24h_flag", "24h_flag","24h_aCD3", "24h_aCD3", "24h_aCD3"),nrow=15,ncol=1))
row.names(CB_colData) <-c("CB_m1_0h", "CB_m2_0h", "CB_m3_0h", "CB_m4_6h_flag", "CB_m5_6h_flag", "CB_m6_6h_flag", "CB_m7_6h_aCD3", "CB_m8_6h_aCD3", "CB_m9_6h_aCD3", "CB_m10_24h_flag", "CB_m11_24h_flag" , "CB_m12_24h_flag" , "CB_m13_24h_aCD3", "CB_m14_24h_aCD3", "CB_m15_24h_aCD3")
colnames(CB_colData)<-c('CB_condition')
CB_countData <-  as.matrix(select(CB_Counts, CB_m1_0h, CB_m2_0h, CB_m3_0h, CB_m4_6h_flag, CB_m5_6h_flag, CB_m6_6h_flag, CB_m7_6h_aCD3, CB_m8_6h_aCD3, CB_m9_6h_aCD3, CB_m10_24h_flag, CB_m11_24h_flag , CB_m12_24h_flag , CB_m13_24h_aCD3, CB_m14_24h_aCD3, CB_m15_24h_aCD3), rownames.force=T )

## run Deseq excluding low count entries
CB_dds <- DESeqDataSetFromMatrix(countData=CB_countData, colData=CB_colData, design= ~ CB_condition)
CB_dds <- CB_dds[rowSums(counts(CB_dds)) > 10, ]
CB_dds <- DESeq(CB_dds, fitType = 'parametric')

CB_table_counts_normalized <- counts(CB_dds, normalized=TRUE)

##### CB PCA Figure 3A and Supplementary Figure 3A ####

CB_rld <- rlog(CB_dds, blind=FALSE)
CB_norm_counts <- as.data.frame(assay(CB_rld))

# visualize PCA using all counts
CB_PCA <- plotPCA(CB_rld, intgroup=c("CB_condition")) + theme_bw(base_size = 18) + geom_point(size=4) + theme(aspect.ratio = 1) + ggtitle('CB') + geom_text_repel(aes(label = name)) 

pdf("CB_PCA.pdf", useDingbats = F)
CB_PCA
dev.off()

##make a custom PCAwith only significantly different genes using autoplot
# select only genes with a significant difference in any comparison
CB_LRT_dds <- DESeq(CB_dds, test="LRT", reduced=~1) 
CB_res <- as.data.frame(results(CB_LRT_dds))
CB_res$gene <- rownames(CB_res)
# make a vector with names of those significant genes
CB_res_sign <- filter(CB_res, padj < .05 & abs(log2FoldChange) > 1) 
CB_res_sign <- CB_res_sign$gene 
CB_norm_counts_LRT <- CB_norm_counts[which(rownames(CB_norm_counts) %in% CB_res_sign), ] 

t_CB_norm_counts_LRT <- as.data.frame(t(CB_norm_counts_LRT)) %>% 
  mutate(Time_point=str_extract(rownames(CB_LRT_prcomp$x), pattern="0h|6h|24h"))  %>% 
  mutate(Treatment=str_extract(rownames(CB_LRT_prcomp$x), pattern="aCD3|flag")) %>% 
  replace_na(list(Treatment="None")) 

pdf("CB_PCA_LRT_autoplot.pdf", useDingbats=FALSE)
CB_autoplot <- autoplot(prcomp(t_CB_norm_counts_LRT[-rev(seq_len(ncol(t_CB_norm_counts_LRT)))[1:2]]), data = t_CB_norm_counts_LRT) +
  geom_point(size=5,aes(shape=Treatment,color=Time_point)) +
  theme_bw(base_size=18) + 
  scale_color_manual(values=c("24h"="steelblue2","6h"="salmon","0h"="gray33")) +
  ggtitle("B. producta PCA") +
  theme(text = element_text(size=16)) +
  theme(aspect.ratio = 2/3) +
  stat_chull(aes(color = Time_point, fill = Time_point), 
             alpha = 0.1, geom = "polygon")
dev.off()


##### generate df containing DES for all CB comparisons ####
CB_flag6h_vs_0h <- as.data.frame(results(CB_dds, contrast=c("CB_condition","6h_flag","0h"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="CB_flag6h_vs_0h") 
CB_aCD36h_vs_0h <- as.data.frame(results(CB_dds, contrast=c("CB_condition","6h_aCD3","0h"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="CB_aCD36h_vs_0h") 
CB_aCD324h_vs_0h <- as.data.frame(results(CB_dds, contrast=c("CB_condition","24h_aCD3", "0h"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="CB_aCD324h_vs_0h") 
CB_flag24h_vs_0h  <- as.data.frame(results(CB_dds, contrast=c("CB_condition","24h_flag","0h"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="CB_flag24h_vs_0h") 
CB_aCD324h_vs_aCD36h <- as.data.frame(results(CB_dds, contrast=c("CB_condition","24h_aCD3", "6h_aCD3"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="CB_aCD324hh_vs_aCD36h") 
CB_flag24h_vs_flag6h <- as.data.frame(results(CB_dds, contrast=c("CB_condition","24h_flag", "6h_flag"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="CB_flag24hh_vs_flag6h") 
CB_aCD36h_vs_flag6h <- as.data.frame(results(CB_dds, contrast=c("CB_condition","6h_aCD3","6h_flag"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="CB_aCD36h_vs_flag6h")
CB_aCD324h_vs_flag24h <- as.data.frame(results(CB_dds, contrast=c("CB_condition","24h_aCD3", "24h_flag"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="CB_aCD324h_vs_flag24h") 

#### make a unique df with all the CB comparisons
CB_complete_df <- rbind(CB_flag6h_vs_0h, CB_aCD36h_vs_0h, CB_aCD324h_vs_0h, CB_flag24h_vs_0h, CB_aCD36h_vs_flag6h, CB_aCD324h_vs_flag24h) %>% 
  mutate(Time_point = str_extract(Comparison, pattern ="6h|24h")) %>% 
  mutate(Treatment = case_when(
    Comparison %in% c("CB_flag6h_vs_0h", "CB_flag24h_vs_0h") ~ "flag",
    Comparison %in% c("CB_aCD36h_vs_0h", "CB_aCD324h_vs_0h") ~ "aCD3",
    Comparison %in% c("CB_aCD36h_vs_flag6h" ,"CB_aCD324h_vs_flag24h") ~ "aCD3vsflag"))  %>% 
  mutate(Gene_name = gsub("(.*);Name=", "", Gene)) %>% mutate(Gene_name = gsub("(\\(E.*)", "", Gene_name)) %>%
  mutate(Gene_number = gsub("(;.*)" , "", Gene)) %>% mutate(Gene_number = gsub("ID=fig\\|208479.6.", "", Gene_number)) %>% mutate(Gene_number = gsub("\\|", "", Gene_number))





##### PD Count ####

# count aligning reads using Feature counts

PD_FtCounts <- featureCounts(PD_list, annot.ext = PD, isGTFAnnotationFile=FALSE, largestOverlap =T, countChimericFragments=TRUE)
PD_Counts <- as.data.frame(PD_FtCounts$counts) %>% select(-1)

colnames(PD_Counts) <- str_replace(colnames(PD_Counts), '_sb93', '')
colnames(PD_Counts) <- str_extract(colnames(PD_Counts), regex(pattern="Sample_m[0-9]{1,2}_[0-9]{1,2}h_flag|Sample_m[0-9]{1,2}_[0-9]{1,2}h_aCD3|Sample_m[0-9]{1}_0h"))
colnames(PD_Counts) <- str_replace(colnames(PD_Counts), 'Sample_', 'PD_') 

# reorder columns
PD_Counts <- PD_Counts[,c("PD_m1_0h", "PD_m2_0h", "PD_m3_0h", "PD_m4_6h_flag", "PD_m5_6h_flag", "PD_m6_6h_flag", "PD_m7_6h_aCD3", "PD_m8_6h_aCD3", "PD_m9_6h_aCD3", "PD_m10_24h_flag", "PD_m11_24h_flag" , "PD_m12_24h_flag" , "PD_m13_24h_aCD3", "PD_m14_24h_aCD3", "PD_m15_24h_aCD3")]

# remove tRNA/rRNA reads
PD_Counts$feature <- row.names(PD_Counts)
PD_Counts <- PD_Counts %>% filter(!grepl(pattern="tRNA-[A-Z, a-z]{3}-[A-Z, a-z]{3}" , PD_Counts$feature))
PD_Counts <- PD_Counts %>% filter(!grepl(pattern="LSU rRNA|SSU rRNA|5S rRNA" , PD_Counts$feature))
rownames(PD_Counts) <- PD_Counts$feature #need to do this as filter drops rownames!
PD_Counts$feature <- NULL


##### Run PD DeSeq2 ####
PD_colData <- as.data.frame(matrix(c("0h", "0h", "0h", "6h_flag","6h_flag", "6h_flag", "6h_aCD3", "6h_aCD3", "6h_aCD3" ,"24h_flag", "24h_flag", "24h_flag","24h_aCD3", "24h_aCD3", "24h_aCD3"),nrow=15,ncol=1))
row.names(PD_colData) <-c("PD_m1_0h", "PD_m2_0h", "PD_m3_0h", "PD_m4_6h_flag", "PD_m5_6h_flag", "PD_m6_6h_flag", "PD_m7_6h_aCD3", "PD_m8_6h_aCD3", "PD_m9_6h_aCD3", "PD_m10_24h_flag", "PD_m11_24h_flag" , "PD_m12_24h_flag" , "PD_m13_24h_aCD3", "PD_m14_24h_aCD3", "PD_m15_24h_aCD3")
colnames(PD_colData)<-c('PD_condition')
PD_countData <-  as.matrix(select(PD_Counts, PD_m1_0h, PD_m2_0h, PD_m3_0h, PD_m4_6h_flag, PD_m5_6h_flag, PD_m6_6h_flag, PD_m7_6h_aCD3, PD_m8_6h_aCD3, PD_m9_6h_aCD3, PD_m10_24h_flag, PD_m11_24h_flag , PD_m12_24h_flag , PD_m13_24h_aCD3, PD_m14_24h_aCD3, PD_m15_24h_aCD3), rownames.force=T )

## run Deseq excluding low count entries
PD_dds <- DESeqDataSetFromMatrix(countData=PD_countData, colData=PD_colData, design= ~ PD_condition)
PD_dds <- PD_dds[rowSums(counts(PD_dds)) > 10, ]
PD_dds <- DESeq(PD_dds, fitType = 'parametric')

PD_table_counts_normalized <- counts(PD_dds, normalized=TRUE)

##### PD PCA Figure 3A and Supplementary Figure 3A ####

PD_rld <- rlog(PD_dds, blind=FALSE)
PD_norm_counts <- as.data.frame(assay(PD_rld))

# visualize PCA using all counts
PD_PCA <- plotPCA(PD_rld, intgroup=c("PD_condition")) + theme_bw(base_size = 18) + geom_point(size=4) + theme(aspect.ratio = 1) + ggtitle('PD') + geom_text_repel(aes(label = name)) 

pdf("PD_PCA.pdf", useDingbats = F)
PD_PCA
dev.off()

##make a custom PCAwith only significantly different genes using autoplot
# select only genes with a significant difference in any comparison
PD_LRT_dds <- DESeq(PD_dds, test="LRT", reduced=~1) 
PD_res <- as.data.frame(results(PD_LRT_dds))
PD_res$gene <- rownames(PD_res)
# make a vector with names of those significant genes
PD_res_sign <- filter(PD_res, padj < .05 & abs(log2FoldChange) > 1) 
PD_res_sign <- PD_res_sign$gene 
PD_norm_counts_LRT <- PD_norm_counts[which(rownames(PD_norm_counts) %in% PD_res_sign), ] 

t_PD_norm_counts_LRT <- as.data.frame(t(PD_norm_counts_LRT)) %>% 
  mutate(Time_point=str_extract(rownames(PD_LRT_prcomp$x), pattern="0h|6h|24h"))  %>% 
  mutate(Treatment=str_extract(rownames(PD_LRT_prcomp$x), pattern="aCD3|flag")) %>% 
  replace_na(list(Treatment="None")) 

pdf("PD_PCA_LRT_autoplot.pdf", useDingbats=FALSE)
PD_autoplot <- autoplot(prcomp(t_PD_norm_counts_LRT[-rev(seq_len(ncol(t_PD_norm_counts_LRT)))[1:2]]), data = t_PD_norm_counts_LRT) +
  geom_point(size=5,aes(shape=Treatment,color=Time_point)) +
  theme_bw(base_size=18) + 
  scale_color_manual(values=c("24h"="steelblue2","6h"="salmon","0h"="gray33")) +
  ggtitle("B. producta PCA") +
  theme(text = element_text(size=16)) +
  theme(aspect.ratio = 2/3) +
  stat_chull(aes(color = Time_point, fill = Time_point), 
             alpha = 0.1, geom = "polygon")
dev.off()


##### generate df containing DES for all PD comparisons ####
PD_flag6h_vs_0h <- as.data.frame(results(PD_dds, contrast=c("PD_condition","6h_flag","0h"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="PD_flag6h_vs_0h") 
PD_aCD36h_vs_0h <- as.data.frame(results(PD_dds, contrast=c("PD_condition","6h_aCD3","0h"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="PD_aCD36h_vs_0h") 
PD_aCD324h_vs_0h <- as.data.frame(results(PD_dds, contrast=c("PD_condition","24h_aCD3", "0h"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="PD_aCD324h_vs_0h") 
PD_flag24h_vs_0h  <- as.data.frame(results(PD_dds, contrast=c("PD_condition","24h_flag","0h"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="PD_flag24h_vs_0h") 
PD_aCD324h_vs_aCD36h <- as.data.frame(results(PD_dds, contrast=c("PD_condition","24h_aCD3", "6h_aCD3"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="PD_aCD324hh_vs_aCD36h") 
PD_flag24h_vs_flag6h <- as.data.frame(results(PD_dds, contrast=c("PD_condition","24h_flag", "6h_flag"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="PD_flag24hh_vs_flag6h") 
PD_aCD36h_vs_flag6h <- as.data.frame(results(PD_dds, contrast=c("PD_condition","6h_aCD3","6h_flag"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="PD_aCD36h_vs_flag6h")
PD_aCD324h_vs_flag24h <- as.data.frame(results(PD_dds, contrast=c("PD_condition","24h_aCD3", "24h_flag"))) %>% rownames_to_column('Gene') %>% mutate(Comparison="PD_aCD324h_vs_flag24h") 

#### make a unique df with all the PD comparisons
PD_complete_df <- rbind(PD_flag6h_vs_0h, PD_aCD36h_vs_0h, PD_aCD324h_vs_0h, PD_flag24h_vs_0h, PD_aCD36h_vs_flag6h, PD_aCD324h_vs_flag24h) %>% 
  mutate(Time_point = str_extract(Comparison, pattern ="6h|24h")) %>% 
  mutate(Treatment = case_when(
    Comparison %in% c("PD_flag6h_vs_0h", "PD_flag24h_vs_0h") ~ "flag",
    Comparison %in% c("PD_aCD36h_vs_0h", "PD_aCD324h_vs_0h") ~ "aCD3",
    Comparison %in% c("PD_aCD36h_vs_flag6h" ,"PD_aCD324h_vs_flag24h") ~ "aCD3vsflag"))  %>% 
  mutate(Gene_name = gsub("(.*);Name=", "", Gene)) %>% mutate(Gene_name = gsub("(\\(E.*)", "", Gene_name)) %>%
  mutate(Gene_number = gsub("(;.*)" , "", Gene)) %>% mutate(Gene_number = gsub("ID=fig\\|823.92.peg.", "", Gene_number)) %>% mutate(Gene_number = gsub("\\|", "", Gene_number))



###########################################
############ CUMULATIVE PLOTS #############
###########################################


##### 16s graphs Figures 1D 2B ####

df <- read.delim('DADA2_SB093_SB172_table.txt')

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

SB093_flag_cecum_gg <- SB093_flag %>%
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

pdf('SB093_flag_16s_cecum_posttreatment.pdf')
SB093_flag_cecum_gg
dev.off()


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

SB093_aCD3_cecum_gg <- SB093_aCD3 %>%
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

pdf('SB093_aCD3_16s_cecum_posttreatment.pdf')
SB093_aCD3_cecum_gg
dev.off()

##### RNA/DNA ratios Supplementary Figures 1E 2E #####


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



##### volcano plots flagellin Figure 1E ####

BP_f6 <- BP_flag6h_vs_0h %>% mutate(Organism='BP')
BS_f6 <- BS_flag6h_vs_0h %>% mutate(Organism='BS')
CB_f6 <- CB_flag6h_vs_0h %>% mutate(Organism='CB')
PD_f6 <- PD_flag6h_vs_0h %>% mutate(Organism='PD')

Flag_6h_all <- bind_rows(BP_f6, BS_f6, CB_f6, PD_f6)

flag_top_5 <- Flag_6h_all %>%
  group_by(Organism) %>%
  arrange(padj) %>%
  top_n(-5, padj) %>% 
  top_n(5, abs(log2FoldChange)) %>%
  ungroup() %>%
  pull(Gene)
  
Flag_6h_all$padj<0.05  
Flag_6h_all$log2FoldChange>1

flag_counts <- Flag_6h_all %>%
  group_by(Organism) %>%
  summarize(n.top=sum(padj<0.05 & log2FoldChange>1,na.rm=TRUE),
            n.bottom=sum(padj<0.05 & log2FoldChange< (-1),na.rm=TRUE)) %>%
  ungroup()


volcanoflag <- Flag_6h_all %>%
  mutate(Significant = case_when(
    padj > 0.05 | (abs(log2FoldChange) < 1) ~ 'Not Sign',
    padj < 0.05 & log2FoldChange > 1 ~ 'UP',
    padj < 0.05 & log2FoldChange <(-1) ~ 'DOWN'
  )) %>% 
  arrange(Significant!="Not Sign") %>%
  ggplot() + 
  geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color=Significant), size=3) + theme_bw(base_size=16) +
  scale_color_manual(values=c("UP"="firebrick2","Not Sign"="darkgrey", 'DOWN'='dodgerblue3')) +   
  geom_vline(xintercept = -1, colour="#990000", linetype="dashed") +
  geom_vline(xintercept = 1, colour="#990000", linetype="dashed") +
  facet_wrap(~factor(Organism, levels=c('CB', 'BP', 'BS', 'PD')), nrow=2) +
  geom_text(data=flag_counts,aes(x=4, y=22,  label=n.top, size=6)) + 
  geom_text(data=flag_counts,aes(x=-3, y=22,  label=n.bottom,  size=6)) +
  theme(strip.background = element_blank(), strip.text.y = element_blank()) 

pdf('flagellin_volcano_6h.pdf', height=6, width=6, useDingbats = F)
volcanoflag
dev.off()

png('flagellin_volcano_6h.png', res = 1200, height=4000, width=10000)
volcanoflag
dev.off()

##### volcano plots aCD3 Figure 2C ####

BP_aCD36 <- BP_aCD36h_vs_0h %>% mutate(Organism='BP')
BS_aCD36 <- BS_aCD36h_vs_0h %>% mutate(Organism='BS')
CB_aCD36 <- CB_aCD36h_vs_0h %>% mutate(Organism='CB')
PD_aCD36 <- PD_aCD36h_vs_0h %>% mutate(Organism='PD')

aCD3_6h_all <- bind_rows(BP_aCD36, BS_aCD36, CB_aCD36, PD_aCD36)

aCD3_top_5 <- aCD3_6h_all %>%
  group_by(Organism) %>%
  arrange(padj) %>%
  top_n(-5, padj) %>% 
  top_n(5, abs(log2FoldChange)) %>%
  ungroup() %>%
  pull(Gene)

aCD3_counts <- aCD3_6h_all %>%
  group_by(Organism) %>%
  summarize(n.top=sum(padj<0.05 & log2FoldChange>1,na.rm=TRUE),
            n.bottom=sum(padj<0.05 & log2FoldChange< (-1),na.rm=TRUE)) %>%
  ungroup()


volcanoaCD3 <- aCD3_6h_all %>%
  mutate(Significant = case_when(
    padj > 0.05 | (abs(log2FoldChange) < 1) ~ 'Not Sign',
    padj < 0.05 & log2FoldChange > 1 ~ 'UP',
    padj < 0.05 & log2FoldChange <(-1) ~ 'DOWN'
  )) %>% 
  arrange(Significant!="Not Sign") %>%
  ggplot() + 
  geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color=Significant), size=3) + theme_bw(base_size=16) +
  scale_color_manual(values=c("UP"="firebrick2","Not Sign"="darkgrey", 'DOWN'='dodgerblue3')) +   
  geom_vline(xintercept = -1, colour="#990000", linetype="dashed") +
  geom_vline(xintercept = 1, colour="#990000", linetype="dashed") +
  facet_wrap(~factor(Organism, levels=c('CB', 'BP', 'BS', 'PD')), nrow=2) +
  geom_text(data=aCD3_counts,aes(x=4, y=22,  label=n.top, size=6)) + 
  geom_text(data=aCD3_counts,aes(x=-3, y=22,  label=n.bottom,  size=6)) +
  theme(strip.background = element_blank(), strip.text.y = element_blank()) 

pdf('aCD3_volcano_6h.pdf', height=6, width=6, useDingbats = F)
volcanoaCD3
dev.off()

png('CD3_volcano_6h.png', res = 1200, height=4000, width=10000)
volcanoaCD3
dev.off()



##### comparison dot-plot aCD3 flagellin Figure 3B #####

all_flag_gg_df <- bind_rows(BP_aCD3_vs_flag_6h_forgg, BS_aCD3_vs_flag_6h_forgg, CB_aCD3_vs_flag_6h_forgg, PD_aCD3_vs_flag_6h_forgg) %>%
  filter(Rel_Sign=='flag') %>% filter(abs(log2FoldChange.y)>1) 

all_aCD3_vs_flag_6h_forgg <- bind_rows(BP_aCD3_vs_flag_6h_forgg, BS_aCD3_vs_flag_6h_forgg, CB_aCD3_vs_flag_6h_forgg, PD_aCD3_vs_flag_6h_forgg) %>%
  mutate(Rel_Sign = case_when(
    # padj < 0.05 ~ "true_dif",
    padj.x < 0.05 & padj.y < 0.05 ~ "both",
    padj.x < 0.05 & (padj.y > 0.05 | padj.y==NA) ~ "aCD3",
    (padj.x > 0.05 | padj.x==NA) & padj.y < 0.05 ~ "flag",
    (padj.x > 0.05 | padj.x==NA) & (padj.y > 0.05 | padj.y==NA) ~ "none")) %>%
  mutate(FC_Sign = case_when(
    log2FoldChange.x > 0 & log2FoldChange.y < 0 ~ "disc, aCD>0", 
    log2FoldChange.x < 0 & log2FoldChange.y > 0 ~ "disc, flag>0",
    log2FoldChange.x > 0 & log2FoldChange.y > 0  ~ "both >0",
    log2FoldChange.x < 0 & log2FoldChange.y < 0 ~ "both <0")) %>%
  mutate(Organism = str_extract(pattern='BS|BP|CB|PD', Comparison)) %>%
  arrange(Rel_Sign!="none",Rel_Sign!="aCD3",Rel_Sign!="flag") %>%
  mutate(alpha=ifelse(Rel_Sign=="none",0.6,1))


pick <- function(condition){
  function(d) d %>% filter_(condition)
}

all_aCD3_vs_flag_6h_gg <- 
  #filter(Rel_Sign != "none") %>%
  ggplot(all_aCD3_vs_flag_6h_forgg, aes(x=log2FoldChange.x, y=log2FoldChange.y, colour=Rel_Sign,alpha=alpha)) +
  stat_cor(inherit.aes = FALSE, mapping = aes(x=log2FoldChange.x, y=log2FoldChange.y), method="pearson", size=6) + #this is necessary as otherwise with inherited aes the value is calculated separately for the color groups!
  # geom_text_repel(data=filter(all_aCD3_vs_flag_24h_forgg, Gene %in% list_glut), segment.size  = 0.4, 
  #                 aes(label=Gene_name)) +
  #geom_hline(yintercept = 0, colour="#990000", linetype="dashed") +
  #geom_vline(xintercept = 0, coloFigure 3Fur="#990000", linetype="dashed") +
  geom_abline(slope=1, colour="#990000", linetype="dashed") +
  geom_point(size=3) +
  scale_alpha_identity() +
  scale_colour_manual(values=c( "none"="darkgrey", "both"= "plum3","aCD3"="red", "flag"="#00B2FF")) +
  theme_bw(base_size=18) + 
  theme(aspect.ratio = 1) +
  #xlim(-7.5, 7.5) +
  ggtitle("all 6h") +
  xlab("log2FC in aCD3") +
  ylab("log2FC in flagellin") +
  theme(legend.text = element_text(size = 20),
        strip.background = element_blank()) +
  facet_wrap(~ Organism, scales = 'free', nrow=1) 

pdf('correlation_aCD3_vs_flag_6h.pdf', useDingbats = F, width=18, height=8)
all_aCD3_vs_flag_6h_gg
dev.off()




##### volcanoplots Figure 3C and Supplementary Figure 3B  #####

all <- bind_rows(BP_complete_df, BS_complete_df, CB_complete_df, PD_complete_df) %>%
  filter(!grepl('aCD36h_vs_flag6h|aCD324h_vs_flag24h', Comparison)) %>%
  mutate(Organism = str_extract(pattern='BP|BS|CB|PD', Comparison)) %>%
  #mutate(padj = coalesce(padj, 1)) %>%
  mutate(Significant = case_when(
    padj>0.05 | abs(log2FoldChange)<1 ~ 'not_sign',
    padj<0.05 & abs(log2FoldChange)>1 ~ 'sign'))

all_6h <-  all %>%
  filter(Time_point=='6h') 

 chaperone_6h <- all_6h %>%
   ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
   geom_point(aes(color = Significant), size=3) +
   geom_point(data = all_6h %>% filter(grepl(("Hsp|Chaperone|chaperone|ClpB"), all_6h$Gene_name)) %>% filter(Significant=='not_sign'), fill="goldenrod1", color="black", pch=21, size=3) +
   geom_point(data = all_6h %>% filter(grepl(("Hsp|Chaperone|chaperone|ClpB"), all_6h$Gene_name)) %>% filter(Significant=='sign'), fill="cornflowerblue", color="black", pch=21, size=3) +
   scale_color_manual(values = c("grey85", "grey70")) +
   theme_bw(base_size = 12) + theme(legend.position = "bottom") +
   ggtitle("Chaperones") +
   facet_wrap( Organism~Treatment, ncol=2) +
   theme(strip.background = element_blank(),
         strip.text.y = element_blank(),
         strip.text.x = element_blank(), 
         aspect.ratio = 1) 
 
pdf('gg_chaperone_6h.pdf', useDingbats = F)
gg_chaperone_6h
dev.off()
 
 
 gg_PUL_6h <- all_6h %>%
   ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
   geom_point(aes(color = Significant), size=3) +
   geom_point(data = all_6h %>% filter(grepl(("(PUL)"), all_6h$Gene_name)) %>% filter(Significant=='not_sign'), fill="goldenrod1", color="black", pch=21, size=3) +
   geom_point(data = all_6h %>% filter(grepl(("(PUL)"), all_6h$Gene_name)) %>% filter(Significant=='sign'), fill="cornflowerblue", color="black", pch=21, size=3) +
   scale_color_manual(values = c("grey85", "grey70")) +
   theme_bw(base_size = 12) + theme(legend.position = "bottom") +
   ggtitle("Polysaccharide Utilization Loci") +
   facet_wrap( Organism~Treatment, ncol=2) +
   theme(strip.background = element_blank(),
         strip.text.y = element_blank(),
         strip.text.x = element_blank(), 
         aspect.ratio = 1) 
 
 pdf("gg_PUL_6h.pdf", useDingbats = F)
 gg_PUL_6h
 dev.off()

 gg_ferric_iron_transporters_6h <- all_6h %>%
   ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
   geom_point(aes(color = Significant), size=3) +
   geom_point(data = all_6h %>% filter(grepl(("Ferric iron | ferric iron"), all_6h$Gene_name)) %>% filter(Significant=='not_sign'), fill="goldenrod1", color="black", pch=21, size=3) +
   geom_point(data = all_6h %>% filter(grepl(("Ferric iron | ferric iron"), all_6h$Gene_name)) %>% filter(Significant=='sign'), fill="cornflowerblue", color="black", pch=21, size=3) +
   scale_color_manual(values = c( "grey85", "grey70")) +
   theme_bw(base_size = 12) + theme(legend.position = "bottom") +
   ggtitle("Ferric iron uptake (ABC transporters and associated two-component systems)") +
   facet_wrap( Treatment~Organism, nrow=2) +
   theme(strip.background = element_blank(),
         strip.text.y = element_blank(),
         strip.text.x = element_blank(),
         aspect.ratio = 1)

 pdf("gg_ferric_iron_transporters_6h.pdf", useDingbats = F)
 gg_ferric_iron_transporters_6h
 dev.off()
 
 
 
 gg_catalases_6h <- all_6h %>%
   ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
   geom_point(aes(color = Significant), size=3) +
   geom_point(data = all_6h %>% filter(grepl(("Catalase|catalase"), all_6h$Gene_name)) %>% filter(Significant=='not_sign'), fill="goldenrod1", color="black", pch=21, size=3) +
   geom_point(data = all_6h %>% filter(grepl(("Catalase|catalase"), all_6h$Gene_name)) %>% filter(Significant=='sign'), fill="cornflowerblue", color="black", pch=21, size=3) +
   scale_color_manual(values = c( "grey85", "grey70")) +
   theme_bw(base_size = 12) + theme(legend.position = "bottom") +
   ggtitle("Catalases") +
   facet_wrap( Organism~Treatment, ncol=2) +
   theme(strip.background = element_blank(),
         strip.text.y = element_blank(),
         strip.text.x = element_blank(), 
         aspect.ratio = 1) 
 
 pdf("gg_catalases_6h.pdf", useDingbats = F)
 gg_catalases_6h
 dev.off()
 
 
 gg_thioredoxin_6h <- all_6h %>%
   ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
   geom_point(aes(color = Significant), size=3) +
   geom_point(data = all_6h %>% filter(grepl(("Thioredoxin|thioredoxin"), all_6h$Gene_name)) %>% filter(Significant=='not_sign'), fill="goldenrod1", color="black", pch=21, size=3) +
   geom_point(data = all_6h %>% filter(grepl(("Thioredoxin|thioredoxin"), all_6h$Gene_name)) %>% filter(Significant=='sign'), fill="cornflowerblue", color="black", pch=21, size=3) +
   scale_color_manual(values = c( "grey85", "grey70")) +
   theme_bw(base_size = 12) + theme(legend.position = "bottom") +
   ggtitle("Thioredoxin") +
   facet_wrap( Treatment~Organism, nrow=2) +
   theme(strip.background = element_blank(),
         strip.text.y = element_blank(),
         strip.text.x = element_blank(),
         aspect.ratio = 1) 
 
 pdf("gg_thioredoxin_6h.pdf", useDingbats = F)
 gg_thioredoxin_6h
 dev.off()
 
 

 
 

 
 
##### housekeeping genes selection for qPCR ####
 
 SB093_BP_HK_genes <- BP_norm_counts %>% mutate(variance=rowVars(BP_norm_counts)) %>% 
   mutate(meancounts=rowMeans(BP_norm_counts)) %>%
   rownames_to_column('Gene') %>% 
   filter(variance < quantile(variance, .02)) %>%# & meancounts > quantile(meancounts, .50)) %>%
   select(Gene, meancounts, variance) %>%
   mutate(Exp='SB093') %>% mutate(Organism='BP') %>% arrange(variance) 
 
 SB093_BS_HK_genes <- BS_norm_counts %>% mutate(variance=rowVars(BS_norm_counts)) %>% 
   mutate(meancounts=rowMeans(BS_norm_counts)) %>%
   rownames_to_column('Gene') %>% 
   filter(variance < quantile(variance, .02)) %>%# & meancounts > quantile(meancounts, .50)) %>%
   select(Gene, meancounts, variance) %>%
   mutate(Exp='SB093') %>% mutate(Organism='BS') %>% arrange(variance) 
 
 SB093_CB_HK_genes <- CB_norm_counts %>% mutate(variance=rowVars(CB_norm_counts)) %>% 
   mutate(meancounts=rowMeans(CB_norm_counts)) %>%
   rownames_to_column('Gene') %>% 
   filter(variance < quantile(variance, .02))%>%# & meancounts > quantile(meancounts, .50)) %>%
   select(Gene, meancounts, variance) %>%
   mutate(Exp='SB093') %>% mutate(Organism='CB') %>% arrange(variance) 
 
 SB093_PD_HK_genes <- PD_norm_counts %>% mutate(variance=rowVars(PD_norm_counts)) %>% 
   mutate(meancounts=rowMeans(PD_norm_counts)) %>%
   rownames_to_column('Gene') %>% 
   filter(variance < quantile(variance, .02)) %>%# & meancounts > quantile(meancounts, .50)) %>%
   select(Gene, meancounts, variance) %>%
   mutate(Exp='SB093') %>% mutate(Organism='PD') %>% arrange(variance) 
 
 bind_rows(SB093_BP_HK_genes, SB093_BS_HK_genes, SB093_CB_HK_genes, SB093_PD_HK_genes) %>%
   write.table('SB093_HK_genes.txt',quote=F, row.names=F, sep='\t')
##### Transcriptional activity levels to be used for Supplementary Figures 1E and 2E #####
 
 
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
 
 pdf('SB093 RNA_reads_overtime.pdf')
 All_reads %>%
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
 
##### get list of top 100 highly expressed genes in FASTA header format for Supplementary Figure 1 #####
 
 BP_counts_baseline <- BP_norm_counts %>% select(1,2,3) %>% mutate(mean_counts = rowMeans(BP_norm_counts)) %>% select(4) %>% mutate(Organism = 'BP') %>% rownames_to_column('Gene') %>% mutate(Gene=gsub('ID=', '>', Gene)) %>% mutate(Gene=gsub("(;.*)" , "", Gene)) %>% arrange( desc(mean_counts)) %>% head(100) %>% select(Gene) #%>% write.table('SB093_BP_counts_baseline_top100.txt', quote=FALSE, row.names=F, col.names=F)
 BS_counts_baseline <- BS_norm_counts %>% select(1,2,3) %>% mutate(mean_counts = rowMeans(BS_norm_counts)) %>% select(4) %>% mutate(Organism = 'BS') %>% rownames_to_column('Gene') %>% mutate(Gene=gsub('ID=', '>', Gene)) %>% mutate(Gene=gsub("(;.*)" , "", Gene)) %>% arrange( desc(mean_counts)) %>% head(100) %>% select(Gene) #%>% write.table('SB093_BS_counts_baseline_top100.txt', quote=FALSE, row.names=F, col.names=F)
 CB_counts_baseline <- CB_norm_counts %>% select(1,2,3) %>% mutate(mean_counts = rowMeans(CB_norm_counts)) %>% select(4) %>% mutate(Organism = 'CB') %>% rownames_to_column('Gene') %>% mutate(Gene=gsub('ID=', '>', Gene)) %>% mutate(Gene=gsub("(;.*)" , "", Gene)) %>% arrange( desc(mean_counts)) %>% head(100) %>% select(Gene) #%>% write.table('SB093_CB_counts_baseline_top100.txt', quote=FALSE, row.names=F, col.names=F)
 PD_counts_baseline <- PD_norm_counts %>% select(1,2,3) %>% mutate(mean_counts = rowMeans(PD_norm_counts)) %>% select(4) %>% mutate(Organism = 'PD') %>% rownames_to_column('Gene') %>% mutate(Gene=gsub('ID=', '>', Gene)) %>% mutate(Gene=gsub("(;.*)" , "", Gene)) %>% arrange( desc(mean_counts)) %>% head(100) %>% select(Gene) #%>% write.table('SB093_PD_counts_baseline_top100.txt', quote=FALSE, row.names=F, col.names=F)
 


 