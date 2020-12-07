### load libraries ####
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(Rsubread)
library(clusterProfiler)
library(stringr)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(grDevices)
library(genefilter)
library(xlsx)
library(tidyr)
library(openxlsx)
library(githubinstall)
library(yingtools2)
library(phyloseq)
library(reshape2)
library(tibble)
library(ggrepel)
library(gplots)
library(VennDiagram)
library(grid)
library(gridExtra)
library(psych)
library(viridis)  
library(ggsci)
library(dplyr)
library(egg)
library(Pigengene)

setwd('/Users/becattis/Desktop/CHM revision/Scripts to be uploaded/CBBP-transcriptomics/3mix vs 4mix/')

### prepare .gff files 

BP <- read.delim("BP.gff", header=T)
BS <- read.delim("BS.gff", header=T)
CB <- read.delim("CB.gff", header=T)
PD <- read.delim("PD.gff", header=T)

#### if wanting to start from original fastq files and performing the alignments, then the following code, adapted for each bacterium, must be run; 
#### otherwise skip this and read in the counts files

# CB_list <- list.files("/Users/becattis/Desktop/SB172_SB151_SB172_SB180/sequences/SAM_BAM_ALIGNMENTS/SB172", pattern = "CB.only_aligned.bam")
# CB_FtCounts <- featureCounts(CB_list, annot.ext = CB, isGTFAnnotationFile=FALSE, largestOverlap =T, countChimericFragments=TRUE)
# CB_Counts <- as.data.frame(CB_FtCounts$counts)
# colnames(CB_Counts) <- str_replace(colnames(CB_Counts), 'SB172.37.CB.only_aligned.bam', 'CB_3mix_0h_1') %>% str_replace('SB172.38.CB.only_aligned.bam', 'CB_3mix_0h_2') %>% str_replace('SB172.39.CB.only_aligned.bam', 'CB_3mix_0h_3') %>%
#  str_replace('SB172.40.CB.only_aligned.bam', 'CB_3mix_6h_4') %>% str_replace('SB172.41.CB.only_aligned.bam', 'CB_3mix_6h_5') %>% str_replace('SB172.42.CB.only_aligned.bam', 'CB_3mix_6h_6') %>% str_replace('SB172.43.CB.only_aligned.bam', 'CB_4mix_6h_7') %>%
#  str_replace('SB172.44.CB.only_aligned.bam', 'CB_4mix_6h_8') %>% str_replace('SB172.45.CB.only_aligned.bam', 'CB_4mix_6h_9') %>% str_replace('SB172.46.CB.only_aligned.bam', 'CB_4mix_0h_10') %>% str_replace('SB172.47.CB.only_aligned.bam', 'CB_4mix_0h_11') %>%
 # str_replace('SB172.48.CB.only_aligned.bam', 'CB_4mix_0h_12')
#CB_Counts$feature <- row.names(CB_Counts)
#CB_Counts <- CB_Counts %>% filter(!grepl(pattern="tRNA-[A-Z, a-z]{3}-[A-Z, a-z]{3}" , CB_Counts$feature))
#CB_Counts <- CB_Counts %>% filter(!grepl(pattern="LSU rRNA|SSU rRNA|5S rRNA" , CB_Counts$feature))
#CB_Counts <- CB_Counts %>% filter(!grepl(pattern="repeat region" , CB_Counts$feature))
#rownames(CB_Counts) <- CB_Counts$feature 
#CB_Counts$feature <- NULL


BP_Counts <- read.delim2('SB172_BP_Counts.txt', header = T, sep=' ') 
BS_Counts <- read.delim2('SB172_BS_Counts.txt', header = T, sep=' ') 
CB_Counts <- read.delim2('SB172_CB_Counts.txt', header = T, sep=' ') 
PD_Counts <- read.delim2('SB172_PD_Counts.txt', header = T, sep=' ') 


#### DeSeq2 and PCA BP ####

BP_colData <- as.data.frame(matrix(c("0h_3mix", "0h_3mix", "0h_3mix","6h_3mix","6h_3mix", "6h_3mix", "6h_4mix", "6h_4mix", "6h_4mix", "0h_4mix","0h_4mix", "0h_4mix" ),nrow=12,ncol=1))
row.names(BP_colData)<-c("BP_3mix_0h_1", "BP_3mix_0h_2", "BP_3mix_0h_3", "BP_3mix_6h_4", "BP_3mix_6h_5",  "BP_3mix_6h_6", "BP_4mix_6h_7", "BP_4mix_6h_8", "BP_4mix_6h_9", "BP_4mix_0h_10", "BP_4mix_0h_11", "BP_4mix_0h_12")
colnames(BP_colData)<-c('BP_condition')
BP_countData <-  as.matrix(BP_Counts, rownames.force=T )

BP_dds <- DESeqDataSetFromMatrix(countData=BP_countData, colData=BP_colData, design= ~ BP_condition)
BP_dds <- BP_dds[rowSums(counts(BP_dds)) > 10, ] 
BP_dds <- DESeq(BP_dds, betaPrior = T, fitType = 'parametric') 

BP_3mix_6h <- as.data.frame(results(BP_dds, contrast=c("BP_condition","6h_3mix","0h_3mix"))) %>% rownames_to_column('Gene') %>% mutate (Gene2=Gene) %>% mutate(Organism="BP", Comparison="3mix_6h") %>% column_to_rownames('Gene2')
BP_4mix_6h <- as.data.frame(results(BP_dds, contrast=c("BP_condition","6h_4mix","0h_4mix"))) %>% rownames_to_column('Gene') %>% mutate (Gene2=Gene) %>% mutate(Organism="BP", Comparison="4mix_6h") %>% column_to_rownames('Gene2')
BP_3mix_vs_4mix_6h <- as.data.frame(results(BP_dds, contrast=c("BP_condition","6h_3mix","6h_4mix"))) %>% rownames_to_column('Gene') %>% mutate (Gene2=Gene) %>% mutate(Organism="BP", Comparison="3mix_vs_4mix_6h") %>% column_to_rownames('Gene2')
BP_3mix_vs_4mix_0h <- as.data.frame(results(BP_dds, contrast=c("BP_condition","0h_3mix","0h_4mix"))) %>% rownames_to_column('Gene') %>% mutate (Gene2=Gene) %>% mutate(Organism="BP", Comparison="3mix_vs_4mix_0h") %>% column_to_rownames('Gene2')

BP_alldata <- bind_rows(BP_3mix_6h, BP_4mix_6h, BP_3mix_vs_4mix_0h, BP_3mix_vs_4mix_6h) %>% 
  mutate(Gene_name = gsub("(.*);Name=", "", Gene)) %>% mutate(Gene_name = gsub("(\\(E.*)", "", Gene_name)) %>%
  mutate(Gene_number = gsub("(;.*)" , "", Gene)) %>% mutate(Gene_number = gsub("ID=fig\\|33035.10.", "", Gene_number)) %>% 
  mutate(Gene_number = gsub("\\|", "", Gene_number)) %>% mutate(Gene_name = gsub( "%2C", ",", Gene_name)) %>%
  mutate(Gene_name = ifelse(grepl("hypothetical protein", Gene), paste0(Gene_name, " ", Gene_number), Gene_name))

BP_rld <- rlog(BP_dds, blind=FALSE)
BP_alone_norm_counts <- as.data.frame(assay(BP_rld))

LRT_BP <- as.data.frame(results(DESeq(BP_dds, test="LRT", reduced=~1)))
LRT_BP$Gene <- rownames(LRT_BP) 
LRT_BP_list <-  LRT_BP %>% filter (padj<0.05, abs(log2FoldChange) >1.5) %>% pull(Gene)
BP_alone_norm_counts$Gene <- rownames(BP_alone_norm_counts)
BP_alone_norm_counts_lrt <- BP_alone_norm_counts %>% filter(Gene %in% LRT_BP_list) ## NOTE: simply removing this will give you AUTOPLOT on ALL VALUES!
BP_alone_norm_counts_lrt$Gene <- NULL
LRT_BP_prcomp <- prcomp(x=t(BP_alone_norm_counts_lrt), center = T, scale=T) #run PCA on the selected table of counts 
LRT_BP_summ <- summary(LRT_BP_prcomp)
LRT_BP_prcomp_graph <- as.data.frame(LRT_BP_prcomp$x[,c(1,2)])
LRT_BP_prcomp_graph$Group <- rep(c("3mix_0h", "3mix_6h", "4mix_6h", "4mix_0h"), c(3,3,3,3))

t_BP_norm_counts_LRT <- as.data.frame(t(BP_alone_norm_counts_lrt)) %>% 
  mutate(Time_point=str_extract(rownames(LRT_BP_prcomp$x), pattern="0h|6h"))  %>% 
  mutate(Community=str_extract(rownames(LRT_BP_prcomp$x), pattern="3mix|4mix")) %>%
  mutate(Group=paste0(Community, ' ', Time_point))

Figure_4C_BP <- autoplot(prcomp(t_BP_norm_counts_LRT[-rev(seq_len(ncol(t_BP_norm_counts_LRT)))[1:3]]), data = t_BP_norm_counts_LRT) +
  geom_point(size=6,aes(shape=Community,color=Time_point)) +
  theme_bw(base_size=18) + 
  scale_color_manual(values=c("24h"="steelblue2","6h"="salmon","0h"="gray33")) +
  ggtitle("SB172 B. producta PCA") +
  theme(text = element_text(size=16)) +
  theme(aspect.ratio = 1) +
  stat_chull(aes(color = Time_point, fill = Group), 
             alpha = 0.1, geom = "polygon")


#### DeSeq2 BS ####

BS_colData <- as.data.frame(matrix(c( "6h_4mix", "6h_4mix", "6h_4mix", "0h_4mix","0h_4mix", "0h_4mix" ),nrow=6,ncol=1))
row.names(BS_colData)<-c( "BS_4mix_6h_7", "BS_4mix_6h_8", "BS_4mix_6h_9", "BS_4mix_0h_10", "BS_4mix_0h_11", "BS_4mix_0h_12")
colnames(BS_colData)<-c('BS_condition')
BS_countData <-  as.matrix(BS_Counts[,c(7,8,9,10,11,12)], rownames.force=T )

BS_dds <- DESeqDataSetFromMatrix(countData=BS_countData, colData=BS_colData, design= ~ BS_condition)
BS_dds <- BS_dds[rowSums(counts(BS_dds)) > 1, ] 
BS_dds <- DESeq(BS_dds,  betaPrior = TRUE,  fitType = 'parametric') 

BS_4mix_6h <- as.data.frame(results(BS_dds, contrast=c("BS_condition","6h_4mix","0h_4mix"))) %>% rownames_to_column('Gene') %>% mutate (Gene2=Gene) %>% mutate(Organism="BS", Comparison="4mix_6h") %>% column_to_rownames('Gene2')


#### DeSeq2 and PCA CB ####

CB_colData <- as.data.frame(matrix(c("0h_3mix", "0h_3mix", "0h_3mix","6h_3mix","6h_3mix", "6h_3mix", "6h_4mix", "6h_4mix", "6h_4mix", "0h_4mix","0h_4mix", "0h_4mix" ),nrow=12,ncol=1))
row.names(CB_colData)<-c("CB_3mix_0h_1", "CB_3mix_0h_2", "CB_3mix_0h_3", "CB_3mix_6h_4", "CB_3mix_6h_5",  "CB_3mix_6h_6", "CB_4mix_6h_7", "CB_4mix_6h_8", "CB_4mix_6h_9", "CB_4mix_0h_10", "CB_4mix_0h_11", "CB_4mix_0h_12")
colnames(CB_colData)<-c('CB_condition')
CB_countData <-  as.matrix(CB_Counts, rownames.force=T )

CB_dds <- DESeqDataSetFromMatrix(countData=CB_countData, colData=CB_colData, design= ~ CB_condition)
CB_dds <- CB_dds[rowSums(counts(CB_dds)) > 10, ] 
CB_dds <- DESeq(CB_dds, betaPrior = TRUE, fitType = 'parametric') 

CB_3mix_6h <- as.data.frame(results(CB_dds, contrast=c("CB_condition","6h_3mix","0h_3mix"))) %>% rownames_to_column('Gene') %>% mutate (Gene2=Gene) %>% mutate(Organism="CB", Comparison="3mix_6h") %>% column_to_rownames('Gene2')
CB_4mix_6h <- as.data.frame(results(CB_dds, contrast=c("CB_condition","6h_4mix","0h_4mix"))) %>% rownames_to_column('Gene') %>% mutate (Gene2=Gene) %>% mutate(Organism="CB", Comparison="4mix_6h") %>% column_to_rownames('Gene2')
CB_3mix_vs_4mix_6h <- as.data.frame(results(CB_dds, contrast=c("CB_condition","6h_3mix","6h_4mix"))) %>% rownames_to_column('Gene') %>% mutate (Gene2=Gene) %>% mutate(Organism="CB", Comparison="3mix_vs_4mix_6h") %>% column_to_rownames('Gene2')
CB_3mix_vs_4mix_0h <- as.data.frame(results(CB_dds, contrast=c("CB_condition","0h_3mix","0h_4mix"))) %>% rownames_to_column('Gene') %>% mutate (Gene2=Gene) %>% mutate(Organism="CB", Comparison="3mix_vs_4mix_0h") %>% column_to_rownames('Gene2')

CB_alldata <- bind_rows(CB_3mix_6h, CB_4mix_6h, CB_3mix_vs_4mix_0h, CB_3mix_vs_4mix_6h) %>% 
  mutate(Gene_name = gsub("(.*);Name=", "", Gene)) %>% mutate(Gene_name = gsub("(\\(E.*)", "", Gene_name)) %>%
  mutate(Gene_number = gsub("(;.*)" , "", Gene)) %>% mutate(Gene_number = gsub("ID=fig\\|33035.10.", "", Gene_number)) %>% 
  mutate(Gene_number = gsub("\\|", "", Gene_number)) %>% mutate(Gene_name = gsub( "%2C", ",", Gene_name)) %>%
  mutate(Gene_name = ifelse(grepl("hypothetical protein", Gene), paste0(Gene_name, " ", Gene_number), Gene_name))

CB_rld <- rlog(CB_dds, blind=FALSE)
CB_alone_norm_counts <- as.data.frame(assay(CB_rld))

LRT_CB <- as.data.frame(results(DESeq(CB_dds, test="LRT", reduced=~1)))
LRT_CB$Gene <- rownames(LRT_CB) 
LRT_CB_list <-  LRT_CB %>% filter (padj<0.05, abs(log2FoldChange) >1.5) %>% pull(Gene)
CB_alone_norm_counts$Gene <- rownames(CB_alone_norm_counts)
CB_alone_norm_counts_lrt <- CB_alone_norm_counts %>% filter(Gene %in% LRT_CB_list)
CB_alone_norm_counts_lrt$Gene <- NULL
LRT_CB_prcomp <- prcomp(x=t(CB_alone_norm_counts_lrt), center = T, scale=T) #run PCA on the selected table of counts 
LRT_CB_summ <- summary(LRT_CB_prcomp)
LRT_CB_prcomp_graph <- as.data.frame(LRT_CB_prcomp$x[,c(1,2)])
LRT_CB_prcomp_graph$Group <- rep(c("3mix_0h", "3mix_6h", "4mix_6h", "4mix_0h"), c(3,3,3,3))

t_CB_norm_counts_LRT <- as.data.frame(t(CB_alone_norm_counts_lrt)) %>% 
  mutate(Time_point=str_extract(rownames(LRT_CB_prcomp$x), pattern="0h|6h"))  %>% 
  mutate(Community=str_extract(rownames(LRT_CB_prcomp$x), pattern="3mix|4mix")) %>%
  mutate(Group=paste0(Community, ' ', Time_point))

Figure_4C_CB <- autoplot(prcomp(t_CB_norm_counts_LRT[-rev(seq_len(ncol(t_CB_norm_counts_LRT)))[1:3]]), data = t_CB_norm_counts_LRT) +
  geom_point(size=6,aes(shape=Community,color=Time_point)) +
  theme_bw(base_size=18) + 
  scale_color_manual(values=c("24h"="steelblue2","6h"="salmon","0h"="gray33")) +
  ggtitle("SB172 C. bolteae PCA") +
  theme(text = element_text(size=16)) +
  theme(aspect.ratio = 1) +
  stat_chull(aes(color = Time_point, fill = Group), 
             alpha = 0.1, geom = "polygon")


#### DeSeq2 and PCA PD ####

PD_colData <- as.data.frame(matrix(c("0h_3mix", "0h_3mix", "0h_3mix","6h_3mix","6h_3mix", "6h_3mix", "6h_4mix", "6h_4mix", "6h_4mix", "0h_4mix","0h_4mix", "0h_4mix" ),nrow=12,ncol=1))
row.names(PD_colData)<-c("PD_3mix_0h_1", "PD_3mix_0h_2", "PD_3mix_0h_3", "PD_3mix_6h_4", "PD_3mix_6h_5",  "PD_3mix_6h_6", "PD_4mix_6h_7", "PD_4mix_6h_8", "PD_4mix_6h_9", "PD_4mix_0h_10", "PD_4mix_0h_11", "PD_4mix_0h_12")
colnames(PD_colData)<-c('PD_condition')
PD_countData <-  as.matrix(PD_Counts, rownames.force=T )

PD_dds <- DESeqDataSetFromMatrix(countData=PD_countData, colData=PD_colData, design= ~ PD_condition)
PD_dds <- PD_dds[rowSums(counts(PD_dds)) > 10, ] 
PD_dds <- DESeq(PD_dds, betaPrior = TRUE, fitType = 'parametric') 

PD_3mix_6h <- as.data.frame(results(PD_dds, contrast=c("PD_condition","6h_3mix","0h_3mix"))) %>% rownames_to_column('Gene') %>% mutate (Gene2=Gene) %>% mutate(Organism="PD", Comparison="3mix_6h") %>% column_to_rownames('Gene2')
PD_4mix_6h <- as.data.frame(results(PD_dds, contrast=c("PD_condition","6h_4mix","0h_4mix"))) %>% rownames_to_column('Gene') %>% mutate (Gene2=Gene) %>% mutate(Organism="PD", Comparison="4mix_6h") %>% column_to_rownames('Gene2')
PD_3mix_vs_4mix_6h <- as.data.frame(results(PD_dds, contrast=c("PD_condition","6h_3mix","6h_4mix"))) %>% rownames_to_column('Gene') %>% mutate (Gene2=Gene) %>% mutate(Organism="PD", Comparison="3mix_vs_4mix_6h") %>% column_to_rownames('Gene2')
PD_3mix_vs_4mix_0h <- as.data.frame(results(PD_dds, contrast=c("PD_condition","0h_3mix","0h_4mix"))) %>% rownames_to_column('Gene') %>% mutate (Gene2=Gene) %>% mutate(Organism="PD", Comparison="3mix_vs_4mix_0h") %>% column_to_rownames('Gene2')

PD_alldata <- bind_rows(PD_3mix_6h, PD_4mix_6h, PD_3mix_vs_4mix_0h, PD_3mix_vs_4mix_6h) %>% 
  mutate(Gene_name = gsub("(.*);Name=", "", Gene)) %>% mutate(Gene_name = gsub("(\\(E.*)", "", Gene_name)) %>%
  mutate(Gene_number = gsub("(;.*)" , "", Gene)) %>% mutate(Gene_number = gsub("ID=fig\\|33035.10.", "", Gene_number)) %>% 
  mutate(Gene_number = gsub("\\|", "", Gene_number)) %>% mutate(Gene_name = gsub( "%2C", ",", Gene_name)) %>%
  mutate(Gene_name = ifelse(grepl("hypothetical protein", Gene), paste0(Gene_name, " ", Gene_number), Gene_name))

PD_rld <- rlog(PD_dds, blind=FALSE)
PD_alone_norm_counts <- as.data.frame(assay(PD_rld))

LRT_PD <- as.data.frame(results(DESeq(PD_dds, test="LRT", reduced=~1)))
LRT_PD$Gene <- rownames(LRT_PD) 
LRT_PD_list <-  LRT_PD %>% filter (padj<0.05, abs(log2FoldChange) >1.5) %>% pull(Gene)
PD_alone_norm_counts$Gene <- rownames(PD_alone_norm_counts)
PD_alone_norm_counts_lrt <- PD_alone_norm_counts %>% filter(Gene %in% LRT_PD_list)
PD_alone_norm_counts_lrt$Gene <- NULL
LRT_PD_prcomp <- prcomp(x=t(PD_alone_norm_counts_lrt), center = T, scale=T) #run PCA on the selected table of counts 
LRT_PD_summ <- summary(LRT_PD_prcomp)
LRT_PD_prcomp_graph <- as.data.frame(LRT_PD_prcomp$x[,c(1,2)])
LRT_PD_prcomp_graph$Group <- rep(c("3mix_0h", "3mix_6h", "4mix_6h", "4mix_0h"), c(3,3,3,3))

t_PD_norm_counts_LRT <- as.data.frame(t(PD_alone_norm_counts_lrt)) %>% 
  mutate(Time_point=str_extract(rownames(LRT_PD_prcomp$x), pattern="0h|6h"))  %>% 
  mutate(Community=str_extract(rownames(LRT_PD_prcomp$x), pattern="3mix|4mix")) %>%
  mutate(Group=paste0(Community, ' ', Time_point))

Figure_4C_PD <- autoplot(prcomp(t_PD_norm_counts_LRT[-rev(seq_len(ncol(t_PD_norm_counts_LRT)))[1:3]]), data = t_PD_norm_counts_LRT) +
  geom_point(size=6,aes(shape=Community,color=Time_point)) +
  theme_bw(base_size=18) + 
  scale_color_manual(values=c("24h"="steelblue2","6h"="salmon","0h"="gray33")) +
  ggtitle("SB172 P. distasonis PCA") +
  theme(text = element_text(size=16)) +
  theme(aspect.ratio = 1) +
  stat_chull(aes(color = Time_point, fill = Group), 
             alpha = 0.1, geom = "polygon")




#### Heatmaps Supplementary Figure 4 ####
my_colour = list(
  Time_point = c('0h'="#868686FF", '6h'="#EFC000FF"),
  Consortium = c('3mix'='bisque2','4mix'='plum4')
)

### BP

select_genes_manual_BP <- LRT_BP %>% filter (padj<0.05, abs(log2FoldChange) >1) %>% pull(Gene)
select_genes_manual_BP <- unique(select_genes_manual_BP)
mat_manual_BP <- assay(BP_rld)[select_genes_manual_BP,]
topVarGenes_BP <- head(order(rowVars(mat_manual_BP), decreasing = TRUE), 200)
mat_manual_BP <- mat_manual_BP[ topVarGenes_BP, ]
mat_manual_BP <- mat_manual_BP - rowMeans(mat_manual_BP)
mat_manual_BP <- as.data.frame(mat_manual_BP)
mat_manual_BP$temp <- rownames(mat_manual_BP)
mat_manual_BP$temp <- gsub("(.*);Name=", "", mat_manual_BP$temp)
mat_manual_BP$temp <- gsub("(\\(E.*)" , "", mat_manual_BP$temp)
mat_manual_BP$n <- seq.int(nrow(mat_manual_BP)) # necessary as R does not let you have rows with identical names, so need to differentiate genes with identical name
mat_manual_BP$temp = paste(mat_manual_BP$n, mat_manual_BP$temp, sep="_") # merge character and n columns
rownames(mat_manual_BP) <- mat_manual_BP$temp
mat_manual_BP$temp <- NULL
mat_manual_BP$n <- NULL
mat_manual_BP <- mat_manual_BP %>% subset(select=c(1,2,3,10,11,12,4,5,6,7,8,9)) #to order samples
mat_manual_BP <- as.matrix(mat_manual_BP)
mat_manual_BP_factors <- data.frame(Sample=colnames(mat_manual_BP)) %>%
  mutate(Consortium=str_extract(Sample, pattern="3mix|4mix"))  %>% 
  mutate(Time_point=str_extract(Sample, pattern="0h|6h"))  %>% 
  column_to_rownames('Sample')
 
pheatmap(mat_manual_BP, cluster_cols = F, show_rownames = F, show_colnames = F, border_color = NA, color = bluered(50), 
         annotation_col = mat_manual_BP_factors, annotation_colors = my_colour, main = "BP 3mix/6mix 6h top DEG" )


### CB

select_genes_manual_CB <- LRT_CB %>% filter (padj<0.05, abs(log2FoldChange) >1) %>% pull(Gene)
select_genes_manual_CB <- unique(select_genes_manual_CB)

mat_manual_CB <- assay(CB_rld)[select_genes_manual_CB,]
topVarGenes_CB <- head(order(rowVars(mat_manual_CB), decreasing = TRUE), 200)
mat_manual_CB <- mat_manual_CB[ topVarGenes_CB, ]
mat_manual_CB <- mat_manual_CB - rowMeans(mat_manual_CB)
mat_manual_CB <- as.data.frame(mat_manual_CB)
mat_manual_CB$temp <- rownames(mat_manual_CB)
mat_manual_CB$temp <- gsub("(.*);Name=", "", mat_manual_CB$temp) 
mat_manual_CB$temp <- gsub("(\\(E.*)" , "", mat_manual_CB$temp) 
mat_manual_CB$n <- seq.int(nrow(mat_manual_CB)) #necessary as R does not let you have rows with identical names, so need to differentiate genes with identical name
mat_manual_CB$temp = paste(mat_manual_CB$n, mat_manual_CB$temp, sep="_") # merge character and n columns
rownames(mat_manual_CB) <- mat_manual_CB$temp
mat_manual_CB$temp <- NULL
mat_manual_CB$n <- NULL
mat_manual_CB <- mat_manual_CB %>% subset(select=c(1,2,3,10,11,12,4,5,6,7,8,9)) # to order samples
mat_manual_CB <- as.matrix(mat_manual_CB)
mat_manual_CB_factors <- data.frame(Sample=colnames(mat_manual_CB)) %>%
  mutate(Consortium=str_extract(Sample, pattern="3mix|4mix"))  %>% 
  mutate(Time_point=str_extract(Sample, pattern="0h|6h"))  %>% 
  column_to_rownames('Sample')

pheatmap(mat_manual_CB, cluster_cols = F, show_rownames = F, show_colnames = F, border_color = NA, color = bluered(50), 
         annotation_col = mat_manual_CB_factors, annotation_colors = my_colour, main = "CB 3mix/6mix 6h top DEG" )


#### PD
select_genes_manual_PD <- LRT_PD %>% filter (padj<0.05, abs(log2FoldChange) >1) %>% pull(Gene)
select_genes_manual_PD <- unique(select_genes_manual_PD)

mat_manual_PD <- assay(PD_rld)[select_genes_manual_PD,]
topVarGenes_PD <- head(order(rowVars(mat_manual_PD), decreasing = TRUE), 200)
mat_manual_PD <- mat_manual_PD[ topVarGenes_PD, ]
mat_manual_PD <- mat_manual_PD - rowMeans(mat_manual_PD)
mat_manual_PD <- as.data.frame(mat_manual_PD)
mat_manual_PD$temp <- rownames(mat_manual_PD)
mat_manual_PD$temp <- gsub("(.*);Name=", "", mat_manual_PD$temp) 
mat_manual_PD$temp <- gsub("(\\(E.*)" , "", mat_manual_PD$temp) 
mat_manual_PD$n <- seq.int(nrow(mat_manual_PD)) # necessary as R does not let you have rows with identical names, so need to differentiate genes with identical name
mat_manual_PD$temp = paste(mat_manual_PD$n, mat_manual_PD$temp, sep="_") # merge character and n columns
rownames(mat_manual_PD) <- mat_manual_PD$temp
mat_manual_PD$temp <- NULL
mat_manual_PD$n <- NULL
mat_manual_PD <- mat_manual_PD %>% subset(select=c(1,2,3,10,11,12,4,5,6,7,8,9)) # to order samples
mat_manual_PD <- as.matrix(mat_manual_PD)
mat_manual_PD_factors <- data.frame(Sample=colnames(mat_manual_PD)) %>%
  mutate(Consortium=str_extract(Sample, pattern="3mix|4mix"))  %>% 
  mutate(Time_point=str_extract(Sample, pattern="0h|6h"))  %>% 
  column_to_rownames('Sample')


pheatmap(mat_manual_PD, cluster_cols = F, show_rownames = F, show_colnames = F, border_color = NA, color = bluered(50), 
         annotation_col = mat_manual_PD_factors, annotation_colors = my_colour, main = "PD 3mix/6mix 6h top DEG" )


#### VENN DIAGRAMS FIG. 4 ####


BP_6h_shared <-   get.venn.partitions(list(BP3mix6h=(BP_3mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene)), 
                                           BP4mix6h =(BP_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene))))

Figure_4B_BP <- venn.diagram(list(BP3mix6h=(BP_3mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene)), 
                             BP4mix6h =(BP_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene))),
                        filename = NULL, output = T, 
                        #  height =2000 , 
                        #  width = 2000, 
                        #  resolution = 600,
                        lwd = 2,
                        lty = 'blank',
                        fill = c("black", "dodgerblue3"),
                        cex = 1,
                        cat.cex = 1,
                        cat.dist = c(0.05, 0.05),
                        cat.default.pos = 'outer',
                        cat.pos= c(-20, 20),
                        fontfamily = 'serif',
                        cat.fontfamily = 'serif',
                        cat.col = c("black", "dodgerblue3"), main= 'CB', sub= 'hyper'
)





#### CB - this gives you DEG genes that are in DE between 3 and 4 mix both at 0 and 6h  
CB_comparison <- intersect(
  (CB_3mix_vs_4mix_0h %>% filter(padj<0.05 & abs(log2FoldChange >1)) %>% pull(Gene)),
  (CB_3mix_vs_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange >1)) %>% pull(Gene))
)

CB_comparison_unique <- c(setdiff(
  (CB_3mix_vs_4mix_0h %>% filter(padj<0.05 & abs(log2FoldChange >1)) %>% pull(Gene)),
  (CB_3mix_vs_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange >1)) %>% pull(Gene))
), setdiff(
  (CB_3mix_vs_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange >1)) %>% pull(Gene)),
  (CB_3mix_vs_4mix_0h %>% filter(padj<0.05 & abs(log2FoldChange >1)) %>% pull(Gene))))


CB_3vs4_shared <-   get.venn.partitions(list(CB_3vs4_0h=(CB_3mix_vs_4mix_0h %>% filter(padj<0.05 & abs(log2FoldChange >1)) %>% pull(Gene)), 
                                             CB_3vs4_6h =(CB_3mix_vs_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange >1)) %>% pull(Gene))))

CB_comparison_6h <- intersect(
  (CB_3mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% arrange(padj)  %>% pull(Gene)),
  (CB_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% arrange(padj)  %>% pull(Gene))
)

CB_comparison_6h_names <- c((CB_3mix_6h %>% arrange(padj) %>% filter(Gene %in%CB_comparison_6h) %>% head(10) %>% pull(Gene)),
                            (CB_4mix_6h %>% arrange(padj) %>% filter(Gene %in%CB_comparison_6h) %>% head(10) %>% pull(Gene)),
                            (CB_3mix_6h %>% arrange(desc(abs(log2FoldChange))) %>% filter(Gene %in%CB_comparison_6h) %>% head(10) %>% pull(Gene)),
                            (CB_4mix_6h %>% arrange(desc(abs(log2FoldChange))) %>% filter(Gene %in%CB_comparison_6h) %>% head(10) %>% pull(Gene))
)

CB_comparison_6h_unique <- c(setdiff(
  (CB_3mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene)),
  (CB_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene))
), setdiff(
  (CB_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene)),
  (CB_3mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene))))



CB_comparison_6h_unique_names <- c((CB_3mix_6h %>% arrange((padj), Gene) %>% filter(Gene %in%  CB_comparison_6h_unique) %>% head(5) %>% pull(Gene)), 
                                   (CB_4mix_6h %>% arrange((padj), Gene) %>% filter(Gene %in%  CB_comparison_6h_unique) %>% head(5) %>% pull(Gene)), 
                                   (CB_3mix_6h %>% arrange(desc(abs(log2FoldChange)), Gene) %>% filter(Gene %in%  CB_comparison_6h_unique) %>% head(5) %>% pull(Gene)), 
                                   (CB_4mix_6h %>% arrange(desc(abs(log2FoldChange)), Gene) %>% filter(Gene %in%  CB_comparison_6h_unique) %>% head(5) %>% pull(Gene)))


CB_6h_shared <-   get.venn.partitions(list(CB3mix6h=(CB_3mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene)), 
                                           CB4mix6h =(CB_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene))))

CB_venn <- venn.diagram(list(CB3mix6h=(CB_3mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene)), 
                  CB4mix6h =(CB_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene))),
             filename = NULL, output = T, 
           #  height =2000 , 
           #  width = 2000, 
           #  resolution = 600,
             lwd = 2,
             lty = 'blank',
             fill = c("black", "dodgerblue3"),
             cex = 1,
             cat.cex = 1,
             cat.dist = c(0.05, 0.05),
             cat.default.pos = 'outer',
             cat.pos= c(-20, 20),
             fontfamily = 'serif',
             cat.fontfamily = 'serif',
             cat.col = c("black", "dodgerblue3"), main= 'CB', sub= 'hyper'
)


#### PD - this gives you DEG genes that are in DE between 3 and 4 mix both at 0 and 6h  

PD_comparison <- intersect(
  (PD_3mix_vs_4mix_0h %>% filter(padj<0.05 & abs(log2FoldChange >1)) %>% pull(Gene)),
  (PD_3mix_vs_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange >1)) %>% pull(Gene))
)

PD_comparison_unique <- c(setdiff(
  (PD_3mix_vs_4mix_0h %>% filter(padj<0.05 & abs(log2FoldChange >1)) %>% pull(Gene)),
  (PD_3mix_vs_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange >1)) %>% pull(Gene))
), setdiff(
  (PD_3mix_vs_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange >1)) %>% pull(Gene)),
  (PD_3mix_vs_4mix_0h %>% filter(padj<0.05 & abs(log2FoldChange >1)) %>% pull(Gene))))


PD_3vs4_shared <-   get.venn.partitions(list(PD_3vs4_0h=(PD_3mix_vs_4mix_0h %>% filter(padj<0.05 & abs(log2FoldChange >1)) %>% pull(Gene)), 
                                             PD_3vs4_6h =(PD_3mix_vs_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange >1)) %>% pull(Gene))))


PD_comparison_6h <- intersect(
  (PD_3mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% arrange(padj)  %>% pull(Gene)),
  (PD_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% arrange(padj)  %>% pull(Gene))
)

PD_comparison_6h_names <- c((PD_3mix_6h %>% arrange(padj) %>% filter(Gene %in%PD_comparison_6h) %>% head(10) %>% pull(Gene)),
                            (PD_4mix_6h %>% arrange(padj) %>% filter(Gene %in%PD_comparison_6h) %>% head(10) %>% pull(Gene)),
                            (PD_3mix_6h %>% arrange(desc(abs(log2FoldChange))) %>% filter(Gene %in%PD_comparison_6h) %>% head(10) %>% pull(Gene)),
                            (PD_4mix_6h %>% arrange(desc(abs(log2FoldChange))) %>% filter(Gene %in%PD_comparison_6h) %>% head(10) %>% pull(Gene))
)

PD_comparison_6h_unique <- c(setdiff(
  (PD_3mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene)),
  (PD_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene))
), setdiff(
  (PD_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene)),
  (PD_3mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene))))



PD_comparison_6h_unique_names <- c((PD_3mix_6h %>% arrange((padj), Gene) %>% filter(Gene %in%  PD_comparison_6h_unique) %>% head(5) %>% pull(Gene)), 
                                   (PD_4mix_6h %>% arrange((padj), Gene) %>% filter(Gene %in%  PD_comparison_6h_unique) %>% head(5) %>% pull(Gene)), 
                                   (PD_3mix_6h %>% arrange(desc(abs(log2FoldChange)), Gene) %>% filter(Gene %in%  PD_comparison_6h_unique) %>% head(5) %>% pull(Gene)), 
                                   (PD_4mix_6h %>% arrange(desc(abs(log2FoldChange)), Gene) %>% filter(Gene %in%  PD_comparison_6h_unique) %>% head(5) %>% pull(Gene)))


PD_6h_shared <-   get.venn.partitions(list(PD3mix6h=(PD_3mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene)), 
                                           PD4mix6h =(PD_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene))))

PD_venn <- venn.diagram(list(PD3mix6h=(PD_3mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene)), 
                             PD4mix6h =(PD_4mix_6h %>% filter(padj<0.05 & abs(log2FoldChange) >1) %>% pull(Gene))),
                        filename = NULL, output = T, 
                        #  height =2000 , 
                        #  width = 2000, 
                        #  resolution = 600,
                        lwd = 2,
                        lty = 'blank',
                        fill = c("black", "dodgerblue3"),
                        cex = 1,
                        cat.cex = 1,
                        cat.dist = c(0.05, 0.05),
                        cat.default.pos = 'outer',
                        cat.pos= c(-20, 20),
                        fontfamily = 'serif',
                        cat.fontfamily = 'serif',
                        cat.col = c("black", "dodgerblue3"), main= 'CB', sub= 'hyper'
)


ggsave(PD_venn, file="PD_venn.pdf", device = "pdf")


##### volcanoplot figure 4F

nb.cols <- 53
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")


Figure_4F_unique <- BP_alldata %>% filter(Comparison == "3mix_6h" |Comparison == "4mix_6h") %>% 
  mutate(Significant = case_when(
    padj < 0.05 & abs(log2FoldChange) >1 & !(Gene %in% BP_comparison_6h) ~ "sign_unique",
    padj > 0.05 | abs(log2FoldChange) <1  ~ "not_sign",
    padj < 0.05 & abs(log2FoldChange) > 1 & Gene %in% BP_comparison_6h ~ "sign_both")) %>%
  arrange(Significant!="not_sign", Significant!="sign_both") %>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) +
  #geom_point(size=4, shape=1, aes(colour=Significant)) +
  geom_point(size=3, shape=16, aes(colour=Significant), alpha = 1/2) +
  scale_colour_manual(values=c("sign_both"="deeppink3", "sign_unique"="#00B2FF","not_sign"="darkgrey")) +
  geom_point(data=pick(~Gene %in% BP_comparison_6h_unique_names_selected), shape=21, size=5, colour="black", aes(fill=as.factor(Gene_name))) +
  scale_fill_jco() +
  geom_hline(yintercept = -log10(0.05), colour="#990000", linetype="dashed") +
  geom_vline(xintercept = 1, colour="#990000", linetype="dashed") +   
  geom_vline(xintercept = -1, colour="#990000", linetype="dashed") +
  theme_bw(base_size = 16) + 
  xlim(-7,7) +
  guides(fill = guide_legend(ncol = 1, title="Gene") )  +
  facet_grid( . ~ Comparison) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank()) 


Figure_4F_common <- BP_alldata %>% filter(Comparison == "3mix_6h" |Comparison == "4mix_6h") %>% 
  mutate(Significant = case_when(
    padj < 0.05 & abs(log2FoldChange) >1 & !(Gene %in% BP_comparison_6h) ~ "sign_unique",
    padj > 0.05 | abs(log2FoldChange) <1  ~ "not_sign",
    padj < 0.05 & abs(log2FoldChange) > 1 & Gene %in% BP_comparison_6h ~ "sign_both")) %>%
  arrange(Significant!="not_sign", Significant!="sign_unique") %>%
  ggplot(aes(x=log2FoldChange, y=-log10(padj))) +
  #geom_point(size=4, shape=1, aes(colour=Significant)) +
  geom_point(size=3, shape=16, aes(colour=Significant), alpha = 1/2) +
  scale_colour_manual(values=c("sign_both"="deeppink3", "sign_unique"="#00B2FF","not_sign"="darkgrey")) +
  geom_point(data=pick(~Gene %in% BP_comparison_6h_names_selected),shape=21, size=5, colour="black", aes(fill=as.factor(Gene_name))) +
  scale_fill_jco() +
  geom_hline(yintercept = -log10(0.05), colour="#990000", linetype="dashed") +
  geom_vline(xintercept = 1, colour="#990000", linetype="dashed") +   
  geom_vline(xintercept = -1, colour="#990000", linetype="dashed") +
  theme_bw(base_size = 16) + 
  xlim(-7,7) +
  guides(fill = guide_legend(ncol = 1, title="Gene") )  +
  facet_grid( . ~ Comparison) +
  theme(strip.background = element_blank(),
        strip.text.y = element_blank()) 


##### Correlation Figure 4E ####

Figure_4E_BP <- BP_4mix_6h %>% full_join(BP_3mix_6h, by='Gene') %>% filter((padj.x <0.05 & abs(log2FoldChange.x) >1) |(padj.y <0.05 & abs(log2FoldChange.y) >1)) %>%
 mutate(padj.x = coalesce(padj.x, 1)) %>%
 mutate(padj.y = coalesce(padj.y, 1)) %>%
  mutate(Significant = case_when(
    (padj.x < 0.05 & abs(log2FoldChange.x) >1) & (padj.y > 0.05 | abs(log2FoldChange.y) <1) ~ "sign_4mix_only",
    (padj.x > 0.05 | abs(log2FoldChange.x) <1) & (padj.y < 0.05 | abs(log2FoldChange.y) >1) ~ "sign_3mix_only",
    (padj.x < 0.05 & abs(log2FoldChange.x) >1) & (padj.y < 0.05 | abs(log2FoldChange.y) >1) ~ "sign_both")) %>%
  arrange(Significant!="not_sign", Significant!="sign_4mix_only", Significant!="sign_3mix_only") %>%
  ggplot(aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
  geom_point(size=4, shape=16, aes(colour=Significant)) +
  scale_colour_manual(values=c("sign_both"="deeppink3", "sign_4mix_only"="#00B2FF","sign_3mix_only"="black", "not_sign"="darkgrey")) +
  theme_bw(base_size = 16) +
  xlab("log2FoldChange in 4mix") +
  ylab("log2FoldChange in 3mix") +
  ggtitle("Correlation BP DEG 4mix vs 3mix (6h)") +
  geom_hline(yintercept = -1, colour="#990000", linetype="dashed") +
  geom_hline(yintercept = +1, colour="#990000", linetype="dashed") +
  geom_vline(xintercept = 1, colour="#990000", linetype="dashed") +   
  geom_vline(xintercept = -1, colour="#990000", linetype="dashed") + 
  stat_cor(method="pearson", size=7) +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        aspect.ratio=1)

Figure_4E_CB <- CB_4mix_6h %>% full_join(CB_3mix_6h, by='Gene') %>% filter((padj.x <0.05 & abs(log2FoldChange.x) >1) |(padj.y <0.05 & abs(log2FoldChange.y) >1)) %>%
  mutate(padj.x = coalesce(padj.x, 1)) %>%
  mutate(padj.y = coalesce(padj.y, 1)) %>%
  mutate(Significant = case_when(
    (padj.x < 0.05 & abs(log2FoldChange.x) >1) & (padj.y > 0.05 | abs(log2FoldChange.y) <1) ~ "sign_4mix_only",
    (padj.x > 0.05 | abs(log2FoldChange.x) <1) & (padj.y < 0.05 | abs(log2FoldChange.y) >1) ~ "sign_3mix_only",
    (padj.x < 0.05 & abs(log2FoldChange.x) >1) & (padj.y < 0.05 | abs(log2FoldChange.y) >1) ~ "sign_both")) %>%
  arrange(Significant!="sign_4mix_only", Significant!="sign_3mix_only") %>%
  ggplot(aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
  geom_point(size=4, shape=16, aes(colour=Significant)) +
  scale_colour_manual(values=c("sign_both"="deeppink3", "sign_4mix_only"="#00B2FF","sign_3mix_only"="black")) +
  theme_bw(base_size = 16) +
  xlab("log2FoldChange in 4mix") +
  ylab("log2FoldChange in 3mix") +
  ggtitle("Correlation CB DEG 4mix vs 3mix (6h)") +
  geom_hline(yintercept = -1, colour="#990000", linetype="dashed") +
  geom_hline(yintercept = +1, colour="#990000", linetype="dashed") +
  geom_vline(xintercept = 1, colour="#990000", linetype="dashed") +   
  geom_vline(xintercept = -1, colour="#990000", linetype="dashed") + 
  stat_cor(method="pearson", size=7) +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        aspect.ratio=1)


Figure_4E_PD <- PD_4mix_6h %>% full_join(PD_3mix_6h, by='Gene') %>% filter((padj.x <0.05 & abs(log2FoldChange.x) >1) |(padj.y <0.05 & abs(log2FoldChange.y) >1)) %>%
  mutate(padj.x = coalesce(padj.x, 1)) %>%
  mutate(padj.y = coalesce(padj.y, 1)) %>%
  mutate(Significant = case_when(
    (padj.x < 0.05 & abs(log2FoldChange.x) >1) & (padj.y > 0.05 | abs(log2FoldChange.y) <1) ~ "sign_4mix_only",
    (padj.x > 0.05 | abs(log2FoldChange.x) <1) & (padj.y < 0.05 | abs(log2FoldChange.y) >1) ~ "sign_3mix_only",
    (padj.x < 0.05 & abs(log2FoldChange.x) >1) & (padj.y < 0.05 | abs(log2FoldChange.y) >1) ~ "sign_both")) %>%
  arrange(Significant!="snot_sign", Significant!="sign_4mix_only", Significant!="sign_3mix_only") %>%
  ggplot(aes(x=log2FoldChange.x, y=log2FoldChange.y)) +
  geom_point(size=4, shape=16, aes(colour=Significant)) +
  scale_colour_manual(values=c("sign_both"="deeppink3", "sign_4mix_only"="#00B2FF","sign_3mix_only"="black")) +
  theme_bw(base_size = 16) +
  xlab("log2FoldChange in 4mix") +
  ylab("log2FoldChange in 3mix") +
  ggtitle("Correlation PD DEG 4mix vs 3mix (6h)") +
  geom_hline(yintercept = -1, colour="#990000", linetype="dashed") +
  geom_hline(yintercept = +1, colour="#990000", linetype="dashed") +
  geom_vline(xintercept = 1, colour="#990000", linetype="dashed") +   
  geom_vline(xintercept = -1, colour="#990000", linetype="dashed") + 
  stat_cor(method="pearson", size=7) +
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        aspect.ratio=1)


BP_alldata %>% filter(Comparison == '4mix_6h') %>% write.table("SB172_BP_4mix_6h_DEG_allsamples.txt")
BS_alldata %>% filter(Comparison == '4mix_6h') %>% write.table("SB172_BS_4mix_6h_DEG_allsamples.txt")
CB_alldata %>% filter(Comparison == '4mix_6h') %>% write.table("SB172_CB_4mix_6h_DEG_allsamples.txt")
PD_alldata %>% filter(Comparison == '4mix_6h') %>% write.table("SB172_PD_4mix_6h_DEG_allsamples.txt")


##### Transcriptional activity levels Supplementary Figure 4A #####
## reads were obtained as follows from featurecounts run

#BP_reads <- data.frame('Reads' = colSums(BP_FtCounts$stat[,-1])) %>%
#  rownames_to_column(var='Sample') %>%
#  mutate(Organism = 'BP')
#BS_reads <- data.frame('Reads' = colSums(BS_FtCounts$stat[,-1])) %>%
#  rownames_to_column(var='Sample') %>%
#  mutate(Organism = 'BS')
#CB_reads <- data.frame('Reads' = colSums(CB_FtCounts$stat[,-1])) %>%
#  rownames_to_column(var='Sample') %>%
#  mutate(Organism = 'CB') 
#PD_reads <- data.frame('Reads' = colSums(PD_FtCounts$stat[,-1])) %>%
#  rownames_to_column(var='Sample') %>%
#  mutate(Organism = 'PD')   

#All_reads <- bind_rows(BP_reads, BS_reads, CB_reads, PD_reads) %>%
#  mutate(Sample = gsub(".BP.only_aligned.bam|.BS.only_aligned.bam|.CB.only_aligned.bam|.PD.only_aligned.bam", "", Sample)) %>%
#  mutate(Sample = gsub('SB172.', '', Sample))
#All_reads$Sample <-  str_replace(All_reads$Sample, '37', '3mix_0h_1') %>% str_replace('38', '3mix_0h_2') %>% str_replace('39', '3mix_0h_3') %>%
#  str_replace('40', '3mix_6h_4') %>% str_replace('41', '3mix_6h_5') %>% str_replace('42', '3mix_6h_6') %>% str_replace('43', '4mix_6h_7') %>%
#  str_replace('44', '4mix_6h_8') %>% str_replace('45', '4mix_6h_9') %>% str_replace('46', '4mix_0h_10') %>% str_replace('47', '4mix_0h_11') %>%
#  str_replace('48', '4mix_0h_12')
#All_reads <- All_reads %>% mutate(Group=str_extract(Sample, pattern="3mix|4mix"))                             
                              
### read-in the above generated txt

All_reads <- read.delim2('SB172_All_reads.txt', sep = '\t')

df <- read.delim('DADA2_SB093_SB172_table.txt')

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

SB172_DNA_df <- SB172 %>%
  filter(grepl('3mix.post|4mix.post|SB172', sample)) %>% 
  mutate(sample=gsub('\\..pool957', '', sample)) %>%
  mutate(sample=gsub('3mix.post.', '', sample)) %>%
  mutate(sample=gsub('4mix.post.', '', sample)) %>%
  mutate(sample=gsub('SB172.INPUT3MIX', '3mix_input', sample)) %>%
  mutate(sample=gsub('SB172.INPUT4MIX', '4mix_input', sample)) %>%
  mutate(Organism = case_when(
    Species == 'Bacteroides sartorii' ~ 'BS',
    Species == 'Blautia producta' ~ 'BP',
    Species == 'Clostridium bolteae' ~ 'CB',
    Species == 'Parabacteroides distasonis' ~ 'PD',
    Species == 'Other' ~ 'Other'
  )) %>% select(-Species) %>%
  filter(!str_detect(sample, 'input')) %>%
  rename(Sample=sample) %>%
  mutate(Sample = gsub('(.*)cecum.', '', Sample)) %>%
  mutate(Group = gsub('-','', Group))

All_reads_for_ratio <- All_reads %>% 
  mutate(Time_point = str_extract(Sample, '0h|6h')) %>%
  mutate(Sample = gsub("(.*)_",  '', Sample)) %>%
  group_by(Time_point, Group, Sample) %>%
  mutate(RNA_pctseq=Reads/sum(Reads)) %>%
  ungroup() %>%
  left_join(SB172_DNA_df, by=c('Sample', 'Organism', 'Group')) %>%
  filter(!is.na(pctseqs))

### according to PATRIC, PD and BS have 7 copies of 16s, BP and CB have 5, 
### take those into account by dividing BS and PD times 7 and BP CB times 5
### then sum/normalize again
All_reads_for_ratio <- All_reads_for_ratio %>% 
  mutate(pctseqs_norm = case_when(
    Organism == 'BS' ~ pctseqs/7,
    Organism == 'PD' ~ pctseqs/7,
    Organism == 'CB'  ~ pctseqs/5,
    Organism == 'BP' ~ pctseqs/5)) %>%
  group_by(Sample) %>%
  mutate(pctseqs_norm=pctseqs_norm/sum(pctseqs_norm)) %>%
  ungroup() %>%
  mutate(RNA_DNA_ratio=RNA_pctseq/pctseqs_norm)

Supplementary_Figure_4A <- All_reads_for_ratio %>%
  ggplot(aes(fill=factor(Organism), 
             y=RNA_DNA_ratio, 
             x=factor(Time_point, levels=c('0h', '6h'))))  + 
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
       title = 'SB172 RNA/DNA ratio') +
  facet_wrap(factor(Organism, levels=c('CB', 'BP', 'BS', 'PD'))~factor(Group, levels=c('3mix', '4mix')), ncol=2 , scales = "fixed", strip.position = 'right') +
  expand_limits(y=0) 






