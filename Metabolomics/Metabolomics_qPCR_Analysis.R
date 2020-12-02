# Metabolomics Analysis Becattini S. et al Cell Host & Microbe

# 0. Prep: Load Packages -----------------------------------------------------------
library(tidyr)
library(dplyr)
library(cowplot)
library(ggplot2)
library(viridis)
library(reshape2)
library(scales)
library("RColorBrewer")

setwd('/Users/becattis/Desktop/CHM revision/Scripts to be uploaded/Matt/')
# 0. Prep: Functions and Colors ----------------------------------------------------


'%!in%' <- function(x,y)!('%in%'(x,y))

# Load Heatmap colors
  bin_values <- read.csv(file="bin_values.csv",header=TRUE,sep=",")
  bin_values$bin_value <- as.character(bin_values$bin_value)
  bin_values$color <- as.character(bin_values$color)
  bin.colors <- bin_values$color
  names(bin.colors) <- bin_values$bin_value

# 0. Prep: Load Compounds/Metadata -------------------------------------------------

# Load compound information
  master_compound_list<- read.csv(file="master_compound_list.csv",header=TRUE,sep=",")%>%
    select(-X)

# Load metadata for mouse experiments
  mouse_meta_data <- read.csv(file="mouse_meta_data.csv",header=TRUE,sep=",")%>%
    select(-X)

# Load metadata for  in vitro CBBP growth 

  in_vitro_meta_data <- read.csv(file="in_vitro_meta_data.csv",header=TRUE,sep=",")%>%
    select(-X)

# 0. Prep: Load GCMS Data --------------------------------------------------

# Load metabolite values from in vivo SPF and CBBP experiments 
  # then join with compound information
  
  mouse_master_clean <- read.csv(file="murine_metabolomics_data.csv",header=TRUE,sep=",")%>%
    select(-X)%>%
    inner_join(master_compound_list)

# Load metabolite values from in vitro CBBP experiments 
  # then join with compound information
  in_vitro_metabolites <- read.csv(file="in_vitro_metabolites.csv") %>%
    select(-X)%>%
    inner_join(master_compound_list)

# 1. In vivo experiments: Normalization to Internal Standards -------------------------------------

mouse_master_clean <- mouse_master_clean %>%
  dplyr::mutate(tube_label=as.numeric(gsub("c","",sample)))

# Generate normalization factors for each sample/injection based on internal standards 
  mouse_normalization_factors <- mouse_master_clean %>%
    inner_join(mouse_meta_data,by=c("run"="run","tube_label"="tube_label"))%>%
    filter(type=="ISTD") %>%
    dplyr::group_by(run,dilution,compound) %>%
    dplyr::mutate(compound_average=mean(value)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(norm_compound = value/compound_average) %>%
    dplyr::mutate(label=paste(run,tube_label,dilution)) %>%
    dplyr::group_by(run,dilution,tube_label)%>%
    dplyr::summarise(normalization_factor=mean(norm_compound))


# 2. Antibiotic Naive (SPF) experiment data processing--------------------------------------------

# Calculate average normalized peak area for each compound in control mouse group
  antibiotic_naive_controls <- mouse_master_clean %>%
    inner_join(mouse_normalization_factors,by=c("run"="run","tube_label"="tube_label","dilution"="dilution"))%>%
    filter(type!="ISTD")%>%
    inner_join(mouse_meta_data,by=c("run"="run","tube_label"="tube_label"))%>%
    dplyr::mutate(norm_value=value/normalization_factor)%>%
    filter(experiment=="SB193")%>%
    filter(Treament=="none")%>%
    dplyr::group_by(compound)%>%
    dplyr::summarise(control_average=mean(norm_value))

# Calculate fold change relative to average control mouse group value
  antibiotic_naive_df <- mouse_master_clean %>%
    inner_join(mouse_normalization_factors,by=c("run"="run","tube_label"="tube_label","dilution"="dilution"))%>%
    filter(type!="ISTD")%>%
    inner_join(mouse_meta_data,by=c("run"="run","tube_label"="tube_label"))%>%
    dplyr::mutate(norm_value=value/normalization_factor)%>%
    filter(experiment=="SB193")%>%
    inner_join(antibiotic_naive_controls,by="compound") %>%
    dplyr::mutate(fold_change = norm_value/control_average)%>%
    mutate(bin_value=cut(fold_change,breaks=c(-1,0.2,.5,(1/1.3),1.3,2,5,10,30,500),
                         labels=c("deplete >5x",
                                  "deplete 2-5x",
                                  "deplete 1.3-2x",
                                  "1.3x deplete-1.3x increase",
                                  "increase 1.3-2x","increase 2-5x",
                                  "increase 5-10x","increase 10-30x",
                                  ">30X")))

#Organize factor levels, clean up compound names 
  antibiotic_naive_df$Treament <- factor(antibiotic_naive_df$Treament,levels=c("none","aCD3"))
  antibiotic_naive_df$compound <- gsub("X","",antibiotic_naive_df$compound)
  antibiotic_naive_df$compound <- gsub("\\.","-",antibiotic_naive_df$compound)
  antibiotic_naive_df$type <- gsub("_"," ",antibiotic_naive_df$type)
  antibiotic_naive_df$type <- factor(antibiotic_naive_df$type,levels=c("SCFA","Branch SCFA","Aminated SCFA","Hydroxylated SCFA","Amino acid","Aromatic","Other"))
  antibiotic_naive_df$compound <- reorder(antibiotic_naive_df$compound,antibiotic_naive_df$molecular_weight)



# 3. Plots of Antibiotic Naive (SPF) Mice ---------------------------------

  
aa_plot <- antibiotic_naive_df %>%
  filter(type=="Amino acid")%>%
  ggplot(aes(x=Treament,y=fold_change,fill=Treament))+
    facet_wrap(~compound,ncol=4,scales="free_y")+
    geom_bar(stat="summary",fun=mean,alpha=.4,color="black")+
    geom_jitter(pch=21,size=4,color="black",width=.2)+
    scale_fill_viridis(discrete = TRUE)+
    ylab(label="Change relative to Average Untreated Mouse")+
    xlab(label="Treatment")+
    theme_minimal()+
    theme(legend.position = "none",
          panel.border = element_rect(fill=NA,color="black",size=.5))
aa_plot



Supplemental_Fig_6c <- ggplot(antibiotic_naive_df,aes(x=compound,y=as.factor(tube_label),fill=bin_value))+
  geom_tile(color="white",size=.15)+
  scale_fill_manual(values=bin.colors)+
  facet_grid(factor(Treament, levels = c('aCD3','none'))~ type, scales="free",space="free",switch = "y")+
  labs(fill="Fold Change to Average \nUntreated Mouse")+
  theme_bw(base_size = 16) + 
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  

pdf('Supplemental_Fig_6c.pdf', width = 10, height = 4)
Supplemental_Fig_6c
dev.off()


# 4. Antibiotic treated (CBBP) experiment data processing --------------------------------------------------

# Calculate average normalized peak area for each compound in control mouse group
  cbbp_controls <- mouse_master_clean %>%
    inner_join(mouse_normalization_factors,by=c("run"="run","tube_label"="tube_label","dilution"="dilution"))%>%
    filter(type!="ISTD")%>%
    inner_join(mouse_meta_data,by=c("run"="run","tube_label"="tube_label"))%>%
    dplyr::mutate(norm_value=value/normalization_factor)%>%
    filter(Group=="4mix")%>%
    filter(Treament=="none")%>%
    filter(group!="antibiotic_naive")%>%
    dplyr::group_by(experiment,compound)%>%
    dplyr::summarise(control_average=mean(norm_value))

# Calculate fold change relative to average control mouse group value
cbbp_df <- mouse_master_clean %>%
  inner_join(mouse_normalization_factors,by=c("run"="run","tube_label"="tube_label","dilution"="dilution"))%>%
  filter(type!="ISTD")%>%
  inner_join(mouse_meta_data,by=c("run"="run","tube_label"="tube_label"))%>%
  dplyr::mutate(norm_value=value/normalization_factor)%>%
  filter(Group=="4mix")%>%
  filter(Mouse_type=="Antibiotic")%>%
  inner_join(cbbp_controls,by=c("compound"="compound","experiment"="experiment")) %>%
  dplyr::mutate(fold_change = norm_value/control_average) 

# Remove compounds not robustly detected in these mice
  cbbp_df$fold_change <- ifelse(cbbp_df$group=="C",NA,cbbp_df$fold_change)

# Cut values into bins for heatmap
  cbbp_df <- cbbp_df %>%
    mutate(bin_value=cut(fold_change,breaks=c(-1,0.2,.5,(1/1.3),1.3,2,5,10,30,max(mouse_master_clean$value)),
                         labels=c("deplete >5x",
                                  "deplete 2-5x",
                                  "deplete 1.3-2x",
                                  "1.3x deplete-1.3x increase",
                                  "increase 1.3-2x",
                                  "increase 2-5x",
                                  "increase 5-10x",
                                  "increase 10-30x",
                                  ">30X increase")))


# Clean up compound names, reorder compounds by molecular weight, order the compound types
  cbbp_df$plot_order <- paste(as.character(cbbp_df$Treament),as.character(cbbp_df$Time))
  cbbp_df$plot_order <- factor(cbbp_df$plot_order,levels=c("none 0","aCD3 6","aCD3 24","flagellin 6","flagellin 24"))
  cbbp_df$compound <- gsub("X","",cbbp_df$compound)
  cbbp_df$compound <- gsub("\\.","-",cbbp_df$compound)
  cbbp_df$type <- gsub("_"," ",cbbp_df$type)
  cbbp_df$type <- factor(cbbp_df$type,levels=c("SCFA","Branch SCFA","Aminated SCFA","Hydroxylated SCFA","Amino acid","Aromatic","Other"))
  cbbp_df$compound <- reorder(cbbp_df$compound,cbbp_df$molecular_weight)

cbbp_df_CD3 <- cbbp_df %>%
  filter(Treament!="flagellin")
cbbp_df_flagellin <- cbbp_df %>%
  filter(Treament!="aCD3")%>%
  filter(experiment=="SB093")


# 5. Antibiotic treated (CBBP) experiment plots ---------------------------


Supplemental_Fig_6a<- ggplot(cbbp_df_CD3,aes(x=compound,y=as.factor(tube_label),fill=bin_value))+
  geom_tile(color="white",size=.15)+
  scale_fill_manual(values=bin.colors,na.value="#D6D6D6")+
  labs(fill="Fold Change to Average\nUntreated Mouse")+
  facet_grid(factor(plot_order, levels = c('aCD3 24','aCD3 6','none 0'))~type,
             scales="free",space="free")+
  #scale_x_discrete(breaks=figure_2A_heat_df$properorder,labels=figure_2A_heat_df$species)+
  theme_bw(base_size = 16) + 
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf('Supplemental_Fig_6a.pdf', width = 10, height = 7)
Supplemental_Fig_6a
dev.off()

Supplemental_Fig_6b<- ggplot(cbbp_df_flagellin,aes(x=compound,y=as.factor(tube_label),fill=bin_value))+
  geom_tile(color="white",size=.15)+
  scale_fill_manual(values=bin.colors,na.value="#D6D6D6")+
  labs(fill="Fold Change to Average\nUntreated Mouse")+
  facet_grid(factor(plot_order, levels = c('flagellin 24','flagellin 6','none 0'))~type,
             scales="free",space="free")+
  theme_bw(base_size = 16) + 
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf('Supplemental_Fig_6b.pdf', width = 10, height = 4)
Supplemental_Fig_6b
dev.off()


# 6. In vitro experiments: Normalization to Internal Standards----------------------------------------------------------------

invitro_normalization_factors <- in_vitro_metabolites %>%
  inner_join(in_vitro_meta_data,by=c("experiment"="experiment","tube_label"="tube_label"))%>%
  filter(type=="ISTD") %>%
  dplyr::group_by(experiment,dilution,compound) %>%
  dplyr::mutate(compound_average=mean(value)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(norm_compound = value/compound_average) %>%
  dplyr::mutate(label=paste(experiment,tube_label,dilution)) %>%
  dplyr::group_by(experiment,dilution,tube_label)%>%
  dplyr::summarise(normalization_factor=mean(norm_compound))


# 7. In vitro experiment data processing ---------------------------------

# Calculate average normalized peak area for each compound in BHI media control
  invitro_media_averages <- in_vitro_metabolites %>%
    inner_join(invitro_normalization_factors,by=c("experiment"="experiment",
                                             "dilution"="dilution",
                                             "tube_label"="tube_label")) %>%
    inner_join(in_vitro_meta_data,by=c("experiment"="experiment","tube_label"="tube_label"))%>%
    dplyr::mutate(final_value=value/normalization_factor) %>%
    filter(type!="ISTD")%>%
    filter(class=="Media")%>%
    dplyr::group_by(experiment,compound,dilution)%>%
    dplyr::summarise(media_average=mean(final_value))    

# Calculate fold change relative to average control mouse group value
  invitro_df <- in_vitro_metabolites %>%
    inner_join(invitro_normalization_factors,by=c("experiment"="experiment",
                                             "dilution"="dilution",
                                             "tube_label"="tube_label")) %>%
    inner_join(in_vitro_meta_data,by=c("experiment"="experiment","tube_label"="tube_label"))%>%
    inner_join(invitro_media_averages,by=c("compound"="compound","dilution"="dilution","experiment"="experiment"))%>%
    dplyr::mutate(final_value=value/normalization_factor) %>%
    dplyr::mutate(fold_change=final_value/media_average)%>%
    filter(class=="Isolate")

# Filter out compounds that are not detected in BHI
  
  invitro_df$fold_change <- ifelse(invitro_df$group=="B",NA,invitro_df$fold_change)

# Cut values into bins for heatmap
  invitro_df<- invitro_df %>%  
    mutate(bin_value=cut(fold_change,breaks=c(-1,.25,.75,1.25,2.5,5.0,7.5,25,50,500),
                         labels=c("0-0.25",".25-.75",
                                  ".75-1.25","1.25-2.5",
                                  "2.5-5.0","5.0-7.5",
                                  "7.5-25","25-50",
                                  "50-500")))

# Label phylum, clean up compound names, order factors
  invitro_df <- invitro_df %>%
    dplyr::mutate(phylum=ifelse(msk_id=="Bacteroides sartorii" | msk_id=="Parabacteroides distasonis",
                                "Bacteroidetes","Firmicutes"))
  invitro_df$compound <- gsub("X","",invitro_df$compound)
  invitro_df$compound <- gsub("\\.","-",invitro_df$compound)
  invitro_df$type <- gsub("_"," ",invitro_df$type)
  invitro_df$type <- factor(invitro_df$type,levels=c("SCFA","Branch SCFA","Aminated SCFA",
                                                     "Hydroxylated SCFA","Amino acid",
                                                     "Aromatic","Other"))
  invitro_df$compound <- reorder(invitro_df$compound,desc(invitro_df$molecular_weight))
  invitro_df$msk_id <- factor(invitro_df$msk_id,levels=c("Parabacteroides distasonis",
                                                         "Bacteroides sartorii",
                                                         "Blautia producta",
                                                         "[Clostridium] bolteae"))


# 8. In vitro CBBP Plots --------------------------------------------------

invitro_df$compound <- reorder(invitro_df$compound,invitro_df$molecular_weight)

Figure_6A <- invitro_df %>%
  filter(time %in% c("48.1","48.2"))%>%
  dplyr::group_by(compound,dilution,type,molecular_weight,class,msk_id,phylum)%>%
  dplyr::summarise(fold_change=mean(fold_change))%>%
  mutate(bin_value=cut(fold_change,breaks=c(-1,.25,.75,1.25,2.5,5.0,7.5,25,50,500),
                       labels=c("0-0.25",".25-.75",".75-1.25","1.25-2.5","2.5-5.0","5.0-7.5", "7.5-25","25-50","50-500")))%>%
  ggplot(aes(x=compound,y=msk_id,fill=bin_value))+
  geom_tile(color="white",size=.15)+
  # scale_fill_brewer()
  scale_fill_manual(values = bin.colors,na.value="#D6D6D6")+
  # scale_fill_viridis(discrete = TRUE)+
  facet_grid(phylum~type,scales="free",space="free",switch = "y")+
  labs(fill="Fold Change to \nMedia")+
  theme_bw(base_size = 16) + 
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf('Figure_6A.pdf', width=10, height = 3.5)
Figure_6A
dev.off()

# 9. Bar plots


Figure_6B <- cbbp_df %>%
  filter(Treament!="flagellin") %>% 
  arrange(compound) %>%
  ggplot(aes(x=plot_order,y=fold_change,fill=plot_order))+ 
  facet_wrap(~factor(compound),ncol=13,scales="free_y")+   #type~compound if you want them grouped by class of molecule
  geom_bar(stat="summary",fun.y="mean",alpha=.4,color="black")+
  geom_jitter(pch=21,size=2,color="black",width=.2)+
  scale_fill_jco()+
  ylab(label="Change relative to Average Untreated Mouse")+
  theme_minimal(base_size = 14)+
  theme(legend.position = "none",
        #     strip.text.x = element_blank(),
        panel.border = element_rect(fill=NA,color="black",size=.5)) +    
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "none 0")

Figure_6B


Figure_6C <- antibiotic_naive_df %>%
  ggplot(aes(x=Treament,y=fold_change,fill=Treament))+
  facet_wrap(.~compound,ncol=13,scales="free_y")+
  geom_bar(stat="summary",fun.y="mean",alpha=.4,color="black")+
  geom_jitter(pch=21,size=2,color="black",width=.2)+
  scale_fill_jco()+
  ylab(label="Change relative to Average Untreated Mouse")+
  theme_minimal(base_size = 14)+
  theme(legend.position = "none",
        #     strip.text.x = element_blank(),
        panel.border = element_rect(fill=NA,color="black",size=.5)) +    
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = "none")

Figure_6C 
