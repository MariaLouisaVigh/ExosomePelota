#### load packages ####

library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)
library(ggrepel)
library(rtracklayer) 
library(systemPipeR)
library(eulerr)
library(GeneOverlap)
library(svglite)
library(reshape2)
library(wesanderson)
library(VennDiagram)
library(readxl)
library(ComplexUpset)
library(openxlsx)

#### miRNA target list ####
###### miRNA target list ######## 

miRNAtargets_rawtable <- readxl::read_xlsx("/binf-isilon/PBgrp/xpj980/TAIR/miRNAtargetlist_Maria_degradome_done.xlsx", col_names = T)


miRNAtargets <- miRNAtargets_rawtable %>% 
  select(miRNA, 
         gene = target_id, 
         target_name)

miRNAtargets <- filter(miRNAtargets, !duplicated(miRNAtargets$gene))


### extract all results from the same library

setwd("/binf-isilon/PBgrp/xpj980/datascratch/2019_exosome_mutants/04.STAR_mapping/sam_files/")

##import my featureCounts table 
raw_featureCount_table <- read.table("featureCountsaraport", header = TRUE, row.names = 1)
glimpse(raw_featureCount_table)

DESeq_table <- raw_featureCount_table %>%
  rownames_to_column() %>%
  select(1,7:15,19:28) %>% 
  column_to_rownames(var = "rowname") %>%
  glimpse()

### with araport11 gff, the exons need to be collapsed 
test_table <- DESeq_table %>% 
  rownames_to_column(var = "gene")  %>%
  as_tibble() %>% 
  separate(col = gene, into = c("gene", "exon", "exonno"), sep = ":") %>%
  select(- c("exon", "exonno")) %>% 
  glimpse()

test_table_summed <- test_table %>%
  group_by(gene) %>%
  summarise(across(everything(), sum)) %>% 
  as.data.frame()

DESeq_table_summed <- data.frame(test_table_summed, row.names = 1) %>% 
  glimpse()

names <- c("wt_R1", "wt_R2", "wt_R3", "GKpelota_R1", "GKpelota_R2", "GKpelota_R3", "hen2.5_R1", "hen2.5_R2", "hen2.5_R3", "rrp45a_R1", "rrp45a_R2", "rrp45a_R3", "rrp4_R1", "rrp4_R2", "rrp4_R3", "SAILpelota_R1", "ski2.5_R1", "ski2.5_R2", "ski2.5_R3")
names(DESeq_table_summed) <- names
glimpse(DESeq_table_summed)

## make metadata matrix
introduce_genotype_to_table  <- data.frame(row.names=colnames(DESeq_table_summed),
                                           "genotype"=c(rep("wt", 3),
                                                        rep("GKpelota", 3), 
                                                        rep("hen2.5", 3),
                                                        rep("rrp45a", 3),
                                                        rep("rrp4", 3),
                                                        rep("SAILpelota", 1),
                                                        rep("ski2.5", 3)))

str(introduce_genotype_to_table)
levels(introduce_genotype_to_table)
introduce_genotype_to_table$genotype  <- as.factor(introduce_genotype_to_table$genotype)
introduce_genotype_to_table$genotype <- relevel(introduce_genotype_to_table$genotype, ref="wt")
levels(introduce_genotype_to_table$genotype)

design  <- model.matrix(~genotype, data=introduce_genotype_to_table)

## convert featurecount table into a matrix
DESeq_matrix <- as.matrix(DESeq_table_summed)
head(DESeq_matrix)

## check that sample order is the same in my metadata matrix and in my featurecount matrix
all(rownames(design) == colnames(DESeq_matrix))

## make the DESeq2 object 
dds <- DESeqDataSetFromMatrix(countData=DESeq_matrix, colData=introduce_genotype_to_table, design = ~ genotype)
head(dds)

##########QA before DESeq analysis ########
## normalize counts & extract them
dds2 <- estimateSizeFactors(dds)
sizeFactors(dds2)
normalized_counts <- counts(dds2, normalized = TRUE)

## log_transformation with vst function 
log_transformeddds2 <- vst(dds2, blind = TRUE)

##hierachial heatmap
vsd_log_transformed <- assay(log_transformeddds2)
pairwise_cor <- cor(vsd_log_transformed)

hierarchial_heat2 <- pheatmap(pairwise_cor)
#ggsave(hierarchial_heat2, file="/binf-isilon/PBgrp/xpj980/figures/20210319hierachialheat2019Set1araport.svg", units = "in", width = 15, height = 15, dpi = 300)

#### Sup. Fig1 #####
####### PCA plot #######
PCAplot_genotype <- plotPCA(log_transformeddds2, intgroup = "genotype", returnData = T)

percentVar <- round(100 * attr(PCAplot_genotype, "percentVar"))

to_plot <- separate(PCAplot_genotype, name, into = c("genotype", "replicate"), sep = "_")

PCAplot <- ggplot(to_plot, aes(x=PC1, y=PC2, color=factor(genotype), label=genotype, shape=replicate)) +
  geom_point(size=3)+
  labs(color="genotype")+
  geom_text_repel(size=4, fontface = 'bold', show.legend = FALSE) +
  #geom_text(aes(label=sample), size=2, hjust=0.9, vjust=1.5)+
  color_palette(palette = "Paired") +
  ggtitle('PC1 and PC2') +
  theme(plot.title = element_text(hjust = 0.5)) +
  cowplot::theme_cowplot() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(aspect.ratio = 1) +
  coord_fixed()

#ggsave(PCAplot, file = "/binf-isilon/PBgrp/xpj980/figures/20210319PCA2019dataset1araport.svg", units = "in", width = 15, height = 15, dpi = 300)

sampleDists <- dist(t(assay(log_transformeddds2)))

sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
test <- pheatmap(sampleDistMatrix,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists,
                 col=colors)

ggsave(test, file = "/binf-isilon/PBgrp/xpj980/figures/20210324distancematdataset2araport.svg", units = "in", width = 15, height = 15, dpi = 300)


########## DESeq Analysis #########
deseqqed <- DESeq(dds2)

resultsNames(deseqqed)

## extract results 
DESeqresultsGKpelota <- DESeq2::results(deseqqed, contrast = c("genotype", "GKpelota", "wt"), alpha = 0.1)
DESeqresultsSAILpelota <- DESeq2::results(deseqqed, contrast = c("genotype", "SAILpelota", "wt"), alpha = 0.1)
DESeqresultsSki2.5 <- DESeq2::results(deseqqed, contrast = c("genotype", "ski2.5", "wt"), alpha = 0.1)
DESeqresultsRrp45a <- DESeq2::results(deseqqed, contrast = c("genotype", "rrp45a", "wt"), alpha = 0.1)
DESeqresultsRrp4 <- DESeq2::results(deseqqed, contrast = c("genotype", "rrp4", "wt"), alpha = 0.1)
DESeqresultsHen2.5 <- DESeq2::results(deseqqed, contrast = c("genotype", "hen2.5", "wt"), alpha = 0.1)

###### gathered DESegresults-table #######
GKpelota_tibble <- DESeqresultsGKpelota %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  mutate(genotype = rep("GKpelota")) %>%
  glimpse()

SAILpelota_tibble <- DESeqresultsSAILpelota %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  mutate(genotype = rep("SAILpelota")) %>%
  glimpse()

ski2.5_tibble <- DESeqresultsSki2.5 %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  mutate(genotype = rep("ski2.5")) %>%
  glimpse()

rrp45a_tibble <- DESeqresultsRrp45a %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  mutate(genotype = rep("rrp45a")) %>%
  glimpse()

rrp4_tibble <- DESeqresultsRrp4 %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  mutate(genotype = rep("rrp4")) %>%
  glimpse()

hen2.5_tibble <- DESeqresultsHen2.5 %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  mutate(genotype = rep("hen2.5")) %>%
  glimpse()

gathered_DESeqresults_rrp4 <- bind_rows(... = rrp4_tibble, ski2.5_tibble, hen2.5_tibble, rrp45a_tibble) %>%
  glimpse()
gathered_DESeqresults_pelota <- bind_rows(... = GKpelota_tibble, SAILpelota_tibble, ski2.5_tibble) %>%
  glimpse()
gathered_DESeqresults_rrp45a <- bind_rows(... = ski2.5_tibble, hen2.5_tibble, rrp45a_tibble) %>%
  glimpse()
gathered_DESeqresults <- bind_rows(... = GKpelota_tibble, SAILpelota_tibble, ski2.5_tibble, rrp45a_tibble, rrp4_tibble,hen2.5_tibble) %>% 
  glimpse()

write.xlsx(x = gathered_DESeqresults, file = "/binf-isilon/PBgrp/xpj980/datascratch/GEO/deseqcontrast2.xlsx")
#save(gathered_DESeqresults_rrp4, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/2019exoDESEQrrp4.RData")
#save(gathered_DESeqresults_pelota, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/2019exoDESEQpelota.RData")
#save(gathered_DESeqresults_rrp45a, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/2019exoDESEQrrp45a.RData")
#save(gathered_DESeqresults, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/2019exoDESEQtable.RData")
#### load DESEQ results ######
load("/binf-isilon/PBgrp/xpj980/R_scripts_collected/2019exoDESEQrrp4.RData")
load("/binf-isilon/PBgrp/xpj980/R_scripts_collected/2019exoDESEQpelota.RData")
load("/binf-isilon/PBgrp/xpj980/R_scripts_collected/2019exoDESEQrrp45a.RData")
load("/binf-isilon/PBgrp/xpj980/R_scripts_collected/2019exoDESEQtable.RData")

## look for deseq miR163 targets for light induction experiment 
gathered_DESeqresults %>% filter(gene %in% c("AT1G66690", "AT1G66700", "AT1G66720", "AT3G44870", "AT3G44860", "AT1G15125")) %>% 
  filter(padj < 0.05)

#### excel sheets for paper ####
test <- left_join(gathered_DESeqresults, to_join, by = "gene") %>% 
  mutate(locus_type = case_when(is.na(locus_type) ~ "intergenic", 
                                T ~ locus_type)) %>% 
  filter(padj < 0.05) %>% 
  mutate(genotype = case_when(genotype == "GKpelota" ~ "pel1.2", 
                              T ~ genotype), 
         genotype = case_when(genotype == "SAILpelota" ~ "pel1.1", 
                              T ~ genotype)) %>% 
  filter(locus_type !="intergenic") %>% 
  left_join(., miRNAtargets, by = "gene", copy = FALSE) %>% 
  mutate(locus_type = case_when(is.na(miRNA) ~ locus_type, 
                                TRUE ~ "miRNA_target")) %>% 
  select(gene, baseMean, log2FoldChange, pvalue, padj, locus_type, genotype) 

#write.xlsx(x = test, file = "/binf-isilon/PBgrp/xpj980/figures/DESeq2ExperimentBWOintergenic.xlsx")

# filter pel1.1 on log2CF instead of padj as it only consists of one replicate (Revision)
part1 <- left_join(gathered_DESeqresults, to_join, by = "gene") %>% 
  mutate(locus_type = case_when(is.na(locus_type) ~ "intergenic", 
                                T ~ locus_type)) %>% 
  filter(padj < 0.05) %>% 
  mutate(genotype = case_when(genotype == "GKpelota" ~ "pel1.2", 
                              T ~ genotype), 
         genotype = case_when(genotype == "SAILpelota" ~ "pel1.1", 
                              T ~ genotype)) %>% 
  filter(locus_type !="intergenic") %>% 
  left_join(., miRNAtargets, by = "gene", copy = FALSE) %>% 
  mutate(locus_type = case_when(is.na(miRNA) ~ locus_type, 
                                TRUE ~ "miRNA_target")) %>% 
  select(gene, baseMean, log2FoldChange, pvalue, padj, locus_type, genotype) %>%
  filter(genotype != "pel1.1")

part2 <- left_join(gathered_DESeqresults, to_join, by = "gene") %>% 
  mutate(locus_type = case_when(is.na(locus_type) ~ "intergenic", 
                                T ~ locus_type)) %>% 
  filter(log2FoldChange > 0.9) %>% 
  mutate(genotype = case_when(genotype == "GKpelota" ~ "pel1.2", 
                              T ~ genotype), 
         genotype = case_when(genotype == "SAILpelota" ~ "pel1.1", 
                              T ~ genotype)) %>% 
  filter(locus_type !="intergenic") %>% 
  left_join(., miRNAtargets, by = "gene", copy = FALSE) %>% 
  mutate(locus_type = case_when(is.na(miRNA) ~ locus_type, 
                                TRUE ~ "miRNA_target")) %>% 
  select(gene, baseMean, log2FoldChange, pvalue, padj, locus_type, genotype) %>% 
  filter(genotype == "pel1.1") %>% 
  filter(locus_type == "miRNA_target")

pel1.1_miRNAtargets_filtered <- c(part2$gene)

new_table_for_paper <- rbind(part1, part2)

#write.xlsx(x = new_table_for_paper, file = "/binf-isilon/PBgrp/xpj980/figures/DESeq2ExperimentBWOintergenicPEL1update.xlsx")

pel2.5_miRNAtargets_filtered <- c(part2$gene)

### Suppl. FigS4 ####
#overlap between pelota alleles 
colors <- c("#986D8E", "#87A8A4")

pel1.2_miRNAtargets <- test %>% 
  filter(genotype == "pel1.2" & locus_type == "miRNA_target" & log2FoldChange > 0) %>% 
  pull(gene)

s4 <- list(pel1.1 = c(pel1.1_miRNAtargets_filtered),
           pel1.2 = c(pel1.2_miRNAtargets))

hits <- plot(euler(s4, shape = "ellipse"), quantities = TRUE, fill = colors)
ggsave(hits, filename = "/binf-isilon/PBgrp/xpj980/figures/pelotaoverlap.svg")


### look for targets in hen2 and rrp4 that overlaps with Voinnet nuclear/cytoplasmic shuttling paper (revision)
look <- test %>% 
  filter(genotype %in% c("ski2.5", "rrp4", "hen2.5")) %>% 
  filter(locus_type == "miRNA_target") %>% 
  filter(gene %in% c("AT1G12820", "AT1G48410", "AT1G53230", "AT1G62670", "AT1G62910", "AT1G62930", "AT1G63070", "AT1G63080", "AT1G63130", "AT1G63150", "AT1G63230", "AT3G23690", "AT3G26810", "AT4G36920", "AT5G16640", "AT5G41610", "AT1G12300", "AT1G12460", "AT1G12620", "AT1G12775", "AT1G62914", "AT1G63630", "AT2G46800")) %>% 
  filter(padj < 0.05) %>% 
  mutate(up_or_down = case_when(log2FoldChange > 0 ~ "up", 
                                log2FoldChange < 0 ~ "down")) %>%
  group_by(genotype,up_or_down) %>%
  add_tally(name = "sum") %>%
  distinct() %>%
  ungroup() %>%
  mutate(sum= case_when(up_or_down == "down" ~ -1*sum,
                        up_or_down == "up" ~ 1*sum)) %>% 
  distinct(genotype, up_or_down, sum)

revision <- left_join(look, miRNAtargets, by = "gene") %>% select(gene, genotype, up_or_down, miRNA)
write.xlsx(revision, file = "/binf-isilon/PBgrp/xpj980/figures/BolognaTable.xlsx")

(barplot_across_geno <- ggplot(look, aes(x= genotype, y=sum, fill = up_or_down)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = sum), vjust=1.6, color="white", size = 3) +
    labs(y= "number of genes\n", x = "genotype") + 
    cowplot::theme_cowplot() + 
    scale_x_discrete(position = "top") + 
    ggtitle("miRNA targets from Bologna et al.") + 
    theme(strip.text.x = element_text(size = 12,
                                      face = "bold",
                                      color = "#22292F"),
          strip.background = element_blank(),
          plot.title = element_text(size = 16, 
                                    face = "bold",
                                    color = "#22292F", 
                                    margin = margin(b = 8)), 
          axis.title.x = element_text(size = 10,
                                      face = "bold",
                                      color = "#22292F", 
                                      margin = margin(t = 15)),
          axis.title.y = element_text(size = 10,
                                      face = "bold",
                                      color = "#22292F"), 
          axis.text.x = element_text(size = 8), 
          axis.text.y = element_text(size = 8),
          legend.title = element_text(size = 9,
                                      face = "bold",
                                      color = "#22292F"), 
          legend.text = element_text(size = 8)))

#ggsave(barplot_across_geno, file = "/binf-isilon/PBgrp/xpj980/figures/Voinnet.svg", device = "svg", units = "in", width = 15, height = 15, dpi = 300)

###### Figure 3A ######
barplot_table <- left_join(x = gathered_DESeqresults_pelota, y = miRNAtargets, by = "gene", copy = FALSE) %>%
  filter(padj < 0.05) 

to_plot <- barplot_table %>% 
  mutate(up_or_down = case_when(log2FoldChange > 0 ~ "up", 
                                log2FoldChange < 0 ~ "down"), 
         miRNA = case_when(is.na(miRNA) ~ "other", 
                           TRUE ~ "miRNA_target")) %>% 
  select(genotype, miRNA, up_or_down)

to_tally <- to_plot %>% 
  group_by(genotype, miRNA,up_or_down) %>%
  add_tally(name = "sum") %>%
  distinct() %>%
  ungroup() %>%
  mutate(sum= case_when(up_or_down == "down" ~ -1*sum,
                        up_or_down == "up" ~ 1*sum),
         plot = paste0(miRNA, up_or_down))

to_tally$genotype <- factor(to_tally$genotype, levels = c("ski2.5", "SAILpelota", "GKpelota"))
to_tally$plot <- factor(to_tally$plot, levels = c("miRNA_targetup", "miRNA_targetdown", "otherdown", "otherup"))

(barplot_across_geno <- ggplot(to_tally, aes(x= miRNA, y=sum, fill = plot)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = sum), vjust=1.6, color="white", size = 3) +
    scale_fill_brewer("more or less sRNAs", 
                      palette = "RdGy",    
                      labels = c("known miRNA target - more sRNAs", "known miRNA target - less sRNAs", "other gene - less sRNAs", "other gene - more sRNAs")) +
    facet_wrap(~genotype, strip.position = "bottom") +
    labs(y= "number of genes\n", x = "gene type\n") + 
    cowplot::theme_cowplot() + 
    scale_x_discrete(position = "top") + 
    ggtitle("Genes with more or less sRNAs compared to WT") + 
    theme(strip.text.x = element_text(size = 12,
                                      face = "bold",
                                      color = "#22292F"),
          strip.background = element_blank(),
          plot.title = element_text(size = 16, 
                                    face = "bold",
                                    color = "#22292F", 
                                    margin = margin(b = 8)), 
          axis.title.x = element_text(size = 10,
                                      face = "bold",
                                      color = "#22292F", 
                                      margin = margin(t = 15)),
          axis.title.y = element_text(size = 10,
                                      face = "bold",
                                      color = "#22292F"), 
          axis.text.x = element_text(size = 8), 
          axis.text.y = element_text(size = 8),
          legend.title = element_text(size = 9,
                                      face = "bold",
                                      color = "#22292F"), 
          legend.text = element_text(size = 8)))

(barplot_across_geno <- ggplot(to_tally, aes(x= genotype, y=sum, fill = plot)) +
    geom_col(alpha = 0.8, width = 0.5) +
    geom_text(aes(label = sum), vjust=1.6, color="black", size = 3) +
    scale_fill_brewer("more or less sRNAs", 
                      palette = "RdGy",    
                      labels = c("known miRNA target - more sRNAs", "known miRNA target - less sRNAs", "other gene - less sRNAs", "other gene - more sRNAs")) +
    facet_wrap(~miRNA, nrow = 1, scales = "free") +
    labs(y= "number of genes\n", x = "gene type\n") + 
    cowplot::theme_cowplot() + 
    geom_hline(yintercept = 0)+
    scale_x_discrete(position = "bottom") + 
    ggtitle("Genes with more or less sRNAs compared to WT") + 
    theme(strip.text.x = element_text(size = 12,
                                      face = "bold",
                                      color = "#22292F"),
          strip.background = element_blank(),
          plot.title = element_text(size = 16, 
                                    face = "bold",
                                    color = "#22292F", 
                                    margin = margin(b = 8)), 
          axis.title.x = element_text(size = 10,
                                      face = "bold",
                                      color = "#22292F", 
                                      margin = margin(t = 15)),
          axis.title.y = element_text(size = 10,
                                      face = "bold",
                                      color = "#22292F"), 
          axis.text.x = element_text(size = 8), 
          axis.text.y = element_text(size = 8),
          legend.title = element_text(size = 9,
                                      face = "bold",
                                      color = "#22292F"), 
          legend.text = element_text(size = 8))) + 
  theme(aspect.ratio = 1)

#ggsave(barplot_across_geno, file = "/binf-isilon/PBgrp/xpj980/figures/20210423Figure3Aaraport.svg", device = "svg", units = "in", width = 15, height = 15, dpi = 300)

#### Figure 3B ######

left_joined <- left_joined <- left_join(x = gathered_DESeqresults_pelota, y = miRNAtargets, by = "gene")

ski2_list <- left_joined[left_joined$genotype=="ski2.5",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(!is.na(miRNA)) %>%
  pull(gene) 

GKpelota_list <- left_joined[left_joined$genotype=="GKpelota",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(!is.na(miRNA)) %>%
  pull(gene) 

SAILpelota_list <- left_joined[left_joined$genotype=="SAILpelota",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(!is.na(miRNA)) %>%
  pull(gene) 

# Object for Carlotta to make Venn Diagrams
mariaslistfig35 <- list("ski2_2019_list" = ski2_list, "GKpelota_list" = GKpelota_list, "SAILpelota_list" = SAILpelota_list, "ski2_2019_nontarg_list" = ski2_other, "GKpelota_nontarg_list" = GKpelota_other,  "SAILpelota_nontarg_list" = SAILpelota_other, "hen2_list" = hen2_list, "rrp4_list" = rrp4_list) 
save(mariaslistfig35, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/2019exosomeVenns.RData")

# and the "other genes overlap"
## and the "other genes" in the mutants, how many do they have in common? 

ski2_other <- left_joined[left_joined$genotype=="ski2.5",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(is.na(miRNA)) %>%
  pull(gene)   

GKpelota_other <- left_joined[left_joined$genotype=="GKpelota",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(is.na(miRNA)) %>%
  pull(gene) 

SAILpelota_other <- left_joined[left_joined$genotype=="SAILpelota",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(is.na(miRNA)) %>%
  pull(gene) 

#### scatterplot ski2-5 & pelota alleles ####

miRNAtargets_toscatter_pel1.1 <- left_join(gathered_DESeqresults, to_join, by = "gene") %>% 
  mutate(locus_type = case_when(is.na(locus_type) ~ "intergenic", 
                                T ~ locus_type)) %>% 
  filter(log2FoldChange > 0.9) %>% 
  mutate(genotype = case_when(genotype == "GKpelota" ~ "pel1.2", 
                              T ~ genotype), 
         genotype = case_when(genotype == "SAILpelota" ~ "pel1.1", 
                              T ~ genotype)) %>% 
  filter(locus_type !="intergenic") %>% 
  left_join(., miRNAtargets, by = "gene", copy = FALSE) %>% 
  mutate(locus_type = case_when(is.na(miRNA) ~ locus_type, 
                                TRUE ~ "miRNA_target")) %>% 
  select(gene, baseMean, log2FoldChange, pvalue, padj, locus_type, genotype) %>% 
  filter(genotype == "pel1.1") %>% 
  filter(locus_type == "miRNA_target") %>%
  dplyr::rename(log2FCSAIL = log2FoldChange) %>% 
  full_join(., y =miRNAtargets) %>%
  mutate(diffexp_pel = case_when(is.na(genotype) ~ "no", 
                                 T ~ "yes"))

miRNAtargets_toscatter_pel1.2 <- left_joined <- left_join(x = gathered_DESeqresults_pelota, y = miRNAtargets, by = "gene") %>% 
  filter(target_name != "NA" & genotype == "GKpelota") %>% 
  dplyr::rename(log2FCGKpel = log2FoldChange) %>% 
  mutate(diffexp_pel = case_when(padj < 0.05 ~ "yes", 
                                 T ~ "no"))

miRNAtargets_toscatter_ski2.5 <- left_joined <- left_join(x = gathered_DESeqresults_pelota, y = miRNAtargets, by = "gene") %>% 
  filter(target_name != "NA" & genotype == "ski2.5") %>% 
  dplyr::rename(log2FCski2 = log2FoldChange) %>% 
  mutate(diffexp_ski2 = case_when(padj < 0.05 ~ "yes", 
                                  T ~ "no"))

to_plot <- left_join(miRNAtargets_toscatter_pel1.2, miRNAtargets_toscatter_ski2.5, by = "gene") %>% 
  mutate(color = case_when(diffexp_ski2 == "yes" & diffexp_pel == "no" ~ "ski2", 
                           T ~ "none"), 
         color = case_when(diffexp_pel == "yes" & diffexp_ski2== "no" ~ "pel1", 
                           T ~ color), 
         color =  case_when(diffexp_pel == "yes" & diffexp_ski2== "yes" ~ "both", 
                            T ~ color)) %>%
  #filter(log2FCGKpel > 0 & log2FCski2 > 0) %>% 
  filter(color %in% c("ski2", "pel1", "both"))

res <- cor.test(to_plot$log2FCGKpel, to_plot$log2FCski2, 
                method = "pearson")
res

cor(to_plot$log2FCGKpel, to_plot$log2FCski2, method = "pearson")

plot <- ggplot(to_plot, aes(x=log2FCGKpel, y = log2FCski2, col = color)) + 
  geom_point() + 
  geom_text_repel(aes(label=ifelse(color == "ski2",as.character(target_name.x),'')),hjust=0,vjust=0, max.overlaps = Inf) + 
  geom_text_repel(aes(label=ifelse(color == "pel1",as.character(target_name.x),'')),hjust=0,vjust=0, max.overlaps = Inf) + 
  #coord_equal() +
  geom_hline(yintercept = 0, size = 0.2, linetype = "dashed") +
  geom_vline(xintercept = 0, size = 0.2, linetype = "dashed") +
  #scale_x_continuous(limits = c(0.01, NA)) +
  #scale_y_continuous(limits = c(0.01, NA)) +
  geom_abline(intercept = 0, slope = 1) +
  cowplot::theme_cowplot() +
  scale_colour_manual(values=c("none"="#aaaaaa","both"="#766161", "ski2" = "#3d84b8", "pel1" = "#ff8474")) +
  ggtitle(label = "log2FC of miRNA targets in ski2-5 and pel1-2 with higher siRNA levels", subtitle = paste0("corr. coeff. r = ",cor(to_plot$log2FCGKpel, to_plot$log2FCski2, method = "pearson"))) +
  #scale_color_manual(values = wes_palette(n=4, name="Darjeeling1"))+
  labs(y= "log2FC ski2-5", x = "log2FC pel1-2") +
  theme(aspect.ratio = 1)
#ggsave(plot, filename = "/binf-isilon/PBgrp/xpj980/figures/scatterski2pel12final.svg", units = "in", width = 15, height = 15, dpi = 300)
#ggsave(plot, filename = "/binf-isilon/PBgrp/xpj980/figures/3scatterski2pel12final.png", units = "in", width = 15, height = 15, dpi = 300)

to_plot <- full_join(miRNAtargets_toscatter_pel1.1, miRNAtargets_toscatter_ski2.5, by = "gene") %>% 
  mutate(color = case_when(diffexp_ski2 == "yes" & diffexp_pel == "no" ~ "ski2", 
                           T ~ "none"), 
         color = case_when(diffexp_pel == "yes" & diffexp_ski2== "no" ~ "pel1", 
                           T ~ color), 
         color =  case_when(diffexp_pel == "yes" & diffexp_ski2== "yes" ~ "both", 
                            T ~ color)) %>%
  #filter(log2FCSAIL > 0 & log2FCski2 > 0) %>%
  filter(color %in% c("ski2", "pel1", "both"))

findlog2SAIL <- filter(to_plot, color == "ski2") %>% pull(gene)
findlog2SAIL2 <- gathered_DESeqresults %>% 
  filter(gene %in% findlog2SAIL & genotype == "SAILpelota") 

test <- left_join(to_plot, findlog2SAIL2, by = "gene") %>% 
  mutate(log2FCSAIL = case_when(is.na(log2FCSAIL) ~ log2FoldChange, 
                                T ~ log2FCSAIL))


cor(test$log2FCSAIL, test$log2FCski2, method = "pearson")

plot <- ggplot(test, aes(x=log2FCSAIL, y = log2FCski2, col = color)) + 
  geom_point() + 
  geom_text_repel(aes(label=ifelse(color == "ski2",as.character(target_name.x),'')),hjust=0,vjust=0, max.overlaps = Inf) + 
  geom_text_repel(aes(label=ifelse(color == "pel1",as.character(target_name.x),'')),hjust=0,vjust=0, max.overlaps = Inf) + 
  coord_equal() +
  #scale_x_continuous(limits = c(0.01, NA)) +
  #scale_y_continuous(limits = c(0.01, NA)) +
  geom_abline(intercept = log2(1), slope = 1) +
  geom_hline(yintercept = 0, size = 0.2, linetype = "dashed") +
  geom_vline(xintercept = 0, size = 0.2, linetype = "dashed") +
  cowplot::theme_cowplot() +
  scale_colour_manual(values=c("none"="#aaaaaa","both"="#766161", "ski2" = "#3d84b8", "pel1" = "#ff8474")) +
  ggtitle(label = "log2FC of miRNA targets in ski2-5 and pel1-1 with higher siRNA levels", subtitle = paste0("corr. coeff. r = ",cor(test$log2FCSAIL, test$log2FCski2, method = "pearson"))) +
  #scale_color_manual(values = wes_palette(n=4, name="Darjeeling1"))+
  labs(y= "log2FC ski2-5)", x = "log2FC pel1-1") +
  theme(aspect.ratio = 1)
ggsave(plot, filename = "/binf-isilon/PBgrp/xpj980/figures/scatterski2pel11log2filtered.svg", units = "in", width = 15, height = 15, dpi = 300)
#ggsave(plot, filename = "/binf-isilon/PBgrp/xpj980/figures/3scatterski2pel11final.png", units = "in", width = 15, height = 15, dpi = 300)

#### Supplm figure S3 #####
### count miRNAs (miRNA genes) ###

all_annotations <- readGFF("/binf-isilon/PBgrp/xpj980/TAIR/Arabidopsis_thaliana.TAIR10.46.gff3")

genes <- all_annotations %>% 
  as_tibble() %>% 
  dplyr::select(biotype, gene=gene_id, Name) %>% 
  filter(!is.na(gene)) 

joined_sign <- left_join(x = gathered_DESeqresults_pelota, y = genes) %>%
  filter(padj < 0.05 , biotype == "miRNA") 

## list all miRNAs (miRNA genes) up- or downregulated in genotypes 

a_list <- lapply(joined_sign$genotype %>% unique(), function(x) {
  joined_sign %>% filter(genotype==x) %>% pull(Name) %>% unique() %>% sort()
})
names(a_list) <- joined_sign$genotype %>% unique()
a_list

### scatterplots ###

joined_sign <- left_join(x = gathered_DESeqresults_pelota, y = genes) %>%
  filter(genotype == "ski2.5", padj < 0.05 , biotype == "miRNA") 

to_gather <- DESeq_table_summed %>% 
  rownames_to_column(var = "gene") %>% 
  glimpse()

keycol <- "genotype"
valcol <- "count"
gathercols <- c(names(to_gather[,2:20]))

test <- gather_(to_gather, keycol, valcol, gathercols)

normalized <- separate(test, col = genotype, into = c("genotyp", "rep"), sep = "_") %>%
  group_by(genotyp,rep) %>%
  mutate(RPM = count*1000000/sum(count)) %>%
  group_by(gene, genotyp) %>% 
  summarise(mean_RPM = mean(RPM))

miRNAs_RPM_table <- normalized %>% 
  left_join(., genes, by = "gene") %>% 
  filter(., biotype %in% "miRNA") %>% 
  spread(., genotyp, mean_RPM) %>% 
  drop_na(., Name) 
save(miRNAs_RPM_table)
sign_miRNAs <- left_join(miRNAs_RPM_table, joined_sign, by = "gene") %>% 
  mutate(sig = case_when(is.na(Name.y) ~ "no", 
                         TRUE ~ "yes")) %>%
  select(gene, biotype.x, Name.x, GKpelota, ski2.5, hen2.5, wt, rrp4, rrp45a, SAILpelota, sig)

left_joined <- left_join(x = gathered_DESeqresults_pelota, y = miRNAtargets, by = "gene", copy = FALSE)

ski2_test <- left_joined[left_joined$genotype=="ski2.5",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(!is.na(miRNA)) 

mir_parse <- function(mir_names, simplify=T) {
  result <- lapply(mir_names, function(mir_name) {
    mir <- mir_name
    mir_base <- gsub("([a-zA-z]*[1-9]*)[a-z]*", "\\1", mir)
    mir_letters <- gsub("[a-zA-z]*[1-9]*([A-Za-z]*)", "\\1", mir) %>% str_split("", simplify = T)
    mir_list <- paste0(mir_base, c("", mir_letters))
    return(mir_list)
  })
  if(simplify) {result <- unlist(result)}
  return(result)
}

miRNA_genes <- mir_parse(ski2_test$miRNA, simplify = T) %>% 
  toupper()

to_merge <- subset(sign_miRNAs, sign_miRNAs$Name.x %in% miRNA_genes) %>% 
  mutate(sign_target = "yes")

to_plot <- left_join(sign_miRNAs, to_merge) %>% 
  mutate(sign_target = case_when(is.na(sign_target) ~ "no", 
                                 T ~ "yes"))
to_plot3 <- to_plot %>% 
  mutate(both = case_when(sign_target == "yes" & sig == "no" ~ "target.miRNA", 
                          sign_target == "no" & sig == "yes" ~ "diff.exp.miRNA", 
                          sign_target == "no" & sig == "no" ~ "miRNA", 
                          sign_target == "yes" & sig == "yes" ~ "diff.exp.miRNA.target")) %>% 
  glimpse()

## this one is correct
ski2_miRNAs <- to_plot3 %>% 
  #mutate(both = factor(both, c("diff.exp.miRNA.target", "diff.exp.miRNA", "target.miRNA", "miRNA"))) %>% 
  #mutate(both = factor(both, c("miRNA", "target.miRNA", "diff.exp.miRNA", "diff.exp.miRNA.target"))) %>%
  ggplot(., aes(x=log2(wt), y=log2(ski2.5), color = both)) +
  #ggplot(.,aes(x=wt, y=rrp45b, color = fct_relevel(both, c("diff.exp.miRNA.target", "diff.exp.miRNA", "target.miRNA", "miRNA")))) +
  #ggplot(.,aes(x=wt, y=rrp45b, color = fct_relevel(both, c("miRNA", "target.miRNA", "diff.exp.miRNA", "diff.exp.miRNA.target")))) +
  geom_point() +
  scale_x_continuous(limits = c(0.01, NA)) +
  scale_y_continuous(limits = c(0.01, NA)) +
  #scale_y_log10(limits = c(0.01,NA), breaks = c(0, 1, 10, 100, 1000)) +
  #scale_x_log10(limits = c(0.01,NA), breaks = c(0, 1, 10, 100, 1000)) +
  coord_equal() +
  #geom_hline(yintercept = 1, colour="#990000", linetype="dashed") + 
  #geom_vline(xintercept = 1, colour="#990000", linetype="dashed") +
  geom_abline(intercept = 0, slope = 1) +
  cowplot::theme_cowplot() +
  geom_text_repel(aes(label=ifelse(both == "diff.exp.miRNA.target",as.character(Name.x),'')),hjust=0,vjust=0) +
  geom_text_repel(aes(label=ifelse(both == "target.miRNA",as.character(Name.x),'')),hjust=0,vjust=0) +
  scale_colour_manual(values=c("diff.exp.miRNA"="#aa4a30","target.miRNA"="#484018", "diff.exp.miRNA.target" = "#e40017", "miRNA" = "#f4c983")) +
  #geom_text_repel(data=subset(to_plot3, rrp45b > 50 & wt > 50), aes(label=Name.x)) +
  ggtitle(label = "ski2-5", subtitle = "scatterplot of miRNA genes in ski2-5 vs WT") +
  #scale_color_manual(values = wes_palette(n=4, name="Darjeeling1"))+
  labs(y= "log2(ski2.5[meanRPM])", x = "log2(WT[meanRPM])") +
  theme(aspect.ratio = 1)
ggsave(ski2_miRNAs, file = "/binf-isilon/PBgrp/xpj980/figures/20210424miRNAgenesscatterski2secondlibaraport.svg", units = "in", width = 15, height = 15, dpi = 300)

## now make the same scatterplot for pel1-1 (SAIL) and pel1-2 (GK) 

joined_sign <- left_join(x = gathered_DESeqresults_pelota, y = genes) %>%
  filter(genotype == "SAILpelota", padj < 0.05 , biotype == "miRNA") 

to_gather <- DESeq_table_summed %>% 
  rownames_to_column(var = "gene") %>% 
  glimpse()

keycol <- "genotype"
valcol <- "count"
gathercols <- c(names(to_gather[,2:20]))

test <- gather_(to_gather, keycol, valcol, gathercols)

normalized <- separate(test, col = genotype, into = c("genotyp", "rep"), sep = "_") %>%
  group_by(genotyp,rep) %>%
  mutate(RPM = count*1000000/sum(count)) %>%
  group_by(gene, genotyp) %>% 
  summarise(mean_RPM = mean(RPM))

miRNAs_RPM_table <- normalized %>% 
  left_join(., genes, by = "gene") %>% 
  filter(., biotype %in% "miRNA") %>% 
  spread(., genotyp, mean_RPM) %>% 
  drop_na(., Name) 

sign_miRNAs <- left_join(miRNAs_RPM_table, joined_sign, by = "gene") %>% 
  mutate(sig = case_when(is.na(Name.y) ~ "no", 
                         TRUE ~ "yes")) %>%
  select(gene, biotype.x, Name.x, GKpelota, ski2.5, hen2.5, wt, rrp4, rrp45a, SAILpelota, sig)

left_joined <- left_join(x = gathered_DESeqresults_pelota, y = miRNAtargets, by = "gene", copy = FALSE)

pel_test <- left_joined[left_joined$genotype=="SAILpelota",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(!is.na(miRNA)) 

mir_parse <- function(mir_names, simplify=T) {
  result <- lapply(mir_names, function(mir_name) {
    mir <- mir_name
    mir_base <- gsub("([a-zA-z]*[1-9]*)[a-z]*", "\\1", mir)
    mir_letters <- gsub("[a-zA-z]*[1-9]*([A-Za-z]*)", "\\1", mir) %>% str_split("", simplify = T)
    mir_list <- paste0(mir_base, c("", mir_letters))
    return(mir_list)
  })
  if(simplify) {result <- unlist(result)}
  return(result)
}

miRNA_genes <- mir_parse(pel_test$miRNA, simplify = T) %>% 
  toupper()

to_merge <- subset(sign_miRNAs, sign_miRNAs$Name.x %in% miRNA_genes) %>% 
  mutate(sign_target = "yes")

sig_targets_pel1.1_filtered <- miRNAtargets_rawtable %>% 
  select(target_id, miRNA_ID) %>% 
  filter(!duplicated(.$target_id)) %>%
  filter(target_id %in% pel1.1_miRNAtargets_filtered) %>% 
  pull(miRNA_ID)

to_plot <- left_join(sign_miRNAs, to_merge) %>% 
  mutate(sign_target = case_when(gene %in% sig_targets_pel1.1_filtered ~ "yes", 
                                 T ~ "no"))

to_plot3 <- to_plot %>% 
  mutate(both = case_when(sign_target == "yes" ~ "target.miRNA", 
                          sign_target == "no"  ~ "miRNA")) %>% 
  glimpse()

## this one is correct
pel1.1_miRNAs <- to_plot3 %>% 
  #mutate(both = factor(both, c("diff.exp.miRNA.target", "diff.exp.miRNA", "target.miRNA", "miRNA"))) %>% 
  #mutate(both = factor(both, c("miRNA", "target.miRNA", "diff.exp.miRNA", "diff.exp.miRNA.target"))) %>%
  ggplot(., aes(x=log2(wt), y=log2(SAILpelota), color = both)) +
  #ggplot(.,aes(x=wt, y=rrp45b, color = fct_relevel(both, c("diff.exp.miRNA.target", "diff.exp.miRNA", "target.miRNA", "miRNA")))) +
  #ggplot(.,aes(x=wt, y=rrp45b, color = fct_relevel(both, c("miRNA", "target.miRNA", "diff.exp.miRNA", "diff.exp.miRNA.target")))) +
  geom_point() +
  scale_x_continuous(limits = c(0.01, NA)) +
  scale_y_continuous(limits = c(0.01, NA)) +
  #scale_y_log10(limits = c(0.01,NA), breaks = c(0, 1, 10, 100, 1000)) +
  #scale_x_log10(limits = c(0.01,NA), breaks = c(0, 1, 10, 100, 1000)) +
  coord_equal() +
  #geom_hline(yintercept = 1, colour="#990000", linetype="dashed") + 
  #geom_vline(xintercept = 1, colour="#990000", linetype="dashed") +
  geom_abline(intercept = 0, slope = 1) +
  cowplot::theme_cowplot() +
  geom_text_repel(aes(label=ifelse(both == "diff.exp.miRNA.target",as.character(Name.x),'')),hjust=0,vjust=0) +
  geom_text_repel(aes(label=ifelse(both == "target.miRNA",as.character(Name.x),'')),hjust=0,vjust=0) +
  scale_colour_manual(values=c("diff.exp.miRNA"="#aa4a30","target.miRNA"="#484018", "diff.exp.miRNA.target" = "#e40017", "miRNA" = "#f4c983")) +
  #geom_text_repel(data=subset(to_plot3, rrp45b > 50 & wt > 50), aes(label=Name.x)) +
  ggtitle(label = "pel1-1", subtitle = "scatterplot of miRNA genes in pel1-1 vs WT") +
  #scale_color_manual(values = wes_palette(n=4, name="Darjeeling1"))+
  labs(y= "log2(pel1.1[meanRPM])", x = "log2(WT[meanRPM])") +
  theme(aspect.ratio = 1)

ggsave(pel1.1_miRNAs, file = "/binf-isilon/PBgrp/xpj980/figures/miRNAgenesscatterSAILpelaraportUPDATEDrevision.svg", units = "in", width = 15, height = 15, dpi = 300)

# GKpelota
joined_sign <- left_join(x = gathered_DESeqresults_pelota, y = genes) %>%
  filter(genotype == "GKpelota", padj < 0.05 , biotype == "miRNA") 

to_gather <- DESeq_table_summed %>% 
  rownames_to_column(var = "gene") %>% 
  glimpse()

keycol <- "genotype"
valcol <- "count"
gathercols <- c(names(to_gather[,2:20]))

test <- gather_(to_gather, keycol, valcol, gathercols)

normalized <- separate(test, col = genotype, into = c("genotyp", "rep"), sep = "_") %>%
  group_by(genotyp,rep) %>%
  mutate(RPM = count*1000000/sum(count)) %>%
  group_by(gene, genotyp) %>% 
  summarise(mean_RPM = mean(RPM))

miRNAs_RPM_table <- normalized %>% 
  left_join(., genes, by = "gene") %>% 
  filter(., biotype %in% "miRNA") %>% 
  spread(., genotyp, mean_RPM) %>% 
  drop_na(., Name) 

sign_miRNAs <- left_join(miRNAs_RPM_table, joined_sign, by = "gene") %>% 
  mutate(sig = case_when(is.na(Name.y) ~ "no", 
                         TRUE ~ "yes")) %>%
  select(gene, biotype.x, Name.x, GKpelota, ski2.5, hen2.5, wt, rrp4, rrp45a, SAILpelota, sig)

left_joined <- left_join(x = gathered_DESeqresults_pelota, y = miRNAtargets, by = "gene", copy = FALSE)

pel_test <- left_joined[left_joined$genotype=="GKpelota",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(!is.na(miRNA)) 

mir_parse <- function(mir_names, simplify=T) {
  result <- lapply(mir_names, function(mir_name) {
    mir <- mir_name
    mir_base <- gsub("([a-zA-z]*[1-9]*)[a-z]*", "\\1", mir)
    mir_letters <- gsub("[a-zA-z]*[1-9]*([A-Za-z]*)", "\\1", mir) %>% str_split("", simplify = T)
    mir_list <- paste0(mir_base, c("", mir_letters))
    return(mir_list)
  })
  if(simplify) {result <- unlist(result)}
  return(result)
}

miRNA_genes <- mir_parse(pel_test$miRNA, simplify = T) %>% 
  toupper()

to_merge <- subset(sign_miRNAs, sign_miRNAs$Name.x %in% miRNA_genes) %>% 
  mutate(sign_target = "yes")

to_plot <- left_join(sign_miRNAs, to_merge) %>% 
  mutate(sign_target = case_when(is.na(sign_target) ~ "no", 
                                 T ~ "yes"))
to_plot3 <- to_plot %>% 
  mutate(both = case_when(sign_target == "yes" & sig == "no" ~ "target.miRNA", 
                          sign_target == "no" & sig == "yes" ~ "diff.exp.miRNA", 
                          sign_target == "no" & sig == "no" ~ "miRNA", 
                          sign_target == "yes" & sig == "yes" ~ "diff.exp.miRNA.target")) %>% 
  glimpse()

## this one is correct
pel1.2_miRNAs <- to_plot3 %>% 
  #mutate(both = factor(both, c("diff.exp.miRNA.target", "diff.exp.miRNA", "target.miRNA", "miRNA"))) %>% 
  #mutate(both = factor(both, c("miRNA", "target.miRNA", "diff.exp.miRNA", "diff.exp.miRNA.target"))) %>%
  ggplot(., aes(x=log2(wt), y=log2(GKpelota), color = both)) +
  #ggplot(.,aes(x=wt, y=rrp45b, color = fct_relevel(both, c("diff.exp.miRNA.target", "diff.exp.miRNA", "target.miRNA", "miRNA")))) +
  #ggplot(.,aes(x=wt, y=rrp45b, color = fct_relevel(both, c("miRNA", "target.miRNA", "diff.exp.miRNA", "diff.exp.miRNA.target")))) +
  geom_point() +
  scale_x_continuous(limits = c(0.01, NA)) +
  scale_y_continuous(limits = c(0.01, NA)) +
  #scale_y_log10(limits = c(0.01,NA), breaks = c(0, 1, 10, 100, 1000)) +
  #scale_x_log10(limits = c(0.01,NA), breaks = c(0, 1, 10, 100, 1000)) +
  coord_equal() +
  #geom_hline(yintercept = 1, colour="#990000", linetype="dashed") + 
  #geom_vline(xintercept = 1, colour="#990000", linetype="dashed") +
  geom_abline(intercept = 0, slope = 1) +
  cowplot::theme_cowplot() +
  geom_text_repel(aes(label=ifelse(both == "diff.exp.miRNA.target",as.character(Name.x),'')),hjust=0,vjust=0) +
  geom_text_repel(aes(label=ifelse(both == "target.miRNA",as.character(Name.x),'')),hjust=0,vjust=0) +
  scale_colour_manual(values=c("diff.exp.miRNA"="#aa4a30","target.miRNA"="#484018", "diff.exp.miRNA.target" = "#e40017", "miRNA" = "#f4c983")) +
  #geom_text_repel(data=subset(to_plot3, rrp45b > 50 & wt > 50), aes(label=Name.x)) +
  ggtitle(label = "pel1-2", subtitle = "scatterplot of miRNA genes in pel1-2 vs WT") +
  #scale_color_manual(values = wes_palette(n=4, name="Darjeeling1"))+
  labs(y= "log2(pel1.2[meanRPM])", x = "log2(WT[meanRPM])") +
  theme(aspect.ratio = 1)
ggsave(pel1.2_miRNAs, file = "/binf-isilon/PBgrp/xpj980/figures/20210424miRNAgenesscatterGKpelaraport.svg", units = "in", width = 15, height = 15, dpi = 300)

## count miRNAs in rrp4 
all_annotations <- readGFF("/binf-isilon/PBgrp/xpj980/TAIR/Arabidopsis_thaliana.TAIR10.46.gff3")

genes <- all_annotations %>% 
  as_tibble() %>% 
  dplyr::select(biotype, gene=gene_id, Name) %>% 
  filter(!is.na(gene)) 

joined_sign <- left_join(x = gathered_DESeqresults_rrp4, y = genes) %>%
  filter(padj < 0.05 , biotype == "miRNA") 

## list all miRNAs (miRNA genes) up- or downregulated in genotypes 

a_list <- lapply(joined_sign$genotype %>% unique(), function(x) {
  joined_sign %>% filter(genotype==x) %>% pull(Name) %>% unique() %>% sort()
})
names(a_list) <- joined_sign$genotype %>% unique()
a_list

# rrp4 

joined_sign <- left_join(x = gathered_DESeqresults_rrp4, y = genes) %>%
  filter(genotype == "rrp4", padj < 0.05 , biotype == "miRNA") 

to_gather <- DESeq_table_summed %>% 
  rownames_to_column(var = "gene") %>% 
  glimpse()

keycol <- "genotype"
valcol <- "count"
gathercols <- c(names(to_gather[,2:20]))

test <- gather_(to_gather, keycol, valcol, gathercols)

normalized <- separate(test, col = genotype, into = c("genotyp", "rep"), sep = "_") %>%
  group_by(genotyp,rep) %>%
  mutate(RPM = count*1000000/sum(count)) %>%
  group_by(gene, genotyp) %>% 
  summarise(mean_RPM = mean(RPM))

miRNAs_RPM_table <- normalized %>% 
  left_join(., genes, by = "gene") %>% 
  filter(., biotype %in% "miRNA") %>% 
  spread(., genotyp, mean_RPM) %>% 
  drop_na(., Name) 

sign_miRNAs <- left_join(miRNAs_RPM_table, joined_sign, by = "gene") %>% 
  mutate(sig = case_when(is.na(Name.y) ~ "no", 
                         TRUE ~ "yes")) %>%
  select(gene, biotype.x, Name.x, GKpelota, ski2.5, hen2.5, wt, rrp4, rrp45a, SAILpelota, sig)

left_joined <- left_join(x = gathered_DESeqresults_rrp4, y = miRNAtargets, by = "gene", copy = FALSE)

rrp4_test <- left_joined[left_joined$genotype=="rrp4",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(!is.na(miRNA)) 

mir_parse <- function(mir_names, simplify=T) {
  result <- lapply(mir_names, function(mir_name) {
    mir <- mir_name
    mir_base <- gsub("([a-zA-z]*[1-9]*)[a-z]*", "\\1", mir)
    mir_letters <- gsub("[a-zA-z]*[1-9]*([A-Za-z]*)", "\\1", mir) %>% str_split("", simplify = T)
    mir_list <- paste0(mir_base, c("", mir_letters))
    return(mir_list)
  })
  if(simplify) {result <- unlist(result)}
  return(result)
}

miRNA_genes <- mir_parse(rrp4_test$miRNA, simplify = T) %>% 
  toupper()

to_merge <- subset(sign_miRNAs, sign_miRNAs$Name.x %in% miRNA_genes) %>% 
  mutate(sign_target = "yes")

to_plot <- left_join(sign_miRNAs, to_merge) %>% 
  mutate(sign_target = case_when(is.na(sign_target) ~ "no", 
                                 T ~ "yes"))
to_plot3 <- to_plot %>% 
  mutate(both = case_when(sign_target == "yes" & sig == "no" ~ "target.miRNA", 
                          sign_target == "no" & sig == "yes" ~ "diff.exp.miRNA", 
                          sign_target == "no" & sig == "no" ~ "miRNA", 
                          sign_target == "yes" & sig == "yes" ~ "diff.exp.miRNA.target")) %>% 
  glimpse()

## this one is correct
rrp4_miRNAs <- to_plot3 %>% 
  #mutate(both = factor(both, c("diff.exp.miRNA.target", "diff.exp.miRNA", "target.miRNA", "miRNA"))) %>% 
  #mutate(both = factor(both, c("miRNA", "target.miRNA", "diff.exp.miRNA", "diff.exp.miRNA.target"))) %>%
  ggplot(., aes(x=log2(wt), y=log2(rrp4), color = both)) +
  #ggplot(.,aes(x=wt, y=rrp45b, color = fct_relevel(both, c("diff.exp.miRNA.target", "diff.exp.miRNA", "target.miRNA", "miRNA")))) +
  #ggplot(.,aes(x=wt, y=rrp45b, color = fct_relevel(both, c("miRNA", "target.miRNA", "diff.exp.miRNA", "diff.exp.miRNA.target")))) +
  geom_point() +
  scale_x_continuous(limits = c(0.01, NA)) +
  scale_y_continuous(limits = c(0.01, NA)) +
  #scale_y_log10(limits = c(0.01,NA), breaks = c(0, 1, 10, 100, 1000)) +
  #scale_x_log10(limits = c(0.01,NA), breaks = c(0, 1, 10, 100, 1000)) +
  coord_equal() +
  #geom_hline(yintercept = 1, colour="#990000", linetype="dashed") + 
  #geom_vline(xintercept = 1, colour="#990000", linetype="dashed") +
  geom_abline(intercept = 0, slope = 1) +
  cowplot::theme_cowplot() +
  geom_text_repel(aes(label=ifelse(both == "diff.exp.miRNA.target",as.character(Name.x),'')),hjust=0,vjust=0) +
  geom_text_repel(aes(label=ifelse(both == "target.miRNA",as.character(Name.x),'')),hjust=0,vjust=0) +
  scale_colour_manual(values=c("diff.exp.miRNA"="#aa4a30","target.miRNA"="#484018", "diff.exp.miRNA.target" = "#e40017", "miRNA" = "#f4c983")) +
  #geom_text_repel(data=subset(to_plot3, rrp45b > 50 & wt > 50), aes(label=Name.x)) +
  ggtitle(label = "rrp4", subtitle = "scatterplot of miRNA genes in rrp4 vs WT") +
  #scale_color_manual(values = wes_palette(n=4, name="Darjeeling1"))+
  labs(y= "log2(rrp4[meanRPM])", x = "log2(WT[meanRPM])") +
  theme(aspect.ratio = 1)
ggsave(rrp4_miRNAs, file = "/binf-isilon/PBgrp/xpj980/figures/20210424miRNAgenesscatterrrp4araport.svg", units = "in", width = 15, height = 15, dpi = 300)


## hen2.5
joined_sign <- left_join(x = gathered_DESeqresults_rrp4, y = genes) %>%
  filter(genotype == "hen2.5", padj < 0.05 , biotype == "miRNA") 

to_gather <- DESeq_table_summed %>% 
  rownames_to_column(var = "gene") %>% 
  glimpse()

keycol <- "genotype"
valcol <- "count"
gathercols <- c(names(to_gather[,2:20]))

test <- gather_(to_gather, keycol, valcol, gathercols)

normalized <- separate(test, col = genotype, into = c("genotyp", "rep"), sep = "_") %>%
  group_by(genotyp,rep) %>%
  mutate(RPM = count*1000000/sum(count)) %>%
  group_by(gene, genotyp) %>% 
  summarise(mean_RPM = mean(RPM))

miRNAs_RPM_table <- normalized %>% 
  left_join(., genes, by = "gene") %>% 
  filter(., biotype %in% "miRNA") %>% 
  spread(., genotyp, mean_RPM) %>% 
  drop_na(., Name) 

sign_miRNAs <- left_join(miRNAs_RPM_table, joined_sign, by = "gene") %>% 
  mutate(sig = case_when(is.na(Name.y) ~ "no", 
                         TRUE ~ "yes")) %>%
  select(gene, biotype.x, Name.x, GKpelota, ski2.5, hen2.5, wt, rrp4, rrp45a, SAILpelota, sig)

left_joined <- left_join(x = gathered_DESeqresults_rrp4, y = miRNAtargets, by = "gene", copy = FALSE)

hen2_test <- left_joined[left_joined$genotype=="hen2.5",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(!is.na(miRNA)) 

mir_parse <- function(mir_names, simplify=T) {
  result <- lapply(mir_names, function(mir_name) {
    mir <- mir_name
    mir_base <- gsub("([a-zA-z]*[1-9]*)[a-z]*", "\\1", mir)
    mir_letters <- gsub("[a-zA-z]*[1-9]*([A-Za-z]*)", "\\1", mir) %>% str_split("", simplify = T)
    mir_list <- paste0(mir_base, c("", mir_letters))
    return(mir_list)
  })
  if(simplify) {result <- unlist(result)}
  return(result)
}

miRNA_genes <- mir_parse(hen2_test$miRNA, simplify = T) %>% 
  toupper()

to_merge <- subset(sign_miRNAs, sign_miRNAs$Name.x %in% miRNA_genes) %>% 
  mutate(sign_target = "yes")

to_plot <- left_join(sign_miRNAs, to_merge) %>% 
  mutate(sign_target = case_when(is.na(sign_target) ~ "no", 
                                 T ~ "yes"))
to_plot3 <- to_plot %>% 
  mutate(both = case_when(sign_target == "yes" & sig == "no" ~ "target.miRNA", 
                          sign_target == "no" & sig == "yes" ~ "diff.exp.miRNA", 
                          sign_target == "no" & sig == "no" ~ "miRNA", 
                          sign_target == "yes" & sig == "yes" ~ "diff.exp.miRNA.target")) %>% 
  glimpse()

## this one is correct
hen2_miRNAs <- to_plot3 %>% 
  #mutate(both = factor(both, c("diff.exp.miRNA.target", "diff.exp.miRNA", "target.miRNA", "miRNA"))) %>% 
  #mutate(both = factor(both, c("miRNA", "target.miRNA", "diff.exp.miRNA", "diff.exp.miRNA.target"))) %>%
  ggplot(., aes(x=log2(wt), y=log2(hen2.5), color = both)) +
  #ggplot(.,aes(x=wt, y=rrp45b, color = fct_relevel(both, c("diff.exp.miRNA.target", "diff.exp.miRNA", "target.miRNA", "miRNA")))) +
  #ggplot(.,aes(x=wt, y=rrp45b, color = fct_relevel(both, c("miRNA", "target.miRNA", "diff.exp.miRNA", "diff.exp.miRNA.target")))) +
  geom_point() +
  scale_x_continuous(limits = c(0.01, NA)) +
  scale_y_continuous(limits = c(0.01, NA)) +
  #scale_y_log10(limits = c(0.01,NA), breaks = c(0, 1, 10, 100, 1000)) +
  #scale_x_log10(limits = c(0.01,NA), breaks = c(0, 1, 10, 100, 1000)) +
  coord_equal() +
  #geom_hline(yintercept = 1, colour="#990000", linetype="dashed") + 
  #geom_vline(xintercept = 1, colour="#990000", linetype="dashed") +
  geom_abline(intercept = 0, slope = 1) +
  cowplot::theme_cowplot() +
  geom_text_repel(aes(label=ifelse(both == "diff.exp.miRNA.target",as.character(Name.x),'')),hjust=0,vjust=0) +
  geom_text_repel(aes(label=ifelse(both == "target.miRNA",as.character(Name.x),'')),hjust=0,vjust=0) +
  #geom_text_repel(data=subset(to_plot3, rrp45b > 50 & wt > 50), aes(label=Name.x)) +
  ggtitle(label = "hen2.5", subtitle = "scatterplot of miRNA genes in hen2-5 vs WT") +
  #scale_color_manual(values = wes_palette(n=4, name="Darjeeling1"))+
  scale_colour_manual(values=c("diff.exp.miRNA"="#aa4a30","target.miRNA"="#484018", "diff.exp.miRNA.target" = "#e40017", "miRNA" = "#f4c983")) +
  labs(y= "log2(hen2.5[meanRPM])", x = "log2(WT[meanRPM])") +
  theme(aspect.ratio = 1)
ggsave(hen2_miRNAs, file = "/binf-isilon/PBgrp/xpj980/figures/20210424miRNAgenesscatterhen2araport.svg", units = "in", width = 15, height = 15, dpi = 300)

## count miRNAs in rrp45a 
all_annotations <- readGFF("/binf-isilon/PBgrp/xpj980/TAIR/Arabidopsis_thaliana.TAIR10.46.gff3")

genes <- all_annotations %>% 
  as_tibble() %>% 
  dplyr::select(biotype, gene=gene_id, Name) %>% 
  filter(!is.na(gene)) 

joined_sign <- left_join(x = gathered_DESeqresults_rrp45a, y = genes) %>%
  filter(padj < 0.05 , biotype == "miRNA") 

## list all miRNAs (miRNA genes) up- or downregulated in genotypes 

a_list <- lapply(joined_sign$genotype %>% unique(), function(x) {
  joined_sign %>% filter(genotype==x) %>% pull(Name) %>% unique() %>% sort()
})
names(a_list) <- joined_sign$genotype %>% unique()
a_list

## rrp45a
joined_sign <- left_join(x = gathered_DESeqresults_rrp45a, y = genes) %>%
  filter(genotype == "rrp45a", padj < 0.05 , biotype == "miRNA") 

to_gather <- DESeq_table_summed %>% 
  rownames_to_column(var = "gene") %>% 
  glimpse()

keycol <- "genotype"
valcol <- "count"
gathercols <- c(names(to_gather[,2:20]))

test <- gather_(to_gather, keycol, valcol, gathercols)

normalized <- separate(test, col = genotype, into = c("genotyp", "rep"), sep = "_") %>%
  group_by(genotyp,rep) %>%
  mutate(RPM = count*1000000/sum(count)) %>%
  group_by(gene, genotyp) %>% 
  summarise(mean_RPM = mean(RPM))

miRNAs_RPM_table <- normalized %>% 
  left_join(., genes, by = "gene") %>% 
  filter(., biotype %in% "miRNA") %>% 
  spread(., genotyp, mean_RPM) %>% 
  drop_na(., Name) 

sign_miRNAs <- left_join(miRNAs_RPM_table, joined_sign, by = "gene") %>% 
  mutate(sig = case_when(is.na(Name.y) ~ "no", 
                         TRUE ~ "yes")) %>%
  select(gene, biotype.x, Name.x, GKpelota, ski2.5, hen2.5, wt, rrp4, rrp45a, SAILpelota, sig)

left_joined <- left_join(x = gathered_DESeqresults_rrp45a, y = miRNAtargets, by = "gene", copy = FALSE)

rrp45a_test <- left_joined[left_joined$genotype=="rrp45a",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(!is.na(miRNA)) 

mir_parse <- function(mir_names, simplify=T) {
  result <- lapply(mir_names, function(mir_name) {
    mir <- mir_name
    mir_base <- gsub("([a-zA-z]*[1-9]*)[a-z]*", "\\1", mir)
    mir_letters <- gsub("[a-zA-z]*[1-9]*([A-Za-z]*)", "\\1", mir) %>% str_split("", simplify = T)
    mir_list <- paste0(mir_base, c("", mir_letters))
    return(mir_list)
  })
  if(simplify) {result <- unlist(result)}
  return(result)
}

miRNA_genes <- mir_parse(rrp45a_test$miRNA, simplify = T) %>% 
  toupper()

to_merge <- subset(sign_miRNAs, sign_miRNAs$Name.x %in% miRNA_genes) %>% 
  mutate(sign_target = "yes")

to_plot <- left_join(sign_miRNAs, to_merge) %>% 
  mutate(sign_target = case_when(is.na(sign_target) ~ "no", 
                                 T ~ "yes"))
to_plot3 <- to_plot %>% 
  mutate(both = case_when(sign_target == "yes" & sig == "no" ~ "target.miRNA", 
                          sign_target == "no" & sig == "yes" ~ "diff.exp.miRNA", 
                          sign_target == "no" & sig == "no" ~ "miRNA", 
                          sign_target == "yes" & sig == "yes" ~ "diff.exp.miRNA.target")) %>% 
  glimpse()

## this one is correct
rrp45a_miRNAs <- to_plot3 %>% 
  #mutate(both = factor(both, c("diff.exp.miRNA.target", "diff.exp.miRNA", "target.miRNA", "miRNA"))) %>% 
  #mutate(both = factor(both, c("miRNA", "target.miRNA", "diff.exp.miRNA", "diff.exp.miRNA.target"))) %>%
  ggplot(., aes(x=log2(wt), y=log2(rrp45a), color = both)) +
  #ggplot(.,aes(x=wt, y=rrp45b, color = fct_relevel(both, c("diff.exp.miRNA.target", "diff.exp.miRNA", "target.miRNA", "miRNA")))) +
  #ggplot(.,aes(x=wt, y=rrp45b, color = fct_relevel(both, c("miRNA", "target.miRNA", "diff.exp.miRNA", "diff.exp.miRNA.target")))) +
  geom_point() +
  scale_x_continuous(limits = c(0.01, NA)) +
  scale_y_continuous(limits = c(0.01, NA)) +
  #scale_y_log10(limits = c(0.01,NA), breaks = c(0, 1, 10, 100, 1000)) +
  #scale_x_log10(limits = c(0.01,NA), breaks = c(0, 1, 10, 100, 1000)) +
  coord_equal() +
  #geom_hline(yintercept = 1, colour="#990000", linetype="dashed") + 
  #geom_vline(xintercept = 1, colour="#990000", linetype="dashed") +
  geom_abline(intercept = 0, slope = 1) +
  cowplot::theme_cowplot() +
  geom_text_repel(aes(label=ifelse(both == "diff.exp.miRNA.target",as.character(Name.x),'')),hjust=0,vjust=0) +
  geom_text_repel(aes(label=ifelse(both == "target.miRNA",as.character(Name.x),'')),hjust=0,vjust=0) +
  scale_colour_manual(values=c("diff.exp.miRNA"="#aa4a30","target.miRNA"="#484018", "diff.exp.miRNA.target" = "#e40017", "miRNA" = "#f4c983")) +
  #geom_text_repel(data=subset(to_plot3, rrp45b > 50 & wt > 50), aes(label=Name.x)) +
  ggtitle(label = "rrp45a", subtitle = "scatterplot of miRNA genes in rrp45a vs WT") +
  #scale_color_manual(values = wes_palette(n=4, name="Darjeeling1"))+
  labs(y= "log2(rrp45a[meanRPM])", x = "log2(WT[meanRPM])") +
  theme(aspect.ratio = 1)
ggsave(rrp45a_miRNAs, file = "/binf-isilon/PBgrp/xpj980/figures/20210424miRNAgenesscatterrrp45aaraport.svg", units = "in", width = 15, height = 15, dpi = 300)

#### Figure 4A #####

barplot_table <- left_join(x = gathered_DESeqresults_rrp4, y = miRNAtargets, by = "gene", copy = FALSE) %>%
  filter(padj < 0.05) 

to_plot <- barplot_table %>% 
  mutate(up_or_down = case_when(log2FoldChange > 0 ~ "up", 
                                log2FoldChange < 0 ~ "down"), 
         miRNA = case_when(is.na(miRNA) ~ "other", 
                           TRUE ~ "miRNA_target")) %>%
  select(genotype, miRNA, up_or_down)

to_tally <- to_plot %>% 
  filter(genotype != "rrp45a") %>%
  group_by(genotype, miRNA,up_or_down) %>%
  add_tally(name = "sum") %>%
  distinct() %>%
  ungroup() %>%
  mutate(sum= case_when(up_or_down == "down" ~ -1*sum,
                        up_or_down == "up" ~ 1*sum),
         plot = paste0(miRNA, up_or_down))

to_tally$genotype <- factor(to_tally$genotype, levels = c("ski2.5", "hen2.5", "rrp4", "rrp45a"))
to_tally$plot <- factor(to_tally$plot, levels = c("miRNA_targetup", "miRNA_targetdown", "otherdown", "otherup"))

(barplot_across_geno <- ggplot(to_tally, aes(x= genotype, y=sum, fill = plot)) +
    geom_col(alpha = 0.8, width = 0.5) +
    geom_text(aes(label = sum), vjust=1.6, color="black", size = 3) +
    scale_fill_brewer("more or less sRNAs", 
                      palette = "RdGy",    
                      labels = c("known miRNA target - more sRNAs", "known miRNA target - less sRNAs", "other gene - less sRNAs", "other gene - more sRNAs")) +
    facet_wrap(~miRNA, nrow = 1, scales = "free") +
    labs(y= "number of genes\n", x = "gene type\n") + 
    cowplot::theme_cowplot() + 
    scale_x_discrete(position = "bottom") + 
    ggtitle("Genes with more or less sRNAs compared to WT") + 
    theme(strip.text.x = element_text(size = 12,
                                      face = "bold",
                                      color = "#22292F"),
          strip.background = element_blank(),
          plot.title = element_text(size = 16, 
                                    face = "bold",
                                    color = "#22292F", 
                                    margin = margin(b = 8)), 
          axis.title.x = element_text(size = 10,
                                      face = "bold",
                                      color = "#22292F", 
                                      margin = margin(t = 15)),
          axis.title.y = element_text(size = 10,
                                      face = "bold",
                                      color = "#22292F"), 
          axis.text.x = element_text(size = 8), 
          axis.text.y = element_text(size = 8),
          legend.title = element_text(size = 9,
                                      face = "bold",
                                      color = "#22292F"), 
          legend.text = element_text(size = 8))) + 
  theme(aspect.ratio = 1)

#ggsave(barplot_across_geno, file = "/binf-isilon/PBgrp/xpj980/figures/20210325Figure4Aaraport.svg", device = "svg", units = "in", width = 15, height = 15, dpi = 300)

###### Figure 4B ##### 
############ count other biotypes ##########
araport_TE_IDs <- import.gff("/binf-isilon/PBgrp/xpj980/TAIR/Araport11_GFF3_genes_transposons.201606.gff")

TE_IDs <- araport_TE_IDs %>% 
  as_tibble() %>% 
  filter(type %in% c("transposable_element", "pseudogene", "transposable_element_gene")) %>% 
  select(type, ID, locus_type) %>% 
  mutate(locus_type = case_when(is.na(locus_type) ~ "transposable_element", 
                                T ~ locus_type)) %>% 
  rename(ID = "gene")


# edit type column from whole shabang to original 
araport_original <- import.gff("/binf-isilon/PBgrp/xpj980/TAIR/Araport11_GFF3_genes_transposons.201606.gff") %>% 
  as.tibble() %>% 
  select(type, ID, locus_type) %>% 
  filter(type =="gene") %>% 
  rename(ID = "gene")

to_join <- bind_rows(TE_IDs, araport_original)

test <- left_join(gathered_DESeqresults_rrp4, to_join, by = "gene") %>% 
  mutate(locus_type = case_when(is.na(locus_type) ~ "intergenic", 
                                T ~ locus_type))

to_plot <- left_join(x = test, y = miRNAtargets) %>%
  filter(padj < 0.05) %>% 
  mutate(up_or_down = case_when(log2FoldChange > 0 ~ "up", 
                                log2FoldChange < 0 ~ "down")) %>%
  mutate(miRNA = case_when(is.na(miRNA) ~ locus_type, 
                           T ~ "miRNAtarget"), 
         miRNA = case_when(miRNA == "antisense_long_noncoding_rna" ~ "long_noncoding_rna", 
                           T ~ miRNA), 
         miRNA = case_when(miRNA == "novel_transcribed_region" ~ "transposable_element", 
                           T ~ miRNA), 
         miRNA = case_when(miRNA == "pseudogene" ~ "transposable_element",
                           T ~ miRNA), 
         miRNA = case_when(miRNA == "transposable_element_gene" ~ "transposable_element",
                           T ~ miRNA), 
         miRNA = case_when(miRNA == "small_nuclear_rna" ~ "snRNA",
                           T ~ miRNA), 
         miRNA = case_when(miRNA == "small_nucleolar_rna" ~ "snRNA",
                           T ~ miRNA)) %>% 
  filter(genotype != "rrp45a",
         miRNA != "antisense_rna", 
         miRNA != "other_rna") %>%
  glimpse()

unique(to_plot$miRNA)

to_plot$genotype = factor(to_plot$genotype, levels = c("ski2.5", "hen2.5", "rrp4"))
to_plot$miRNA= factor(to_plot$miRNA, levels = c("miRNAtarget", "mirna", "long_noncoding_rna", "protein_coding", "ribosomal_rna", "pre_trna", "snRNA",  "transposable_element", "intergenic"))

(loveplot <- to_plot %>% 
    ggplot(aes(x=miRNA, y=log2FoldChange, col=miRNA)) +
    geom_point(shape=21,
               col="black",
               size = 0.5,
               fill="transparent",
               position=position_dodge2(width=0.5),
               alpha=.5) +
    scale_y_continuous(position = "right") +
    geom_boxplot(fill="transparent", size=0.5, outlier.shape=1, outlier.size=1, width=0.5, position = position_dodge(), coef=10) +
    geom_hline(yintercept = 0, linetype = 2) +
    facet_wrap(~genotype, nrow = 3, strip.position="left") +
    #coord_flip() +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
    scale_color_manual(values = c(wes_palette("Darjeeling1", n = 9, type = ("continuous")))) 
  # theme(aspect.ratio = 1)
)

ggsave(loveplot, file = "/binf-isilon/PBgrp/xpj980/figures/Figure4boxplotwithTEandintergenicNSincluded.svg", units = "in", width = 10, height = 10, dpi = 300)


### barplot version ### 
to_plot %>% 
  group_by(genotype, miRNA, up_or_down) %>% 
  add_tally() %>%
  filter(up_or_down == "up") %>%
  ggplot(aes(x= miRNA, y = n, fill = miRNA)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~genotype, scales = "free_y") +
  xlab("miRNA targets vs other locustypes") +
  ylab("sum") +
  cowplot::theme_cowplot() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c(wes_palette("Darjeeling1", n = 9, type = ("continuous")))) 
ggsave(loveplot, file = "/binf-isilon/PBgrp/xpj980/figures/Figure4barplotversion.svg", units = "in", width = 10, height = 10, dpi = 300)


##### Figure 4C #####
########## try with sizes ########
# next time just load this one 
load("/binf-isilon/PBgrp/xpj980/R_scripts_collected/2019exosomedata1_sizes_normalizedreads.RData")

genes <- all_annotations %>% 
  as_tibble() %>% 
  dplyr::select(biotype, gene=gene_id, Name) %>% 
  filter(!is.na(gene))

#### sizes on genes from Deseq analysis 
# edit type column from whole shabang to original 
araport_original <- import.gff("/binf-isilon/PBgrp/xpj980/TAIR/Araport11_GFF3_genes_transposons.201606.gff") %>% 
  as.tibble() %>% 
  select(type, ID, locus_type) %>% 
  filter(type =="gene") %>% 
  rename(ID = "gene")

to_join <- bind_rows(TE_IDs, araport_original)

test <- left_join(gathered_DESeqresults_rrp4, to_join, by = "gene") %>% 
  mutate(locus_type = case_when(is.na(locus_type) ~ "intergenic", 
                                T ~ locus_type))

all_genes_table <- left_join(x = test, y = miRNAtargets, by = "gene", copy = FALSE) %>%
  filter(padj < 0.05) %>%
  glimpse()

all_genes_table1 <- left_join(x = all_genes_table, y = genes, by = "gene", copy = FALSE) %>% 
  mutate(up_or_down = case_when(log2FoldChange > 0 ~ "up", 
                                log2FoldChange < 0 ~ "down"), 
         locustype = case_when(is.na(miRNA) ~ locus_type, 
                               TRUE ~ "known miRNA target")) %>%
  glimpse()


test <- filter(normalized, gene %in% all_genes_table1$gene) 
test_errorbar <- filter(normalized_for_errorbars, gene %in% all_genes_table1$gene)

all_genes_table1$up_or_down = factor(all_genes_table1$up_or_down, levels = c("up", "down"))


### new one ###
#### now the mean_RPM of sizes of sRNA is grouped by genotype also ######
to_level <- filter(all_genes_table, genotype == "ski2.5" & log2FoldChange > 0) 

ski2_genes_meanRPM <- test_errorbar %>% 
  filter(., gene %in% c(to_level$gene)) %>% 
  left_join(., genes) %>% 
  left_join(., to_level) %>% 
  mutate(locustype = case_when(is.na(miRNA) ~ locus_type, 
                               TRUE ~ "known miRNA target")) %>%
  filter(genotyp == "ski2.5") %>%
  group_by(genotyp, size, locustype, rep) %>% 
  summarise(mean_RPM = mean(RPM))

to_level <- filter(all_genes_table, genotype == "rrp4" & log2FoldChange > 0) 

rrp4_genes_meanRPM <- test_errorbar %>% 
  filter(., gene %in% c(to_level$gene)) %>% 
  left_join(., genes) %>% 
  left_join(., to_level) %>% 
  mutate(locustype = case_when(is.na(miRNA) ~ locus_type, 
                               TRUE ~ "known miRNA target")) %>%
  filter(genotyp == "rrp4") %>%
  group_by(genotyp, size, locustype, rep) %>% 
  summarise(mean_RPM = mean(RPM))

to_level <- filter(all_genes_table, genotype == "hen2.5" & log2FoldChange > 0) 

hen2_genes_meanRPM <- test_errorbar %>% 
  filter(., gene %in% c(to_level$gene)) %>% 
  left_join(., genes) %>% 
  left_join(., to_level) %>% 
  mutate(locustype = case_when(is.na(miRNA) ~ locus_type, 
                               TRUE ~ "known miRNA target")) %>%
  filter(genotyp == "hen2.5") %>%
  group_by(genotyp, size, locustype, rep) %>% 
  summarise(mean_RPM = mean(RPM))

to_plot <- bind_rows(ski2_genes_meanRPM, rrp4_genes_meanRPM, hen2_genes_meanRPM) %>% 
  filter(locustype == "known miRNA target" | locustype == "protein_coding" | locustype == "intergenic" |locustype == "transposable_element")

to_plot$genotyp = factor(to_plot$genotyp, levels = c("ski2.5", "hen2.5", "rrp4" ))
to_plot$locustype =factor(to_plot$locustype, levels = c("known miRNA target","protein_coding","transposable_element", "intergenic"))


solution4_errorbars <- ggplot(to_plot, aes(x = size, y = mean_RPM, fill = genotyp)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  facet_grid(genotyp ~locustype, scales = "free_y") + 
  scale_fill_manual(values = wes_palette("Darjeeling2", n = 3)) +
  cowplot::theme_cowplot() + 
  theme(aspect.ratio = 1)
#ggsave(solution4_errorbars, file = "/binf-isilon/PBgrp/xpj980/figures/figure4Dsizesaraport_correctothergrid.svg", units = "in", width = 15, height = 10, dpi = 300)

#### Figure 5A ####
geneset <- gathered_DESeqresults_rrp4 %>% 
  semi_join(miRNAtargets, by = "gene") %>% 
  filter(padj < 0.05) %>% 
  pull(gene) %>% 
  unique()

log2toplot <- gathered_DESeqresults_rrp4 %>%
  select(gene, log2FoldChange, genotype) %>% 
  filter(gene %in% geneset) %>%
  spread(genotype, log2FoldChange) %>%
  column_to_rownames("gene")

join_with <- rownames_to_column(log2toplot, var = "gene") 

norm_counts_matrix <- as.data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  select(1:4,8:10,14:16,18:20) %>% 
  as_tibble() %>% 
  semi_join(y = join_with, by = "gene", copy = FALSE) %>% 
  mutate(wt_avg = ((wt_R1 + wt_R2 + wt_R3)/3), 
         ski2.5_avg = ((ski2.5_R1 + ski2.5_R2 + ski2.5_R3)/3), 
         hen2.5_avg = ((hen2.5_R1 + hen2.5_R2 + hen2.5_R3)/3),
         rrp4_avg = ((rrp4_R1 + rrp4_R2 + rrp4_R3)/3)) %>% 
  select(gene, wt_avg, ski2.5_avg, hen2.5_avg, rrp4_avg) %>% 
  column_to_rownames(var = "gene")
## this is just a suggestion to Peter, if he likes it then work more on it. If not go back to orignal by deleting annotation_row argument 

annot_cols <- norm_counts_matrix %>% 
  rownames_to_column(var = "gene") %>% 
  left_join(., miRNAtargets) %>% 
  select(gene, miRNA) %>% 
  column_to_rownames(var = "gene")

normalized_heatmap_rowclust <- pheatmap(norm_counts_matrix,
                                        #annotation_row = annot_cols,
                                        #annotation_colors = annot_colors, 
                                        breaks = c(-3:3),
                                        kmeans_k = NA, 
                                        cluster_rows = 1, 
                                        cluster_cols = 1,
                                        cellwidth = 25, 
                                        cellheight = 7,
                                        fontsize_row = 8,
                                        scale = "row",
                                        #cutree_rows = 2,
                                        color = rev(RColorBrewer::brewer.pal(name = "RdBu", n=6)),
                                        main = "significant DE miRNA targets RPM")
ggsave(normalized_heatmap_rowclust, file = "/binf-isilon/PBgrp/xpj980/figures/heatmapsuggestion.svg", units = "in", width = 15, height = 15, dpi = 300)

### try and filter this one also 
geneset_filtered <- gathered_DESeqresults_rrp4 %>% 
  semi_join(miRNAtargets, by = "gene") %>% 
  filter(padj < 0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  pull(gene) %>% 
  unique()

log2toplot_filtered <- gathered_DESeqresults_rrp4 %>%
  select(gene, log2FoldChange, genotype) %>% 
  filter(gene %in% geneset_filtered) %>%
  spread(genotype, log2FoldChange) %>%
  column_to_rownames("gene")

join_with <- rownames_to_column(log2toplot_filtered, var = "gene") 

norm_counts_matrix <- as.data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  select(1:4,8:10,14:16,18:20) %>% 
  as_tibble() %>% 
  semi_join(y = join_with, by = "gene", copy = FALSE) %>% 
  mutate(wt_avg = ((wt_R1 + wt_R2 + wt_R3)/3), 
         ski2.5_avg = ((ski2.5_R1 + ski2.5_R2 + ski2.5_R3)/3), 
         hen2.5_avg = ((hen2.5_R1 + hen2.5_R2 + hen2.5_R3)/3),
         rrp4_avg = ((rrp4_R1 + rrp4_R2 + rrp4_R3)/3)) %>% 
  select(gene, wt_avg, ski2.5_avg, hen2.5_avg, rrp4_avg) %>% 
  column_to_rownames(var = "gene")

annot_cols <- norm_counts_matrix %>% 
  rownames_to_column(var = "gene") %>% 
  left_join(., miRNAtargets) %>% 
  select(gene, miRNA) %>% 
  column_to_rownames(var = "gene")

normalized_heatmap_rowclust <- pheatmap(norm_counts_matrix,
                                        annotation_row = annot_cols,
                                        #annotation_colors = annot_colors, 
                                        breaks = c(-3:3),
                                        kmeans_k = NA, 
                                        cluster_rows = 1, 
                                        cluster_cols = 0,
                                        cellwidth = 25, 
                                        cellheight = 7,
                                        fontsize_row = 8,
                                        scale = "row",
                                        color = rev(RColorBrewer::brewer.pal(name = "RdBu", n=6)),
                                        main = "significant DE miRNA targets RPM")
#ggsave(normalized_heatmap_rowclust, file = "/binf-isilon/PBgrp/xpj980/figures/20210419normheatrrp42019dataFILTEREDaraport.svg", units = "in", width = 15, height = 15, dpi = 300)

#### Figure 5B #####

left_joined <- left_join(x = gathered_DESeqresults_rrp4, y = miRNAtargets, by = "gene")

ski2_list <- left_joined[left_joined$genotype=="ski2.5",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(!is.na(miRNA)) %>%
  pull(gene) 

hen2_list <- left_joined[left_joined$genotype=="hen2.5",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(!is.na(miRNA)) %>%
  pull(gene) 

rrp4_list <- left_joined[left_joined$genotype=="rrp4",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(!is.na(miRNA)) %>%
  pull(gene) 

#### suppl. fig 6 ####
###### Figure 6A ######
barplot_table <- left_join(x = gathered_DESeqresults_rrp4, y = miRNAtargets, by = "gene", copy = FALSE) %>%
  filter(padj < 0.05) 

to_plot <- barplot_table %>% 
  mutate(up_or_down = case_when(log2FoldChange > 0 ~ "up", 
                                log2FoldChange < 0 ~ "down"), 
         miRNA = case_when(is.na(miRNA) ~ "other", 
                           TRUE ~ "miRNA_target")) %>%
  select(genotype, miRNA, up_or_down)

to_tally <- to_plot %>% 
  #filter(genotype != "rrp45a") %>%
  group_by(genotype, miRNA,up_or_down) %>%
  add_tally(name = "sum") %>%
  distinct() %>%
  ungroup() %>%
  mutate(sum= case_when(up_or_down == "down" ~ -1*sum,
                        up_or_down == "up" ~ 1*sum),
         plot = paste0(miRNA, up_or_down))

to_tally$genotype <- factor(to_tally$genotype, levels = c("ski2.5", "hen2.5", "rrp4", "rrp45a"))
to_tally$plot <- factor(to_tally$plot, levels = c("miRNA_targetup", "miRNA_targetdown", "otherdown", "otherup"))

(barplot_across_geno <- ggplot(to_tally, aes(x= genotype, y=sum, fill = plot)) +
    geom_col(alpha = 0.8, width = 0.5) +
    geom_text(aes(label = sum), vjust=1.6, color="black", size = 3) +
    scale_fill_brewer("more or less sRNAs", 
                      palette = "RdGy",    
                      labels = c("known miRNA target - more sRNAs", "known miRNA target - less sRNAs", "other gene - less sRNAs", "other gene - more sRNAs")) +
    facet_wrap(~miRNA, nrow = 1, scales = "free") +
    labs(y= "number of genes\n", x = "gene type\n") + 
    cowplot::theme_cowplot() + 
    scale_x_discrete(position = "bottom") + 
    ggtitle("Genes with more or less sRNAs compared to WT") + 
    theme(strip.text.x = element_text(size = 12,
                                      face = "bold",
                                      color = "#22292F"),
          strip.background = element_blank(),
          plot.title = element_text(size = 16, 
                                    face = "bold",
                                    color = "#22292F", 
                                    margin = margin(b = 8)), 
          axis.title.x = element_text(size = 10,
                                      face = "bold",
                                      color = "#22292F", 
                                      margin = margin(t = 15)),
          axis.title.y = element_text(size = 10,
                                      face = "bold",
                                      color = "#22292F"), 
          axis.text.x = element_text(size = 8), 
          axis.text.y = element_text(size = 8),
          legend.title = element_text(size = 9,
                                      face = "bold",
                                      color = "#22292F"), 
          legend.text = element_text(size = 8))) + 
  theme(aspect.ratio = 1)
ggsave(barplot_across_geno, file = "/binf-isilon/PBgrp/xpj980/figures/rrp45a.svg", units = "in", width = 15, height = 15, dpi = 300)

#### Supplemental Figure #### 
#ski2-5 overlaps 

# load 2015 dataset and also sop1 dataset 

ski2_setB <- gathered_DESeqresults %>% 
  filter(genotype == "ski2.5") %>% 
  filter(padj < 0.05) %>% 
  left_join(x = ., y = miRNAtargets, by = "gene", copy = FALSE) %>% 
  mutate(up_or_down = case_when(log2FoldChange > 0 ~ "up", 
                                log2FoldChange < 0 ~ "down"), 
         miRNA = case_when(is.na(miRNA) ~ "other", 
                           TRUE ~ "miRNA_target")) 

load(file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/2015_gatheredDESeqresults.RData")

ski2_setA <- gathered_DESeqresults %>% 
  filter(genotype == "ski2.5") %>% 
  filter(padj < 0.05) %>% 
  left_join(x = ., y = miRNAtargets, by = "gene", copy = FALSE) %>% 
  mutate(up_or_down = case_when(log2FoldChange > 0 ~ "up", 
                                log2FoldChange < 0 ~ "down"), 
         miRNA = case_when(is.na(miRNA) ~ "other", 
                           TRUE ~ "miRNA_target")) 

look <- filter(ski2_setA, up_or_down == "up" & miRNA == "miRNA_target")

# what if I had filtered on certain log2FoldChange instead of padj for sample with one replicate 

ski2_setA_filtered <- gathered_DESeqresults %>%
  filter(genotype == "ski2.5") %>% 
  #filter(padj < 0.05) %>%
  filter(log2FoldChange > 1) %>%
  left_join(x = ., y = miRNAtargets, by = "gene", copy = FALSE) %>% 
  mutate(up_or_down = case_when(log2FoldChange > 0 ~ "up", 
                                log2FoldChange < 0 ~ "down"), 
         miRNA = case_when(is.na(miRNA) ~ "other", 
                           TRUE ~ "miRNA_target")) %>%
  filter(miRNA == "miRNA_target")

miRNAs_logFCfiltered <- left_join(ski2_setA_filtered, miRNAtargets, by = "gene")

load(file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/2019exovol2DESEQtable.RData")

ski2_setC <- gathered_DESeqresults %>% 
  filter(genotype == "ski2.5") %>% 
  filter(padj < 0.05) %>% 
  left_join(x = ., y = miRNAtargets, by = "gene", copy = FALSE) %>% 
  mutate(up_or_down = case_when(log2FoldChange > 0 ~ "up", 
                                log2FoldChange < 0 ~ "down"), 
         miRNA = case_when(is.na(miRNA) ~ "other", 
                           TRUE ~ "miRNA_target")) 

filter(ski2_setC, up_or_down == "up" & miRNA == "miRNA_target")

ski2_setA_miRNAtargets_up <- ski2_setA_filtered %>% filter(up_or_down == "up" & miRNA == "miRNA_target") %>% pull(gene)
ski2_setB_miRNAtargets_up <- ski2_setB %>% filter(up_or_down == "up" & miRNA == "miRNA_target") %>% pull(gene)
ski2_setC_miRNAtargets_up <- ski2_setC %>% filter(up_or_down == "up" & miRNA == "miRNA_target") %>% pull(gene)

#ski2_setA_miRNAtargets_down <- ski2_setA %>% filter(up_or_down == "down" & miRNA == "miRNA_target") %>% pull(gene)
ski2_setB_miRNAtargets_down <- ski2_setB %>% filter(up_or_down == "down" & miRNA == "miRNA_target") %>% pull(gene)
ski2_setC_miRNAtargets_down <- ski2_setC %>% filter(up_or_down == "down" & miRNA == "miRNA_target") %>% pull(gene)

#ski2_setA_other_up <- ski2_setA %>% filter(up_or_down == "up" & miRNA == "other") %>% pull(gene)
ski2_setB_other_up <- ski2_setB %>% filter(up_or_down == "up" & miRNA == "other") %>% pull(gene)
ski2_setC_other_up <- ski2_setC %>% filter(up_or_down == "up" & miRNA == "other") %>% pull(gene)

#ski2_setA_other_down <- ski2_setA %>% filter(up_or_down == "down" & miRNA == "other") %>% pull(gene)
ski2_setB_other_down <- ski2_setB %>% filter(up_or_down == "down" & miRNA == "other") %>% pull(gene)
ski2_setC_other_down <- ski2_setC %>% filter(up_or_down == "down" & miRNA == "other") %>% pull(gene)

colors <- c("#FFA6D5", "#8E0505", "#9CC094")

# hits overlap
s4 <- list(setA = c(ski2_setA_miRNAtargets_up),
           setB = c(ski2_setB_miRNAtargets_up), 
           setC = c(ski2_setC_miRNAtargets_up))

hits <- plot(euler(s4, shape = "ellipse"), quantities = TRUE, fill = colors)
ggsave(hits, filename = "/binf-isilon/PBgrp/xpj980/figures/sets_euler_miRNAtargetuplog2FC09.svg")

colors <- c("#8E0505", "#9CC094")
s4 <- list(setB = c(ski2_setB_miRNAtargets_down), 
           setC = c(ski2_setC_miRNAtargets_down))

hits <- plot(euler(s4, shape = "ellipse"), quantities = TRUE, fill = colors)
ggsave(hits, filename = "/binf-isilon/PBgrp/xpj980/figures/sets_euler_miRNAtargetdownfilt.svg")

s4 <- list(
  setB = c(ski2_setB_other_up), 
  setC = c(ski2_setC_other_up))

hits <- plot(euler(s4, shape = "ellipse"), quantities = TRUE, fill = colors)
ggsave(hits, filename = "/binf-isilon/PBgrp/xpj980/figures/sets_euler_otherupfilt.svg")

s4 <- list(
  setB = c(ski2_setB_other_down), 
  setC = c(ski2_setC_other_down))

hits <- plot(euler(s4, shape = "ellipse"), quantities = TRUE, fill = colors)
ggsave(hits, filename = "/binf-isilon/PBgrp/xpj980/figures/sets_euler_otherdownfilt.svg")

# overlap with Anjas paper 

ski2.4_miRNAtargets <- c("AT1G62670", "AT1G62930", "AT1G62910", "AT1G63150", "AT1G63080", "AT1G63130", "AT1G63400", "AT1G63330", "AT1G62590", "AT1G63070", "AT1G12775", "AT1G12620", "AT5G43270", "AT2G28350", "AT1G77850", "AT5G39610", "AT2G34710", "AT1G30490", "AT5G60690", "AT1G30330", "AT1G48410", "AT2G45160", "AT4G36920", "AT3G15030", "AT4G18390", "AT1G53230", "AT3G23690", "AT1G12820", "AT2G33770", "AT1G12290", "AT5G43740", "AT1G02860", "AT3G21170")

# all four ski2 samples
s4 <- list(setA = c(ski2_setA_miRNAtargets_up),
           setB = c(ski2_setB_miRNAtargets_up), 
           setC = c(ski2_setC_miRNAtargets_up), 
           setD = c(ski2.4_miRNAtargets))
colors <- c("#FFA6D5", "#8E0505", "#9CC094", "#7F7C82")

hits <- plot(euler(s4, shape = "ellipse"), quantities = TRUE, fill = colors)
ggsave(hits, filename = "/binf-isilon/PBgrp/xpj980/figures/sets_euler_miRNAtargetuplog2FCAnjas.svg")


#### qPCR ####

FL_qPCR <- readxl::read_xlsx("/binf-isilon/PBgrp/xpj980/qPCR_tables/NAR_FLResultssaved.xlsx", col_names = T) %>% 
  select(Well, Cq, target, genotype) %>% 
  filter(target != "empty")  

# remove WT outliers
FL_qPCR %>% filter(genotype == "WT")

FL_qPCR <- filter(FL_qPCR, Well != "B06") %>% 
  filter(., Well != "C11")

control <- filter(FL_qPCR, target == "ACT2") %>% 
  select(genotype, Cq) %>% 
  rename(Cq = "CTRLcq") %>%
  group_by(genotype) %>%
  summarise(CTRLcq = mean(CTRLcq)) %>% 
  mutate(CF = rep("FL"))

deltaCq <- FL_qPCR %>% 
  filter(target != "ACT2") %>%
  left_join(., control) %>% 
  mutate(delta = Cq - CTRLcq, 
         CF = rep("FL"))

average_WT <- deltaCq %>% 
  filter(genotype == "WT") %>% 
  group_by(target) %>% 
  summarise(avg_delta = mean(delta))

to_plot <- left_join(deltaCq, average_WT) %>% 
  mutate(deltadelta = delta - avg_delta) %>%
  mutate(deltadeltapow = 2^(-deltadelta))

FL_toplot <- to_plot  

plot <- FL_toplot %>% 
  ggplot(., aes(x = genotype, y = log2(deltadeltapow), fill = genotype)) +
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2) +
  geom_point() +
  facet_wrap(~target, nrow = 1) +
  cowplot::theme_cowplot()

# good now add 5CF and 3CF data #
#OBS make individual RSG2 plot for 3CF as it is the only sample with measurable Ct
# 5CF 
CF5_qPCR <- readxl::read_xlsx("/binf-isilon/PBgrp/xpj980/qPCR_tables/NAR_5CF.xlsx", col_names = T) %>% 
  select(Well, Cq, target, genotype) %>% 
  filter(target != "empty")


# remove WT outliers
CF5_qPCR %>% filter(genotype == "WT")
# none to remove

control <- filter(CF5_qPCR, target == "ACT2") %>% 
  select(genotype, Cq) %>% 
  rename(Cq = "CTRLcq") %>%
  group_by(genotype) %>%
  summarise(CTRLcq = mean(CTRLcq)) %>% 
  mutate(CF = rep("5CF"))

deltaCq <- CF5_qPCR %>% 
  filter(target != "ACT2") %>%
  left_join(., control) %>% 
  mutate(delta = Cq - CTRLcq, 
         CF = rep("5CF"))

average_WT <- deltaCq %>% 
  filter(genotype == "WT") %>% 
  group_by(target) %>% 
  summarise(avg_delta = mean(delta))

to_plot <- left_join(deltaCq, average_WT) %>% 
  mutate(deltadelta = delta - avg_delta) %>%
  mutate(deltadeltapow = 2^(-deltadelta))

CF5_toplot <- to_plot 
# filter(genotype != "WT") 

# 3CF 
CF3_qPCR <- readxl::read_xlsx("/binf-isilon/PBgrp/xpj980/qPCR_tables/NAR_3CF.xlsx", col_names = T) %>% 
  select(Well, Cq, target, genotype) %>% 
  filter(target != "empty") 

# remove WT outliers
CF3_qPCR %>% filter(genotype == "WT")

CF3_qPCR <- filter(CF3_qPCR, Well != "C09")

control <- filter(CF3_qPCR, target == "ACT2") %>% 
  select(genotype, Cq) %>% 
  rename(Cq = "CTRLcq") %>%
  group_by(genotype) %>%
  summarise(CTRLcq = mean(CTRLcq)) %>% 
  mutate(CF = rep("3CF"))

deltaCq <- CF3_qPCR %>% 
  filter(target != "ACT2") %>%
  left_join(., control) %>% 
  mutate(delta = Cq - CTRLcq, 
         CF = rep("3CF"))

average_WT <- deltaCq %>% 
  filter(genotype == "WT") %>% 
  group_by(target) %>% 
  summarise(avg_delta = mean(delta))

to_plot <- left_join(deltaCq, average_WT) %>% 
  mutate(deltadelta = delta - avg_delta) %>%
  mutate(deltadeltapow = 2^(-deltadelta))

CF3_toplot <- to_plot 

# merge all three tables 

gathered_qPCR <- rbind(... = FL_toplot, CF5_toplot, CF3_toplot)

gathered_qPCR$genotype <- factor(gathered_qPCR$genotype, levels = c("ski2.5", "ski2.2", "pel1.1", "pel1.2", "WT"))
gathered_qPCR$CF <- factor(gathered_qPCR$CF, levels = c("5CF", "FL", "3CF"))

plot <- gathered_qPCR %>% 
  group_by(genotype, CF, target) %>%
  ggplot(., aes(x = CF, y = log2(deltadeltapow), fill = CF)) +
  stat_summary(geom = "bar", fun = "mean", position = "dodge") + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2, color = "grey") +
  geom_point() +
  ggtitle("facet by genotype") +
  scale_fill_brewer() +
  facet_grid(target~genotype, scales = "free") +
  cowplot::theme_cowplot()
ggsave(plot, filename = "/binf-isilon/PBgrp/xpj980/figures/qPCRFCgeno.svg")

my_pal <- c("#142F43", "#32C1CD", "#E9A6A6","#C85C5C")

plot <- gathered_qPCR %>% 
  filter(genotype != "WT") %>%
  group_by(genotype, CF, target) %>%
  ggplot(., aes(x = genotype, y = log2(deltadeltapow), fill = genotype,  alpha = 0.7)) +
  stat_summary(geom = "bar", fun = "mean", position = "dodge", width = 0.5) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.2, width = 0.5) +
  geom_point(size = 0.4) +
  ggtitle("log2(FC) value and facet by CF") +
  geom_hline(yintercept = 0, size =0.2) +
  scale_fill_manual(values = my_pal) +
  facet_grid(target~CF, scales = "free") +
  cowplot::theme_cowplot()
ggsave(plot, filename = "/binf-isilon/PBgrp/xpj980/figures/qPCRlog2FCCF.svg")

#### Tukey test ####

# test on one target 

PHO2_Tukey <- filter(gathered_qPCR, target == "PHO2" & CF == "FL")

model <- aov(log2(deltadeltapow)~genotype, data=PHO2_Tukey)
summary(model)

TukeyHSD(model, conf.level=.95)

RSG2_Tukey <- filter(gathered_qPCR, target == "RSG2" & CF == "3CF")

model <- aov(log2(deltadeltapow)~genotype, data=RSG2_Tukey)
summary(model)

TukeyHSD(model, conf.level=.95)

SPL2_Tukey <- filter(gathered_qPCR, target == "SPL2" & CF == "FL")

model <- aov(log2(deltadeltapow)~genotype, data=SPL2_Tukey)
summary(model)

TukeyHSD(model, conf.level=.95)

#### seed count ####

seeds_count <- readxl::read_xlsx("/binf-isilon/PBgrp/xpj980/qPCR_tables/counts.xlsx", col_names = T) %>% 
  column_to_rownames(var = "genotype")

seed_frq <- prop.table(data.matrix(seeds_count), 1) %>%
  melt() %>%
  as.data.frame()

colors <- c("#7E370C", "#FFCE45")

seeds <- ggplot(seed_frq, aes(x=Var1, y= value, fill = Var2)) +
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(values = colors)

ggsave(seeds, filename = "/binf-isilon/PBgrp/xpj980/figures/seedcount.svg")

