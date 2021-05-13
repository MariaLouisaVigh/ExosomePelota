##sop1 araport ##
#20-4-21#

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

setwd("/binf-isilon/PBgrp/xpj980/datascratch/2019_exosome_mutants_vol2/04.STAR_mapping/samfiles/")

##import my featureCounts table 
raw_featureCount_table <- read.table("featureCountsaraport", header = TRUE, row.names = 1)
glimpse(raw_featureCount_table)
## make the table DESeq2-usable by removing useless columns (chr, start, end, strand and length) and also hen2.sop3
DESeq_table <- raw_featureCount_table %>%
  rownames_to_column() %>%
  select(1,7:9,13:15,22:27) %>% 
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

names <- c("wt_R1", "wt_R2", "wt_R3", "hen2.5_R1", "hen2.5_R2", "hen2.5_R3", "ski2.5_R1", "ski2.5_R2", "ski2.5_R3", "sop1.5_R1", "sop1.5_R2", "sop1.5_R3")
names(DESeq_table_summed) <- names
glimpse(DESeq_table_summed)

## make metadata matrix
introduce_genotype_to_table  <- data.frame(row.names=colnames(DESeq_table_summed),
                                           "genotype"=c(rep("wt", 3), 
                                                        rep("hen2.5", 3),
                                                        rep("ski2.5", 3), 
                                                        rep("sop1.5", 3)))

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
#ggsave(hierarchial_heat2, file="/binf-isilon/PBgrp/xpj980/figures/20210420hierachialheatsop1araport.svg", units = "in", width = 15, height = 15, dpi = 300)

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

#ggsave(PCAplot, file = "/binf-isilon/PBgrp/xpj980/figures/20210420PCAsop1araport.svg", units = "in", width = 15, height = 15, dpi = 300)

sampleDists <- dist(t(assay(log_transformeddds2)))

sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
test <- pheatmap(sampleDistMatrix,
                 clustering_distance_rows=sampleDists,
                 clustering_distance_cols=sampleDists,
                 col=colors)

#ggsave(test, file = "/binf-isilon/PBgrp/xpj980/figures/20210420distancematsop1araport.svg", units = "in", width = 15, height = 15, dpi = 300)


########## DESeq Analysis #########
deseqqed <- DESeq(dds2)

resultsNames(deseqqed)

## extract results 
DESeqresultsSki2.5 <- DESeq2::results(deseqqed, contrast = c("genotype", "ski2.5", "wt"), alpha = 0.1)
DESeqresultsHen2.5 <- DESeq2::results(deseqqed, contrast = c("genotype", "hen2.5", "wt"), alpha = 0.1)
DESeqresultsSop1.5 <- DESeq2::results(deseqqed, contrast = c("genotype", "sop1.5", "wt"), alpha = 0.1)

###### gathered DESegresults-table #######

ski2.5_tibble <- DESeqresultsSki2.5 %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  mutate(genotype = rep("ski2.5")) %>%
  glimpse()

hen2.5_tibble <- DESeqresultsHen2.5 %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  mutate(genotype = rep("hen2.5")) %>%
  glimpse()

sop1_tibble <- DESeqresultsSop1.5 %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  mutate(genotype = rep("sop1.5")) %>%
  glimpse()

gathered_DESeqresults <- bind_rows(... = ski2.5_tibble, hen2.5_tibble, sop1_tibble) %>%
  glimpse()

write.xlsx(x = gathered_DESeqresults, file = "/binf-isilon/PBgrp/xpj980/datascratch/GEO/deseqcontrast3.xlsx")
save(gathered_DESeqresults, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/2019exovol2DESEQtable.RData")
load(file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/2019exovol2DESEQtable.RData")


#### excel sheet for paper without intergenic regions### 
test <- left_join(gathered_DESeqresults, to_join, by = "gene") %>% 
  mutate(locus_type = case_when(is.na(locus_type) ~ "intergenic", 
                                T ~ locus_type)) %>% 
  filter(padj < 0.05) %>% 
  filter(locus_type !="intergenic") %>% 
  left_join(., miRNAtargets, by = "gene", copy = FALSE) %>% 
  mutate(locus_type = case_when(is.na(miRNA) ~ locus_type, 
                                TRUE ~ "miRNA_target")) %>% 
  select(gene, baseMean, log2FoldChange, pvalue, padj, locus_type, genotype)
write.xlsx(x = test, file = "/binf-isilon/PBgrp/xpj980/figures/DESeq2ExperimentCWOintergenic.xlsx")


##### Suppl. Fig2 #####

# scatterplot miRNA genes ski2-5
all_annotations <- readGFF("/binf-isilon/PBgrp/xpj980/TAIR/Arabidopsis_thaliana.TAIR10.46.gff3")

genes <- all_annotations %>% 
  as_tibble() %>% 
  dplyr::select(biotype, gene=gene_id, Name) %>% 
  filter(!is.na(gene)) 

joined_sign <- left_join(x = gathered_DESeqresults, y = genes) %>%
  filter(padj < 0.05 , biotype == "miRNA") 

## list all miRNAs (miRNA genes) up- or downregulated in genotypes 

a_list <- lapply(joined_sign$genotype %>% unique(), function(x) {
  joined_sign %>% filter(genotype==x) %>% pull(Name) %>% unique() %>% sort()
})
names(a_list) <- joined_sign$genotype %>% unique()
a_list

### scatterplots ###
# ski2-5
joined_sign <- left_join(x = gathered_DESeqresults, y = genes) %>%
  filter(genotype == "ski2.5", padj < 0.05 , biotype == "miRNA") 

to_gather <- DESeq_table_summed %>% 
  rownames_to_column(var = "gene") %>% 
  glimpse()

keycol <- "genotype"
valcol <- "count"
gathercols <- c(names(to_gather[,2:13]))

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
  select(gene, biotype.x, Name.x, ski2.5, hen2.5, wt, sop1.5, sig)

left_joined <- left_join(x = gathered_DESeqresults, y = miRNAtargets, by = "gene", copy = FALSE)

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
ggsave(ski2_miRNAs, file = "/binf-isilon/PBgrp/xpj980/figures/20210424miRNAgenesscatterski2thirdlibaraport.svg", units = "in", width = 15, height = 15, dpi = 300)


# hen2.5
joined_sign <- left_join(x = gathered_DESeqresults, y = genes) %>%
  filter(genotype == "hen2.5", padj < 0.05 , biotype == "miRNA") 

to_gather <- DESeq_table_summed %>% 
  rownames_to_column(var = "gene") %>% 
  glimpse()

keycol <- "genotype"
valcol <- "count"
gathercols <- c(names(to_gather[,2:13]))

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
  select(gene, biotype.x, Name.x, ski2.5, hen2.5, wt, sop1.5, sig)

left_joined <- left_join(x = gathered_DESeqresults, y = miRNAtargets, by = "gene", copy = FALSE)

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
  scale_colour_manual(values=c("diff.exp.miRNA"="#aa4a30","target.miRNA"="#484018", "diff.exp.miRNA.target" = "#e40017", "miRNA" = "#f4c983")) +
  #geom_text_repel(data=subset(to_plot3, rrp45b > 50 & wt > 50), aes(label=Name.x)) +
  ggtitle(label = "hen2-5", subtitle = "scatterplot of miRNA genes in ski2-5 vs WT") +
  #scale_color_manual(values = wes_palette(n=4, name="Darjeeling1"))+
  labs(y= "log2(hen2.5[meanRPM])", x = "log2(WT[meanRPM])") +
  theme(aspect.ratio = 1)
ggsave(hen2_miRNAs, file = "/binf-isilon/PBgrp/xpj980/figures/20210424miRNAgenesscatterhen2thirdlibaraport.svg", units = "in", width = 15, height = 15, dpi = 300)

# sop1
joined_sign <- left_join(x = gathered_DESeqresults, y = genes) %>%
  filter(genotype == "sop1.5", padj < 0.05 , biotype == "miRNA") 

to_gather <- DESeq_table_summed %>% 
  rownames_to_column(var = "gene") %>% 
  glimpse()

keycol <- "genotype"
valcol <- "count"
gathercols <- c(names(to_gather[,2:13]))

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
  select(gene, biotype.x, Name.x, ski2.5, hen2.5, wt, sop1.5, sig)

left_joined <- left_join(x = gathered_DESeqresults, y = miRNAtargets, by = "gene", copy = FALSE)

sop1_test <- left_joined[left_joined$genotype=="sop1.5",] %>% 
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

miRNA_genes <- mir_parse(sop1_test$miRNA, simplify = T) %>% 
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

sop1_miRNAs <- to_plot3 %>% 
  #mutate(both = factor(both, c("diff.exp.miRNA.target", "diff.exp.miRNA", "target.miRNA", "miRNA"))) %>% 
  #mutate(both = factor(both, c("miRNA", "target.miRNA", "diff.exp.miRNA", "diff.exp.miRNA.target"))) %>%
  ggplot(., aes(x=log2(wt), y=log2(sop1.5), color = both)) +
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
  ggtitle(label = "sop1-5", subtitle = "scatterplot of miRNA genes in ski2-5 vs WT") +
  #scale_color_manual(values = wes_palette(n=4, name="Darjeeling1"))+
  labs(y= "log2(sop1.5[meanRPM])", x = "log2(WT[meanRPM])") +
  theme(aspect.ratio = 1)
ggsave(sop1_miRNAs, file = "/binf-isilon/PBgrp/xpj980/figures/20210424miRNAgenesscattersop1thirdlibaraport.svg", units = "in", width = 15, height = 15, dpi = 300)

###### Suppl. Fig3 ######
barplot_table <- left_join(x = gathered_DESeqresults, y = miRNAtargets, by = "gene", copy = FALSE) %>%
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

to_tally$genotype <- factor(to_tally$genotype, levels = c("ski2.5", "hen2.5", "sop1.5"))
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

#ggsave(barplot_across_geno, file = "/binf-isilon/PBgrp/xpj980/figures/20210319Figure3Aaraport.svg", device = "svg", units = "in", width = 15, height = 15, dpi = 300)

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
ggsave(barplot_across_geno, file = "/binf-isilon/PBgrp/xpj980/figures/20210420sop1Aaraport.svg", device = "svg", units = "in", width = 15, height = 15, dpi = 300)
