###### figure 2 #####

###### set wd and load packages ######
setwd("/binf-isilon/PBgrp/xpj980/datascratch/2015_exosome_mutants/03.STAR_mapping/sam_files")

library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)
library(ggrepel)
library(rtracklayer) 
library(VennDiagram)
library(png)
library(systemPipeR)
library(eulerr)
library(GeneOverlap)
library(svglite)
library(reshape2)
library(wesanderson)
library(GenomicRanges)
library(GenomicAlignments)
library(readxl)
library(magrittr)

#### DESeqresults table ####
load("/binf-isilon/PBgrp/xpj980/R_scripts_collected/2015_gatheredDESeqresults.RData")

###### miRNA target list ######## 

miRNAtargets_rawtable <- readxl::read_xlsx("/binf-isilon/PBgrp/xpj980/TAIR/miRNAtargetlist_Maria_degradome_done.xlsx", col_names = T)

miRNAtargets <- miRNAtargets_rawtable %>% 
  select(miRNA, 
         gene = target_id, 
         target_name, 
         miRNA_ID)

miRNAtargets <- filter(miRNAtargets, !duplicated(miRNAtargets$gene))

##import my featureCounts table 
raw_featureCount_table <- read.table("/binf-isilon/PBgrp/xpj980/datascratch/2015_exosome_mutants/03.STAR_mapping/sam_files/featureCountsaraport", header = TRUE, row.names = 1)
glimpse(raw_featureCount_table)

## make the table DESeq2-usable by removing useless columns (chr, start, end, strand and length)
DESeq_table <- raw_featureCount_table %>%
  rownames_to_column() %>%
  select(1,7:12,15:20) %>% 
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

names <- c("rrp45b_R1", "rrp45b_R2", "wt_R1", "wt_R2", "hen2.5_R1", "hen2.5_R2", "rrp4_R1", "rrp4_R2", "ski2.5_R1", "ski2.5_R2", "ski3_R1", "ski3_R2")
names(DESeq_table_summed) <- names
glimpse(DESeq_table_summed)

## reorder columns so that WT replicates are the first two columns 
colorder <- c("wt_R1", "wt_R2", "rrp45b_R1", "rrp45b_R2", "hen2.5_R1", "hen2.5_R2", "rrp4_R1", "rrp4_R2", "ski2.5_R1", "ski2.5_R2", "ski3_R1", "ski3_R2")
reordered_deseq_table <- DESeq_table_summed[,colorder]

######### make deseq design table - remember to take out bad replicate ski2-5 R2 ########

deseq_table <- reordered_deseq_table %>% 
  rownames_to_column(var = "gene") %>% 
  select(1:5,10,12:13) %>% 
  column_to_rownames(var = "gene") %>% 
  glimpse()

## make metadata matrix
introduce_genotype_to_table  <- data.frame(row.names=colnames(deseq_table),
                                           "genotype"=c(rep("wt", 2),
                                                        rep("rrp45b", 2), 
                                                        rep("ski2.5", 1),
                                                        rep("ski3", 2)))


str(introduce_genotype_to_table)
levels(introduce_genotype_to_table)
introduce_genotype_to_table$genotype  <- as.factor(introduce_genotype_to_table$genotype)
introduce_genotype_to_table$genotype <- relevel(introduce_genotype_to_table$genotype, ref="wt")
levels(introduce_genotype_to_table$genotype)

design  <- model.matrix(~genotype, data=introduce_genotype_to_table)

## convert featurecount table into a matrix
DESeq_matrix <- as.matrix(deseq_table)
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

##### average count calculated by myself for later use #####

to_gather <- deseq_table %>% 
  rownames_to_column(var = "gene") 

keycol <- "genotype"
valcol <- "count"
gathercols <- c(names(to_gather[,2:8]))

test <- gather_(to_gather, keycol, valcol, gathercols)

normalized <- separate(test, col = genotype, into = c("genotyp", "rep"), sep = "_") %>%
  group_by(genotyp,rep) %>%
  mutate(RPM = count*1000000/sum(count)) %>%
  group_by(gene, genotyp) %>% 
  summarise(mean_RPM = mean(RPM))

## log_transformation with vst function 
log_transformeddds2 <- vst(dds2, blind = TRUE)

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

#ggsave(PCAplot, file = "/binf-isilon/PBgrp/xpj980/figures/20210302PCAski2araport.svg", units = "in", width = 15, height = 15, dpi = 300)

#### Sup. Fig1 #####
# only run it once because it needs the additional ski2 R2 in it ##

#sampleDists <- dist(t(assay(log_transformeddds2)))

#sampleDistMatrix <- as.matrix(sampleDists)
#colnames(sampleDistMatrix) <- NULL
#colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
#test <- pheatmap(sampleDistMatrix,
         #clustering_distance_rows=sampleDists,
         #clustering_distance_cols=sampleDists,
         #col=colors)

#ggsave(test, file = "/binf-isilon/PBgrp/xpj980/figures/20210324distancematski2araport.png", units = "in", width = 15, height = 15, dpi = 300)


########## DESeq Analysis #########
deseqqed <- DESeq(dds2)


resultsNames(deseqqed)

## plot estimated dispersions 
svg(file = "/binf-isilon/PBgrp/xpj980/figures/20210302dispski2araport.svg")
(dispersionEst <- plotDispEsts(deseqqed))
dev.off()

## extract results 
DESeqresultsRrp45b <- DESeq2::results(deseqqed, contrast = c("genotype", "rrp45b", "wt"), alpha = 0.1)
DESeqresultsSki2.5 <- DESeq2::results(deseqqed, contrast = c("genotype", "ski2.5", "wt"), alpha = 0.1)
DESeqresultsSki3 <- DESeq2::results(deseqqed, contrast = c("genotype", "ski3", "wt"), alpha = 0.1)

###### gathered DESegresults-table #######
rrp45b_tibble <- DESeqresultsRrp45b %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  mutate(genotype = rep("rrp45b")) %>%
  glimpse()

ski2.5_tibble <- DESeqresultsSki2.5 %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  mutate(genotype = rep("ski2.5")) %>%
  glimpse()

ski3_tibble <- DESeqresultsSki3 %>% 
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  mutate(genotype = rep("ski3")) %>%
  glimpse()

gathered_DESeqresults <- bind_rows(... = rrp45b_tibble, ski2.5_tibble, ski3_tibble) 

#write.xlsx(x = gathered_DESeqresults, file = "/binf-isilon/PBgrp/xpj980/datascratch/GEO/deseqcontrast1.xlsx")
#save(gathered_DESeqresults, 
  #   file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/2015_gatheredDESeqresults.RData")

### excel sheet for paper without intergenic regions ## 
test <- left_join(gathered_DESeqresults, to_join, by = "gene") %>% 
  mutate(locus_type = case_when(is.na(locus_type) ~ "intergenic", 
                                T ~ locus_type)) %>% 
  filter(padj < 0.05) %>% 
  filter(locus_type !="intergenic") %>% 
  left_join(., miRNAtargets, by = "gene", copy = FALSE) %>% 
  mutate(locus_type = case_when(is.na(miRNA) ~ locus_type, 
                                TRUE ~ "miRNA_target")) %>% 
  select(gene, baseMean, log2FoldChange, pvalue, padj, locus_type, genotype)
#write.xlsx(x = test, file = "/binf-isilon/PBgrp/xpj980/figures/DESeq2ExperimentAWOintergenic.xlsx")

# Supplementary Table S2 (2015 dataset - need to filter ski2-5 on log2FC instead of padj (revision))

part1 <- left_join(gathered_DESeqresults, to_join, by = "gene") %>% 
  mutate(locus_type = case_when(is.na(locus_type) ~ "intergenic", 
                                T ~ locus_type)) %>% 
  filter(padj < 0.05) %>% 
  filter(locus_type !="intergenic") %>% 
  left_join(., miRNAtargets, by = "gene", copy = FALSE) %>% 
  mutate(locus_type = case_when(is.na(miRNA) ~ locus_type, 
                                TRUE ~ "miRNA_target")) %>% 
  select(gene, baseMean, log2FoldChange, pvalue, padj, locus_type, genotype) %>%
  filter(genotype != "ski2.5")

part2 <- left_join(gathered_DESeqresults, to_join, by = "gene") %>% 
  mutate(locus_type = case_when(is.na(locus_type) ~ "intergenic", 
                                T ~ locus_type)) %>% 
  filter(log2FoldChange > 1) %>% 
  filter(locus_type !="intergenic") %>% 
  left_join(., miRNAtargets, by = "gene", copy = FALSE) %>% 
  mutate(locus_type = case_when(is.na(miRNA) ~ locus_type, 
                                TRUE ~ "miRNA_target")) %>% 
  select(gene, baseMean, log2FoldChange, pvalue, padj, locus_type, genotype) %>% 
  filter(genotype == "ski2.5") %>% 
  filter(locus_type == "miRNA_target")

new_table_for_paper <- rbind(part1, part2)

write.xlsx(x = new_table_for_paper, file = "/binf-isilon/PBgrp/xpj980/figures/DESeq2ExperimentBWOintergenicSKI2update.xlsx")

#### HEATMAPS ####
## log2fold heatmap
geneset <- gathered_DESeqresults %>% 
  semi_join(miRNAtargets, by = "gene") %>% 
  filter(padj < 0.05) %>% 
  pull(gene) %>% 
  unique()

log2toplot <- gathered_DESeqresults %>%
  select(gene, log2FoldChange, genotype) %>% 
  filter(gene %in% geneset) %>%
  spread(genotype, log2FoldChange) %>%
  column_to_rownames("gene")

heatlog2scale <- pheatmap(log2toplot, 
                          breaks = c(-4:4), 
                          color = rev(RColorBrewer::brewer.pal(name = "RdYlBu", n=8)), 
                          main = "significant DE miRNA targets log2FoldChange",
                          cluster_cols = 1,
                          cluster_rows = 1,
                          cutree_rows = 2,
                          cellwidth = 25, 
                          cellheight = 7,
                          fontsize_row = 8)  
#ggsave(heatlog2scale, file = "/binf-isilon/PBgrp/xpj980/figures/20210302heatlog2figure2araport.svg", units = "cm", width = 20, height = 25, dpi = 300)

## normalized counts heatmap

## normalized counts heatmap - they don't work as well as log2foldchange
join_with <- rownames_to_column(log2toplot, var = "gene") 

norm_counts_matrix <- as.data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>% 
  semi_join(y = join_with, by = "gene", copy = FALSE) %>% 
  mutate(wt_avg = ((wt_R1 + wt_R2)/2), 
         ski2.5_avg = ((ski2.5_R1)/1), 
         rrp45b_avg = ((rrp45b_R1 + rrp45b_R2)/2),
         ski3_avg = ((ski3_R1 + ski3_R2)/2)) %>% 
  select(gene, wt_avg, ski2.5_avg, rrp45b_avg, ski3_avg) %>% 
  column_to_rownames(var = "gene")
## this is just a suggestion to Peter, if he likes it then work more on it. If not go back to orignal by deleting annotation_row argument 

annot_cols <- norm_counts_matrix %>% 
  rownames_to_column(var = "gene") %>% 
  left_join(., miRNAtargets) %>% 
  select(gene, miRNA) %>% 
  column_to_rownames(var = "gene")

normalized_heatmap_rowclust <- pheatmap(norm_counts_matrix,
                                        annotation_row = annot_cols,
                                        #annotation_colors = ann_colors, 
                                        breaks = c(-3:3),
                                        kmeans_k = NA, 
                                        cluster_rows = 1, 
                                        cluster_cols = 1,
                                        cellwidth = 25, 
                                        cellheight = 7,
                                        fontsize_row = 8,
                                        scale = "row",
                                        cutree_rows = 2,
                                        color = rev(RColorBrewer::brewer.pal(name = "RdYlBu", n=6)),
                                        main = "significant DE miRNA targets RPM")
#ggsave(normalized_heatmap_rowclust, file = "/binf-isilon/PBgrp/xpj980/figures/20210301heatmapski2Araport.svg", units = "in", width = 15, height = 15, dpi = 300)

######### Figure 2A ################

barplot_table <- left_join(x = gathered_DESeqresults, y = miRNAtargets, by = "gene", copy = FALSE) %>%
  filter(padj < 0.05) 

to_plot <- barplot_table %>% 
  mutate(up_or_down = case_when(log2FoldChange > 0 ~ "up", 
                                log2FoldChange < 0 ~ "down"), 
         miRNA = case_when(is.na(miRNA) ~ "other", 
                           TRUE ~ "known\n miRNA target")) %>%
  dplyr::select(gene, genotype, miRNA, up_or_down)

to_tally <- to_plot %>% 
  group_by(genotype, miRNA,up_or_down) %>%
  add_tally(name = "sum") %>%
  distinct() %>%
  ungroup() %>%
  mutate(sum= case_when(up_or_down == "down" ~ -1*sum,
                        up_or_down == "up" ~ 1*sum), 
         plot = paste0(miRNA, up_or_down))

to_tally$genotype <- factor(to_tally$genotype, levels = c("ski2.5", "ski3", "rrp45b"))
to_tally$plot <- factor(to_tally$plot, levels = c("known\n miRNA targetup", "known\n miRNA targetdown", "otherdown", "otherup"))

for2A <- to_tally %>% 
  select(genotype, miRNA, up_or_down, sum, plot) %>% 
  distinct()

(barplot_across_geno <- ggplot(for2A, aes(x= miRNA, y=sum, fill = plot)) +
    geom_bar(stat= "identity", alpha = 0.8) +
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

(barplot_across_geno <- ggplot(for2A, aes(x= genotype, y=sum, fill = plot)) +
    geom_col(alpha = 0.8, width = 0.5) +
    geom_text(aes(label = sum), vjust=1.6, color="black", size = 3) +
    scale_fill_brewer("more or less sRNAs", 
                      palette = "RdGy",    
                      labels = c("known miRNA target - more sRNAs", "known miRNA target - less sRNAs", "other gene - less sRNAs", "other gene - more sRNAs")) +
    facet_wrap(~miRNA, nrow = 1, scales = "free") +
    labs(y= "number of genes\n", x = "gene type\n") + 
    cowplot::theme_cowplot() + 
    geom_hline(yintercept = 0) +
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

ggsave(barplot_across_geno, file = "/binf-isilon/PBgrp/xpj980/figures/2021042302Figure2Aaraport.svg", device = "svg", units = "in", width = 15, height = 15, dpi = 300)

##### Sup. Fig2 #####
## scatterplot of miRNAs (miRNA genes) ##
all_annotations <- readGFF("/binf-isilon/PBgrp/xpj980/TAIR/Arabidopsis_thaliana.TAIR10.46.gff3")

genes <- all_annotations %>% 
  as_tibble() %>% 
  dplyr::select(biotype, gene=gene_id, Name) %>% 
  filter(!is.na(gene)) 


##### miRNA genes scatterplot ####

joined_sign <- left_join(x = gathered_DESeqresults, y = genes) %>%
  filter(genotype == "rrp45b", padj < 0.05 , biotype == "miRNA") 

miRNAs_RPM_table <- normalized %>% 
  left_join(., genes, by = "gene") %>% 
  filter(., biotype %in% "miRNA") %>% 
  spread(., genotyp, mean_RPM) %>% 
  drop_na(., Name) 

sign_miRNAs <- left_join(miRNAs_RPM_table, joined_sign, by = "gene") %>% 
  mutate(sig = case_when(is.na(Name.y) ~ "no", 
                         TRUE ~ "yes")) %>%
  select(gene, biotype.x, Name.x, rrp45b, ski2.5, ski3, wt, sig)

left_joined <- left_join(x = gathered_DESeqresults, y = miRNAtargets, by = "gene", copy = FALSE)

rrp45b_test <- left_joined[left_joined$genotype=="rrp45b",] %>% 
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

miRNA_genes <- mir_parse(rrp45b_test$miRNA, simplify = T) %>% 
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

rrp45b_miRNAs <- to_plot3 %>% 
  #mutate(both = factor(both, c("diff.exp.miRNA.target", "diff.exp.miRNA", "target.miRNA", "miRNA"))) %>% 
  #mutate(both = factor(both, c("miRNA", "target.miRNA", "diff.exp.miRNA", "diff.exp.miRNA.target"))) %>%
  ggplot(., aes(x=log2(wt), y=log2(rrp45b), color = both)) +
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
  geom_text_repel(aes(label=ifelse(both == "target.miRNA",as.character(Name.x),'')),hjust=0,vjust=0) +
  #geom_text_repel(aes(label=ifelse(both == "diff.exp.miRNA",as.character(Name.x),'')),hjust=0,vjust=0) +
  #geom_text_repel(data=subset(to_plot3, rrp45b > 50 & wt > 50), aes(label=Name.x)) +
  ggtitle(label = "rrp45b", subtitle = "scatterplot of miRNA genes in rrp45b vs WT") +
  #scale_color_manual(values = wes_palette(n=4, name="Darjeeling1"))+
  scale_colour_manual(values=c("diff.exp.miRNA.target"="#5b6d5b","target.miRNA"="#484018", "diff.exp.miRNA" = "#e40017", "miRNA" = "#f4c983")) +
 labs(y= "log2(rrp45b[meanRPM])", x = "log2(WT[meanRPM])") +
  theme(aspect.ratio = 1)
ggsave(rrp45b_miRNAs, file = "/binf-isilon/PBgrp/xpj980/figures/20210423miRNAgenesscatterrrp45baraport.svg", units = "in", width = 15, height = 15, dpi = 300)

## now make the same scatterplot for ski2 and ski3 

joined_sign <- left_join(x = gathered_DESeqresults, y = genes) %>%
  filter(genotype == "ski2.5", padj < 0.05 , biotype == "miRNA") 

miRNAs_RPM_table <- normalized %>% 
  left_join(., genes, by = "gene") %>% 
  filter(., biotype %in% "miRNA") %>% 
  spread(., genotyp, mean_RPM) %>% 
  drop_na(., Name) 

sign_miRNAs <- left_join(miRNAs_RPM_table, joined_sign, by = "gene") %>% 
  mutate(sig = case_when(is.na(Name.y) ~ "no", 
                         TRUE ~ "yes")) %>%
  select(gene, biotype.x, Name.x, rrp45b, ski2.5, ski3, wt, sig)

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

ski2_miRNAs <- ggplot(to_plot3, aes(x=log2(wt), y=log2(ski2.5), color = both)) +
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
  geom_text_repel(aes(label=ifelse(both == "target.miRNA",as.character(Name.x),'')),hjust=0,vjust=0) +
  geom_text_repel(aes(label=ifelse(both == "diff.exp.miRNA",as.character(Name.x),'')),hjust=0,vjust=0) +
  #geom_text_repel(data=subset(to_plot3, rrp45b > 50 & wt > 50), aes(label=Name.x)) +
  ggtitle(label = "ski2.5", subtitle = "scatterplot of miRNA genes in ski2.5 vs WT") +
  #scale_color_manual(values = wes_palette(n=4, name="Darjeeling1"))+
  scale_colour_manual(values=c("diff.exp.miRNA.target"="#5b6d5b","target.miRNA"="#484018", "diff.exp.miRNA" = "#e40017", "miRNA" = "#f4c983")) +
  labs(y= "log2(ski2.5[meanRPM])", x = "log2(WT[meanRPM])") +
  theme(aspect.ratio = 1)
  
ggsave(ski2_miRNAs, file = "/binf-isilon/PBgrp/xpj980/figures/20210423miRNAgenesscatterski2araport.svg", units = "in", width = 15, height = 15, dpi = 300)



# ski3
joined_sign <- left_join(x = gathered_DESeqresults, y = genes) %>%
  filter(genotype == "ski3", padj < 0.05 , biotype == "miRNA") 

miRNAs_RPM_table <- normalized %>% 
  left_join(., genes, by = "gene") %>% 
  filter(., biotype %in% "miRNA") %>% 
  spread(., genotyp, mean_RPM) %>% 
  drop_na(., Name) 

sign_miRNAs <- left_join(miRNAs_RPM_table, joined_sign, by = "gene") %>% 
  mutate(sig = case_when(is.na(Name.y) ~ "no", 
                         TRUE ~ "yes")) %>%
  select(gene, biotype.x, Name.x, rrp45b, ski2.5, ski3, wt, sig)

left_joined <- left_join(x = gathered_DESeqresults, y = miRNAtargets, by = "gene", copy = FALSE)

ski3_test <- left_joined[left_joined$genotype=="ski3",] %>% 
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

ski3_miRNAs <- ggplot(to_plot3, aes(x=log2(wt), y=log2(ski3), color = both)) +
  geom_point() +
  scale_x_continuous(limits = c(0.01, NA)) +
  scale_y_continuous(limits = c(0.01, NA)) +
  #scale_y_log10(limits = c(0.01,NA), breaks = c(0, 1, 10, 100, 1000)) +
  #scale_x_log10(limits = c(0.01,NA), breaks = c(0, 1, 10, 100, 1000)) +
  coord_equal() +
  geom_abline(intercept = 0, slope = 1) +
  cowplot::theme_cowplot() +
  geom_text_repel(aes(label=ifelse(both == "target.miRNA",as.character(Name.x),'')),hjust=0,vjust=0) +
  geom_text_repel(aes(label=ifelse(both == "diff.exp.miRNA",as.character(Name.x),'')),hjust=0,vjust=0) +
  #geom_text_repel(data=subset(to_plot3, rrp45b > 50 & wt > 50), aes(label=Name.x)) +
  ggtitle(label = "ski3", subtitle = "scatterplot of miRNA genes in ski3 vs WT") +
  scale_colour_manual(values=c("diff.exp.miRNA.target"="#5b6d5b","target.miRNA"="#484018", "diff.exp.miRNA" = "#e40017", "miRNA" = "#f4c983")) +
  labs(y= "log2(ski3[meanRPM])", x = "log2(WT[meanRPM])") +
  theme(aspect.ratio = 1)
ggsave(ski3_miRNAs, file = "/binf-isilon/PBgrp/xpj980/figures/20210423miRNAgenesscatterski3araport.svg", units = "in", width = 15, height = 15, dpi = 300)

######Figure 2B #####

left_joined <- left_join(x = gathered_DESeqresults, y = miRNAtargets, by = "gene", copy = FALSE)

ski2_list <- left_joined[left_joined$genotype=="ski2.5",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(!is.na(miRNA)) %>%
  pull(gene) 

ski3_list <- left_joined[left_joined$genotype=="ski3",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(!is.na(miRNA)) %>%
  pull(gene) 

rrp45b_list <- left_joined[left_joined$genotype=="rrp45b",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(!is.na(miRNA)) %>%
  pull(gene) 

## eulers were made with these lists by Carlotta 
mariaslistfig2 <- list("ski2_list" = ski2_list, "ski3_list" = ski3_list, "rrp45b_list" = rrp45b_list, "ski2_nontarg_list" = ski2_other, "ski3_nontarg_list" = ski3_other,  "rrp45b_nontarg_list" = rrp45b_other)
save(mariaslistfig2, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/2015Vennlist.RData")

test <- venn.diagram(
  x = list(ski2_list, ski3_list, rrp45b_list), 
  category.names = c("ski2" , "ski3" , "rrp45b"),
  output = TRUE ,
  filename = NULL, 
  imagetype="svg" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#283148", "#913535", "#bbbbbb"),
  fill = c(alpha("#283148"), alpha("#913535"), alpha("#bbbbbb")),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("black", 'black', 'black'),
  rotation = 1
)
ggsave(test, filename = "/binf-isilon/PBgrp/xpj980/figures/202103022Bvennaraport.svg")


## and the "other genes" in the mutants, how many do they have in common? 

ski2_other <- left_joined[left_joined$genotype=="ski2.5",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(is.na(miRNA)) %>%
  pull(gene)   

ski3_other <- left_joined[left_joined$genotype=="ski3",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(is.na(miRNA)) %>%
  pull(gene) 

rrp45b_other <- left_joined[left_joined$genotype=="rrp45b",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(is.na(miRNA)) %>%
  pull(gene) 

mariaslistfig2 <- list("ski2_list" = ski2_list, "ski3_list" = ski3_list, "rrp45b_list" = rrp45b_list, "ski2_nontarg_list" = ski2_other, "ski3_nontarg_list" = ski3_other,  "rrp45b_nontarg_list" = rrp45b_other)


test <- venn.diagram(
  x = list(ski2_other, ski3_other, rrp45b_other), 
  category.names = c("ski2" , "ski3" , "rrp45b"),
  output = TRUE ,
  filename = NULL, 
  imagetype="svg" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#283148", "#913535", "#bbbbbb"),
  fill = c(alpha("#283148"), alpha("#913535"), alpha("#bbbbbb")),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("black", 'black', 'black'),
  rotation = 1
)
ggsave(test, filename = "/binf-isilon/PBgrp/xpj980/figures/202103022Bvenn_otheraraport.svg")
