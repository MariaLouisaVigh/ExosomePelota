####### figure 6 ###### 
#### phasing clocks #####

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


# custum functions
rc_string <- function(string) {
  string %>% Biostrings::DNAStringSet() %>% 
    Biostrings::reverseComplement() %>% 
    as.character()
}

c("ATCG", "AGGT") %>% rc_string()

## load gtf and target list 
mygtf <- import("/binf-isilon/PBgrp/xpj980/TAIR/Arabidopsis_thaliana.TAIR10.44.gtf") 

targets <- read.table("/binf-isilon/PBgrp/xpj980/datascratch/2015_exosome_mutants/03.STAR_mapping/sam_files/annotated_ASRP_miRNA_target_sites.bed.txt", header = T, sep = "\t", quote = "", strip.white = T, stringsAsFactors = F) %>%
  as_tibble()

target_GR <- makeGRangesFromDataFrame(df = targets, keep.extra.columns = T, seqnames.field = "chr", start.field = "start", end.field = "stop", strand.field = "strand")
target_GR <- target_GR[!duplicated(target_GR)]

targets <- read_xlsx("/binf-isilon/PBgrp/xpj980/TAIR/miRNAtargetlist_Maria_degradome_done.xlsx", col_names = T)

### loop ###
target_vector <- c("AT2G39675", "AT1G62910", "AT2G34710", "AT4G18390", "AT1G48410")
genotypes_sub <- c("Col.0.WT", "ski2.5", "hen2.5", "rrp4")

bam_paths <- list.files("/binf-isilon/PBgrp/xpj980/datascratch/2019_exosome_mutants/04.STAR_mapping/sam_files/by_size/bam", full.names = T, pattern = "Aligned.out.21.sorted.bam")
index_paths <- list.files("/binf-isilon/PBgrp/xpj980/datascratch/2019_exosome_mutants/04.STAR_mapping/sam_files/by_size/bam/bai", full.names = T, pattern = "Aligned.out.21.sorted.bam")

bam_objects_list <- BamFileList(bam_paths, index = index_paths)

### LOOP
for (target in target_vector) {
  
  #GRange for target of interest
  to_subset <- subset(x = mygtf, 
                      mygtf$gene_id == target & mygtf$type %in% "gene")
  
  #GRange for targetsite of target of interest
  targetsite <- subsetByOverlaps(target_GR, to_subset)
  targetsite <- targetsite[1,]
  
  #find cleavage site of target of interest 
  target_midpoint <- targetsite %>% 
    promoters(., downstream = 11) %>% 
    resize(., width = 1, fix = "end")
  
  target_midpoint_df <- target_midpoint %>% as.data.frame()
  
  ### create GRanges of bam file on target
  what <- c("seq", "cigar")
  which <- GRanges(to_subset)
  param <- ScanBamParam(which=which, what=what)
  
  test <- lapply(bam_objects_list, function(x) {readGAlignments(x, param = param)})
  target_GRanges_list <- lapply(test, function(x) {as(x, "GRanges")})
  
  # rename GRange lists 
  file_names <- list.files("/binf-isilon/PBgrp/xpj980/datascratch/2019_exosome_mutants/04.STAR_mapping/sam_files/by_size/bam", full.names = T, pattern = "21.sorted.bam")
  names(target_GRanges_list)
  good_names <- file_names %>%
    str_remove(pattern = ".*by_size/bam/") %>%
    str_remove(pattern = ".Aligned.out.21.sorted.bam") %>%
    str_replace_all(pattern = "-", replacement = "\\.")
  
  names(target_GRanges_list) <- good_names
  
  ## make a column with genotype and a column with replicate 
  new_list <- lapply(seq_along(target_GRanges_list), function(x) {
    test_name <- names(target_GRanges_list[x])
    test <- target_GRanges_list[[x]]
    test$sample <- test_name
    df <- test %>% as.data.frame() %>% tidyr::separate(sample, c("genotype", "replicate"), sep="_") %>%
    mutate(target = rep(target))
    return(df)
  })
  
  assign(paste0("gathered_df_",target), purrr::reduce(new_list,.f = bind_rows))
  }
  
  

gathered_gr_AT1G62910 <- as(gathered_df_AT1G62910, "GRanges")
gathered_gr_AT2G34710 <- as(gathered_df_AT2G34710, "GRanges")
gathered_gr_AT2G39675 <- as(gathered_df_AT2G39675, "GRanges")  
gathered_gr_AT4G18390 <- as(gathered_df_AT4G18390, "GRanges")
gathered_gr_AT1G48410 <- as(gathered_df_AT1G48410, "GRanges")

for (target in target_vector) {

to_subset <- subset(x = mygtf, 
                      mygtf$gene_id == target & mygtf$type %in% "gene")

#GRange for targetsite of target of interest
targetsite <- subsetByOverlaps(target_GR, to_subset)
targetsite <- targetsite[1,]

#find cleavage site of target of interest 
assign(paste0("target_midpoint",target),targetsite %>% 
  promoters(., downstream = 11) %>% 
  resize(., width = 1, fix = "end"))

assign(paste0("target_midpoint_df", target), targetsite %>% 
         promoters(., downstream = 11) %>% 
         resize(., width = 1, fix = "end") %>% as_data_frame())
}

to_subset_m <- subset(x = mygtf, 
                              mygtf$gene_id == "AT1G62910" & mygtf$type %in% "gene") 

#GRange for targetsite of target of interest
targetsite <- subsetByOverlaps(target_GR, to_subset_m)
targetsite <- targetsite[1,]

test <- targetsite %>% 
  promoters(., downstream = 11) %>% 
  resize(., width = 1, fix = "end") 

##### make phasing columns ####
gathered_gr_AT1G62910$phase <- 1 + (start(gathered_gr_AT1G62910) - start(test)-1) %% 21
gathered_gr_AT1G62910$side <- case_when(start(gathered_gr_AT1G62910) < start(test) ~ "left", T ~ "right")
if (target_midpoint_dfAT1G62910$strand == "-") {gathered_gr_AT1G62910$CF <- case_when(gathered_gr_AT1G62910$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT1G62910$CF <- case_when(gathered_gr_AT1G62910$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT1G62910 <- as.data.frame(gathered_gr_AT1G62910)

gathered_gr_AT2G34710$phase <- 1 + (start(gathered_gr_AT2G34710) - start(target_midpointAT2G34710)-1) %% 21
gathered_gr_AT2G34710$side <- case_when(start(gathered_gr_AT2G34710) < start(target_midpointAT2G34710) ~ "left", T ~ "right")
if (target_midpoint_dfAT2G34710$strand == "-") {gathered_gr_AT2G34710$CF <- case_when(gathered_gr_AT2G34710$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT2G34710$CF <- case_when(gathered_gr_AT2G34710$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT2G34710 <- as.data.frame(gathered_gr_AT2G34710)

gathered_gr_AT2G39675$phase <- 1 + (start(gathered_gr_AT2G39675) - start(target_midpointAT2G39675)-1) %% 21
gathered_gr_AT2G39675$side <- case_when(start(gathered_gr_AT2G39675) < start(target_midpointAT2G39675) ~ "left", T ~ "right")
if (target_midpoint_dfAT2G39675$strand == "-") {gathered_gr_AT2G39675$CF <- case_when(gathered_gr_AT2G39675$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT2G39675$CF <- case_when(gathered_gr_AT2G39675$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT2G39675 <- data.frame(gathered_gr_AT2G39675)

gathered_gr_AT4G18390$phase <- 1 + (start(gathered_gr_AT4G18390) - start(target_midpointAT4G18390)-1) %% 21
gathered_gr_AT4G18390$side <- case_when(start(gathered_gr_AT4G18390) < start(target_midpointAT4G18390) ~ "left", T ~ "right")
if (target_midpoint_dfAT4G18390$strand == "-") {gathered_gr_AT4G18390$CF <- case_when(gathered_gr_AT4G18390$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT4G18390$CF <- case_when(gathered_gr_AT4G18390$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT4G18390 <- data.frame(gathered_gr_AT4G18390)
  

#### make a dataframe to plot & count phases ####
dt_plot_AT1G62910 <- gathered_AT1G62910 %>%
    filter(cigar=="21M") %>%
    group_by(phase, strand, side, CF, genotype, replicate, target) %>%
    summarise(count=n()) %>%
    ungroup()

to_plot_AT1G62910 <- if (dt_plot_AT1G62910[which.max(dt_plot_AT1G62910$count),]$CF == "3CF") {filter(dt_plot_AT1G62910, CF %in% "3CF")} else {filter(dt_plot_AT1G62910, CF %in% "5CF")}

summed_to_plot_AT1G62910 <- to_plot_AT1G62910 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>%
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT1G62910"))

dt_plot_AT2G34710 <- gathered_AT2G34710 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  ungroup()

to_plot_AT2G34710 <- if (dt_plot_AT2G34710[which.max(dt_plot_AT2G34710$count),]$CF == "3CF") {filter(dt_plot_AT2G34710, CF %in% "3CF")} else {filter(dt_plot_AT2G34710, CF %in% "5CF")}

summed_to_plot_AT2G34710 <- to_plot_AT2G34710 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT2G34710"))
   
dt_plot_AT2G39675 <- gathered_AT2G39675 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  ungroup()

to_plot_AT2G39675 <- if (dt_plot_AT2G39675[which.max(dt_plot_AT2G39675$count),]$CF == "3CF") {filter(dt_plot_AT2G39675, CF %in% "3CF")} else {filter(dt_plot_AT2G39675, CF %in% "5CF")}

summed_to_plot_AT2G39675 <- to_plot_AT2G39675 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT2G39675"))

dt_plot_AT4G18390 <- gathered_AT4G18390 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  ungroup()

to_plot_AT4G18390 <- if (dt_plot_AT4G18390[which.max(dt_plot_AT4G18390$count),]$CF == "3CF") {filter(dt_plot_AT4G18390, CF %in% "3CF")} else {filter(dt_plot_AT4G18390, CF %in% "5CF")}

summed_to_plot_AT4G18390 <- to_plot_AT4G18390 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT4G18390"))
  
all_gathered_df <- bind_rows(to_plot_AT1G62910, to_plot_AT2G34710, to_plot_AT2G39675, to_plot_AT4G18390)
all_gathered_sum_reads <- bind_rows(summed_to_plot_AT1G62910, summed_to_plot_AT2G34710, summed_to_plot_AT2G39675, summed_to_plot_AT4G18390)
all_gathered_df_sub <-subset(all_gathered_df, genotype %in% genotypes_sub) 

all_gathered_df_sub$genotype <- factor(all_gathered_df_sub$genotype, levels = c("Col.0.WT", "ski2.5", "hen2.5", "rrp4"))
all_gathered_sum_reads$genotype <- factor(all_gathered_sum_reads$genotype, levels = c("Col.0.WT", "ski2.5", "hen2.5", "rrp4"))
all_gathered_df_sub$target <- factor(all_gathered_df_sub$target, levels = c("AT1G62910","AT2G39675", "AT2G34710", "AT4G18390"))
all_gathered_sum_reads$target <- factor(all_gathered_sum_reads$target, levels = c("AT1G62910","AT2G39675", "AT2G34710", "AT4G18390")) 


read_intensity_plot <- all_gathered_sum_reads %>% 
  filter(., genotype != "Col.0.WT") %>%
  ggplot(., aes(x=genotype, y=sum, fill = target, alpha = 0.5)) +
  stat_summary(geom = "bar", fun = "mean", position = "dodge", aes(width=0.5)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.5, aes(width=0.5))+
  #geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~target, scales = "free", nrow = 4) +
  scale_fill_brewer(palette = "Dark2") + 
  cowplot::theme_cowplot()
ggsave(read_intensity_plot, filename = paste0("/isdata/PBgrp/xpj980/figures/readsbarplotfigure6b.svg"), units = "cm", width = 18, height = 18)


to_plot_AT1G62910$genotype <- factor(to_plot_AT1G62910$genotype, levels = c("Col.0.WT", "ski2.5", "hen2.5", "rrp4"))
to_plot_AT1G62910$replicate <- factor(to_plot_AT1G62910$replicate, levels = c("R3", "R2", "R1"))
together <- to_plot_AT1G62910 %>% 
  filter(genotype %in% genotypes_sub) %>% 
  group_by(genotype, replicate, target) %>%
  mutate(ncount = count/sum(count)) %>%
  ggplot(aes(x=factor(phase, levels = 1:21), y=ncount)) +
  geom_bar(stat="identity", position= "identity", aes(group=replicate, fill=replicate), col="black", alpha=.8) + 
  coord_polar() +
  facet_grid(target~genotype, switch="both") + 
  #ggtitle(label = my_title) +
  theme(panel.grid = element_line(
    size = 0.25, linetype = 'solid', colour = "grey"),
    strip.background = element_rect(color="black"), 
    axis.title.x = element_text(margin = margin(t = 15))) +
  scale_fill_brewer(palette = rev("Reds")) +
  scale_y_continuous(limits = c(0,1))
  

# plot with wt, ski2, hen2, rrp4 
  plot <- all_gathered_df_sub %>% 
    group_by(genotype, replicate, target) %>%
    mutate(ncount = count/sum(count)) %>%
    ggplot(aes(x=factor(phase, levels = 1:21), y=ncount)) +
    geom_bar(stat="identity", position= "identity", aes(group=target, fill=target), col="black", alpha=.8) + 
    coord_polar() +
    facet_grid(target~genotype, switch="both") + 
    #ggtitle(label = my_title) +
    theme(panel.grid = element_line(
          size = 0.25, linetype = 'solid', colour = "grey"),
          strip.background = element_rect(color="black"), 
          axis.title.x = element_text(margin = margin(t = 15))) +
    scale_fill_brewer(palette = "Dark2") +
    scale_y_continuous(limits = c(0,1))
  
  ggsave(plot, filename = paste0("/isdata/PBgrp/xpj980/figures/figure6A.svg"), units = "cm", width = 18, height = 18)
  


#### investigating modifications on the end ####
  dt_plot_AT2G34710_invest <- gathered_AT2G34710 %>%
    filter(cigar!="21M" & strand %in% "-" & phase %in% 2)  
  
test_grep_PHB <- gathered_AT2G34710[grep("CATCCCAATCATCTGA", gathered_AT2G34710$seq),] %>%
  filter(genotype %in% c("Col.0.WT", "ski2.5", "hen2.5", "rrp4")) 


test_grep_PHB$new_c <- with(test_grep_PHB, ifelse(seq == "CTTCATCCCAATCATCTGAAC" & cigar == "21M",0, 
                                                  ifelse(cigar == "2S19M",2, 
                                                         ifelse (cigar =="1S20M", 1, 
                                                                 ifelse(cigar =="3S18M", 3, 
                                                                        ifelse(seq == "TCATCCCAATCATCTGAACCC" & cigar == "21M",-2,
                                                                               ifelse(seq == "TCATCCCAATCATCTGAACCT" & cigar == "20M1S", -2, 
                                                                                      ifelse(seq =="CATCCCAATCATCTGAACCCA" & cigar == "21M", -3, 
                                                                                             ifelse(seq == "TTCATCCCAATCATCTGAACC" & cigar == "21M", -3, "no"))))))))) %>%
  as.numeric()

test_grep_PHB %>% 
  group_by(new_c) %>%
summarise(sum(new_c))

to_plot <- test_grep_PHB %>% 
  group_by(genotype, new_c) %>%
  add_tally() %>%
  ggplot(., aes(x = genotype, y= new_c, size = n, color = target)) + 
  geom_point() +  
  cowplot::theme_cowplot() +
  scale_color_manual(values = "#5d54a4") + 
  scale_y_continuous(limits = c(-3,3), breaks = c(-3,-2,-1,0,1,2,3))

test_grep_TAS1c <- gathered_AT2G39675[grep("GAGAAAAATCA", gathered_AT2G39675$seq),] %>%
  filter(genotype %in% c("Col.0.WT", "ski2.5", "hen2.5", "rrp4"))
 

test_grep_TAS1c$new_c <- with(test_grep_TAS1c, ifelse(seq == "TAAATGGTCTATTCGCTTGTA" & cigar == "21M", 0,
                                                      ifelse(cigar == "16M5S", 1, 
                                                      ifelse(cigar == "18M3S", 1, "no")))) %>%
  as.numeric()

unique(test_grep_TAS1c$cigar)
subset(test_grep_TAS1c, cigar == "16M5S")
together <- bind_rows(test_grep_PHB, test_grep_TAS1c, test_grep_TCP2) 

to_plot <- together %>%
  filter(., new_c !="NA") %>%
  group_by(genotype, new_c, target) %>%
  add_tally() %>%
  ggplot(., aes(x = genotype, y= new_c, size = n, color = target)) + 
  geom_point() +  
  facet_wrap(~target) +
  cowplot::theme_cowplot() +
  scale_color_manual(values = c("#158467","#5d54a4", "#318fb5")) + 
  scale_y_continuous(limits = c(-3,7), breaks = c(-3,-2,-1,0,1,2,3,4,5,6,7))
ggsave(to_plot, filename = paste0("/isdata/PBgrp/xpj980/figures/trimtail.svg"), units = "cm", width = 18, height = 18)


test_grep_TCP2 <- gathered_AT4G18390[grep("GGGTCCCCCT", gathered_AT4G18390$seq),] %>%
  filter(genotype %in% c("Col.0.WT",  "ski2.5", "hen2.5", "rrp4")) 

test_grep_TCP2$new_c <- with(test_grep_TCP2, ifelse(seq == "AAAGGGTCCCCCTGTTGAAAC" & cigar == "1S20M", 3, 
                                                    ifelse(cigar=="6S15M", 7, 
                                                           ifelse(cigar == "4S17M", 6, "no")))) %>%
  as.numeric()

test_grep_Rgene <- gathered_AT1G62910[grep("TCCAAATGTAGTCA", gathered_AT1G62910$seq),]
  
#### individual targets for figure #####

#### Tas1C ####

### target site
to_subset_m <- subset(x = mygtf, 
                      mygtf$gene_id == "AT2G39675" & mygtf$type %in% "gene") 

#GRange for targetsite of target of interest
targetsite <- subsetByOverlaps(target_GR, to_subset_m)
targetsite <- targetsite[1,]

test <- targetsite %>% 
  promoters(., downstream = 13) %>% 
  resize(., width = 1, fix = "end") 

gathered_gr_AT2G39675$phase <- 1 + (start(gathered_gr_AT2G39675) - start(test)-1) %% 21
gathered_gr_AT2G39675$side <- case_when(start(gathered_gr_AT2G39675) < start(test) ~ "left", T ~ "right")
if (target_midpoint_dfAT2G39675$strand == "-") {gathered_gr_AT2G39675$CF <- case_when(gathered_gr_AT2G39675$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT2G39675$CF <- case_when(gathered_gr_AT2G39675$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT2G39675 <- data.frame(gathered_gr_AT2G39675)

dt_plot_AT2G39675 <- gathered_AT2G39675 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  filter(genotype %in% genotypes_sub) %>% 
  ungroup()

to_plot_AT2G39675 <- if (dt_plot_AT2G39675[which.max(dt_plot_AT2G39675$count),]$CF == "3CF") {filter(dt_plot_AT2G39675, CF %in% "3CF")} else {filter(dt_plot_AT2G39675, CF %in% "5CF")}

summed_to_plot_AT2G39675 <- to_plot_AT2G39675 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT2G39675"))

to_plot_AT2G39675$genotype <- factor(to_plot_AT2G39675$genotype, levels = c("Col.0.WT", "ski2.5", "hen2.5", "rrp4"))
to_plot_AT2G39675$replicate <- factor(to_plot_AT2G39675$replicate, levels = c("R3", "R2", "R1"))

together <- to_plot_AT2G39675 %>% 
  group_by(phase, genotype, target, strand) %>%
  summarise(meancount = mean(count)) %>%
  group_by(genotype, target, strand) %>%
  mutate(ncount = meancount/sum(meancount)) 

tas1 <- together %>% 
  ggplot(aes(x=factor(phase, levels = 1:21), y=ncount)) +
  geom_bar(stat="identity", position= "identity", aes(fill=strand), col="black", size=0.25, alpha=.9) +
  geom_hline(yintercept = seq(0, 4, by = 1), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.25), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.50), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.75), colour = "darkgrey", size = 0.1) +
  geom_vline(xintercept = 1:21, colour = "darkgrey", size = 0.1) + 
  coord_polar() +
  facet_grid(target~genotype, switch="both") + 
  theme(panel.grid = element_blank()) +
  #ggtitle(label = my_title) +
  #theme(panel.grid = element_line(
  #size = 0.25, linetype = 'solid', colour = "grey"),
  #strip.background = element_rect(color="black"), 
  #axis.title.x = element_text(margin = margin(t = 15))) +
  scale_fill_manual(values = c("#7ea04d", "#335d2d")) +
  scale_y_continuous(limits = c(0,1))
ggsave(tas1, filename = paste0("/isdata/PBgrp/xpj980/figures/figure6ATas1geneavg.svg"), units = "cm", width = 18, height = 18)

summed_to_plot_AT2G39675$genotype <- factor(summed_to_plot_AT2G39675$genotype, levels = c("Col.0.WT", "ski2.5", "hen2.5", "rrp4"))

read_intensity_plot <- summed_to_plot_AT2G39675 %>% 
  filter(target %in% "AT2G39675") %>%
  ggplot(., aes(x=genotype, y=sum, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge", aes(width=0.5)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.5, aes(width=0.5))+
  #geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~target, scales = "free", nrow = 4) +
  scale_fill_brewer(palette = "Greens") + 
  cowplot::theme_cowplot()
ggsave(read_intensity_plot, filename = paste0("/isdata/PBgrp/xpj980/figures/figure6ATas1genereadsavg.svg"), units = "cm", width = 18, height = 18)


#### PHB #### 
to_subset_m <- subset(x = mygtf, 
                      mygtf$gene_id == "AT2G34710" & mygtf$type %in% "gene") 

#GRange for targetsite of target of interest
targetsite <- subsetByOverlaps(target_GR, to_subset_m)
targetsite <- targetsite[1,]

test <- targetsite %>% 
  promoters(., downstream = 11) %>% 
  resize(., width = 1, fix = "end") 

gathered_gr_AT2G34710$phase <- 1 + (start(gathered_gr_AT2G34710) - start(test)-1) %% 21
gathered_gr_AT2G34710$side <- case_when(start(gathered_gr_AT2G34710) < start(test) ~ "left", T ~ "right")
if (target_midpoint_dfAT2G34710$strand == "-") {gathered_gr_AT2G34710$CF <- case_when(gathered_gr_AT2G34710$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT2G34710$CF <- case_when(gathered_gr_AT2G34710$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT2G34710 <- as.data.frame(gathered_gr_AT2G34710)

dt_plot_AT2G34710 <- gathered_AT2G34710 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  filter(genotype %in% genotypes_sub) %>%
  ungroup()

to_plot_AT2G34710 <- if (dt_plot_AT2G34710[which.max(dt_plot_AT2G34710$count),]$CF == "3CF") {filter(dt_plot_AT2G34710, CF %in% "3CF")} else {filter(dt_plot_AT2G34710, CF %in% "5CF")}

summed_to_plot_AT2G34710 <- to_plot_AT2G34710 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT2G34710"))

to_plot_AT2G34710$genotype <- factor(to_plot_AT2G34710$genotype, levels = c("Col.0.WT", "ski2.5", "hen2.5", "rrp4"))
to_plot_AT2G34710$replicate <- factor(to_plot_AT2G34710$replicate, levels = c("R3", "R2", "R1"))
together <- to_plot_AT2G34710 %>% 
  group_by(phase, genotype, target, strand) %>%
  summarise(meancount = mean(count)) %>%
  group_by(genotype, target, strand) %>%
  mutate(ncount = meancount/sum(meancount)) 

test <- together %>% 
  group_by(genotype, strand) %>% 
  #mutate(avg = sum(ncount)/3)
  summarise(sum = sum(ncount))

PHB <- together %>% 
  group_by(genotype, target, strand) %>%
  ggplot(aes(x=factor(phase, levels = 1:21), y=ncount)) +
  #stat_summary(geom = "bar", fun = "mean", position = "dodge", aes(fill=strand)) +
  geom_bar(stat="identity", position= "identity", aes(fill=strand), col="black", size=0.25, alpha=.9) +
  geom_hline(yintercept = seq(0, 4, by = 1), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.25), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.50), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.75), colour = "darkgrey", size = 0.1) +
  geom_vline(xintercept = 1:21, colour = "darkgrey", size = 0.1) + 
  coord_polar() +
  facet_grid(target~genotype, switch="both") + 
  theme(panel.grid = element_blank()) +
  #ggtitle(label = my_title) +
  #theme(panel.grid = element_line(
  #size = 0.25, linetype = 'solid', colour = "grey"),
  #strip.background = element_rect(color="black"), 
  #axis.title.x = element_text(margin = margin(t = 15))) +
  scale_fill_manual(values = c("#31112c", "#7d0633")) +
  scale_y_continuous(limits = c(0,1))
ggsave(PHB, filename = paste0("/isdata/PBgrp/xpj980/figures/figure6APHBavg.svg"), units = "cm", width = 18, height = 18)

summed_to_plot_AT2G34710$genotype <- factor(summed_to_plot_AT2G34710$genotype, levels = c("Col.0.WT", "ski2.5", "hen2.5", "rrp4"))

read_intensity_plot <- summed_to_plot_AT2G34710 %>% 
  filter(target %in% "AT2G34710") %>%
  ggplot(., aes(x=genotype, y=sum, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge", aes(width=0.5)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.5, aes(width=0.5))+
  #geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~target, scales = "free", nrow = 4) +
  scale_fill_brewer(palette = "Purples") + 
  cowplot::theme_cowplot()
ggsave(read_intensity_plot, filename = paste0("/isdata/PBgrp/xpj980/figures/figure6APHBreadsavg.svg"), units = "cm", width = 18, height = 18)

#### TCP2 ####

to_subset_m <- subset(x = mygtf, 
                      mygtf$gene_id == "AT4G18390" & mygtf$type %in% "gene") 

#GRange for targetsite of target of interest
targetsite <- subsetByOverlaps(target_GR, to_subset_m)
targetsite <- targetsite[1,]

test <- targetsite %>% 
  promoters(., downstream = 11) %>% 
  resize(., width = 1, fix = "end") 

gathered_gr_AT4G18390$phase <- 1 + (start(gathered_gr_AT4G18390) - start(test)-1) %% 21
gathered_gr_AT4G18390$side <- case_when(start(gathered_gr_AT4G18390) < start(test) ~ "left", T ~ "right")
if (target_midpoint_dfAT4G18390$strand == "-") {gathered_gr_AT4G18390$CF <- case_when(gathered_gr_AT4G18390$side =="right" ~ "5CF", T ~ "3CF")} else {gathered_gr_AT4G18390$CF <- case_when(gathered_gr_AT4G18390$side =="right" ~ "3CF", T ~ "5CF")}
gathered_AT4G18390 <- data.frame(gathered_gr_AT4G18390)

dt_plot_AT4G18390 <- gathered_AT4G18390 %>%
  filter(cigar=="21M") %>%
  group_by(phase, strand, side, CF, genotype, replicate, target) %>%
  summarise(count=n()) %>%
  filter(genotype %in% genotypes_sub) %>% 
  ungroup()

to_plot_AT4G18390 <- if (dt_plot_AT4G18390[which.max(dt_plot_AT4G18390$count),]$CF == "3CF") {filter(dt_plot_AT4G18390, CF %in% "3CF")} else {filter(dt_plot_AT4G18390, CF %in% "5CF")}

summed_to_plot_AT4G18390 <- to_plot_AT4G18390 %>% 
  filter(genotype %in% genotypes_sub) %>%
  group_by(genotype, replicate) %>% 
  summarise(sum = sum(count)) %>% 
  #summarise(mean_sum = mean(sum)) %>%
  mutate(target = rep("AT4G18390"))

to_plot_AT4G18390$genotype <- factor(to_plot_AT4G18390$genotype, levels = c("Col.0.WT", "ski2.5", "hen2.5", "rrp4"))
to_plot_AT4G18390$replicate <- factor(to_plot_AT4G18390$replicate, levels = c("R3", "R2", "R1"))
together <- to_plot_AT4G18390 %>% 
  group_by(phase, genotype, strand, target) %>%
  summarise(meancount = mean(count)) %>%
  group_by(genotype, strand, target) %>%
  mutate(ncount = meancount/sum(meancount)) 

test <- together %>% 
  group_by(genotype, strand) %>% 
  #mutate(avg = sum(ncount)/3)
  summarise(sum = sum(ncount))  
  
  tcp2 <- together %>% 
  ggplot(aes(x=factor(phase, levels = 1:21), y=ncount)) +
  geom_bar(stat="identity", position= "identity", aes(fill=strand), col="black", size=0.25) +
  geom_hline(yintercept = seq(0, 4, by = 1), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.25), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.50), colour = "darkgrey", size = 0.1) +
  geom_hline(yintercept = seq(0, 4, by = 0.75), colour = "darkgrey", size = 0.1) +
  geom_vline(xintercept = 1:21, colour = "darkgrey", size = 0.1) + 
  coord_polar() +
  facet_grid(target~genotype, switch="both") + 
  theme(panel.grid = element_blank()) +
  #ggtitle(label = my_title) +
  #theme(panel.grid = element_line(
  #size = 0.25, linetype = 'solid', colour = "grey"),
  #strip.background = element_rect(color="black"), 
  #axis.title.x = element_text(margin = margin(t = 15))) +
  scale_fill_manual(values = c("#00334e", "#145374")) +
  scale_y_continuous(limits = c(0,1))
ggsave(tcp2, filename = paste0("/isdata/PBgrp/xpj980/figures/figure6ATCP2avg.svg"), units = "cm", width = 18, height = 18)

summed_to_plot_AT4G18390$genotype <- factor(summed_to_plot_AT4G18390$genotype, levels = c("Col.0.WT", "ski2.5", "hen2.5", "rrp4"))

read_intensity_plot <- summed_to_plot_AT4G18390 %>% 
  filter(target %in% "AT4G18390") %>%
  ggplot(., aes(x=genotype, y=sum, fill = genotype)) + 
  stat_summary(geom = "bar", fun = "mean", position = "dodge", aes(width=0.5)) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", position = "dodge", size = 0.5, aes(width=0.5))+
  #geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~target, scales = "free", nrow = 4) +
  scale_fill_brewer(palette = "Blues") + 
  cowplot::theme_cowplot()
ggsave(read_intensity_plot, filename = paste0("/isdata/PBgrp/xpj980/figures/figure6ATCP2readsavg.svg"), units = "cm", width = 18, height = 18)

