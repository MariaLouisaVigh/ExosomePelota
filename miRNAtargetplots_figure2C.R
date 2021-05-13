##### miRNAtarget plots figure aC ####

## set working directory and load packages 
setwd("/binf-isilon/PBgrp/xpj980/datascratch/2015_exosome_mutants/04.individual_target_plots/new_target_plots/")

library(magrittr)
library(tidyverse)
library(reshape2)
library(stringr)
library(rtracklayer)
library(GenomicRanges)
library(devtools)
library(RColorBrewer)

###### load target list and make it into GRanges #######
targets <- read.table("/binf-isilon/PBgrp/xpj980/datascratch/2015_exosome_mutants/03.STAR_mapping/sam_files/annotated_ASRP_miRNA_target_sites.bed.txt", header = T, sep = "\t", quote = "", strip.white = T, stringsAsFactors = F) %>%
  as_tibble()

target_GR <- makeGRangesFromDataFrame(df = targets, keep.extra.columns = T, seqnames.field = "chr", start.field = "start", end.field = "stop", strand.field = "strand")
target_GR <- target_GR[!duplicated(target_GR)]

# list with 1bp GRange corresponding to cleavage site  
target_midpoint <- resize(target_GR, width = 1, fix = "start")
target_midpoint_df <- as.data.frame(target_midpoint)
target_midpoint_end_fix <- resize(target_GR, width = 1, fix = "end") 
target_midpoint_end_fix_df <- as.data.frame(target_midpoint_end_fix)

# gtf file subsetted for IRanges in targetGR 
target_id_list <- c(unique(target_midpoint$target_id))

mygtf <- import("/binf-isilon/PBgrp/xpj980/TAIR/Arabidopsis_thaliana.TAIR10.44.gtf")

gtf_sub <- subset(x = mygtf, 
                  mygtf$gene_id %in% target_id_list & 
                    mygtf$type %in% c("five_prime_utr", "CDS", "three_prime_utr"))

####### load bw files and rename them ######

## prepare names for bw forward and reverse files 
files_names_forward <- list.files("/binf-isilon/PBgrp/xpj980/datascratch/2015_exosome_mutants/03.STAR_mapping/sam_files/by_size/bw/", full.names = T, pattern = "Forward")
good_names_fw <- files_names_forward %>%
  str_remove(pattern = ".*by_size/bw//") %>%
  str_remove(pattern = ".Aligned.out") %>%
  str_remove(pattern = ".sorted.Forward.bw") %>%
  str_replace(pattern = "\\.", replacement = "_") %>%
  str_replace_all(pattern = "-", replacement = "\\.") %>% 
  str_replace(pattern = "heso1.3.", replacement = "heso1.3_") %>% # exceptional code for this library 
  str_replace(pattern = "cer7.3", replacement = "rrp45b") # exceptional code for this library 
names(files_names_forward) <- good_names_fw

files_names_reverse <- list.files("/binf-isilon/PBgrp/xpj980/datascratch/2015_exosome_mutants/03.STAR_mapping/sam_files/by_size/bw/", full.names = T, pattern = "Reverse")
good_names_rv <- files_names_reverse %>%
  str_remove(pattern = ".*by_size/bw//") %>%
  str_remove(pattern = ".Aligned.out") %>%
  str_remove(pattern = ".sorted.Reverse.bw") %>%
  str_replace(pattern = "\\.", replacement = "_") %>%
  str_replace_all(pattern = "-", replacement = "\\.") %>% 
  str_replace(pattern = "heso1.3.", replacement = "heso1.3_") %>% # exceptional code for this library
  str_replace(pattern = "cer7.3", replacement = "rrp45b") # exceptional code for this library
names(files_names_reverse) <- good_names_rv

fw_select <- files_names_forward %>%
  str_remove(pattern = ".*Aligned.out\\.") %>%
  str_remove(pattern = "\\.sorted.Forward.bw") %in% 21:24
rv_select <- files_names_reverse %>%
  str_remove(pattern = ".*Aligned.out\\.") %>%
  str_remove(pattern = "\\.sorted.Reverse.bw") %in% 21:24
all(fw_select==rv_select)

## import and subset bw files 

bw_list_fw <- lapply(files_names_forward[fw_select], function(file_name) {import(file_name, which=gtf_sub)})
bw_list_rv <- lapply(files_names_reverse[rv_select], function(file_name) {import(file_name, which=gtf_sub)})

df_list_fw <- lapply(seq_along(bw_list_fw), function(x) {
  test_name <- names(bw_list_fw[x])
  test <- bw_list_fw[[x]]
  test$sample <- test_name
  df <- test %>% as.data.frame() %>% tidyr::separate(sample, c("genotype", "replicate", "size"), sep="_") %>%
    mutate(strand="+")
  return(df)
})

df_list_rv <- lapply(seq_along(bw_list_rv), function(x) {
  test_name <- names(bw_list_rv[x])
  test <- bw_list_rv[[x]]
  test$sample <- test_name
  df <- test %>% as.data.frame() %>% tidyr::separate(sample, c("genotype", "replicate", "size"), sep="_") %>%
    mutate(strand="-")
  return(df)
})

gathered_bw <- purrr::reduce(c(df_list_fw,df_list_rv),.f = bind_rows)
gathered_bw_gr <- as(gathered_bw, "GRanges")

##### import featureCounts summary #######
#for total reads per genotype to merge with bw dataframes later 
featureCounts_sum <- read.table("/binf-isilon/PBgrp/xpj980/datascratch/2015_exosome_mutants/03.STAR_mapping/sam_files/featureCounts.summary", header = T, row.names = 1)  %>%
  as.data.frame()
new_names <- names(featureCounts_sum) %>% str_remove("\\.Aligned\\.out\\.sam") %>% 
  str_remove("set1\\.") %>% 
  str_replace(pattern = "heso1.3.", replacement = "heso1.3_") %>% # exceptional code for this library 
  str_replace(pattern = "cer7.3", replacement = "rrp45b") # exceptional code for this library 

featureCounts_renamed <- featureCounts_sum
names(featureCounts_renamed) <- new_names

total_alignments <- featureCounts_renamed[1,] %>%
  gather(key = "genotype", value = "totalread") %>% 
  separate(genotype, c("genotype", "replicate"), sep="_")  

#### one plot to test #### 
a_gene <- subset(gtf_sub, gtf_sub$gene_id == "AT2G28190")

a_gene_bw <- subsetByOverlaps(gathered_bw_gr, a_gene, ignore.strand = T)
a_bw_df <- data.frame(a_gene_bw) %>% 
  dplyr::filter(., genotype %in% c("rrp45b", "ski2.5", "Col.0", "ski3")) %>% # exceptional code for this library
  dplyr::filter(., genotype != "ski2.5" | replicate != "R2") # exceptional code for this library

a_gene_df <- data.frame(a_gene) %>%
  mutate(mean_RP10M=0,
         genotype="gene model")
a_target_midpoint <- subset(target_midpoint_df, target_midpoint_df$target_id == "AT2G28190")

a_target_midpoint <- a_target_midpoint[1,]  

RP10M <- a_bw_df %>% 
  left_join(., total_alignments, by = c("genotype" = "genotype", "replicate" = "replicate"), copy = F) %>%
  mutate(., RP10M=10000000*.$score/.$totalread) 

Check <- RP10M %>% 
  select(genotype, seqnames, start, end, strand, width, size) %>% 
  unique() %>% 
  group_by(genotype) %>%
  summarise(n = n())

new <- RP10M %>% complete(nesting(seqnames, start, end, width, strand, genotype, size), replicate, fill = list(RP10M = 0)) %>% 
  filter(., !(genotype == "ski2.5" & replicate == "R2"))

new %>% group_by(genotype, replicate) %>% summarise(n = n())

group_by(new, seqnames, start, end, width, strand, genotype, size) %>% 
  summarise(n = n()) %>% 
  pull(n) %>% 
  unique()

meanRP10M <- new %>% 
  group_by(genotype, size, strand, start, end) %>% 
  summarise("mean_RP10M" = mean(RP10M)) %>% 
  ungroup()

a_gene_df$genotype <- factor(a_gene_df$genotype, levels = c("Col.0", "ski2.5", "ski3", "rrp45b", "gene model")) # exceptional code for this library 
meanRP10M$genotype <- factor(meanRP10M$genotype, levels = c("Col.0", "ski2.5", "ski3", "rrp45b", "gene model")) # exceptional code for this library 

mytitle <- paste(a_target_midpoint$miRNA,
                 a_target_midpoint$target_id, 
                 a_target_midpoint$target_name,
                 sep=" - ") 

(plot <- ggplot(data=meanRP10M, aes(x=start, y=mean_RP10M, fill=size)) +
    scale_fill_brewer(palette = "Set1") +
    geom_bar(stat = "identity", position = "stack") +
    cowplot::theme_cowplot() +
    ggtitle(label = mytitle) +
    labs(y= "Reads per 10 million", x = "Gene coordinates", fill= "nt size") +
    #scale_y_continuous(position = "right", 
    #breaks = c(max(abs(meanRP10M$mean_RP10M))*-1, 0, max(abs(meanRP10M$mean_RP10M))),
    #limits = c(max(abs(meanRP10M$mean_RP10M))*-1, max(abs(meanRP10M$mean_RP10M))),
    #labels = scales::number_format(accuracy = 1,
    #decimal.mark = '.')) +
    scale_y_continuous(position = "right", 
                       breaks = c(-8, 0, 8),
                       limits = c(-8, 8),
                       labels = scales::number_format(accuracy = 1,
                                                      decimal.mark = '.')) +
    facet_grid(genotype~., switch = "both") + 
    geom_vline(xintercept = a_target_midpoint$start, lty=2) +
    geom_hline(yintercept = 0, size = 0.2) +
    geom_rect(data=a_gene_df, aes(ymin=min(meanRP10M$mean_RP10M)/2, ymax=(max(meanRP10M$mean_RP10M)/2), xmin=start,
                                  xmax=end, fill=type), alpha =0.5, inherit.aes = F) + 
    theme(strip.text.y = element_text(angle = 180, hjust = 1), 
          strip.background = element_blank(), 
          plot.title = element_text(size = 24, 
                                    face = "bold",
                                    color = "#22292F", 
                                    margin = margin(b = 8)), 
          axis.title.x = element_text(margin = margin(t = 15)),
          axis.title.y = element_text(margin = margin(r = 15))) +
    if (a_gene_df$strand[1] == "-") {scale_x_continuous(trans = "reverse", 
                                                        breaks = c(min(meanRP10M$start), a_target_midpoint$start, max(meanRP10M$start)))} else 
                                                        {scale_x_continuous(breaks = c(min(meanRP10M$start), a_target_midpoint$start, max(meanRP10M$start)))}
)

ggsave(plot, filename = paste0("/isdata/PBgrp/xpj980/figures/miRNAtargets/CSD2new2015.svg"), device = "svg", units = "cm", width = 18, height = 18)

######## targets described in paper ###### 
figure_targets <- c("AT2G28350", "AT2G34710", "AT1G30490", "AT2G45160", "AT3G15030", "AT1G48410", "AT5G43740", "AT2G28190")

##### LOOP #####

loop_vector <- unique(gtf_sub$transcript_id)

test_loop_vector <- loop_vector[1:10]

for (i in loop_vector) {
  a_gene <- subset(gtf_sub, gtf_sub$transcript_id == i)
  a_gene_bw <- subsetByOverlaps(gathered_bw_gr, a_gene, ignore.strand = T)
  a_bw_df <- data.frame(a_gene_bw) %>% 
    dplyr::filter(., genotype %in% c("Col.0.WT", "ski2.5", "ski3", "rrp45b")) 
  a_gene_df <- data.frame(a_gene) %>%
    mutate(mean_RP10M=0,
           genotype="gene model")
  a_target_midpoint <- subset(target_midpoint_df, target_midpoint_df$gene_id == unique(a_gene_df$gene_id))
  
  RP10M <- a_bw_df %>% 
    group_by(... = genotype, replicate) %>% 
    left_join(., total_alignments, by = c("genotype" = "genotype", "replicate" = "replicate"), copy = F) %>%
    ungroup() %>%
    mutate(., RP10M=10000000*.$score/.$totalread) 
  
  meanRP10M <- group_by(.data = RP10M, genotype, size, start, end) %>% 
    summarise("mean_RP10M" = mean(RP10M)) %>% 
    ungroup() 
  
  # filter away low reads 
  if (meanRP10M$mean_RP10M %>% abs() %>% max() < 2) {next} 
  
  # stuff for plotting 
  a_gene_df$genotype <- factor(a_gene_df$genotype, levels = c("Col.0.WT", "ski2.5", "ski3", "rrp45b", "gene model"))
  meanRP10M$genotype <- factor(meanRP10M$genotype, levels = c("Col.0.WT", "ski2.5", "ski3", "rrp45b", "gene model"))
  
  mytitle <- paste(a_target_midpoint$miRNA,
                   i, 
                   a_target_midpoint$target_name,
                   sep=" - ") 
  
  plot <- ggplot(data=meanRP10M, aes(x=start, y=mean_RP10M, fill=size)) +
    scale_fill_brewer(palette = "Set1") +
    geom_bar(stat = "identity", position = "stack") +
    cowplot::theme_cowplot() +
    ggtitle(label = mytitle) +
    labs(y= "Reads per 10 million", x = "Gene coordinates", fill= "nt size") +
    scale_y_continuous(position = "right", 
                       breaks = c(max(abs(meanRP10M$mean_RP10M))*-1, 0, max(abs(meanRP10M$mean_RP10M))),
                       limits = c(max(abs(meanRP10M$mean_RP10M))*-1, max(abs(meanRP10M$mean_RP10M))),
                       labels = scales::number_format(accuracy = 1,
                                                      decimal.mark = '.')) +
    facet_grid(genotype~., switch = "both") + 
    geom_vline(xintercept = a_target_midpoint$start, lty=2) +
    geom_rect(data=a_gene_df, aes(ymin=min(meanRP10M$mean_RP10M)/2, ymax=(max(meanRP10M$mean_RP10M)/2), xmin=start,
                                  xmax=end, fill=type), alpha =0.5, inherit.aes = F) + 
    theme(strip.text.y = element_text(angle = 180, hjust = 1), 
          strip.background = element_blank(), 
          plot.title = element_text(size = 24, 
                                    face = "bold",
                                    color = "#22292F", 
                                    margin = margin(b = 8)), 
          axis.title.x = element_text(margin = margin(t = 15)),
          axis.title.y = element_text(margin = margin(r = 15))) +
    if (a_gene_df$strand[1] == "-") {scale_x_continuous(trans = "reverse", 
                                                        breaks = c(min(meanRP10M$start), a_target_midpoint$start, max(meanRP10M$start)))} else 
                                                        {scale_x_continuous(breaks = c(min(meanRP10M$start), a_target_midpoint$start, max(meanRP10M$start)))} # next time use + or - to revert x-scale
  
  
  filename <- paste(a_target_midpoint$miRNA[1],
                    i,
                    a_target_midpoint$target_name[1],
                    ".svg", sep="_") %>%
    str_remove_all("/")
  
  ggsave(plot, filename = paste0("/isdata/PBgrp/xpj980/datascratch/figures/miRNAtargetplots", filename), device = "svg", units = "cm", width = 18, height = 18)
}  
