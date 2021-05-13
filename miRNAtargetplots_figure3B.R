#### target plots for figure 3 - pelota figure 

setwd("/binf-isilon/PBgrp/xpj980/datascratch/2019_exosome_mutants/05.target_plots/")

library(magrittr)
library(tidyverse)
library(reshape2)
library(stringr)
library(rtracklayer)
library(GenomicRanges)
library(devtools)
library(RColorBrewer)
library(forcats)
library(readxl)
###### load target list and make it into GRanges #######
miRNAtargets_rawtable <- readxl::read_xlsx("/binf-isilon/PBgrp/xpj980/TAIR/miRNAtargetlist_Maria_degradome_done.xlsx", col_names = T)

miRNAtargets <- filter(miRNAtargets_rawtable, !duplicated(miRNAtargets_rawtable$target_id))

target_GR <- makeGRangesFromDataFrame(df = miRNAtargets, keep.extra.columns = T, seqnames.field = "chr", start.field = "start", end.field = "stop", strand.field = "strand")
target_GR <- target_GR[!duplicated(target_GR$target_id)]
target_GR_dataframe <- as.data.frame(target_GR) %>% 
  rename(target_id = "gene_id")

# list with 1bp GRange corresponding to cleavage site  
target_midpoint <- target_GR %>% 
  promoters(., downstream = 11) %>% 
  resize(., width = 1, fix = "end")
target_midpoint_df <- as.data.frame(target_midpoint) %>% 
  rename(target_id = "gene_id")

# gtf file subsetted for transcript_ids corresponding to gene_ids in targetGR 

mygtf <- import("/binf-isilon/PBgrp/xpj980/TAIR/Arabidopsis_thaliana.TAIR10.44.gtf") 
mygtf_df <- import("/binf-isilon/PBgrp/xpj980/TAIR/Arabidopsis_thaliana.TAIR10.44.gtf") %>% 
  as.data.frame()

#first put transcript_id in target_GR 
target_transcript_id <- left_join(mygtf_df, target_GR_dataframe, by = "gene_id", copy = FALSE) %>% 
  select(gene_id, transcript_id, miRNA) %>% 
  drop_na()

# list to subset with 
transcript_id_list <- c(unique(target_transcript_id$transcript_id))

gtf_sub <- subset(x = mygtf, 
                  mygtf$transcript_id %in% transcript_id_list & 
                    mygtf$type %in% c("five_prime_utr", "CDS", "three_prime_utr"))

######## import bw files ########

##### load bw files and rename them 
files_names_forward <- list.files("/binf-isilon/PBgrp/xpj980/datascratch/2019_exosome_mutants/04.STAR_mapping/sam_files/by_size/bw", full.names = T, pattern = "Forward")
good_names_fw <- files_names_forward %>%
  str_remove(pattern = ".*by_size/bw/") %>%
  str_remove(pattern = ".Aligned.out") %>%
  str_remove(pattern = ".sorted.Forward.bw") %>%
  str_replace(pattern = "\\.", replacement = "_") %>%
  str_replace_all(pattern = "-", replacement = "\\.")
names(files_names_forward) <- good_names_fw


files_names_reverse <- list.files("/binf-isilon/PBgrp/xpj980/datascratch/2019_exosome_mutants/04.STAR_mapping/sam_files/by_size/bw", full.names = T, pattern = "Reverse")
good_names_rv <- files_names_reverse %>%
  str_remove(pattern = ".*by_size/bw/") %>%
  str_remove(pattern = ".Aligned.out") %>%
  str_remove(pattern = ".sorted.Reverse.bw") %>%
  str_replace(pattern = "\\.", replacement = "_") %>%
  str_replace_all(pattern = "-", replacement = "\\.")
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
## import featureCounts summary for total reads per genotype to merge with tidyMetaProfile later 
featureCounts_sum <- read.table("/binf-isilon/PBgrp/xpj980/datascratch/2019_exosome_mutants/04.STAR_mapping/sam_files/featureCounts.summary", header = T, row.names = 1)  %>%
  as.data.frame()
new_names <- names(featureCounts_sum) %>% str_remove("\\.Aligned\\.out\\.sam")

featureCounts_renamed <- featureCounts_sum
names(featureCounts_renamed) <- new_names

total_alignments <- featureCounts_renamed[1,] %>%
  gather(key = "genotype", value = "totalread") %>% 
  separate(genotype, c("genotype", "replicate"), sep="_")                                       

load("/binf-isilon/PBgrp/xpj980/R_scripts_collected/2019exoDESEQpelota.RData")

left_joined <- left_joined <- left_join(x = gathered_DESeqresults_pelota, y = miRNAtargets, by = c("gene"="target_id"))

GKpelota_list <- left_joined[left_joined$genotype=="GKpelota",] %>% 
  filter(padj<0.05) %>% 
  filter(log2FoldChange > 0) %>% 
  filter(!is.na(miRNA)) %>%
  pull(gene) 

#### find specific gene to plot to test if script works 
a_gene <- subset(gtf_sub, gtf_sub$transcript_id == "AT5G43730.5")
a_gene_bw <- subsetByOverlaps(gathered_bw_gr, a_gene, ignore.strand = T)
a_bw_df <- data.frame(a_gene_bw) %>% 
  dplyr::filter(., genotype %in% c("Col.0.WT", "ski2.5", "GKpelota", "SAILpelota"))
a_gene_df <- data.frame(a_gene) %>%
  mutate(mean_RP10M=0,
         genotype="gene model")
a_target_midpoint <- target_midpoint_df %>% subset(., target_midpoint_df$gene_id == unique(a_gene_df$gene_id))

#a_target_midpoint <- a_target_midpoint[1,]  

RP10M <- a_bw_df %>% 
  dplyr::left_join(., total_alignments, by = c("genotype" = "genotype", "replicate" = "replicate"), copy = F) %>%
  mutate(., RP10M=10000000*.$score/.$totalread) %>% glimpse() 

Check <- RP10M %>% 
  select(genotype, seqnames, start, end, strand, width, size) %>% 
  unique() %>% 
  group_by(genotype) %>%
  summarise(n = n())

new <- RP10M %>% complete(nesting(seqnames, start, end, width, strand, genotype, size), replicate, fill = list(RP10M = 0)) %>% 
  filter(., !(genotype == "SAILpelota" & replicate %in% c("R2", "R3")))

new %>% group_by(genotype, replicate) %>% summarise(n = n())

group_by(new, seqnames, start, end, width, strand, genotype, size) %>% 
  summarise(n = n()) %>% 
  pull(n) %>% 
  unique()

meanRP10M <- new %>% 
  group_by(genotype, size, strand, start, end) %>% 
  summarise("mean_RP10M" = mean(RP10M)) %>% 
  ungroup() 

a_gene_df$genotype <- factor(a_gene_df$genotype, levels = c("Col.0.WT", "ski2.5","SAILpelota", "GKpelota", "gene model"))
meanRP10M$genotype <- factor(meanRP10M$genotype, levels = c("Col.0.WT", "ski2.5", "SAILpelota", "GKpelota", "gene model"))

mytitle <- paste(a_target_midpoint$miRNA,
                 a_gene_df$transcript_id, 
                 a_target_midpoint$target_name,
                 sep=" - ") 

(plot <- ggplot(data=meanRP10M, aes(x=start, y=mean_RP10M, fill=size)) +
    scale_fill_brewer(palette = "Set1") +
    geom_bar(stat = "identity", position = "stack") +
    cowplot::theme_cowplot() +
    ggtitle(label = mytitle) +
    labs(y= "Reads per 10 million", x = "Gene coordinates", fill= "nt size") +
     #scale_y_continuous(position = "right", 
      #                  breaks = c(max(abs(meanRP10M$mean_RP10M))*-1, 0, max(abs(meanRP10M$mean_RP10M))),
       #                 limits = c(max(abs(meanRP10M$mean_RP10M))*-1, max(abs(meanRP10M$mean_RP10M))),
        #                labels = scales::number_format(accuracy = 1,
         #                                              decimal.mark = '.')) +
     scale_y_continuous(position = "right", 
                        breaks = c(-4, 0, 4),
                        limits = c(-4, 4),
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
ggsave(plot, filename = paste0("/isdata/PBgrp/xpj980/figures/miRNAtargets/RSG2iso5.svg"), device = "svg", units = "cm", width = 18, height = 18)

####### LOOP over transcript_IDs ######

loop_vector <- unique(gtf_sub$transcript_id)

test_loop_vector <- loop_vector[1:10]

for (i in loop_vector) {
  a_gene <- subset(gtf_sub, gtf_sub$transcript_id == i)
  a_gene_bw <- subsetByOverlaps(gathered_bw_gr, a_gene, ignore.strand = T)
  a_bw_df <- data.frame(a_gene_bw) %>% 
    dplyr::filter(., genotype %in% c("Col.0.WT", "ski2.5", "GKpelota", "SAILpelota")) 
  a_gene_df <- data.frame(a_gene) %>%
    mutate(mean_RP10M=0,
           genotype="gene model")
  a_target_midpoint <- subset(target_midpoint_df, target_midpoint_df$gene_id == unique(a_gene_df$gene_id))
  
  RP10M <- a_bw_df %>% 
    group_by(... = genotype, replicate) %>% 
    left_join(., total_alignments, by = c("genotype" = "genotype", "replicate" = "replicate"), copy = F) %>%
    ungroup() %>%
    mutate(., RP10M=10000000*.$score/.$totalread) 

  new <- RP10M %>% complete(nesting(seqnames, start, end, width, strand, genotype, size), replicate, fill = list(RP10M = 0)) %>% 
    filter(., !(genotype == "SAILpelota" & replicate %in% c("R2", "R3")))
    
  
  meanRP10M <- new %>% 
    group_by(genotype, size, strand, start, end) %>% 
    summarise("mean_RP10M" = mean(RP10M)) %>% 
    ungroup() 
  
  # filter away low reads 
  if (meanRP10M$mean_RP10M %>% abs() %>% max() < 2) {next} 
  
  # stuff for plotting 
  a_gene_df$genotype <- factor(a_gene_df$genotype, levels = c("Col.0.WT", "ski2.5", "GKpelota", "SAILpelota", "gene model"))
  meanRP10M$genotype <- factor(meanRP10M$genotype, levels = c("Col.0.WT", "ski2.5", "GKpelota", "SAILpelota", "gene model"))
  
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
  
  ggsave(plot, filename = paste0("/isdata/PBgrp/xpj980/datascratch/2019_exosome_mutants/06.target_plots/", filename), device = "svg", units = "cm", width = 18, height = 18)
}  


