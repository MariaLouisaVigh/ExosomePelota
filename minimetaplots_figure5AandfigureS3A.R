#### mini-meta target plot #### 
## make it work with ski2-5 21 targets in Figure2A - use it on other datasets afterwards ##

setwd("/binf-isilon/PBgrp/xpj980/datascratch/2015_exosome_mutants/03.STAR_mapping/sam_files")
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(ggrepel)
library(rtracklayer) 
library(reshape2)
library(wesanderson)
library(GenomicRanges)
library(GenomicAlignments)
library(readxl)
library(magrittr)
library(ensembldb)

load("/binf-isilon/PBgrp/xpj980/R_scripts_collected/2015_gatheredDESeqresults.RData")

#### miRNA target midpoint ####
#(both genomic and transcript coordinates)

miRNA_targets <- openxlsx::read.xlsx("/binf-isilon/PBgrp/xpj980/TAIR/miRNAtargetlist_Maria_degradome_done.xlsx", 1) %>%
  as_tibble() %>% 
  dplyr::filter(!duplicated(.$target_id)) %>% 
  dplyr::rename(gene = target_id)

targets_GR <- makeGRangesFromDataFrame(df = miRNA_targets, keep.extra.columns = T, seqnames.field = "chr", start.field = "start", end.field = "stop", strand.field = "strand")

target_midpoint <- targets_GR %>% 
  promoters(., downstream = 11) %>% 
  resize(., width = 1, fix = "end") %>% 
  as.data.frame()

target_midpoint_GR <- makeGRangesFromDataFrame(target_midpoint, keep.extra.columns = T, seqnames.field = "seqnames", start.field = "start", end.field = "end", strand.field = "strand")

## mRNA coordinates of midpoints ##
edb <- EnsDb("/binf-isilon/PBgrp/xpj980/TAIR/EnsDb.Athaliana.v94.sqlite")

#transcriptmidpoint <- genomeToTranscript(target_midpoint_GR, edb) %>% 
  #endoapply(., function(x) x[1]) %>% 
  #unlist() %>% 
  #as.data.frame()
#save(transcriptmidpoint, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/targetmidpointsTranscriptcoordinates.RData")
load(file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/targetmidpointsTranscriptcoordinates.RData")
# still a few 6 missing but I can maybe add them manually later 

filtered_midpoints <- transcriptmidpoint %>% 
  group_by(names) %>% 
  dplyr::filter(n() == 1, !is.na(tx_id))

#### genome coordinates of all miRNA targets converted to transcriptome coordinates ####
# mygtf <- import("/binf-isilon/PBgrp/xpj980/TAIR/Arabidopsis_thaliana.TAIR10.44.gtf")
Carlottagtf <- import("/binf-isilon/PBgrp/qbp693/genome_data/Araport11_GTF_genes_transposons.Mar202021_relev.gtf")

gtf_sub_Carlotta <- subset(x = Carlottagtf, 
                            Carlottagtf$gene_id %in% target_midpoint$gene &
                              Carlottagtf$type == "exon")

gtf_sub <- gtf_sub_Carlotta[str_detect(gtf_sub_Carlotta$transcript_id, pattern = ".1$"), ] 


#gtf_sub <- subset(x = mygtf, 
 #                 mygtf$gene_id %in% target_midpoint$gene &
  #                  mygtf$type == "exon")

# load(file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/miRNAtargetsTranscriptcoordinates.RData")

# test <- genomeToTranscript(gtf_sub, edb) %>% 
#   endoapply(., function(x) x[1]) %>% 
#   unlist() %>% 
#   as.data.frame()

#test_filtered <- dplyr::filter(test, !duplicated(test$exon_id))

# test_Carlotta <- genomeToTranscript(gtf_sub, edb) %>% 
#   endoapply(., function(x) x[1]) %>% 
#   unlist() %>% 
#   as.data.frame()
# 
# save(test_Carlotta, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/miRNAtargetsTranscriptcoordinates.RData")

load(file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/miRNAtargetsTranscriptcoordinates.RData")


test_filtered <- dplyr::filter(test_Carlotta, !duplicated(test_Carlotta$exon_id)) %>% glimpse()

#save(test_filtered, file = "/binf-isilon/PBgrp/xpj980/R_scripts_collected/miRNAtargetsTranscriptcoordinates.RData")

##### expand the transcript intervals #####

expand_transcript <- function(transcript_id, miRNAtargetsTranscriptcoordinates = test_filtered) {
  
  expanded <- function(x) {  
    
    if (x[["seq_strand"]]=="-") {
      rowdf <- data.frame(seq_name = rep(x$seq_name),
                          "transcript" = seq(x$start, x$end), 
                          "genomic" = seq(x$seq_end, x$seq_start), 
                          name = rep(x$names), 
                          exon_id = rep(x$exon_id), 
                          seq_strand = rep(x$seq_strand))
    } else {
      rowdf <- data.frame(seq_name = rep(x$seq_name),
                          "transcript" = seq(x$start, x$end), 
                          "genomic" = seq(x$seq_start, x$seq_end), 
                          name = rep(x$names), 
                          exon_id = rep(x$exon_id), 
                          seq_strand = rep(x$seq_strand))
    }
    
    return(rowdf)
  }
  
  a_trx_df <- dplyr::filter(miRNAtargetsTranscriptcoordinates, tx_id == transcript_id)
  
  df_list <- lapply(1:nrow(a_trx_df), function(row_N) {
    expanded(a_trx_df[row_N,])
  }) %>%
    purrr::reduce(bind_rows)
  
  return(df_list)
  
}


##### import featureCounts summary #######
## import featureCounts summary for total reads per genotype to merge with tidyMetaProfile later 
featureCounts_sum <- read.table("/binf-isilon/PBgrp/xpj980/datascratch/2015_exosome_mutants/03.STAR_mapping/sam_files/featureCounts.summary", header = T, row.names = 1)  %>%
  as.data.frame()
new_names <- names(featureCounts_sum) %>% str_remove("\\.Aligned\\.out\\.sam") %>% 
  str_remove(pattern = "set1.") %>%
  str_replace(pattern = ".R([12])$", "_R\\1") %>%
  str_replace(pattern = "cer7.3", replacement = "rrp45b") 

featureCounts_renamed <- featureCounts_sum
names(featureCounts_renamed) <- new_names

total_alignments <- featureCounts_renamed[1,] %>%
  gather(key = "genotype", value = "totalread") %>% 
  separate(genotype, c("genotype", "replicate"), sep="_")
  
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

##### REMEMBER TO PUT SOME miRNA target midpoints manually!!! #### 
to_filter <- dplyr::filter(transcriptmidpoint, is.na(tx_id))
manual4 <- dplyr::filter(target_midpoint, start %in% to_filter$seq_start) %>%
  dplyr::filter(., target_type != "pseudogene" & target_biotype != "other_rna") %>% 
  pull(gene) %>% 
  str_c(".1") 

first <- expand_transcript(transcript_id = "AT1G07770.1") 
extract <- dplyr::filter(target_midpoint,target_midpoint$gene == "AT1G07770")[["start"]]
first <- first[first$genomic == extract,] %>% 
  dplyr::rename(names = "name")

second <- expand_transcript(transcript_id = "AT2G34710.1") 
extract <- dplyr::filter(target_midpoint,target_midpoint$gene == "AT2G34710")[["start"]]+1
second <- second[second$genomic == extract,] %>% 
  dplyr::rename(names = "name")

third <- expand_transcript(transcript_id = "AT5G18100.1") 
extract <- dplyr::filter(target_midpoint,target_midpoint$gene == "AT5G18100")[["start"]]-10
third <- third[third$genomic == extract,] %>% 
  dplyr::rename(names = "name")

fourth <- expand_transcript(transcript_id = "AT5G60690.1") 
extract <- dplyr::filter(target_midpoint,target_midpoint$gene == "AT5G60690")[["start"]]-1
fourth <- fourth[fourth$genomic == extract,] %>% 
  dplyr::rename(names = "name")

selected_transcriptmid <- filtered_midpoints %>% dplyr::select(seq_name, start, seq_start, names, exon_id, seq_strand) %>% 
  dplyr::rename(transcript = "start", genomic = "seq_start") 

complete_midpoint_transcript <- rbind(selected_transcriptmid, first, second, third, fourth)


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

get_meta_cov <- function(selected_mutant = "ski3",
                         selected_transcript = "AT1G48410.1") {
  
  selected_mut_bw <- gathered_bw %>% 
    dplyr::filter(., genotype == selected_mutant) %>% 
    as("GRanges")
  
  transcriptrange  <- subset(gtf_sub, transcript_id == selected_transcript) 
  
  selected_transcript_df <- subsetByOverlaps(selected_mut_bw, transcriptrange, ignore.strand = T) %>% 
    as.data.frame()
  
  expanded_transcript <- expand_transcript(selected_transcript)
  
 #midpoint <- expanded_transcript[expanded_transcript$genomic == dplyr::filter(complete_midpoint_transcript, names == selected_transcript)[["seq_start"]], ]
  midpoint2 <- complete_midpoint_transcript[complete_midpoint_transcript$names == selected_transcript, ]$transcript 
 #midpoint <- dplyr::filter(filtered_midpoints, names == selected_transcript)
  # dplyr::filter(filtered_midpoints, names == selected_transcript)
  # target_midpoint_GR %>%
  #   subset(target_name=="AGO1") %>%
  #   genomeToTranscript( db = edb)  
      
  joined <- left_join(selected_transcript_df, expanded_transcript, 
                      by = c("seqnames" = "seq_name", "start"="genomic")) %>% 
    mutate(CS = midpoint2, 
           distance = transcript - CS)
  
  RP10M <-  joined %>% 
    dplyr::left_join(., total_alignments, by = c("genotype" = "genotype", "replicate" = "replicate"), copy = F) %>%
    mutate(., RP10M=10000000*.$score/.$totalread) 
  
  new <- RP10M %>% complete(nesting(distance, strand, genotype, size), replicate, fill = list(RP10M = 0)) %>% 
    dplyr::select(distance, strand, genotype, size, replicate, RP10M)
  
  meanRP10M <- new %>% 
    group_by(genotype, size, strand, distance) %>% 
    summarise("mean_RP10M" = mean(RP10M)) %>% 
    ungroup() %>%
    mutate(tx_id =selected_transcript, 
           meanRP10M_final = mean_RP10M/sum(mean_RP10M))  
  
  return(meanRP10M)
}

ski2_list <- left_joined[left_joined$genotype=="ski2.5",] %>% 
  dplyr::filter(padj<0.05) %>% 
  dplyr::filter(log2FoldChange > 0) %>% 
  dplyr::filter(!is.na(miRNA)) %>%
  pull(gene) 
ski2 <- ski2_list[str_c(ski2_list,".1") %in% complete_midpoint_transcript$names] %>% str_c(".1")
rrp45_list <- left_joined[left_joined$genotype=="rrp45b",] %>% 
  dplyr::filter(padj<0.05) %>% 
  dplyr::filter(log2FoldChange > 0) %>% 
  dplyr::filter(!is.na(miRNA)) %>%
  pull(gene) 
rrp45 <- rrp45_list[str_c(rrp45_list,".1") %in% complete_midpoint_transcript$names] %>% str_c(".1")
ski3_list <- left_joined[left_joined$genotype=="ski3",] %>% 
  dplyr::filter(padj<0.05) %>% 
  dplyr::filter(log2FoldChange > 0) %>% 
  dplyr::filter(!is.na(miRNA)) %>%
  pull(gene) 
ski3 <- ski3_list[str_c(ski3_list,".1") %in% complete_midpoint_transcript$names] %>% str_c(".1")


listone <- lapply(ski3, get_meta_cov, selected_mutant = "ski3")

plot <- purrr::reduce(listone, bind_rows) %>%
  ggplot(aes(x=distance, y = meanRP10M_final, fill = size)) +
  geom_bar(stat = "summary", fun="sum", position = "stack", aes(group=interaction(strand, size))) +
  scale_fill_brewer(palette = "Set1") +
  cowplot::theme_cowplot() +
  ggtitle(label = "minimetaplot",paste0(selected_mutant)) + 
  labs(y= "mean(Reads per 10 million)", x = "absolute distance from cleavage site (nt)", fill= "nt size") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.3) +
  geom_hline(yintercept = 0, size = 0.1) +
  theme(strip.text.y = element_text(angle = 180, hjust = 1), 
        strip.background = element_blank(), 
        plot.title = element_text(size = 24, 
                                  face = "bold",
                                  color = "#22292F", 
                                  margin = margin(b = 8)), 
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

ggsave(plot, filename = paste0("/isdata/PBgrp/xpj980/figures/minimetaski32015.svg"), device = "svg", units = "cm", width = 18, height = 18)


### now do minimeta plots for ski2, hen2, rrp4 in dataset 2019 ###

# load deseq genes #
load("/binf-isilon/PBgrp/xpj980/R_scripts_collected/2019exoDESEQrrp4.RData")
load("/binf-isilon/PBgrp/xpj980/R_scripts_collected/2019exoDESEQpelota.RData")
left_joined <- left_join(x = gathered_DESeqresults_rrp4, y = miRNA_targets, by = "gene", copy = FALSE)
left_joined <- left_join(x = gathered_DESeqresults_pelota, y = miRNA_targets, by = "gene", copy = FALSE)
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

####### load bw files and rename them ######

## prepare names for bw forward and reverse files 
files_names_forward <- list.files("/binf-isilon/PBgrp/xpj980/datascratch/2019_exosome_mutants/04.STAR_mapping/sam_files/by_size/bw/", full.names = T, pattern = "Forward")
good_names_fw <- files_names_forward %>%
  str_remove(pattern = ".*by_size/bw//") %>%
  str_remove(pattern = ".Aligned.out") %>%
  str_remove(pattern = ".sorted.Forward.bw") %>%
  str_replace(pattern = "\\.", replacement = "_") %>%
  str_replace_all(pattern = "-", replacement = "\\.") 
names(files_names_forward) <- good_names_fw

files_names_reverse <- list.files("/binf-isilon/PBgrp/xpj980/datascratch/2019_exosome_mutants/04.STAR_mapping/sam_files/by_size/bw/", full.names = T, pattern = "Reverse")
good_names_rv <- files_names_reverse %>%
  str_remove(pattern = ".*by_size/bw//") %>%
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

get_meta_cov <- function(selected_mutant = "rrp4",
                         selected_transcript = "AT1G48410.1") {
  
  selected_mut_bw <- gathered_bw %>% 
    dplyr::filter(., genotype == selected_mutant) %>% 
    as("GRanges")
  
  transcriptrange  <- subset(gtf_sub, transcript_id == selected_transcript) 
  
  selected_transcript_df <- subsetByOverlaps(selected_mut_bw, transcriptrange, ignore.strand = T) %>% 
    as.data.frame()
  
  expanded_transcript <- expand_transcript(selected_transcript)
  
  #midpoint <- expanded_transcript[expanded_transcript$genomic == dplyr::filter(complete_midpoint_transcript, names == selected_transcript)[["seq_start"]], ]
  midpoint2 <- complete_midpoint_transcript[complete_midpoint_transcript$names == selected_transcript, ]$transcript 
  #midpoint <- dplyr::filter(filtered_midpoints, names == selected_transcript)
  # dplyr::filter(filtered_midpoints, names == selected_transcript)
  # target_midpoint_GR %>%
  #   subset(target_name=="AGO1") %>%
  #   genomeToTranscript( db = edb)  
  
  joined <- left_join(selected_transcript_df, expanded_transcript, 
                      by = c("seqnames" = "seq_name", "start"="genomic")) %>% 
    mutate(CS = midpoint2, 
           distance = transcript - CS)
  
  RP10M <-  joined %>% 
    dplyr::left_join(., total_alignments, by = c("genotype" = "genotype", "replicate" = "replicate"), copy = F) %>%
    mutate(., RP10M=10000000*.$score/.$totalread) 
  
  new <- RP10M %>% complete(nesting(distance, strand, genotype, size), replicate, fill = list(RP10M = 0)) %>% 
    dplyr::select(distance, strand, genotype, size, replicate, RP10M)
  
  meanRP10M <- new %>% 
    group_by(genotype, size, strand, distance) %>% 
    summarise("mean_RP10M" = mean(RP10M)) %>% 
    ungroup() %>%
    mutate(tx_id =selected_transcript, 
           meanRP10M_final = mean_RP10M/sum(mean_RP10M))  
  
  return(meanRP10M)
}

ski2_list <- left_joined[left_joined$genotype=="ski2.5",] %>% 
  dplyr::filter(padj<0.05) %>% 
  dplyr::filter(log2FoldChange > 0) %>% 
  dplyr::filter(!is.na(miRNA)) %>%
  pull(gene) 
ski2 <- ski2_list[str_c(ski2_list,".1") %in% complete_midpoint_transcript$names] %>% str_c(".1")
rrp4_list <- left_joined[left_joined$genotype=="rrp4",] %>% 
  dplyr::filter(padj<0.05) %>% 
  dplyr::filter(log2FoldChange > 0) %>% 
  dplyr::filter(!is.na(miRNA)) %>%
  pull(gene) 
rrp4 <- rrp4_list[str_c(rrp4_list,".1") %in% complete_midpoint_transcript$names] %>% str_c(".1")
hen2_list <- left_joined[left_joined$genotype=="hen2.5",] %>% 
  dplyr::filter(padj<0.05) %>% 
  dplyr::filter(log2FoldChange > 0) %>% 
  dplyr::filter(!is.na(miRNA)) %>%
  pull(gene) 
hen2 <- hen2_list[str_c(hen2_list,".1") %in% complete_midpoint_transcript$names] %>% str_c(".1")
pel11_list <- left_joined[left_joined$genotype=="SAILpelota",] %>% 
  dplyr::filter(padj<0.05) %>% 
  dplyr::filter(log2FoldChange > 0) %>% 
  dplyr::filter(!is.na(miRNA)) %>%
  pull(gene) 
pel11 <- pel11_list[str_c(pel11_list,".1") %in% complete_midpoint_transcript$names] %>% str_c(".1")
pel12_list <- left_joined[left_joined$genotype=="GKpelota",] %>% 
  dplyr::filter(padj<0.05) %>% 
  dplyr::filter(log2FoldChange > 0) %>% 
  dplyr::filter(!is.na(miRNA)) %>%
  pull(gene) 
pel12 <- pel12_list[str_c(pel12_list,".1") %in% complete_midpoint_transcript$names] %>% str_c(".1")

listone <- lapply(pel12, get_meta_cov, selected_mutant = "GKpelota")

plot <- purrr::reduce(listone, bind_rows) %>%
  ggplot(aes(x=distance, y = meanRP10M_final, fill = size)) +
  geom_bar(stat = "summary", fun="sum", position = "stack", aes(group=interaction(strand, size))) +
  scale_fill_brewer(palette = "Set1") +
  cowplot::theme_cowplot() +
  ggtitle(label = "minimetaplot",paste0(selected_mutant)) + 
  labs(y= "mean(Reads per 10 million)", x = "absolute distance from cleavage site (nt)", fill= "nt size") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.3) +
  geom_hline(yintercept = 0, size = 0.1) +
  theme(strip.text.y = element_text(angle = 180, hjust = 1), 
        strip.background = element_blank(), 
        plot.title = element_text(size = 24, 
                                  face = "bold",
                                  color = "#22292F", 
                                  margin = margin(b = 8)), 
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)))

ggsave(plot, filename = paste0("/isdata/PBgrp/xpj980/figures/minimetaGKpel122019.svg"), device = "svg", units = "cm", width = 18, height = 18)


