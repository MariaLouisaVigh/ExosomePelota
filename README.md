# ExosomePelota
R code for figures in Exosome/Pelota article

#### Small RNA libraries were mapped with STAR with the following settings: 
--outFilterMultimapNmax 20 --alignIntronMax 1 --outFilterMismatchNmax 1 --outFilterMismatchNoverLmax 0 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 15 --alignSJDBoverhangMin 100 --outFileNamePrefix $file --outMultimapperOrder Random 

#### Sam files were quantified with featureCounts:  
-a "/binf-isilon/PBgrp/xpj980/TAIR/exonTEintergenicmiRNA.gff3" -F GFF -g ID -t whole_shabang -s 0 -T 48 --minOverlap 3 -M -f -O --largestOverlap -o featureCountsaraport *.sam

#### Sam files were furthermore converted to bam and divided by size for miRNA target plots using samtools, awk and wigToBigWig
parallel "samtools view -H {2} > {2.}.{1}.sam" ::: {15..35} ::: *.out.sam ;
parallel "awk 'length(\$10) ~ /{1}/ {print\$0}' {2} >> {2.}.{1}.sam" ::: {15..35} ::: *.out.sam ;

parallel "samtools view -b {} -o {.}.bam" ::: *sam

parallel "samtools sort {} -o {.}.sorted.bam" ::: *bam

parallel "samtools index -b {}" ::: *sorted.bam 

parallel "bam2wig.py -i {} -s /isdata/PBgrp/kld287/scratch/TAIR/MyMess/Chr.size -o {.} -d '++,--'" ::: *bam 

parallel "/isdata/PBgrp/kld287/scratch/my_tools/UCSC_binaries/wigToBigWig -clip {} /isdata/PBgrp/kld287/scratch/TAIR/MyMess/Chr.size {.}.bw" ::: *wig
 
