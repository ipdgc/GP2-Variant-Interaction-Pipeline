#load libraries
library(tidyverse)
library(data.table)

#setwd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#read nominated loci
Nominated_loci <- read.table("nalls2019_nominated_loci_hg38.txt", sep = "\t", header = T)

# Reformat nominated loci to avinput format (using lifted nalls2019_nominated_loci_hg38.txt file)
avinput <- Nominated_loci %>%
  select(hg38_chr = hg38_chromosome, 
         chr_start = hg38_position, 
         chr_end = hg38_position, 
         oa, 
         ea)
avinput$hg38_chr <- avinput$hg38_chr <- gsub("^.{0,3}", "", avinput$hg38_chr)

write.table(avinput[1:90,], "90_loci.avinput", sep = "\t", quote = F, row.names=F, col.names=F)

#annotate using ANNOVAR
input_file <- "90_loci.avinput"
annovar_folder <- "/Users/oswaldo/Documents/annovar/"
build <- "-buildver hg38"
output_file <- "-out 90_loci.annoavinput"
arguments <- "-dot2underline -remove -protocol refGene,knownGene,ensGene,gnomad211_exome,clinvar_20160302,dbnsfp30a -arg '-splicing_threshold=6','-splicing_threshold=6','-splicing_threshold=6',,, -operation g,g,g,f,f,f -otherinfo"
cmd <- paste("perl", paste0(annovar_folder,"table_annovar.pl"), input_file, paste0(annovar_folder,"humandb"), build, output_file, arguments)
system(cmd)


