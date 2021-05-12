#load libraries
library(tidyverse)
library(data.table)

#setwd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#read nominated loci
Nominated_loci <- read.table("example_data/GWAS_input.bim", sep = "\t", header = F)

# Reformat nominated loci to avinput format (using BIM file)
colnames(Nominated_loci) <- c("chr", "ID", "cM", "position", "A1", "A2")
avinput <- Nominated_loci %>%
  select(hg38_chr = chr, 
         chr_start = position, 
         chr_end = position, 
         A1, 
         A2)

write.table(avinput, "90_loci.avinput", sep = "\t", quote = F, row.names=F, col.names=F)

#annotate using ANNOVAR
input_file <- "90_loci.avinput"
annovar_folder <- "/Users/oswaldo/Documents/annovar/"
build <- "-buildver hg38"
output_file <- "-out 90_loci.annoavinput"
arguments <- "-dot2underline -remove -protocol refGene,knownGene,ensGene,gnomad211_exome,clinvar_20160302,dbnsfp30a -arg '-splicing_threshold=6','-splicing_threshold=6','-splicing_threshold=6',,, -operation g,g,g,f,f,f -otherinfo"
cmd <- paste("perl", paste0(annovar_folder,"table_annovar.pl"), input_file, paste0(annovar_folder,"humandb"), build, output_file, arguments)
system(cmd)


