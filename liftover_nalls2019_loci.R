
library(data.table)
library(tidyverse)

# make bed file for crossmap
dat <- fread('nalls2019_nominated_loci.txt')
dat$pos <- gsub(',', '', dat$BP)
dat$pos <- as.numeric(dat$pos)
dat$chr <- paste('chr', dat$CHR, sep='')
dat$posMinus1 <- as.integer(dat$pos-1)

# can i get CrossMap to grab this file without writing it out like this?
fwrite(dat[,c('chr','posMinus1','pos')], 'nalls2019_nominated_loci_bed4liftover.txt', sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

# download hg19tohg38 liftover chain file and run CrossMap
cmd='wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

.local/bin/CrossMap.py bed hg19ToHg38.over.chain.gz nalls2019_nominated_loci_bed4liftover.txt > nalls2019_nominated_loci_crossmap_hg19ToHg38.txt'
cat('cmd:', cmd, '\n')
system(cmd)

# edit to make a snp list
hg19_bed <- fread('nalls2019_nominated_loci_crossmap_hg19ToHg38.txt')

hg19_bed$V4 <- NULL
hg19_bed$V2 <- NULL
hg19_bed$V6 <- NULL
names(hg19_bed)[1:4] <- c('hg19_chr', 'hg19_pos', 'hg38_chromosome', 'hg38_position')

hg19_bed$chrpos_hg19 <- paste(hg19_bed$hg19_chr, hg19_bed$hg19_pos, sep=':')
hg19_bed$chrpos_hg38 <- paste(hg19_bed$hg38_chromosome, hg19_bed$hg38_position, sep=':')

dat$chrpos_hg19 <- paste(dat$chr, dat$pos, sep=':')
hg19_bed_with_alleles <- left_join(hg19_bed, dat[,c('chrpos_hg19','Effect allele', 'Other allele', 'Effect allele frequency', 'Beta, all studies','P, all studies')], by='chrpos_hg19')

names(hg19_bed_with_alleles)[7:11] <- c('ea', 'oa', 'eaf', 'beta','p')

# write it out
fwrite(hg19_bed_with_alleles, 'nalls2019_nominated_loci_hg38.txt', sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
fwrite(hg19_bed_with_alleles[,'chrpos_hg38'], 'nalls2019_nominated_loci_hg38_snplist4plink.txt', sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)
