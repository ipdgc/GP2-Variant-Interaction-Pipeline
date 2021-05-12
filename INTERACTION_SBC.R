#load libraries
library(tidyverse)
library(data.table)
library(pheatmap)

#get user line arguments; change this
#args = commandArgs(trailingOnly = TRUE)
#cat('1:', args[1], '\n')

#set up input arguments after connecting with everyone else on the team to final.
#i.e. GWAS input file, also another argument for potentially focusing on or skipping the analysis stratified by annotation types.

#set up the plink files/dataset prefix. Should have .bed, .bim, amd .fam files
plinkFilePrefix='GWAS_input'

#read variants of interest
user_variants <- read.table("plinkFilePrefix", sep = "\t", header = F)

#reformat user variants of interest to avinput format (using BIM file)
colnames(user_variants) <- c("chr", "ID", "cM", "position", "A1", "A2")
avinput <- user_variants %>%
  select(hg38_chr = chr, 
         chr_start = position, 
         chr_end = position, 
         A1, 
         A2)

write.table(avinput, "user_variants.avinput", sep = "\t", quote = F, row.names=F, col.names=F)

#annotate using ANNOVAR
input_file <- "user_variants.avinput"
annovar_folder <- "/Users/annovar/"
build <- "-buildver hg38"
output_file <- "-out user_variants.annoavinput"
arguments <- "-dot2underline -remove -protocol refGene,knownGene,ensGene,gnomad211_exome,clinvar_20160302,dbnsfp30a -arg '-splicing_threshold=6','-splicing_threshold=6','-splicing_threshold=6',,, -operation g,g,g,f,f,f -otherinfo"
cmd <- paste("perl", paste0(annovar_folder,"table_annovar.pl"), input_file, paste0(annovar_folder,"humandb"), build, output_file, arguments)
system(cmd)

#OPTIONAL: subset by class of variants and MAF, these can be input later if you just want to use them as snpExtractFile options
intergenic_SNPs <- subset(output_file, Func_refGene == "intergenic" & AF_nfe >= 0.005 & AF_nfe <= 0.995)
intergenic_SNPs$chr_pos <- paste(intergenic_SNPs$Chr, ":", intergenic_SNPs$Start, sep="")
write.table(intergenic_SNPs, "intergenic_SNPs.csv", quote = F, sep = ",", row.names = F)

intronic_SNPs <- subset(output_file, Func_refGene == "intronic" & AF_nfe >= 0.005 & AF_nfe <= 0.995)
intronic_SNPs$chr_pos <- paste(intronic_SNPs$Chr, ":", intronic_SNPs$Start, sep="")
write.table(intronic_SNPs, "intronic_SNPs.csv", quote = F, sep = ",", row.names = F)

exonic_SNPs <- subset(output_file, Func_refGene == "exonic" & ExonicFunc_refGene = 'nonsynonymous SNV' & AF_nfe >= 0.005 & AF_nfe <= 0.995)
exonic_SNPs$chr_pos <- paste(exonic_SNPs$Chr, ":", exonic_SNPs$Start, sep="")
write.table(exonic_nonsynonymous_SNPs, "exonic_nonsynonymous_SNPs.csv", quote = F, sep = ",", row.names = F)

#name of file storing the SNPs to be extracted

#for loop here that runs through different classes of variants.
class_in =c("intergenic_SNPs","intronic_SNPs","exonic_SNPs")
#start loop.
for(i in class_in){
snpExtractFile=class_in
if(file.exists(snpExtractFile)==F)
{
  cat('snpExtractFile file not found\n')
  q()
}

#prefix of output file 
outFilePrefix='Test2'


#path to plink cmd
plinkCmd='/Users/plink'
if(file.exists(plinkCmd)==F)
{
  cat('plink cmd not found\n')
  q()
}

#read the covariate file
covariateFile='covariates.txt'
covar=fread(covariateFile)
head(covar)

#set up the command
cmd=str_glue('{plinkCmd} --bfile {plinkFilePrefix}  --extract {snpExtractFile} --recode A --out {outFilePrefix}')
cat('cmd:', cmd, '\n')
system(cmd)

#read the raw file
genotypeRawFile=str_glue('{outFilePrefix}.raw')
geno=fread(genotypeRawFile)
head(geno)
ncols=ncol(geno)

#number of varaints
numVarsInput=ncols-6
numVarsInput

varNameVec=colnames(geno)[7:ncol(geno)]
varNameVec

#merge the geno and covariate file
final_covar=merge(covar, geno, by.x = "FID", by.y = "FID")
head(final_covar)

#scale to deal with zero values
final_covar[,20:(19+numVarsInput)] = final_covar[,20:(19+numVarsInput)] +1

colnames(final_covar)
#get the 2 pair combinations
varPairList=combn(colnames(final_covar)[20:(20+numVarsInput-1)], 2, simplify = F)
varPairList


#function to call regression
callRegress = function(inVarList,covarDf, inPhenoColumn)
{
  #inVarList='4_17968811_A 4_77110365_G'
  #covarDf=final_covar
  cat('\nvars are:', inVarList, '\n')
  s=unlist(strsplit(inVarList, ' '))
  firstVar=s[1]
  scndVar=s[2]
  cat('firstVar:', firstVar, ' scndVar:', scndVar, '\n')
  #cat('nrow df:', nrow(covarDf), '\n')
  
  #change column names
  colVec=c()
  for(colName in colnames(covarDf))
  {
    if(colName==firstVar)
    {
      colName='firstVar'
    }else if(colName==scndVar)
    {
      colName='scndVar'
    }
    colVec=c(colVec, colName)
  }
  colnames(covarDf)=colVec
  
  
  #create inetarction term for regression model
  covarDf$interactionTerm = covarDf$firstVar*covarDf$scndVar
  
  #assign the index of the two varaints for re assigning the actual names
  firstVarIndex=9
  scndVarIndex=10
  
  #call the model
  fmla=str_glue('{inPhenoColumn} ~  AGE + SEX + PC1 + PC2 + PC3 + PC4 + PC5 + firstVar + scndVar + interactionTerm')
  fmla
  fmla=as.formula(fmla)
  interactModel = glm(fmla, data = covarDf, family="binomial")
  summaryStats = summary(interactModel)
  statDf=data.frame('Variable'=row.names(coef(summaryStats)), coef(summaryStats))
  statDf
  colnames(statDf)=c('Variable', 'Estimate', 'StdError', 'Zvalue', 'Pvalue')
  #assign back the variant names
  statDf$Variable[firstVarIndex]=firstVar
  statDf$Variable[scndVarIndex]=scndVar
  return(statDf)
}

#call the regression model per each two variants
resList=lapply(varPairList, callRegress, covarDf=final_covar, inPhenoColumn='PHENO')
length(resList)
resList


#loop thru the list and write the data frames ; one file per df
#name each output with the two varaints names

#loop thru the list and write the data frames ; one file per df
#name each output with the two varaints names

matDf = data.frame(matrix(0,ncol=numVarsInput,nrow=numVarsInput))
matDf
colnames(matDf)=varNameVec
rownames(matDf)=varNameVec
matDf


for (num in 1:length(resList))
{
  output_tab <- resList[[num]]
  outTabName <- as.character(paste("chr",sub(":","_",output_tab[9,1]),"_chr",sub(":","_",output_tab[10,1]),".tab", sep=""))
  write.table(output_tab, file=outTabName, sep='\t', row.names = F, quote = F)
  #fill the df
  matDf[output_tab[9,1],output_tab[10,1]]=output_tab[11,5]
  matDf[output_tab[10,1], output_tab[9,1]]=output_tab[11,5]
}

#matDf


#plot
matInput=as.matrix(matDf)
fileName=str_glue('{outFilePrefix}_heatmap.png')
fileName
png(fileName)
corrplot(matInput, is.corr = FALSE, type="upper")
dev.off()
# remember to flag outputs with suffix from class_in list.
}