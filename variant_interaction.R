#load libraries

suppressMessages( library(tidyverse) )
suppressMessages( library(data.table) )
suppressMessages( library(corrplot) )
cat('\nFinsihed loading libraries\n')
#call script from command line as: Rscript variant_interaction.R [plinkFilePrefix] [snpFile] [outFilePrefix] [covarFile]
#ex:Rscript variant_interaction.R GWAS_input chr1_4SNPs.txt covariates.txt Test3 
args = commandArgs(trailingOnly = T)
cat('\nSetting varaibles\n')
#set up the plink files/dataset prefix. Should have .bed, .bim, amd .fam files
plinkFilePrefix=args[1]
#name of file storing the SNPs to be extracted
#snpExtractFile='chr1_2SNPs.txt'
snpExtractFile=args[2]
if(file.exists(snpExtractFile)==F)
{
  cat('snpExtractFile file not found\n')
  q()
}

#get the covariate file name
covariateFile=args[3]

#prefix of output file 
outFilePrefix=args[4]


cat('\nUser arguments are:\n')
cat('plink file prefix:', plinkFilePrefix, '\n')
cat('snpExtract file:', snpExtractFile, '\n')
cat('covariate file:', covariateFile, '\n')
cat('output file prefix:', outFilePrefix, '\n')

#########################
########Methods##########
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
########End of Methods##########
#########################


#read the covar file
covar=fread(covariateFile)
head(covar)

#path to plink cmd
plinkCmd='/Users/rami/Documents/tools/plink/plink_mac_20210416/plink'
if(file.exists(plinkCmd)==F)
{
  cat('plink cmd not found\n')
  q()
}


#set up the plink command
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
cat('Number of input variants:', numVarsInput, '\n')

#make a vector of the input variants
varNameVec=colnames(geno)[7:ncol(geno)]
varNameVec

#merge the geno and covariate file
final_covar=merge(covar, geno, by.x = "FID", by.y = "FID")
head(final_covar)

#scale to deal with zero values
final_covar[,20:(19+numVarsInput)] = final_covar[,20:(19+numVarsInput)] +1

#colnames(final_covar)
#get the 2 pair combinations of variants
varPairList=combn(colnames(final_covar)[20:(20+numVarsInput-1)], 2, simplify = F)
varPairList

#call the regression model per each two variants
resList=lapply(varPairList, callRegress, covarDf=final_covar, inPhenoColumn='PHENO')
#length(resList)
#resList

#initialize matrix to fill in p-values for plotting
matDf = data.frame(matrix(0,ncol=numVarsInput,nrow=numVarsInput))
#matDf
colnames(matDf)=varNameVec
rownames(matDf)=varNameVec
#matDf

#loop thru the list and write the data frames ; one file per df name each output with the two varaints used
for (num in 1:length(resList))
{
  output_tab = resList[[num]]
  outTabName = as.character(paste(outFilePrefix,"_chr",sub(":","_",output_tab[9,1]),"_chr",sub(":","_",output_tab[10,1]),".tab", sep=""))
  write.table(output_tab, file=outTabName, sep='\t', row.names = F, quote = F)
  #fill the df
  matDf[output_tab[9,1],output_tab[10,1]]=output_tab[11,5]
  matDf[output_tab[10,1], output_tab[9,1]]=output_tab[11,5]
}


#plot
matInput=as.matrix(matDf)
fileName=str_glue('{outFilePrefix}_heatmap.png')
fileName
png(fileName)
corrplot(matInput, is.corr = FALSE, type="upper")
dev.off()


