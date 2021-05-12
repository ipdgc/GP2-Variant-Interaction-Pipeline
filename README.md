# GP2-Variant-Interaction-Pipeline

## Contributors
Jeff Kim, Sara Bandres Ciga, Geoffrey Rateau, Catherine Storm, Rami Al-Ouran, Oswaldo Lorenzo Betancor

## Introduction / Motivation
This is a variant interaction pipeline built as part of the IPDGC x GP2 Hackathon 2021 Project. This pipeline can be used to investigate epistasis between pairs of variants.

Genome-wide association studies (GWAS) have identified many common variants associated with complex diseases like Parkinson's disease. As GWAS sample sizes continue to grow, genetic variants with weaker effects on the phenotype are being discovered. It has become clear that GWAS-identified variants do not fully account for the heritability of many common diseases ([Manolio et al. 2009](https://www.nature.com/articles/nature08494)).

Epistasis has been proposed to underlie some of this "missing heritability" identified in the post-GWAS era, though this topic remains controversial ([Ritchie et al. 2015](https://pubmed.ncbi.nlm.nih.gov/25403525/); [Ritchie and Van Steen 2018](https://pubmed.ncbi.nlm.nih.gov/29862246/); [Bandres-Ciga et al. 2020](https://pubmed.ncbi.nlm.nih.gov/31991247/)).Epistasis occurs when there is an interaction between genetic variants. For example, an allele an locus 1 may prevent an alleles at locus 2 from manifesting ([Bateson 1909](https://scholar.google.com/scholar?q=Bateson,+W.+(1909)+Mendel%27s+Principles+of+Heredity.+Cambridge+University+Press,+Cambridge.)), a phenomenon known as masking.

Genetics of complex traits can not be attributed to single susceptibility genes identified with Genome-Wide Association studies (GWAS). The "missing heritability" of complex traits can be further explained with genetic interactions. The concept of epistasis has been extensively studied in model systems, but its importance in humans continues to be a matter of debate. One of the reasons for this discussion is that there is a complete lack of consensus when performing epistasis analyses and replication studies are required.

There is an abundance of methods for the exploration of genetic interactions, including different statistical methods as well as machine learning and data mining techniques, each of them leading to a different result even when analyzing the same dataset.

The current project is a variant interaction pipeline built as part of the IPDGC x GP2 Hackathon 2021 Project to start a collaborative environment to generate a user friendly tool to investigate epistasis between pairs of variants provided by the user.

## Description of files/data


## How to use
Define some variables here and run the variant interaction script:
e.g.

```
export ANNOVAR_FOLDER="file path to annovar application"
export PLINK_FOLDER="file path to plink folder"

Rscript final_variant_interaction_script.R
```

## Example
As an example, we provide data to test for epistasis between candidate variants identified by the most recent GWAS meta-analysis in European-ancestry Parkinson's disease patients ([Nalls et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31701892/)). We test for epistasis between these variants in openly-available 1000 data from the Genomes Project (citation).

