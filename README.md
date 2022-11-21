<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/github_username/repo_name">
    <img src="images/GP2_logo.png" alt="Logo" width="300" height="70">
  </a>

<h3 align="center">Variant interaction pipeline</h3>

  <p align="center">
    One of the projects from the 2021 GP2/IPDGC Hackathon. The related manuscript can be found on [biorxiv](https://www.biorxiv.org/content/10.1101/2022.05.04.490670v1) 
    <br />
    Contributers: Jeff Kim, Sara Bandres Ciga, Geoffrey Rateau, Catherine Storm, Rami Al-Ouran, Oswaldo Lorenzo Betancor
    <br />
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#quick-description">Quick Description</a></li>
        <li><a href="#background/motivation">Background/motivation</a></li>
        <li><a href="#workflow-summary">Workflow Summary</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

![Project Screen Shot][project-screenshot]

### Quick Description

The goal for this project was to develop a pipeline that investigates variant interaction.

### Background/motivation

Genome-wide association studies (GWAS) have identified many common variants associated with complex diseases like Parkinson's disease. As GWAS sample sizes continue to grow, genetic variants with weaker effects on the phenotype are being discovered. It has become clear that GWAS-identified variants do not fully account for the heritability of many common diseases ([Manolio et al. 2009](https://www.nature.com/articles/nature08494)). Indeed, genetics of complex traits can not be attributed to single susceptibility variants identified with GWAS.

Epistasis has been proposed to underlie some of this "missing heritability" identified in the post-GWAS era, though this topic remains controversial ([Ritchie et al. 2015](https://pubmed.ncbi.nlm.nih.gov/25403525/); [Ritchie and Van Steen 2018](https://pubmed.ncbi.nlm.nih.gov/29862246/); [Bandres-Ciga et al. 2020](https://pubmed.ncbi.nlm.nih.gov/31991247/)). Epistasis occurs when there is an interaction between genetic variants. For example, an allele an locus 1 may prevent an alleles at locus 2 from manifesting ([Bateson 1909](https://scholar.google.com/scholar?q=Bateson,+W.+(1909)+Mendel%27s+Principles+of+Heredity.+Cambridge+University+Press,+Cambridge.)), a phenomenon known as masking.

The concept of epistasis has been extensively studied in model systems, but its importance in humans continues to be a matter of debate ([Carlborg and Haley 2004](https://pubmed.ncbi.nlm.nih.gov/15266344/); [Ritchie and Van Steen 2018](https://pubmed.ncbi.nlm.nih.gov/29862246/)). One of the reasons for this discussion is that there is a no of consensus on the methods used to perform epistasis analyses and replication studies are required. 

The current project is a variant interaction pipeline built as part of the IPDGC x GP2 Hackathon 2021 Project to start a collaborative environment to generate a user friendly tool to investigate epistasis between pairs of variants provided by the user.
   

### Workflow Summary

1. Variant annotation with ANNOVAR using chromosome, position, REF and ALT alleles from the BIM file. 
2. Data harmonization to ensure that risk allele is being treated as such for every single variant.
3. Filtering step according to minor allele frequency (MAF) and functional category. Eventually, filtering step will include also p-value from GWAS summary statistics.
4. Pair to pair interaction analysis using glm() and adjusted for age, sex, PC1, PC2, PC3, PC4 and PC5.
5. Results plot.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

* [PLINK 1.9](https://www.cog-genomics.org/plink/)
* [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) 
* [R 3.6](https://www.r-project.org/)

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/ipdgc/GP2-Variant-Interaction-Pipeline.git
   ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Usage

As an example, we test for epistasis between candidate variants identified by the most recent GWAS meta-analysis in European-ancestry Parkinson's disease patients ([Nalls et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31701892/)) and test for epistasis between these variants in openly-available 1000 data from the Genomes Project. However any variants can be used.

Input files:
* List of variants in PLINK binary format.
* Co-variates file in tab delimited format.

Define some variables here and run the variant interaction script:
e.g.

```
export ANNOVAR_FOLDER="file path to annovar application"
export PLINK_FOLDER="file path to plink folder"

Rscript final_variant_interaction_script.R
```

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* [Nalls et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31701892/)


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[project-screenshot]: images/project_screenshot.png
