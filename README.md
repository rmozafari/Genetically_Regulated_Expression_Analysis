# GREX analysis pipeline

## Overview

This repo contains a part of a comprehensive Snakemake pipelines designed for processing genomic data, specifically focusing on enhancing the accuracy of gene expression imputation and identifying gene-trait associations. The pipeline integrates various steps, including data preparation, variant filtering, expression analysis, and the imputation of genetically regulated expression (GReX) using a modified version of the AffiXcan package.

By leveraging transcriptome-wide association studies (TWAS), this project aims to prioritize trait-associated genes by correlating complex traits with the genetically regulated components of gene expression. Instead of relying solely on individual genetic variants, our approach utilizes the sequence composition of proximal gene regulatory regions, represented by transcription factor binding affinities. This innovative methodology has demonstrated improved imputation accuracy compared to traditional models that regress expression directly on genotype.

The pipeline is structured to facilitate reproducible research in genomics, enabling researchers to analyze and interpret genomic data consistently. By applying this pipeline to datasets such as GTEx and the Alzheimer’s Disease Neuroimaging Initiative (ADNI), we aim to uncover novel gene-trait associations and enhance our understanding of the genetic underpinnings of complex traits and diseases. This work not only advances the methodologies used in TWAS but also provides a robust framework for integrating genomic and transcriptomic data, ultimately contributing to the identification of potential therapeutic targets in complex diseases.


## Directory Structure

```
.
├── config
│   └── snakemake_config.yml
├── src
│   ├── MultiAssayExp.constructor.R
│   ├── SumExp.constructor.R
│   └── transpose.pl
│   └── glm.R
│   └── Kfold_CrossValidation.R
└── snakefiles
    └── crossvalidation.smk
```