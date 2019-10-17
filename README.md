# DRIAD: Drug Repurposing in Alzheimer's Disease

## Funding
We gratefully acknowledge support by NIA grant R01 AG058063: Harnessing Diverse BioInformatic Approaches to Repurpose Drugs for Alzheimers Disease.

## Installation

The R package can be installed directly from GitHub using

``` r
if(!require(devtools)) install.packages("devtools")
devtools::install_github( "labsyspharm/DRIAD" )
```

Once installed, the package can be loaded using the standard `library()` command:
``` r
library( DRIAD)
```

## Wrangling AMP-AD datasets

The package provides two functions to wrangle AMP-AD datasets, one for ROSMAP and another for MSBB. Prior to wrangling the datasets, please ensure that you have sufficient access to these data. Data usage agreements can be found at the following URLs:

* ROSMAP: https://www.synapse.org/#!AccessRequirements:ID=syn3219045&TYPE=ENTITY
* MSBB: https://www.synapse.org/#!AccessRequirements:ID=syn3159438&TYPE=ENTITY

Once the data usage agreements have been accepted, each function can be called with a directory path where the wrangled data should live:

``` r
fnROSMAP <- wrangleROSMAP("~/data/ROSMAP")
fnMSBB   <- wrangleMSBB("~/data/MSBB")
```

Both functions return the corresponding filenames of wrangled datasets in tab-delimited format. In case of MSBB, there are four such files: one for each region of the brain profiled by the project.

``` r
fnROSMAP
# [1] "~/data/ROSMAP/rosmap.tsv.gz"

fnMSBB
# [1] "~/data/MSBB/msbb10.tsv.gz" "~/data/MSBB/msbb22.tsv.gz"
# [3] "~/data/MSBB/msbb36.tsv.gz" "~/data/MSBB/msbb44.tsv.gz"
```

## Setting up a prediction task



## Evaluating gene sets

Individual gene sets can be assembled by-hand or parsed from a .gmt file using the provided `read_gmt()`. For example, here's the gene set proposed by Mariet Allen as the most differentially-expressed contrast of the temporal cortex in AD vs. non-AD subjects ([PMID: 29107053](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5866744/), Suppl. Table 3, genes with Bonferroni-corrected p-value < 0.05):

``` r
gs <- c("ADAMTS2", "ADCYAP1", "AFG3L1", "ANGEL2", "ANGPT1", "ART3", 
        "BOLA3", "C1orf89", "C21orf119", "C7orf52", "C9orf24", "COL21A1", 
        "CREB3L1", "CRH", "DACT1", "EFCAB5", "EMG1", "EP400", "FLJ10986", 
        "FNIP1", "FRMPD2", "FUNDC1", "GNA15", "GPATCH8", "HES4", "HSPB3", 
        "JAK2", "KATNAL2", "KHDC1", "KRT5", "LRRC17", "LSM14A", "LYST", 
        "MAS1", "MAX", "MLLT10", "MRPL27", "MRPS18A", "MSC", "MTMR12", 
        "NDUFA13", "NDUFA8", "NIN", "NOLA2", "PIK3C2B", "PLAC8L1", "POMC", 
        "PPEF1", "PPL", "PRKAB2", "PRKG2", "PROC", "PSCD1", "RAPGEF4", 
        "RERG", "RPH3A", "SERTAD4", "SGCG", "SLN", "SLTM", "SNAI2", "SST", 
        "TFPT", "TXNDC17", "USP47", "VWC2", "XRN1", "ZBTB5")
 ```
 
