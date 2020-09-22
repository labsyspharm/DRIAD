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
library(DRIAD)
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

The package considers binary classification tasks that contrast early (A), intermediate (B) and late (C) disease stages, as defined by the Braak staging through neuropathological assessment. To set up a classification task, we need to specify a dataset wrangled by `wrangleROSMAP()` or `wrangleMSBB()` as well as which pair of classes we want to compare. For example, the early-vs-late classification task in MSBB (Brodmann Area 36) can be constructed via

``` r
# Recall that fnMSBB contains four filenames, the third of which points to msbb36.tsv.gz
task1 <- prepateTask( fnMSBB[3], "AC" )
```

Likewise, the intermediate-vs-late classification task in ROSMAP is constructed via

``` r
task2 <- prepareTask( fnROSMAP, "BC" )
```

To evaluate predictor peformance for a given binary classification task, we consider [leave-pair-out cross-validation (LPOCV)](http://proceedings.mlr.press/v8/airola10a/airola10a.pdf), where one sample from each of the two classes is withheld in every cross-validation fold. The samples are age-matched to address the potential bias of age as a confounder in performance estimation. The package can automatically construct the pairing of age-matched samples, using a task definition produced by `prepareTask()`. Continuing the example above,

``` r
pairs1 <- preparePairs( task1 )
# $`1`
# [1] "AMPAD_MSSM_0000003566" "AMPAD_MSSM_0000016747"
# 
# $`2`
# [1] "AMPAD_MSSM_0000004030" "AMPAD_MSSM_0000077958"
# 
# $`3`
# [1] "AMPAD_MSSM_0000004310" "AMPAD_MSSM_0000041763"
# ...and 138 more pairs
```

## Evaluating gene sets

Individual gene sets can be assembled by-hand or parsed from a .gmt file using the provided `read_gmt()`. For example, here's the gene set proposed by Mariet Allen and colleagues as the most differentially-expressed contrast of the temporal cortex in AD vs. non-AD subjects ([PMID: 29107053](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5866744/), Suppl. Table 3, genes with Bonferroni-corrected p-value < 0.05):

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
 
A list of one or more gene sets can be provided to `evalGeneSets()` together with the classification task and the correpsonding pairing of samples for cross-validation. For example, let's give the gene set above a name and evaluate it on our first classification task (A-vs-C in MSBB36) against 100 background sets:

``` r
evalGeneSets( list("Allen, et al."=gs), task1, pairs1, 100 )
#   Set           Feats        AUC BK           pval
#   <chr>         <list>     <dbl> <list>      <dbl>
# 1 Allen, et al. <chr [61]> 0.816 <dbl [100]>  0.05
```

The function returns the set of gene names that are found in the dataset (61 out of the 68 we specified), the area under the ROC curve (AUC) estimated through LPOCV, the 100 performance values associated with background sets, and the resulting empricial p-value.

By default, the function runs logistic regression. The user is able to run other methods by specifying the appropriate `method=` parameter values: `"svm"` for support vector machines, `"xgb"` for xgboost (random forest), and `"nn"` for neural networks.

``` r
pairs2 <- preparePairs(task2)
evalGeneSets( list("Allen, et al."=gs), task2, pairs2, 100, method="svm" )
#   Set           Feats        AUC BK           pval
#   <chr>         <list>     <dbl> <list>      <dbl>
# 1 Allen, et al. <chr [61]> 0.669 <dbl [100]>  0.16
```
