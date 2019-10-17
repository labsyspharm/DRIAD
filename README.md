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
