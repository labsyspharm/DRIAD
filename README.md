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

