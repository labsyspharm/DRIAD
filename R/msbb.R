## Wrangling of MSBB RNAseq data and matching clinical covariates
##   Reduces gene space to the protein coding regions
##
## by Artem Sokolov

library( magrittr )

destDir <- "~/data/amp-ad/msbb"
syn <- function(s) {synapser::synGet(s, downloadLocation=destDir)$path}

## Create directory if it doesn't exist
dir.create( destDir, recursive=TRUE, showWarnings=FALSE )
cat( "Wrangling MSBB dataset to", destDir, "\n" )

## Login to Synapse and download/wrangle data
cat( "Logging in to Synapse...\n" )
synapser::synLogin()

## Load biotype annotations and retrieve names of protein-coding genes
cat( "Downloading biotype annotations...\n" )
fnBT <- syn( "syn14236139" )
BT <- readr::read_csv( fnBT, col_types=readr::cols() ) %>%
    dplyr::filter( gene_biotype=="protein_coding" ) %>%
    dplyr::select( ENSEMBL = gene_id, HUGO = gene_name )

## Read raw expression matrix
## PINX1 is the only gene with duplicate entries, but one of the entries has
##   a higher total count, so we keep it and discard the other entry.
cat( "Downloading expression data...\n" )
synX <- c( BM10="syn16796116", BM22="syn16796117", BM36="syn16796121", BM44="syn16796123" )
X <- purrr::map( synX, function(s) {
    read.delim( syn(s), check.names=FALSE ) %>%
        tibble::rownames_to_column( "ENSEMBL" ) %>%
        dplyr::inner_join( BT, by="ENSEMBL" ) %>%
        dplyr::filter( ENSEMBL != "ENSG00000258724" ) %>%
        dplyr::select( -ENSEMBL ) %>%
        tidyr::gather( RNA_ID, Value, -HUGO ) %>%
        dplyr::mutate( barcode = stringr::str_split(RNA_ID, "_", simplify=TRUE)[,3],
                      RNA_ID = NULL ) }) %>% dplyr::bind_rows()

## Match sample barcodes against individual IDs and brain region information
cat( "Annotating samples with brain region...\n" )
XZ <- syn( "syn6100548" ) %>% readr::read_csv(col_types=readr::cols()) %>%
    dplyr::select( BrodmannArea, barcode, individualIdentifier ) %>%
    dplyr::distinct() %>% dplyr::mutate_at( "barcode", as.character ) %>%
    dplyr::inner_join( X, by = "barcode" )

## Load clinical covariates and combine with the expression matrix
cat( "Downloading clinical covariates...\n" )
XY <- syn( "syn6101474" ) %>% readr::read_csv(col_types=readr::cols()) %>%
    dplyr::select( individualIdentifier, PMI, AOD, CDR, bbscore ) %>%
    dplyr::inner_join( XZ, by = "individualIdentifier" )

## Flatten the matrix to samples-by-(clin+genes) and save to file
cat( "Finalizing...\n" )
RR <- tidyr::spread( XY, HUGO, Value ) %>%
    dplyr::rename( ID = individualIdentifier, Braak = bbscore, Barcode = barcode )
fnOut <- file.path( destDir, "msbb.tsv.gz" )
cat( "Writing output to", fnOut, "\n" )
readr::write_tsv( RR, fnOut )
