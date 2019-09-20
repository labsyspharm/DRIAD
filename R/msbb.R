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
BT <- readr::read_csv( syn("syn14236139"), col_types=readr::cols() ) %>%
    dplyr::filter( gene_biotype=="protein_coding" ) %>%
    dplyr::select( ENSEMBL = gene_id, HUGO = gene_name )

## Match sample barcodes against individual IDs
cat( "Downloading ID map...\n" )
Z <- syn("syn6100548") %>% readr::read_csv(col_types=readr::cols()) %>%
    dplyr::select( Barcode=barcode, individualIdentifier ) %>%
    dplyr::distinct() %>% dplyr::mutate_at( "Barcode", as.character )

## Load clinical covariates and combine with the expression matrix
cat( "Downloading clinical covariates...\n" )
Y <- syn("syn6101474") %>% readr::read_csv(col_types=readr::cols()) %>%
    dplyr::select( individualIdentifier, PMI, AOD, CDR, Braak=bbscore ) %>%
    dplyr::inner_join( Z, by = "individualIdentifier" ) %>%
    dplyr::rename( ID=individualIdentifier )

## Read raw expression matrix
## PINX1 is the only gene with duplicate entries, but one of the entries has
##   a higher total count, so we keep it and discard the other entry.
cat( "Downloading and wrangling expression data...\n" )
synX <- c( BM10="syn16795931", BM22="syn16795934", BM36="syn16795937", BM44="syn16795940" )
X <- purrr::map( synX, function(s) {
    read.delim( syn(s), check.names=FALSE ) %>%
        tibble::rownames_to_column( "ENSEMBL" ) %>%
        dplyr::inner_join( BT, by="ENSEMBL" ) %>%
        dplyr::filter( ENSEMBL != "ENSG00000258724" ) %>%
        dplyr::select( -ENSEMBL ) %>%
        tidyr::gather( RNA_ID, Value, -HUGO ) %>%
        dplyr::mutate( Barcode = stringr::str_split(RNA_ID, "_", simplify=TRUE)[,3],
                      RNA_ID = NULL ) %>%
        dplyr::inner_join( Y, by="Barcode" ) %>%
        tidyr::spread( HUGO, Value ) })

## Compose filenames
## Save each region into a separate file
cat( "Writing output to", destDir, "\n" )
fnOut <- names(X) %>% stringr::str_sub( 3 ) %>%
    stringr::str_c( "msbb", ., ".tsv.gz" ) %>% purrr::map( ~file.path(destDir,.x) )
R <- purrr::map2( X, fnOut, readr::write_tsv )
