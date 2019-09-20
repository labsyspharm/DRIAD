## A variant of rosmap.R that reduces the space of molecular
##   features to protein-coding (pc) genes
##
## by Artem Sokolov

library( magrittr )

wrangleROSMAP <- function( destDir = "~/data/amp-ad/rosmap" )
{
    ## Create directory if it doesn't exist
    dir.create( destDir, recursive=TRUE, showWarnings=FALSE )
    cat( "Wrangling ROS/MAP dataset to", destDir, "\n" )

    ## Login to Synapse and download/wrangle data
    cat( "Logging in to Synapse... " )
    synapser::synLogin()

    ## Read raw expression matrix
    cat( "Downloading expression data...\n" )
    fnX <- synapser::synGet( "syn3505720", downloadLocation = destDir )$path
    Xraw <- readr::read_tsv( fnX, col_types = readr::cols() )

    ## Load biotype annotations and retrieve names of protein-coding genes
    cat( "Downloading biotype annotations...\n" )
    fnBT <- synapser::synGet( "syn14236139", downloadLocation = destDir )$path
    BT <- readr::read_csv( fnBT, col_types=readr::cols() ) %>%
        dplyr::filter( gene_biotype=="protein_coding" ) %>%
        dplyr::select( ENSEMBL = gene_id, HUGO = gene_name )

    ## Map ENSEMBL Gene IDs to HUGO
    ## Reduce feature space to protein-coding genes
    ## There is a single duplicate: ENSG00000254093.3 and ENSG00000258724.1 map to the
    ##   same HUGO ID. However, ENSG00000258724.1 is almost all 0s, so we drop it.
    cat( "Mapping gene IDs to HUGO...\n" )
    X <- Xraw %>% dplyr::filter( gene_id != "ENSG00000258724.1" ) %>%
        dplyr::mutate( ENSEMBL = stringr::str_split( gene_id, "\\.", simplify=TRUE )[,1] ) %>%
        dplyr::inner_join( BT, by="ENSEMBL" ) %>%
        dplyr::select( -tracking_id, -gene_id, -ENSEMBL )

    ## Log-transform the data and combine the replicates
    cat( "Additional processing...\n" )
    fmed <- function(x) {x %>% as.matrix %>% apply( 1, median )}
    XX <- X %>% dplyr::mutate_at( dplyr::vars(-HUGO), ~log2(.x+1) ) %>%
        dplyr::mutate( `492_120515_j` = fmed(dplyr::select(., dplyr::contains("492_120515"))) ) %>%
        dplyr::select( -`492_120515_0`, -`492_120515_6`, -`492_120515_7` ) %>%
        tidyr::gather( rnaseq_id, Value, -HUGO ) %>%
        dplyr::mutate( rnaseq_id = stringr::str_sub( rnaseq_id, 0, -3 ) )

    ## Match sample IDs against individual identifiers
    cat( "Matching sample and individual IDs...\n" )
    fnZ <- synapser::synGet( "syn3382527", downloadLocation = destDir )$path
    XZ <- readr::read_csv(fnZ, col_types=readr::cols()) %>%
        dplyr::select( projid, rnaseq_id ) %>% na.omit %>%
        dplyr::distinct() %>% dplyr::inner_join( XX, ., by="rnaseq_id" )

    ## Match expression data up against the following clinical covariates:
    ## ID, PMI, AOD, CDR, Braak
    cat( "Matching against clinical covariates...\n" )
    fnY <- synapser::synGet( "syn3191087", version=6, downloadLocation = destDir )$path
    Y <- suppressWarnings( readr::read_csv(fnY, col_types=readr::cols()) )%>%
        dplyr::select( projid, PMI = pmi, AOD = age_death, CDR = cogdx, Braak = braaksc )

    ## Combining everything into a common data frame
    cat( "Finalizing...\n" )
    XY <- dplyr::inner_join( Y, XZ, by="projid" ) %>%
        dplyr::rename( ID = projid, Barcode = rnaseq_id ) %>%
        tidyr::spread( HUGO, Value )

    ## Write out wrangled dataset to file
    fnOut <- file.path( destDir, "rosmap.tsv.gz" )
    cat( "Writing output to", fnOut, "\n" )
    readr::write_tsv( XY, fnOut )
}
