## A variant of rosmap.R that reduces the space of molecular
##   features to protein-coding (pc) genes
##
## by Artem Sokolov

library( tidyverse )

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
    Xraw <- read_tsv( fnX, col_types = cols() )

    ## Load biotype annotations and retrieve names of protein-coding genes
    cat( "Downloading biotype annotations...\n" )
    fnBT <- synapser::synGet( "syn14236139", downloadLocation = destDir )$path
    BT <- read_csv( fnBT, col_types=cols() ) %>% filter( gene_biotype=="protein_coding" ) %>%
        select( ENSEMBL = gene_id, HUGO = gene_name )

    ## Map ENSEMBL Gene IDs to HUGO
    ## Reduce feature space to protein-coding genes
    ## There is a single duplicate: ENSG00000254093.3 and ENSG00000258724.1 map to the
    ##   same HUGO ID. However, ENSG00000258724.1 is almost all 0s, so we drop it.
    cat( "Mapping gene IDs to HUGO...\n" )
    X <- Xraw %>% filter( gene_id != "ENSG00000258724.1" ) %>%
        mutate( ENSEMBL = str_split( gene_id, "\\.", simplify=TRUE )[,1] ) %>%
        inner_join( BT, by="ENSEMBL" ) %>% select( -tracking_id, -gene_id, -ENSEMBL )

    ## Log-transform the data and combine the replicates
    cat( "Additional processing...\n" )
    fmed <- function(x) {x %>% as.matrix %>% apply( 1, median )}
    XX <- X %>% mutate_at( vars(-HUGO), ~log2(.x+1) ) %>%
        mutate( `492_120515_j` = fmed(select( ., contains("492_120515") )) ) %>%
        select( -`492_120515_0`, -`492_120515_6`, -`492_120515_7` ) %>%
        gather( rnaseq_id, Value, -HUGO ) %>%
        mutate( rnaseq_id = str_sub( rnaseq_id, 0, -3 ) )

    ## Match sample IDs against individual identifiers
    cat( "Matching sample and individual IDs...\n" )
    fnZ <- synapser::synGet( "syn3382527", downloadLocation = destDir )$path
    XZ <- read_csv(fnZ, col_types=cols()) %>%
        select( projid, rnaseq_id ) %>% na.omit %>%
        distinct %>% inner_join( XX, ., by="rnaseq_id" )

    ## Match expression data up against the following clinical covariates:
    ## ID, PMI, AOD, CDR, Braak, BrodmannArea
    cat( "Matching against clinical covariates...\n" )
    fnY <- synapser::synGet( "syn3191087", version=6, downloadLocation = destDir )$path
    Y <- suppressWarnings( read_csv(fnY, col_types=cols()) )%>%
        select( projid, PMI = pmi, AOD = age_death, CDR = cogdx, Braak = braaksc ) %>%
        mutate( BrodmannArea = "BM9,BM46" )

    ## Combining everything into a common data frame
    cat( "Finalizing...\n" )
    XY <- inner_join( Y, XZ, by="projid" ) %>% rename( ID = projid, Barcode = rnaseq_id ) %>%
        spread( HUGO, Value )

    ## Write out wrangled dataset to file
    fnOut <- file.path( destDir, "rosmap.tsv.gz" )
    cat( "Writing output to", fnOut, "\n" )
    write_tsv( XY, fnOut )
}
