## Quick-and-dirty plots for one-off gene set runs on ROSMAP
## For large-scale analyses, use the BTR framework: https://github.com/pvtodorov/btr
##
## by Artem Sokolov

library( magrittr )

## Reduces a dataset to the binary task of interest
## fn - filename of a dataset wrangled by wrangleROSMAP() or wrangleMSBB()
## task - one of {"AB", "AC", "BC"}
prepareTask <- function( fn, task )
{
    ct <- readr::cols(ID=readr::col_character(),
                      Braak=readr::col_integer())
    X <- readr::read_tsv( fn, col_types=ct )

    ## Define the task mapping
    ## Note that the mapping has 1-based indexing
    ##   while Braak score ranges from 0 to 6
    taskMap <- list(
        AB = c("neg", "neg", "neg", "pos", "pos", NA, NA),
        AC = c("neg", "neg", "neg", NA, NA, "pos", "pos"),
        BC = c(NA, NA, NA, "neg", "neg", "pos", "pos")
    )
    
    ## Argument verification
    stopifnot( task %in% names(taskMap) )
    stopifnot( is.integer( X$Braak ) )

    ## Reduce rows accordingly
    ## Note the +1 to align 0-based and 1-based indexing
    X %>% dplyr::mutate( Label = taskMap[[task]][Braak+1] ) %>%
        dplyr::filter( !is.na(Label) )
}

## Prepares a set of pairs for leave-pair-out cross-validation
## XY - dataset, as loaded by prepareTask()
## rs - random seed for reproducibility
preparePairs <- function( XY, rs=100 )
{
    ## Argument verification
    stopifnot( all(c("ID","AOD","Label") %in% colnames(XY)) )
    stopifnot( all(sort(unique(XY$Label)) == c("neg","pos")) )

    ## Select the relevant columns, standardize age
    S <- XY %>% dplyr::select( ID, AOD, Label ) %>%
        dplyr::mutate_at( "AOD", dplyr::recode, "90+" = "90" ) %>%
        dplyr::mutate_at( "AOD", as.numeric )

    ## For each sample, identify potential pair candidates from the other class
    ## Compute the distance in age space and select the closest candidate
    ## Break ties randomly
    set.seed(rs)
    S1 <- S %>% dplyr::mutate( Cand = purrr::map(Label, ~dplyr::filter(S, Label != .x)) ) %>%
        dplyr::mutate_at( "Cand", purrr::map, dplyr::rename_all, stringr::str_c, "c" ) %>%
        tidyr::unnest() %>% dplyr::group_by( ID ) %>%
        dplyr::mutate( Dist = abs(AODc-AOD) ) %>%
        dplyr::filter( Dist == min(Dist) ) %>%
        dplyr::slice( sample(1:(dplyr::n()),1) ) %>% dplyr::ungroup()

    ## Finalize the format
    S1 %>% dplyr::select( ID, IDc ) %>% dplyr::mutate( Index = 1:(dplyr::n()) ) %>%
        with(split(., Index)) %>% purrr::map(dplyr::select, -Index) %>%
            purrr::map(unlist) %>% purrr::map(unname)
}

## Reduces a dataset to the gene set of interest
## gs - vector character containing the gene set of interest
## X - Dataset, as loaded by prepareTask()
reduceToGenes <- function( gs, X )
{
    ## Argument verification
    stopifnot( all(gs %in% colnames(X)) )

    ## Reduce columns accordingly
    X %>% dplyr::select( ID, Label, dplyr::one_of(gs) )
}

## Train-test for a single pair using liblinear implementation
liblinear <- function( XY, vTest )
{
    ## Argument verification
    stopifnot( all( sort(unique(XY$Label)) == c("neg","pos") ) )
    
    ## Split the data into train and test
    Xte <- XY %>% dplyr::filter( ID %in% vTest ) %>%
        dplyr::select( -ID, -Label ) %>% as.matrix
    Xtr <- XY %>% dplyr::filter( !(ID %in% vTest) ) %>%
        dplyr::select( -ID, -Label ) %>% as.matrix
    ytr <- dplyr::filter( XY, !(ID %in% vTest) )$Label

    ## Train a model and apply it to test data
    m <- LiblineaR::LiblineaR( Xtr, ytr, type=0 )
    ypred <- predict( m, Xte, proba=TRUE )$probabilities[,"pos"]
    XY %>% dplyr::filter( ID %in% vTest ) %>%
        dplyr::select( ID, Label ) %>% dplyr::mutate( Pred = ypred )
}

## Leave pair out cross-validation for a given dataset
## XY - dataset must contain columns ID (denoting sample identifiers) and Label
## lPairs - list of vectors-of-length-2 specifying IDs to withhold
lpocv <- function( XY, lPairs )
{
    ## Computes AUC from LPOCV
    ## hi - prediction values associated with high-label example in the pair
    ## lo - prediction values associated with low-label example in the pair
    AUC_LPOCV <- function( hi, lo )
    { (0.5*sum(hi==lo) + sum(hi>lo)) / length(hi) }
    
    cat( "." )
    
    ## Ensure that only pairs of samples are withheld
    stopifnot( all( range(purrr::map_int( lPairs, length )) == c(2,2) ) )
    stopifnot( all( purrr::flatten_chr(lPairs) %in% XY$ID ) )

    ## Traverse the pairs and perform cross-validation
    RR <- purrr::map( lPairs, ~liblinear(XY, .x) )

    ## Compute AUC
    dplyr::bind_rows( RR, .id="index" ) %>% dplyr::select( -ID ) %>%
        tidyr::spread( Label, Pred ) %>% with( AUC_LPOCV(pos, neg) )
}

## Given a dataset and a gene set of interest, generates background sets of equal size
## gsi - gene set of interest
## X - Dataset, as loaded by prepareTask()
## nBK - number of background sets to generate
## vExclude - set of identifiers to exclude from sampling
genBK <- function( gsi, X, nBK, vExclude=c("ID", "PMI", "AOD", "CDR",
                                           "Braak", "Barcode", "Label") )
{
    ## Intersect the gene set of interest against dataset's feature space
    vFullSet <- setdiff( colnames(X), vExclude )
    vGSI <- intersect( gsi, vFullSet )
    nGSI <- length(vGSI)

    ## Sample random sets of genes of matching size
    seq(1, length.out=nBK) %>% purrr::map( ~sample(vFullSet, nGSI) ) %>%
        rlang::set_names( rep("BK", nBK) )
}

## Evaluates a gene set of interest in the context of a given dataset
## gsi - gene set of interest as a character vector 
## XY - Dataset, as loaded by prepareTask()
## lP - List of pairs for cross-validation
## nBK - number of background sets that should be generated alongside
evalGeneSet <- function( gsi, XY, lP, nBK=0 )
{
    lgs <- list( GSI=gsi )
    
    ## Generate background sets, if requested
    if( nBK > 0 )
        lgs <- c(lgs, genBK( gsi, XY, nBK ))

    ## Downsample the data according to the requested gene sets
    cat( "Generating data slices...\n" )
    SS <- tibble::enframe( lgs, "Name", "Feats" ) %>%
        dplyr::mutate( Data = purrr::map(Feats, reduceToGenes, XY) )
    
    ## Run LPOCV on each slice of data
    cat( "Running LPOCV" )
    RR <- SS %>% dplyr::mutate( AUC = purrr::map_dbl(Data, lpocv, lP) )
    cat( "\n" )

    RR %>% dplyr::select( -Data )
}

## Collapses the result of evalGeneSet into a single-row data frame
summarizeBK <- function( .df )
{
    vBK <- .df %>% dplyr::filter( Name == "BK" ) %>% dplyr::pull(AUC)
    .df %>% dplyr::filter( Name == "GSI" ) %>%
        dplyr::mutate( BK = list(vBK), Name=NULL )
}

## Evaluates multiple gene sets after matching them up against a given dataset
## lGSI - list of character vectors, each encapsulating a gene set of interest (GSI)
## XY - Dataset, as loaded by prepareTask()
## lP - list of pairs for cross-validation
## nBK - number of background sets to generate for each GSI
## rs - random seed to allow for reproducibility
evalGeneSets <- function( lGSI, XY, lP, nBK=0, rs=100 )
{
    set.seed(rs)
    
    ## Isolate genes that exist in the data
    lGSI <- purrr::map( lGSI, intersect, colnames(XY))

    ## Fill in empty names
    vn <- stringr::str_c("GeneSet",1:length(lGSI))
    if( is.null(names(lGSI)) ) names(lGSI) <- vn
    names(lGSI) <- purrr::map2( names(lGSI), vn, ~`if`(.x=="",.y,.x) )

    ## Evaluate individual gene sets and combine results into a single data frame
    R <- purrr::map( lGSI, evalGeneSet, XY, lP, nBK ) %>%
        purrr::map( summarizeBK ) %>% dplyr::bind_rows( .id = "Set" )

    ## Compute empirical p-values
    R %>% dplyr::mutate( pval = purrr::map2_dbl(AUC, BK,
                                                ~`if`(length(.y)==0, NA, mean(.x <= .y))) )
}

