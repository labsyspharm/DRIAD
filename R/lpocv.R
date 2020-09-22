#' Generate a binary prediction task
#'
#' Loads a previously-wrangled dataset and composes a binary prediction task
#' 
#' @param fn Filename of a dataset wrangled by wrangleROSMAP() or wrangleMSBB()
#' @param task One of {"AB", "AC", "BC"}, specifying the binary task of interest
#' @importFrom magrittr %>%
#' @export
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

#' Pairs for leave-pair-out cross-validation
#'
#' Generates a set of age-matched pairs for leave-pair-out cross-validation
#' 
#' @param XY dataset, as loaded by prepareTask()
#' @param rs random seed for reproducibility
#' @importFrom magrittr %>%
#' @export
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
        tidyr::unnest( "Cand" ) %>% dplyr::group_by( ID ) %>%
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
    X[,c("ID","Label",gs)]
}

## Leave pair out cross-validation for a given dataset
## XY - dataset must contain columns ID (denoting sample identifiers) and Label
## lPairs - list of vectors-of-length-2 specifying IDs to withhold
## method - one of {"lgr", "svm", "xgb", "nn"}
lpocv <- function( XY, lPairs, method="lgr" )
{
    if( !(method %in% c("lgr", "svm", "xgb", "nn")) )
        stop( "method must be one of lgr, svm, xgb, nn" )
    
    ## Computes AUC from LPOCV
    ## hi - prediction values associated with high-label example in the pair
    ## lo - prediction values associated with low-label example in the pair
    AUC_LPOCV <- function( hi, lo )
    { (0.5*sum(hi==lo) + sum(hi>lo)) / length(hi) }
    
    ## Ensure that only pairs of samples are withheld
    stopifnot( all( range(purrr::map_int( lPairs, length )) == c(2,2) ) )
    stopifnot( all( purrr::flatten_chr(lPairs) %in% XY$ID ) )

    ## Split the XY frame into features and labels
    X <- XY %>% as.data.frame() %>% tibble::column_to_rownames("ID") %>%
        dplyr::select( -Label ) %>% as.matrix
    y <- with( XY, rlang::set_names(Label, ID) )
    
    ## Traverse the pairs and perform cross-validation
    if( method == "lgr" ) {
        RR <- purrr::map( lPairs, ~liblinear_lgr(X, y, .x) )
    } else if( method == "svm" ) {
        RR <- purrr::map( lPairs, ~liblinear_svm(X, y, .x) )
    } else if( method == "xgb" ) {
        RR <- purrr::map( lPairs, ~xgboost(X, y, .x) )
    } else {
        RR <- purrr::map( lPairs, ~nnet(X, y, .x) )
    }

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

#' Evaluate a single gene set
#' 
#' Evaluates a gene set of interest in the context of a given dataset
#' 
#' @param gsi gene set of interest, provided as a character vector 
#' @param XY dataset, as loaded by prepareTask()
#' @param lP list of pairs for cross-validation, as generated by preparePairs()
#' @param nBK number of background sets that should be generated alongside the gene set of interest
#' @param method one of {"lgr", "svm, xgb", "nn"} for logistics regression,
#'    support vector machines, xgboost and neural nets, respectively
#' @return A data frame of results, where the first row corresponds to the gene set of interest,
#'    and all subsequent rows correspond to the generated background sets
#' @importFrom magrittr %>%
#' @export
evalGeneSet <- function( gsi, XY, lP, nBK=0, method="lgr" )
{
    lgs <- list( GSI=gsi )
    
    ## Generate background sets, if requested
    if( nBK > 0 )
        lgs <- c(lgs, genBK( gsi, XY, nBK ))

    ## Downsample the data according to the requested gene sets
    SS <- tibble::enframe( lgs, "Name", "Feats" ) %>%
        dplyr::mutate( Data = purrr::map(Feats, reduceToGenes, XY) )

    ## Run LPOCV on each slice of data
    RR <- SS %>% dplyr::mutate( AUC = purrr::map_dbl(Data, lpocv, lP, method) )

    RR %>% dplyr::select( -Data )
}

## Collapses the result of evalGeneSet into a single-row data frame
summarizeBK <- function( .df )
{
    vBK <- .df %>% dplyr::filter( Name == "BK" ) %>% dplyr::pull(AUC)
    .df %>% dplyr::filter( Name == "GSI" ) %>%
        dplyr::mutate( BK = list(vBK), Name=NULL )
}

#' Evaluate multiple gene sets
#' 
#' Evaluates multiple gene sets after matching them up against a given dataset
#' 
#' @param lGSI list of character vectors, each encapsulating a gene set of interest (GSI)
#' @param XY dataset, as loaded by prepareTask()
#' @param lP list of pairs for cross-validation, as generated by preparePairs()
#' @param nBK number of background sets to generate for each GSI
#' @param method one of {"lgr", "svm, xgb", "nn"} for logistics regression,
#'    support vector machines, xgboost and neural nets, respectively
#' @param rs random seed to allow for reproducibility
#' @return A data frame of results, with one row per GSI
#' @importFrom magrittr %>%
#' @export
evalGeneSets <- function( lGSI, XY, lP, nBK=0, method="lgr", rs=100 )
{
    set.seed(rs)
    
    ## Isolate genes that exist in the data
    lGSI <- purrr::map( lGSI, intersect, colnames(XY))

    ## Fill in empty names
    vn <- stringr::str_c("GeneSet",1:length(lGSI))
    if( is.null(names(lGSI)) ) names(lGSI) <- vn
    names(lGSI) <- purrr::map2( names(lGSI), vn, ~`if`(.x=="",.y,.x) )

    ## Evaluate individual gene sets and combine results into a single data frame
    fo <- furrr::future_options( seed=TRUE )
    R <- furrr::future_map( lGSI, evalGeneSet, XY, lP, nBK, method, .options=fo ) %>%
        purrr::map( summarizeBK ) %>% dplyr::bind_rows( .id = "Set" )

    ## Compute empirical p-values
    R %>% dplyr::mutate( pval = purrr::map2_dbl(AUC, BK,
                                                ~`if`(length(.y)==0, NA, mean(.x <= .y))) )
}

#' Load gene sets from a .gmt file
#' 
#' Parses a .gmt file and puts it into the list format
#' 
#' @param iName index of the token containing gene set names
#'    (This is typically 1 for the Broad's MSigDB sets, and 2 for PathwayCommons sets)
#' @return A named list of character vectors, one per gene set
#' @importFrom magrittr %>%
#' @export
read_gmt <- function( fn, iName=1 )
{
    readr::read_lines(fn) %>% stringr::str_split( "\\t" ) %>%
        rlang::set_names( purrr::map_chr(., dplyr::nth, iName) ) %>%
        purrr::map( ~.x[-2:-1] )
}

