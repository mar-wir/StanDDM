#' Automatic DDM Parameter Extraction
#' 
#' Automatically extracts parameter estimates from a stanfit object.
#' @param stanfit_object Output of a previous stanfit.
#' @param num_subjects Amount of fitted subjects. Required for automatic extraction.
#' @param names Subject names (Vector of strings) for correct plotting. Is not required, but recommended.
#' @export
#' @return Returns a CSV file with all parameter estimates.
parameter_extraction <- function(fit, numsub, names=NULL){
    suj_params <- rstan::extract(fit, permuted = TRUE)
    suj_params %<>% map(., function(x) colMeans(as.matrix(x), dims = 1)) %>% map(unname)
    suj_params <- suj_params[as.vector(map(suj_params, length)==numsub)];
    suj_params %<>% as.data.frame %>% select(-contains("_pr"), -contains("_lik")) 
    row.names(suj_params) <- names
    suj_params
}