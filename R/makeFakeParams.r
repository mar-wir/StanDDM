#' Pseudo-Random DDM Parameter Generation
#' 
#' Generates a data_frame with pseudo-randomly generated DDM parameters for model
#' validation or similar. 
#' @export
#' @param nsub Number of "simulated" subjects.
#' @param include A vector with the desired inter-trial variabilities.
#' 'sv' for Drift Rate inter-trial variability.
#' 'sz' for Bias inter-trial variability.
#' 'st' for Non-Decision Time inter-trial variability.
#' @return Returns a Data Frame with the desired parameters.
#' @examples 
#' 
#' makeFakeParams(10, include = c('sv', 'st'))
makeFakeParams <- function(nsub=1, include=NULL){
    
    fakeparams <- lapply(1:nsub, function(x) {
        list(
            a     = runif(1, 0.5, 2),
            v    = runif(1, 0.5, 2),
            t     = runif(1, 0.15, 1),
            z    = runif(1, 0.45, 0.55)
        )
    })
    
    if('sv' %in% include){
        fakeparams %<>% map(., function(x) {c(x, sv = runif(1, 0.05, 0.2))})
    }
    
    if('sz' %in% include){
        fakeparams %<>% map(., function(x) {c(x, sz = runif(1, 0.05, 0.2))})
    }
    
    if('st' %in% include){
        fakeparams %<>% map(., function(x) {c(x, st = runif(1, 0.05, 0.2))})}
    
    fakeparams <- as.data.frame(data.table::rbindlist(fakeparams))
    
}