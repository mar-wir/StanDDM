#' Data Generation From DDM Parameters
#' 
#' A wrapper around the function \code{\link{drift_simuls}} which simulates as many data sets as 
#' data frames with parameters are provided in a list. This is a convenience function and 
#' the function \code{\link{drift_simuls}} can also be used directly on just one set on parameters.
#' @export
#' @param parameter_data_frames A list of data frames that contain DDM parameters in
#' the format returned by the 'makeFakeParams' function: 1 column per parameter 
#' and as many rows as subjects. Can also be applied to just one data frame!
#' The rows can be named (in a data frame, in tibbles it is deprecated).
#' @return Returns a list of Data Frames with simulated data. Or alternatively, a list of 1
#' data frame.
#' @examples 
#' 
#' parameters_1 #parameters extracted from a fitting procedure
#' parameters_2 #simulated parameters
#' 
#' sims <- simulDat(list(parameters_1, parameters_2))
#' 
#' sim_1 <- sims[[1]]
#' sim_2 <- sims[[2]]
simulDat <- function(..., n = 1000){
    
    x <- list(...)
    
    simul <- function(params, .samples = n){
        sims <- split(params, seq(nrow(params))) %>% 
            map(., ~drift_simuls(.x, samples = .samples))
        
        names(sims) <- rownames(params) 
        
        rts <- sims %>% map('rt') %>% data.table::melt()
        
        cor <- sims %>% map('response') %>% data.table::melt()
        
        sims <- data.frame('suj'=as.factor(rts$L1), 'rt'=rts$value, 'cor'= as.factor(cor$value), 'cond'='Sim')
        
        onsiders <- sims %>% 
            group_by(suj) %>% summarise(mean = 
                                            mean(as.numeric(as.character(cor)))) %>% 
            filter(mean == 0 | mean == 1) %>% 
            rownames() %>% length()
        
        if(onsiders > 0){
            warning('\n\nOne or more subjects in the provided data only features one type of responses
(only 0s or 1s). This data set will not be able to be processed by "experimental_data_processing"
and cannot be fit.\n')
        }
        
        list(params=params, data=sims) 
    }
    
    map(x, simul)
    
}