#' Precheck for Underflow Values Imputation
#' 
#' Checks Stanfit objects for underflow values before applying the \code{\link{model_judgement}} function.
#' 
#' @param arg A Stanfit object fitted on synthetic data.
#' @param current_model Name of the current model (Currently not in use.)
#' @param lik_name Name under which the log likelihoods have been saved in the models. Needs to be identical 
#' @return A list with information regarding the underflow values for each Stanfit object.
precheck <- function(arg, lik_name) {
    
    nms <- names(arg)
    
    for(.. in 1:length(arg)){
        
        lik <- as.matrix(arg[[..]], pars = lik_name)
        colnames(lik) <- NULL
        current_model <- nms[..]
        check <- sum(apply(log(exp(lik)), 2, FUN = sum))
        
        if(!is.finite(check)==TRUE | is.nan(check)==TRUE | is.na(check)==TRUE){
            status <- "Problematic!"
            
        }
        else{ status <- 'OK'
        # cat('\n',current_model, 'is', status)
        }
        
        arg[[..]] <- list("Status" = status)
        
    }
    
    df <- rbindlist(arg, fill=TRUE)
    df <- cbind(Model = names(arg), df)
    gc()
    
    if("Problematic!" %in% df$Status){
        
        return(list("status" = 0, 
                    "troublemakers" = subset.data.frame(as.data.frame(df), Status=="Problematic!") ))
    }
    
    else{
        return(list("status" = 1))
    }
}