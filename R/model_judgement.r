#' Model Comparison via LOO-CV, WAIC and Raw LPPD
#' 
#' Applies log-likelihood based model comparison to any number of stanfit objects and extracts LOO-CV, 
#' WAIC and Raw LPPD measures and plots the models on a absolute scale ranging from 0 (best models) to 
#' 1 (worst models). WAIC and LOO-CV code has been adapted directly from the "LOO" package. Because
#' most methods to prevent log-exp calculation-underflow failed, this function imputes underflow values
#' with the mean of sampled values. When this happens, the function will report the amount of imputed values
#' to inform the user. More than 5% of imputed values should not be accepted.
#' Mind how Stan saves log-likelihood values! Consult the Stan manual to check how they are saved correctly!
#' Because all models need to be loaded into memory, be wary if loading to big stanfits. Lighten the Stanfit
#' objects by discarding unnecesary iterations is advised. 
#' 
#' @export
#' @param stanfits At least 2 stanfit objects.
#' @param lik_name Name under which the log likelihoods have been saved in the models. Needs to be identical 
#' across all Stanfit objects.
#' @param impute_inf A boolean which regulates if underflow values should be automatically imputed or not. If
#' \code{FALSE}, the models with such values will just be ignored. If \code{TRUE}, a report will be generated
#' on how many values were imputed and for which models.
#' @examples 
#' 
#' load(fit1.Rdata)
#' load(fit2.Rdata)
#' 
#' model_judgement(fit1, fit2, impute_inf = TRUE)
#' @return Prints a table and generates a plot with the model ranked according to LOO-CV, WAIC and Raw LPPD.
model_judgement <- function(..., lik_name = "log_lik", impute_inf = TRUE) {
    #compare all models in environment: compare <- ls(pattern = 'StanDDM')
    if (!require("data.table")) {install.packages("data.table", dependencies = TRUE); library(data.table)}
    if (!require("rstan")) {install.packages("rstan", dependencies = TRUE); library(rstan)}
    
    arg <- list(...)
    nms <- substitute(list(...))[-1]
    names(arg) <- as.list(nms)
    
    if(impute_inf){
        imps <- arg} 
    else{
        pre <- precheck(arg, lik_name)
        if(pre[[1]]==1){
            cat('\nAll models are OK to go.\n')
        }
        else{
            cat('\nWarning: The following models feature Log Exp underflow NaNs:\n\n')
            print(pre$troublemakers)
            cat('___________________________________\n\n')
        }
    }
    
    
    
    
    for(.. in 1:length(arg)){
        
        current_model <- nms[..]
        
        judgement <- (waic(arg[[..]], current_model, lik_name, impute_inf))
        
        ic <- judgement['ic']
        arg[[..]] <- list("lppd" = as.numeric(ic$ic$total['lpd']),
                          "PSIS-LOO" = as.numeric(ic$ic$total['elpd_loo']),
                          "WAIC" = as.numeric(ic$ic$total['waic']))
        
        if(impute_inf){
            imp <- judgement['imps']
            imps[[..]] <- list('Model Name' = as.character(current_model), 
                               "%Likelihood imp." = imp$imps$lik_imp,
                               "% LOO CV imp." = imp$imps$loo_imp) 
        }
        
        
    }
    if(impute_inf){
        cat('\n\n___________________________________')
        cat('\n Warning: The following models had -Inf values imputed (shows % of data imputed):\n\n')
        print(rbindlist(imps, fill=TRUE))
        cat('___________________________________\n\n')
    }
    
    #plotting:
    gc()
    df <- rbindlist(arg, fill=TRUE)
    
    madmax <- apply(abs(df), 2, FUN = max, na.rm=TRUE)
    
    df <- sweep(abs(df), 2, madmax, "/")
    df <- cbind(Model = names(arg), df)
    molten_df <- melt(df)
    # pd <- position_dodge(0.05)  #The errorbars overlapped, so use position_dodge to move them horizontally
    plot <- ggplot(molten_df, aes(x=variable, y=value, colour=Model,group=Model)) +
        # geom_errorbar(aes(ymin=value-se, ymax=mean+se), width=.1,position=pd) +
        geom_line() +
        geom_point()+
        geom_text(aes(label=Model),hjust=0, vjust=0, 
                  position = position_dodge(width=0.3),  size=2.5) +
        ylab("Absolute ranking of model fit") +
        xlab("IC")+
        scale_colour_hue(name="Models:") +
        # ylim(0, 1) +
        theme(legend.position="bottom")
    
    gc()
    return(list(df, plot))
    
    
}