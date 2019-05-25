#' LOO-CV, WAIC and Raw LPPD Calculations
#' 
#' Implements LOO-CV, WAIC and Raw LPPD for the \code{\link{model_judgement}} function. Contains
#' the helper function "colVars" which calculates row-wise variances efficiently.
#' 
#' @param stanfit A Stanfit object fitted on synthetic data.
#' @param current_model Name of the current model (Currently not in use.)
#' @param lik_name Name under which the log likelihoods have been saved in the models. Needs to be identical 
#' across all Stanfit objects.
#' @param impute_inf A boolean which regulates if underflow values should be automatically imputed or not. If
#' \code{FALSE}, the models with such values will just be ignored. If \code{TRUE}, a report will be generated
#' on how many values were imputed and for which models.
#' @return A list with LOO-CV, WAIC and Raw LPPD calculations.
waic <- function(stanfit, current_model, lik_name, impute_inf){
    #http://kylehardman.com/BlogPosts/View/6 DIC code also from Gelman
    #Modified code of www.stat.columbia.edu/~gelman/research/unpublished/waic_stan.pdf
    #from gist.github.com/ihrke for underflow probs
    
    colVars <- function(a) {
        n <- dim(a)[[1]]; 
        c <- dim(a)[[2]]; 
        result <- (.colMeans(((a - matrix(.colMeans(a, n, c), 
                                          nrow = n, ncol = c, byrow = TRUE)) ^ 2), n, c) * n / (n - 1))
        return(result)
    }
    
    log_lik <- rstan::extract(stanfit, lik_name)$log_lik
    
    
    if(impute_inf){
        lik_imp <- sum(is.infinite(log_lik))/length(log_lik)
        log_lik[is.infinite(log_lik)] <- mean(log_lik[!is.infinite(log_lik)])
    }
    dim(log_lik) <- if (length(dim(log_lik))==1) c(length(log_lik),1) else
        c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))
    S <- nrow(log_lik)
    n <- ncol(log_lik)
    #log pointwise 
    lpd <- log(colMeans(exp(log_lik))) #only when posterior simulations in Stan are correctly made (and possible)
    
    #waic
    p_waic <- colVars(log_lik)
    elpd_waic <- lpd - p_waic
    waic <- -2*elpd_waic
    
    #loo
    loo_weights_raw <- 1/exp(log_lik-max(log_lik))
    
    if(impute_inf){
        loo_imp <- sum(is.infinite(loo_weights_raw))/length(loo_weights_raw)
        loo_weights_raw[is.infinite(loo_weights_raw)] <- mean(loo_weights_raw[!is.infinite(loo_weights_raw)])
    }
    
    loo_weights_normalized <- loo_weights_raw/matrix(colMeans(loo_weights_raw),nrow=S,ncol=n,byrow=TRUE)
    loo_weights_regularized <- pmin (loo_weights_normalized, sqrt(S))
    
    elpd_loo <- log(colMeans(exp(log_lik)*loo_weights_regularized)/colMeans(loo_weights_regularized))
    p_loo <- lpd - elpd_loo
    elpd_loo <- elpd_loo*-2
    pointwise <- cbind(waic,lpd,p_waic,elpd_waic,p_loo,elpd_loo)
    total <- colSums(pointwise)
    se <- sqrt(n*colVars(pointwise))
    ic <- list(waic=total["waic"], elpd_waic=total["elpd_waic"],
               p_waic=total["p_waic"], elpd_loo=total["elpd_loo"], p_loo=total["p_loo"],
               pointwise=pointwise, total=total, se=se)
    
    if(impute_inf){
        imps <- list('lik_imp' = lik_imp, 'loo_imp' = loo_imp)
        return(list('ic' = ic, 'imps' = imps))
    } else{
        return(list('ic' = ic))
    }
}