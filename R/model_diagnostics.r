#' Stanfit MCMC Sampling Diagnostics
#' 
#' A wrapper function around several model diagnostics provided by \code{rstan}. It prints
#' the amount of divergent transitions, treedepth limit hits, energy, effective transitions 
#' and R-hat. 
#' 
#' @export
#' @param fit Any stanfit object.
#' @return Prints model diagnostics to console.
model_diagnostics <- function(fit){
    
    cat('\nDivergencies:\n')
    check_divergences(fit)
    cat('\nTreedepth:\n')
    rstan::check_treedepth(fit)
    cat('\nE-BFMI:\n')
    check_energy(fit)
}