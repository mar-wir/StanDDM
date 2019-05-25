#' StanDDM: Multi-level Bayesian Fitting Procedure for Decision-Making Data
#' 
#' This is the main function, which calles all models defined in "include" and applies all Stan sampler
#' arguments to the fitting.
#' 
#' @export
#' @param data If simulations are not done (\code{simulation = FALSE}), then a data frame with decision-making
#' data has to be provided. Required columns are 'suj' (subject ID), 'rt' (reaction times in seconds), 'crit' for the 
#' diffusion criteria (accuracy coding or stimuli coding, more info on the github readme) and 'cor' (correct or 
#' incorrect answer). See the github readme for a example.
#' @param include_models A vector of strings to indicate which models should be fitted to the data or used in 
#' the simulations. 'Pure', 'st', 'sv', 'sz', 'sv_sz', 'st_sv', 'sz_st' and 'sv_sz_st' can be included in any order.
#' It is also possible to load in custom models. See the github readme for instructions.
#' @return Returns a list of Data Frames with simulated data. Or alternatively, a list of 1
#' @param num_cores Number of CPU cores to be used for the fitting procedure. Note that for each core, a MCMC
#' chain will be assigned automatically. The minimum amount of available cores on your machine should be 4.
#' @param simulation A boolean that if \code{simulation = TRUE}, a parameter recovery study for the current model is launched.
#' In that mode, a set of parameters will be randomly generated and also the corresponding data for 10 simulated 
#' subjects. That synthetic data will then be fitted with the indicated model. Plots and fitting quality will be
#' saved in your working directory. Interesting are the by-iteration plots for each parameter.
#' If \code{simulation = FALSE}, a data frame has to be provided in the correct format. The models will fit
#' the provided data.
#' @param warmup Warmup iterations which will be discarded. 500 to 1000 are recommended, but should be adjusted
#' to accomodate rig capacities and model properties.
#' @param iter MCMC iterations per chain. 1000 to 5000 are recommended but should be adjusted
#' to accomodate rig capacities and model properties.
#' @param stepsize Initial NUTS step-size. Will adapt after a some iterations. Adjust to own parameter space
#' if necessary. See Stan manual for more information.
#' @param adapt_delta See Stan manual for more information. 
#' @param max_treedepth See Stan manual for more information.
StanDDM <- function(data = NULL,
                    include_models = c(),
                    num_cores = 4, 
                    simulation = FALSE, 
                    warmup = 500, 
                    iter = 2000, 
                    stepsize = 0.15, 
                    adapt_delta = 0.85, 
                    max_treedepth = 30){
    
    seed <- 45100 
    set.seed(seed = seed)
    
    cat('\n/////////////////////////////////////////////////////////////////////////')
    cat('\n//////////////////// HIERARCHICAL STAN DDM FITTING //////////////////////')
    cat('\n/////////////////////////////////////////////////////////////////////////\n')

    if(!simulation){
        if(is.null(data)){
            stop('\n\nNo data found. Please provide some data.\n')
         }
        modal <- '_EXPERIMENTAL'
        cat('\nAttention: Experimental data will be fitted ...\n')
    }else{
        modal <- '_SIMULATION'
        cat('\nAttention: Simulating data for parameter recovery and model testing...\n')
    }
    
    if(parallel::detectCores() >= num_cores){
        cat('\n4 or more cores are available: Proceeding...\n')
        options(mc.cores = num_cores) #for not occupying to many cpus
        rstan::rstan_options(auto_write = TRUE) #better memory management


        control <- lapply(1:num_cores, function(i) {
            list(
                stepsize = stepsize,
                adapt_delta = adapt_delta,
                max_treedepth = max_treedepth
            )
        })

        if(!simulation){
            cat('\nProcessing EXPERIMENTAL data...\n')
        }
        
        args <- list(simulation, modal, data, control, 
                     seed, warmup, num_cores, iter)        

        i <- include_models
        
        if (is.null(i) | any(i == 'Pure') | any(i == 'pure' | any(i == 'none'))) {
            do.call(StanDDM_NCEN_Pure, args)
        } else if (any(i == 'st')) {
            do.call(StanDDM_NCEN_St, args)
        } else if (any(i == 'sv')) {
            do.call(StanDDM_NCEN_Sv, args)
        } else if (any(i == 'sz')) {
            do.call(StanDDM_NCEN_Sz, args)
        } else if (any(i == 'st_sz' | i == 'sz_st')) {
            do.call(StanDDM_NCEN_St_Sz, args)
        } else if (any(i == 'st_sv' | i == 'sv_st')) {
            do.call(StanDDM_NCEN_Sv_St, args)
        } else if (any(i == 'sv_sz' | i == 'sz_sv')) {
            do.call(StanDDM_NCEN_Sv_Sz, args)
        } else if (any(i == 'sv_sz_st')) {
            do.call(StanDDM_NCEN_Sv_Sz_St, args)
        }else {
        tryCatch(do.call(paste0(i[1]), args),
                error = stop('\n\nSomething went wrong with the "include" argument.\n 
                             Check which models you want to fit and how to inlude them.\n'))    
        }        
  
    }else{
        stop('\n\nNot enough cores available: Aborting...\n')
    }

}# end function