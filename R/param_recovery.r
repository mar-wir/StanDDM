#' Parameter Recovery for Model Testing without Experimental Data
#' 
#' This function extracts DDM parameters from a Stanfit object run on simulated data. Then,
#' it and compares those parameters fitted on simulated data with the parameters used to simulate
#' the data in the first place. This function is intended to check if the model is capable to
#' "recover" parameters which are known. Plots for visual inspection are made automatically.
#' Automatic fit diagnostics yet to be implementd with this function. The comparison is done
#' for each parameter over all iterations to detect convergence problems or to spot difficult parameter
#' spaces.
#' 
#' @export
#' @param fit Stanfit object fitted on synthetic data.
#' @param n_subjects Number of "Simulated Subjects".
#' @param auth_params "Ground-truth" parameters used to simulate the data. In the form of the
#' output from the function \code{\link{makeFakeParams}}.
#' @param model_name Optional name for the model to distinguish the plots and created
#' folder easier. Default is 'NAME_UNDEFINED'.
#' @return Creates a folder and saves plots which compare recovered parameters and "ground-truth" parameters
#' for each iteration of the MCMC sampler.
param_recovery <- function(fit, n_subjects, auth_params, model_name='NAME_UNDEFINED'){
    # plot divergencies by parameter:
    mat <- rstan::extract(fit, permuted=FALSE)
    
    n_iter <- dim(mat)[1]
    n_chains <- dim(mat)[2]
    
    param <- rstan::extract(fit, permuted=TRUE, inc_warmup=FALSE) %>% map(unname)
    
    param <- param[as.vector(map(param, length)==(n_iter*n_subjects*n_chains))]
    #filter elements in list with the amount of data only subject-level params have
    
    param %<>% map(function(x){as.data.frame(x)}) %>%
        map(function(x){ apply(x, 2, function(y) cumsum(y)/seq_along(y))}) %>%
        map(unname) %>%
        map(function(x){data.table::melt(x)}) %>%
        map(function(x){x[,c('Var2', 'value')]}) %>%
        map(function(x){x$Var2 <- as.factor(x$Var2); return(x)}) %>%
        melt()
    
    param %<>% select('Var2', 'value', 'L1')
    
    param$L1 <- as.factor(param$L1)
    
    names(param) <- c('suj', 'value', 'param')
    
    ##auth param
    tmpa <- c('a', 'v', 't', 'z', 'st', 'sz', 'sv')
    tmpb <- c('alpha_mu', 'delta_mu', 'tau_mu', 'beta_mu', 'tau_sigma', 'beta_sigma', 'delta_sigma')
    for(.. in 1:length(tmpa)){
        if(tmpa[..] %in% names(auth_params)){
            names(auth_params)[names(auth_params) == tmpa[..]] <- tmpb[..]  
        }    
    }
    
    auth_params$suj <- as.factor(seq.int(nrow(auth_params)))
    
    auth_params %<>% data.table::melt(id.vars='suj')
    names(auth_params) <- c('suj', 'param', 'value')
    param %<>% left_join(., auth_params , by=c('suj', 'param') )
    
    dd <- param %>% filter(param!='log_lik' & param!="lp__") %>% group_by(param) %>%  nest() %>%
        mutate(plot=map2(data, as.list(.$param), ~ggplot(data = .x ) +
                             geom_line(aes(y=value.x, x=1:length(.x$value.x))) +
                             geom_line(aes(y=value.y, x=1:length(.x$value.x)), linetype = "dashed", color='red') +
                             theme(legend.position="bottom") +
                             scale_colour_hue(name="Data Types:")  +
                             theme(legend.position="bottom") + facet_wrap(~suj, scales='free') +
                             labs(title = paste('Parameter:', .y, sep=' '),
                                  subtitle = "Parameter recovery per simulated subject",
                                  x = "Iterations",
                                  y = 'Parameter Values',
                                  caption = 'Note: Warmup iterations not included. The dashed line represents the
                                  parmeter value which generated the data.')))
    
    walk2(.x = dd$plot, .y = dd$param, 
          ~ ggsave(device = 'pdf',
                   scale = 1, width = 12, height = 8, units = c("in"),
                   filename = paste0(model_name, '/', model_name, '_param_recovery_', .y, ".pdf"), 
                   plot = .x))
    
}