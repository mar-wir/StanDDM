#' Initialize All Required Packages
#' 
#' A function initializing the needed functions. Only required when not installing the package.
init_packages <- function(){
    suppressPackageStartupMessages(if (!require("rstan")) {install.packages("rstan", dependencies = TRUE); library(rstan)})
    suppressPackageStartupMessages(if (!require("StanHeaders")) {install.packages("StanHeaders", dependencies = TRUE); library(StanHeaders)})
    suppressPackageStartupMessages(if (!require("rstantools")) {install.packages("rstantools", dependencies = TRUE); library(rstantools)})
    suppressPackageStartupMessages(if (!require("magrittr")) {install.packages("magrittr", dependencies = TRUE); library(magrittr)})
    suppressPackageStartupMessages(if (!require("data.table")) {install.packages("data.table", dependencies = TRUE); library(data.table)})
    suppressPackageStartupMessages(if (!require("plyr")) {install.packages("plyr", dependencies = TRUE); library(plyr)})
    suppressPackageStartupMessages(if (!require("reshape2")) {install.packages("reshape2", dependencies = TRUE); library(reshape2)})
    suppressWarnings(suppressMessages(if (!require("devtools")) {install.packages("devtools", dependencies = TRUE); library(devtools)}))
    suppressWarnings(suppressMessages(if (!require("tidyverse")) {install.packages("tidyverse", dependencies = TRUE); library(tidyverse)}))
    ggplot2::theme_set(theme_bw())
    }

#' Experimental Data Processing
#' 
#' Read in and process experimental data for model fitting.
#' The column structure of the data has to be following:
#' SUBJECT_ID, REACTION TIME, STIMULI SIDE (alternatively ANSWER SIDE), CORRECTNESS OF RESPONSE
#' Mind that the data needs to be in a long and tidy format.
#' As of now, this function has to be manually modified to change to stimuli fitting.
#' @export
#' @param experimental_data Data frame or tibble. Mind the required columns and names.
#' @return Returns a list with two entries: "forstan" is the data in the correct format to be used in fitting. "forsim" is for data simulation.
experimental_data_processing <- function(a){
#ADAPT THIS FUNCTION SO IT FITS WITH YOUR DATA    
 
    a$suj <- as.factor(a$suj)
    SUB <- length(levels(a$suj))#amount of subjects
    names <- levels(a$suj) #<- 1:SUB
 
    
    aggreg <- ddply(a, c('suj'), summarise, freq = length(rt), .drop = FALSE )
    dat_amt <- aggreg$freq #vector with amount of data per subject
    
    N_min <- min(aggreg$freq) # min amount of data
    N_max <- max(aggreg$freq) # max amount of data
    
    
-----#ADAPT so that the result is a data frame with following columns:
    
    #SUBJECT_ID, REACTION TIME, STIMULI SIDE (alternatively ANSWER SIDE), CORRECTNESS OF RESPONSE
    
    a <- a[c('suj', 'rt', 'click', 'cor')]  
    # a$rt <- a$rt/1000 # only required if RT not in milliseconds
    a$rt[a$rt<0.2] <- 0.2 #lowest amount of RT that works
    # a$click <- a$click-1 #now its 0s and 1s

-----#ADAPT: DEFINE THE CODING OF THE MODEL:
        #STIMULI CODING:  value = sum(STIMULI_SIDE==1 OR 0)
        #ACCURACY CODING:  value = sum(CORECTNESS==1 OR 0)
     
    Nu <- ddply(a, c('suj'), summarise, value = sum(cor==1), .drop = FALSE )$value 
    Nl <- ddply(a, c('suj'), summarise, value = sum(cor==0), .drop = FALSE )$value 

    minRT <- ddply(a, c('suj'), summarise, value = min(rt), .drop = FALSE )$value #min RT per subject
    
    # Reaction times for upper and lower boundary responses, PADDED matrices
    RTu <- array(-1, c(SUB, max(Nu)))
    RTl <- array(-1, c(SUB, max(Nl)))
    
    # Store each subjects' reaction time data
    
-----#ADAPT: DEFINE THE CODING OF THE MODEL:
        #STIMULI CODING:  tmp$rt[tmp$SIMULI_SIDE==1 or 0]
        #ACCURACY CODING:  tmp$rt[tmp$CORRECTNESS==1 or 0]
    for (i in 1:SUB){
        curSubj         <- names[i]
        tmp             <- subset(a, a$suj == curSubj)
        RTu[i, 1:Nu[i]] <- tmp$rt[tmp$cor==1] # (Nu/Nl[i]+1):Nu/Nl_max will be padded with 0's
        RTl[i, 1:Nl[i]] <- tmp$rt[tmp$cor==0] # 0 padding is skipped in likelihood calculation
    }
    
    RTbound <- 0.1
    
    # List of data sent to Stan
    forstan <- list(
        N       = SUB, # Number of subjects
        Nu_max  = max(Nu),  # Max (across subjects) number of upper boundary responses
        Nl_max  = max(Nl),  # Max (across subjects) number of lower boundary responses
        Nu      = Nu,       # Number of upper boundary responses for each subj
        Nl      = Nl,       # Number of lower boundary responses for each subj
        RTu     = RTu,      # upper boundary response times
        RTl     = RTl,      # lower boundary response times
        minRT   = minRT,    # minimum RT for each subject of the observed data
        RTbound = RTbound   # lower bound or RT across all subjects (e.g., 0.1 second)
    )
    
    a$cor <- as.factor(a$cor)
    forsim <- list(
        names = names,
        rawdat = a[, c('suj', 'rt', 'cor')]
    )

        return(list(forstan = forstan, forsim = forsim))
}

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

#' Data Processing for Simulations
#' 
#' Processes output from data simulation functions for posterior use in model validation procedures. Only package internal use.
#' @param simulated_data Data frame or tibble. 
#' @return Returns a list with two entries: "forstan" is the data in the correct format to be used in fitting. "forsim" is for data simulation. 
fake_data_processing <- function(data){
    
    a <- data
    
    a$suj <- as.factor(a$suj)
    SUB <- length(levels(a$suj))#amount of subjects
    names <- levels(a$suj) <- 1:SUB
    
    aggreg <- ddply(a, c('suj'), summarise, freq = length(rt), .drop = FALSE )
    dat_amt <- aggreg$freq #vector with amount of data per subject
    
    N_min <- min(aggreg$freq) # min amount of data
    N_max <- max(aggreg$freq) # max amount of data

    Nu <- ddply(a, c('suj'), summarise, value = sum(cor==1), .drop = FALSE )$value
    Nl <- ddply(a, c('suj'), summarise, value = sum(cor==0), .drop = FALSE )$value 

    
    minRT <- ddply(a, c('suj'), summarise, value = min(rt), .drop = FALSE ) #min RT per subject
    minRT$suj <- as.numeric(minRT$suj)
    minRT <- minRT$value
    
    # Reaction times for upper and lower boundary responses, PADDED matrices
    RTu <- array(-1, c(SUB, max(Nu)))
    RTl <- array(-1, c(SUB, max(Nl)))
    
    # Store each subjects' reaction time data
    for (i in 1:SUB){
        tmp             <- subset(a, a$suj == i)
        RTu[i, 1:Nu[i]] <- tmp$rt[tmp$cor == 1] # (Nu/Nl[i]+1):Nu/Nl_max will be padded with 0's
        RTl[i, 1:Nl[i]] <- tmp$rt[tmp$cor == 0] # 0 padding is skipped in likelihood calculation
    }
    
    RTbound <- 0.1
    
    # List of data sent to Stan
    forstan <- list(
        N       = SUB, # Number of subjects
        Nu_max  = max(Nu),  # Max (across subjects) number of upper boundary responses
        Nl_max  = max(Nl),  # Max (across subjects) number of lower boundary responses
        Nu      = Nu,       # Number of upper boundary responses for each subj
        Nl      = Nl,       # Number of lower boundary responses for each subj
        RTu     = RTu,      # upper boundary response times
        RTl     = RTl,      # lower boundary response times
        minRT   = minRT,    # minimum RT for each subject of the observed data
        RTbound = RTbound   # lower bound or RT across all subjects (e.g., 0.1 second)
    )
    
    
    forsim <- list(
        names = 1:SUB,
        rawdat = a[c('suj', 'rt', 'cor')]
    )
    
    return(list(forstan = forstan, forsim = forsim))
}

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
            z    = runif(1, 0.3, 0.7)
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


#' Data Generation From DDM Parameters
#' 
#' A wrapper around the function 'drift_simuls' which simulates as many data sets as 
#' data frames are provided in a list. This is a convenience function and 
#' the function 'drift_simuls' can also be used directly on just one set on parameters.
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
simulDat <- function(...){
    
    x <- list(...)
    
    simul <- function(params){
        sims <- split(params, seq(nrow(params))) %>% 
            map(drift_simuls)
        
        names(sims) <- rownames(params) 

        rts <- sims %>% map('rt') %>% data.table::melt()
        
        cor <- sims %>% map('response') %>% data.table::melt()
        
        sims <- data.frame('suj'=as.factor(rts$L1), 'rt'=rts$value, 'cor'= as.factor(cor$value), 'cond'='Sim')
        list(params=params, data=sims) 
    }
    
    map(x, simul)
    
}

#' Save R Objects With External Name
#' 
#' A wrapper function around \code{\link[base]{save}}. Assign a string as filename for a R object, 
#' instead the environment name. Will save file to disk.
#' @export
#' @param file Any R object in the environment of any type or class.
#' @param string A string used as filename for the saved R object
saveit <- function(..., string, file) {

    x <- list(...)
    names(x) <- string
    save(list=names(x), file=file, envir=list2env(x))
}

#' Root Mean Squared Error (RMSE)
#' 
#' Calculates the RMSE between simulated and experimental data, stored in a data frame.
#' This is a helper function for \code{\link{fit_quality}}.
#' @param data_frame A data frame that contains simulated and experimental data, one in 
#' each column.
#' @return Returns the RMSE value as double.
rmse <- function(df){ 
    a <- df$Data
    b <- df$Sim
    sqrt(sum((a-b)**2)/length(a))
} 

#' Data Fit Assessment
#' 
#' Calculates \eqn{R^2} values for correct and incorrect mean reaction times and RMSE values
#' for the proportions of correct answers per reaction times cuantiles. 
#' @export
#' @param experim_dat Experimental data in the form produced by \code{\link{experimental_data_processing}}.
#' @param experim_dat Simulated data in the form produced by \code{\link{simulDat}}.
#' @param model_name Optional name for the model to distinguish to which data/model the
#' function was applied to. Default is 'NAME_UNDEFINED'.
#' 
#' @return Prints a table with RMSE and \eqn{R^2} values.
fit_quality <- function(experim_dat, simul_dat, model_name='NAME_UNDEFINED'){
    
    experim_dat$cond <- "Data"
    
    dat <- rbind(experim_dat, simul_dat$data)
    
    ##########
    nbins <- 6

    dat_cuants <- quantile(experim_dat$rt, probs = seq(0, 1, length.out = nbins+1)) 
    dat_cuants[length(dat_cuants)] <- tail(dat_cuants, n=1)+0.1

    sim_cuants <- quantile(simul_dat$data$rt, probs = seq(0, 1, length.out = nbins+1)) 
    sim_cuants[length(sim_cuants)] <- tail(sim_cuants, n=1)+0.1

    dat$cuants[dat$cond=='Data'] <-  experim_dat %>%
        select(-cond) %>%
        split(experim_dat$suj) %>%  map(c('rt')) %>%
        map(function(x) {cut(x, breaks = dat_cuants, right=FALSE, na.rm = TRUE)}) %>%
        melt() %>%  .$value

    dat$cuants[dat$cond=='Sim'] <-  simul_dat$data %>%
        split(simul_dat$data$suj) %>%  map(c('rt')) %>%
        map(function(x) {cut(x, breaks = sim_cuants, right=FALSE, na.rm = TRUE)}) %>%
        melt()  %>%  .$value

    dat$cuants <- unlist(dat$cuants) %>% as.factor()
    dat$cor <- as.numeric(levels(dat$cor))[dat$cor]
    #########
    
   a <- dat %>% group_by(cond, cuants) %>%
        dplyr::summarise(mean = mean(cor)) %>% 
        spread(cond,mean) %>%
        select(-cuants) %>%
        rmse()
   
  
   b <- dat %>% 
        group_by(suj, cond) %>% 
        dplyr::summarize(mean_rt = mean(rt)) %>%
        spread(cond, mean_rt) %>%
        ungroup() %>%
        dplyr::summarise(rmse = rmse(.)) %>%
        pull()
    
   c <- dat %>% 
        group_by(suj, cond) %>% 
        dplyr::summarize(mean_rt = mean(rt)) %>%
        spread(cond, mean_rt) %>%
        ungroup() %>%
        dplyr::summarise(r2 = summary(lm(data=., Data~Sim))$adj.r.squared) %>%
        pull() 
    
   d <- dat %>% group_by(cor, cond, suj) %>%
        dplyr::summarise(mean = mean(rt)) %>% 
        ungroup() %>%
        spread(cond,mean) %>%
        nest(-cor) %>%
        mutate(data = map(data, ~ rmse(as.data.frame(.x)))) %>%
        mutate(data = unlist(data)) %>%
        pull(data)
    
   e <- dat %>% group_by(cor, cond, suj) %>%
        dplyr::summarise(mean = mean(rt)) %>% 
        ungroup() %>%
        spread(cond,mean) %>%
        nest(-cor) %>%
        mutate(data = map(data, ~ summary(lm(data=as.data.frame(.x), Data~Sim))$adj.r.squared)) %>%
        mutate(data = unlist(data)) %>%
        pull(data) 
    
    stats <- list("RMSE (Correct answers per RT cuantiles)" = a,
         "RMSE (Overall Means of RTs aggr. by Subj.)" = b,
         "R2 (Overall Means of RTs aggr. by Subj.)" = c,
         "RMSE (For Incorrect Mean RTs aggr. by Subj.)" = d[1],
         "RMSE (For Correct Mean RTs aggr. by Subj.)" = d[2], 
         "R2 (For Incorrect Mean RTs aggr. by Subj.)" = e[1],
         "R2 (For Correct Mean RTs aggr. by Subj.)" = e[2]
         )

    print(data.frame(stats=unlist(stats)))
    
}

#' Data Simulation Plots
#' 
#' Plots simulated and experimental data. The created plots are: Reaction times densities for 
#' correct and incorrect responses and the proportion of correct answeres for each one of 
#' six reaction time bins. A folder with the model name is created and the plots saved 
#' automatically in the current working directory.
#' @export
#' @param experim_dat Experimental data in the form produced by \code{\link{experimental_data_processing}}.
#' @param experim_dat Simulated data in the form produced by \code{\link{simulDat}}.
#' @param model_name Optional name for the model to distinguish to which data/model the
#' function was applied to. Default is 'NAME_UNDEFINED'.
#' @return Generates and saves plots in a directory with name of \code{model_name}.
models_plots <- function(experim_dat, simul_dat, model_name='NAME_UNDEFINED'){
    
    if(is.null(model_name) | nchar(model_name) < 4 | model_name=='NAME_UNDEFINED'){
        model_name <- "NAME_UNDEFINED"
        dir.create(paste(model_name), showWarnings = FALSE)
    }
    
    
    the_combine <- rbind(data.frame(experim_dat, 'cond'='Data'), simul_dat$data)

    #------------------------------------------------------------
    plot1 <-  ggplot() +
        geom_density(aes(x=rt, group=cor, fill=cor),data=simul_dat$data, alpha=.3) +
        theme(legend.position="bottom") +
        labs(title = paste(model_name),
             subtitle = "Comparison of Simulated Correct vs. Incorrect RTs",
             x = "Reaction Times",
             y = 'Density',
             fill='Correct') +
        scale_colour_hue(name="Data Types:")
    ggsave(plot1, filename = paste(model_name, '/',model_name, '_simRTcorvsncorr.pdf',sep = ''), 
           device = 'pdf',
           scale = 1, width = 12, height = 8, units = "in")
    #------------------------------------------------------------
    
    
    #------------------------------------------------------------
    plot2tmp <- the_combine
    plot2tmp$rt[plot2tmp$cor==0] <- plot2tmp$rt[plot2tmp$cor==0] *-1
    plot2 <-  ggplot() +
        geom_density(aes(x=rt,y=..density..,color=cond),data=plot2tmp) +
        # geom_density(aes(x=value,y=..density..,color=cond),data=c) +
        facet_wrap(~plot2tmp$suj) +
        labs(title = paste(model_name),
             subtitle = "Comparison of Data and Simulated RTs",
             x = "Reaction Times",
             y = 'Density',
             caption = 'Note: Negative/Positive RTs belong to incorrect/correct trials respectively') +
        theme(legend.position="bottom") +
        scale_colour_hue(name="Data Types:")
    ggsave(plot2, filename = paste(model_name, '/', model_name, '_simRTcomparison.pdf',sep = ''), device = 'pdf',
           scale = 1, width = 12, height = 8, units = c("in"))
    #------------------------------------------------------------
    
    #PLOT3: Check in what range the error rates are Sim vs. Data
    #------------------------------------------------------------
    nbins <- 6
    
    dat_cuants <- quantile(experim_dat$rt, probs = seq(0, 1, length.out = nbins+1)) ; dat_cuants[length(dat_cuants)] <- tail(dat_cuants, n=1)+0.1
    
    sim_cuants <- quantile(simul_dat$data$rt, probs = seq(0, 1, length.out = nbins+1)) ; sim_cuants[length(sim_cuants)] <- tail(sim_cuants, n=1)+0.1
    
    the_combine$cuants[the_combine$cond=='Data'] <-  experim_dat %>%  
        split(experim_dat$suj) %>%  map(c('rt')) %>%
        map(function(x) {cut(x, breaks = dat_cuants, right=FALSE, na.rm = TRUE)}) %>% 
        melt() %>%  .$value   
    
    the_combine$cuants[the_combine$cond=='Sim'] <-  simul_dat$data %>%  
        split(simul_dat$data$suj) %>%  map(c('rt')) %>%
        map(function(x) {cut(x, breaks = sim_cuants, right=FALSE, na.rm = TRUE)}) %>% 
        melt()  %>%  .$value 
    
    the_combine$cuants <- unlist(the_combine$cuants) %>% as.factor()
    the_combine$cor <- as.numeric(levels(the_combine$cor))[the_combine$cor]
    
    smrzd <- ddply(the_combine, c("cond","cuants"), summarise,
                         N    = length(cor),
                         mean = mean(cor),
                         sd   = sd(cor),
                         se   = sd / sqrt(N))
    
    pd <- position_dodge(0.1)  #The errorbars overlapped, so use position_dodge to move them horizontally
    plot3 <- ggplot(smrzd, aes(x=cuants, y=mean, colour=cond, group=cond)) +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1,position=pd) +
        geom_line(position=pd) +
        geom_point(position=pd) +
        labs(title = paste(model_name),
             subtitle = "Comparison of Response correctness per RT Quantiles",
             x = "RT Quantile",
             y = 'Proportion Correct Answers',
             caption = 'Note: Bars represent standard errors') +
        theme(legend.position="bottom") +
        scale_colour_hue(name="Data Types:") + ylim(0, 1) +
        theme(legend.position="bottom")

    ggsave(plot3, filename = paste(model_name, '/', model_name, '_cuan_comp.pdf',sep = ''), 
           device = 'pdf',
           scale = 1, width = 12, height = 8, units = "in")
    #------------------------------------------------------------
    

    #------------------------------------------------------------
    smrzd <- ddply(the_combine, c("suj","cond","cuants"), summarise,
                         N    = length(cor),
                         mean = mean(cor),
                         sd   = sd(cor),
                         se   = sd / sqrt(N))
    pd <- position_dodge(0.1)  #The errorbars overlapped, so use position_dodge to move them horizontally
    plot4 <- ggplot(smrzd, aes(x=cuants, y=mean, colour=cond, group=cond)) +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1,position=pd) +
        geom_line(position=pd) +
        geom_point(position=pd) +
        labs(title = paste(model_name),
             subtitle = "Comparison of Response correctness per RT Quantiles and Subject",
             x = "RT Quantile",
             y = 'Proportion Correct Answers',
             caption = 'Note: Bars represent standard errors') +
        theme(legend.position="bottom") +
        scale_colour_hue(name="Data Types:") + ylim(0, 1) +
        theme(legend.position="bottom") + facet_wrap(~smrzd$suj)
    ggsave(plot4, filename = paste(model_name, '/', model_name, '_cuan_comp_suj.pdf',sep = ''),
           device = 'pdf',
           scale = 1, width = 12, height = 8, units = "in")
    
    
    #POSTERIOR PLOT AND OTHERS
    # plot <- rstan::stan_dens(fit)

    #------------------------------------------------------------

}

#' Stanfit MCMC Sampling Diagnostics
#' 
#' A wrapper function around several model diagnostics provided by \code{rstan}. It prints
#' the amount of divergent transitions, treedepth limit hits, energy, 
#' 
#' @export
#' @param experim_dat Experimental data in the form produced by \code{\link{experimental_data_processing}}.
#' @param experim_dat Simulated data in the form produced by \code{\link{simulDat}}.
#' @param model_name Optional name for the model to distinguish to which data/model the
#' function was applied to. Default is 'NAME_UNDEFINED'.
#' 
#' @return Prints a table with RMSE and \eqn{R^2} values.
model_diagnostics <- function(fit){
    
    cat('\nDivergencies:\n')
    check_divergences(fit)
    cat('\nTreedepth:\n')
    rstan::check_treedepth(fit)
    cat('\nE-BFMI:\n')
    check_energy(fit)
    cat('\nEffectice sample size:\n') 
    check_n_eff(fit)
    cat('\nR-hat:\n')  
    check_rhat(fit)
    # plot <- rstan::stan_rhat(fit)
}

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
#MODEL VALIDATION LOO WAIC AND RAW POINTWISE


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


runif2 <- function(n, loc, scale){ # wrap runif in location scale variant
    out <- loc + scale * runif(n)
    as.vector(out)
}

drift_simuls <- function(params, samples = 1000, dt = 1e-04, intra_sv = 1){
    #Translated function from HDDM
    params <- as.list(params)
    nn <- 1000
    a <- params[['a']]
    v <- params[['v']]

    if('st' %in% names(params)){
        start_delay <- runif2(n = samples, 
                              loc = params[['t']],
                              scale = params[['st']]) -params[['st']]/2
    }else{
        start_delay <- rep(1, samples)*params[['t']]
    }
    
    if('sz' %in% names(params)){
        starting_points <- runif2(n = samples,
                                  loc = params[['z']], 
                                  scale = params[['sz']]) -params[['sz']]/2*a
    } else {
        starting_points <- rep(1, samples)*params[['z']]*a
    }
    
    rts <- c()
    step_size <- sqrt(dt)*intra_sv
    drifts <- c()
    
    for (i_sample in 1:samples){
        
       drift <- c()         
       crossed <- FALSE
       iter <- 0
       y_0 <- starting_points[i_sample]
       #drifting:
       
       if('sv' %in% names(params)){
           if(params[['sv']] != 0){
               drift_rate <- rnorm(n = 1, 
                                   mean = v, 
                                   sd = params[['sv']])
           } else {
               drift_rate <- v
           }
       }else {
           drift_rate <- v
       } 
       
       prob_up = 0.5*(1+sqrt(dt)/intra_sv*drift_rate)
       
       while(!crossed){
           # Generate nn steps
           iter <- iter + 1
        
           steps <- runif(nn, min = 0, max = 1)
           
           position <- ((steps < prob_up)*2 - 1) * step_size
           
           position[1] <- position[1] + y_0
           
           position <- cumsum(position)
           
           # Find boundary crossings
           cross_idx <- position[(position < 0) | (position > a)][1]
           cross_idx <- which(cross_idx == position)
           drift <- c(drift, position)
           if(length(cross_idx) > 0){
               crossed <- TRUE
           } else {
               # If not crossed, set last position as starting point
               # for next nn steps to continue drift
               y_0 <- tail(position, n = 1)
           }
           
       }#end while loop for drifts
       
       #find the boundary interception
       y2 = position[cross_idx[1]]
       
       if(cross_idx[1] != 1){
           y1 <- position[cross_idx[1]-1]
       } else {
           y1 <- y_0
       }
       m <- (y2 - y1)/dt #slope
       b <- y2 - m*((iter-1)*nn+cross_idx[1])*dt # intercept 
       if(y2 < 0){
           rt <- ((0 - b) / m)
       } else {
           rt <- ((a - b) / m)
       }

       rts[i_sample] <- (rt + start_delay[i_sample])*sign(y2)
       # delay <- start_delay[i_sample]/dt
       # drifts <- c(drifts, c(rep(1,as.integer(delay))*starting_points[i_sample], 
       # drift[1:as.integer(abs(rt)/dt)]))
       data <- data.frame('rt' = rts)
       data$response <- 1
       data$response[data$rt<0] <- 0
       data$rt <- abs(data$rt)
       
    }#end i_sample
    data
}#end function
