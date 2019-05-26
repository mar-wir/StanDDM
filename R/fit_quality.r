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
        data.table::melt() %>%  .$value
    
    dat$cuants[dat$cond=='Sim'] <-  simul_dat$data %>%
        split(simul_dat$data$suj) %>%  map(c('rt')) %>%
        map(function(x) {cut(x, breaks = sim_cuants, right=FALSE, na.rm = TRUE)}) %>%
        data.table::melt()  %>%  .$value
    
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