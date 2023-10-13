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
    
    ggplot2::theme_set(theme_bw())
    
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
        geom_density(aes(x=rt,y=after_stat(!!str2lang("density")),color=cond),data=plot2tmp) +
         # geom_density(aes(x=value,y=after_stat(!!str2lang("density")),color=cond),data=c) +
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
        as_tibble() %>% tidyr::pivot_longer(cols=everything(), names_to = 'variable', values_to = 'value')  %>%  .$value   
    
    the_combine$cuants[the_combine$cond=='Sim'] <-  simul_dat$data %>%  
        split(simul_dat$data$suj) %>%  map(c('rt')) %>%
        map(function(x) {cut(x, breaks = sim_cuants, right=FALSE, na.rm = TRUE)}) %>% 
        as_tibble() %>% tidyr::pivot_longer(cols=everything(), names_to = 'variable', values_to = 'value')  %>%  .$value 
    
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
    
}