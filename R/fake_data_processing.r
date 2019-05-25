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