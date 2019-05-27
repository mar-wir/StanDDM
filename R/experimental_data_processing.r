#' Experimental Data Processing
#' 
#' Read in and process experimental data for model fitting. It is assumed that there is at least 1 
#' response per side (left/right)! Else, this function will not work!
#' The column structure of the data has to be following:
#' SUBJECT_ID, REACTION TIME, STIMULI SIDE (alternatively ANSWER SIDE), CORRECTNESS OF RESPONSE
#' Mind that the data needs to be in a long and tidy format. The data processing is inspired by Ahn et al.'s 
#' procedure.  
#' @export
#' @references hBayesDM Package and paper: \url{https://www.ncbi.nlm.nih.gov/pubmed/29601060}
#' @param a Data frame or tibble. Mind the required columns and names. 
#' See the example data set with \code{StanDDM::example_data}
#' @return Returns a list with two entries: "forstan" is the data in the correct format to be used in fitting. 
#' "forsim" is for data simulation.
experimental_data_processing <- function(a){
    
     if(!('suj' %in% names(a) & 'rt' %in% names(a) & 'cor' %in% names(a) & 
         'crit' %in% names(a))){
        stop('\n\nThe provided data has not the correct column names or lacks one or more columns. 
Please run "StanDDM:example_data" to see how the correct format looks like.\n')
     }
    
    onsiders <- a %>% 
        group_by(suj) %>% summarise(mean = 
                                        mean(as.numeric(as.character(crit)))) %>% 
        filter(mean == 0 | mean == 1) %>% 
        rownames() %>% length()
    
    if(onsiders > 0){
        stop('\n\nOne or more subjects in the provided data only features one type of responses
(only 0s or 1s). If you simulated this data with "simulDat", try to increase the amount of trials
(increase the argument "n") or alternatively, re-run the simulations.\n')
    }
    
    a$suj <- as.factor(a$suj)
    SUB <- length(levels(a$suj))#amount of subjects
    names <- levels(a$suj) #<- 1:SUB
    
    aggreg <- plyr::ddply(a, c('suj'), summarise, freq = length(rt), .drop = FALSE )
    dat_amt <- aggreg$freq #vector with amount of data per subject
    
    N_min <- min(aggreg$freq) # min amount of data
    N_max <- max(aggreg$freq) # max amount of data
    
    
    #A data frame with following columns:
    
    #SUBJECT_ID, REACTION TIME, STIMULI SIDE (alternatively ANSWER SIDE), CORRECTNESS OF RESPONSE
    
    a <- a[c('suj', 'rt', 'crit', 'cor')]  
    a$rt[a$rt<0.2] <- 0.2 #lowest amount of RT that works
    
    Nu <- plyr::ddply(a, c('suj'), summarise, value = sum(crit==1), .drop = FALSE )$value 
    Nl <- plyr::ddply(a, c('suj'), summarise, value = sum(crit==0), .drop = FALSE )$value 
    
    minRT <- plyr::ddply(a, c('suj'), summarise, value = min(rt), .drop = FALSE )$value #min RT per subject
    
    # Reaction times for upper and lower boundary responses, PADDED matrices
    RTu <- array(-1, c(SUB, max(Nu)))
    RTl <- array(-1, c(SUB, max(Nl)))
    
    # Store each subjects' reaction time data
    
    #ADAPT: DEFINE THE CODING OF THE MODEL:
    #STIMULI CODING:  tmp$rt[tmp$SIMULI_SIDE==1 or 0]
    #ACCURACY CODING:  tmp$rt[tmp$CORRECTNESS==1 or 0]
    for (i in 1:SUB){
        curSubj         <- names[i]
        tmp             <- subset(a, a$suj == curSubj)
        RTu[i, 1:Nu[i]] <- tmp$rt[tmp$crit==1] # (Nu/Nl[i]+1):Nu/Nl_max will be padded with 0's
        RTl[i, 1:Nl[i]] <- tmp$rt[tmp$crit==0] # 0 padding is skipped in likelihood calculation
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