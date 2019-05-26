#' Drift Diffusion Data Simulation
#' 
#' Simulates reaction times and binary decisions (0/1) from a set of DDM parameters. This function is not likelihood
#' based, but executes all steps of the diffusion process. Can be used directly, but intended use is within the 
#' wrapper function \code{\link{simulDat}}. This function has been adapted from the python package "HDDM". 
#' 
#' @export
#' @param params A list of parameters which should at least contain 'a', 'v', 'z' and 't'. Inter trial variabilities 
#' can also be included: 'sv', 'sz' and/or 'st', in any combination.
#' @param samples Amount of decisions to be simulated. Can be conceptualized as 'trials'. 
#' @param dt Function steps or 'resolution'. Should not be altered, but can.
#' @param intra_sv Intra-trial variability. When simulating data with parameters fitted via Stan, then the value of
#' \code{1} should be used. The implemented likelihood in Stan is calibrated to that value. When using parameters 
#' fitted with a Ratcliff procedure, then this value should be set to \code{0.1}.
#' @return A data frame with two columns, 'rt' and 'response', as long as 'samples'.
drift_simuls <- function(params, samples = 500, dt = 1e-04, intra_sv = 1){
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
