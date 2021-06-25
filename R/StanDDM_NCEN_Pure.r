StanDDM_NCEN_Pure <- function(simulation, modal, data, control, seed, warmup, num_cores, iter){
    
    model_name <- paste('StanDDM_NCEN_Pure', modal, sep="") #modal assigned in main
    dir.create(paste(model_name, sep=''), showWarnings = FALSE)
    cat('\n///////////////Computing model: ', model_name, '///////////////\n')
    
    if(simulation){
        cat('\nCreating random parameters and simulating Data...\n')
        fakeParams <- makeFakeParams(nsub = 5, include = c())    
        fakedat <- simulDat(fakeParams, n = 30)
        fakedat <- fakedat[[1]]
        data <- fake_data_processing(fakedat$data)
        
    }
    
    model_definition <-
    "
        data {
            int<lower=1> N;      // Number of subjects
            int<lower=0> Nu_max; // Max (across subjects) number of upper boundary responses
            int<lower=0> Nl_max; // Max (across subjects) number of lower boundary responses
            int<lower=0> Nu[N];  // Number of upper boundary responses for each subj
            int<lower=0> Nl[N];  // Number of lower boundary responses for each subj
            real RTu[N, Nu_max];  // upper boundary response times
            real RTl[N, Nl_max];  // lower boundary response times
            real minRT[N];       // minimum RT for each subject of the observed data
            real RTbound;        // lower bound or RT across all subjects (e.g., 0.1 second)
        }
    
    parameters {
        vector[4] mu_p;
        vector<lower=0>[4] sigma;
        
        // Subject-level raw parameters (for Matt trick)
        vector[N] alpha_pr;
        vector[N] beta_pr;
        vector[N] delta_pr;
        vector[N] tau_pr;
    }
    
    transformed parameters {
        // Transform subject-level raw parameters
        vector<lower=0>[N]         alpha; // boundary separation
        vector<lower=0, upper=1>[N] beta;  // initial bias
        vector<lower=0>[N]         delta; // drift rate
        vector<lower=RTbound, upper=max(minRT)>[N] tau; // nondecision time
        
        for (i in 1:N) {
            beta[i] = Phi_approx(mu_p[2] + sigma[2] * beta_pr[i]); //Phi approx so bounded between 0 and 1
            tau[i]  = Phi_approx(mu_p[4] + sigma[4] * tau_pr[i]) * (minRT[i]-RTbound) + RTbound; //non decision time //needs to be necessarily smaller than the RT. The bound makes sure that the sampled min non dec time stays //below the smallest RT by the order of the padding, which gets added after the sampling is done(tis sis why //there is a sum at the end. Phi approx = inverse probit.)
        }
        alpha = exp(mu_p[1] + sigma[1] * alpha_pr); //reparametrization as in Gelman manual second ed pg 313 and Kruschkes manual pg. 281
        delta = exp(mu_p[3] + sigma[3] * delta_pr); //exponential link so only positive values, hard limit on 0
    }
    
    model {
        // Hyperparameters
        mu_p  ~ normal(0, 1);
        sigma ~ cauchy(0, 5);
        
        // Individual parameters for non-centered parameterization
        alpha_pr ~ normal(0, 1);
        beta_pr  ~ normal(0, 1);
        delta_pr ~ normal(0, 1);
        tau_pr   ~ normal(0, 1);
        
        // Begin subject loop
        for (i in 1:N) {
            // Response time distributed along wiener first passage time distribution
            RTu[i, :Nu[i]] ~ wiener(alpha[i], tau[i], beta[i], delta[i]);
            RTl[i, :Nl[i]] ~ wiener(alpha[i], tau[i], 1-beta[i], -delta[i]);
            
        } // end of subject loop
    }
    
    generated quantities {
    
        // For log likelihood calculation
        real log_lik[N];
    
        
        { // local section, this saves time and space
            // Begin subject loop
            for (i in 1:N) {
                log_lik[i] = wiener_lpdf(RTu[i, :Nu[i]] | alpha[i], tau[i], beta[i], delta[i]);
                log_lik[i] = log_lik[i] + wiener_lpdf(RTl[i, :Nl[i]] | alpha[i], tau[i], 1-beta[i], -delta[i]);
            }
        }
    }
    "
    
    inits_fixed <- c(0.5, 0.5, 0.5, 0.15)
    init <- lapply(1:num_cores, function(i) {
        list(
            mu_p     = c(log(inits_fixed[1]), qnorm(inits_fixed[2]), log(inits_fixed[3]), 
                         qnorm(inits_fixed[4])),
            sigma    = c(1.0, 1.0, 1.0, 1.0),
            alpha_pr = rep(log(inits_fixed[1]), data$forstan$N),
            beta_pr  = rep(qnorm(inits_fixed[2]), data$forstan$N),
            delta_pr = rep(log(inits_fixed[3]), data$forstan$N),
            tau_pr   = rep(qnorm(inits_fixed[4]), data$forstan$N)
        )
    })
    
    
    pars <-  c("alpha", "beta", "delta", "tau", "log_lik")
            
            
    cat('\nCompiling and fitting the model...\n')
    #-----------------------------------------------------------------        
        fit <- stan(model_code = model_definition,
                    seed = seed,
                    data = data$forstan,
                    pars = pars,
                    warmup = warmup,
                    cores = num_cores,
                    control = control,
                    init = init,
                    iter = iter,
                    chains = num_cores)
            
    
    cat('\nSaving Stanfit object...\n')
    #-----------------------------------------------------------------        
        saveit(fit = fit, 
               string = model_name, 
               file=paste(model_name, '/', model_name, '.RData', sep = ''))
            
            
    cat('\nExtracting and saving parameters...\n')
    #----------------------------------------------------------------- 
        extracted_params <- parameter_extraction(
            fit, 
            numsub=data$forstan$N, 
            names=data$forsim$names)
    
        names(extracted_params) <- c('a', 'z', 'v', 't') #dependent on model
        write.csv(extracted_params, file = paste(model_name, '/',model_name,".csv", sep = ''), 
                  row.names=FALSE)
        
        if(simulation){ #needs model_name, else that would be in main
            write.csv(fakeParams, file = paste(model_name, '/',model_name,"_template_params.csv", 
                                               sep = ''), row.names=FALSE)
            
        }
             
    cat(paste('\nSimulating', data$forstan$N, 'subjects with fitted parameters:\n'))
    #-----------------------------------------------------------------
        sims <- simulDat(extracted_params)
        sims <- sims[[1]]

    cat(paste('\nSaving Plots...\n'))
    #-----------------------------------------------------------------
        models_plots(experim_dat = data$forsim$rawdat, 
                     simul_dat = sims, 
                     model_name = model_name)
                     
        fit_quality(experim_dat = data$forsim$rawdat, 
                 simul_dat = sims, 
                 model_name = model_name)   
                     
        if(simulation){
            param_recovery(fit, data$forstan$N, fakeParams, model_name = model_name)
        }  
    
    cat(paste('\nComputing/plotting model diagnostics:\n'))
    #-----------------------------------------------------------------
        model_diagnostics(fit)
    
    
    cat(paste('\nFor Reproducibility:\n'))
    #-----------------------------------------------------------------
        print(devtools::session_info())
        cat(paste('Seed:', seed, sep = ' '))
    
    
    cat('\n///////////////', model_name, 'has been computed.///////////////\n')
}#end function