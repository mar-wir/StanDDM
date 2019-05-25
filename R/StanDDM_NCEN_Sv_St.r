StanDDM_NCEN_Sv_St <- function(simulation, modal, data, control, seed, warmup, num_cores, iter){
    
    model_name <- paste('StanDDM_NCEN_Sv_St', modal, sep="") #modal assigned in main
    dir.create(paste(model_name, '/', sep=''), showWarnings = FALSE)
    cat('\n///////////////Computing model: ', model_name, '///////////////\n')
    
    if(simulation){
        cat('\nCreating random parameters and simulating Data...\n')
        fakeParams <- makeFakeParams(nsub = 10, include = c('sv', 'st'))     
        fakedat <- simulDat(fakeParams)
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
    //Group level parameters in vecotized form, one for each parameter.
      vector[4] mu_mu; //group level mean means
      vector<lower=0>[4] mu_sigma; //group level mean variances
    
      vector[4] sigma_mu; //group level variance means
      vector<lower=0>[4] sigma_sigma; //group level variance variances
    
      // Subject-level raw parameters (for Matt trick)
      vector[N] alpha_mu_pr;
      vector[N] beta_mu_pr;
    
      vector[N] delta_mu_pr;
      vector[N] delta_sigma_pr;
    
      vector[N] tau_mu_pr;
      vector[N] tau_sigma_pr;
    
    	
      // Trial-level raw parameters (for Matt trick)
      vector[(sum(Nu)+sum(Nl))] tau_pr;
      vector[(sum(Nu)+sum(Nl))] delta_pr;
    }
    
    transformed parameters {
      // Transform  raw parameters (primes)
    
      // Subject-level declarations: MEANS
      vector<lower=0>[N]         alpha_mu; // boundary separation mean subject level
      vector<lower=0, upper=1>[N] beta_mu;  // initial bias mean subject level
      vector<lower=0>[N]         delta_mu; // drift rate mean subject level
      vector<lower=RTbound, upper=max(minRT)>[N] tau_mu; // nondecision time mean subject level
    
      // Subject-level declarations: VARIANCES
      vector<lower=0>[N]         tau_sigma; // drift rate variance subject level
      vector<lower=0>[N]         delta_sigma; // drift rate variance subject level
    
      // Trial-level declarations
      vector[(sum(Nu)+sum(Nl))] tau;
      vector[(sum(Nu)+sum(Nl))] delta;
    
    //-------------------------------------
    
      //Subject-level MEANS transformations
      for (i in 1:N) {
        beta_mu[i] = Phi_approx(mu_mu[2] + mu_sigma[2] * beta_mu_pr[i]);
        tau_mu[i]  = Phi_approx(mu_mu[4] + mu_sigma[4] * tau_mu_pr[i]) * (minRT[N]-RTbound) + RTbound;
      }
    
      alpha_mu = exp(mu_mu[1] + mu_sigma[1] * alpha_mu_pr); //reparametrization as in  Gelman manual second ed pg 313 and Kruschkes manual pg. 281
      delta_mu = exp(mu_mu[3] + mu_sigma[3] * delta_mu_pr);
    
      //Subject-level VARIANCES transformations
      tau_sigma = exp(sigma_mu[4] + sigma_sigma[4] * tau_sigma_pr);//already vectorized for subject level in declaration
      delta_sigma = exp(sigma_mu[3] + sigma_sigma[3] * delta_sigma_pr);//already vectorized for subject level in declaration
    
    //Trial-level transformations: MUST BE SLICED AS IN MODEL BLOCK
    
        for (i in 1:N) {//begin subject loop
    
            //loop-scoped definitions for index calcs
            int NuL; int NuU; int NlL; int NlU;
    
            NuL = sum(Nu[1:i]) -Nu[i] +1;
            NuU = sum(Nu[1:i]);
            NlL = sum(Nl[1:i]) -Nl[i] +1 + sum(Nu);
            NlU = sum(Nl[1:i]) + sum(Nu);
    
    		tau[NuL:NuU] = Phi_approx(tau_mu[i] + tau_sigma[i] * tau_pr[NuL:NuU])* (minRT[i]-RTbound) + RTbound;
    		tau[NlL:NlU] = Phi_approx(tau_mu[i] + tau_sigma[i] * tau_pr[NlL:NlU])* (minRT[i]-RTbound) + RTbound;
    
    		delta[NuL:NuU] = exp(delta_mu[i] + delta_sigma[i] * delta_pr[NuL:NuU]);
    		delta[NlL:NlU] = exp(delta_mu[i] + delta_sigma[i] * delta_pr[NlL:NlU]);
    
        }//end subject loop
    
    }
    
    model {
        // Group level parameters (all vectorized!)
        mu_mu  ~ normal(0, 1); //Group mean mean
        mu_sigma ~ cauchy(0, 5); //Group mean variance
    
        sigma_mu ~ cauchy(0, 5); //Group variance mean
        sigma_sigma ~ cauchy(0, 5); //Group variance variance
    
        // Subject-level sampling statements: MEANS
        alpha_mu_pr ~ normal(0, 1);
        beta_mu_pr ~ normal(0, 1);
        delta_mu_pr ~ normal(0, 1);
        tau_mu_pr ~ normal(0, 1);
    
        // Subject-level sampling statements: VARIANCES
        tau_sigma_pr ~ cauchy(0, 5);
        delta_sigma_pr ~ cauchy(0, 5);
    
        // Trial-level sampling statements: PRIMES
        //alpha_pr ~ normal(0, 1);
        //beta_pr  ~ normal(0, 1);
        delta_pr ~ normal(0, 1);
        tau_pr   ~ normal(0, 1);
    
        // Begin subject loop
        for (i in 1:N) {
        
            int NuL; int NuU; int NlL; int NlU;
        
            NuL = sum(Nu[1:i]) -Nu[i] +1;
            NuU = sum(Nu[1:i]);
        
            RTu[i, 1:Nu[i]] ~ wiener(alpha_mu[i], tau[NuL:NuU], beta_mu[i], delta[NuL:NuU]);
        
            NlL = sum(Nl[1:i]) -Nl[i] +1 + sum(Nu);
            NlU = sum(Nl[1:i]) + sum(Nu);
        
            RTl[i, 1:Nl[i]] ~ wiener(alpha_mu[i], tau[NlL:NlU], 1-beta_mu[i], -delta[NlL:NlU]);
    
        } // end of subject loop
    }// end model block
    
    generated quantities {
        
        vector[(sum(Nu)+sum(Nl))] log_lik;
        
        for (i in 1:N) {
            
            int NuL; int NuU; int NlL; int NlU;
            
            NuL = sum(Nu[1:i]) -Nu[i] +1;
            NuU = sum(Nu[1:i]);
            
            for (j in NuL:NuU){
                log_lik[j] = wiener_lpdf(RTu[i, 1:Nu[i]] | alpha_mu[i], tau[j], beta_mu[i], delta[j]);
            }
            
            NlL = sum(Nl[1:i]) -Nl[i] +1 + sum(Nu);
            NlU = sum(Nl[1:i]) + sum(Nu);
            
            for (j in NlL:NlU){
                log_lik[j] = wiener_lpdf(RTl[i, 1:Nl[i]] | alpha_mu[i], tau[j], 1-beta_mu[i], -delta[j]);
            }
        }
    }
    "
    
    inits_fixed <- c(0.5, 0.5, 0.5, 0.15)
    init <- lapply(1:num_cores, function(i) {
        list(
            mu_mu     = c(log(inits_fixed[1]), qnorm(inits_fixed[2]), log(inits_fixed[3]), qnorm(inits_fixed[4])),
            mu_sigma    = c(1.0, 1.0, 1.0, 1.0),
            sigma_mu     = c((inits_fixed[1]), qnorm(inits_fixed[2]), log(inits_fixed[3]), qnorm(inits_fixed[4])),
            sigma_sigma    = c(0.1, 0.1, 0.1, 0.1),
            alpha_mu_pr = rep(log(inits_fixed[1]), data$forstan$N),
            beta_mu_pr  = rep(qnorm(inits_fixed[2]), data$forstan$N),
            tau_mu_pr   = rep(qnorm(inits_fixed[4]), data$forstan$N),
            tau_sigma_pr   = rep(0.1, data$forstan$N),
            delta_mu_pr   = rep(qnorm(inits_fixed[3]), data$forstan$N),
            tau_pr = rep(log(inits_fixed[4]), sum(data$forstan$Nl)+sum(data$forstan$Nu)),
            delta_pr = rep(log(inits_fixed[3]), sum(data$forstan$Nl)+sum(data$forstan$Nu))
            
        )
    })
    
    pars <- c('delta_sigma', 'tau_sigma', "alpha_mu", "beta_mu", "delta_mu", "tau_mu", 'log_lik')
            
            
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
    
        names(extracted_params) <- c('sv', 'st','a', 'z', 'v', 't')  
        write.csv(extracted_params, file = paste(model_name, '/',model_name,".csv", sep = ''), 
                  row.names=FALSE)
        
        if(simulation){ #needs model_name, else that would be in main
            write.csv(fakeParams, file = paste(model_name, '/',model_name,"_template_params.csv", 
                                               sep = ''), row.names=FALSE)
            
        }
        #TODO: name parsing?
            
    cat(paste('\nSimulating', data$forstan$N, 'subjects with fitted parameters:\n'))
    #-----------------------------------------------------------------
        sims <- simulDat(extracted_params)
        sims <- sims[[1]]
 
    cat(paste('\nSaving Plots...\n'))
    #-----------------------------------------------------------------
        models_plots(experim_dat = data$forsim$rawdat, 
                     simul_dat = sims, 
                     model_name = model_name)
    
    cat(paste('\nComputing/plotting model diagnostics:\n'))
    #-----------------------------------------------------------------
        model_diagnostics(fit)
    
    
    cat(paste('\nFor Reproducibility:\n'))
    #-----------------------------------------------------------------
        print(devtools::session_info())
        cat(paste('Seed:', seed, sep = ' '))
    
    
    cat('\n///////////////', model_name, 'has been computed.///////////////\n')
}#end function