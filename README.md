# StanDDM

[![forthebadge](https://forthebadge.com/images/badges/built-with-science.svg)](https://forthebadge.com)
[![forthebadge](https://forthebadge.com/images/badges/gluten-free.svg)](https://forthebadge.com)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/Seneketh/StanDDM/issues)

This R package contains a Multi-Level Bayesian fitting procedure for (Perceptual) Decision-Making Data and a collection of several Drift Diffusion Models, implemented in the probabilistic programming language Stan. A set of convenience functions for data simulation, model comparison and plotting are also supplied. The aim was to write and test "non-centered" parametrizations of Multi-Level Bayesian models which incorporate inter-trial variabilities and are able to process different amounts of data per subject. The package can be installed to use several utility functions (see below). Also includes random generation of parameters, data simulation and automated recovery for model testing. 

**In its current state, to run the models on a server or cluster it is recommended to download the code directly from this repository and load the functions in directly.**

#### Installation

``` r
devtools::install_github('https://github.com/Seneketh/StanDDM.git', ref = 'master')
```
#### Un-install

``` r
remove.packages("StanDDM")
```
## Decision-Making Data Simulation and Comparison with Experimental Data

Simulates and structures perceptual decision-making data. These steps create simualted data from a set
of DDM parameter estimates, plots and RMSE/R2 for comparison/interpretation.
Should this package not have been installed, but the R scripts downloaded, then the provided functions have to be loaded in like this:

#### Preamble
```r
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Locate working directory where the script is.
source('./StanDDM/functions.r', local=TRUE) #Source all functions of StanDDM
init_packages() # Load the required packages. Use only if package is not installed.
```
Else, the functions can just be called after loading in the package `library('StanDDM')` if it has been installed.

#### Loading and Structuring Data 

```r
load("./example_name.RData")

patients <- a %>% filter(group == "patient") #separate the data if necessary by group/condition

#structure the data (adapt that function if necessary). See "functions.r" (second function)
data <- experimental_data_processing(patients) 
```
The data should be in long format. Reaction times and diffusion criteria should be indexed by Subject-ID. One line per trial:

| Subject-ID (Str.) | Reaction Times (sec.) | Diffusion Criterion|
| ------------------|-----------------------| --------|
| Subject1          | 1.53                  | 1       |
| Subject1          | 0.8                   | 0       |
| Subject1          | 1.2                   | 1       |
| ...               | ...                   | ...     |

The Diffusion Criterion must either stand for the position (left/right) the stimuli has been presented during that trial (stimuli coding) or if the subject decided correctly or not (accuracy coding).

Now the data has been loaded and formatted, the DDM parameters need also to be read in:

#### Loading Parameters 
```r 
params <- read_csv("./Example_Data_Params/StanDDM_NCEN_Pure_sz.csv")
rownames(params) <- data$forsim$names #giving simulation names of subjects (important!!)
```
With those parameters, data can be simulated now:

#### Simulate Data 
```r
sims <- simulDat(params)
sims <- sims[[1]]

#Rename Vars for Convenience---------------------
experim_dat <- data$forsim$rawdat #Output of 'experimental_data_processing()'
simul_dat <- sims
```
Now the simulated and experimental data are both in the working environment, they can be fed to the comparison functions:

#### Plot and Compare
```r
fit_quality(experim_dat, simul_dat) #calculate RMSE and R2
models_plots(experim_dat, simul_dat) #will create new folder with plots comparing data and simuls
```
This will create a folder with several plots for visual inspection of how well the parameters are capable of reproducing the loaded data. RMSE and R-suquared statistics will be printed to the console.

## Model Fitting (section under construction)

```r
source('functions.r', local=TRUE)

main <- function(num_cores = 4, simulation = FALSE){
    
    seed <- 45100 
    set.seed(seed = seed)

    cat('\n/////////////////////////////////////////////////////////////////////////')
    cat('\n//////////////////// HIERARCHICAL STAN DDM FITTING //////////////////////')
    cat('\n/////////////////////////////////////////////////////////////////////////\n')

    if(!simulation){
        modal <- '_EXPERIMENTAL'
        cat('\nAttention: Experimental data will be fitted ...\n')
    }else{
        modal <- '_SIMULATION'
        cat('\nAttention: Simulating data for parameter recovery and model testing...\n')
    }
    
    warmup <- 20
    iter <- 50

    if(parallel::detectCores() >= num_cores){
        cat('\n4 or more cores are available: Proceeding...\n')
        options(mc.cores = num_cores) #for not occupying to many cpus
        rstan::rstan_options(auto_write = TRUE) #better memory management


        control <- lapply(1:num_cores, function(i) {
            list(
                stepsize = 0.15,
                adapt_delta = 0.85,
                max_treedepth = 30
            )
        })

        cat('\nInitializing packages...\n')
        init_packages()
        

        if(!simulation){
        cat('\nLoading and processing EXPERIMENTAL data...\n')
        load("rt_confs.RData")
        data <- experimental_data_processing(a)
        }else{
            #pass
        }
        
        #source('StanDDM_NCEN_Pure.r', local=TRUE)
        source('StanDDM_NCEN_St.r', local=TRUE)
        #source('StanDDM_NCEN_Sv_Sz_St.r', local=TRUE)
        
    }else{
        cat('\nNot enough cores available: Aborting...\n')
    }
```
