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