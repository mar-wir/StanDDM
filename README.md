# StanDDM

[![forthebadge](https://forthebadge.com/images/badges/built-with-science.svg)](https://forthebadge.com)
[![forthebadge](https://forthebadge.com/images/badges/gluten-free.svg)](https://forthebadge.com)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/Seneketh/StanDDM/issues)

This R package contains a Multi-Level Bayesian fitting procedure for (Perceptual) Decision Making Data and a collection of several Drift Diffusion Models, implemented in the probabilistic programming language Stan. A set of convenience functions for data simulation, model comparison and plotting are also supplied. The aim was to write and test "non-centered" parametrizations of Multi-Level Bayesian models which incorporate inter-trial variabilities and are able to process different amounts of data per subject. The package can be installed to use several utility functions (see below). Also includes random generation of parameters, data simulation and automated recovery for model testing. 

**In its current state, to run the models on a server or cluster it is recommended to download the code directly from this repository and load the functions in directly.**

#### Installation

``` r
devtools::install_github('https://github.com/Seneketh/StanDDM.git', ref = 'master')
```
#### Un-install

``` r
remove.packages("StanDDM")
```
## Usage

### Decision Making Data Simulation and Comparison with Experimental Data
