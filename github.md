
<!-- README.md is generated from README.Rmd. Please edit that file -->

# epcTools

<!-- badges: start -->

[![status](https://img.shields.io/badge/status:-v0.1c-green)]()
[![status: support](https://img.shields.io/badge/support:-yes-green)]()
<!-- badges: end -->

epcTools (**e**co**p**hysiological **c**onstrained trait **Tools**), is
a R package designed to identify how the environment influences
parameters of Ornstein-Uhlenbeck models of trait evolution. This
approach is currently oriented around ecophysiologically constrained
traits, or traits whose expression is modulated by a certain level of
environmental resource availability, but it can be used for other time
varying processes. The current implementation is for single traits.

## Installation

‘epcTools’ can be installed from this github repository. Simply enter
the following code:

``` r
# install.packages("devtools")
devtools::install_github("mason-linscott/epcTools")
```

## General usage

The general workflow for empirical datasets is as follows: create an EPC
cache object -\> run epc.max.lik search

While it may sound simple, all the data needs to be in the right format.
We have included an example empirical dataset with the package to help
users understand the required input format.

``` r
data(diatom_data)

#diatom_cache<-cache.epc(mtree=diatom_data$diatom_tree,epc_params = "alpha",environment_data = #diatom_data$diatom_environment,n_slice = 10, trait_data = diatom_data$diatom_traits[,1,drop=FALSE])

#diatom_epc_search<-epc.max.lik(diatom_cache)
```

## A simulated example

We will need to simulate a phylogeny and environmental vector. We can
then use these two to simulate an EPC cache object which will vary with
the trait of interest.

The sim.epc function requires several arguments. First, the tree and
environmental values but also including other parameters documented
elsewhere such as the type of EPC relationship and what variable(s) are
affected by the environment. Here we simulate a cache object with alpha
and theta varying with the environment and ten time slices.

``` r
#Tree and environment
library(epcTools)
tree_example<-epcTools::sim.tree(100,200,100,1)
env_example<-(c(6,8,11,13,9,8,5,4,8,10))

#Parameters to simulate under
base_sig2=0.5
base_b0_t=1
base_b1_t=2
base_b0_a=0.01
base_b1_a=0.03
base_slice=10


sim0.5f_200b_at_linear<-sim.epc(tree_example,epc_params=c("alpha","theta"),m_type="linear",X=env_example,n_slice=base_slice,1,b0_a=base_b0_a,b1_a=base_b1_a,sig2=base_sig2,b0_t=base_b0_t,b1_t=base_b1_t)
```

## Visualize how the traits are varying with the environment

``` r
#EPC phenogram
sim.epc.phenogram(sim0.5f_200b_at_linear[[1]],10)
```

<img src="man/figures/README-phenogram-1.png" width="100%" />

## epcTools maximum likelihood search

Now, lets see if we can recover the parameters we simulated our model
under. This requires invoking the epc.max.lik() function. Note that this
function can take a vector of starting parameters fed by the user or
find one automatically using the start.searcher function.

``` r
#EPC phenogram
example_results<-epc.max.lik(sim0.5f_200b_at_linear[[1]])
#> Starting initial iterative maximum likelihood search...
#> First search     AIC:344.62
#> Starting second iterative maximum likelihood search...
#> Second search    AIC:344.62
#> Starting third iterative maximum likelihood search...
#> Third search     AIC:691.33
#> Starting last iterative maximum likelihood search...
#> Last search  AIC:344.62
#> [1] "Calculating intervals at a confidence level of 95%"
#> [1] "Done replicate 500"
#> [1] "CI of values (the 34 replicates within 5.53524884675818 neglnL of the optimum)"
#>        neglnL parameter_1 parameter_2 parameter_3 parameter_4 parameter_5
#> [1,] 167.3115   0.4632926    2.023406   0.4051995  0.01828872  0.02284182
#> [2,] 172.5429   0.6444879    2.068737   0.5214812  0.02714467  0.02934762
#> [1] "Rough volume of good region is 5.50290734920944e-08"
#> [1] "Done replicate 1000"
#> [1] "CI of values (the 114 replicates within 5.53524884675818 neglnL of the optimum)"
#>        neglnL parameter_1 parameter_2 parameter_3 parameter_4 parameter_5
#> [1,] 167.3115   0.4632926    2.023020   0.4051995  0.01828872  0.02284182
#> [2,] 172.8132   0.6746954    2.068737   0.5214812  0.02714467  0.02934762
#> [1] "Rough volume of good region is 6.47497688816286e-08"
#> [1] "Done replicate 1500"
#> [1] "CI of values (the 210 replicates within 5.53524884675818 neglnL of the optimum)"
#>        neglnL parameter_1 parameter_2 parameter_3 parameter_4 parameter_5
#> [1,] 167.3115   0.4632926    2.023020   0.3921361  0.01219756  0.02284182
#> [2,] 172.8141   0.7026640    2.068737   0.5214812  0.02714467  0.02934762
#> [1] "Rough volume of good region is 1.37645025339221e-07"
#> [1] "Done replicate 2000"
#> [1] "CI of values (the 285 replicates within 5.53524884675818 neglnL of the optimum)"
#>        neglnL parameter_1 parameter_2 parameter_3 parameter_4 parameter_5
#> [1,] 167.3115   0.4632926    2.023020   0.3921361  0.01170811  0.02284182
#> [2,] 172.8259   0.7026640    2.068737   0.5214812  0.02714467  0.02934762
#> [1] "Rough volume of good region is 1.42152258412326e-07"
#> [1] "Done replicate 2500"
#> [1] "CI of values (the 372 replicates within 5.53524884675818 neglnL of the optimum)"
#>        neglnL parameter_1 parameter_2 parameter_3 parameter_4 parameter_5
#> [1,] 167.3115   0.4632926    2.023020   0.3921361 0.009373159  0.02284182
#> [2,] 172.8259   0.7026640    2.068737   0.5214812 0.027144672  0.02934762
#> [1] "Rough volume of good region is 1.63654373893493e-07"

summary(example_results)
#> An EPC - alpha model was ran where those paremeters had an linearrelationship with the environment
#> 
#> Maximum likelihood estimate:
#> 167.3114692715
#> 
#> Number of Parameters:
#> 5
#> 
#> AIC:
#> 344.622938543 An EPC - theta model was ran where those paremeters had an linearrelationship with the environment
#> 
#> Maximum likelihood estimate:
#> 167.3114692715
#> 
#> Number of Parameters:
#> 5
#> 
#> AIC:
#> 344.622938543
#> 
#> Parameters Estimates: 
#>                       b0_t     b1_t      sig2        b0_a       b1_a
#> best             0.5526942 2.045722 0.4574184 0.022423777 0.02678971
#> lower.CI         0.4632926 2.023020 0.3921361 0.009373159 0.02284182
#> upper.CI         0.7026640 2.068737 0.5214812 0.027144672 0.02934762
#> lowest.examined  0.4047647 1.505148 0.3455309 0.008933104 0.02005322
#> highest.examined 0.7129060 2.565388 0.6783171 0.028742192 0.03762531
```

Looking at the summarized output, we can see that the maximum likelihood
estimate did indeed recover the linear relationship of the environment
on alpha and theta. However, how confident can we be in these paremeter
estimates? To solve this issue, ‘epcTools’ uses ‘dentist’ to generate
the 95% CI around each of the ML parameter estimates which we can see in
the summary function. However, how do these 95% CI estimates look on a
likelihood surface? Are the parameter estimates on a ridge (not
desirable) or on a well defined peak (desirable). Let’s see what these
look like:

## Dentist parameter visualization

``` r
#Dentist output
plot(example_results)
```

<img src="man/figures/README-visualization-1.png" width="100%" />

Well, some parameters certainly look more peak like than others but that
is expected. Sig2 and alpha are known to have a ridge like relationship
with each other and this is true for the constituent EPC parameters that
go into either of these (b0a and b1a here). Notice, however, that the
slope parameter (b1a) has a much more defined peak than the intercept
parameter (b0a) which suggests support for a linear relationship. To
really understand if this is real however, we need to compare the AIC of
the EPC-AT model to a single peak OU or BM model. Here, we will just do
OU to save time.

## EPC-AT vs. OU comparison

``` r
#First, we have to change the cache epc params to be "OU" to run an OU search. Then we will skip dentist as we only want the AIC output.
sim0.5f_200b_at_linear[[1]]$epc_params<-"OU"
OU_example<-epc.max.lik(sim0.5f_200b_at_linear[[1]],skip_dentist = TRUE)
#> Starting initial iterative maximum likelihood search...
#> First search     AIC:767.26
#> Starting second iterative maximum likelihood search...
#> Second search    AIC:767.26
#> Starting third iterative maximum likelihood search...
#> Third search     AIC:767.26
#> Starting last iterative maximum likelihood search...
#> Last search  AIC:767.26


cat(paste0("EPC-AT: ",example_results$fit$AIC, "OU: ",OU_example$fit$AIC))
#> EPC-AT: 344.622938543OU:  EPC-AT: 344.622938543OU:
```

As you can see, the EPC-AT model has a much lower AIC than the OU model,
indiciating that there is strong support for an EPC-environment
relationship on the trait.
