
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
#> First search     AIC:350.43
#> Starting second iterative maximum likelihood search...
#> Second search    AIC:350.43
#> Starting third iterative maximum likelihood search...
#> Third search     AIC:350.43
#> Starting last iterative maximum likelihood search...
#> Last search  AIC:350.43
#> [1] "Calculating intervals at a confidence level of 95%"
#> [1] "Done replicate 500"
#> [1] "CI of values (the 29 replicates within 5.53524884675818 neglnL of the optimum)"
#>        neglnL parameter_1 parameter_2 parameter_3 parameter_4 parameter_5
#> [1,] 170.2129   0.8391762    1.976150   0.4294995  0.02183338  0.02946018
#> [2,] 175.6873   1.0248757    2.018695   0.5279465  0.03164273  0.03506429
#> [1] "Rough volume of good region is 4.27577142350907e-08"
#> [1] "Done replicate 1000"
#> [1] "CI of values (the 143 replicates within 5.53524884675818 neglnL of the optimum)"
#>        neglnL parameter_1 parameter_2 parameter_3 parameter_4 parameter_5
#> [1,] 170.2129   0.8391762    1.965736   0.4294995  0.01476198  0.02862064
#> [2,] 175.6873   1.0924705    2.019601   0.5351442  0.03164273  0.03535046
#> [1] "Rough volume of good region is 1.6374662653123e-07"
#> [1] "Done replicate 1500"
#> [1] "CI of values (the 369 replicates within 5.53524884675818 neglnL of the optimum)"
#>        neglnL parameter_1 parameter_2 parameter_3 parameter_4 parameter_5
#> [1,] 170.2129   0.8391762    1.965736   0.4294995  0.01177903  0.02862064
#> [2,] 175.7282   1.0924705    2.019601   0.5351442  0.03164273  0.03668197
#> [1] "Rough volume of good region is 2.30804415278014e-07"
#> [1] "Done replicate 2000"
#> [1] "CI of values (the 528 replicates within 5.53524884675818 neglnL of the optimum)"
#>        neglnL parameter_1 parameter_2 parameter_3 parameter_4 parameter_5
#> [1,] 170.2129   0.8391762    1.965736   0.4285947  0.01133985  0.02862064
#> [2,] 175.7403   1.0924705    2.019601   0.5351442  0.03164273  0.03668197
#> [1] "Rough volume of good region is 2.37927846133854e-07"
#> [1] "Done replicate 2500"
#> [1] "CI of values (the 676 replicates within 5.53524884675818 neglnL of the optimum)"
#>        neglnL parameter_1 parameter_2 parameter_3 parameter_4 parameter_5
#> [1,] 170.2129   0.8391762    1.965216   0.4285947 0.002576456  0.02862064
#> [2,] 175.7403   1.0924705    2.019601   0.5351442 0.031642727  0.03778539
#> [1] "Rough volume of good region is 3.9099059819795e-07"

summary(example_results)
#> An EPC - alpha model was ran where those paremeters had an linearrelationship with the environment
#> 
#> Maximum likelihood estimate:
#> 170.212942386153
#> 
#> Number of Parameters:
#> 5
#> 
#> AIC:
#> 350.425884772307 An EPC - theta model was ran where those paremeters had an linearrelationship with the environment
#> 
#> Maximum likelihood estimate:
#> 170.212942386153
#> 
#> Number of Parameters:
#> 5
#> 
#> AIC:
#> 350.425884772307
#> 
#> Parameters Estimates: 
#>                       b0_t     b1_t      sig2        b0_a       b1_a
#> best             0.8560128 2.007563 0.4949117 0.031123088 0.03063689
#> lower.CI         0.8391762 1.965216 0.4285947 0.002576456 0.02862064
#> upper.CI         1.0924705 2.019601 0.5351442 0.031642727 0.03778539
#> lowest.examined  0.7011635 1.412827 0.3081205 0.002215806 0.01969982
#> highest.examined 1.1557557 2.441308 0.6098790 0.037081521 0.04244923
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
#EPC phenogram
plot(example_results)
```

<img src="man/figures/README-visualization-1.png" width="100%" />

As you can see, there is some nuance to the relationships each variable
has with each other and the roughness of the likelihood surface.
However, each variable exists in a clearly defined peak.
