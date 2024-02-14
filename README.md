# Homework-7440

Homework for Spring 2024 PUBH 7440 at the University of Minnesota.

# HW1

* Statistical Concepts: 
  + review of Bayes Rule and application of the rule to work with conventional conditional probability examples like blindly drawing 
    balls of red/blue colors from two urns 
  + visualization of prior, posterior, and observed data likelihood distributions over some distributional parameters $\theta$. 
    Concepts applied to simple models with conventional exponential family distributions for visualizations. More complex MCMC examples in later HW assignments 
  + hypothesis testing for differences between two groups using assumed binomial model with beta conjugate prior distribution. 
    Main examples uses data driven Bayes method to find an appropriate set of parameters for the prior distribution. 
  + More examples in appendix show how the choice of more and more informative beta distribution prior yields different results 
  + Application of the poisson model for the difference in rates using gamma conjugate prior for the parameter $\lambda$. 
  
* Implementation via R code: 
  + All code to do data wrangling, hypothesis testing, and visualization of results was done using a mix of base R and tidyverse code 
  
* Packages: 
  + No new packages for Bayesian data analysis were introduced yet 
  
# HW3 

* Statistical Concepts: 
  + Gibbs sampling technique: impute suppressed values of small death counts for counties in PA and estimate posterior distribution of 
    parameter $\lambda$ for poisson distribution 
  + derivation of a full hierarchical model given data likelihood, prior distribution of the parameter, and model for suppressed/censored values of the data (which we need for the posterior distribution)
  + interpretation of the prior and posterior parameters of poisson and gamma distributions 
  
* Implementation via R code: 
  + iterative Gibbs sampling procedure for sampling suppressed/censored data and sampling parameters from the posterior distribution based 
    on imputed data 
  + mapping of age-adjusted death rates by county using `maps` objects 

* Packages: 
  + `maps`: update to `maptools`: easy framework for importing maps to R and plotting shapes/outlines 
    - example contains map of PA counties
