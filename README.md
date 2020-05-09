# Package: 'danalyze'

The goal of 'danalyze' is to provide a suite of functions to allow semi-automated causal analysis.

## Installation

You can install 'danalyze' from github using the following code:

  ``` r
# if devtools package is not installed run: install.packages("devtools")
devtools::install_github("jacobaro/danalyze")

```

## Dependencies

Most of the underlying models run come from different packages. Make sure the following packages are installed 
for full funcitonality: lme4, MASS, rstanarm, survival, tidyverse, and timereg.

In addition, to run a Bayesian version of survival analysis it is necessary to install a development version of rstanarm.
To do this run: 

  ``` r
  devtools::install_github("stan-dev/rstanarm", ref = "feature/survival", build_vignettes = F)
```

## Example

This is a basic example to show the workflow:

  ``` r
  # first load the data
  data(lalonde.psid, package = "causalsens")

  # set the formula
  f = re78 ~ treat
  
  # select data and include only complete cases (optional)
  dt = dplyr::filter(dplyr::select(lalonde.psid, all.vars(f)), complete.cases(lalonde.psid))
  
  # create a hypothesis that we want to test -- in this case the impact of moving from a value of '0' to '1' for "treat"
  predictions = pr_list(treat = c(1, 0))
  
  # run the analysis -- the function will figure out what models to run
  out = analysis(runs = 1000, formula = f, data = dt)
  
  # examine the results -- the funciton will test the hypothesis in "predictions"
  results(m.out, predictions)
  
```
