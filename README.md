# danalyze

The goal of 'danalyze' is to provide a suite of functions to allow semi-automated causal analysis.

## Installation

You can install 'danalyze' from github with:


  ``` r
# install.packages("devtools")
devtools::install_github("jacobaro/danalyze")
```

## Example

This is a basic example to show the workflow:

  ``` r
  # load data
  data(lalonde.psid, package = "causalsens")

  # set formulas
  f = re78 ~ treat
  
  # select data
  dt = dplyr::select(lalonde.test, all.vars(f.f))
  
  # create prediction list
  predictions = pr_list(treat = c(1, 0))
  
  # run analysis
  out = analysis(runs = 1000, formula = f, data = dt)
  
  # get results
  results(m.out, predictions)
```
