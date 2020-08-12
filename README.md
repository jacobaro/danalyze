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
for full functionality: lme4, MASS, rstanarm, survival, tidyverse, and timereg.

In addition, to run a Bayesian version of survival analysis it is necessary to install a development version of rstanarm.
To do this run: 

  ``` r
  devtools::install_github("stan-dev/rstanarm", ref = "feature/survival", build_vignettes = F)
```

## Package architecture

### Purpose

The ultimate goal of the 'danalyze' package is to provide a suite of functions that make it easier to turn qualitative information about 
a case, situation, or issue into quantitative results produced by best-practice causal methods.

Currently, all of the main components of  the package are in place. Work will continue to improve the quantitative methodologies and the 
generalizability of the included methods.

### Main components

A user will generally have qualitative information about a case that they would like to assess. For example, actions whose effect is plausible 
but unknown, other factors that have some influence, beliefs about the conditions under which an action may be more or less effective, and 
a theory of how an action works. (Data on these factors must, of course, already be collected!)

This qualitative knowledge maps to four different categories of variable: treatment, control, interaction, and mediation. Lets take a situation where 
the is interested in personal earnings and the impact of job training (the classic 'lalonde' data).

First load the package and the data:

```r
  # load library
  library(danalyze)

  # load data
  data(lalonde.psid, package = "causalsens")

  # set the data
  dt = lalonde.psid

  # make race a categorical variable -- doesnt need to happen just makes things a little neater
  dt$race = dplyr::case_when(
    dt$black == 1 ~ "Black",
    dt$hispanic == 1 ~ "Hispanic",
    T ~ "Other"
  )
```

Next, the user enters their qualitative knowledge by category and provides labels:

```r
  # identify main variables
  treatment = c(
    "Job training" = "treat"
  )

  # identify interactions
  interaction = c(
    "Race" = "race",
    "Married" = "married",
    "No Degree" = "nodegree"
  )

  # other controls
  control = c(
    "Age" = "age",
    "Education" = "education",
    "Race" = "race",
    "Married" = "married",
    "No Degree" = "nodegree",
    "Earnings in 1975" = "re75",
    "Earnings in 1974" = "re74"
  )
```

With this information, the package can create a research plan:

```r
  # create a plan -- not terribly interesting since the data is so simple
  plan = research_plan(
    treatment = treatment,
    interaction = interaction,
    control = control,
    data = dt
  )

  # set outcome as extra for now -- will not be needed eventually
  plan$outcome = .outcome ~ .
  plan$outcome.occurrence = list("Earnings" = dt$re78)
```

The qualitative knowledge is entered (although no knowledge about possible mediating factors was provided). Eventually, 
the outcome will also be entered in a more convenient manner.

The 'research_plan' function creates all the necessary formulas for testing, proper transformations for variables, and 
predictions to test various effects. All of this is determined by the category of the variable and the variable's 
underlying data.

With the plan created, the user can know analyze this plan. The following command is used with the plan passed to the 
function and the return from the function saved to the variable 'results':

```r
  # now analyze the plan
  results = analyze_plan(research.plan = plan)
```

The results can take all the results from statistical assessment of the models identified in the research plan. The 
function automatically identifies the correct model to use, provides useful model defaults (when necessary), selects 
necessary variables in a formula to include, runs the model, and aggregates the results in an easy to understand manner.

For instance, the predictions and contrasts returned for the main effect are as follows:

```r
# PREDICTIONS:
#   .outcome           .main.variable .prediction.id treat      c  c.low c.high p.value draws
#   <chr>              <chr>                   <dbl> <dbl>  <dbl>  <dbl>  <dbl>   <dbl> <dbl>
# 1 Earnings           treat                       1     1 23056. 21146. 24816.       0   500
# 2 Earnings           treat                       2     0 20329. 19955. 20744.       0   500

# CONTRASTS:
#   .outcome           .main.variable .contrast .prediction.id     c c.low c.high p.value draws v1.high.p1 v1.low.p1
#   <chr>              <chr>          <chr>     <chr>          <dbl> <dbl>  <dbl>   <dbl> <dbl> <chr>      <chr>    
# 1 Earnings           treat          1 vs. 0   1, 2           2761.  767.  4638.   0.008   500 1          0        
```

These results show the predicted impact of treatment at a value of '1' and '0' (median effect, low and high confidence 
interval, and statistical significance). The contrasts show the statistical difference between the two predictions. In this 
case the job training program increases wages by $2,761 and the effect is statistically significant at conventional levels.

Understanding these statistical results may be a challenge or there may be so many statistical results that it is nearly 
impossible to make sense of them.

To address this problem, the package also provides tools to provide a human-interpretable understanding of the results. To 
do so the user would pass these results to the function 'qualitative_assessment' as follows:

```r
  # produce qualitative assessment -- will return NULL if no variables or conditional effects are significant
  assessment = qualitative_assessment(research.plan = plan, all.results = results)
```

The returned assessment provides a textual descriptio of all of the effects. In this case, the following text is returned;

```r
  `Earnings, positive`
  [1] "The outcome 'Earnings' is increased by one variable."                                                                                                                          
  [2] "Variable 'Job training' has a positive unconditional effect on 'Change in earnings' and a 'positive' conditional effect 
       when 'No Degree' is 'high' and/or 'Married' or 'Race' are 'low'."
```

This function can also handle time and mediating variables if present. If there are many possible treatments, the function can also rank 
order these effects to provide some understanding of which variables are most (or least important) and under what conditions they have 
an effect.

Eventually this system will also allow the user to specify the situation under which the ranking should occur (e.g., the variable that 
has the largest impact on earnings for married individuals).


## Example

This is a basic example to show how the analysis functions are used:

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
  
  # examine the results -- the function will test the hypothesis in "predictions"
  results(out, predictions)
  
```
