# run to recreate roxygen stuff
# devtools::document()

# run tests for various glm models
test_glm = function() {
  ## SETUP DATA

  # load library
  library(danalyze)

  # load data
  data(lalonde.psid, package = "causalsens")

  # set the formula
  m = treat ~ education * nodegree + black + hispanic + log1p(age) + married + log1p(re75) + log1p(re74) * nodegree
  f = re78 ~ treat
  f.full = update(m, re78 ~ treat + .)
  f.nb = round(re78) ~ treat

  # select data and include only complete cases (optional)
  dt = dplyr::filter(dplyr::select(lalonde.psid, all.vars(m), all.vars(f)), complete.cases(lalonde.psid))

  # get weights
  dt$.weights = CBPS::CBPS(m, dt)$weights


  ## RUN STANDARD VERSION

  # run various glm versions with and without weights
  summary(glm(f, dt, family = gaussian))
  summary(glm(f, dt, family = gaussian, weights = dt$.weights))
  summary(glm(f.full, dt, family = gaussian))
  summary(MASS::glm.nb(formula = f.nb, data = dt, weights = dt$.weights))


  ## RUN THIS VERSION

  # run the analysis -- both with and without weights
  out.unweighted = analysis(runs = 1000, formula = f, data = dt)
  out.full = analysis(runs = 1000, formula = f.full, data = dt)
  out.weighted = analysis(runs = 1000, formula = f, data = dt, weights = dt$.weights)

  # also automatically supports over-dispersed count data -- determined based on the nature of the data for dependent variable (integer, non-negative)
  out.nb = analysis(runs = 1000, formula = f.nb, data = dt, weights = dt$.weights)

  # create a hypothesis that we want to test -- in this case the impact of moving from a value of '0' to '1' for "treat"
  predictions = pr_list(treat = c(1, 0))

  # examine the results, which should be very close to the summaries above
  results(object = out.unweighted, predictions = predictions)
  results(object = out.weighted, predictions = predictions)
  results(object = out.full, predictions = predictions)
  results(object = out.nb, predictions = predictions)

  # results should be comparable to the above (contrasts show the treatment effect)


  ## MODEL CAN ALSO BE RUN USING BAYESIAN INFERENCE

  # run -- currently weights don't work as expected for glm models in the rstanarm package so run without
  out.bayes.full = analysis(runs = 1000, formula = f.full, data = dt, inference = "bayesian")

  # get results for bayesian version
  results(object = out.bayes.full, predictions = predictions)

  # results should match the full model run using lm and using the non-bayesian version of the model


  ## ADDITIONAL OPTIONS

  # also supports multiple predictions
  results(object = out.full, predictions = c(predictions, pr_list(age = c(44, 25))))

  # to make the analysis faster, select fewer draws

  # analyze with a lot of runs
  out.weighted = analysis(runs = 5000, formula = f, data = dt, weights = dt$.weights)

  # produce results faster by using a random sample of draws
  results(object = out.full, predictions = c(predictions, pr_list(age = c(44, 25))), draws = 500)
}

# run tests for survival analysis
test_survival = function() {
  ## LOAD DATA

  # load library
  library(danalyze)

  # try survival
  data(mela.pop, package = "timereg")

  # formula
  f = survival::Surv(start, stop, status) ~ sex + log1p(age) + rate


  ## RUN COXPH MODEL

  # run survival model
  m.surv = survival::coxph(f, mela.pop, cluster = id, x = T)
  summary(m.surv)


  ## RUN PACKAGE VERSION

  # run
  out = analysis(runs = 1000, formula = f, data = mela.pop, cluster = ~ id)
  out.bayes = analysis(runs = 500, formula = f, data = mela.pop, inference = "bayesian")

  # also try with an exponential baseline hazard
  out.bayes.exp = analysis(runs = 500, formula = f, data = mela.pop, inference = "bayesian", model.extra.args = list(basehaz = "exp"))

  # create prediction list
  predictions = pr_list(sex = c(2, 1))

  # get results -- should be comparable across methods and the coefs should be similar to the above
  results(object = out, predictions = predictions, times = 1:4)
  results(object = out.bayes, predictions = predictions, times = 1:4)
  results(object = out.bayes.exp, predictions = predictions, times = 1:4)

  # the bayesian version uses a different approach for the baseline hazard (and has priors) so results will vary slightly
}

# runs a test for mediation analysis -- compares output from this package to 'mediation'
test_mediation = function() {
  ## SETUP DATA

  # load library
  library(danalyze)

  # load data
  data(jobs, package = "mediation")


  ## RUN MEDIATION PACKAGE VERSION

  # mediation models
  m.m = lm(job_seek ~ treat + econ_hard + sex + age, data = jobs)
  m.y = lm(depress2 ~ treat + job_seek + econ_hard + sex + age, data = jobs)

  # run -- indirect effect: first difference for mediator caused by treatment holding treatment constant; direct effect: first difference for treatment holding mediator constant
  m.mediate = mediation::mediate(m.m, m.y, sims = 1000, treat = "treat", mediator = "job_seek")

  # show results
  summary(m.mediate)


  ## RUN DANALYZE PACKAGE VERSION

  # mediation models
  md.m = analysis(runs = 2000, formula = job_seek ~ treat + econ_hard + sex + log1p(age), data = jobs)
  md.y = analysis(runs = 2000, formula = depress2 ~ treat + job_seek + econ_hard + sex + log1p(age), data = jobs)

  # set predictions
  m.predictions = pr_list(treat = c(1, 0))

  # get mediation results
  out.med = results_mediation(m.mediator = md.m, m.outcome = md.y, predictions = m.predictions, draws = 500)

  # show results
  out.med

  # p-values and effect sizes should be very similar (direct.effect == ADE, indirect.effect == ACME, etc.) -- in particular the "indirect.effect" should fully capture uncertainty in both models


  ## ALSO COMPARE BAYESIAN VERSION

  # bayesian mediation models
  md.bayes.m = analysis(runs = 2000, formula = job_seek ~ treat + econ_hard + sex + log1p(age), data = jobs, inference = "bayesian")
  md.bayes.y = analysis(runs = 2000, formula = depress2 ~ treat + job_seek + econ_hard + sex + log1p(age), data = jobs, inference = "bayesian")

  # show effect
  out.med.bayes = results_mediation(m.mediator = md.bayes.m, m.outcome = md.bayes.y, predictions = m.predictions, draws = 500)
  out.med.bayes

  # results should essentially be the same as the non-bayesian version
}
