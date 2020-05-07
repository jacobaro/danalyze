# run to recreate roxygen stuff
# devtools::document()

# test the functions

# test
test_analysis = function() {
  # # library
  # needs::needs(tidyverse)

  # get matching functions
  source("../foreign_aid/cbps_utilities.r")

  # set data
  data(lalonde.psid, package = "causalsens")
  data.test = lalonde.psid

  # set formulas
  f = re78 ~ treat
  f.full = re78 ~ treat + education * nodegree + black + hispanic + log1p(age) + married + log1p(re75) + log1p(re74) * nodegree

  # get weights
  cbps = cbps_match(update(f.full, treat ~ . - treat), data.test, lasso = F)

  # select data
  dt = dplyr::select(data.test, all.vars(f.full))

  # set weights
  dt$.weights = cbps$weights

  # set weights
  # filter out zero weights
  dt = dplyr::filter(dt, .weights > 0)

  # check rstan version
  # cores = 4, warmup = 250,
  system.time(m.stan <- rstanarm::stan_glm(f.full, data = dt, prior = rstanarm::normal(scale = 2.5), prior_intercept = rstanarm::normal(scale = 2.5), cores = 8, iter = 1000, warmup = 250))
  summary(m.stan, probs = c(0.025, 0.5, 0.975))
  summary(lm(f.full, dt))

  # create prediction list
  predictions = pr_list(treat = c(1, 0))

  # run
  system.time(out.stan <- analysis(runs = 1500, formula = f.full, data = dt, inference = "bayesian", model.extra.args = list(cores = 6, warmup = 250)))
  out.nonweight = analysis(runs = 1000, formula = f.full, data = dt)
  out.weight = analysis(runs = 1000, formula = f, data = dt, weights = dt$.weights)

  # get predictions
  results(m.stan, predictions)
  results(out.stan, predictions)
  results(out.nonweight, predictions)
  results(out.weight, predictions)

  # compare

  # unweighted
  summary(lm(f, dt))
  summary(lm(f.full, dt))
  summary(MASS::glm.nb(f.full, dt))

  # weighted
  summary(lm(f, dt, weights = dt$.weights))
  summary(lm(f.full, dt, weights = dt$.weights))
  summary(MASS::glm.nb(f, dt, weights = dt$.weights))




  # try survival
  data(mela.pop, package = "timereg")
  data.test = mela.pop

  # formula
  f = survival::Surv(start, stop, status) ~ sex + log1p(age) + rate

  # run survival model
  m.surv = survival::coxph(f, data.test, cluster = id)
  summary(m.surv)
  summary(exp(-predict(m.surv, type = "expected"))) # survival rate

  # stan version
  m.stan = rstanarm::stan_surv(f, data = data.test, prior = rstanarm::normal(scale = 1), prior_intercept = rstanarm::normal(scale = 1), warmup = 250, iter = 1000, cores = 4, basehaz = "exp")

  # create prediction list
  predictions = pr_list(rate = c(0.07, 0))

  # # run
  # out = analysis(runs = 1000, formula = f, data = data.test, cluster = ~ id, weights = NULL)

  # # check quick summary
  # out$stan_summary

  # get results
  results(m.stan, predictions, times = 1:4)
  # results(out, predictions)
}

# runs a test for mediation analysis -- compares output from this package to 'mediation'
test_mediation = function() {
  # load data
  data(jobs, package = "mediation")

  ## RUN MEDIATION PACKAGE VERSION

  # mediation models
  m.m = lm(job_seek ~ treat + econ_hard + sex + age, data = jobs)
  m.y = lm(depress2 ~ treat + job_seek + econ_hard + sex + age, data = jobs)

  # run -- indirect effect: first difference for mediator caused by treatment holding treatment constant; direct effect: first difference for treatment holding mediator constant
  m.mediate = mediation::mediate(m.m, m.y, sims = 2000, treat = "treat", mediator = "job_seek", boot = T)

  # show results
  summary(m.mediate)


  ## RUN DANALYZE PACKAGE VERSION

  # mediation models
  md.m = danalyze::analysis(runs = 2000, formula = job_seek ~ treat + econ_hard + sex + age, data = jobs)
  md.y = danalyze::analysis(runs = 2000, formula = depress2 ~ treat + job_seek + econ_hard + sex + age, data = jobs)

  # set predictions
  m.predictions = danalyze::pr_list(treat = c(1, 0))

  # get mediation results
  out.med = danalyze::results_mediation(m.mediator = md.m, m.outcome = md.y, predictions = m.predictions)
}
