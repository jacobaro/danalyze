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

# runs a test to compare the cumulative incidence of a survival model to a series of logits
test_cumulative_logit = function() {
  ## LOAD DATA

  # load library
  library(danalyze)

  # get survival data
  data(mela.pop, package = "timereg")

  # formula -- the data has fractional stop times, which wont work for our purposes so just use start + 1
  mela.pop$stop2 = mela.pop$start + 1

  # and set the actual formula
  f = survival::Surv(start, stop2, status) ~ sex + log1p(age) + rate


  ## RUN COXPH MODEL

  # run survival model
  m.surv = survival::coxph(formula = f, data = mela.pop, cluster = id, x = T)
  summary(m.surv)

  m.surv.timereg = timereg::timecox(formula = f, data = mela.pop)

  # check proportional hazards
  cox.zph(m.surv) # standard test says that they are all fine
  summary(m.surv.timereg) # this test says that 'sex' is fine but that rate is not, which matches what we find with our linear model
  plot(m.surv.timereg) # plot shows that hazards are not proporitonal -- sex looks good but rate shows a rise and then a huge dip as does our linear model

  # run package version
  out = analysis(runs = 1000, formula = f, data = mela.pop, cluster = ~ id)

  # create a basic set of predictions for the normal survival model
  pr.basic = bind_rows(tibble(sex = 2, age = mean(mela.pop$age), rate = mean(mela.pop$rate), start = 0, stop2 = 1:10, status = 0),
                       tibble(sex = 1, age = mean(mela.pop$age), rate = mean(mela.pop$rate), start = 0, stop2 = 1:10, status = 0))

  # get the survival rate over time
  r.basic = predict(m.surv, newdata = pr.basic, type = "survival", se.fit = T)

  # turn it into cumulative incidence (the CIF common in interpreting competing risks models)
  tibble(
    .time = c(1:10, 1:10),
    c = 1 - r.basic$fit,
    c.low = c - qnorm(0.975) * r.basic$se.fit,
    c.high = c + qnorm(0.975) * r.basic$se.fit
  )

  # these basic predictions (which don't handle contrasts in a straightforward way can be compared to package predictions)

  # predictions to run
  main.pred = pr(sex = c(2, 1)) # pr(rate = quantile(mela.pop$rate, c(0.9, 0.1)))
  # predictions from the package
  res.surv = results(object = out, predictions = main.pred, times = 1:10)

  # setup data for logit models

  # helper function
  within_time = function(x, t, outcome = 1) {
    maxl = length(x)
    sapply(1:maxl, function(z) {
      maxt = min(maxl, z + t)
      if(any(na.omit(x[z:maxt]) == outcome)) return(outcome) # we know it happened
      if(z + t > maxt || any(is.na(x[1:maxt]))) return(NA) # these are outcomes we dont know -- they havent happened yet but could in a time we dont record data for
      return(0) # we know it didnt happen
    })
  }

  # we will create 10 dependent variables that record the occurrence of the outcome within X time ('status_c0' is identical to 'status')
  mela.pop = mela.pop %>% arrange(id, start, stop2) %>% group_by(id) %>%
    mutate(across(status, .fns = list(c0 = ~ within_time(.x, 0), c1 = ~ within_time(.x, 1), c2 = ~ within_time(.x, 2), c3 = ~ within_time(.x, 3), c4 = ~ within_time(.x, 4),
                                      c5 = ~ within_time(.x, 5), c6 = ~ within_time(.x, 6), c7 = ~ within_time(.x, 7), c8 = ~ within_time(.x, 8), c9= ~ within_time(.x, 9))))

  # run a linear model (on binomial data)
  res.logit = lapply(0:9, function(i) {
    # set the formula
    f.t = update(f, as.formula(paste0("status_c", i, " ~ .")))

    # run model -- binomial outcome but we are running it using a gaussian family -- could also just use lm but this makes it easier to check what happens when you swap family
    out.t = glm(formula = f.t, data = mela.pop, family = gaussian)

    # get prediction -- the SE is pretty large when using just id, but that is likely the more accurate assessment
    r.t = get_prediction_frequentist(out.t, cluster = ~ id, predictions = main.pred)

    # save
    r.t
  })

  # check
  out.logit = map_dfr(res.logit, "contrasts", .id = ".time")
  out.logit$.time = as.numeric(out.logit$.time)
  out.logit

  # compare to the survival model
  res.surv$contrasts

  # set plot data
  out.logit$.type = "Cumulative Linear"
  res.surv$contrasts$.type = "Survival"
  out.all = bind_rows(out.logit %>% select(.type, .time, c:c.high), res.surv$contrasts %>% select(.type, .time, c:c.high))

  # plot
  ggplot(out.all, aes(x = .time, y = c, ymin = c.low, ymax = c.high)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey80", size = 1.5) +
    geom_ribbon(aes(fill = .type), alpha = 0.4) + geom_line(aes(color = .type), size = 1) +
    scale_fill_manual("Model Type", values = c("orangered", "steelblue4")) + scale_color_manual("Model Type", values = c("orangered", "steelblue4")) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal() + labs(x = "Time", y = "Change in Cumulative Survival, Male vs. Female")

  # the lines are very similar -- this only holds as long as proportional hazards is not violated
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
