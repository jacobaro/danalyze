# utility functions

# internal function to remove backticks from a string
remove_backticks = function(str) {
  r = sapply(str, function(x) { if(substring(x, 1, 1) == "`") x = substring(x, 2); l = nchar(x); if(substring(x, l, l) == "`") x = substring(x, 1, l - 1); x })
  r = as.character(r)
  r
}

# identify time varying covariates for survival model
# more info: https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf + https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6015946/
identify_tve = function(formula, data) {
  data$time.step = dplyr::if_else(data$stop < 2.172603, 0, 1)
  formula = survival::Surv(start, stop, status) ~ sex + log1p(age) + rate

  # run model
  mx = survival::coxph(survival::Surv(start, stop, status) ~ sex + pspline(log1p(age)) + rate, data, cluster = id)
  summary(mx)
  mxt = timecox(survival::Surv(start, stop, status) ~ const(sex) + log1p(age) + rate, data, clusters = data$id)
  summary(mxt)

  # check cox.zph
  mx.zph = survival::cox.zph(mx)

  # plot
  plot(mx.zph, resid = F)
  plot(mxt)

  # temp data frame
  dx = tibble::as_tibble(mx.zph$y)
  dx$.time = mx.zph$time

  # get breakpoints
  br.time = dx$.time[strucchange::Fstats(formula = sex ~ .time, data = dx)$breakpoint]
}

# fast identification of baseline hazard
fast_bh = function(model) {
  # get survfit
  sfit = survival::survfit(model, se.fit = F)

  # remove means
  chaz = sfit$cumhaz * exp(-sum(model$means * model$coefficients))

  # return
  return(list(hazard = chaz, time = sfit$time))
}

# TODO: bring in some functions from performance package to test model assumptions
