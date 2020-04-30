# utility functions

# internal function to remove backticks from a string
remove_backticks = function(str) {
  r = sapply(str, function(x) { if(substring(x, 1, 1) == "`") x = substring(x, 2); l = nchar(x); if(substring(x, l, l) == "`") x = substring(x, 1, l - 1); x })
  r = as.character(r)
  r
}

# identify time varying covariates for survival model
# more info: https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
identify_tve = function(formula, data) {
  data$time.step = dplyr::if_else(data$stop < 2.172603, 0, 1)
  formula = survival::Surv(start, stop, status) ~ sex:strata(time.step) + log1p(age) + rate

  # run model
  mx = survival::coxph(formula, data)

  # check cox.zph
  mx.zph = survival::cox.zph(mx)

  # plot
  plot(mx.zph, resid = F)

  # temp data frame
  dx = tibble::as_tibble(mx.zph$y)
  dx$.time = mx.zph$time

  # get breakpoints
  br.time = dx$.time[strucchange::Fstats(formula = sex ~ .time, data = dx)$breakpoint]
}
