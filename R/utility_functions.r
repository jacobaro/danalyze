# utility functions

# internal function to remove backticks from a string
remove_backticks = function(str) {
  r = sapply(str, function(x) { if(substring(x, 1, 1) == "`") x = substring(x, 2); l = nchar(x); if(substring(x, l, l) == "`") x = substring(x, 1, l - 1); x })
  r = as.character(r)
  r
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

# turn a formula into a character - surprisingly annoying!
as.character.formula = function(f) {
  paste(trimws(deparse(f)), collapse = " ")
}

# transformation functions

#' log1p for positive and negative values
#'
#' @export
#'
symlog = function(x) { sign(x) * log1p(abs(x)) }

#' sqrt for positive and negative values
#'
#' @export
#'
symsqrt = function(x) { sign(x) * sqrt(abs(x)) }

# TODO: bring in some functions from performance package to test model assumptions

# intersection of vectors
intersect_all = function(...) {
  Reduce(intersect, list(...))
}
