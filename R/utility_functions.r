# utility functions

#' Remove backticks.
#'
#' This function takes in a string and outputs a string with
#' backticks removed.
#'
#' @param str The string to remove backticks from.
#' @return A string with backticks removed.
#' @export
#'
remove_backticks = function(str) {
  r = sapply(str, function(x) { if(substring(x, 1, 1) == "`") x = substring(x, 2); l = nchar(x); if(substring(x, l, l) == "`") x = substring(x, 1, l - 1); x })
  r = as.character(r)
  r
}

#' Fast identification of baseline hazard.
#'
#' This function takes in a survival model and outputs a list with the
#' baseline hazard and time.
#'
#' @param model A survival model running using the "survival" package.
#' @return A list with the baseline hazard and time.
#' @export
#'
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

#' Symmetric 'log1p'.
#'
#' This function takes in a positive or negative number and returns correctly
#' transformed values.
#'
#' @param x A numeric vector.
#' @return The transformed numeric vector.
#' @export
#'
symlog = function(x) { sign(x) * log1p(abs(x)) }

#' @rdname symlog
#' @export
#'
symsqrt = function(x) { sign(x) * sqrt(abs(x)) }

# TODO: bring in some functions from performance package to test model assumptions

# intersection of vectors
intersect_all = function(...) {
  Reduce(intersect, list(...))
}

#' Collapse time.
#'
#' Summarize a dataframe to different units of analysis. The 'time.var' is
#' divided by 'time.factor' and turned into an integer.
#'
#' @param data The data to summarize.
#' @param time.var Symbols identifying the time variable in data.
#' @param time.factor A number that 'time.var' will be divided by.
#' @param group.by A vector of symbols that identify the unit of analysis.
#' @param variables A list of variables to apply special functions to.
#' @return A summarized dataframe.
#' @export
#'
collapse_time = function(data, time.var, time.factor, group.by, variables) {
  # function to collapse time -- could really recode this to do start and end time using blocks of X time that would be neat

  # first make sure it is not grouped
  data = dplyr::ungroup(data)

  # set :=
  `:=` = rlang::`:=`

  # enquo
  time.var = rlang::enquo(time.var)
  group.by = rlang::enquo(group.by)

  # create the new time group
  data = dplyr::mutate(data, .time.factor = ceiling(!!time.var / time.factor))

  # create time per group
  data.group = dplyr::mutate(dplyr::group_by_at(data, dplyr::vars(!!group.by)), .time.group = !!time.var - min(!!time.var, na.rm = T) + 1)

  # grouped data
  data.group = dplyr::group_by_at(data.group, dplyr::vars(!!group.by, .time.factor))

  # summarize
  data.base = dplyr::summarize(data.group,
                               !!time.var := min(!!time.var, na.rm = T),
                               .time.start = min(.time.group, na.rm = T) - 1,
                               .time.end = max(.time.group, na.rm = T), .groups = "keep")

  # get column order
  column.order = c(colnames(data.base), colnames(data)[!colnames(data) %in% colnames(data.base)])

  # summarize data

  # identify user-supplied variable names
  user.var.names = unique(unlist(lapply(variables, function(x) colnames(dplyr::select(data, !!!x$vars)))))

  # make sure we are not trying to transform grouping or time variables
  user.var.names = user.var.names[!user.var.names %in% c(dplyr::group_vars(data.group), ".time.factor", ".time.group", ".time.start", ".time.end")]

  # columns to remove
  column.remove = colnames(dplyr::select(dplyr::ungroup(data.group), !!time.var, !!group.by, .time.factor, .time.group))

  # apply our special functions
  data.special = lapply(variables, function(x) {
    dplyr::select(dplyr::ungroup(dplyr::summarize_at(data.group, .vars = colnames(dplyr::select(data, !!!x$vars)), .funs = x$func)), -tidyselect::any_of(column.remove))
  })

  # bind
  data.special = dplyr::bind_cols(data.special)

  # remaining variables to select
  data.group.remain = dplyr::select(data.group, -!!user.var.names)

  # numeric data
  data.numeric = dplyr::select(dplyr::ungroup(dplyr::summarize_if(data.group.remain, is.numeric, mean, na.rm = T)), -tidyselect::any_of(column.remove))

  # categorical data
  data.categorical = dplyr::select(dplyr::ungroup(dplyr::summarize_if(data.group.remain, function(x) !is.numeric(x), function(x) dplyr::first(na.omit(x)))), -tidyselect::any_of(column.remove))

  # combine
  data.full = dplyr::bind_cols(lapply(list(data.base, data.special, data.numeric, data.categorical), dplyr::ungroup))

  # select columns
  data.full = dplyr::select(data.full, !!column.order)

  # return
  return(data.full)
}

#' Create high and low values
#'
#' This function takes in a vector (numeric or character) and outputs nicely formatted
#' test values.
#'
#' @param x A numeric vector.
#' @param .quantile The high and low quantiles for a numeric vector.
#' @param .places The number of places to round a numeric vector to.
#' @return For a numeric vector this function return high and low values.
#'   For a character vector this returns unique values.
#' @export
#'
create_values = function(x, .quantile = c(0.975, 0.025), .places = 2) {
  # make sure quantile is okay
  if(!is.numeric(.quantile) | !length(.quantile) == 2) {
    warning("Quantile is not structured correctly. Using default.")
    .quantile = c(0.975, 0.025)
  }

  # check
  if(is.numeric(x)) {
    # if there are a large number of zeros than we want to produce values absent zero
    if((length(x == 0) / length(x)) > 0.75) {
      x.t = x[x != 0]
    } else {
      x.t = x
    }

    # get quantile
    r = quantile(x.t, .quantile, na.rm = T)

    # make sure the values are okay
    if(diff(range(r)) < .Machine$double.eps ^ 0.5) {
      r = c(max(x), min(x))
    }

    # round
    r = round(r, .places)
  } else {
    r = unique(x)
  }

  # return
  return(r)
}
