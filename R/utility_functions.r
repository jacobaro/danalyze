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

# function to collapse time -- could really recode this to do start and end time using blocks of X time that would be neat

#' Summarize a dataframe to a different unit of analysis.
#'
#' @export
#'
collapse_time = function(data, time.var, time.factor, group.by, variables) {
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
