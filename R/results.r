# functions to get results from analysis

#' Function to get summary statistics from a dataframe.
#'
#' This function produces a table of summary statistics based on a formula and a dataframe.
#'
#' @family post analysis exploration
#' @param formula.list A formula or list of formulas that indicate which variables are included in summary statistics.
#' @param data Dataframe with variables to produce summary statistics from.
#' @param labels An optional vector of labels: c("variable" = "label).
#' @return A dataframe of summary statistics.
#' @keywords summary_stats, summary
#' @export
#' @examples
#' get_summary_statistics(formula.list = ~ variable, data = df, labels = c("variable" = "label"))
#'
get_summary_statistics = function(formula.list, data, labels = NULL) {
  # process formula
  f = lapply(formula.list, function(x) delete.response(terms(lme4:::nobars(x))))

  # get all variables
  f.all.vars = unique(unlist(lapply(f, all.vars)))

  # select frame
  d.f = dplyr::select(dplyr::ungroup(data), tidyselect::all_of(f.all.vars))

  # create model matrix preserving NAs
  d.mm = model.matrix.lm(~ 0 + ., d.f, na.action = na.pass)

  # get colnames -- do we want to do some reordering?
  if(!is.null(labels)) {
    suppressWarnings(col.names <- dplyr::left_join(tibble::tibble(colnames = colnames(d.mm)), tibble::tibble(values = names(labels), colnames = labels), by = "colnames"))
    col.names$values[is.na(col.names$values)] = col.names$colnames[is.na(col.names$values)]
    col.names = col.names$values
  } else {
    col.names = colnames(d.mm)
  }


  # convert to tibble
  d.mm = tibble::as_tibble(d.mm)

  # create the data frame
  d.mm = tibble::tibble(
    "Variable" = col.names,
    "N" = sapply(d.mm, function(x) length(x[!is.na(x)])),
    "Mean" = sapply(d.mm, mean, na.rm = T),
    "Std. Dev." = sapply(d.mm, sd, na.rm = T),
    "Min." = sapply(d.mm, min, na.rm = T),
    "Max." = sapply(d.mm, max, na.rm = T)
  )

  # format
  d.mm = dplyr::mutate_if(d.mm, is.numeric, function(x) format(round(x, 3), scientific = F, big.mark = ",", trim = T))

  # return
  return(d.mm)
}

#' Create predicitons for all variables in a formula
#'
#' This function generates a list of predictions based on the data for all variables in the 'variables' formula. This is an
#' alternative to showing coefficient estimates. Coefficient estimates can be hard to interpret or compare. This function
#' produces standardized effect estimates for all variables.
#'
#' For character (factor) variables, all unique values of the variable are compared against each other. For numeric variables, the
#' impact of moving from the 10th to the 90th percentile is estimated.
#'
#' @family post analysis exploration
#' @param variables A formula or list of formulas that indicate which variables are included in summary statistics.
#' @param data Dataframe with variables to produce summary statistics from.
#' @param labels An optional vector of labels: c("variable" = "label).
#' @param .variable.set An optional subset of variables to examine.
#' @param .perc For numeric variables this sets the low and high quantile values. The default is to examine the
#'   effect of moving from the 10th percentile (0.1) to the 90th percentile (0.9).
#' @return A dataframe that shows the impact of all variables (median, confidence interval, and P-value). If labels
#'   are provided this dataframe will also contain the labels.
#' @keywords coefficients, summary, effects
#'
get_variable_effects = function(variables, data, labels = NULL, .variable.set = NULL, .perc = c(0.1, 0.9)) {
  # function to get predictions for assessing variable effects

  # can accept either a formula or a vector
  if(rlang::is_formula(variables)) {
    # process formula
    f = all.vars(delete.response(terms(lme4:::nobars(variables))))
  } else {
    f = as.character(variables)
  }

  # select frame -- if the character vector is named it automatically renames too -- pretty cool
  d.f = dplyr::select(dplyr::ungroup(data), tidyselect::any_of(f))

  # reoder
  if(!is.null(labels)) {
    # select present
    labels.present = labels[labels %in% colnames(d.f)]
  } else {
    labels.present = colnames(d.f)
    names(labels.present) = colnames(d.f)
  }

  # shrink to the variables we want to examine
  if(!is.null(.variable.set)) {
    labels.present = labels.present[stringr::str_detect(labels.present, .variable.set)]
  }

  # create predictions
  r = lapply(1:length(labels.present), function(x) {
    # make sure we dont double name
    x = labels.present[x]

    # get values
    if(is.factor(d.f[[x]]) | is.character(d.f[[x]])) {
      # we have a factor

      # get levels
      if(is.factor(d.f[[x]])) {
        levels = rev(levels(d.f[[x]]))
      } else {
        levels = unique(na.omit(d.f[[x]]))
      }

      # cant seem to get rlang to work here so just do it the hard way
      args = list(x = levels)
      names(args) = x

      # return
      r = do.call(pr_list, args = args)
    } else {
      # we have a number

      # get low and high
      l = quantile(d.f[[x]], .perc[1], na.rm = T)
      h = quantile(d.f[[x]], .perc[2], na.rm = T)

      # are they the same? if so then just get low and high
      if(l == h) {
        l = min(d.f[[x]])
        h = max(d.f[[x]])
      }

      # cant seem to get rlang to work here so just do it the hard way
      args = list(x = c(round(h, 3), round(l, 3)))
      names(args) = x

      # return
      r = do.call(pr_list, args = args)
    }

    # set names since many high and low values can be the same across variables
    names(r) = paste0(names(labels.present)[labels.present == x], ": ", names(r))

    # and send it back
    return(r)
  })

  # combine
  r = c(unlist(r, recursive = F))

  # return
  return(r)
}

# function to identify observations in our data that match our results

#' Identify observations that support our results.
#'
#' This functions takes a set of results (main and interaction but not mediation) and identifies observations
#' in the dataset used that support the results. Supportive observations are only identified for statistically
#' significant results.
#'
#' For most normal regression models, an observation is supportive when (1) it occurs in the same unit of analysis
#' as a positive outcome and (2) the variable has a high value. For interactions, both the main variable and the interaction
#' need to have the specified values (conditional effects can occur when the interaction variable has a low value--this is
#' accounted for). For event history models, a supportive variable does not need to occur in the same observation as the outcome.
#' Observations that are close in time to the outcome are also considered supportive.
#'
#'
#' @family post analysis exploration
#' @param results The output from 'analyze_plan'.
#' @param research.plan The research plan used in the produce of 'results'.
#' @param unit The unit of analysis to aggregate examples within.
#' @param time The tine variable used in the results. Only necessary for time-series analysis.
#' @param .time.distance If a results was from a survival model, this parameter identifies the quantile
#'   of time that indicates a "close" event. If time to the outcome is below this quantile of the time
#'   variable, the observation is considered close.
#' @param .quantile The quantile cutoffs for low and high values of a variable. A value above the first
#'   cutoff is considered a "high" value while a value below the second cutoff is considered a "low' value.
#' @param .variables Names of the variable columns in the results. Defaults to main and interaction effects.
#'   If results contains more than two variables, additional values need to be added here.
#' @return Function returns a dataframe that contains a supportive observation that matches one of the statistically
#'   significant results.
#' @export
#'
identify_examples = function(results, research.plan, unit = NULL, time = NULL, .time.distance = 0.1, .quantile = c(0.8, 0.2), .variables = c(".main.variable", ".main.interaction")) {
  # basic logic is to find observations in the data where the outcome occurs and our main variable has a higher value
  # TODO: make this work for non numeric variables -- use a generic function that makes predictions for variables from danalyze
  # the results could include the reseaarch plan -- this might make sense since the plan was used for the results

  # vars to group by
  .variables = .variables[.variables %in% colnames(results)]
  var.to.group = c(".outcome", .variables)
  var.to.group = var.to.group[var.to.group %in% colnames(results)]
  results$.id = 1:nrow(results)

  # need to identify high and low values
  var.value.labels = dplyr::group_map(dplyr::group_by_at(results, var.to.group), ~ {
    r = sapply(1:length(.variables), function(x) {
      v = paste0("v", x, ".high")
      dplyr::if_else(.x[[v]] >= max(.x[[v]], na.rm = T), "High", "Low")
    })
    r = matrix(r, ncol = length(.variables))
    colnames(r) = paste0(.variables, ".status")
    r = tibble::as_tibble(r)
    r$.id = .x$.id
    r
  }, .keep = T)

  # bind
  var.value.labels = dplyr::bind_rows(var.value.labels)

  # add to results
  results = dplyr::left_join(dplyr::ungroup(results), var.value.labels, by = ".id")

  # identify significant results that we want to identify cases for
  results = dplyr::summarize(dplyr::group_by_at(dplyr::ungroup(results), c(var.to.group, colnames(dplyr::select(var.value.labels, -.id)))),
                             coefficient = mean(c, na.rm = T),
                             direction = dplyr::if_else(coefficient < 0, "negative", "positive"),
                             is.significant = dplyr::if_else(any(p.value < 0.05), 1, 0),
                             .groups = "keep")

  # filter
  results = dplyr::filter(results, is.significant == 1)

  # set data
  data = research.plan$data

  # arrange data
  data = dplyr::arrange(data, dplyr::across(tidyselect::any_of(c(unit, time))))

  # loop through results and find plausible supportive observations
  r = lapply(1:nrow(results), function(i) {
    # we may not have all the variables so subset
    .variables.t = .variables[sapply(.variables, function(x) !is.na(results[i, x]))]

    # values
    var.values = lapply(.variables.t, function(x) {
      danalyze:::create_values(data[[results[[i, x]]]], .quantile = .quantile)[if(results[i, paste0(x, ".status")] == "High") 1 else 2]
    })
    names(var.values) = .variables.t

    # need to deal with character values

    # identify low and high values of our variable
    variable.high = sapply(names(var.values), function(x) {
      if(results[i, paste0(x, ".status")] == "High") {
        dplyr::if_else(data[[results[[i, x]]]] >= var.values[[x]], T, F)
      } else {
        dplyr::if_else(data[[results[[i, x]]]] <= var.values[[x]], T, F)
      }
    })
    variable.high = rowSums(variable.high, na.rm = T) == ncol(variable.high)

    # identify observations where the outcome occurred
    outcome.occur = dplyr::if_else(research.plan$outcome.occurrence[[results$.outcome[i]]] == T, T, F)

    # set "near to end"
    if(!is.null(time)) {
      # identify outcome periods

      # create data with time and outcome stuff
      data.t = dplyr::select_at(data, c(unit, time))
      data.t$.outcome = outcome.occur
      data.t$.period = cumsum(data.t$.outcome)

      # first create time to outcome
      outcome.near = unlist(dplyr::group_map(dplyr::group_by_at(data.t, c(unit, ".period")), ~ {
        max(.x[[time]], na.rm = T) - .x[[time]]
      }))

      # set in data.t so we can see
      data.t$.near = outcome.near

      # identify observations that are near (or far if the effect is negative) to the outcome
      if(results$direction[i] == "negative") {
        outcome.near = dplyr::if_else(outcome.near >= quantile(outcome.near, 1 - .time.distance, na.rm = T), T, F)
      } else {
        outcome.near = dplyr::if_else(outcome.near <= quantile(outcome.near, .time.distance, na.rm = T), T, F)
      }
    } else {
      # identify outcomes (or not outcomes if effect is negative)
      if(results$direction[i] == "negative") {
        outcome.near = !outcome.near
      } else {
        outcome.near = outcome.occur
      }
    }

    # values to pull
    obs.to.keep = variable.high & outcome.near

    # not all observations are always present so drop nas
    obs.to.keep = obs.to.keep & !is.na(obs.to.keep)

    # do we have any
    if(all(obs.to.keep == F)) {
      return(NULL)
    }

    # find plausible observations
    support = data[obs.to.keep, ]

    # select columns
    support = dplyr::select(support, tidyselect::any_of(c(unit, time)))

    # set the variable name
    variables.name = plan$variables$label$var.label[match(results[i, .variables.t], plan$variables$label$var.name)]

    # combine variable into a name
    variable.label = paste0(variables.name, " (", results[i, paste0(.variables.t, ".status")], ")")

    # set additional values
    support$variable.checked = paste(variable.label, collapse = " + ")
    support$outcome = paste0(results$.outcome[i], " (", results$direction[i], ")")
    support$time = support[[time]]

    # return
    return(support)
  })

  # bind
  r = dplyr::bind_rows(r)

  # now summarize by our grouping vars
  r = dplyr::summarize(dplyr::group_by_at(r, c(unit, "outcome")),
                       time = paste(unique(time), collapse = ", "),
                       variable.checked = paste(unique(variable.checked), collapse = ", "),
                       num.obs = dplyr::n(),
                       .groups = "keep")

  # return
  return(r)
}

#' Function to summarize a dataframe and produce confidence intervals.
#'
#' This function summarizes predictions to create confidence intervals. This works for both the results of
#' a bootstrap and the posterior distribution of a Bayesian analysis.
#'
#' @family results
#' @param data Dataframe with variables to be summarized.
#' @param .ci Confidence interval. Defaults to 95 percent, which is the 2.5th to the 97.5th percentile of the distribution.
#' @param .grouping Variables that are excluded from the group and summarize. The variable "prob" should contain the estimates.
#' @param .round The number of digits to round the results to. Can be set to 'NULL' to prevent rounding.
#' @return A dataframe that contains the variable name, the median value, the lower confidence interval, the upper confidence
#'   interval, the simulated P-value, and the number of "draws" used to produce the confidence interval.
#' @keywords prediction confidence_interval ci summarize
#' @export
#' @examples
#' summarize_interval(results)
#'
summarize_interval = function(data, .ci = 0.95, .grouping = c(".prediction.id", "prob"), .round = 3) {
  # summarize interval -- could use HDI -- should be equivalent to http://mc-stan.org/rstanarm/reference/posterior_interval.stanreg.html

  # set return
  r = data

  # if a variable is a list than take the median
  r = dplyr::mutate_if(r, is.list, function(x) sapply(x, median))

  # identify variables to group by
  vars.group = colnames(r)[!colnames(r) %in% .grouping]

  # group -- suppresses implicit NA warning
  suppressWarnings(r <- dplyr::group_by_at(r, vars.group))

  # drop na probabilities
  r = dplyr::filter(r, !is.na(prob))

  # summarize
  r = dplyr::summarize(
    r,
    q = NA, #list(hdi(prob, .ci)),
    c = median(prob),
    c.low = quantile(prob, (1 - .ci) / 2), # q[[1]][1], # quantile(prob, (1 - .ci) / 2),
    c.high = quantile(prob, 1 - ((1 - .ci) / 2)), # q[[1]][2], # quantile(prob, 1 - ((1 - .ci) / 2)),
    p.value = min(dplyr::if_else(c < 0, 1 - ecdf(prob)(0), ecdf(prob)(0)) * 2, 1),
    draws = length(prob),
    .groups = "keep"
  )

  # drop q
  r = dplyr::select(r, -q)

  # ungroup
  r = dplyr::ungroup(r)

  # round
  if(!is.null(.round)) {
    r = dplyr::mutate_if(r, is.numeric, round, digits = .round)
  }

  # return
  return(r)
}


#' Structure predictions to make testing easier.
#'
#' This function structures data for testing a prediction. Deals with predictions that are a matrix for mediation analysis. Also
#' properly format factor variables.
#'
#' @family internal functions
#' @param formula The formula used to identify all variables that are in the model that will be run.
#' @param data The data used to create predictions.
#' @param predictions The predictions to test. This should be the object returned from 'pr_list'.
#' @param method A string indicating how variables not in 'predictions' should be treated. There are two options; "observed values" and "mean."
#' @return A dataframe
#'
structure_predictions = function(formula, data, predictions, method) {
  # internal function to make it easy to run one predict with observed values with factors and nicely setup contrasts

  ## GET DATA TO USE

  # identify variable
  vars.to.include = all.vars(formula)

  # subset the data
  if(any(vars.to.include == ".")) {
    data.all = dplyr::ungroup(data)
  } else {
    data.all = dplyr::select(dplyr::ungroup(data), all.vars(formula))
  }

  # select complete cases
  data.all = dplyr::filter(data.all, complete.cases(data.all))


  ## CREATE THE PREDICTION FRAME

  # if we have predictions make the frame and the data
  if(!is.null(predictions)) {
    # get variable classes
    classes = sapply(data.all, class)

    # fix factors in predictions frame -- loop allows predictions to be updated
    for(x in 1:length(predictions)) {
      for(x2 in 1:length(predictions[[x]])) {
        for(x3 in colnames(predictions[[x]][[x2]])) {
          if(any(x3 == names(classes)) & classes[x3] == "factor") {
            predictions[[x]][[x2]][[x3]] = factor(predictions[[x]][[x2]][[x3]], levels = levels(data.all[[names(classes[x3])]]))
          }
        }
      }
    }

    # combine -- gets rid of "Outer names are only allowed for unnamed scalar atomic inputs" warning
    suppressWarnings(prediction.frame <- dplyr::bind_rows(predictions))

    # remove ".name" column if present
    if(rlang::has_name(prediction.frame, ".name")) {
      prediction.frame = dplyr::select(prediction.frame, -.name)
    }

    # remove ".model" column if present
    if(rlang::has_name(prediction.frame, ".model")) {
      prediction.frame = dplyr::select(prediction.frame, -.model)
    }

    # make distinct (after name is removed) -- throws a warning for type 'list' which is what we use to pass a sampling distribution -- fix eventually but okay for now
    suppressWarnings(prediction.frame <- dplyr::distinct(prediction.frame))

    # use mean values or observed values
    if(dplyr::first(method) == "observed values") {
      data.prediction = data.all
    } else {
      # use mean for numeric values and mode otherwise
      data.prediction = dplyr::summarize_all(data.all, function(x) if (is.numeric(x)) mean(x, na.rm = T) else { r = unique(x); r[which.max(tabulate(match(x, r)))] })
    }

    # expand so that we have all observed values for each prediction
    data.prediction = data.prediction[rep(1:nrow(data.prediction), times = nrow(prediction.frame)), ]

    # check to make sure prediction vars are in data frame
    if(!all(colnames(prediction.frame) %in% colnames(data.prediction))) {
      stop("Prediction variables are missing from data frame.")
    }

    # add prediction id
    prediction.frame$.prediction.id = 1:nrow(prediction.frame)

    ## create the full data frame of observed values used for predictions

    # the prediction values to overlay on data.prediction
    prediction.temp = prediction.frame[rep(1:nrow(prediction.frame), each = nrow(data.prediction) / dplyr::n_distinct(prediction.frame$.prediction.id)), ]

    # set missing from master
    for(var in colnames(prediction.temp)) {
      prediction.temp[[var]][is.na(prediction.temp[[var]])] = data.prediction[[var]][is.na(prediction.temp[[var]])]
    }

    # overlay prediction values on main data frame

    # get colnames
    colnames.overlap = colnames(prediction.temp)

    # merge
    data.merge =
      lapply(colnames.overlap,
             function(x) {
               # identify which values need to be replaced
               add = as.vector(is.na(prediction.temp[, x]) | apply(prediction.temp[, x], 2, function(y) sapply(y, is.null)))

               # create return vector and replace needed values
               t = prediction.temp[[x]]
               t[add] = data.prediction[[x]][add]

               # return
               return(t)
             })
    data.merge = dplyr::as_tibble(data.merge, .name_repair = "minimal")
    colnames(data.merge) = colnames.overlap
    data.prediction[, colnames.overlap] = data.merge[, colnames.overlap]

    # this samples for each observed value -- we want to just take mean values and the inflate by the distribution size
    # if we passed a distribution instead of a discrete value, then go through and sample a random pair for each observation
    data.prediction = dplyr::mutate_if(data.prediction, is.list, function(x) sapply(x, function(y) if(!is.null(y)) as.list(y) else NULL))
    data.prediction = tidyr::unnest(data.prediction, cols = colnames(data.prediction)[sapply(data.prediction, is.list)], keep_empty = T)
    data.prediction = dplyr::mutate_if(data.prediction, is.list, function(x) sapply(x, function(y) { y[is.null(y)] = NA_real_; y }))

    # # transform to model matrix
    #
    # # save prediction id
    # .prediction.id = matrix(data = data.prediction$.prediction.id, ncol = 1, dimnames = list(NULL, ".prediction.id"))
    #
    # # transform
    # data.prediction = model.matrix(obj.formula, data.prediction)
    # data.prediction = cbind(data.prediction, .prediction.id)
    #
    # # turn into data frame
    # data.prediction = tibble::as_tibble(data.prediction)
  } else {
    prediction.frame = NULL
    data.prediction = NULL
  }


  ## CREATE THE CONTRAST DATA

  # if we have predictions create contrasts
  if(!is.null(predictions)) {
    # factor function
    create_blank_factor = function(t) {
      if(class(prediction.frame[[t]]) == "factor") {
        return(tibble::tibble(!!t := factor(NA_character_, levels = levels(prediction.frame[[t]]))))
      } else if(class(prediction.frame[[t]]) == "character" || class(prediction.frame[[t]]) == "list") {
        return(tibble::tibble(!!t := NA_character_))
      } else {
        return(tibble::tibble(!!t := NA_real_))
      }
    }

    # create a blank NA column that has correct factors
    all.na = dplyr::bind_cols(lapply(colnames(prediction.frame)[!colnames(prediction.frame) %in% c(".prediction.id", ".add")], create_blank_factor))

    # set contrast -- because we can have NAs in our prediction
    set_contrast = function(x) {
      # create a blank one row dataframe to stuff our values in
      cols = colnames(x)[colnames(x) %in% colnames(all.na)]
      r = all.na

      # stuff our values
      r[cols] = x[cols]

      # set reference
      r = dplyr::mutate_if(r, is.list, function(x) sapply(x, data.table::address))
      prediction.frame.temp = dplyr::mutate_if(prediction.frame, is.list, function(x) { y = sapply(x, data.table::address); y[sapply(x, is.null)] = NA_character_; y })

      # this mimics functionality in dplyr unique -- matches by memory reference when the dataframes have a list
      r = dplyr::left_join(r, prediction.frame.temp, by = colnames(r)[colnames(r) %in% colnames(prediction.frame.temp)])$.prediction.id

      # # identify the row in prediction.frame that matches r, which allows us to pull the proper .prediction.id
      # r = apply(sapply(colnames(r), function(x) sapply(1:nrow(prediction.frame), function(y) all(r[[x]][[1]] == prediction.frame[[x]][[y]]) || (is.na(r[[x]][[1]]) && is.na(prediction.frame[[x]][[y]])))), 1, all)
      # r[is.na(r)] = F
      #
      # # pull the prediction id
      # r = dplyr::first(prediction.frame$.prediction.id[r])

      # return
      return(r)
    }

    # set contrasts -- need to manually set the NAs in the predictions so we can left join correctly
    contrast.frame = lapply(predictions, function(y) sapply(y, set_contrast))
  } else {
    contrast.frame = NULL
  }


  ## DONE NOW RETURN

  # return the values
  return(list(predictions = prediction.frame, contrasts = contrast.frame, data = data.prediction))
}

#' Get draws from a frequentist or Bayesian object.
#'
#' This function selects the draws dataframe from a model object produced by 'analyze'. This function also
#' works on model objects produced by 'rstanarm', e.g., from a run of 'stan_glm'.
#'
#' @family internal functions
#' @param object The model object.
#' @return The dataframe containing all of the draws (model runs).
#'
get_draws = function(object) {
  # internal function to get draws (model results) for each run

  # get the values for each coefficient draw
  if(is.list(object)) {
    if(rlang::has_name(object, "stanfit")) {
      if(is.matrix(object$stanfit)) {
        draws = object$stanfit
      } else {
        draws = as.matrix(object)
      }
    } else if(rlang::has_name(object, "fit")) {
      draws = as.matrix(object$fit)
    }

    # colnames
    col.names = remove_backticks(colnames(draws))

    # turn into tibble
    draws = tibble::as_tibble(draws)
    colnames(draws) = col.names
  } else {
    draws = NULL
  }

  # return
  return(draws)
}

#' Get formula from a frequentist or Bayesian object.
#'
#' This function selects the formula from a model object produced by 'analyze'. This function also
#' works on model objects produced by 'rstanarm', e.g., from a run of 'stan_glm'.
#'
#' @family internal functions
#' @param object The model object.
#' @return The formula used in the model.
#'
get_formula = function(object) {
  # internal function to get formula from object

  # check
  if(is.null(object) || typeof(object) != "list" || !rlang::has_name(object, "formula")) {
    return(NULL)
  }

  # object is good get the formula

  # get formula
  if(purrr::is_formula(object$formula) || class(object$formula) %in% "brmsformula") {
    obj.formula = terms(as.formula(object$formula)) # good formula
  } else if(purrr::is_formula(object$formula$formula)) {
    obj.formula = terms(object$formula$formula) # nested formula
  } else if(is.list(object$formula)) {
    obj.formula  = lapply(object$formula, terms) # multiple formula
  } else {
    obj.formula = NULL
  }

  # return
  return(obj.formula)
}

#' Get data from a frequentist or Bayesian object.
#'
#' This function selects the data from a model object produced by 'analyze'. This function also
#' works on model objects produced by 'rstanarm', e.g., from a run of 'stan_glm'.
#'
#' @family internal functions
#' @param object The model object.
#' @return The data used in the model.
#'
get_data = function(object, data = NULL) {
  # internal function to get data from object

  # check
  if(is.null(object) || typeof(object) != "list" || !rlang::has_name(object, "data")) {
    return(data)
  }

  # return
  return(object$data)
}

#' Produce coefficients from a frequentist or Bayesian object.
#'
#' This function takes the draws result from a model and produces a coefficient table with confidence
#' intervals and P-values.
#'
#' @family internal functions
#' @param obj.formula The formula used in the model. A return from 'get_formula'.
#' @param obj.draws The draws produced from the model. A return from 'get_draws'.
#' @param obj.data The data used in the model. A return from 'get_data'.
#' @param times An optional vector identifying the times to produce coefficients for. Only relevant for
#'   survival models.
#' @return The coefficient table.
#'
results_coefficients = function(obj.formula, obj.draws, obj.data, times = NULL) {
  # internal function to process and return coefficients

  # TODO: make it get the time-varying coefficients as well

  # delete response to make it easier to work with -- especially with more complicated brms formulas
  obj.formula = formula(delete.response(obj.formula))

  # remove bars if needed
  if(lme4:::anyBars(obj.formula)) {
    # strip bars
    obj.formula = formula(lme4:::nobars(obj.formula))
  }

  # break into coefficients, predictions, contrasts functions

  # fix backticks in obj.draws
  colnames(obj.draws) = remove_backticks(colnames(obj.draws))

  # formatted formula without backticks
  formula.labels = remove_backticks(colnames(model.matrix(obj.formula, obj.data)))
  formula.labels = formula.labels[formula.labels %in% colnames(obj.draws)]

  # get coefficients
  coefs = tidyr::pivot_longer(dplyr::select(obj.draws, all_of(formula.labels)), everything(), names_to = "coefficient", values_to = "prob")

  # return
  return(coefs)
}

#' Produce predictions from a frequentist or Bayesian object.
#'
#' This function takes the draws result from a model and produces a prediction table with confidence
#' intervals and P-values based on structured predictions.
#'
#' @family internal functions
#' @param object The model object.
#' @param pr The structured predictions. A return from 'structure_predictions'.
#' @param obj.draws The draws produced from the model. A return from 'get_draws'.
#' @param method A string indicating how variables not in 'predictions' should be treated. There are two options; "observed values" and "mean."
#' @param times An optional vector identifying the times to produce coefficients for. Only relevant for
#'   survival models.
#' @param m The model to get results for. Only relevated for a multivariate Bayesian model (a model with
#'   multiple dependent variables).
#' @return The predictions table.
#'
results_predictions = function(object, pr, obj.draws, method = c("observed values", "mean"), draws = 1000, times = NULL, m = NULL) {
  # internal function to process and return predictions

  # do we have predictions
  if(is.null(pr)) {
    return(NULL)
  }

  # do we want to sample from the draws?
  if(!is.null(draws) && draws <= nrow(obj.draws)) {
    obj.draws = obj.draws[sample.int(nrow(obj.draws), size = draws), ]
  }

  # set m to null if no list of formula
  if(!is.list(object$formula)) m = NULL

  # type of predict
  if(any(class(object) %in% c("stansurv", "survival"))) {
    # for survival using rstanarm: for each time period (.pp_data_surv -> .pp_predict_surv [evaluate_log_surv]) send this package to (.pp_summarise_surv) = exp(lp) * bh

    # set times
    if(is.null(times) || !is.numeric(times)) {
      times = seq(1:as.integer(max(object$eventtime)))
    }

    # # do we need to supply our own baseline hazard
    # overwrite.baseline = rlang::has_name(object$basehaz, "raw.data")
    #
    # # get baseline hazard for times identified
    # if(overwrite.baseline) {
    #   # interpolate hazard
    #   bh.time = approx(x = object$basehaz$raw.data$time, y = object$basehaz$raw.data$hazard, xout = times)
    #
    #   # limit time to only what is available in the data
    #   bh.time$y[is.na(bh.time$y)] = max(bh.time$y, na.rm = T)
    #
    #   # format more nicely
    #   bh.time = tibble::tibble(time = bh.time$x, haz = bh.time$y)
    #
    #   # save original call
    #   original.call = rstanarm:::evaluate_log_basesurv
    #
    #   # hook into rstanarm:::evaluate_log_basesurv to pass our own baseline hazard when needed
    #   assignInNamespace("evaluate_log_basesurv", function(times, basehaz, aux, intercept = NULL) {
    #     # do we want our own baseline hazard or pass back to rstanarm
    #     if(overwrite.baseline) {
    #       return(-bh.time$haz[times])
    #     } else {
    #       return(original.call(times = times, basehaz = basehaz, aux = aux, intercept = intercept))
    #     }
    #   }, ns = "rstanarm")
    # }

    # # function to get the curve
    # predict_survival = function(time, object, pr, obj.draws, type = "cdf") {
    #   if(rlang::has_name(object$basehaz, "raw.data")) {
    #     # generate linear predictor
    #     pp.args = as.matrix(obj.draws) %*% t(pr$data[, colnames(obj.draws)])
    #
    #     # do exp(-exp(lp) * basehaz)
    #     # pp.args = exp(-bh.time$haz[bh.time$time == time]) ^ exp(pp.args)
    #     # pp.args = exp(-exp(pp.args) * bh.time$haz[bh.time$time == time])
    #     pp.args = exp(pp.args) * bh.time$haz[bh.time$time == time]
    #
    #     # get the cdf
    #     pp.args = 1 - exp(-pp.args)
    #
    #     # format the return
    #     predict.df = tibble::as_tibble(t(pp.args), .name_repair = "minimal")
    #     colnames(predict.df) = 1:ncol(predict.df)
    #
    #     # add prediction id
    #     predict.df$.prediction.id = pr$data$.prediction.id
    #     predict.df$.time = time
    #
    #     # return
    #     return(predict.df)
    #   }
    # }
    #
    # # combine
    # predict.df = dplyr::bind_rows(lapply(times, predict_survival, object = object, pr = pr, obj.draws = obj.draws))

    # # add intercept if needed
    # if(!rlang::has_name(pr$data, "(Intercept)") & !rlang::has_name(obj.draws, "(Intercept)")) {
    #   pr$data[, "(Intercept)"] = 0
    #   obj.draws[, "(Intercept)"] = 0
    # }

    # predict survival function -- default to cdf which is 1 - surv
    predict_survival = function(time, object, newdata, obj.draws, type = "cdf") {
      # to get the CDF it is basehaz * exp(beta)

      # get parameters
      pp.pars = rstanarm:::extract_pars.stansurv(object = object, stanmat = obj.draws, means = F)
      pp.pars = lapply(pp.pars, as.matrix) # saves as a dataframe so convert to a matrix

      # format the data
      pp.data = rstanarm:::.pp_data_surv(object = object, newdata = newdata, times = rep(time, nrow(newdata)), at_quadpoints = T)

      # get the prediction -- rows are draws, columns are newdata observations
      pp.args = rstanarm:::.pp_predict_surv.stansurv(object = object, data = pp.data, pars = pp.pars, type = type)

      # format the return
      predict.df = tibble::as_tibble(t(pp.args), .name_repair = "minimal")
      colnames(predict.df) = 1:ncol(predict.df)

      # add prediction id
      predict.df$.prediction.id = newdata$.prediction.id
      predict.df$.time = time

      # return
      return(predict.df)
    }

    # combine
    predict.df = dplyr::bind_rows(lapply(times, predict_survival, object = object, newdata = pr$data, obj.draws = obj.draws))

    # # return
    # if(overwrite.baseline) {
    #   # return the original function to the namespace
    #   assignInNamespace("evaluate_log_basesurv", original.call, ns = "rstanarm")
    # }
  } else {
    # using rstanarm built-in function -- each row is a draw, each column is a row in the prediction frame -- pp_eta gets the linear predictor, pp_args gets the inverse link of the LP

    # TODO: need to update to make it work with random effects -- requires setting is.mer

    # first get formatted data -- basically model.matrix -- average over random effects
    pp.data = rstanarm:::pp_data(object, pr$data, re.form = NA, offset = rep(0, nrow(pr$data)), m = m)

    # summarize if desired
    # if(dplyr::first(method) == "observed values") {
      # save prediction id
      .prediction.id = pr$data$.prediction.id
    # } else {
    #   # add in prediction id
    #   pp.data$x = tibble::as_tibble(pp.data$x)
    #   pp.data$x$.prediction.id = pr$data$.prediction.id
    #
    #   # summarize -- get means instead of observed values for predictions
    #   pp.data$x = dplyr::summarize_all(dplyr::group_by(pp.data$x, .prediction.id), mean)
    #
    #   # save prediction id
    #   .prediction.id = pp.data$x$.prediction.id
    #
    #   # save data back
    #   pp.data$x = as.matrix(dplyr::select(pp.data$x, -.prediction.id))
    # }

      # sometimes the object data (pp.data) can have more columns than actually ran so fix that -- occurs with fixed effects in frequentist models
      if("frequentist" %in% class(object)) pp.data$x = as.matrix(pp.data$x[, colnames(pp.data$x)[colnames(pp.data$x) %in% colnames(obj.draws)]])

      # obj.draws can also have missing values for the frequentist model -- if that is the case just make it zero
      obj.draws[is.na(obj.draws)] = 0


    # get eta -- the built-in function needs stanmat to be an S4 object, which is a pain so just send it manually
    # if(any(object$algorithm == "bootstrap")) {
      pp.eta = rstanarm:::pp_eta(object = object, data = pp.data, stanmat = as.matrix(obj.draws), m = m) # would like to fix this if possible
    # } else {
    #   pp.eta = rstanarm:::pp_eta(object, pp.data)
    # }

    # get the linear predictors and then transform by the inverse link function -- all pretty basic
    pp.args = rstanarm:::pp_args(object, pp.eta, m = m)

    # turn into usable data frame
    predict.df = tibble::as_tibble(t(pp.args$mu), .name_repair = "none")
    names(predict.df) = as.character(1:ncol(predict.df))
    predict.df$.prediction.id = .prediction.id
  }

  # group vars
  group.vars = c(".prediction.id", ".time")[c(".prediction.id", ".time") %in% colnames(predict.df)]

  # get mean observed values -- this basically synchronizes mean and observed value approaches
  if(dplyr::first(method) == "observed values") {
    predict.df = dplyr::summarize_all(dplyr::group_by_at(predict.df, group.vars), mean)
  }

  # pivot longer
  predict.df = tidyr::pivot_longer(predict.df, -tidyselect::all_of(group.vars), names_to = ".extra.name", values_to = "prob")

  # nicer version of predictions to add -- make a sent distribution look nicer
  nice.predictions = dplyr::mutate_if(pr$predictions, is.list, function(x) sapply(x, function(y) if(is.null(y)) NA else median(y)))

  # add variable names
  predict.df = dplyr::left_join(dplyr::select(dplyr::ungroup(predict.df), -.extra.name), nice.predictions, by = ".prediction.id")

  # return
  return(predict.df)
}

#' Produce contrasts from a frequentist or Bayesian object.
#'
#' This function takes the draws result from a model and produces a contrast table with confidence
#' intervals and P-values based on structured predictions. Contrasts are comparisons between
#' predictions.
#'
#' Supports a two-way and a four-way comparison. A two-way comparison compares two different predictions
#' while a four-way comparison compares the difference between two two-way comparisons. For instance, identifying
#' a treatment effect is a two-way comparison (high vs. low value). Comparing the treatment effect in two different
#' situations is a four way comparison (high vs. low in situation one compared to high vs. low in situation two).
#'
#'
#' @family internal functions
#' @param pr The structured predictions. A return from 'structure_predictions'.
#' @param predict.df The predictions produced from the model. A return from 'results_predictions'.
#' @return The contrasts table.
#'
results_contrasts = function(pr, predict.df) {
  # internal function to process and return contrasts

  # check if we have contrasts
  if(is.null(pr) | is.null(predict.df)) {
    return(NULL)
  }

  # function to get contrast
  get_contrasts = function(x, predict.df) {
    # length
    l = length(x)

    # get the contrast
    if(l == 2) { # 1 - 2
      r = predict.df$prob[predict.df$.prediction.id == x[1]] - predict.df$prob[predict.df$.prediction.id == x[2]]
    } else if(l == 4) { # (1 - 2) - (3 - 4)
      r = (predict.df$prob[predict.df$.prediction.id == x[1]] - predict.df$prob[predict.df$.prediction.id == x[2]]) -
        (predict.df$prob[predict.df$.prediction.id == x[3]] - predict.df$prob[predict.df$.prediction.id == x[4]])
    } else {
      r = NA
    }

    # could easily just make a new tibble that combines portions of predict.df -- and change the variable names to see what predictions where used (1, 2, 3, 4)

    # create tibble
    r = tibble::tibble(
      .contrast = dplyr::first(names(pr$contrasts)[sapply(pr$contrasts, function(c) all(c == x))]),
      .prediction.id = paste(x, collapse = ", "),
      prob = r)

    # if we have time add it
    if(rlang::has_name(predict.df, ".time")) {
      r$.time = predict.df$.time[predict.df$.prediction.id == x[1]]
    }

    # return
    return(r)
  }

  # get raw contrasts
  contrast.df = lapply(pr$contrasts, get_contrasts, predict.df = predict.df)
  contrast.df = dplyr::bind_rows(contrast.df)

  # return
  return(contrast.df)
}

#' Function to get results from an object returned by analysis.
#'
#' This function allows you to get results based on specified predictions from the return object produced by the
#' 'analyze' function. Results are produced on the scale of the response variable.
#'
#' Predictions estimate the value of the response variable when one or more independent variables is set to a desired
#' value (as indicated by 'predictions'). Contrasts compare two or more predictions.
#'
#' @family results
#' @param object Object returned from 'analysis' function.
#' @param predictions The predictions returned from "pr_list."
#' @param method Method for producing predictions Either uses observed values (the default) or mean values.
#' @param times Vector of times to produce predictions for. Only used for survival analysis.
#' @param .full.matrix Whether to return structured predictions (the default) or the full matrix.
#' @return A list with coefficient values, predictions, and contrasts. If '.full.matrix' is set, the unsummarized matrix
#'   for both predictions and contrasts is returned as well.
#' @keywords bootstrap results prediction
#' @export
#' @examples
#' results(object = output, predictions = main.predictions)
#'
results = function(object, predictions = NULL, method = c("observed values", "mean"), times = NULL, draws = 1000, .full.matrix = F) {
  ## get needed variables

  # get results from the model run
  obj.draws = get_draws(object)

  # get formula used
  obj.formula = get_formula(object)

  # get the data used for running the model
  obj.data = get_data(object)

  # are we dealing with variable distributions
  has.distributions = any(unlist(sapply(predictions, function(x) lapply(x[[1]], is.list))))

  # we cant do observed values and a variable with a distribution (memory problems) so just set to mean and go from there
  if(has.distributions) method = "mean"

  # need to deal with the possibility that formula/data are lists and that draws contains values for each list item -- turn the draws into a list too if that is the case
  if("list" %in% class(obj.formula)) {
    # check to make sure data is also a list
    if(!"list" %in% class(obj.data) | length(obj.data) != length(obj.formula)) {
      stop("Expecting both formula and data to be the same length.")
    }

    # both formula and data are lists of the same length now turn draws into a list too if it isnt already
    if(!"list" %in% class(obj.draws)) {
      obj.draws.t = sapply(names(obj.formula), function(nm) {
        # string to search for
        str.draw = paste0(nm, "\\|")

        # pull the values we need
        t.draws = obj.draws[, colnames(obj.draws)[stringr::str_detect(colnames(obj.draws), str.draw)]]

        # replace column names
        colnames(t.draws) = stringr::str_replace_all(colnames(t.draws), str.draw, "")

        # return
        return(t.draws)
      }, simplify = F)
    }

    # assemble
    obj.info = lapply(1:length(obj.formula), function(i) {
      list(obj.draws = obj.draws.t[[i]], obj.formula = obj.formula[[i]], obj.data = obj.data[[i]])
    })

    # set names
    names(obj.info) = names(obj.formula)
  } else {
    # create a list of one
    obj.info = list("y1" = list(obj.draws = obj.draws, obj.formula = obj.formula, obj.data = obj.data))
  }

  ## create coefficient frames

  # run through list
  r.coefs = lapply(obj.info, function(obj) {
    # get coefficients
    r.coefs = results_coefficients(obj.formula = obj$obj.formula, obj.draws = obj$obj.draws, obj.data = obj$obj.data)

    # summarize to get coefficient values
    r.coefs = summarize_interval(r.coefs)

    # return
    return(r.coefs)
  })

  # combine
  if(length(r.coefs) > 1) {
    r.coefs = dplyr::bind_rows(r.coefs, .id = ".model")
  } else {
    r.coefs = r.coefs[[1]]
  }


  ## create predictions and contrastsframe

  # get structured predictions
  if(!is.null(predictions)) {
    # structure predictions data frame
    pr = lapply(obj.info, function(x) structure_predictions(formula = x$obj.formula, data = x$obj.data, predictions = predictions, method = method))

    ## A PROBLEM IS THAT pr$data can have negative values for a variable that is being logged, etc. -- this breaks "make_model_frame"

    # run through the list for predictions - get predictions -- mvmer uses the full draws object so just pass that even though we made it nice above
    predict.df = lapply(1:length(obj.info), function(m) results_predictions(object = object, pr = pr[[m]], obj.draws = obj.draws, method = method, draws = draws, times = times, m = m))
    names(predict.df) = names(obj.info)

    # produce predictions that include the prediction id
    r.preds = lapply(predict.df, function(x) summarize_interval(x, .grouping = c("prob")))

    # get contrasts
    contrast.df = lapply(names(predict.df), function(x) results_contrasts(pr[[x]], predict.df[[x]]))
    names(contrast.df) = names(predict.df)

    # summarize
    r.contr = lapply(names(contrast.df), function(x) {
      # create contrast
      t.contr = summarize_interval(contrast.df[[x]], .grouping = c("prob"))

      # add prediction info
      pred.vals = lapply(t.contr$.prediction.id, function(pred.id) {
        # the id
        pred.id = as.integer(unlist(stringr::str_split(pred.id, ", ")))

        # collect prediction info
        r = lapply(1:length(pred.id), function(id) {
          # select
          r = dplyr::select(dplyr::filter(pr[[x]]$predictions, .prediction.id == pred.id[id]), -.prediction.id)
          r = r[, sapply(r, function(x) !all(is.na(x)))]

          # set column names
          colnames(r) = paste0("v", if(id %in% c(1, 2)) 1 else 2, ".", if(id %in% c(1, 3)) "high" else "low", ".p", 1:ncol(r))

          # make it a character so we can stack the return
          dplyr::mutate_all(r, as.character)
        })
      })

      # bind
      t.contr = dplyr::bind_cols(t.contr, dplyr::bind_rows(lapply(pred.vals, dplyr::bind_cols)))

      # return
      return(t.contr)
    })
    names(r.contr) = names(contrast.df)

    # combine all lists
    if(length(r.preds) > 1) {
      predict.df = dplyr::bind_rows(predict.df, .id = ".model")
      r.preds = dplyr::bind_rows(r.preds, .id = ".model")
      contrast.df = dplyr::bind_rows(contrast.df, .id = ".model")
      r.contr = dplyr::bind_rows(r.contr, .id = ".model")
    } else {
      predict.df = predict.df[[1]]
      r.preds = r.preds[[1]]
      contrast.df = contrast.df[[1]]
      r.contr = r.contr[[1]]
    }
  } else {
    # set everything to null
    predict.df = NULL
    r.preds = NULL
    contrast.df = NULL
    r.contr = NULL
  }


  ## return stuff

  # set return
  r = list(coefficients = r.coefs, predictions = r.preds, contrasts = r.contr)

  # add in full matrix if desired
  if(.full.matrix) {
    r[["predictions.matrix"]] = predict.df
    r[["contrasts.matrix"]] = contrast.df
  }

  # clean memory
  rm(r.coefs, r.preds, r.contr, predict.df, contrast.df)

  # return
  return(r)
}

#' Function to run mediation analysis.
#'
#' This function allows you to run mediation analysis based on two sets of returns from the 'analysis' function. The dependent variable
#' from 'm.mediator' must be present as an independent variable in 'm.outcome.' Transformations in the outcome model are allowed but transformations
#' in the mediator are not.
#'
#' All estimated effects are on the scale of the response variable.
#'
#' @family results
#' @param m.mediator Object returned from 'analysis' function that identifies the effect of the treatment on the mediator.
#' @param m.outcome Object returned from 'analysis' function that identifies the effect of the treatment and the mediator on the outcome.
#' @param predictions A 'pr_list' object with the desired change in the treatment. The treatment should be in both equations.
#' @param times An optional vector indicating the time(s) to conduct the mediation analysis if the outcome is a survival model.
#' @param draws An optional integer indicating the number of draws (randomly sampled) from the results to use. Smaller values allow a
#'   faster run.
#' @param .outcome Optional name for the outcome. If this is missing the name of the dependent variable in 'm.outcome' is used.
#' @return A datafame showing the direct effect, the indirect effect, the total effect, the portion of the total effect mediated,
#'   and the effect of the treatment on the mediator.
#' @keywords bootstrap results mediation
#' @export
#' @examples
#' results_mediation(m.mediator, m.outcome, predictions = predictions.mediation)
#'
results_mediation = function(m.mediator, m.outcome, predictions, times = NULL, draws = 1000, .outcome = NULL) {
  # mediation analysis using posterior distributions
  # the name of the mediator in m.med should be the same as in m.out + predictions should be in both -- add checks
  # needs a non-survival model upfront but can work with anything at the back
  # based on: https://imai.fas.harvard.edu/research/files/BaronKenny.pdf + https://imai.fas.harvard.edu/research/files/mediationP.pdf

  # use imai's updated baron-kenney approach to do mediation analysis

  # basic idea:
  #   model the mediator as a function of the treatment and the pre-treatment variables (causal effect of treatment on mediator)
  #   model the outcome as a function of the treatment, the mediator, and the pre-treatment/mediator variables (causal effect of treatment/mediator on outcome)
  #   identify change in outcome due to change in mediator caused by treatment

  # get the effect of predictions on the mediator -- should it just pass back the full matrix and the summaries
  effect.mediator = results(object = m.mediator, predictions = predictions, times = times, draws = draws, .full.matrix = T)

  # currently we are using point predictions for the effect of the treatment on the mediator
  # we should take into account uncertainty in this when carrying through our estimate of the indirect effect -- so if the treatment on the mediator is insignificant this needs to be accounted for
  # this is probably why our CIs are a little too low

  # create a prediction list that allows sampling (instead of a single value we pass it the distribution of values identified)
  mediation.var.name = all.vars(rlang::f_lhs(m.mediator$formula)) # allows us to have a transformation on the DV for the mediator -- still doesnt take the transformation into account -- should fix!!
  pr.mat = matrix(effect.mediator$predictions.matrix$prob, ncol = dplyr::n_distinct(effect.mediator$predictions.matrix$.prediction.id))
  pr.mediator = if(ncol(pr.mat) == 4) list(list(pr.mat[, 1] - pr.mat[, 2], pr.mat[, 3] - pr.mat[, 4])) else list(list(pr.mat[, 1], pr.mat[, 2])) # we either have a diff or a diff-in-diff
  names(pr.mediator) = mediation.var.name
  pr.mediator =  do.call(pr_list, pr.mediator)

  # # create predictions (using pr_list) for mediator based on the above results (the change in the mediator that results from the treatment predictions)
  # pr.mediator = list(effect.mediator$predictions$c)
  # names(pr.mediator) = as.character(as.formula(m.mediator$formula)[2]) # should probably have an option to set the mediator name or at least make sure it is identified consistently
  # pr.mediator =  do.call(pr_list, pr.mediator)

  # set outcome name if none provided
  if(is.null(.outcome)) .outcome = as.character(as.formula(m.outcome$formula)[2])

  # when getting results we could try sample from the treatment effect distribution instead of taking point estimates -- effect.mediator$predictions.matrix (the prediction ids for the contrast)


  # get raw data on mediation effect -- we need to bring in uncertainty about the effect of the treatment on mediator when producing predictions
  effect.outcome = results(object = m.outcome, predictions = c(predictions, pr.mediator), times = times, draws = draws, .full.matrix = T)

  # identify pairs -- this is the prediction for the treatment paired with the prediction for the mediator caused by the treatment
  # could modify to rely on the structure predictions function instead of hand coding here
  # pairs = lapply(1:nrow(effect.mediator$contrasts), function(x) {
  #   # get the predictions that correspond to this contrast
  #   r = as.numeric(unlist(stringr::str_split(effect.mediator$contrasts.matrix$.prediction.id[x], ", ")))
  #
  #   # make the pairs
  #   if(length(r) == 2) {
  #     r = paste0(effect.mediator$predictions$c[r[1]], " vs. ", effect.mediator$predictions$c[r[2]])
  #   } else {
  #     r = paste0(effect.mediator$predictions$c[r[1]], " vs. ", effect.mediator$predictions$c[r[2]], " <- ", effect.mediator$predictions$c[r[3]], " vs. ", effect.mediator$predictions$c[r[4]])
  #   }
  #
  #   # return
  #   return(c(effect.mediator$contrasts.matrix$.contrast[x], r))
  # })

  # create a vector of the correct length
  vec_length = function(v, length = length(v)) {
    r = rep_len(NA, length)
    r[1:length(v)] = v
    return(r)
  }

  # identify the mediation effect -- currently it is averaging over time -- should account for time -- the matrix from the result needs to return time
  r = lapply(length(names(predictions)), function(x) {
    # get the contrasts
    med.contrast = names(pr.mediator)[x]
    treat.contrast = names(predictions)[x]

    # get the effect of the treatment on the mediator
    mediator = effect.mediator$contrasts.matrix$prob[effect.mediator$contrasts.matrix$.contrast == treat.contrast]

    # get the effect of the treatment on the outcome independent of the mediator
    direct = effect.outcome$contrasts.matrix$prob[effect.outcome$contrasts.matrix$.contrast == treat.contrast]

    # get the effect of the mediator on the outcome due to the effect of the treatment
    indirect = effect.outcome$contrasts.matrix$prob[effect.outcome$contrasts.matrix$.contrast == med.contrast]

    # get the total effect of the treatment and the mediator
    total = direct + indirect

    # get the portion of the total effect explained by the mediator
    mediated = indirect / total

    # make sure everything is the correct length
    length = max(length(mediator), length(direct), length(indirect), length(total))

    # return
    return(tibble::tibble(.treatment = treat.contrast, .mediator = as.character(m.mediator$formula[2]), .mediator.contrast = med.contrast, .outcome = .outcome,
                          direct.effect = vec_length(direct, length), indirect.effect = vec_length(indirect, length), total.effect = vec_length(total, length),
                          prop.mediated = vec_length(mediated, length), treat.on.mediator = vec_length(mediator, length)))
  })

  # bind
  r = dplyr::bind_rows(r)

  # take the raw data from above and process -- kept in stages so what is done is transparent and easy to follow

  # set name
  r$.mediation = paste(r$.treatment, "->", r$.mediator.contrast)

  # pivot
  r = tidyr::pivot_longer(r, direct.effect:treat.on.mediator, names_to = ".effect", values_to = "prob")

  # select
  r = dplyr::select(r, .treatment, .effect, .mediator, .mediator.contrast, .outcome, prob)

  # summarize
  r = summarize_interval(r)

  # send back
  return(r)
}

#' Helper function to parse the 'contrasts' string in a result.
#'
#' Currently this function should not be used as it is built in to to the 'results' function.
#'
#' @export
#'
parse_contrast = function(df) {
  # parse effects -- not robust to factors/characters with " vs." or ", "
  r = tibble::as_tibble(t(sapply(df$.contrast, function(str) {
    r = unlist(stringr::str_split(str, " vs. "))
    if(stringr::str_detect(str, ", ")) {
      r = unlist(stringr::str_split(r, ", "))
      c(v1.high = r[1], v1.low = r[3], v2.high = r[2], v2.low = r[4])
    } else {
      c(v1.high = r[1], v1.low = r[2])
    }
  })))
  r
}

#' Provide a qualitative assessment for a set of results (main, interaction, and mediation effects).
#'
#' This function provides a human interpretable understanding of a research plan and the results from analyzing
#' the research plan.
#'
#' @family results
#' @param research.plan Object returned from 'research_plan' function.
#' @param all.results Object returned from 'analyze_plan' function.
#' @return A list providing an English-language description of the results for each outcome/direction. Currently, the
#'   a description of all of the variables that impact an outcome in a direction (and when and how), as well as a rank
#'   ordering of the magnitude of the effect of the variables for an outcome in a direction. Only results that are
#'   statistically significant at the 0.1 level or better are described.
#' @keywords results, qualitative
#' @export
#'
qualitative_assessment = function(research.plan, all.results) {
  # TODO: allow user to select criteria to scope the assessment (e.g., effect when source has a large economy and target has a large army, etc.)

  # set the parts
  main.res = purrr::map_dfr(all.results, "main.contrasts")
  interaction.res = purrr::map_dfr(all.results, "interaction.contrasts", .id = ".main.variable")
  mediation.res = purrr::map_dfr(all.results, "mediation.contrasts", .id = ".main.variable")

  # # duplicates in mediation.res because our plan had duplicates -- should get the plan to check duplicates
  # mediation.res = mediation.res[!duplicated(dplyr::select(mediation.res, .main.variable, .main.mediation, .outcome, .effect)), ]

  # how many significant values would you need to look through
  # table(c(main.res$p.value, interaction.res$p.value, mediation.res$p.value) < 0.1)

  # function to easily paste a list -- add below to variable_effect
  paste_list = function(items, wrap = NULL, .join = "and") {
    # wrap if needed
    if(is.character(wrap)) {
      items = paste0(wrap, items, wrap)
    }

    # drop NA
    items = unique(items[!is.na(items)])

    # return if empty
    if(length(items) < 1) {
      return("")
    }

    # commas or commas and
    if(length(items) < 3) {
      str = paste(items, collapse = paste0(" ", .join, " "))
    } else {
      str = paste(paste(items[-length(items)], collapse = ", "), items[length(items)], sep = paste0(", ", .join, " "))
    }

    # return
    return(str)
  }

  # provide qualitative indicators of the direction, significance, and size of an effect
  qualitative_effect = function(df) {
    # make sure it has the right columns
    if(!all(c("c", "c.low", "c.high", "p.value", ".compare") %in% colnames(df))) {
      stop("Dataframe does not have correct columns.")
    }

    # add blank time if missing
    if(!rlang::has_name(df, ".time")) df$.time = NA

    # set additional variable name
    .additional = if(rlang::has_name(df, ".main.interaction")) ".main.interaction" else ".main.mediation"

    # create the return
    r = tibble::tibble(
      # add labels
      .main.label = plan$variables$label$var.label[match(df$.main.variable, plan$variables$label$var.name)],
      .additional.label = plan$variables$label$var.label[match(df[[.additional]], plan$variables$label$var.name)],

      # set time
      time = factor(dplyr::case_when(
        is.na(df$.time) ~ "",
        df$.time < quantile(df$.time, 0.33, na.rm = T) ~ "short",
        df$.time < quantile(df$.time, 0.66, na.rm = T) ~ "medium",
        df$.time >= quantile(df$.time, 0.66, na.rm = T) ~ "long"
      ), levels = c("short", "medium", "long")),

      # how significant the effect is
      significance = dplyr::case_when(
        df$p.value < 0.01 ~ "very significant",
        df$p.value < 0.05 ~ "significant",
        df$p.value < 0.1 ~ "possible",
        T ~ "insignificant"
      ),

      # the direction of the effect
      direction = dplyr::case_when(
        df$c < 0 ~ "negative",
        T ~ "positive"
      ),

      # size of effect relative to comparison -- if compare is of zero size just default to the highest value
      percent = dplyr::if_else(df$.compare == 0, 999, abs(df$c) / abs(df$.compare)),

      # formatted
      percent.label = dplyr::case_when(
        percent > 10 ~ "(>1,000%)",
        percent < 0.1 ~ "(<10%)",
        T ~ paste0("(", scales::percent(percent, big.mark = ",", accuracy = 1), ")")
      ),

      # turn the percentage into a
      size = factor(dplyr::case_when(
        # percent > 3 ~ "huge",
        percent > 2.5 ~ "very large",
        percent > 1.5 ~ "large",
        percent < 1/1.5 ~ "small",
        percent < 1/2.5 ~ "very small",
        # percent < 1/3 ~ "tiny",
        T ~ "average"
      ), levels = c("very small", "small", "average", "large", "very large")),

      # formatted
      size.label = paste(size, percent.label)
    )


    # return
    return(r)
  }

  # add NA time column if not present
  if(!rlang::has_name(main.res, ".time")) main.res$.time = NA
  if(!rlang::has_name(interaction.res, ".time")) interaction.res$.time = NA

  # combine and select
  all.effects = dplyr::bind_rows(
    dplyr::bind_cols(dplyr::select(dplyr::ungroup(main.res), tidyselect::any_of(c(".outcome", ".time", ".main.variable", "c", "c.low", "c.high", "p.value"))), parse_contrast(main.res)),
    dplyr::bind_cols(dplyr::select(dplyr::ungroup(interaction.res), tidyselect::any_of(c(".outcome", ".time", ".main.variable", ".main.interaction", "c", "c.low", "c.high", "p.value"))), parse_contrast(interaction.res))
  )

  # set mediation too
  med.effects = mediation.res

  # add compare
  all.effects = dplyr::mutate(dplyr::group_by(all.effects, .outcome, .time), .compare = mean(abs(c[p.value < 0.1])))
  if(!is.null(med.effects) && nrow(med.effects) > 0) med.effects = dplyr::mutate(dplyr::group_by(med.effects, .outcome), .compare = mean(abs(c[.effect == "total.effect"])))

  # identify low and high values
  all.effects = dplyr::mutate(dplyr::group_by(all.effects, .outcome, .time, .main.variable, .main.interaction),
                              v1.size = dplyr::case_when(v1.high != v1.low ~ "both", v1.high == max(v1.high) ~ "high", T ~ "low"),
                              v2.size = dplyr::case_when(is.na(v2.high) ~ NA_character_, v2.high != v2.low ~ "both", v2.high == max(v2.high) ~ "high", T ~ "low"))

  # describe effects
  all.effects = dplyr::bind_cols(all.effects, qualitative_effect(all.effects))
  if(!is.null(med.effects) && nrow(med.effects) > 0) med.effects = dplyr::bind_cols(med.effects, qualitative_effect(med.effects))

  # first filter out insignificant results
  all.signif = dplyr::filter(all.effects, significance != "insignificant")
  if(!is.null(med.effects) && nrow(med.effects) > 0) {
    med.signif = dplyr::select(dplyr::filter(
      dplyr::mutate(dplyr::group_by(med.effects, .outcome, .main.variable, .main.mediation),
                    med.direction = if(any(c[.effect == "treat.on.mediator"] * c[.effect == "indirect.effect"] >= 0)) "positive" else "negative",
                    .to.keep = any(significance[.effect %in% c("indirect.effect", "treat.on.mediator")] != "insignificant")),
      .to.keep == T), -.to.keep)
  } else {
    med.signif = NULL
  }

  # summarize to get unique effects
  all.signif.outcome = dplyr::summarize(
    dplyr::group_by(all.signif, .outcome, .main.variable, .main.label, .main.interaction, .additional.label, v2.size, significance, direction, size),
    time = paste_list(time),
    label = size.label[median(dplyr::n())],
    .groups = "keep"
  )

  # make sure we have significant outcomes
  if(nrow(all.signif.outcome) == 0) {
    return(list(all = NULL, ranked = NULL))
  }

  # two ways to summarize -- by variable and by outcome

  # summarize by outcome
  all.signif.outcome = dplyr::group_map(dplyr::group_by(all.signif.outcome, .outcome, direction), ~ {
    # determine if there are multiple variables
    .plural = if(dplyr::n_distinct(.x$.main.label) > 1) T else F

    # start string
    str = paste0("The outcome '", unique(.x$.outcome), "' is ", if(unique(.x$direction) == "positive") "increased" else "decreased", " by ",
                 if(.plural) "several variables." else "one variable.")

    ## ADD IN SIZE OF EFFECT, TIMING OF EFFECT, COMPARISON OF EFFECT SIZE, AND ABILITY TO PULL OUT AN EFFECT FOR A SPECIFIC CASE

    # create variable strings
    var.string = dplyr::group_map(dplyr::group_by(.x, .main.variable), ~ {
      # set unconditional
      .unconditional = if(any(is.na(.x$.main.interaction))) T else F

      # # and if we have a vowel before size
      # .vowel = if(grepl("^(a|e|i|o|u)", unique(.x$size))) T else F

      # create string
      str = paste0("Variable '", unique(.x$.main.label), "' has a ", unique(.x$direction), " ", if(.unconditional) "unconditional ",
            "effect on '", unique(.x$.outcome), "'")

      # set when (interaction)
      if(any(!is.na(.x$.main.interaction))) {
        when.str = dplyr::group_map(dplyr::group_by(dplyr::filter(.x, !is.na(.main.interaction)), v2.size), ~ {
          # determine if there are multiple variables
          .plural = if(dplyr::n_distinct(.x$.main.interaction) > 1) T else F

          # string
          paste0(paste_list(.x$.additional.label, wrap = "'", .join = "or"), if(.plural) " are '" else " is '", unique(.x$v2.size), "'")
        }, .keep = T)

        # combine
        when.str = paste0(if(.unconditional) paste0(" and a '", unique(.x$direction), "' conditional effect"), " when ", paste_list(when.str, .join = "and/or"))
      } else {
        when.str = ""
      }

      # now deal with mediation when it is present
      if(!is.null(med.signif)) {
        # pull mediation data -- only for a certain direction so also filter that
        .x.med = dplyr::filter(med.signif, .main.variable == unique(.x$.main.variable), .outcome == unique(.x$.outcome), med.direction == unique(.x$direction))
        .x.med = dplyr::mutate(dplyr::group_by(.x.med, .main.mediation),
                               .treat.direction = dplyr::if_else(c[.effect == "treat.on.mediator"] >= 0, "increasing", "decreasing"))

        # describe mediation effect if present
        if(nrow(.x.med) > 0) {
          # get the individual effects
          med.dir.effect = dplyr::group_map(dplyr::group_by(.x.med, .treat.direction), ~ {
            paste0(unique(.x.med$.treat.direction), " ", paste_list(.x.med$.additional.label, wrap = "'"))
          })

          med.str = paste0(". Variable '", unique(.x$.main.label), "' also has a mediating impact by ", paste_list(med.dir.effect))
        } else {
          med.str = ""
        }
      } else {
        med.str = ""
      }

      # combine
      str = paste0(str, when.str, med.str, ".")

      # return
      return(str)
    }, .keep = T)

    # combine
    str = c(unlist(str), unlist(var.string))

    # return
    return(str)
  }, .keep = T)

  # expand add in the comparison of results to identify what is best for a particular outcome -- should allow general guidance and also guidance under specific circumstances -- all effects or only significant?

  # if we look at conditionality we get information overload; if we don't we can get misleading estimates; really this will work best when a user can select some conditions under which they want to estimate effects

  # rank order effects
  all.signif.ranked = dplyr::group_map(dplyr::group_by(dplyr::filter(all.signif, significance != "possible"), .outcome, direction), ~ {
    # arrange by largest effect first
    signif.arrange = dplyr::arrange(.x, dplyr::desc(percent))

    # reverse factor ordering just so we have highest first
    signif.arrange$size = forcats::fct_rev(signif.arrange$size)

    # and group by time and size
    str = dplyr::group_map(dplyr::group_by(signif.arrange, time, size), ~ {
      # and if we have a vowel before size
      .vowel = if(grepl("^(a|e|i|o|u)", unique(.x$size))) T else F

      # add condition
      .x$.full.label = dplyr::if_else(is.na(.x$.additional.label), paste0(.x$.main.label, " (unconditional)"), paste0(.x$.main.label, " (conditional)"))
      # .x$.full.label = dplyr::if_else(is.na(.x$.additional.label), .x$.main.label, paste0(.x$.main.label, " (when ", .x$.additional.label, " is '", .x$v2.size, "')"))

      # determine if there are multiple variables
      .plural = if(dplyr::n_distinct(.x$.full.label) > 1) T else F

      # create string
      paste0(paste_list(unique(.x$.full.label)), " ", if(.plural) "have" else "has", " a", if(.vowel) "n " else " ", unique(.x$size), " impact")
    }, .keep = T)

    # create string -- also deal with no time being present
    paste0(if(is.na(unique(.x$time))) "In general " else paste0("In the ", unique(.x$time), " run "), paste_list(str))
  }, .keep = T)


  # set the names
  names(all.signif.outcome) = apply(unique(dplyr::select(dplyr::ungroup(all.signif), .outcome, direction)), 1, paste, collapse = ", ")
  names(all.signif.ranked) = apply(unique(dplyr::select(dplyr::ungroup(all.signif), .outcome, direction)), 1, paste, collapse = ", ")

  # return -- has main effects, interaction effects, and mediation effects sorted by outcome
  return(list(all = all.signif.outcome, ranked = all.signif.ranked))
}

