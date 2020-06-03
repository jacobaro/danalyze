# functions to get results from analysis

# summary stats

#' Function to get summary statistics from a dataframe.
#'
#' This function produces a table of summary statistics based on a formula and a dataframe.
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


# summarize interval -- could use HDI -- should be equivalent to http://mc-stan.org/rstanarm/reference/posterior_interval.stanreg.html

#' Function to summarize a dataframe and produce confidence intervals.
#'
#' This function summarizes predictions to create confidence intervals.
#' @param data Dataframe with variables to be summarized.
#' @param .ci Confidence interval. Defaults to 95 percent, which is the 2.5th to the 97.5th percentile of the distribution.
#' @param .grouping Variables that are excluded from the group and summarize. The variable "prob" should contain the estimates.
#' @param .round The number of digits to round the results to. Can be set to 'NULL' to prevent rounding.
#' @keywords prediction confidence_interval ci summarize
#' @export
#' @examples
#' summarize_interval(results)
#'

summarize_interval = function(data, .ci = 0.95, .grouping = c(".prediction.id", "prob"), .round = 3) {
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
    draws = length(prob)
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

# internal function
# structure predictions to make it easy to run one predict with observed values with factors and nicely setup contrasts
structure_predictions = function(formula, data, predictions, method) {

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

    # combine -- move to this instead: vctrs::vec_c()
    prediction.frame = dplyr::bind_rows(dplyr::combine(predictions))

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

# internal function to get draws (model results) for each run
get_draws = function(object) {
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
  } else {
    draws = NULL
  }

  # colnames
  col.names = remove_backticks(colnames(draws))

  # turn into tibble
  draws = tibble::as_tibble(draws)
  colnames(draws) = col.names

  # return
  return(draws)
}

# internal function to get formula from object
get_formula = function(object) {
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

# internal function to get data from object
get_data = function(object, data = NULL) {
  # check
  if(is.null(object) || typeof(object) != "list" || !rlang::has_name(object, "data")) {
    return(data)
  }

  # return
  return(object$data)
}

# internal function to process and return coefficients
# TODO: make it get the time-varying coefficients as well
results_coefficients = function(obj.formula, obj.draws, obj.data, times = NULL) {
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

# internal function to process and return predictions
results_predictions = function(object, pr, obj.draws, method = c("observed values", "mean"), draws = 1000, times = NULL) {
  # do we have predictions
  if(is.null(pr)) {
    return(NULL)
  }

  # do we want to sample from the draws?
  if(!is.null(draws) && draws <= nrow(obj.draws)) {
    obj.draws = obj.draws[sample.int(nrow(obj.draws), size = draws), ]
  }

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
    pp.data = rstanarm:::pp_data(object, pr$data, re.form = NA, offset = rep(0, nrow(pr$data)))

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

    # get eta -- the built-in function needs stanmat to be an S4 object, which is a pain so just send it manually
    # if(any(object$algorithm == "bootstrap")) {
      pp.eta = rstanarm:::pp_eta(object, pp.data, stanmat = as.matrix(obj.draws)) # would like to fix this if possible
    # } else {
    #   pp.eta = rstanarm:::pp_eta(object, pp.data)
    # }

    # get the linear predictors and then transform by the inverse link function -- all pretty basic
    pp.args = rstanarm:::pp_args(object, pp.eta)

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

# internal function to process and return contrasts
results_contrasts = function(pr, predict.df) {
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

# function to get predictions from an "out" object produced by analyze

#' Function to get results from an object returned by analysis.
#'
#' This function allows you to get results from an analysis.
#' @param object Object returned from 'analysis' function.
#' @param predictions The predictions returned from "pr_list."
#' @param method Method for producing predictions Either uses observed values (the default) or mean values.
#' @param times Vector of times to produce predictions for. Only used for survival analysis.
#' @param .full.matrix Whether to return structured predictions (the default) or the full matrix.
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


  ## create coefficient frames

  # get coefficients
  r.coefs = results_coefficients(obj.formula, obj.draws, obj.data)

  # summarize to get coefficient values
  r.coefs = summarize_interval(r.coefs)

  ## create predictions and contrastsframe

  # get structured predictions
  if(!is.null(predictions)) {
    # structure predictions data frame
    pr = structure_predictions(formula = obj.formula, data = obj.data, predictions = predictions, method = method)

    # get predictions
    predict.df = results_predictions(object = object, pr = pr, obj.draws = obj.draws, method = method, draws = draws, times = times)

    # produce predictions that include the prediction id
    r.preds = summarize_interval(predict.df, .grouping = c("prob"))

    # get contrasts
    contrast.df = results_contrasts(pr, predict.df)

    # summarize
    r.contr = summarize_interval(contrast.df, .grouping = c("prob"))
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


# mediation analysis using posterior distributions
# the name of the mediator in m.med should be the same as in m.out + predictions should be in both -- add checks
# needs a non-survival model upfront but can work with anything at the back
# based on: https://imai.fas.harvard.edu/research/files/BaronKenny.pdf + https://imai.fas.harvard.edu/research/files/mediationP.pdf

#' Function to run a mediation analysis.
#'
#' This function allows you to run mediation analysis based on two sets of returns from the 'analysis' function.
#' @param m.mediator Object returned from 'analysis' function that identifies the effect of the treatment on the mediator.
#' @param m.outcome Object returned from 'analysis' function that identifies the effect of the treatment and the mediator on the outcome.
#' @param predictions A 'pr_list' object with the desired change in the treatment. The treatment should be in both equations.
#' @param .outcome Optional name for the outcome. If this is missing the name of the variable is used automatically.
#' @keywords bootstrap results mediation
#' @export
#' @examples
#' results_mediation(m.mediator, m.outcome, predictions = predictions.mediation)
#'

results_mediation = function(m.mediator, m.outcome, predictions, times = NULL, .outcome = NULL) {
  # use imai's updated baron-kenney approach to do mediation analysis

  # basic idea:
  #   model the mediator as a function of the treatment and the pre-treatment variables (causal effect of treatment on mediator)
  #   model the outcome as a function of the treatment, the mediator, and the pre-treatment/mediator variables (causal effect of treatment/mediator on outcome)
  #   identify change in outcome due to change in mediator caused by treatment

  # get the effect of predictions on the mediator -- should it just pass back the full matrix and the summaries
  effect.mediator = results(object = m.mediator, predictions = predictions, times = times, .full.matrix = T)

  # currently we are using point predictions for the effect of the treatment on the mediator
  # we should take into account uncertainty in this when carrying through our estimate of the indirect effect -- so if the treatment on the mediator is insignificant this needs to be accounted for
  # this is probably why our CIs are a little too low

  # create a prediction list that allows sampling (instead of a single value we pass it the distribution of values identified)
  mediation.var.name = as.character(as.formula(m.mediator$formula)[2])
  pr.mat = matrix(effect.mediator$predictions.matrix$prob, ncol = dplyr::n_distinct(effect.mediator$predictions.matrix$.prediction.id))
  pr.mediator = list(list(pr.mat[, 1], pr.mat[, 2]))
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
  effect.outcome = results(object = m.outcome, predictions = c(predictions, pr.mediator), times = times, .full.matrix = T)

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
