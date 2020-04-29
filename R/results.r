# functions to get results from analysis

#' Function to summarize a dataframe and produce confidence intervals.
#'
#' This function summarizes predictions to create confidence intervals.
#' @param data Dataframe with variables to be summarized.
#' @param .ci Confidence interval. Defaults to 95% (2.5% to 97.5% of distribution).
#' @param .grouping Variables that are excluded from the group and summarize. The variable "prob" should contain the estimates.
#' @param .round The number of digits to round the results to. Can be set to 'NULL' to prevent rounding.
#' @keywords prediction, confidence interval, summarize
#' @export
#' @examples
#' summarize_interval(results)
#'


# summarize interval -- could use HDI -- should be equivalent to http://mc-stan.org/rstanarm/reference/posterior_interval.stanreg.html
summarize_interval = function(data, .ci = 0.95, .grouping = c(".prediction.id", "prob"), .round = 3) {
  # set return
  r = data

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
    r = dplyr::mutate_if(r, is.numeric, round, digits = round)
  }

  # return
  return(r)
}

# internal function
# structure predictions to make it easy to run one predict with observed values with factors and nicely setup contrasts
structure_predictions = function(formula, data, predictions) {

  ## GET DATA TO USE

  # subset the data
  data.all = dplyr::select(dplyr::ungroup(data), all.vars(formula))

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

    # combine
    prediction.frame = dplyr::bind_rows(dplyr::combine(predictions))

    # remove ".name" column if present
    if(rlang::has_name(prediction.frame, ".name")) {
      prediction.frame = dplyr::select(prediction.frame, -.name)
    }

    # remove ".model" column if present
    if(rlang::has_name(prediction.frame, ".model")) {
      prediction.frame = dplyr::select(prediction.frame, -.model)
    }

    # make distinct (after name is removed)
    prediction.frame = dplyr::distinct(prediction.frame)

    # expand so that we have all observed values for each prediction
    data.prediction = data.all[rep(1:nrow(data.all), times = nrow(prediction.frame)), ]

    # check to make sure prediction vars are in data frame
    if(!all(colnames(prediction.frame) %in% colnames(data.prediction))) {
      stop("Prediction variables are missing from data frame.")
    }

    # add prediction id
    prediction.frame$.prediction.id = 1:nrow(prediction.frame)

    ## create the full data frame of observed values used for predictions

    # the prediction values to overlay on data.prediction
    prediction.temp = prediction.frame[rep(1:nrow(prediction.frame), each = nrow(data.all)), ]

    # set missing from master
    for(var in colnames(prediction.temp)) {
      prediction.temp[[var]][is.na(prediction.temp[[var]])] = data.prediction[[var]][is.na(prediction.temp[[var]])]
    }

    # overlay prediction values on main data frame
    data.prediction[, colnames(prediction.temp)] = prediction.temp[, colnames(prediction.temp)]
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
      } else if(class(prediction.frame[[t]]) == "character") {
        return(tibble::tibble(!!t := NA_character_))
      } else {
        return(tibble::tibble(!!t := NA_real_))
      }
    }

    # create a blank NA column that has correct factors
    all.na = dplyr::bind_cols(lapply(colnames(prediction.frame)[!colnames(prediction.frame) %in% c(".prediction.id", ".add")], create_blank_factor))

    # set contrast -- because we can have NAs in our prediction
    set_contrast = function(x) {
      cols = colnames(x)[colnames(x) %in% colnames(all.na)]
      r = all.na
      r[cols] = x[cols]
      r = suppressWarnings(dplyr::left_join(r, prediction.frame)$.prediction.id)
      return(r)
    }

    # set contrasts -- need to manually set the NAs in the predictions so we can left join correctly
    contrast.frame = lapply(predictions, function(y) suppressMessages(sapply(y, set_contrast)))
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
  } else if(purrr::is_formula(object$formula$formul)) {
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
results_coefficients = function(obj.formula, obj.draws, obj.data) {
  # delete response to make it easier to work with -- especially with more complicated brms formulas
  obj.formula = formula(delete.response(obj.formula))

  # remove bars if needed
  if(lme4:::anyBars(obj.formula)) {
    # strip bars
    obj.formula = formula(lme4:::nobars(obj.formula))
  }

  # break into coefficients, predictions, contrasts functions

  # formatted formula without backticks
  formula.labels = colnames(model.matrix(obj.formula, obj.data))

  # get coefficients
  coefs = tidyr::pivot_longer(dplyr::select(obj.draws, all_of(formula.labels)), everything(), names_to = "coefficient", values_to = "prob")

  # return
  return(coefs)
}

# internal function to process and return predictions
results_predictions = function(object, pr, obj.draws, method = c("observed values", "mean"), times = NULL) {
  # do we have predictions
  if(is.null(pr)) {
    return(NULL)
  }

  # type of predict
  if(any(class(object) %in% c("stansurv", "survival"))) {
    # for survival using rstanarm: for each time period (.pp_data_surv -> .pp_predict_surv [evaluate_log_surv]) send this package to (.pp_summarise_surv) = exp(lp) * bh

    # set times
    if(is.null(times) || !is.numeric(times)) {
      times = seq(1:as.integer(max(object$eventtime)))
    }

    # predict survival function -- default to cdf which is 1 - surv
    predict_survival = function(time, object, pr, stanmat, type = "cdf") {
      # first format the data
      pp.data = rstanarm:::.pp_data_surv(object = object, newdata = pr$data, times = rep(time, nrow(pr$data)), at_quadpoints = T)

      # summarize if desired
      if(dplyr::first(method) == "observed values") {
        # save prediction id
        .prediction.id = pr$data$.prediction.id
      } else {
        # add in prediction id
        pp.data$x = tibble::as_tibble(pp.data$x)
        pp.data$x$.prediction.id = pr$data$.prediction.id

        # set other vars
        pp.data$pts = sapply(unique(pr$data$.prediction.id), function(x) mean(pp.data$pts[pp.data$x$.prediction.id == x], na.rm = T))
        pp.data$wts = sapply(unique(pr$data$.prediction.id), function(x) mean(pp.data$wts[pp.data$x$.prediction.id == x], na.rm = T))
        pp.data$ids = sapply(unique(pr$data$.prediction.id), function(x) mean(pp.data$ids[pp.data$x$.prediction.id == x], na.rm = T))
        pp.data$s = matrix(0, length(pp.data$pts), 0) #### NEES TO BE FIXED for time varying covariates

        # summarize -- get means instead of observed values for predictions
        pp.data$x = dplyr::summarize_all(dplyr::group_by(pp.data$x, .prediction.id), mean)

        # save prediction id
        .prediction.id = pp.data$x$.prediction.id

        # save data back
        pp.data$x = as.matrix(dplyr::select(pp.data$x, -.prediction.id))
      }

      # now get the full package for prediction
      pp.pars = rstanarm:::extract_pars.stansurv(object = object, stanmat = stanmat, means = F)
      pp.pars = lapply(pp.pars, as.matrix) # saves as a dataframe so convert to a matrix

      # get the prediction -- rows are draws, columns are newdata observations
      pp.args = rstanarm:::.pp_predict_surv.stansurv(object = object, data = pp.data, pars = pp.pars, type)

      # format the return
      predict.df = tibble::as_tibble(t(pp.args), .name_repair = "minimal")
      colnames(predict.df) = 1:ncol(predict.df)

      # add prediction id
      predict.df$.prediction.id = .prediction.id
      predict.df$.time = time

      # return
      return(predict.df)
    }

    # combine
    predict.df = dplyr::bind_rows(lapply(times, predict_survival, object = object, pr = pr, stanmat = obj.draws))
  } else {
    # using rstanarm built-in function -- each row is a draw, each column is a row in the prediction frame -- pp_eta gets the linear predictor, pp_args gets the inverse link of the LP

    # first get formatted data -- basically model.matrix -- average over random effects
    pp.data = rstanarm:::pp_data(object, pr$data, re.form = NA)

    # summarize if desired
    if(dplyr::first(method) == "observed values") {
      # save prediction id
      .prediction.id = pr$data$.prediction.id
    } else {
      # add in prediction id
      pp.data$x = tibble::as_tibble(pp.data$x)
      pp.data$x$.prediction.id = pr$data$.prediction.id

      # summarize -- get means instead of observed values for predictions
      pp.data$x = dplyr::summarize_all(dplyr::group_by(pp.data$x, .prediction.id), mean)

      # save prediction id
      .prediction.id = pp.data$x$.prediction.id

      # save data back
      pp.data$x = as.matrix(dplyr::select(pp.data$x, -.prediction.id))
    }

    # get eta -- the built-in function needs stanmat to be an S4 object, which is a pain so just send it manually
    if(object$algorithm == "bootstrap") {
      pp.eta = rstanarm:::pp_eta(object, pp.data, stanmat = as.matrix(obj.draws)) # would like to fix this if possible
    } else {
      pp.eta = rstanarm:::pp_eta(object, pp.data)
    }

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

  # add variable names
  predict.df = dplyr::left_join(dplyr::select(dplyr::ungroup(predict.df), -.extra.name), pr$predictions, by = ".prediction.id")

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
    }
    else if(l == 4) { # (1 - 2) - (3 - 4)
      r = (predict.df$prob[predict.df$.prediction.id == x[1]] - predict.df$prob[predict.df$.prediction.id == x[2]]) -
        (predict.df$prob[predict.df$.prediction.id == x[3]] - predict.df$prob[predict.df$.prediction.id == x[4]])
    }
    else {
      r = NA
    }

    # return
    return(r)
  }

  # times
  times = nrow(predict.df) / nrow(pr$predictions)

  # get contrasts
  contrast.df =
    tibble::tibble(
      .contrast = rep(names(pr$contrasts), each = times),
      .prediction.id = rep(format(pr$contrasts), each = times),
      prob = as.numeric(sapply(pr$contrasts, get_contrasts, predict.df = predict.df))
    )

  # if we have time add it
  if(rlang::has_name(predict.df, ".time")) {
    contrast.df$.time = predict.df$.time[predict.df$.prediction.id == 1]
  }

  # return
  return(contrast.df)
}

#' Function to get results from an object returned by analysis.
#'
#' This function allows you to get results from an analysis.
#' @param object Object returned from 'analysis' function.
#' @param predictions The predictions returned from "pr_list."
#' @param method Method for producing predictions Either uses observed values (the default) or mean values.
#' @param times Vector of times to produce predictions for. Only used for survival analysis.
#' @param .full.matrix Whether to return structured predictions (the default) or the full matrix.
#' @keywords bootstrap, results, prediction
#' @export
#' @examples
#' results(object = output, predictions = main.predictions)
#'

# function to get predictions from an "out" object produced by analyze
results = function(object, predictions = NULL, method = c("observed values", "mean"), times = NULL, .full.matrix = F) {
  ## get needed variables

  # get results from the model run
  obj.draws = get_draws(object)

  # get formula used
  obj.formula = get_formula(object)

  # get the data used for running the model
  obj.data = get_data(object)

  # get structured predictions
  if(!is.null(predictions)) {
    pr = structure_predictions(obj.formula, obj.data, predictions)
  } else {
    pr = NULL
  }


  ## create coefficient frames

  # get coefficients
  r.coefs = results_coefficients(obj.formula, obj.draws, obj.data)

  # summarize to get coefficient values
  r.coefs = summarize_interval(r.coefs)


  ## create predictions frame

  # get predictions
  predict.df = results_predictions(object, pr, obj.draws, method, times)

  # produce predictions that include the prediction id
  r.preds = summarize_interval(predict.df, .grouping = c("prob"))


  ## create contrasts frame

  # get contrasts
  contrast.df = results_contrasts(pr, predict.df)

  # summarize
  r.contr = summarize_interval(contrast.df, .grouping = c("prob"))


  ## return stuff

  # set return
  if(.full.matrix) {
    r = list(coefficients = r.coefs, predictions = predict.df, contrasts = contrast.df)
  } else {
    r = list(coefficients = r.coefs, predictions = r.preds, contrasts = r.contr)
  }

  # clean memory
  rm(r.coefs, r.preds, r.contr, predict.df, contrast.df)

  # return
  return(r)
}
