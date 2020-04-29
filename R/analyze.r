# analyze functions

# function to easily create prediction list
pr_list = function(..., .constant = NULL, .diff.in.diff = NULL, .add = NULL) {
  # expand to create combinations -- moving to expand_grid (or removing strings as factors) fixes a problem with predictions being reversed
  pr = dplyr::as_tibble(tidyr::expand_grid(...))

  # check to make sure that all the vars are present
  if(!is.null(.diff.in.diff) & !all(.diff.in.diff %in% colnames(pr))) {
    stop("Missing vars in 'diff-in-diff'")
  }

  # check to make sure that all the vars are present
  if(!is.null(.constant) & !all(.constant %in% colnames(pr))) {
    stop("Missing vars in 'constant'")
  }

  # sort to make prettier predictions
  pr = dplyr::select(pr, !!.diff.in.diff, dplyr::everything(), !!.constant)

  # set the names
  pr$.name = sapply(1:nrow(pr), function(x) return(paste(sapply(pr[x, ], format, scientific = F, big.mark = ",", trim = T), collapse = ", ")))

  # extra model info
  if(!is.null(.add)) pr$.model = .add

  # extra variables to hold constant
  .extra = colnames(pr)[!colnames(pr) %in% c(.constant, .diff.in.diff, ".name", ".model")]

  # are we doing a two-way or four-way comparison
  ways = dplyr::if_else(!is.null(.diff.in.diff), 4, 2)

  # check diff-in-diff and constant
  if(ways == 4) {
    # create a data frame of all combinations
    pr.list = tidyr::expand_grid(x1 = 1:nrow(pr), x2 = 1:nrow(pr), x3 = 1:nrow(pr), x4 = 1:nrow(pr))

    # find contrasts that meet our constant requirement (all extra variables are the same within the category and the constant variable is the same across everything)
    pr.list$.constant =
      apply(pr[pr.list$x1, .extra] == pr[pr.list$x2, .extra], 1, all) &
      apply(pr[pr.list$x3, .extra] == pr[pr.list$x4, .extra], 1, all) &
      apply(pr[pr.list$x1, .constant] == pr[pr.list$x2, .constant], 1, all) &
      apply(pr[pr.list$x2, .constant] == pr[pr.list$x3, .constant], 1, all) &
      apply(pr[pr.list$x3, .constant] == pr[pr.list$x4, .constant], 1, all)

    # find contrasts that meet our diff-in-diff requirement
    pr.list$.diff.in.diff =
      apply(pr[pr.list$x1, .diff.in.diff] == pr[pr.list$x3, .diff.in.diff], 1, all) &
      apply(pr[pr.list$x2, .diff.in.diff] == pr[pr.list$x4, .diff.in.diff], 1, all) &
      apply(pr[pr.list$x1, .diff.in.diff] != pr[pr.list$x2, .diff.in.diff], 1, all)

    # make sure they are not all the same
    pr.list$.not.same = !(pr.list$x1 == pr.list$x2 | pr.list$x3 == pr.list$x4 | (pr.list$x1 == pr.list$x3 & pr.list$x2 == pr.list$x4) | (pr.list$x1 == pr.list$x4 & pr.list$x2 == pr.list$x3))

    # identify entries that test the same thing
    pr.list$.unique = paste(pr.list$x1, pr.list$x2, pr.list$x3, pr.list$x4) == mapply(
      function(x1, x2, x3, x4) {
        r = t(matrix(c(x1, x2, x3, x4, x2, x1, x4, x3, x3, x4, x1, x2, x4, x3, x2, x1, x1, x2, x4, x3, x2, x1, x3, x4), nrow = 4))
        r = r[order(r[, 1], r[, 2], r[, 3], r[, 4]), ]
        return(paste(r[1, ], collapse = " "))
      }, pr.list$x1, pr.list$x2, pr.list$x3, pr.list$x4)
  } else {
    # create a data frame of all combinations
    pr.list = tidyr::expand_grid(x1 = 1:nrow(pr), x2 = 1:nrow(pr))

    # find contrasts that meet our constant requirement
    pr.list$.constant = apply(pr[pr.list$x1, .constant] == pr[pr.list$x2, .constant], 1, all)

    # find contrasts that meet our diff-in-diff requirement
    pr.list$.diff.in.diff = T

    # make sure they are not all the same
    pr.list$.not.same = !(pr.list$x1 == pr.list$x2)

    # identify entries that test the same thing
    pr.list$.unique = paste(pr.list$x1, pr.list$x2) == mapply(function(x1, x2) if(x1 < x2) return(paste(x1, x2)) else return(paste(x2, x1)), pr.list$x1, pr.list$x2)
  }

  # select only good contrasts
  pr.list = dplyr::filter(pr.list, pr.list$.constant & pr.list$.diff.in.diff & pr.list$.not.same & pr.list$.unique)

  # now create the list of contrast lists
  pr.full = lapply(1:nrow(pr.list), function(x) {
    if(ways == 4) {
      return(list(pr[pr.list$x1[x], ], pr[pr.list$x2[x], ], pr[pr.list$x3[x], ], pr[pr.list$x4[x], ]))
    } else {
      return(list(pr[pr.list$x1[x], ], pr[pr.list$x2[x], ]))
    }
  })

  # make names vector
  names = sapply(1:nrow(pr.list), function(x) {
    if(ways == 4) {
      return(paste0(pr[pr.list$x1[x], ".name"], " <- ", pr[pr.list$x2[x], ".name"], " vs. ", pr[pr.list$x3[x], ".name"], " <- ", pr[pr.list$x4[x], ".name"], if(is.null(.add)) "" else paste0(" [", .add, "]")))
    } else {
      return(paste0(pr[pr.list$x1[x], ".name"], " vs. ", pr[pr.list$x2[x], ".name"], if(is.null(.add)) "" else paste0(" [", .add, "]")))
    }
  })

  # add names
  names(pr.full) = names

  # return
  return(pr.full)
}

# summarize interval -- could use HDI -- should be equivalent to http://mc-stan.org/rstanarm/reference/posterior_interval.stanreg.html
summarize_interval = function(data, .ci = 0.95, .grouping = c(".prediction.id", "prob"), .round = T) {
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
  if(.round == T) {
    r = dplyr::mutate_if(r, is.numeric, round, 3)
  }

  # return
  return(r)
}

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

# generate resamples -- http://jee3.web.rice.edu/cluster-paper.pdf
# resample hierarchical data: http://biostat.mc.vanderbilt.edu/wiki/Main/HowToBootstrapCorrelatedData + https://stats.idre.ucla.edu/r/dae/mixed-effects-logistic-regression/
resample_df = function(data, n = 1000, cluster = NULL, weights = NULL) {
  # the basic idea is to resample clusters based on average weight and then resample observations within a clsuter based on individual weights

  # add row id
  data = tibble::rowid_to_column(data, ".row.id")

  # set weights
  data$.weights = if(is.null(weights)) { 1 } else { weights }

  # set cluster
  if(is.null(cluster)) cluster = ~ .row.id

  # select vars
  clust = dplyr::group_by_at(dplyr::select(data, c(".row.id", ".weights", all.vars(cluster))), dplyr::vars(all.vars(cluster)))

  # add group id
  clust$.group.id = dplyr::group_indices(clust)

  # sort weights to help with sample
  clust = dplyr::arrange(clust, .group.id, .weights)

  # set unique groups with weights equal to the average weight in a group
  groups = dplyr::summarize(dplyr::group_by(clust, .group.id), .group.length = dplyr::n(), .group.weights = sum(.weights) / dplyr::n())

  # get the indices
  if(nrow(groups) == nrow(data)) {
    # no clusters
    i.smp = sample.int(n = nrow(groups), size = n * nrow(groups), replace = T, prob = groups$.group.weights)
  } else {
    # sample groups
    g.smp = sample.int(n = nrow(groups), size = n * nrow(groups), replace = T)

    # sample within groups
    i.smp = lapply(g.smp, function(x) {
      # get length
      length = groups$.group.length[groups$.group.id == x]

      # get indices for a group -- use the internal function as it is faster than calling with overhead
      id = 1:length
      # id = .Internal(sample(length, length, T, clust$.weights[clust$.group.id == x]))
      # id = sample.int(length, replace = T)

      # return
      return(clust$.row.id[clust$.group.id == x][id])
    })
  }

  # get actual index
  idx.full = lapply(1:n, function(x) unlist(i.smp[((x-1) * nrow(groups) + 1):(x * nrow(groups))]))

  # check lengths
  # summary(sapply(idx.full, length))

  # return
  return(idx.full)
}

# remove backticks function
remove_backticks = function(str) {
  r = sapply(str, function(x) { if(substring(x, 1, 1) == "`") x = substring(x, 2); l = nchar(x); if(substring(x, l, l) == "`") x = substring(x, 1, l - 1); x })
  r = as.character(r)
  r
}

# get draws (model results) for each run
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

# get formula from object
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

# get data
get_data = function(object, data = NULL) {
  # check
  if(is.null(object) || typeof(object) != "list" || !rlang::has_name(object, "data")) {
    return(data)
  }

  # return
  return(object$data)
}

# set coefficients
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

# set predictions
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

# set contrasts
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

# new basic model function that just runs the model
# what is needed: libraries, data, formula, model, model args, converge check, cleaner
model_run.frequentist = function(idx, formula.parsed, model.functions) {
  # set data
  df = formula.parsed$data[idx, ]

  # do we require any libraries
  if(!is.null(model.functions$libraries)) {
    require(model.functions$libraries, character.only = T)
  }

  # run the call
  m = NULL
  try(m <- do.call(model.functions$model.run, c(list(formula = formula.parsed$formula, data = df), model.functions$model.args)), T)

  # did the model run?
  if(is.null(m) | length(m) == 0) {
    return(NULL)
  }

  # check if it converged
  if(!is.null(model.functions$converged)) {
    if(!model.functions$converged(m)) {
      return(NULL)
    }
  }

  # get the model results
  coefs = model.functions$coefs(m)

  # set alpha and beta
  names = names(coefs)
  alpha = if(any(names == "(Intercept)")) coefs[names == "(Intercept)"] else NULL
  beta = coefs[names != "(Intercept)"]

  # get auxiliary var
  aux = model.functions$special(m)

  # get mean fitted values
  mean_PPD = model.functions$fitted(m)

  # get the log likelihood -- the lower the better
  lp_ = model.functions$performance(m)

  # return
  return(list(alpha = alpha, beta = beta, aux = aux, mean_PPD = mean_PPD, lp_ = lp_))
}

# the model function
model_run.bayesian = function(chain.id, formula.parsed, model.functions) {
  # do we require any libraries
  if(!is.null(model.functions$libraries)) {
    lapply(model.functions$libraries, require, character.only = T)
  }

  # run the call
  m = NULL
  try(m <- do.call(model.functions$model.run, c(list(formula = formula.parsed$formula, data = formula.parsed$data, chain_id = chain.id), model.functions$model.args)), T)

  # get the fit from the object
  if(rlang::has_name(m, "stanfit")) {
    m.fit <- m$stanfit
  } else if(rlang::has_name(m, "fit")) {
    m.fit <- m$fit
  } else {
    m.fit <- NULL
  }

  # structure return
  r = list(model = m, fit = m.fit)

  # return
  return(r)
}

# run -- replicates is a list of observations in data to run on, boot_function is what is run, should return NULL for a bad run
run_all_replicates = function(replicates, boot_function, N = length(replicates), args = NULL, max.threads = 20, parallel.libraries = NULL) {
  # makes sure replicates is a list
  if(!is.list(replicates)) {
    stop("Replicates must be a list of data frames.")
  }

  # check to make sure its a function
  if(!is.function(boot_function)) {
    stop("Boot function is not a function.")
  }

  # check to make sure N is okay
  if(N < 1 | N > length(replicates)) {
    stop(paste0("Bad number of replicates generated to run ", N, " times."))
  }

  # setup parallel
  `%dopar%` <- foreach:::`%dopar%`

  # insert serial backend, otherwise error in repetetive tasks
  foreach::registerDoSEQ()

  # just to be safe
  closeAllConnections()

  # set number of threads
  num.threads = max(1, min(max.threads, parallel::detectCores(T, T), na.rm = T), na.rm = T)

  # create the cluster
  cl = parallel::makeCluster(num.threads, type = "PSOCK", methods = F, outfile = "")

  # add libraries if necessary
  if(!is.null(parallel.libraries)) {
    # parallel::clusterEvalQ(cl, { library(parallel.libraries, character.only = T) }) # doesnt work
  }

  # register the cluster
  doSNOW::registerDoSNOW(cl)

  # create progress bar
  pb = txtProgressBar(1, N, style = 3, char = "-", width = 25)
  progress = function(n) setTxtProgressBar(pb, n)
  foreach.options = list(progress = progress)

  # set time and out
  time = NULL
  out = NULL

  # error handling function
  has_error = function(e) {
    # say error
    cat(e)

    # return null
    return(NULL)
  }

  # run the bootstrap function
  tryCatch(
    time <- system.time(
      out <- foreach::foreach(i = replicates, .errorhandling = "pass", .options.snow = foreach.options) %dopar% {
        do.call(boot_function, c(list(i), args))
      }
    ), error = has_error)

  # did we run
  if(!is.null(time)) {
    # set the time
    time = time[3]

    # print how long it took
    cat(paste0(" (", round(time, 2), " seconds)"))
  }

  # close the progress bar
  close(pb)

  # stop the cluster
  parallel::stopCluster(cl)

  # insert serial backend, otherwise error in repetetive tasks
  foreach::registerDoSEQ()

  # just to be safe
  closeAllConnections()

  # memory cleanup
  rm(replicates)
  gc()

  # filter bad results
  if(length(out) > 0) {
    out[sapply(out, is.null)] = NULL
  }

  # return
  return(list(result = out, time = time))
}

# create stan fit object
create_stan_fit = function() {
  # set stuff
  model.name = "glm"

  # the categories returned
  sample.pars = c("alpha", "beta", "lp__")

  # the number of items in each category
  sample.dims = list(NULL, 3, NULL)
  names(sample.dims) = sample.pars

  # simulations
  sim = list(samples = samples, # need to fix
             iter = runs,
             thin = 1,
             warmup = 0,
             chains = 1,
             n_save = runs,
             warmup2 = 0,
             permutation = list(sample.int(runs)),
             pars_oi = sample.pars,
             dims_oi = dims_oi,
             fnames_oi = dotfnames_to_sqrfnames(par_fnames),
             n_flatnames = length(par_fnames))

  # create null args
  null.args = list(chain_id = 1, iter = runs, thin = 1, seed = seed, warmup = 0, init = "random", algorithm = "bootstrap", save_warmup = F, method = "bootstrap", control = list())

  # create a null dso to pass to the null model
  null.dso = new("cxxdso", sig = list(character(0)), dso_saved = F, dso_filename = character(0), modulename = character(0),
                  system = R.version$system, cxxflags = character(0), .CXXDSOMISC = new.env(parent = emptyenv()))

  # create a null model
  null.model = new("stanmodel", model_name = model.name, model_code = character(0), model_cpp = list(), dso = null.dso)

  # create the fit in the proper format
  boot.fit = new("stanfit",
      model_name = model.name,
      model_pars = sample.pars,
      par_dims = sample.dims,
      mode = mode,
      sim = sim,
      inits = list(),
      stan_args = null.args,
      stanmodel = null.model,
      date = sdate,
      .MISC = new.env(parent = emptyenv()))

  # return
  return(boot.fit)
}

# create a stan summary object
create_stan_summary = function(full.matrix, full.names) {
  # create a long version
  summary.long = tidyr::pivot_longer(tibble::as_tibble(full.matrix[, c(full.names$alpha, full.names$beta, full.names$aux, full.names$mean_PPD, full.names$lp_)]), dplyr::everything())

  # filter out NAs
  summary.long = dplyr::filter(summary.long, !is.na(value))


  summary.long = dplyr::summarize_all(
    dplyr::group_by(summary.long, name),
    list(mean = mean,
         se_mean = function(x) sd(x)/sqrt(length(x)),
         sd = sd,
         `2.5%` = function(x) quantile(x, 0.025),
         `10%` = function(x) quantile(x, 0.1),
         `25%` = function(x) quantile(x, 0.25),
         `50%` = function(x) quantile(x, 0.5),
         `75%` = function(x) quantile(x, 0.75),
         `90%` = function(x) quantile(x, 0.9),
         `95%` = function(x) quantile(x, 0.95),
         `97.5%` = function(x) quantile(x, 0.975),
         n_eff = function(x) length(x)))

  # save row names
  row.names = summary.long$name

  # get rid of names column and reorder
  summary.long = as.matrix(summary.long[, -1])
  summary.long = summary.long[match(c(full.names$alpha, full.names$beta, full.names$aux, full.names$mean_PPD, full.names$lp_), row.names), ]

  # set rownames
  rownames(summary.long) = row.names[match(c(full.names$alpha, full.names$beta, full.names$aux, full.names$mean_PPD, full.names$lp_), row.names)]

  # return
  return(summary.long)
}

# create model string summary
create_model_string = function(model.type, formula.parsed, cluster, weights, inference) {
  # create model string
  model.string = dplyr::case_when(
    formula.parsed$random.effects && formula.parsed$fixed.effects ~ "Random and Fixed Effects",
    formula.parsed$random.effects ~ "Random Effects",
    formula.parsed$fixed.effects ~ "Fixed Effects",
    T ~ NA_character_
  )

  # create model string
  se.string = dplyr::case_when(
    !is.null(cluster) && !is.null(weights) ~ "Clustered and Weighted",
    !is.null(cluster) && is.null(weights) ~ "Clustered",
    is.null(cluster) && !is.null(weights) ~ "Weighted",
    T ~ NA_character_
  )

  # create inference string
  inf.string = dplyr::case_when(
    inference == "bayesian" ~ "Bayesian",
    T ~ "Frequentist"
  )

  # assemble string
  return.string =
    paste0("Model: ", model.type,
           if(!is.na(model.string)) paste0(" with ", model.string) else "",
           if(!is.na(se.string)) paste0(", Data: ", se.string) else "",
           ", Method: ", inf.string)

  # return
  return(return.string)
}

# basic logic:
# create N dataframes that incorporate weights and clusters (accounts for effect of weights and clusters)
# run the model N times on all dataframes and save the cleaned model object
# run prediction and contrast on base data for each model run
# return results

# the function to run the model and save results + model info
analysis = function(runs, formula, data, cluster = NULL, weights = NULL, model.type = NULL, model.extra.args = NULL, inference = c("frequentist", "bayesian")) {
  # select only relevant variabls from data
  data = dplyr::select(data, all.vars(formula))
  data = dplyr::filter(data, complete.cases(data))

  # parse the formula
  formula.parsed = parse_formula(formula = formula, data = data)

  # set model type if we need to
  if(is.null(model.type)) {
    # identify model type automatically
    model.type = determine_model(formula.parsed = formula.parsed, data = data)
  } else {
    # check to make sure it is specified correctly
    if(!model.type %in% unlist(model.type.list)) {
      stop(paste("Incorrect model type --", model.type, "- specified"))
    }
  }

  # set inference
  if(dplyr::first(inference) %in% c("Bayesian", "bayesian", "Bayes", "bayes")) {
    inference = "bayesian"
  } else {
    inference = "frequentist"
  }

  # print
  cat(paste("Running --", create_model_string(model.type = model.type, formula.parsed = formula.parsed, cluster = cluster, weights = weights, inference = inference), "-- analysis\n"))

  # get model extra
  model.functions = produce_model_function(model.type = model.type, formula.parsed = formula.parsed, inference = inference, model.extra.args = model.extra.args)

  # run the model -- either bayesian or frequentist

  # main division
  if(inference == "bayesian") {
    # print
    cat("Running model\n")

    # set resamples
    resamples = as.list(1:model.functions$model.args$cores)
    model.functions$model.args$cores = 1

    # set iterations
    model.functions$model.args$iter = runs

    # run multicore
    out = run_all_replicates(
      replicates = resamples,
      boot_function = model_run.bayesian,
      args = list(formula.parsed = formula.parsed, model.functions = model.functions),
      parallel.libraries = model.functions$libraries
    )

    # combine the chains together
    if(length(out$result) > 0) {
      # first our return frame
      out.result = out$result[[1]]$model

      # combine the fits
      all.fits = rstan::sflist2stanfit(purrr::map(out$result, "fit"))

      # set the combined model fit -- different for rstanarm and brms
      if(rlang::has_name(out.result, "stanfit")) out.result$stanfit = all.fits else out.result$fit = all.fits
    } else {
      out.result = NULL
    }

    # add back to out
    out = out.result
  } else {
    # do we have a number of runs or the actual resamples
    if(is.list(runs)) {
      resamples = runs
    } else {
      # print
      cat(paste("Generating --", runs, "-- replicates of data"))

      # generate data resamples
      time.resample = system.time(resamples <- resample_df(data = data, n = runs, cluster = cluster, weights = weights))

      # print time
      cat(paste0(" (", round(time.resample[3], 2), " seconds)\n"))
    }

    # print
    cat("Running model\n")

    # run model on full data initially
    base.model = do.call(model.functions$model.run, c(list(formula = formula.parsed$formula, data = formula.parsed$data), model.functions$model.args))

    # run the replicates
    out =
      run_all_replicates(
        replicates = resamples,
        boot_function = model_run.frequentist,
        args = list(formula.parsed = formula.parsed, model.functions = model.functions),
        parallel.libraries = model.functions$libraries
      )

    # check if something went wrong
    if(!is.list(out$result) && !is.list(out$result[[1]]) && length(out$result[[1]]) != 5) {
      print(out$result[[1]])
      stop("Model did not run. Please respecify and try again.")
    }

    # save formula and family
    out$formula = formula.parsed$formula
    out$terms = terms(formula.parsed$formula)
    out$family = model.functions$family

    # full matrix
    full.matrix = t(sapply(out$result, function(x) { .Internal(unlist(c(x$alpha, x$beta, x$aux, x$mean_PPD, x$lp_), F, F)) }))

    # fix colnames by getting rid of backticks
    colnames(full.matrix) = remove_backticks(colnames(full.matrix))

    # get all names
    full.names = lapply(out$result[[1]], function(x) remove_backticks(names(x)))

    # save pseudo draws object based on coefs from runs -- just alpha, beta, and aux
    out$stanfit = full.matrix[, c(full.names$alpha, full.names$beta, full.names$aux)]

    # set coefficients
    out$coefficients = apply(out$stanfit[, c(full.names$alpha, full.names$beta)], 2, median, na.rm = T)
    out$ses = apply(out$stanfit[, c(full.names$alpha, full.names$beta)], 2, sd, na.rm = T)

    # create a stan_summary
    out$stan_summary = create_stan_summary(full.matrix, full.names)

    # save data
    out$data = formula.parsed$data

    # save the model frame
    out$model = formula.parsed$data

    # set offset
    out$offset = rep(0, nrow(formula.parsed$data))

    # set algorithm -- not sure if we want to do this
    out$algorithm = c("optimizing", "bootstrap")

    # set draws
    out$asymptotic_sampling_dist = out$stanfit

    # set model specific stuff
    if(model.type == model.type.list$survival) {
      # get baseline hazard
      out$basehaz = survival::basehaz(base.model, centered = F)

      # get info about Y
      model.y = as.matrix(base.model$y)
      colnames(model.y) = as.character(base.model$terms[[2]])[-1]

      # set entrytime, eventtime, event, and delayed
      out$entrytime = model.y[, "start"]
      out$eventtime = model.y[, "stop"]
      out$event = as.logical(model.y[, "status"])
      out$delayed = model.y[, "start"] != 0
    }

    # set class
    class(out) = c("stanreg", "bootstrap")
  }

  # free memory
  rm(formula.parsed, model.functions, resamples)

  # do main garbage collection after each full model run
  gc()

  # return
  return(out)
}

