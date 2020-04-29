# analyze functions

#' Function to easily create a nicely formatted prediction list.
#'
#' This function allows you to create formatted predictions and contrasts.
#' @param ... Variables that will be used to create the predictions.
#' @param .constant Variable(s) to hold constant within a contrast.
#' @param .diff.in.diff Variable(s) to use to create a "difference-in-difference" contrast.
#' @param .add Text to add to a contrast.
#' @keywords prediction contrast hypothesis
#' @export
#' @examples
#' pr_list(treatment = c(1, 0), condition = c("High", "Low"), .constant = "condition")
#'

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

# generate resamples -- http://jee3.web.rice.edu/cluster-paper.pdf
# resample hierarchical data: http://biostat.mc.vanderbilt.edu/wiki/Main/HowToBootstrapCorrelatedData + https://stats.idre.ucla.edu/r/dae/mixed-effects-logistic-regression/

#' Function to resample a dataframe that takes into account user-specified weights and clusters.
#'
#' This function allows you to resample a dataframe for analysis.
#' @param data Dataframe with variables to be resampled
#' @param n Number of resamples to generate.
#' @param cluster Formula specifying the variable in data that identifies clusters within the data.
#' @param weights Vector of weights to use.
#' @keywords bootstrap resample cluster weight
#' @export
#' @examples
#' resample_df(data, n = 2000, cluster = ~ cluster, weights = data$.weights)
#'

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

#' Main entry function to conduct analysis of a dataset.
#'
#' This function runs either a frequentist or bayesian analysis of a dataset.
#' @param runs The number of runs to conduct. Optionally, an object returned from "resample_df."
#' @param formula The formula to run.
#' @param data The data to use for analysis.
#' @param cluster Formula specifying the variable in data that identifies clusters within the data.
#' @param weights Vector of weights to use.
#' @param model.type The model type is determined automatically from the formula and data. Can optionally be specified directly.
#' @param model.extra.args Extra arguments to the model that will be run.
#' @param inference The type of inference to use: frequentist (the default) or bayesian.
#' @keywords boostrap analysis frequentist bayesian
#' @export
#' @examples
#' analysis(runs = 1000, formula = out ~ treat, data = dt, inference = "bayesian")
#'

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

