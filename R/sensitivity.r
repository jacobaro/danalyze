# libraries
needs::needs(tidyverse)

# sensitivity tests

# try the de-confounder approach
# basically include a variable that captures the "distribution of assigned causes" which is a stand-in for unobserved factors influencing our treatments
# a posterior predictive check makes sure that the latent variable actually does a good job of capturing the assigned causes (must pass to use)
# if the latent variable is good, this means that the causes are conditionally independent given the latent factor (assumption is no single-cause variables included)

# article: https://arxiv.org/pdf/1805.06826.pdf

# python code: https://github.com/blei-lab/deconfounder_tutorial/blob/master/deconfounder_tutorial.ipynb

# pca methods: https://github.com/hredestig/pcaMethods

# the logic:
# (1) remove causes that are highly correlated with other causes (merge highly correlated causes)
# (2) standardize data used for principal component analysis
# (3) holdout some data for posterior predictive checks used to assess factor model
# (4) run a random factor or probabilistic PCA model on the standardized data to infer a latent variable that is the deconfounder
#     (continuous for continuous variables and factor for factor variables to enforce overlap)
# (5)

# load the data
test_load_data = function(drop.cor = T) {
  # column names
  cnames = c(
    "id",
    "outcome",
    "radius_mean",
    "texture_mean",
    "perimeter_mean",
    "area_mean",
    "smoothness_mean",
    "compactness_mean",
    "concavity_mean",
    "concave_points_mean",
    "symmetry_mean",
    "fractal_dimension_mean",
    "radius_standard_error",
    "texture_standard_error",
    "perimeter_standard_error",
    "area_standard_error",
    "smoothness_standard_error",
    "compactness_standard_error",
    "concavity_standard_error",
    "concave_points_standard_error",
    "symmetry_standard_error",
    "fractal_dimension_standard_error",
    "radius_worst",
    "texture_worst",
    "perimeter_worst",
    "area_worst",
    "smoothness_worst",
    "compactness_worst",
    "concavity_worst",
    "concave_points_worst",
    "symmetry_worst",
    "fractal_dimension_worst"
  )

  # load
  d = readr::read_csv("test/wdbc.data.csv", col_names = cnames, col_types = paste(c("dc", rep("d", 30)), collapse = ""))

  # set outcome
  d$outcome = factor(d$outcome, levels = c("B", "M"), labels = c("Benign", "Malignant"))

  # set outcome binary
  d$outcome.binary = dplyr::if_else(d$outcome == "Benign", 1, 0)

  # select
  d = dplyr::select(d, outcome, outcome.binary, radius_mean:fractal_dimension_mean)

  # check
  summary(d)

  # drop correlated variables as per blei's code
  if(drop.cor) {
    # check correlation
    cor = cor(dplyr::select(d, -outcome, -outcome.binary))

    # find highly correlated variables -- this matches the result in the example code
    # to.drop = caret::findCorrelation(cor, cutoff = 0.95, names = T)
    to.drop = c("perimeter_mean", "area_mean")

    # select out
    d = dplyr::select(d, -dplyr::all_of(to.drop))

    # check
    summary(d)
  }

  # return
  return(d)
}

# create holdout mask
create_holdout_mask = function(n.col, n.row, portion = 0.2) {
  # create training and validation data -- the example code sets 20% of data points (rows and columns) to zero
  holdout.mask = matrix(sample(c(0, 1), size = n.col * n.row, replace = T, prob = c(1 - portion, portion)), nrow = n.row, ncol = n.col)

  # return
  return(holdout.mask)
}

# function to holdout data
holdout_data = function(data, portion = 0.20) {
  # rows and columns
  n.col = ncol(data)
  n.row = nrow(data)

  # get mask
  holdout.mask = create_holdout_mask(n.col = n.col, n.row = n.row, portion = portion)

  # identify rows chosen
  rows.chosen = which(rowSums(holdout.mask) > 0)

  # training and validation data
  d.val = data * holdout.mask
  d.train = data * (1 - holdout.mask)

  # return
  return(list(train = d.train, validation = d.val, mask = holdout.mask, rows = rows.chosen))
}

# run factor model
factor_model = function(formula, data, latent = 2, model.type = c("ppca", "linear", "quadratic"), seed = as.integer(Sys.time()), data.std.dev = 0.1, z.std.dev = 2.0) {
  # get model type
  model.type = match.arg(arg = model.type, choices = c("ppca", "linear", "quadratic"))

  # select the data
  data.std = dplyr::select(data, all.vars(formula))

  # parse the formula
  f.parts = stringr::str_split(attr(terms(formula), "term.labels"), "\\|")[[1]]

  # do we have a unit
  if(length(f.parts) > 1) {
    formula.use = as.formula(paste0("~", f.parts[1]))
    formula.unit = all.vars(as.formula(paste0("~", f.parts[2])))
  } else {
    formula.use = formula
    formula.unit = NULL
  }

  # identify complete cases
  data.rows = complete.cases(data.std)

  # turn into a model matrix and remove intercept
  data.std = tibble::as_tibble(model.matrix(formula.use, data.std))[, -1]

  # standardize
  data.std =  dplyr::mutate_all(data.std, function(x) as.numeric(scale(x)))

  # remove columns that are all NA
  data.std = data.std[, apply(data.std, 2, function(x) !all(is.na(x)))]

  # make sure we still have data
  if(nrow(data.std) == 0 || ncol(data.std) < 2) {
    return(NULL)
  }

  # get data
  holdout = holdout_data(data.std)

  # helper functions -- not needed currently
  normal_lpdf = function(x, mean, sd) { return(sum(dnorm(x, mean, sd, log = T), na.rm = T)) }

  # bootstrap a ppca model
  if(model.type == "ppca") {
    # the logic for each run: get train/validation split, create latent variable using training data, check fit for validation data

    # boot function
    factor_boot = function(data, latent, portion, create_holdout_mask) {
      # create mask
      # mask = create_holdout_mask(n.col = ncol(data), n.row = nrow(data), portion = portion)
      mask = sample.int(nrow(data), as.integer(nrow(data) * portion), replace = F)

      # training and validation data -- just use observation-based out-of-sample
      # data.train = as.matrix(data * mask)
      # data.validation = as.matrix(data * (1 - mask))
      data.train = as.matrix(data[-mask, ])
      data.validation = as.matrix(data[mask, ])

      # need to install BiocManager, cate + BiocManager::install("sva")

      # fit the factor
      factor.fit = cate::factor.analysis(Y = data.train, r = latent, method = "ml")

      # get scores
      scores = as.matrix(factor.fit$Gamma)
      colnames(scores) = paste0("z", seq(latent))

      # get validation prediction
      validation.scores = data.validation %*% scores

      # check model fit by examining variation explained (r^2) in holdout data
      r = apply(data.validation, 2, function(Y) summary(lm(Y ~ validation.scores))$r.squared)

      # score dataframe
      scores.df = tibble::tibble(variable = rep(rownames(scores), times = ncol(scores)), confounder = rep(colnames(scores), each = nrow(scores)), prob = as.vector(scores))

      # return
      return(list(scores = scores.df, q2 = r))
    }

    # boot it
    r = danalyze::run_all_replicates(seq(500), boot_function = factor_boot, args = list(data = data.std, latent = latent, portion = 0.2, create_holdout_mask = create_holdout_mask))

    # get scores
    W = purrr::map_dfr(r$result, "scores")
    W = danalyze::summarize_interval(W)

    # get Z -- the deconfounder
    Z = tidyr::pivot_wider(dplyr::select(W, variable, confounder, c), names_from = confounder, values_from = c)
    Z = Z[match(colnames(data.std), Z$variable), -1]
    confounder = as.matrix(data.std) %*% as.matrix(Z)

    # check correlation
    confounder.cor = as.matrix(cor(data.std, confounder))

    # get q2 for factor model -- not p-value but explaind variation based on out-of-sample
    p.value.train = tidyr:::pivot_longer(purrr::map_dfr(r$result, "q2"), everything(), names_to = "variable", values_to = "prob")
    p.value.train = danalyze::summarize_interval(p.value.train)

    # create data frame to return -- this allows data.std to have missing rows
    confounder = dplyr::left_join(tibble::tibble(.row = seq(nrow(data))),
                                  dplyr::bind_cols(tibble::tibble(.row = seq(nrow(data))[data.rows]), tibble::as_tibble(confounder)),
                                  by = ".row")
    confounder = dplyr::select(confounder, -.row)

    # set remaining
    a.hat = a.hat.cor = NULL
  } else {
    # setup stan model -- preloaded and saved as RDS
    if(model.type == "linear") {
      # model = rstan::stan_model(file = "test/linearfactor_knownvar.stan")
      # saveRDS(model, "test/linear_factor_stan.rds", compress = T)
      model = readRDS("test/linear_factor_stan.rds")
    } else if(model.type == "quadratic") {
      # model = rstan::stan_model(file = "test/quadraticfactor_knownvar.stan")
      # saveRDS(model, "test/quadratic_factor_stan.rds", compress = T)
      model = readRDS("test/quadratic_factor_stan.rds")
    }

    # model data
    args = list(N = nrow(data.std), D = ncol(data.std), K = latent, X_train = holdout$train, X_vad = holdout$validation, holdout_mask = holdout$mask, X_all = data.std, data_std = data.std.dev, Z_std = z.std.dev)

    # fit the model
    factor.map = rstan::optimizing(object = model, data = args, as_vector = F, iter = 1000, seed = seed)

    # get the point estimate
    factor.init = factor.map$par

    # draw from the posterior using variational inference
    # old: factor.fit = rstan::vb(object = model, data = args, init = factor.init, iter = 10000, eta = 0.25, adapt_engaged = 0, tol_rel_obj = 1e-2, output_samples = 5, seed = seed)

    # set arguments
    args.vb = list(object = model, data = args, init = factor.init, iter = 10000, tol_rel_obj = 1e-2, output_samples = 100, seed = seed, eta = 0.25, adapt_engaged = 0)

    # run it multiple times
    # factor.fit = danalyze:::run_all_replicates(as.list(1:6), function(i, ...) do.call(rstan::vb, args = list(...)), args = args.vb)
    factor.fit = danalyze:::run_all_replicates(seq(6), rstan::vb, args = args.vb)
    # factor.fit = do.call(rstan::vb, args = args.vb)

    # combine posterior
    factor.fit = rstan::sflist2stanfit(factor.fit$result)

    # get samples from fitted object -- Z (latent variable), X_pred (reconstructed causes), rep_lp (log posterior for prediction), and vad_lp (log posterior for validation)
    la = rstan::extract(factor.fit)

    # next do posterior predictive check

    # what portion of the predictions for an observation have below mean log-posterior values
    p.value.train = mean(sapply(holdout$rows, function(i) mean(la$rep_lp[, i] <= mean(la$vad_lp[, i]))))
    # p.value.full = mean(apply(la$rep_lp, 1, sum) < mean(apply(la$vad_lp, 1, sum)))

    # we want a p-value of >= 0.1

    # get the reconstructed confounder (full return is in la$Z)
    confounder = apply(la$Z, c(2, 3), mean, na.rm = T)
    colnames(confounder) = paste0("z", seq(latent))
    confounder = tibble::as_tibble(confounder)

    # check correlation
    confounder.cor = as.matrix(cor(data.std, confounder))

    # get reconstructed causes
    a.hat = apply(la$X_pred, c(2, 3), mean, na.rm = T)
    colnames(a.hat) = paste0(colnames(data.std), "_x")
    a.hat = tibble::as_tibble(a.hat)

    # check correlation
    a.hat.cor = as.matrix(cor(data.std, a.hat))
  }

  # return value
  r = list(z = confounder, z.cor = confounder.cor, a = a.hat, a.cor = a.hat.cor, p.train = p.value.train)

  # return
  return(r)
}

# test the deconfounder approach
test_deconfounder = function() {
  # load data
  data = test_load_data()

  # main formula
  f = ~ radius_mean + texture_mean + smoothness_mean + compactness_mean + concavity_mean + concave_points_mean + symmetry_mean + fractal_dimension_mean

  # produce substitute confounder
  cf = factor_model(formula = f, data = data, latent = 2, model.type = "linear")

  # check p-value and correlation

  # p-value should be >= 0.1
  cf$p.train

  # correlation should be ~= 0.5 (lower means it does a poor job deconfounding and higher means high variance)
  cf$z.cor

  # add to data
  data.run = dplyr::bind_cols(dplyr::mutate_at(data, all.vars(f), function(x) as.numeric(scale(x))), cf$z)

  # run -- run just on the training data -- could bootstrap across samples of the deconfounder to get sensitivity
  m.base = glm(formula = update(f, outcome.binary ~ .), data = data.run, family = binomial)
  m.conf = glm(formula = update(f, outcome.binary ~ . + z1 + z2), data = data.run, family = binomial)

  # check
  summary(m.base)
  summary(m.conf)
}

