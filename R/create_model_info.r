# functions to set model info

# supported model types
model.type.list =
  list(
    linear = "Gaussian",
    beta = "Beta",
    binary = "Binomial",
    categorical = "Categorical",
    ordered = "Categorical (ordered)",
    poisson = "Count",
    neg.binomial = "Count (over-dispersed)",
    survival = "Event History",
    comp.risk = "Competing Risks"
  )

# function to parse the formula and data and pull out relevant components
# input: raw formula and raw data
# output: parsed formula, identificaiton of random/fixed effects, formula and data to run on
parse_formula = function(formula, data) {
  ## setup a few variables

  # get the sides of the formula
  formula.explanatory = formula(delete.response(terms(formula)))
  formula.response = update(formula, . ~ 1)
  formula.add = NULL

  ## find and deal with random effects

  # get random effects terms
  formula.random = lme4::findbars(formula.explanatory) # by hand: unlist(stringr::str_extract_all(as.character(formula[3]), "(?<=\\()[^\\)]*[\\|][^\\)|]*"))

  # remove random effects if needed
  if(!is.null(formula.random)) {
    # remove bars
    formula.explanatory = lme4::nobars(formula.explanatory)

    # add to extra
    formula.random = paste0("(", formula.random, ")")
  }

  ## find and deal wit special terms

  # special terms
  specials = c("strata", "tt", "frailty")

  # get the special terms for event history models
  time.terms = terms(formula.explanatory, specials = specials)

  # list of items to drop
  special.names = lapply(specials, survival:::untangle.specials, tt = time.terms)
  formula.survival = unlist(purrr::map(special.names, "vars"))
  to.drop = unlist(purrr::map(special.names, "terms"))

  # modify main formula
  if(!is.null(to.drop)) {
    formula.explanatory = formula(drop.terms(time.terms, to.drop, keep.response = F), env = environment(formula))
  }

  ## find and deal with fixed effects

  # first get the model frame so we can identify fixed effect
  data.frame = model.frame(formula.explanatory, data)

  # loop through and identify factors or characters with more than five categories
  fixed.effects = sapply(colnames(data.frame), function(x) (is.character(data.frame[[x]]) || is.factor(data.frame[[x]]) && dplyr::n_distinct(data.frame[[x]]) > 4))

  ## assemble the additional formulas

  # add a formula additional?
  if(!is.null(formula.random) || !is.null(formula.survival)) {
    formula.additional = formula(paste("~", paste(c(formula.random, formula.survival), collapse = " + ")), env = environment(formula))
  } else {
    formula.additional = NULL
  }

  ## assemble the full data and the full formula to test on

  # get the model matrix and the test data -- this needs to go here so that we can pull the column names
  data.full = as.matrix(model.matrix(update(formula.explanatory, ~ . + 0), data))

  # make the full formula -- the full formula is always response ~ explanatory + added
  formula.full = formula(paste(as.character(formula.response)[2], "~", paste(c(sapply(colnames(data.full), glue::backtick), formula.random, formula.survival), collapse = " + ")), env = environment(formula))
  # formula.full = formula(paste(as.character(formula.response)[2], "~", paste(c(".", formula.random, formula.survival), collapse = " + ")), env = environment(formula))

  # add response and additional to full data -- also converts to data_frame
  data.full = dplyr::bind_cols(list(dplyr::select(data, all.vars(formula.response)), as.data.frame(data.full), dplyr::select(data, all.vars(formula.add))))

  # the list that we will return

  # return list
  r = list(
    # the formula and data for running the model
    formula = formula,
    data = data,
    data.raw = data,

    # formula components
    formula.original = formula,
    response = formula.response,
    explanatory = formula.explanatory,
    added = formula.additional,

    # flags for random and fixed effects
    random.effects = !is.null(formula.random),
    fixed.effects = any(fixed.effects)
  )

  # return
  return(r)
}

# function to determine the type of model to run
# input: parsed formula and raw data
# output: model type
determine_model = function(formula.parsed, data) {
  # make character
  response.char = as.character(formula.parsed$response[2])

  # set model type
  model.type = NULL
  model.effects = NULL

  # get the dependent variable
  response.length = all.vars(formula.parsed$response)

  # no response variable
  if(length(response.length) == 1 & all(response.length == ".")) {
    response.length = NULL
  }

  # start identifying what type of model to run
  if(length(response.length) > 1) {
    # time series
    if(stringr::str_detect(response.char, "Surv")) {
      # set model type
      model.type = model.type.list$survival

      # fix formula
      if(!stringr::str_detect(response.char, "survival::")) {
        formula.parsed$response = as.formula(paste(stringr::str_replace(as.character(formula.parsed$response[2]), stringr::fixed("Surv"), "Surv"), "~ 1"), env = .GlobalEnv) # add package reference?
      }
    }
    else if(stringr::str_detect(response.char, "Event")) {
      # set model type
      model.type = model.type.list$comp.risk

      # fix formula
      if(!stringr::str_detect(response.char, "timereg::")) {
        formula.parsed$response = as.formula(paste(stringr::str_replace(as.character(formula.parsed$response[2]), stringr::fixed("Event"), "Event"), "~ 1"), env = .GlobalEnv) # add package reference?
      }
    }
    else {
      warning("Unsuported response variable.")
    }
  }
  else if(length(response.length) > 0) {
    # check the dv
    data.response = model.frame(formula.parsed$response, data)[, 1]

    # check types
    if(is.factor(data.response)) {
      # we have a binomial, categorical, or ordered model
      if(n_distinct(data.response) > 2) {
        if(is.ordered(data.response)) {
          model.type = model.type.list$ordered
        } else {
          model.type = model.type.list$categorical
        }
      } else {
        model.type = model.type.list$binary
      }
    }
    else if(all(data.response >= 0) & all(data.response %% 1 == 0)) {
      # count data

      # check for over dispersion

      # first run a poisson to see if it fits well or not
      m.t = suppressWarnings(glm(formula.parsed$formula, data, family = poisson))

      # check if residual deviance/df is much bigger than one -- also: summary(m.t)$deviance / m.t$df.residual > 1.5
      if(1 - pchisq(summary(m.t)$deviance, m.t$df.residual) < 0.05) {
        # over-dispersed when accounting for covaraiates
        model.type = model.type.list$neg.binomial
      } else {
        # not over-dispersed so run poisson
        model.type = model.type.list$poisson
      }
    } else {
      model.type = model.type.list$linear
    }
  }
  else {
    # set to null
    model.type = NULL

    # we have a problem!
    warning("No response variables.")
  }

  # return the model type
  return(model.type)
}

# function to assemble the model calls for running
# input: model type, parsed formula, extra arguments
# output: list to allow a call to model functions (run, predict, residual) and extra model libraries
produce_model_function = function(model.type, formula.parsed, inference, model.extra.args = NULL) {
  # the return list -- could make it an s4 class
  model.extra = list(
    # main model calls
    model.run = NULL,
    model.args = NULL,
    family = NULL,
    class = NULL,

    # access functions
    coefs = NULL,
    residuals = NULL,

    # additional info
    converged = NULL,
    fitted = NULL,
    performance = NULL,
    special = NULL,
    libraries = NULL,

    # additional rstanarm stuff
    cores = NULL,
    chains = NULL,
    iter = NULL,
    seed = NULL
  )

  ## set the model info list -- these are all just families for glm (negative binomial requires an extra estimation of theta)

  # set base model calls -- defaults to gaussian generalized linear model
  model.extra$model.run = stats::glm
  model.extra$model.args = list(family = stats::gaussian(link = "identity"), model = F, y = F, trace = F)
  model.extra$family = model.extra$model.args$family
  model.extra$class = c("frequentist", "glm")

  # set functions to get data
  model.extra$coefs = stats::coef
  model.extra$residuals = stats::residuals

  # set additional utility functions
  model.extra$converged = function(x) x$boundary == F & x$converged == T
  model.extra$fitted = function(x) c("mean_PPD" = mean(x$fitted.values, na.rm = T))
  model.extra$performance = function(x) c("log-posterior" = x$rank - (x$aic / 2))
  model.extra$special = function(x) c("sigma" = stats::sigma(x))

  # additions/changes for binomial
  if(model.type == model.type.list$binary) {
    model.extra$model.args$family = stats::binomial(link = "logit")
    model.extra$family = model.extra$model.args$family
  }

  # additions/changes for poisson
  if(model.type == model.type.list$poisson) {
    model.extra$model.args$family = stats::poisson(link = "log")
    model.extra$family = model.extra$model.args$family
    model.extra$class = c("frequentist", "poisson",  "glm")
  }

  # additions/changes for negative binomial
  if(model.type == model.type.list$neg.binomial) {
    if(inference == "frequentist") {
      # get theta
      theta = do.call(MASS::glm.nb, c(list(formula = formula.parsed$formula, data = formula.parsed$data, trace = 0), model.extra.args))

      # set
      model.extra$model.args$family = MASS::negative.binomial(theta = theta$theta, link = "log")
      model.extra$family = model.extra$model.args$family
      model.extra$class = c("frequentist", "negbin",  "glm")
    } else {
      # baysian neg binomial family
      model.extra$family = rstanarm::neg_binomial_2(link = "log")
    }
  }

  # additions/changes for beta regression
  if(model.type == model.type.list$beta) {
    model.extra$model.run = betareg::betareg
    model.extra$model.args = list(link = "logit", model = F, y = F)
    model.extra$family = make.link(model.extra$model.args$link)
    model.extra$class = c("frequentist", "betareg")
  }

  # additions/changes for ordered regression
  if(model.type == model.type.list$ordered) {
    model.extra$model.run = MASS::polr
    model.extra$model.args = list(method = "logistic")
    model.extra$family = make.link("logit")
    model.extra$class = c("frequentist", "polr")
    model.extra$converged = function(x) x$converged == 1
  }

  # additions/changes for categorical regression
  if(model.type == model.type.list$categorical) {
    model.extra$model.run = nnet::multinom
    model.extra$model.args = list(model = F, trace = F)
    model.extra$family = make.link("logit")
    model.extra$class = c("frequentist", "nnet")
    model.extra$converged = function(x) x$converged == 1

    # Residuals for Multinomial Models: https://www.jstor.org/stable/2673571
    # predictions for MNL https://cran.r-project.org/web/packages/MNLpred/vignettes/OVA_Predictions_For_MNL.html
  }

  # additions/changes for mixed effects
  if(formula.parsed$random.effects & inference == "frequentist") {
    # we cant do some models with random effects
    if(!model.type %in% c(model.type.list$linear, model.type.list$binary)) {
      stop(paste("Unable to run --", model.type, "-- with Random Effects"))
    }

    # model changes
    model.extra$model.run = lme4::glmer
    model.extra$model.args = c(model.extra$model.args, list(control = lme4::lmerControl(calc.derivs = F, optimizer = "nloptwrap")))
    model.extra$class = c(model.extra$class, "merMod")

    # function changes
    model.extra$coefs = lme4::fixef
    model.extra$residuals = lme4:::residuals
  }

  # for survival
  if(model.type == model.type.list$survival) {
    model.extra$model.run = survival::coxph
    model.extra$model.args = list(ties = "efron", y = T, x = T, model = T)
    model.extra$class = c("frequentist", "survival")
    model.extra$fitted = function(x) c("mean_PPD" = mean(x$y[, ncol(x$y)] - x$residuals, na.rm = T))
    model.extra$performance = function(x) c("log-posterior" = x$loglik[2]) # could also use concordance: x$concordance["concordance"]
    model.extra$residuals = function(x) x$residuals # martingale residuals = true events - predict(fit, type = 'expected') events
    model.extra$converged = function(x) x$info["convergence"] == 0

    # fast_bh -- this returns our full baseline hazard instead of a shape parameter as would occur in a stan model -- right now just assume an exponential baseline
    model.extra$special = function(x) { bh = danalyze:::fast_bh(x); intercept = log(sum(bh$hazard) / sum(bh$time)); c("(Intercept)" = intercept) }
  }

  # additional stuff for bayesian models
  if(inference == "bayesian") {
    # set model
    model.extra$model.run = rstanarm::stan_glm

    # set mixed effects version
    if(formula.parsed$random.effects) {
      model.extra$model.run = rstanarm::stan_glmer
    }

    # set class
    model.extra$class[1] = "bayesian"

    # set priors
    model.extra$model.args = NULL
    model.extra$model.args$prior = rstanarm::normal(0, 2.5)
    model.extra$model.args$prior_intercept = rstanarm::normal(0, 2.5)
    model.extra$model.args$prior_aux = rstanarm::exponential(1)

    # set multicore stuff
    model.extra$model.args$cores = 6
    model.extra$model.args$chains = 6
    model.extra$model.args$warmup = 250
    model.extra$model.args$iter = 1000
    model.extra$model.args$seed = 24021985

    # additions/changes for ordered regression
    if(model.type == model.type.list$ordered) {
      model.extra$model.run = rstanarm::stan_polr
      model.extra$model.args$prior = rstanarm::R2(NULL)
      model.extra$model.args$prior_counts = rstanarm::dirichlet(1)
    }

    # additions/changes for survival
    if(model.type == model.type.list$survival) {
      model.extra$model.run = rstanarm::stan_surv
      model.extra$model.args$prior_aux = NULL
      model.extra$model.args$basehaz = "exp"
    }

    ## need to add in brms zero inflated
  }

  # if(model.type == model.type.list$survival && !formula.parsed$random.effects) {
  #   # from terry thernau
  #   # A good rule for Cox models is to have 10-20 events for each coefficient. When models get below 2 events/coef the
  #   # results can be unreliable, both numerically and biologically, and some version of a shrinkage model is called for.
  #
  #   # set model info
  #   model.extra$model.run = survival::coxph
  #   model.extra$model.args = list(ties = "efron", y = T, x = T, model = T)
  #   model.extra$coefs = stats::coef
  #   model.extra$fitted = function(x) c("mean_PPD" = mean(x$y[, ncol(x$y)] - x$residuals, na.rm = T))
  #   model.extra$performance = function(x) c("log-posterior" = x$loglik[2]) # could also use concordance: x$concordance["concordance"]
  #   model.extra$residuals = function(x) x$residuals # martingale residuals = true events - predict(fit, type = 'expected') events
  #   model.extra$converged = function(x) x$info["convergence"] == 0
  #   model.extra$special = function(x) NULL # no shape parameter for the baseline hazard estimated
  #   # model.extra$libraries = "survival"
  # }
  # else if(model.type == model.type.list$survival && formula.parsed$random.effects) {
  #   # # get extras
  #   # .extra = as.list(as.call(model.info$response[[2]]))[-1]
  #   # .args = c("time", "time2", "event", "type", "origin", "cens.code")
  #   # .rep = if(is.null(names(.extra))) 1:length(.extra) else 1:sum(names(.extra) == "")
  #   # names(.extra)[.rep] = .args[!.args %in% names(.extra)][.rep]
  #   #
  #   # # set model info
  #   # model.extra$model.run = survival_wrap
  #   # model.extra$model.args = list(ties = "efron", y = T, x = T, refine.n = 0, .response.var = as.character(.extra$event), .cens.code = .extra$cens.code, .re = T)
  #   # model.extra$clean = clean_glm
  #   # model.extra$coefs = coef_survival
  #   # model.extra$predict = predict_survival
  #   # model.extra$predict.args = NULL
  #   # model.extra$residuals = NULL
  #   # model.extra$libraries = "coxme"
  # }
  # else {
  #   stop(paste("Incorrect model type --", model.info$model, "-- selected."))
  # }

  # add in extra model args
  if(!is.null(model.extra.args)) {
    if(is.null(model.extra$model.args)) {
      model.extra$model.args = model.extra.args
    } else {
      model.extra$model.args[names(model.extra.args)] = model.extra.args
    }
  }

  # return
  return(model.extra)
}
