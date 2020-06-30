# creating a "research plan" by taking some mapping of treatments/controls and using that to create a set of test formulas and new variables

# add transform to a var in a formula
add_transform = function(x, transforms = NULL, lag = NULL) {
  # strip lags
  x.nolag = stringr::str_replace_all(x, paste(paste0("_l", lag), collapse = "|"), "")

  # identify transforms
  x.transform = transforms$transform[match(x.nolag, transforms$var.name)]

  # add the transforms
  r = sapply(1:length(x), function(t) if(!is.na(x.transform[t]) & x.transform[t] != "") paste0(x.transform[t], "(", x[t], ")") else x[t])

  # return
  return(r)
}

# create a new formula object
create_new_formula = function(dv = NULL, iv, controls, type = "*", transforms = NULL, lag = NUL) {
  # the return object
  r = list(
    dv = dv,
    iv = add_transform(iv, transforms),
    iv.transform = type,
    controls = add_transform(controls, transforms)
  )

  # add lengths
  r$dv.length = length(r$dv)
  r$iv.length = length(r$iv)
  r$controls.length = length(controls)

  # set class
  class(r) = "planformula"

  # return
  return(r)
}

#' Print and show methods for planformula
#'
#' @method print planformula
#' @export
#'
print.planformula = function(x, useS4 = F) { print(formula(x)) }

#' @rdname print.planformula
#' @method show planformula
#' @export
#'
show.planformula = function(x) { formula(x) }

#' Formula method for planformula
#'
#' @export
#'
formula.planformula = function(x, to.drop = NULL, transform.iv = NULL) {
  # if we need to drop vars drop them from controls
  if(!is.null(to.drop)) {
    controls = x$controls[!x$controls %in% to.drop]
  } else {
    controls = x$controls
  }

  # do we need to transform the iv
  if(!is.null(transform.iv)) {
    iv = paste0(transform.iv, "(", x$iv, ")")
  } else {
    iv = x$iv
  }

  # create the string
  r = paste("~", paste(c(paste(iv, collapse = x$iv.transform), controls), collapse = " + "))

  # add the dv if we have it
  if(!is.null(x$dv)) r = paste(x$dv, r)

  # turn it into a formula
  as.formula(r, env = .GlobalEnv)
}

#' Terms method for planformula
#'
#' @method terms planformula
#' @export
#'
terms.planformula = function(x) {
  # get formula
  terms(formula(x))
}

#' As.character method for planformula
#'
#' @method as.character planformula
#' @export
#'
as.character.planformula = function(x) {
  # create a the string
  r = paste("~", paste(c(paste(x$iv, collapse = x$iv.transform), x$controls), collapse = " + "))
  if(!is.null(x$dv)) r = paste(x$dv, r)

  # return
  r
}

#' Reformulates the formula so that the IV is the DV.
#'
#' This function is a helper function to allow easy identification of the controls that correlate with the IV.
#' @param x Object of class 'planformula'.
#' @keywords formulas, variables
#' @export
#'
formula_iv_as_dv = function(x) {
  # just loops through and creates a new formula for each iv with the iv as the response and the controls as the explanatory variables
  r = lapply(x$iv, function(y) as.formula(paste(y, "~", paste(x$controls, collapse = " + ")), env = .GlobalEnv))

  # if length is one just pass back
  if(length(r) == 1) r[[1]] else r
}

# this function identifies "other" values for a unit (i.e., the "func" of the variables for each observation of a unit minus the current values for the unit)
build_other_values = function(.data, ..., unit, prefix = "other", func = base::sum, .func.args = NULL, .other = T) {
  # to make this error proof could check to make sure all variables are present and that no variables with "_other" appended exist

  # identify vars that we want to run a function on -- dplyr 1.0 uses the name in summarize_at so remove it
  vars = as.character(tidyselect::vars_select(dplyr::tbl_vars(.data), !!!rlang::enquos(...)))

  # group
  .data = dplyr::group_by_at(.data, .vars = unit)

  # pull function from rlang so we don't need to load the package
  `:=` = rlang:::`:=`

  # summarize
  r = dplyr:::summarize_at(.data, .vars = vars, .funs = rlang::list2(!!prefix := func), !!!.func.args)

  # need to deal with situation where there is only one var
  if(length(vars) == 1) colnames(r)[length(colnames(r))] = paste0(vars, "_", prefix)

  # link back to full data that is "distinct" by unit
  r = dplyr::left_join(dplyr:::select_at(.data, vars), dplyr::ungroup(r), by = unit)

  # subtract out the current value if requested
  if(.other) r[, paste0(vars, "_", prefix)] = dplyr::ungroup(r)[, paste0(vars, "_", prefix)] - dplyr::mutate_at(dplyr::ungroup(r)[, vars], .vars = vars, .funs = function(x) sapply(x, func))

  # select just the relevant variables
  r = dplyr:::select_at(dplyr::ungroup(r), paste0(vars, "_", prefix))

  # rename to a dot-based naming convention
  colnames(r)[match(paste0(vars, "_", prefix), colnames(r))] = paste0(prefix, ".", vars)

  # return
  return(r)
}

#' Identify variables that correlate with our response and modify the formula.
#'
#' This function allows easy identification of the controls that correlate with the IV.
#' @param formula Formula to examine.
#' @param data Data to use.
#' @param cluster Optional data cluster.
#' @param inference Type of inference to use.
#' @param threshold P-value threshold.
#' @param save.re Whether to save the random effects in a formula.
#' @keywords formula, variable, selection
#' @export
#'
trim_formula = function(formula, data, cluster = NULL, inference = "frequentist", threshold = 0.4, save.re = F) {
  # ungroup the data just to be sure
  data = dplyr::ungroup(data)

  # select just the relevant data and remove incomplete cases
  data = dplyr::select(data, all.vars(formula), all.vars(cluster))
  data = dplyr::filter(data, complete.cases(data))

  # run the model
  out.identify = danalyze::analysis(runs = 1000, formula = formula, data = data, cluster = cluster, inference = inference)

  # get the results
  results.identify = danalyze::results(object = out.identify, predictions = NULL, draws = NULL)

  # identify variables in formula
  form.matrix = model.matrix(formula, data)

  # values for current formula
  formula.terms = terms(formula)
  formula.length = 1:length(attr(formula.terms, "term.labels"))

  # identify variables to drop given our threshold -- currently this includes factor components -- need to deal with this
  variables.drop = results.identify$coefficients$coefficient[results.identify$coefficients$p.value >= threshold]

  # identify formula elements to drop
  formula.drop = formula.length[!formula.length %in% (attr(form.matrix, "assign")[!colnames(form.matrix) %in% variables.drop])]

  # create the new formula
  formula.new = formula(drop.terms(formula.terms, dropx = formula.drop, keep.response = T))

  # save random effects if desired
  if(save.re) {
    # determine if we have a mixed effects model
    formula.bars = lme4::findbars(formula)

    # if we do, add back in the RE terms
    if(length(formula.bars) > 0) {
      formula.new = lasso2:::merge.formula(formula.new, formula(paste("~", paste(paste0("(", formula.bars, ")"), collapse = " + "))))
      environment(formula.new) = environment(formula.terms)
    }
  }

  # return
  return(list(trimmed.formula = formula.new, dropped = variables.drop))
}

# internal functions for the research plan

# create variable table
create_variable_table = function(data, treatment, control, interaction, mediation, factor.transform, lag, prefix) {
  # variables that exist and that are created
  # exist: treatment, control, interaction, and mediation (source and target)
  # created: relative, other source, other target, lags

  ## CREATE VARIABLE TABLE

  # the treatment can be a vector or a list -- if a vector than make it an action, if list then also add alternative
  if(is.list(treatment)) {
    treatment.action = if(rlang::has_name(treatment, "action")) treatment$action else treatment[[1]]
    treatment.factor = if(rlang::has_name(treatment, "factor")) treatment$factor else treatment[[2]]
  } else {
    treatment.action = treatment
    treatment.factor = NULL
  }

  # make a table of all variables
  all.vars = list("treatment" = treatment.action, "alternative" = treatment.factor, "control" = control, "interaction" = interaction, "mediation" = mediation)

  # bind
  all.vars = dplyr::bind_rows(lapply(names(all.vars), function(x) if(!is.null(all.vars[[x]])) tibble::tibble(type = x, var.label = names(all.vars[[x]]), var.name = all.vars[[x]])))

  # data var names
  data.var.names = colnames(data)

  # set source, target, relative, source other, target other, and lags if necessary
  all.vars.st = lapply(1:nrow(all.vars), function(x) {
    # get variable name
    var = all.vars$var.name[x]

    # get the base return table
    r = tibble::tibble(
      source = dplyr::if_else(paste(prefix$source, var, sep = ".") %in% data.var.names, paste(prefix$source, var, sep = "."), var), # if there is no source + target it just gets shoved into source
      target = dplyr::if_else(paste(prefix$target, var, sep = ".") %in% data.var.names, paste(prefix$target, var, sep = "."), NA_character_),
      relative = dplyr::if_else(!is.na(source) & !is.na(target), paste(prefix$relative, var, sep = "."), NA_character_), # to make a relative var need both source and target
      source.other = dplyr::if_else(all.vars$type[x] == "treatment", paste(prefix$other, prefix$source, var, sep = "."), NA_character_), # other source only for treatment "action"
      target.other = dplyr::if_else(all.vars$type[x] == "treatment", paste(prefix$other, prefix$target, var, sep = "."), NA_character_) # other target only for treatment "action"
    )

    # set lag names
    if(!is.null(lag) | (is.numeric(lag) && !all(lag == 0))) {
      r = dplyr::mutate(r, dplyr::across(.fns = lapply(lag[lag != 0], function(l) formula(paste0("~ dplyr::if_else(!is.na(.), paste(., 'l", l, "', sep = '.'), NA_character_)"))), .names = "{col}.l{fn}"))
    }

    # return
    return(r)
  })

  # bind
  all.vars.st = dplyr::bind_rows(all.vars.st)

  # combine
  all.vars = dplyr::bind_cols(all.vars, all.vars.st)

  # set the class
  all.vars = dplyr::mutate(all.vars, .factor = sapply(source, function(x) if(!is.na(x) & !is.numeric(data[[x]])) T else F))

  # set default transform
  all.vars$var.transform = sapply(all.vars$.factor, function(x) if(x) function(y) as.numeric(as.factor(y)) else NA)

  # add functions
  if(length(factor.transform) > 0) {
    # only one transform function per variable
    factor.transform = factor.transform[unique(names(factor.transform))]

    # link it
    all.vars$var.transform[which(all.vars$var.name %in% names(factor.transform))] = factor.transform[na.omit(match(all.vars$var.name, names(factor.transform)))]
  }

  # check for missing variables
  all.vars$.missing = is.na(all.vars$source) & is.na(all.vars$target)

  # send a warning if variables are missing
  if(any(all.vars$.missing)) {
    warning(paste("The following variables are missing both source and target names and will be excluded from the plan:", paste(all.vars$var.name[all.vars$.missing], collapse = ", ")))
  }

  # select the columns we need
  all.vars = dplyr::select(dplyr::filter(all.vars, .missing == F), type, var.transform, var.label, var.name, source, target, relative, dplyr::everything(), -.missing, -.factor)


  ## CREATE NAME TABLE

  # create name table
  name.table = sapply(1:nrow(all.vars), function(x) {
    # set the columns to run through
    name.columns = c(if(all.vars$var.name[x] == all.vars$source[x]) "var.name" else c("source", "target"), "relative", "source.other", "target.other")

    # loop through
    lapply(name.columns, function(c) {
      if(!is.na(all.vars[[c]][x])) {
        lapply(c(0, lag), function(l) {
          # set lag
          var.lag = dplyr::if_else(l > 0, paste0(", Lag: ", l, ")"), ")")
          var.lag.name = dplyr::if_else(l > 0, paste0(all.vars[[c]][x], ".l", l), all.vars[[c]][x])

          # set the label
          var.label =
            dplyr::case_when(
              c == "source" ~ paste0(all.vars$var.label[x], " (Source", var.lag),
              c == "target" ~ paste0(all.vars$var.label[x], " (Target", var.lag),
              c == "relative" ~ paste0(all.vars$var.label[x], " (Relative", var.lag),
              c == "source.other" ~ paste0(all.vars$var.label[x], " (Other to Source", var.lag),
              c == "target.other" ~ paste0(all.vars$var.label[x], " (Other to Target", var.lag),
              T ~ paste0(all.vars$var.label[x], if(var.lag == ")") "" else var.lag)
            )

          # return the tibble
          tibble::tibble(type = all.vars$type[x], var.label = var.label, var.name = var.lag.name)
        })
      }
    })
  })

  # bind and make distinct
  name.table = dplyr::distinct(dplyr::filter(dplyr::bind_rows(unlist(name.table, recursive = F)), !is.na(var.name)), type, var.name, .keep_all = T)


  ## RETURN

  # return the variable table
  return(list(variables = all.vars, names = name.table))
}

# select the data
create_plan_data = function(data, all.vars, lag, prefix, unit, keep.vars) {
  # select the data

  # vars to select
  vars.to.select = unique(na.omit(c(keep.vars, unlist(unit), all.vars$source, all.vars$target)))

  # make sure we have the variables
  if(!all(vars.to.select %in% colnames(data))) {
    warning(paste("The following variables are not in the data:", paste(vars.to.select[!vars.to.select %in% colnames(data)], collapse = ", ")))
  }

  # select and subset the relevant data to allow us to easily check variance, etc.
  .data = dplyr::select(dplyr::ungroup(data), tidyselect::any_of(vars.to.select))

  # the extent of what we want to create: lags of treatment, control, interaction and mediation; other values of treatment against source/target

  # arrange
  .data = dplyr::arrange_at(.tbl = .data, .vars = c(unit$unit, unit$source, unit$target, unit$time))

  #  create other values and add to our data frame -- done first so that we can create lags later -- need to deal with FACTOR VARS
  if(!is.null(unit$source)) {
    other.source.vars = build_other_values(.data, !!all.vars$var.name[all.vars$type == "treatment"], unit = c(unit$unit, unit$source, unit$time), prefix = paste(prefix$other, prefix$source, sep = "."))
    .data = dplyr::bind_cols(.data, other.source.vars)
  }

  # do also for target if present
  if(!is.null(unit$target)) {
    other.target.vars = build_other_values(.data, !!all.vars$var.name[all.vars$type == "treatment"], unit = c(unit$unit, unit$target, unit$time), prefix = paste(prefix$other, prefix$target, sep = "."))
    .data = dplyr::bind_cols(.data, other.target.vars)
  }

  # create relative variables and lags
  for(i in 1:nrow(all.vars)) {
    # we have a source and target so we can make a relative variable
    if(!is.na(all.vars$relative[i])) {
      # set source and target
      t.source = .data[[all.vars$source[i]]]
      t.target = .data[[all.vars$target[i]]]

      # if not numeric then process
      if(!is.na(all.vars$var.transform[i])) {
        t.source = all.vars$var.transform[[i]](t.source)
        t.target = all.vars$var.transform[[i]](t.target)
      }

      # set relative var and if both are zero set relative var to zero -- could create relative variables in this way or in the formula: ~ source + target + I(source / (source + target))
      .data[[paste("rel", all.vars$var.name[i], sep = ".")]] = dplyr::if_else(t.source == 0 & t.target == 0, 0, t.source / (t.source + t.target))
    }

    # we need to make lags so make lags
    if(!is.null(lag) | (is.numeric(lag) && !all(lag == 0))) {
      # variables to create lags for
      col.other.names = na.omit(c(all.vars$source[i], all.vars$target[i], all.vars$relative[i], all.vars$source.other[i], all.vars$target.other[i]))

      # create the lags
      .data = dplyr::mutate(.data, dplyr::across(.cols = tidyselect::any_of(col.other.names), .fns = lapply(lag[lag != 0], function(l) formula(paste0("~ lag(., n = ", l, ")"))), .names = "{col}.l{fn}"))

      # zeros should be okay instead of NA for actions but not for control or factors
    }
  }

  # ungroup data just to make sure
  .data = dplyr::ungroup(.data)

  # return
  return(.data)
}

# function to identify variables that should not be included in the analysis
identify_problematic_variables = function(.data, unit.main, variables.to.check) {
  # select -- ordering so we keep source and target as highest priority
  .data.col = dplyr::select(.data, tidyselect::any_of(c(unit.main, variables.to.check)))

  # remove entirely problematic variables
  bad.var = colnames(.data.col)[unlist(dplyr::summarise_all(.data.col, function(x) any(is.infinite(x)) | any(is.nan(x))))]

  # first find variables with little variance -- exclude unit and treatment from check
  no.var = caret::nearZeroVar(.data.col[!colnames(.data.col) %in% unit.main], freqCut = 95/5, uniqueCut = 0, foreach = T, names = T)

  # find variables that dont vary within our groups -- exclude treatment (unit is part of group_by)
  .data.group = dplyr::group_by_at(.tbl = .data.col, .vars = unit.main)
  group.variation = dplyr::group_map(.data.group, ~ caret::nearZeroVar(.x, freqCut = 99/1, uniqueCut = 0, foreach = T, names = T)) # less stringent criteria
  group.variation = 1 - unlist(as.list(table(unlist(group.variation)))) / dplyr::n_groups(.data.group) # what portion of groups have variation
  no.group.var = names(group.variation)[group.variation < 0.01]

  # now identify collinear vars -- drop factors could eventually look for linear combinations for factors -- also drops no variation vars since we already exclude those
  data.for.cor = dplyr::select_if(.data.col[!colnames(.data.col) %in% c(unit.main, no.var, bad.var)], is.numeric)
  col.var = caret::findCorrelation(cor(data.for.cor, use = "pairwise.complete"), cutoff = 0.95, names = T, exact = T)

  # set final no var -- these are good candidates to drop from our analysis since there is so little variation -- tell the user though
  all.prob.vars = unique(c(no.var, no.group.var, col.var))

  # we could also try to find linear combinations but probably not needed at this step -- can also just rely on priors to address this

  # return
  return(all.prob.vars)
}

# function to identify variable transformations
identify_variable_transformations = function(.data, variables.to.transform) {
  # helpful: https://www.ibm.com/support/pages/transforming-variable-normality-parametric-statistics

  # need to take care of source, target, other source, other target, and source/target/other lags
  all.transform = tibble::tibble(var.name = variables.to.transform)

  # dont transform factor variables so just drop them for now
  all.transform = all.transform[sapply(all.transform$var.name, function(x) is.numeric(.data[[x]])), ]

  # set sign
  all.transform$sign = sapply(all.transform$var.name, function(x) {
    dat = na.omit(.data[[x]])
    dplyr::case_when(dplyr::n_distinct(dat) == 2 ~ "dummy", all(dat >= 0) & all(dat <= 0) ~ "portion", all(dat <= 0) ~ "negative", all(dat >= 0) ~ "positive", T ~ "both")
  })

  # identify portion min or max
  all.transform$portion.min = sapply(all.transform$var.name, function(x) mean(.data[[x]] == min(.data[[x]], na.rm = T), na.rm = T))
  all.transform$portion.max = sapply(all.transform$var.name, function(x) mean(.data[[x]] == max(.data[[x]], na.rm = T), na.rm = T))

  # identify skewness: < 1 = normal to half-normal; < 2 = exponential; > 2 = lognormal
  all.transform$skewness = sapply(all.transform$var.name, function(x) e1071::skewness(.data[[x]], na.rm = T, type = 2))

  # identify kurtosis: < 0 = fat tails; > 0 = thin tails; ~0 = normal tails
  all.transform$kurtosis = sapply(all.transform$var.name, function(x) e1071::kurtosis(.data[[x]], na.rm = T, type = 2))

  ## TODO:
  # when testing mediation analysis its possible for the mediation equation to produce negative values of an otherwise positive variable
  # variables on very different scales

  # set the transform
  all.transform$transform =
    dplyr::case_when(
      # all.transform$skewness > 1.5 & all.transform$portion.min < 0.5 & all.transform$sign == "positive" ~ "log1p",
      # all.transform$skewness > 1.5 & all.transform$sign == "positive" ~ "sqrt",
      all.transform$sign == "dummy" ~ "",
      all.transform$sign == "portion" ~ "asin",
      # all.transform$skewness > 1.5 & all.transform$portion.min < 0.5 ~ "danalyze::symsqrt",
      all.transform$skewness > 1.5 ~ "danalyze::symlog",
      T ~ "scale"
    )

  # return
  return(all.transform)
}

# function to create formulas to test for the research plan
create_research_formulas = function(all.vars, all.transform, lag) {
  # idea is to create a formula: each treatment, each treatment interacted with a condition, each treatment as a cause of a mediator (+ the treatment and the mediator)

  # ultimately we will want to identify variables that don't vary enough over time, within a unit, or between units

  # treatments to loop through
  main.ivs = all.vars$source[all.vars$type %in% c("treatment", "alternative")]

  # loop through each treatment
  all.formulas = lapply(main.ivs, function(t) {
    # a lag of "0" means just the normal variable
    if(is.null(lag)) {
      lag = 0
    }

    # columns we want to include
    t.vars = lapply(lag, function(l) {
      # column select
      c.extra = dplyr::if_else(l == 0, "", paste0(".l", l))

      # return var names
      c(all.vars[all.vars$type == "control", paste0("source", c.extra)],
        all.vars[all.vars$type == "control", paste0("target", c.extra)],
        all.vars[all.vars$type == "control", paste0("relative", c.extra)],
        all.vars[, paste0("source.other", c.extra)],
        all.vars[, paste0("target.other", c.extra)])
    })

    # get character vectoe
    t.vars = as.character(na.omit(unlist(t.vars)))

    # pull from rlang
    `:=` = rlang:::`:=`

    # for each formula block we need: (1) list of treatment + interaction + mediator -- tried this using an r object

    # old: rlang::list2(!!t := as.formula(paste("~", paste(add_transform(c(t, t.vars, t.lags)), collapse = " + ")), env = .GlobalEnv))

    # combine to create the main formula
    f.base = create_new_formula(iv = t, controls = t.vars[!t.vars %in% c(t)], transforms = all.transform, lag = lag)

    # create formulas for interactions
    f.int.vars = na.omit(c(all.vars$source[all.vars$type == "interaction"], all.vars$target[all.vars$type == "interaction"], all.vars$relative[all.vars$type == "interaction"]))
    f.int = lapply(f.int.vars, function(i) create_new_formula(iv = c(t, i), controls = t.vars[!t.vars %in% c(t, i)], transforms = all.transform, lag = lag))
    names(f.int) = f.int.vars

    # create formulas for mediating effects -- requires two "med ~ iv" and "outcome ~ iv + med"
    f.med.vars = dplyr::if_else(sapply(all.vars$target[all.vars$type == "mediation"], is.null), all.vars$source[all.vars$type == "mediation"], all.vars$target[all.vars$type == "mediation"])
    f.med = lapply(f.med.vars, function(m) list(mediator = create_new_formula(dv = m, iv = t, controls = t.vars[!t.vars %in% c(t, m)], transforms = all.transform, lag = lag),
                                                outcome = create_new_formula(iv = c(t, m), type = "+", controls = t.vars[!t.vars %in% c(t, m)], transforms = all.transform, lag = lag)))
    names(f.med) = f.med.vars

    # create full return
    f.full = list(main = f.base, interaction = f.int, mediation = f.med)

    # return
    return(f.full) #rlang::list2(!!t := f.full))
  })

  # set names
  names(all.formulas) = all.vars$source[all.vars$type %in% c("treatment", "alternative")]

  # return
  return(all.formulas)
}

# create a research plan for the given variables

#' Function to get create new variables and formulas for testing a set of variable relationships.
#'
#' This function produces a set of formulas and new data for testing a research plan.
#' @param treatment Variable containing names of treatments to test. Can be a vector for just actions or a list that contains both 'action' and 'factor' to test.
#' @param control Optional vector containing variable names for controls to include.
#' @param interaction Optional vector containing interactions for the treatment.
#' @param mediation Optional vector containing mediating variables for the treatment.
#' @param data Data containing variables in treatment, control, interaction, and mediation.
#' @param keep.vars Optional vector of variable names to keep in the data frame returned.
#' @param prefix List of prefixes for source and target in data as well as for "other" treatments used against source/target. Source and target are optional.
#' @param unit List of variables that correspond to units, times, and optional source/target in data.
#' @param lag Optional vector containing the number of lags to use.
#' @keywords formulas, variables
#' @export
#' @examples
#' research_plan(treatment = c("use.force", "show.force"),
#' control = c("wdi.economy", "wdi.terrain"), data = data,
#' unit = list(unit = "dispute.id", time = "time"), lag = 1:2)
#'

research_plan = function(treatment, control = NULL, interaction = NULL, mediation = NULL, data, keep.vars = NULL, factor.transform = NULL,
                         prefix = list(source = "sv", target = "tv", other = "oth", relative = "rel"), unit = list(unit = NULL, source = NULL, target = NULL, time = NULL), lag = NULL) {
  # below, we identify the variables we need to use, create new data to correspond to these variables, create formulas to use the variables, and then return everything

  # check to make sure arguments are good
  if(!(is.list(treatment) | is.character(treatment)) | !is.data.frame(data)) {
    stop("To create a research plan it is necessary to have a treatment and data.")
  }


  ## CREATE VARIABLE TABLE

  # TODO: could add a consecutive transformation (past or within X * lag)

  # create a table with the variables
  variable.table = create_variable_table(data = data, treatment = treatment, control = control, interaction = interaction, mediation = mediation, factor.transform = factor.transform, lag = lag, prefix = prefix)

  # set all vars
  all.vars = variable.table$variables


  ## SELECT DATA AND CREATE NEW VARIABLES

  # select the data and create new variables
  .data = create_plan_data(data = data, all.vars = variable.table$variables, lag = lag, prefix = prefix, unit = unit, keep.vars = keep.vars)


  ## IDENTIFY PROBLEMATIC VARIABLES

  # list of variables to check -- all variables except the main treatment action/alternative variables
  variables.to.check = na.omit(unlist(dplyr::select(variable.table$variables, dplyr::everything(), -type, -var.transform, -var.label, -var.name)))
  variables.to.check = as.character(variables.to.check[!variables.to.check %in% c(variable.table$variables$source[variable.table$variables$type %in% c("treatment", "alternative")],
                                                                                  variable.table$variables$target[variable.table$variables$type %in% c("treatment", "alternative")])])

  # run the function to identify problematic variables
  all.prob.vars = identify_problematic_variables(.data = .data, unit.main = unit$unit, variables.to.check = variables.to.check)

  # remove vars that are problematic and drop vars that no longer have any values to test
  variable.table$variables = dplyr::mutate_all(variable.table$variables, .funs = function(x) { x[x %in% all.prob.vars] = NA_character_; x })
  variable.table$variables = dplyr::filter(variable.table$variables, !is.na(source) | !is.na(target) | !is.na(relative) | !is.na(source.other) | !is.na(target.other))


  ## IDENTIFY VARIABLE TRANSFORMATIONS

  # set variables that we want to identify transformations for
  variables.to.transform = unique(na.omit(unlist(dplyr::select(variable.table$variables, -c(type:var.name)))))

  # identify variable transformations
  all.transform = identify_variable_transformations(.data = .data, variables.to.transform = variables.to.transform)


  ## CREATE FORMULAS FOR TESTING

  # create formulas
  all.formulas = create_research_formulas(all.vars = variable.table$variables, all.transform = all.transform, lag = lag)

  # the returned data frame is missing the outcome and any fixed/random effects (i.e., all variables that are not IVs in some way) -- need to think about how to deal with this


  ## RETURN

  # return all
  return(list(formulas = all.formulas, data = .data, variables = list(table = variable.table$variables, transform = all.transform, label = variable.table$names)))
}

