# creating a "research plan" by taking some mapping of treatments/controls and using that to create a set of test formulas and new variables

#' Update a formula with provided transform functions.
#'
#' This internal function takes a formula and modifies it to include the
#' variable transforms provided as an input. This function is used in the 'create_new_formula'
#' function.
#'
#' @family research plan functions
#' @param x The variables to add transformations to.
#' @param transforms A dataframe identifying variables and transforms to use for each variable.
#' @param lag The lag level to remove from the formula.
#' @return An updated formula with appropriate variable transformation functions added.
#'
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

#' Create a new 'planformula' object.
#'
#' This internal function takes components of a formula and creates a new formula object
#' that is used by other research plan functions. The 'planformula' object also has its
#' own print and update functions.
#'
#' @family research plan functions
#' @param dv A character string naming the dependent variable.
#' @param iv A vector of character strings naming the main independent variables.
#' @param controls A vector of character strings naming the control variables.
#' @param type A character indicating the relationship between the 'iv's.
#' @param transforms A list providing variable names and transformation for those variables.
#' @param lag Lags to include (currently unused).
#' @return An object that includes all the components (including variable transformations) for a
#' testable formula. This formula object is passed to other functions to analyze.
#'
create_new_formula = function(dv = NULL, iv, controls, type = "*", transforms = NULL, lag = NULL) {
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
#' These functions provide easy ways to manpulate or show the contents of a
#' 'planformula' object.
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

#' @rdname print.planformula
#' @method formula planformula
#' @export
#' @param x A return from the function 'create_new_formula'.
#' @param to.drop A list of variable names to remove from the formula.
#' @param transform.iv A string that identifies a function name to wrap the main IV.
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

#' @rdname print.planformula
#' @method terms planformula
#' @export
#'
terms.planformula = function(x) {
  # get formula
  terms(formula(x))
}

#' @rdname print.planformula
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
#'
#' @family research plan functions
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

#' Automates creation of variables in formula.
#'
#' This function is a helper function to automatically create variables (in both a formula and the data)
#' that identify other "actions" taken against a target in a unit of analysis. Specifically, it identifies
#' the "func" (defaults to "sum") of the variables for each observation of a unit minus the current values
#' for that unit.
#'
#' The goal is to make this process as generalizable as possible so these functions are usable in a variety
#' of contexts.
#'
#' @family research plan functions
#' @param .data The dataframe containing all variables.
#' @param ... Variables for which to create "other" values. Should be entered a symbols and not strings.
#' @param unit A vector listing the variables that describe the unit of analysis in the data.
#' @param prefix A string indicating the string to be added as a prefix to the variable name to indicate that
#'   it is an "other" variable.
#' @param func The function applied to the vector of valuables in each unit of analysis.
#' @param .func.args Arguments to pass to 'func'.
#' @param .other Whether or not to subtract observations in each row from the "other" observations in each unit.
#'   The default is to do so, i.e., if you want to sum the number of other attacks against a target the number of
#'   attacks initiated by the source in the current observation is subtracted.
#' @return A dataframe with the other variables created. The column names are the names of the variables. The number of
#'   of rows should match the numbr of rows in'.data'.
#' @keywords formulas, variables
#' @export
#'
build_other_values = function(.data, ..., unit, prefix = "other", func = base::sum, .func.args = NULL, .other = T) {
  # this function identifies "other" values for a unit (i.e., the "func" of the variables for each observation of a unit minus the current values for the unit)

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
#' This function allows easy identification of the control variables that correlate with the IV. This function
#' can be used in conjunction with 'formula_iv_as_dv' to simplify the process of taking a formula, identifying
#' variables in the formula that correlate with the main independent variable(s) and trimming the formula to
#' only include these variables.
#'
#' @family research plan functions
#' @param formula Formula to examine.
#' @param data Data to use.
#' @param cluster Optional data cluster.
#' @param inference Type of inference to use.
#' @param threshold P-value threshold for inclusion or exclusion of control variables. Variables with p-values above
#'   this threshold will be dropped.
#' @param save.re Whether to save the random effects in the returned formula.
#' @return Returns a list with two values. The first is a new formula that only includes the variables that correlate
#'   with the dependent variable. The second is a list of the variables that were dropped from the formula.
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
  out.identify = analysis(runs = 1000, formula = formula, data = data, cluster = cluster, inference = inference)

  # get the results
  results.identify = results(object = out.identify, predictions = NULL, draws = NULL)

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
  if(length(formula.drop) > 0) {
    formula.new = formula(drop.terms(formula.terms, dropx = formula.drop, keep.response = T))
  } else {
    formula.new = formula(formula.terms)
  }


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

#' Create all the new variables to include in each formula.
#'
#' This function creates a large number of new variables to potentially include in a formula. These variables
#' include source values, target values, relative values, and other values.
#'
#' For included variable lists (treatment, control, interaction, and mediation) the name provides a human-understandable
#' label while the value gives the name as found in 'data'.
#'
#' @family research plan functions
#' @param data Data to use.
#' @param treatment The list of treatment variables to check. May have two names in the list: action and factor.
#'   An 'action' is a variable that can have source/target/other transforms. A 'factor' is not transformed.
#' @param control The list of control variables.
#' @param interaction The list of interaction variables.
#' @param mediation The list of mediation variables.
#' @param factor.transform A provided function that allows a factor/character variable to be transformed.
#' @param lag The number of lags to include.
#' @param unit A list identifying the unit of analysis in the data. May be NULL if no unit is present.
#' @param prefix A list that provides the prefixes for source, target, relative, and other.
#' @return A list with two values. The first value is a dataframe that lists all of the variables created and their
#'   transformations. The second provides labels for all included variables.
#' @keywords formula, variable, selection
#'
create_variable_table = function(data, treatment, control, interaction, mediation, factor.transform, lag, unit, prefix) {
  # variables that exist and that are created
  # exist: treatment, control, interaction, and mediation (source and target)
  # created: relative, other source, other target, lags

  ## CREATE VARIABLE TABLE

  # the treatment can be a vector or a list -- if a vector than make it an action, if list then also add alternative
  if(is.list(treatment)) {
    treatment.action = if(rlang::has_name(treatment, "action")) treatment$action else treatment$action
    treatment.factor = if(rlang::has_name(treatment, "factor")) treatment$factor else treatment$factor
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

    # we only want to create "other" varaibles if we have a unit

    # get the base return table
    r = tibble::tibble(
      source = dplyr::if_else(paste(prefix$source, var, sep = ".") %in% data.var.names, paste(prefix$source, var, sep = "."), var), # if there is no source + target it just gets shoved into source
      target = dplyr::if_else(paste(prefix$target, var, sep = ".") %in% data.var.names, paste(prefix$target, var, sep = "."), NA_character_),
      relative = dplyr::if_else(!is.na(source) & !is.na(target), paste(prefix$relative, var, sep = "."), NA_character_), # to make a relative var need both source and target
      source.other = dplyr::if_else(all.vars$type[x] == "treatment" & !is.null(unit$unit), paste(prefix$other, prefix$source, var, sep = "."), NA_character_), # other source only for treatment "action"
      target.other = dplyr::if_else(all.vars$type[x] == "treatment" & !is.null(unit$unit), paste(prefix$other, prefix$target, var, sep = "."), NA_character_) # other target only for treatment "action"
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

#' Create all data that corresponds to newly created variables.
#'
#' Given variables specified in 'create_variable_table' this function creates the new corresponding variables
#' in the data.
#'
#' @family research plan functions
#' @param data Data to use.
#' @param all.vars A list of new variables to create. This is the first list item returned from 'create_variable_table'.
#' @param lag The number of lags to include.
#' @param prefix A list that provides the prefixes for source, target, relative, and other.
#' @param unit A vector of strings identifying the name of the unit of analysis.
#' @param keep.vars Additional variables to keep in 'data' when returned.
#' @return A dataframe that includes all of the new variables created. This dataframe will also include the variables
#'   identified by 'keep.vars'.
#' @keywords formula, variable, selection
#'
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

#' Exclude variables with insufficient variation to allow inclusion.
#'
#' This function takes as input a dataframe, a unit of analysis, and a list of variables and identifies which variables have
#' sufficient variation to allow inclusion in a regression function. Variables that have insufficient variation in the data,
#' insufficient variation on average in each unit of analysis, or that are collinear are flagged for exclusion.
#'
#' This processing step is only applied to numeric variables.
#'
#' @family research plan functions
#' @param .data Data to use.
#' @param unit.main A list of variables that identifies the main unit within the data. This is only applicable if there are
#'   multiple observations for each unit (e.g., observations over time for a single actor).
#' @param variables.to.check The list of variables to check for sufficient variation.
#' @return This function returns a list of variable names that are problematic to include in regression analysis. These are
#'   variables that have no variation, no variation within-group variaiton, or that are collinear.
#' @keywords formula, variable, selection
#'
identify_problematic_variables = function(.data, unit.main, variables.to.check) {
  # function to identify variables that should not be included in the analysis

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

#' Based on the data, identify how to transform a variable.
#'
#' This function takes as input a dataframe and a list of variables and identifies what transformations
#' should be applied to each variable to make it more normal.
#'
#' This processing step is only applied to numeric variables.
#'
#' @family research plan functions
#' @param .data Data to use.
#' @param variables.to.transform The list of variables to identify the appropriate transformation for.
#' @return This function returns a dataframe where each row is a variable and each column characterizes the data of that variable
#'   and the appropriate transformation function to use. Columns include description of the sign, the portion of minimum and maximum values,
#'   the skewness, and the kurtosis. Based on these characteristics, an appropriate transformation function for the variable is suggested.
#' @keywords formula, variable, transformation
#'
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
      all.transform$skewness > 1.25 ~ "danalyze::symlog",
      T ~ "scale"
    )

  # return
  return(all.transform)
}

#' Create formulas for testing.
#'
#' Given suggested variables to include and how to transform them, this function creates the actual formulas that
#' can be tested in a regression model.
#'
#' For each variable identified as a treatment this function creates a formula for the main effect, a formula for each
#' interaction, and formulas for each possible mediator. For each formula, a prediction to test the formula is also created.
#'
#' @family research plan functions
#' @param all.vars The return value from the function 'create_variable_table'.
#' @param all.transform The return value from the function 'identify_variable_transformations'.
#' @param lag The number of lags for each control variable to include.
#' @return A list of formulas and predictions that allow main effects, interaction effects, and mediating effects to be
#'   tested.
#' @keywords formula, variable, transformation
#'
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

    # check if we have variables
    if(length(all.vars$target[all.vars$type == "interaction"]) > 0) {
      # create formulas for interactions
      f.int.vars = na.omit(c(all.vars$source[all.vars$type == "interaction"], all.vars$target[all.vars$type == "interaction"], all.vars$relative[all.vars$type == "interaction"]))
      f.int = lapply(f.int.vars, function(i) create_new_formula(iv = c(t, i), controls = t.vars[!t.vars %in% c(t, i)], transforms = all.transform, lag = lag))
      names(f.int) = f.int.vars
    } else {
      f.int = NULL
    }

    # check if we have variables
    if(length(all.vars$target[all.vars$type == "mediation"]) > 0) {
      # create formulas for mediating effects -- requires two "med ~ iv" and "outcome ~ iv + med"
      f.med.vars = dplyr::if_else(sapply(all.vars$target[all.vars$type == "mediation"], is.null), all.vars$source[all.vars$type == "mediation"], all.vars$target[all.vars$type == "mediation"])
      f.med = lapply(f.med.vars, function(m) list(mediator = create_new_formula(dv = m, iv = t, controls = t.vars[!t.vars %in% c(t, m)], transforms = all.transform, lag = lag),
                                                  outcome = create_new_formula(iv = c(t, m), type = "+", controls = t.vars[!t.vars %in% c(t, m)], transforms = all.transform, lag = lag)))
      names(f.med) = f.med.vars
    } else {
      f.med = NULL
    }

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

#' Function to get create new variables and formulas for testing a set of variable relationships.
#'
#' This function produces a set of formulas and new data for testing a research plan.
#'
#' @family research plan functions
#' @param treatment Variable containing names of treatments to test. Can be a vector for just actions or a list that contains both 'action' and 'factor' to test.
#' @param control Optional vector containing variable names for controls to include.
#' @param interaction Optional vector containing interactions for the treatment.
#' @param mediation Optional vector containing mediating variables for the treatment.
#' @param data Data containing variables in treatment, control, interaction, and mediation.
#' @param keep.vars Optional vector of variable names to keep in the data frame returned.
#' @param prefix List of prefixes for source and target in data as well as for "other" treatments used against source/target. Source and target are optional.
#' @param unit List of variables that correspond to units, times, and optional source/target in data.
#' @param lag Optional vector containing the number of lags to use.
#' @return A research plan object that contains all the formulas and data to allow automated analysis.
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

  # update to incorporate selection based on dependent variable

  # should suggest a research design (differenced if unit has pre-post, etc.)

  # check to make sure arguments are good
  if(!(is.list(treatment) | is.character(treatment)) | !is.data.frame(data)) {
    stop("To create a research plan it is necessary to have a treatment and data.")
  }


  ## CREATE VARIABLE TABLE

  # TODO: could add a consecutive transformation (past or within X * lag)

  # create a table with the variables
  variable.table = create_variable_table(data = data, treatment = treatment, control = control, interaction = interaction, mediation = mediation, factor.transform = factor.transform, lag = lag, unit = unit, prefix = prefix)

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

#' Function to automatically analyze a research plan in a systematic and efficient manner.
#'
#' Takes a research plan created from the function 'research_plan' and performs all necessary analysis. This analysis is then
#' returned to the user to allow various types of assessment.
#'
#' @family research plan functions
#' @param research.plan A research plan to analyze. This should be the return value from the function 'research_plan'.
#' @param .threshold The p-value used to identify which variables correlate with the treatment and thus which potential
#'   confounding variables are included in the analysis.
#' @param only.run The name of a variable in the research plan to analyze. If this is provided only the named variable will
#'   be analyzed and returned,
#' @param fill.in.plan An optional results from a previous 'analyze_plan' run. If provided, the function will ignore values
#'   already present and just run variables, interactions, or mediators that are not present. Using this option is an easy way
#'   to add a new variable to a plan and run just this new variable.
#' @param run.basic Logical indicator of whether to run the main effects for a variable.
#' @param run.mediation Logical indicator of whether to run the interaction effects for a variable.
#' @param run.interaction Logical indicator of whether to run the mediation effects for a variable.
#' @export
#' @return The statistical results of the analysis. Each list item corresponds to a variable in the research plan. Within each
#'   variable separate dataframes contain the results for the main analysis, the interaction effects, and the mediation analysis.
#'   If there are multiple outcomes examined these are added to the dataframe. If time is a component (e.g., when conducting
#'   survival anlysis) time is also added to the dataframe.
#'
#'   Mediation analysis is only performed if the main result for a variable/outcome is significant.
#' @keywords automated, analysis, research plan
#'
analyze_plan = function(research.plan, .threshold = 0.3, only.run = NULL, fill.in.plan = NULL, run.basic = T, run.mediation = T, run.interaction = T) {
  # first check to make sure research.plan is correctly formatted
  if(!is.list(research.plan) | !rlang::has_name(research.plan, "formulas") | !rlang::has_name(research.plan, "data") | !rlang::has_name(research.plan, "variables")) {
    stop("The research plan does not have the correct fields please check to make sure the variable is correct.")
  }

  # get list of main variables to run
  main.variables = names(research.plan$formulas)

  # subset of main variables to run when analyzing the plan
  if(!is.null(only.run)) {
    main.variables = main.variables[main.variables %in% only.run]
  }

  # make sure we still have variables to run
  if(is.null(main.variables) || length(main.variables) < 1) {
    stop("No variables in research pla nto run analysis on.")
  }

  # console output
  cat(paste0("\n\n** Running research plan with -- ", length(main.variables), " -- variable(s).\n"))

  # loop structure is outcome --> main variable (basic, interaction, mediation)

  # loop through variables
  all.variables = lapply(main.variables, function(variable) {
    # for each outcome, loop through the variables and run the main, interaction, and mediator

    # set data
    temp.data = research.plan$data

    # console output
    cat(paste0("\n** Identifying correlates of -- ", variable, ".\n"))

    # add random/fixed effects if needed
    if(!is.null(research.plan$effects)) {
      formula.effects = lasso2:::merge.formula(formula_iv_as_dv(research.plan$formulas[[variable]]$main), research.plan$effects)
    } else {
      formula.effects = formula_iv_as_dv(research.plan$formulas[[variable]]$main)
    }

    # get the list of variables to drop for our main effect
    main.drop = trim_formula(formula = formula.effects, data = temp.data, cluster = research.plan$cluster, threshold = .threshold)

    # loop through variables
    r.outcome = lapply(names(research.plan$outcome.occurrence), function(dv) {
      # console output
      cat(paste0("\n** Working on outcome -- ", dv, ".\n"))

      # set outcome -- needs to be integrated with setting the outcome in the plan -- right now a little artificial
      temp.data$.outcome = research.plan$outcome.occurrence[[dv]]

      # run main models if needed
      if(run.basic) {
        # check to see if we already have information or if we need to fill it in
        if(is.list(fill.in.plan) && rlang::has_name(fill.in.plan[[variable]], "main.contrasts") &&
           dv %in% fill.in.plan[[variable]]$main.contrasts$.outcome) {
          # we already have information in a plan that we are updating so tell the user and skip the additional analysis
          cat(paste0("\n** Already have main model results for -- ", variable, ".\n"))

          # save recycle our existing results
          r.main = list(list(variable = variable, model = NULL,
                             coefficients = dplyr::filter(fill.in.plan[[variable]]$main.coefficients, .outcome == dv),
                             predictions = dplyr::filter(fill.in.plan[[variable]]$main.predictions, .outcome == dv),
                             contrasts = dplyr::filter(fill.in.plan[[variable]]$main.contrasts, .outcome == dv)))
        } else {
          # console output
          cat(paste0("\n** Running main model -- ", variable, ".\n"))

          # add the outcome to our formula so we now have a testable formula -- can suggest one component of a factor which of course is problematic
          new.formula = update(formula(research.plan$formulas[[variable]]$main, to.drop = main.drop$dropped), research.plan$outcome)

          # merge random effects if present
          if(!is.null(research.plan$effects)) {
            new.formula = lasso2:::merge.formula(new.formula, research.plan$effects)
          }

          # analyze our formula and data
          model = NULL
          try(model <- analysis(
            runs = 500,
            formula = new.formula,
            main.ivs = variable,
            data = temp.data,
            inference = "bayesian",
            model.extra.args = list(prior = rstanarm::normal(0, 1), prior_intercept = rstanarm::normal(0, 1), adapt_delta = 0.9, warmup = 250, iter = 500)
          ), T)

          if(!is.null(model)) {
            # create prediction for treatment
            prediction = pr_list(!!variable := create_values(temp.data[[variable]]))

            # set times -- this needs to be fixed so that it is not hard coded
            times = unique(quantile(temp.data$.time.end, seq(1, 9, by = 2) / 10, na.rm = T))

            # get results
            results.m = results(object = model, predictions = prediction, draws = 500, times = times)

            # save -- double list so it is easier to use purrr
            r.main = list(list(variable = variable, model = model, coefficients = results.m$coefficients, predictions = results.m$predictions, contrasts = results.m$contrasts))
          } else {
            r.main = list(list(variable = variable, model = NULL, coefficients = NULL, predictions = NULL, contrasts = NULL))
          }
        }

        # remove nulls from results and set names
        r.main[sapply(r.main, is.null)] = NULL
        if(!is.null(r.main)) { n = purrr::map_chr(r.main, "variable"); if(!is.null(n)) names(r.main) = n }
      } else {
        r.main = NULL
      }

      # get list of interaction variables to run
      interaction.variables = names(research.plan$formulas[[variable]]$interaction)

      # now run the interaction
      if(run.interaction & !is.null(interaction.variables)) {
        # loop through interaction variables
        r.interaction = lapply(interaction.variables, function(interaction) {
          # make sure the interaction is not the same as the main variable
          if(interaction == variable) {
            return(NULL)
          }

          # check to see if we already have information or if we need to fill it in
          if(is.list(fill.in.plan) && rlang::has_name(fill.in.plan[[variable]], "interaction.contrasts") &&
             dv %in% fill.in.plan[[variable]]$interaction.contrasts$.outcome && interaction %in% fill.in.plan[[variable]]$interaction.contrasts$.main.interaction) {
            # we already have information in a plan that we are updating so tell the user and skip the additional analysis
            cat(paste0("\n** Already have interaction results for -- ", variable, " X ", interaction, ".\n"))

            # save recycle our existing results
            r.inte = list(variable = variable, interaction = interaction, model = NULL,
                          coefficients = dplyr::filter(fill.in.plan[[variable]]$interaction.coefficients, .outcome == dv, .main.interaction == interaction),
                          predictions = dplyr::filter(fill.in.plan[[variable]]$interaction.predictions, .outcome == dv, .main.interaction == interaction),
                          contrasts = dplyr::filter(fill.in.plan[[variable]]$interaction.contrasts, .outcome == dv, .main.interaction == interaction))
          } else {
            # console output
            cat(paste0("\n** Running interaction -- ", variable, " X ", interaction, ".\n"))

            # trim the formula -- now we just have variables that correlate with our treatment -- should also make it work for multiple variables

            # add the outcome to our formula so we now have a testable formula
            new.formula = update(formula(research.plan$formulas[[variable]]$interaction[[interaction]], to.drop = main.drop$dropped), research.plan$outcome)

            # merge random effects if present
            if(!is.null(research.plan$effects)) {
              new.formula = lasso2:::merge.formula(new.formula, research.plan$effects)
            }

            # analyze our formula and data
            model.i = NULL
            try(model.i <- analysis(
              runs = 500,
              formula = new.formula,
              main.ivs = variable,
              data = temp.data,
              inference = "bayesian",
              model.extra.args = list(prior = rstanarm::normal(0, 1), prior_intercept = rstanarm::normal(0, 1), adapt_delta = 0.9, warmup = 250, iter = 500)
            ), T)

            # save if possible
            if(!is.null(model.i)) {
              # create prediction for treatment
              prediction = pr_list(!!variable := create_values(temp.data[[variable]]), !!interaction := create_values(temp.data[[interaction]]), .constant = interaction)

              # set times
              times = unique(quantile(temp.data$.time.end, seq(1, 9, by = 2) / 10, na.rm = T))

              # get results
              results.i = results(object = model.i, predictions = prediction, draws = 500, times = times)

              # save
              r.inte = list(variable = variable, interaction = interaction, model = model.i, coefficients = results.i$coefficients, predictions = results.i$predictions, contrasts = results.i$contrasts)
            } else {
              r.inte = list(variable = variable, interaction = interaction, model = NULL, coefficients = NULL, predictions = NULL, contrasts = NULL)
            }
          }

          # return
          return(r.inte)
        })

        # remove nulls from results and set names
        r.interaction[sapply(r.interaction, is.null)] = NULL
        if(!is.null(r.interaction)) { n = purrr::map_chr(r.interaction, "interaction"); if(!is.null(n)) names(r.interaction) = n }
      } else {
        r.interaction = NULL
      }

      # get list of interaction variables to run
      mediation.variables = names(research.plan$formulas[[variable]]$mediation)

      # now run the mediation if we have any significant effects
      if(run.mediation && !is.null(mediation.variables) && is.list(r.main) && (is.null(r.main[[1]]$contrasts) || any(r.main[[1]]$contrasts$p.value < 0.1))) {
        # loop through mediation variables
        r.mediation = lapply(mediation.variables, function(mediation) {
          # make sure the interaction is not the same as the main variable
          if(mediation == variable) {
            return(NULL)
          }

          # check to see if we already have information or if we need to fill it in
          if(is.list(fill.in.plan) && rlang::has_name(fill.in.plan[[variable]], "mediation.contrasts") && !is.null(fill.in.plan[[variable]]$mediation.contrasts) &&
             dv %in% fill.in.plan[[variable]]$interaction.contrasts$.outcome && mediation %in% fill.in.plan[[variable]]$mediation.contrasts$.main.mediation) {
            # we already have information in a plan that we are updating so tell the user and skip the additional analysis
            cat(paste0("\n** Already have mediation results for -- ", variable, " -> ", mediation, ".\n"))

            # save recycle our existing results
            r.medi = list(variable = variable, mediation = mediation, model.med = NULL, model.out = NULL,
                          contrasts = dplyr::filter(fill.in.plan[[variable]]$mediation.contrasts, .outcome == dv, .main.mediation == mediation))
          } else {
            # need to deal with both formulas in turn -- first mediator then outcome

            # console output
            cat(paste0("\n** Running mediator -- ", variable, " -> ", mediation, ".\n"))

            # mediator formula
            new.formula.med = formula(research.plan$formulas[[variable]]$mediation[[mediation]]$mediator, to.drop = main.drop$dropped)

            # don't check thee mediator too -- just check the treatment for now
            new.formula.out = update(formula(research.plan$formulas[[variable]]$mediation[[mediation]]$outcome, to.drop = main.drop$dropped), research.plan$outcome)

            # merge random effects if present
            if(!is.null(research.plan$effects)) {
              new.formula.med = lasso2:::merge.formula(new.formula.med, research.plan$effects)
              new.formula.out = lasso2:::merge.formula(new.formula.out, research.plan$effects)
            }

            # first part of mediation
            model.med = NULL
            try(model.med <- analysis(
              runs = 500,
              formula = new.formula.med,
              main.ivs = variable,
              data = temp.data,
              inference = "bayesian",
              model.extra.args = list(prior = rstanarm::normal(0, 1), prior_intercept = rstanarm::normal(0, 1), adapt_delta = 0.9, warmup = 250, iter = 500)
            ), T)

            # second part of mediation
            model.out = NULL
            try(model.out <- analysis(
              runs = 500,
              formula = new.formula.out,
              main.ivs = variable,
              data = temp.data,
              inference = "bayesian",
              model.extra.args = list(prior = rstanarm::normal(0, 1), prior_intercept = rstanarm::normal(0, 1), adapt_delta = 0.9, warmup = 250, iter = 500)
            ), T)


            if(!is.null(model.med) & !is.null(model.out)) {
              # create prediction for treatment
              prediction = pr_list(!!variable := create_values(temp.data[[variable]]))

              # set times -- this needs to be fixed so that it is not hard coded
              times = unique(quantile(temp.data$.time.end, seq(1, 9, by = 2) / 10, na.rm = T))

              # get mediation results
              results.m = results_mediation(m.mediator = model.med, m.outcome = model.out, predictions = prediction, times = as.integer(median(times)), .outcome = dv)

              # save
              r.medi = list(variable = variable, mediation = mediation, model.med = model.med, model.out = model.out, contrasts = results.m)
            } else {
              r.medi = list(variable = variable, mediation = mediation, model.med = NULL, model.out = NULL, contrasts = NULL)
            }
          }

          # return
          return(r.medi)
        })

        # remove nulls from results and set names
        r.mediation[sapply(r.mediation, is.null)] = NULL
        if(!is.null(r.mediation)) { n = purrr::map_chr(r.mediation, "mediation"); if(!is.null(n)) names(r.mediation) = n }
      } else {
        r.mediation = NULL
      }

      # data to return
      r = list(variable = variable, outcome = dv, main = r.main, interaction = r.interaction, mediation = r.mediation)

      # return
      return(r)
    })

    # set the names
    if(rlang::has_name(r.outcome, "outcome")) {
      r.outcome = list(r.outcome)
    }
    names(r.outcome) = purrr::map_chr(r.outcome, "outcome")

    # add coefficients and predictions

    # get the main coefficients
    f.main.coefficients = dplyr::bind_rows(lapply(names(r.outcome), function(x) {
      if(!is.null(r.outcome[[x]]$main)) {
        df = purrr::map_dfr(r.outcome[[x]]$main, "coefficients", .id = ".main.variable")
        df = dplyr::mutate(df, .outcome = x, .before = ".main.variable")
      } else {
        return(NULL)
      }
    }))

    # get the main predictions
    f.main.predictions = dplyr::bind_rows(lapply(names(r.outcome), function(x) {
      if(!is.null(r.outcome[[x]]$main)) {
        df = purrr::map_dfr(r.outcome[[x]]$main, "predictions", .id = ".main.variable")
        df = dplyr::mutate(df, .outcome = x, .before = ".main.variable")
      } else {
        return(NULL)
      }
    }))

    # get the main effects
    f.main.contrasts = dplyr::bind_rows(lapply(names(r.outcome), function(x) {
      if(!is.null(r.outcome[[x]]$main)) {
        df = purrr::map_dfr(r.outcome[[x]]$main, "contrasts", .id = ".main.variable")
        df = dplyr::mutate(df, .outcome = x, .before = ".main.variable")
      } else {
        return(NULL)
      }
    }))

    # get the interaction coefficients -- the main variable is not named correctly -- fix!!!
    f.interaction.coefficients = dplyr::bind_rows(lapply(names(r.outcome), function(x) {
      if(!is.null(r.outcome[[x]]$interaction)) {
        df = dplyr::mutate(purrr::map_dfr(r.outcome[[x]]$interaction, "coefficients", .id = ".main.interaction"), .main.variable = x, .before = ".main.interaction")
        df = dplyr::mutate(df, .outcome = x, .before = ".main.variable")
      } else {
        return(NULL)
      }
    }))

    # get the interaction predictions
    f.interaction.predictions = dplyr::bind_rows(lapply(names(r.outcome), function(x) {
      if(!is.null(r.outcome[[x]]$interaction)) {
        df = dplyr::mutate(purrr::map_dfr(r.outcome[[x]]$interaction, "predictions", .id = ".main.interaction"), .main.variable = x, .before = ".main.interaction")
        df = dplyr::mutate(df, .outcome = x, .before = ".main.variable")
      } else {
        return(NULL)
      }
    }))

    # get the interaction contrasts
    f.interaction.contrasts = dplyr::bind_rows(lapply(names(r.outcome), function(x) {
      if(!is.null(r.outcome[[x]]$interaction)) {
        df = dplyr::mutate(purrr::map_dfr(r.outcome[[x]]$interaction, "contrasts", .id = ".main.interaction"), .main.variable = x, .before = ".main.interaction")
        df = dplyr::mutate(df, .outcome = x, .before = ".main.variable")
      } else {
        return(NULL)
      }
    }))

    # get the mediation
    f.mediation.contrasts = dplyr::bind_rows(lapply(names(r.outcome), function(x) {
      if(!is.null(r.outcome[[x]]$mediation)) {
        df = dplyr::mutate(purrr::map_dfr(r.outcome[[x]]$mediation, "contrasts", .id = ".main.mediation"), .main.variable = x, .before = ".main.mediation")
        df = dplyr::mutate(df, .outcome = x, .before = ".main.variable")
      } else {
        return(NULL)
      }
    }))

    # the full results for an outcome
    r = list(variable = variable,
             main.coefficients = f.main.coefficients, main.predictions = f.main.predictions, main.contrasts = f.main.contrasts,
             interaction.coefficients = f.interaction.coefficients, interaction.predictions = f.interaction.predictions, interaction.contrasts = f.interaction.contrasts,
             mediation.contrasts = f.mediation.contrasts)

    # return
    return(r)
  })

  # set names
  if(length(all.variables) > 0) {
    if(rlang::has_name(all.variables, "variable")) {
      all.variables = list(all.variables)
    }
    names(all.variables) = purrr::map_chr(all.variables, "variable")
  } else {
    all.variables = NULL
  }

  # final return
  return(all.variables)
}

