# creating a "research plan" by taking some mapping of treatments/controls and using that to create a set of test formulas and new variables

# add transform to a var in a formula
add_transform = function(x, transforms = NULL, lag = NULL) {
  # strip lags
  x.nolag = stringr::str_replace_all(x, paste(paste0("_l", lag), collapse = "|"), "")

  # identify transforms
  x.transform = transforms$transform[match(x.nolag, transforms$var.name)]

  # add the transforms
  r = sapply(1:length(x), function(t) if(!is.na(x.transform[t])) paste0(x.transform[t], "(", x[t], ")") else x[t])

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

# create a research plan for the given variables

#' Function to get create new variables and formulas for testing a set of variable relationships.
#'
#' This function produces a set of formulas and new data for testing a research plan.
#' @param treatment Vector containing variable names for treatments to test.
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
                         prefix = list(source = "sv", target = "tv", other = "oth"), unit = list(unit = NULL, source = NULL, target = NULL, time = NULL), lag = NULL) {
  # below, we identify the variables we need to use, create new data to correspond to these variables, create formulas to use the variables, and then return everything

  ## CREATE A TABLE WITH THE VARAIBLES

  # check to make sure arguments are good
  if(!is.character(treatment) || !is.data.frame(data)) {
    stop("The treatment must be a character string and the data must be a dataframe.")
  }

  # make a table of all variables
  all.vars = list("treatment" = treatment, "control" = control, "interaction" = interaction, "mediation" = mediation)
  all.vars = dplyr::bind_rows(lapply(names(all.vars), function(x) if(!is.null(all.vars[[x]])) tibble::tibble(type = x, var.label = names(all.vars[[x]]), var.name = all.vars[[x]])))

  # set source, target, and relative
  all.vars.st = dplyr::bind_rows(lapply(all.vars$var.name, function(x) {
    cn = colnames(data)
    tibble::tibble(
      source = if(paste(prefix$source, x, sep = ".") %in% cn) paste(prefix$source, x, sep = ".") else x, # if there is no source + target it just gets shoved into source
      target = if(paste(prefix$target, x, sep = ".") %in% cn) paste(prefix$target, x, sep = ".") else NA_character_,
      relative = if(!is.na(source) & !is.na(target)) paste("rel", x, sep = ".") else NA_character_
    )
  }))

  # combine
  all.vars = dplyr::bind_cols(all.vars, all.vars.st)

  # set the class
  all.vars = dplyr::mutate(all.vars, factor = sapply(source, function(x) if(!is.na(x) & !is.numeric(data[[x]])) T else F))

  # set default transform
  all.vars$var.transform = sapply(all.vars$factor, function(x) if(x) function(y) as.numeric(as.factor(y)) else NA)

  # add functions
  if(length(factor.transform) > 0) {
    # only one transform function per variable
    factor.transform = factor.transform[unique(names(factor.transform))]

    # link it
    all.vars$var.transform[which(all.vars$var.name %in% names(factor.transform))] = factor.transform[na.omit(match(all.vars$var.name, names(factor.transform)))]
  }

  # set default var.transform


  # check for missing variables
  all.vars$missing = is.na(all.vars$source) & is.na(all.vars$target)

  # send a warning if variables are missing
  if(any(all.vars$missing)) {
    warning(paste("The following variables are missing both source and target names and will be excluded from the plan:", paste(all.vars$var.name[all.vars$missing], collapse = ", ")))
  }

  # select the columns we need
  all.vars = dplyr::select(dplyr::filter(all.vars, missing == F), type, factor, var.transform, var.label, var.name, source, target, relative)



  ## SELECT THE DATA

  # vars to select
  vars.to.select = unique(na.omit(c(keep.vars, unlist(unit), all.vars$source, all.vars$target)))

  # make sure we have the variables -- shoud modify to allow it to continue without these variables
  if(!all(vars.to.select %in% colnames(data))) {
    warning(paste("The following variables are not in the data:", paste(vars.to.select[!vars.to.select %in% colnames(data)], collapse = ", ")))
  }

  # select
  .data = dplyr::select(data, tidyselect::any_of(vars.to.select))


  ## CREATE NEW VARIABLES

  # the extent of what we want to create: lags of treatment, control, interaction and mediation; other values of treatment against source/target

  # first create relative variables
  for(i in 1:nrow(all.vars)) {
    # we have a source and target so we can make a relative variable -- also probably just need to do unique values
    if(!is.na(all.vars$relative[i])) {
      # set source and target
      t.source = .data[[all.vars$source[i]]]
      t.target = .data[[all.vars$target[i]]]

      # if not numeric then process
      if(all.vars$factor[i]) {
        t.source = all.vars$var.transform[[i]](t.source)
        t.target = all.vars$var.transform[[i]](t.target)
      }

      # set relative var
      .data[[paste("rel", all.vars$var.name[i], sep = ".")]] = dplyr::if_else(t.source == 0 & t.target == 0, 0, t.source / (t.source + t.target))
    }

    # could create relative variables in this way or in the formula: ~ source + target + I(source / (source + target))
  }

  # create raw vector for other values -- these are just created for the treatment
  all.vars$source.other = NA_character_
  all.vars$target.other = NA_character_

  # now create other values and add to our data frame -- need to deal with FACTOR VARS
  other.source.vars = build_other_values(.data, !!all.vars$var.name[all.vars$type == "treatment"], unit = c(unit$unit, unit$source, unit$time), prefix = paste(prefix$other, prefix$source, sep = "."))
  .data = dplyr::bind_cols(.data, other.source.vars)
  all.vars$source.other[all.vars$type == "treatment"] = colnames(other.source.vars)

  # do also for target if present
  if(!is.null(unit$target)) {
    other.target.vars = build_other_values(.data, !!all.vars$var.name[all.vars$type == "treatment"], unit = c(unit$unit, unit$target, unit$time), prefix = paste(prefix$other, prefix$target, sep = "."))
    .data = dplyr::bind_cols(.data, other.target.vars)
    all.vars$target.other[all.vars$type == "treatment"] = colnames(other.target.vars)
  }

  # arrange
  .data = dplyr::arrange_at(.tbl = .data, .vars = c(unit$unit, unit$source, unit$target, unit$time))

  # add time lags
  if(is.numeric(lag)) {
    # now create lags
    .data = dplyr::mutate_at(
      dplyr::group_by_at(.tbl = .data, .vars = c(unit$unit, unit$source, unit$target)),
      .vars = unique(na.omit(c(all.vars$source, all.vars$target, all.vars$relative, all.vars$source.other, all.vars$target.other))),
      .funs = sapply(lag, function(x) rlang::list2(!!paste0("l", x) := formula(paste0("~ lag(., ", x, ")"))))
    )

    # TODO: need to deal with just one .vars
  }

  # ungroup data
  .data = dplyr::ungroup(.data)


  ## IDENTIFY PROBLEMATIC VARIABLES

  # select -- ordering so we keep source and target as highest priority
  .data.col = dplyr::select(.data, unique(na.omit(c(unit$unit, all.vars$relative, all.vars$source.other, all.vars$target.other, all.vars$source, all.vars$target))))

  # remove entirely problematic variables
  bad.var = colnames(.data.col)[unlist(dplyr::summarise_all(.data.col, function(x) any(is.infinite(x)) | any(is.nan(x))))]

  # first find variables with little variance -- exclude unit and treatment from check
  no.var = caret::nearZeroVar(.data.col[!colnames(.data.col) %in% c(unit$unit, all.vars$source[all.vars$type == "treatment"])], freqCut = 95/5, uniqueCut = 0, foreach = T, names = T)

  # find variables that dont vary within our groups -- exclude treatment (unit is part of group_by)
  .data.group = dplyr::group_by_at(.tbl = .data.col[!colnames(.data.col) %in% all.vars$source[all.vars$type == "treatment"]], .vars = unit$unit)
  group.variation = dplyr::group_map(.data.group, ~caret::nearZeroVar(.x, freqCut = 99/1, uniqueCut = 0, foreach = T, names = T)) # less stringent criteria
  group.variation = 1 - unlist(as.list(table(unlist(group.variation)))) / dplyr::n_groups(.data.group) # what portion of groups have variation
  no.group.var = names(group.variation)[group.variation < 0.01]

  # now identify collinear vars -- drop factors could eventually look for linear combinations for factors -- also drops no variation vars since we already exclude those
  data.for.cor = dplyr::select_if(.data.col[!colnames(.data.col) %in% c(unit$unit, all.vars$source[all.vars$type == "treatment"], no.var, bad.var)], is.numeric)
  col.var = caret::findCorrelation(cor(data.for.cor, use = "pairwise.complete"), cutoff = 0.95, names = T, exact = T)

  # set final no var -- these are good candidates to drop from our analysis since there is so little variation -- tell the user though
  all.prob.vars = unique(c(no.var, no.group.var, col.var))

  # we could also try to find linear combinations but probably not needed at this step -- can also just rely on priors to address this

  # remove vars that are problematic and drop vars that no longer have any values to test
  all.vars = dplyr::mutate_at(all.vars, .vars = dplyr::vars(source:target.other), function(x) { x[x %in% all.prob.vars] = NA_character_; x })
  all.vars = dplyr::filter(all.vars, !is.na(source) | !is.na(target) | !is.na(relative) | !is.na(source.other) | !is.na(target.other))


  ## IDENTIFY VARIABLE TRANSFORMATIONS

  # helpful: https://www.ibm.com/support/pages/transforming-variable-normality-parametric-statistics

  # need to take care of source, target, other source, other target, and source/target/other lags
  all.transform = tibble::tibble(var.name = unique(na.omit(c(all.vars$source, all.vars$target, all.vars$source.other, all.vars$target.other))))

  # dont transform factor variables so just drop them for now
  all.transform = all.transform[sapply(all.transform$var.name, function(x) is.numeric(.data[[x]])), ]

  # set sign
  all.transform$sign = sapply(all.transform$var.name, function(x) {
    dat = na.omit(.data[[x]])
    dplyr::case_when(all(dat >= 0) & all(dat <= 0) ~ "portion", all(dat <= 0) ~ "negative", all(dat >= 0) ~ "positive", T ~ "both")
  })

  # identify portion min or max
  all.transform$portion.min = sapply(all.transform$var.name, function(x) mean(.data[[x]] == min(.data[[x]], na.rm = T), na.rm = T))
  all.transform$portion.max = sapply(all.transform$var.name, function(x) mean(.data[[x]] == max(.data[[x]], na.rm = T), na.rm = T))

  # identify skewness: < 1 = normal to half-normal; < 2 = exponential; > 2 = lognormal
  all.transform$skewness = sapply(all.transform$var.name, function(x) e1071::skewness(.data[[x]], na.rm = T, type = 2))

  # identify kurtosis: < 0 = fat tails; > 0 = thin tails; ~0 = normal tails
  all.transform$kurtosis = sapply(all.transform$var.name, function(x) e1071::kurtosis(.data[[x]], na.rm = T, type = 2))

  # set the transform -- when testing mediation analysis its possible for the mediation equation to produce negative values of an otherwise positive variable -- NEED TO ADDRESS!!
  all.transform$transform =
    dplyr::case_when(
      # all.transform$skewness > 1.5 & all.transform$portion.min < 0.5 & all.transform$sign == "positive" ~ "log1p",
      # all.transform$skewness > 1.5 & all.transform$sign == "positive" ~ "sqrt",
      all.transform$skewness > 1.5 & all.transform$portion.min < 0.5 ~ "danalyze::symlog",
      all.transform$skewness > 1.5 ~ "danalyze::symsqrt",
      all.transform$sign == "portion" ~ "asin",
      T ~ "scale"
    )


  ## CREATE FORMULAS FOR TESTING

  # idea is to create a formula: each treatment, each treatment interacted with a condition, each treatment as a cause of a mediator (+ the treatment and the mediator)

  # ultimtely we will want to identify variables that don't vary enough over time, within a unit, or between units

  # loop through each treatment
  all.formulas = lapply(all.vars$source[all.vars$type == "treatment"], function(t) {
    # assemble list of variables to include
    t.vars = as.character(na.omit(c(
      # other treatments -- keep for current treatment too: [!c(colnames(other.source.vars), colnames(other.target.vars)) %in% c(paste(paste(prefix$other, prefix$source, sep = "."), t, sep = "."), paste(paste(prefix$other, prefix$target, sep = "."), t, sep = "."))]
      all.vars$source.other,
      all.vars$target.other,

      # source controls
      all.vars$source[all.vars$type == "control"],

      # target controls
      all.vars$target[all.vars$type == "control"],

      # relative controls
      all.vars$relative[all.vars$type == "control"]
    )))

    # add lags
    if(is.numeric(lag)) {
      t.lags = as.character(sapply(paste0("l", lag), function(x) paste(c(t, t.vars), x, sep = "_")))
    } else {
      t.lags = NULL
    }

    # pull from rlang
    `:=` = rlang:::`:=`

    # for each formula block we need: (1) list of treatment + interaction + mediator -- tried this using an r object

    # old: rlang::list2(!!t := as.formula(paste("~", paste(add_transform(c(t, t.vars, t.lags)), collapse = " + ")), env = .GlobalEnv))

    # combine to create the main formula
    f.base = create_new_formula(iv = t, controls = c(t.vars[!t.vars %in% c(t)], t.lags), transforms = all.transform, lag = lag)

    # create formulas for interactions
    f.int.vars = na.omit(c(all.vars$source[all.vars$type == "interaction"], all.vars$target[all.vars$type == "interaction"], all.vars$relative[all.vars$type == "interaction"]))
    f.int = lapply(f.int.vars, function(i) create_new_formula(iv = c(t, i), controls = c(t.vars[!t.vars %in% c(i, t)], t.lags), transforms = all.transform, lag = lag))

    # create formulas for mediating effects -- requires two "med ~ iv" and "outcome ~ iv + med"
    f.med.vars = dplyr::if_else(sapply(all.vars$target[all.vars$type == "mediation"], is.null), all.vars$source[all.vars$type == "mediation"], all.vars$target[all.vars$type == "mediation"])
    f.med = lapply(f.med.vars, function(m) list(mediator = create_new_formula(dv = m, iv = t, controls = c(t.vars[!t.vars %in% c(m, t)], t.lags), transforms = all.transform, lag = lag),
                                                outcome = create_new_formula(iv = c(t, m), type = "+", controls = c(t.vars[!t.vars %in% c(m, t)], t.lags), transforms = all.transform, lag = lag)))

    # create full return
    f.full = list(main = f.base, interaction = f.int, mediation = f.med)

    # return
    return(f.full) #rlang::list2(!!t := f.full))
  })
  names(all.formulas) = all.vars$source[all.vars$type == "treatment"]

  # the returned data frame is missing the outcome and any fixed/random effects (i.e., all variables that are not IVs in some way) -- need to think about how to deal with this


  ## RETURN

  return(list(formulas = all.formulas, data = .data, variables = list(name = all.vars, transform = all.transform)))
}