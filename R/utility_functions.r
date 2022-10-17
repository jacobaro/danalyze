# utility functions

#' Remove backticks.
#'
#' This function takes in a string and outputs a string with
#' backticks removed.
#'
#' @param str The string to remove backticks from.
#' @return A string with backticks removed.
#' @export
#'
remove_backticks = function(str) {
  r = sapply(str, function(x) { if(substring(x, 1, 1) == "`") x = substring(x, 2); l = nchar(x); if(substring(x, l, l) == "`") x = substring(x, 1, l - 1); x })
  r = as.character(r)
  r
}

#' Fast identification of baseline hazard.
#'
#' This function takes in a survival model and outputs a list with the
#' baseline hazard and time.
#'
#' @param model A survival model running using the "survival" package.
#' @return A list with the baseline hazard and time.
#' @export
#'
fast_bh = function(model) {
  # get survfit
  sfit = survival::survfit(model, se.fit = F)

  # remove means
  chaz = sfit$cumhaz * exp(-sum(model$means * model$coefficients))

  # return
  return(list(hazard = chaz, time = sfit$time))
}

# turn a formula into a character - surprisingly annoying!
as.character.formula = function(f) {
  paste(trimws(deparse(f)), collapse = " ")
}

#' Symmetric 'log1p'.
#'
#' This function takes in a positive or negative number and returns correctly
#' transformed values.
#'
#' @param x A numeric vector.
#' @return The transformed numeric vector.
#' @export
#'
symlog = function(x) { sign(x) * log1p(abs(x)) }

#' @rdname symlog
#' @export
#'
symsqrt = function(x) { sign(x) * sqrt(abs(x)) }

# TODO: bring in some functions from performance package to test model assumptions

# intersection of vectors
intersect_all = function(...) {
  Reduce(intersect, list(...))
}

#' Collapse time.
#'
#' Summarize a dataframe to different units of analysis. The 'time.var' is
#' divided by 'time.factor' and turned into an integer.
#'
#' @param data The data to summarize.
#' @param time.var Symbols identifying the time variable in data.
#' @param time.factor A number that 'time.var' will be divided by.
#' @param group.by A vector of symbols that identify the unit of analysis.
#' @param variables A list of variables to apply special functions to.
#' @return A summarized dataframe.
#' @export
#'
collapse_time = function(data, time.var, time.factor, group.by, variables) {
  # function to collapse time -- could really recode this to do start and end time using blocks of X time that would be neat

  # first make sure it is not grouped
  data = dplyr::ungroup(data)

  # set :=
  `:=` = rlang::`:=`

  # enquo
  time.var = rlang::enquo(time.var)
  group.by = rlang::enquo(group.by)

  # create the new time group
  data = dplyr::mutate(data, .time.factor = ceiling(!!time.var / time.factor))

  # create time per group
  data.group = dplyr::mutate(dplyr::group_by_at(data, dplyr::vars(!!group.by)), .time.group = !!time.var - min(!!time.var, na.rm = T) + 1)

  # grouped data
  data.group = dplyr::group_by_at(data.group, dplyr::vars(!!group.by, .time.factor))

  # summarize
  data.base = dplyr::summarize(data.group,
                               !!time.var := min(!!time.var, na.rm = T),
                               .time.start = min(.time.group, na.rm = T) - 1,
                               .time.end = max(.time.group, na.rm = T), .groups = "keep")

  # get column order
  column.order = c(colnames(data.base), colnames(data)[!colnames(data) %in% colnames(data.base)])

  # summarize data

  # identify user-supplied variable names
  user.var.names = unique(unlist(lapply(variables, function(x) colnames(dplyr::select(data, !!!x$vars)))))

  # make sure we are not trying to transform grouping or time variables
  user.var.names = user.var.names[!user.var.names %in% c(dplyr::group_vars(data.group), ".time.factor", ".time.group", ".time.start", ".time.end")]

  # columns to remove
  column.remove = colnames(dplyr::select(dplyr::ungroup(data.group), !!time.var, !!group.by, .time.factor, .time.group))

  # apply our special functions
  data.special = lapply(variables, function(x) {
    dplyr::select(dplyr::ungroup(dplyr::summarize_at(data.group, .vars = colnames(dplyr::select(data, !!!x$vars)), .funs = x$func)), -tidyselect::any_of(column.remove))
  })

  # bind
  data.special = dplyr::bind_cols(data.special)

  # remaining variables to select
  data.group.remain = dplyr::select(data.group, -!!user.var.names)

  # numeric data
  data.numeric = dplyr::select(dplyr::ungroup(dplyr::summarize_if(data.group.remain, is.numeric, mean, na.rm = T)), -tidyselect::any_of(column.remove))

  # categorical data
  data.categorical = dplyr::select(dplyr::ungroup(dplyr::summarize_if(data.group.remain, function(x) !is.numeric(x), function(x) dplyr::first(na.omit(x)))), -tidyselect::any_of(column.remove))

  # combine
  data.full = dplyr::bind_cols(lapply(list(data.base, data.special, data.numeric, data.categorical), dplyr::ungroup))

  # select columns
  data.full = dplyr::select(data.full, !!column.order)

  # return
  return(data.full)
}

#' Create high and low values
#'
#' This function takes in a vector (numeric or character) and outputs nicely formatted
#' test values.
#'
#' @param x A numeric vector.
#' @param .quantile The high and low quantiles for a numeric vector.
#' @param .places The number of places to round a numeric vector to.
#' @return For a numeric vector this function return high and low values.
#'   For a character vector this returns unique values.
#' @export
#'
create_values = function(x, .quantile = c(0.975, 0.025), .places = 2) {
  # make sure quantile is okay
  if(!is.numeric(.quantile) | !length(.quantile) == 2) {
    warning("Quantile is not structured correctly. Using default.")
    .quantile = c(0.975, 0.025)
  }

  # check
  if(is.numeric(x)) {
    # if there are a large number of zeros than we want to produce values absent zero
    if((length(x == 0) / length(x)) > 0.75) {
      x.t = x[x != 0]
    } else {
      x.t = x
    }

    # get quantile
    r = quantile(x.t, .quantile, na.rm = T)

    # make sure the values are okay
    if(diff(range(r)) < .Machine$double.eps ^ 0.5) {
      r = c(max(x, na.rm = T), min(x, na.rm = T))
    }

    # round
    r = round(r, .places)
  } else {
    r = unique(x)
  }

  # return
  return(r)
}

#' Nice model output.
#'
#' This function prints out a nice models summary and optionally incorporates clustered standard errors.
#'
#' @param m The model to print a summary for.
#' @param cluster An optional formula identifying variables to use for clustered standard errors.
#' @param drop.factor Drop factor variables from print output or keep them. Defaults to drop.
#' @param .level Level of statistical significance for confidence intervals. Defaults to P < 0.05.
#' @param .round Number of digits to round results to. Defaults to three.
#' @return For a numeric vector this function return high and low values.
#'   For a character vector this returns unique values.
#' @export
#'
ct_to_out = function(m, cluster = NULL, vcov = NULL, drop.factor = T, .level = 0.95, .round = 3) {
  # get coeftest
  if(any(c("lmerMod", "glmerMod", "lmerModLmerTest", "glmerModLmerTest") %in% class(m))) {
    mt = as.data.frame(summary(m)$coefficients)
    mt$p.value = (1 - pnorm(abs(mt[, 3]))) * 2
    colnames(mt) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
    mt = as.matrix(mt)
  } else {
    if(!is.null(vcov)) {
      mt = lmtest::coeftest(m, vcov)
    } else {
      mt = lmtest::coeftest(m, get_vcov(m, cluster = cluster, iv.fix = F))
    }
  }

  # drop factors if desired
  if(drop.factor) {
    mt = mt[!stringr::str_detect(rownames(mt), "factor\\("), , drop = F] # make sure it stays a matrix
    class(mt) = "coeftest"
  }

  # get data -- returns the full data frame with missing values for lmer models
  d.data = get_data(m)

  # t
  t = qnorm(0.5 + .level / 2)

  # create the output dataframe
  r = tibble::tibble(coefficient = rownames(mt), c = mt[, 1], c.low = mt[, 1] - t * mt[, 2], c.high = mt[, 1] + t * mt[, 2], p.value = mt[, 4], draws = nrow(d.data))

  # round
  if(!is.null(.round)) {
    r = dplyr::mutate_if(r, is.numeric, round, digits = .round)
  }

  # return
  return(list(coeftest = mt, output = r))
}

# output a table for a fixest model
create_table = function(mlist, title = "Coefficient Table", labels = NULL) {
  # # allow times font
  # windowsFonts(Times = windowsFont("Times New Roman"))

  # make list
  if(!is.list(mlist)) mlist = list("Model" = mlist)

  # t
  tval = qnorm(0.5 + 0.95 / 2)

  # function to format coefficients
  format_coefs = function(ms, name, group = "Variables", tval, keep = NULL) {
    # no coeftable so create one
    if(!tibble::has_name(ms, "coeftable")) {
      ms$coeftable = matrix(data = c(ms$coefficients, se = ms$ses, zs = ms$coefficients / ms$ses, p_value = 1 - pnorm(abs(ms$coefficients / ms$ses))), ncol = 4)
      rownames(ms$coeftable) = names(ms$coefficients)
    }

    # format coef table
    ct = tibble::tibble(
      var = rownames(ms$coeftable),
      value =
        paste0(scales::number(ms$coeftable[, 1], 0.001, big.mark = ","),
               dplyr::case_when(ms$coeftable[, 4] < 0.001 ~ "***", ms$coeftable[, 4] < 0.01 ~ "**", ms$coeftable[, 4] < 0.05 ~ "*", ms$coeftable[, 4] < 0.1 ~ "†", T ~ ""),
               "<br>(", scales::number(ms$coeftable[, 1] - tval * ms$coeftable[, 2], 0.001, big.mark = ","), " to ", scales::number(ms$coeftable[, 1] + tval * ms$coeftable[, 2], 0.001, big.mark = ","), ")"),
      group = group,
      model = name
    )

    # filter if needed
    if(!is.null(keep)) {
      ct = dplyr::filter(ct, var %in% keep)
    }

    # drop some vars
    to.drop = c("m-splines-coef")
    ct = dplyr::filter(ct, !str_detect(var, to.drop))

    # return
    return(ct)
  }

  # set name and order
  name_and_order = function(x, labels) {
    # if null return
    if(is.null(labels)) {
      return(x)
    }

    # length
    x.length = 1:nrow(x)

    # get ordering
    fe.order = na.omit(match(names(labels), x$var))

    # set name -- could also limit to variables & x$group == "Variables"
    new.names = labels[match(x$var, names(labels))]
    new.names[is.na(new.names)] = x$var[is.na(new.names)]
    x$var = new.names

    # add non-present back in
    fe.order = c(fe.order, x.length[!x.length %in% fe.order])

    # set ordering
    x = x[fe.order, ]

    # return
    return(x)
  }

  # loop through
  r = lapply(names(mlist), function(m) {
    # get summary
    ms = mlist[[m]]

    # model name
    model.name = paste0("(", which(m == names(mlist)), ")<br>", m)

    # format coef table
    ct = format_coefs(ms, name = model.name, tval = tval)

    # get first stage if iv analysis
    if(!is.null(ms$iv) && ms$iv == T) {
      # get ivs
      iv = lapply(ms$iv_first_stage, function(x) {
        dplyr::bind_rows(format_coefs(x, name = model.name, tval = tval, keep = ms$iv_inst_names_xpd),
                         tibble::tibble(var = "F-stat", value = as.character(scales::number(fixest::fitstat(x, "ivf")[[1]]$stat, 0.001, big.mark = ",")), group = "", model = model.name))
      })

      # bind
      iv = dplyr::bind_rows(iv, .id = "group")

      # set lavel if present
      if(!is.null(labels)) {
        iv$group = labels[match(iv$group, names(labels))]
      }

      # set group name
      iv$group = paste0("First Stage: ", iv$group)
    } else {
      iv = NULL
    }

    # random effects
    if(!is.null(ms$has_bars) && ms$has_bars) {
      re = tibble::tibble(
        var = sapply(lme4::findbars(ms$formula$formula), function(x) dplyr::last(all.vars(x))),
        value = "T",
        group = "Random Effects",
        model = model.name
      )
    } else {
      re = NULL
    }


    # fixed effects
    if(!is.null(ms$fixef_vars)) {
      fe = tibble::tibble(
        var = ms$fixef_vars,
        value = "T",
        group = "Fixed Effects",
        model = model.name
      )
    } else {
      fe = NULL
    }

    # identify standard error
    se = attr(ms$se, "type")
    if("formula" %in% class(se)) if(!is.null(labels)) se = paste(labels[match(all.vars(se), names(labels))], collapse = ", ") else se = paste(all.vars(se), collapse = ", ")

    # set r2
    if(c("fixest") %in% class(ms)) {
      if(is.null(ms$family) || ms$family == "gaussian") {
        r2 = c("Adj. R2" = as.numeric(fixest::fitstat(ms, "ar2")$ar2))
        fstat = c("F-stat" = fixest::fitstat(ms, "wf")[[1]]$stat)
      } else {
        r2 = c("Adj. Pseudo R2" = as.numeric(fixest::fitstat(ms, "apr2")$apr2))
        fstat = c("Log-likelihood" = fixest::fitstat(ms, "ll")[[1]])
      }

      # extra -- nobs vs. nobs_origin
      extra = tibble::tibble(
        var = c("Observations", names(r2), "SE", names(fstat)),
        value = c(scales::number(ms$nobs_origin, big.mark = ","), round(r2, 3), se, as.character(scales::number(fstat, 0.001, big.mark = ","))),
        group = "Model Information",
        model = model.name
      )
    } else {
      extra = tibble::tibble(
        var = c("Observations"),
        value = c(scales::number(nrow(ms$data), big.mark = ",")),
        group = "Model Information",
        model = model.name
      )
    }



    # combine
    return(dplyr::bind_rows(list(ct, iv, re, fe, extra)))
  })

  # bind
  r = dplyr::bind_rows(r)

  # arrange the group variables, iv variables, fixed effects, model information
  r$group = factor(r$group, levels = unique(r$group))

  # set labels
  r = name_and_order(r, labels)

  # arrange back to correct group order
  r = dplyr::arrange(r, group)

  # widen
  r = tidyr::pivot_wider(r, names_from = "model", values_from = "value", values_fill = "")

  # column names
  tcoln = lapply(colnames(r)[!colnames(r) %in% c("var", "group")], gt::md)
  names(tcoln) = colnames(r)[!colnames(r) %in% c("var", "group")]

  # create table -- should get rid of magittr reference
  rt = r %>%
    gt::gt(rowname_col = "var", groupname_col = "group") %>%
    gt::tab_header(
      title = gt::md(title)
      # subtitle = gt::md("")
    ) %>%
    gt::tab_source_note(gt::md("P-values: *** < 0.001; ** < 0.01; * < 0.05; † < 0.1.<br>Confidence intervals (95%) shown in parentheses below coefficients.")) %>%
    gt::fmt_markdown(columns = c(-var, -group)) %>% gt::cols_label(!!!tcoln) %>%
    gt::cols_align("center", columns = c(-var, -group)) %>%
    gt::tab_options(table.font.names = c("Times New Roman", "Times", "times"), table.font.size = gt::px(11))

  # return
  return(rt)
}
