% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze.r
\name{model_run.bayesian}
\alias{model_run.bayesian}
\title{Run a Bayesian model}
\usage{
model_run.bayesian(chain.id, formula.parsed, model.functions)
}
\arguments{
\item{chain.id}{A unique ID to assign to the chain.}

\item{formula.parsed}{The parsed formula object.}

\item{model.functions}{The model object to run on the data.}
}
\value{
A return list that includes both the full model and the 'stanfit' component.
}
\description{
This internal function runs a Bayesian model and returns the results. Each run uses a single chain/core.
}
