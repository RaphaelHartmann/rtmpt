% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_params.R
\name{nu2const}
\alias{nu2const}
\alias{set_nu_const}
\title{Set process drift rate to constants}
\usage{
nu2const(model, names, constants = NA)

set_nu_const(model, names, constants = NA)
}
\arguments{
\item{model}{An object of the class \code{rtmpt_model}.}

\item{names}{Character vector with process names.}

\item{constants}{Numerical vector of length one or \code{length(names)}. You have the following options for the elements of the numeric vector:
\itemize{
 \item \code{-Inf < constants < Inf}: set the named drift rate parameter(s) to constant value(s)
 \item \code{NA}: estimate the named drift rate parameter(s)
}}
}
\value{
An object of the class \code{drtmpt_model}.
}
\description{
Setting process drif rate (parameter nu) to constants or change it back to be estimated.
}
\examples{
####################################################################################
# Detect-Guess variant of the Two-High Threshold model.
# The encoding and motor execution times are assumed to be equal for each category.
# The process drift rate for guessing (g) will be set to 1.0.
####################################################################################

mdl_2HTM <- "
# targets
do+(1-do)*g
(1-do)*(1-g)

# lures
(1-dn)*g
dn+(1-dn)*(1-g)

# do: detect old; dn: detect new; g: guess
"

model <- to_drtmpt_model(mdl_file = mdl_2HTM)

## setting drift rate for g to a constant (1.0):
new_model <- nu2const(model = model, names = c("g"), constants = c(1.0))
new_model


## setting drift rate of g to a constant (1.0):
new_model <- set_nu_const(model = model, names = c("g"), constants = c(1.0))
new_model
}
\seealso{
\code{\link{delta2delta}}, \code{\link{a2const}}, \code{\link{a2a}}, \code{\link{nu2nu}}, \code{\link{omega2const}} and \code{\link{omega2omega}}
}
\author{
Raphael Hartmann
}
