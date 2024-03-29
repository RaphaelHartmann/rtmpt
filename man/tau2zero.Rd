% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_params.R
\name{tau2zero}
\alias{tau2zero}
\alias{set_tau_zero}
\title{Set process completion times to zero}
\usage{
tau2zero(model, names, outcomes, values = 0)

set_tau_zero(model, names, outcomes, values = 0)
}
\arguments{
\item{model}{A list of the class \code{ertmpt_model}.}

\item{names}{Character vector with process names.}

\item{outcomes}{Character vector of length \code{length(names)} indicating for which process outcome the process completion time should 
be zero or changed back to be estimated. Allowed characters are:
\itemize{
 \item \code{"minus"}: the negative outcome of the process.
 \item \code{"plus"}: the positive outcome of the process.
}}

\item{values}{Numerical vector of length one or \code{length(names)}. You have the following options for the elements of the numeric vector:
\itemize{
 \item \code{0}: suppress the process time/rate, i.e., set the process completion time (tau) with the specified output to zero.
 \item \code{NA}: estimate the process time (tau)
}}
}
\value{
A list of the class \code{ertmpt_model}.
}
\description{
Setting process completion times (taus) to zero or change it back to be estimated.
}
\examples{
####################################################################################
# Detect-Guess variant of the Two-High Threshold model.
# The encoding and motor execution times are assumed to be equal for each category.
# The process completion times for both failed detections will be suppressed.
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

model <- to_ertmpt_model(mdl_file = mdl_2HTM)

## removing the process times (tau) for the failed (="minus") detection ("do" and "dn") 
new_model <- tau2zero(model = model, names = c("dn", "do"),
                      outcomes = c("minus", "minus"), values = 0)
new_model


## removing the process times (tau) for the failed (="minus") detection ("do" and "dn") 
new_model <- set_tau_zero(model = model, names = c("dn", "do"),
                          outcomes = c("minus", "minus"), values = 0)
new_model
}
\seealso{
\code{\link{delta2delta}}, \code{\link{theta2const}}, \code{\link{theta2theta}} and \code{\link{tau2tau}}
}
\author{
Raphael Hartmann
}
