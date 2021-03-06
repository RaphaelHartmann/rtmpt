% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SBC_fit_rtmpt.R
\name{fit_rtmpt_SBC}
\alias{fit_rtmpt_SBC}
\title{Simulate data from RT-MPT models}
\usage{
fit_rtmpt_SBC(model, seed, n.eff_samples = 99, n.chains = 4,
  n.iter = 5000, n.burnin = 200, n.thin = 1, Rhat_max = 1.05,
  Irep = 1000, n.subj = 40, n.trials = 30, prior_params = NULL)
}
\arguments{
\item{model}{A list of the class \code{rtmpt_model}.}

\item{seed}{Random seed number.}

\item{n.eff_samples}{Number of effective samples. Default is 99, leading to 100 possible ranks (from 0 to 99).}

\item{n.chains}{Number of chains to use. Default is 4. Must be larger than 1 and smaller or equal to 16.}

\item{n.iter}{Number of samples per chain. Default is 5000. Must be larger or equal to \code{n.eff_samples}.}

\item{n.burnin}{Number of warm-up samples. Default is 200.}

\item{n.thin}{Thinning factor. Default is 1.}

\item{Rhat_max}{Maximal Potential scale reduction factor: A lower threshold that needs to be reached before the actual sampling starts. Default is 1.05}

\item{Irep}{Every \code{Irep} samples an interim state with the current maximal potential scale reduction
factor is shown. Default is 1000. The following statements must hold true for \code{Irep}:
\itemize{
  \item \code{n.burnin} is smaller than or equal to \code{Irep},
  \item \code{Irep} is a multiple of \code{n.thin} and
  \item \code{n.iter} is a multiple of \code{Irep / n.thin}.
}}

\item{n.subj}{Number of subjects. Default is 40.}

\item{n.trials}{Number of trials per tree. Default is 30.}

\item{prior_params}{Named list of parameters from which the data will be generated. This must be the same named list as \code{prior_params} from 
\code{\link{fit_rtmpt}} and has the same defaults. It is not recommended to use the defaults since they lead to many probabilities close or 
equal to \code{0} and/or \code{1} and to RTs close or equal to \code{0}. Allowed parameters are:
\itemize{
  \item \code{mean_of_exp_mu_beta}: This is the expected exponential rate (\code{E(exp(beta)) = E(lambda)}) and 
        \code{1/mean_of_exp_mu_beta} is the expected process time (\code{1/E(exp(beta)) = E(tau)}). The default
        mean is set to \code{10}, such that the expected process time is \code{0.1} seconds.
  \item \code{var_of_exp_mu_beta}: The group-specific variance of the exponential rates. Since
        \code{exp(mu_beta)} is Gamma distributed, the rate of the distribution is just mean divided by variance and
        the shape is the mean times the rate. The default is set to \code{100}.
  \item \code{mean_of_mu_gamma}: This is the expected \emph{mean parameter} of the encoding and response execution times,
        which follow a normal distribution truncated from below at zero, so \code{E(mu_gamma) < E(gamma)}. The default is \code{0}.
  \item \code{var_of_mu_gamma}: The group-specific variance of the \emph{mean parameter}. Its default is \code{10}.
  \item \code{mean_of_omega_sqr}: This is the expected residual variance (\code{E(omega^2)}). The default is \code{0.005}.
  \item \code{var_of_omega_sqr}: The variance of the residual variance (\code{Var(omega^2)}). The default is
        \code{0.01}. The default of the mean and variance is equivalent to a shape and rate of \code{0.0025} and 
        \code{0.5}, respectivly.
  \item \code{df_of_sigma_sqr}: degrees of freedom for the individual variance of the response executions. The
        individual variance follows a scaled inverse chi-squared distribution with \code{df_of_sigma_sqr} degrees of freedom and
        \code{omega^2} as scale. \code{2} is the default and it should be an integer.
  \item \code{sf_of_scale_matrix_SIGMA}: The original scaling matrix (S) of the (scaled) inverse Wishart distribution for the process 
        related parameters is an identity matrix \code{S=I}. \code{sf_of_scale_matrix_SIGMA} is a scaling factor, that scales this 
        matrix (\code{S=sf_of_scale_matrix_SIGMA*I}). Its default is \code{1}.
  \item \code{sf_of_scale_matrix_GAMMA}: The original scaling matrix (S) of the (scaled) inverse Wishart distribution for the encoding and
        motor execution parameters is an identity matrix \code{S=I}. \code{sf_of_scale_matrix_GAMMA} is a scaling factor that scales 
        this matrix (\code{S=sf_of_scale_matrix_GAMMA*I}). Its default is \code{1}.
  \item \code{prec_epsilon}: This is epsilon in the paper. It is the precision of mu_alpha and all xi (scaling parameter
        in the scaled inverse Wishart distribution). Its default is also \code{1}.
  \item \code{add_df_to_invWish}: If \code{P} is the number of parameters or rather the size of the scale matrix used in the (scaled)
        inverse Wishart distribution then \code{add_df_to_invWish} is the number of degrees of freedom that can be added to it. So
        \code{DF = P + add_df_to_invWish}. The default for \code{add_df_to_invWish} is \code{1}, such that the correlations are uniformly 
        distributed within \code{[-1, 1]}.
}}
}
\value{
A list of the class \code{rtmpt_sbc} containing 
  \itemize{
    \item \code{ranks}: the rank statistic for all parameters,
    \item \code{sim_list}: an object of the class \code{rtmpt_sim},
    \item \code{fit_list}: an object of the class \code{rtmpt_fit},
    \item \code{specs}: some specifications like the model, seed number, etc.,
  }
}
\description{
Simulate data from RT-MPT models using \code{rtmpt_model} objects. The difference to \code{\link{sim_rtmpt_data}} is that here only scalars are allowed. This makes it usable for 
  simulation-based calibration (SBC; Talts et al., 2018). You can specify the random seed, number of subjects, number of trials, and some
  parameters (same as \code{prior_params} from \code{\link{fit_rtmpt}}).
}
\examples{
########################################################################################
# Detect-Guess variant of the Two-High Threshold model.
# The encoding and motor execution times are assumed to be different for each response.
########################################################################################

mdl_2HTM <- "
# targets
d+(1-d)*g     ; 0
(1-d)*(1-g)   ; 1

# lures
(1-d)*g       ; 0
d+(1-d)*(1-g) ; 1

# d: detect; g: guess
"

model <- to_rtmpt_model(mdl_file = mdl_2HTM)

params <- list(mean_of_exp_mu_beta = 10, 
               var_of_exp_mu_beta = 10, 
               mean_of_mu_gamma = 0.5, 
               var_of_mu_gamma = 0.0025, 
               mean_of_omega_sqr = 0.005, 
               var_of_omega_sqr = 0.000025,
               df_of_sigma_sqr = 10, 
               sf_of_scale_matrix_SIGMA = 0.1, 
               sf_of_scale_matrix_GAMMA = 0.01, 
               prec_epsilon = 10,
               add_df_to_invWish = 5)
\donttest{
SBC_out <- fit_rtmpt_SBC(model, seed = 123, prior_params = params)
SBC_out$ranks
}
\donttest{
# For 2000 replications
## This takes too long to run and in addition Rhat should always be 
## checked as well as the effective sample size.
R = 2000
rank_mat <- data.frame()

for (r in 1:R) {
  SBC_out <- fit_rtmpt_SBC(model, seed = r*123, prior_params = params, n.eff_samples = 99)
  rank_mat <- rbind(rank_mat, SBC_out$ranks)
}

## pearson chi square for testing uniformity
x <- apply(rank_mat[1:R,], 2, table)
expect <- R/100 # 100 = number of bins/cells (0:99)
pearson <- apply(X = x, MARGIN = 2, FUN = function(x) {sum((x-expect)^2/expect)})
z95 <- qchisq(0.95, 99) # 99 = degrees of freedom
sum(pearson>z95) / length(pearson)
}
}
\references{
Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A. (2018). Validating Bayesian inference algorithms with simulation-based calibration. \emph{arXiv preprint arXiv:1804.06788}.
}
\author{
Raphael Hartmann
}
