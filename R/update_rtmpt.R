
#' Update Diffusion RT-MPT Models
#' 
#' Given a diffusion RT-MPT fit object and the data, this function updates the fit. All information is taken from the fit object, except
#'   the number of iterations and thinning as well as some flags and control options.
#' @param fit A list of the class \code{drtmpt_model}.
#' @param data Optimally, a list of class \code{drtmpt_data}. Also possible is a \code{data.frame} or a
#'   path to the text file. Both, \code{data.frame} and the text file must contain the column names "subj",
#'   "group", "tree", "cat", and "rt" preferably but not necessarily in this order. The values of the latter must
#'   be in milliseconds. It is always advised to use \code{\link{to_drtmpt_data}} first, which gives back a
#'   \code{drtmpt_data} list with information about the changes in the data, that were needed.
#' @param n.iter Number of samples per chain. Default is 1000.
#' @param n.thin Thinning factor. Default is 1.
#' @param Irep Every \code{Irep} samples an interim state with the current maximal potential scale reduction
#'   factor is shown. Default is 1000. The following statements must hold true for \code{Irep}:
#'   \itemize{
#'     \item \code{Irep} is a multiple of \code{n.thin} and
#'     \item \code{n.iter} is a multiple of \code{Irep / n.thin}.
#'   }
#' @param flags Either NULL or a list of
#'   \itemize{
#'     \item \code{old_label} If set to \code{TRUE} the old labels of "subj" and "group" of the data will be used in the elements of the output list.
#'       Default is \code{FALSE}.
#'     \item \code{indices} Model selection indices. If set to \code{TRUE} the log-likelihood for each iteration and trial will be
#'       stored temporarily and with that the DIC, as well as the WAIC and LOOIC will be calculated via the \code{loo} package. If you want to have the
#'       log-likelihood matrix stored in the output of this function, you can set \code{loglik} to \code{TRUE}. Default for
#'       \code{indices} is \code{FALSE}.
#'     \item \code{loglik} If set to \code{TRUE} and \code{indices = TRUE} the log-likelihood matrix for each iteration and trial will
#'       be saved in the output as a matrix. Default is \code{FALSE}.
#'   }
#' @param control Either NULL or a list of
#'   \itemize{
#'     \item \code{maxthreads} for the ML estimation of the initial values and the calculation of the DIC values one can use more than
#'       \code{n.chains} threads for parallelization. Default is 4 like \code{n.chians}. \code{maxthreads} must be larger or equal to
#'       \code{n.chains}.
#'   }
#' @return A list of the class \code{drtmpt_fit} containing
#'   \itemize{
#'     \item \code{samples}: the posterior samples as an \code{mcmc.list} object,
#'     \item \code{diags}: some diagnostics like deviance information criterion, posterior predictive checks for the frequencies and latencies,
#'                         potential scale reduction factors, and also the 99\% and 95\% HDIs and medians for the group-level parameters,
#'     \item \code{specs}: some model specifications like the model, arguments of the model call, and information about the data transformation,
#'     \item \code{indices} (optional): if enabled, WAIC and LOO,
#'     \item \code{LogLik} (optional): if enabled, the log-likelihood matrix used for WAIC and LOO.
#'     \item \code{summary} includes posterior mean and median of the main parameters.
#'   }
#' @references
#' Klauer, K. C. (2010). Hierarchical multinomial processing tree models: A latent-trait approach. \emph{Psychometrika, 75(1)}, 70-98.
#'
#' Spiegelhalter, D. J., Best, N. G., Carlin, B. P., & Van Der Linde, A. (2002). Bayesian measures of model complexity and fit.
#'   \emph{Journal of the royal statistical society: Series b (statistical methodology), 64(4)}, 583-639.
#'
#' Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC.
#'   \emph{Statistics and Computing, 27(5)}, 1413-1432.
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory.
#'   \emph{Journal of Machine Learning Research, 11(Dec)}, 3571-3594.
#'
#' @examples
#' ####################################################################################
#' # Detect-Guess variant of the Two-High Threshold model.
#' # The encoding and motor execution times are assumed to be equal for each response.
#' ####################################################################################
#'
#' mdl_2HTM <- "
#' # targets
#' do+(1-do)*g
#' (1-do)*(1-g)
#'
#' # lures
#' (1-dn)*g
#' dn+(1-dn)*(1-g)
#'
#' # do: detect old; dn: detect new; g: guess
#' "
#'
#' model <- to_drtmpt_model(mdl_file = mdl_2HTM)
#'
#' data_file <- system.file("extdata/data.txt", package="rtmpt")
#' data <- read.table(file = data_file, header = TRUE)
#' data_list <- to_drtmpt_data(raw_data = data, model = model)
#' \donttest{
#' # This might take some time
#' drtmpt_out <- fit_drtmpt(model = model, data = data_list, Rhat_max = 1.1)
#' 
#' drtmpt_updt <- update_drtmpt(fit = drtmpt_out, data = data_list, 
#'                              n.iter = 2000, n.thin = 1, Irep = 2000)
#' }
#' @author Raphael Hartmann
#' @useDynLib "rtmpt", .registration=TRUE
#' @export
#' @importFrom coda gelman.diag mcmc.list mcmc
#' @importFrom methods as callGeneric new is
#' @importFrom utils read.table write.table
update_drtmpt <- function(fit, 
                          data,
                          n.iter = 1000,
                          n.thin = 1,
                          Irep = 1000, 
                          flags = NULL,
                          control = NULL) {

  # CHECKS
  if (Irep < 1) stop("\"Irep\" must be larger than or equal to one.")
  if (n.thin < 1) stop("\"n.thin\" must be larger than or equal to one.")
  
  if (Irep %% n.thin != 0) stop("\"Irep\" must be a multiple of \"n.thin\".")
  if (n.iter %% (Irep/n.thin) != 0) stop("\"n.iter\" must be a multiple of \"Irep\" / \"n.thin\".")
  if (n.iter < Irep/n.thin) stop("\"n.iter\" must be greater or equal to \"Irep\" / \"n.thin\" = ", Irep*n.thin, ".")
  
  # if (Rhat_max < 1) stop("\"Rhat_max\" must be larger than or equal to one.")
  
  if (!is.null(flags) && !is.list(flags)) stop("\"flags\" must be a list or NULL.")
  if (!is.null(control) && !is.list(control)) stop("\"control\" must be a list or NULL.")
  
  if (!is.null(flags)) {
    if ("old_label" %in% names(flags)) if (!is.logical(flags$old_label)) stop("\"flags$old_label\" must either be TRUE or FALSE.")
    if ("indices" %in% names(flags)) if (!is.logical(flags$indices)) stop("\"flags$indices\" must either be TRUE or FALSE.")
    if ("loglik" %in% names(flags)) if (!is.logical(flags$loglik)) stop("\"flags$loglik\" must either be TRUE or FALSE.")
  }
  
  if (!is.null(control)) {
    if ("maxthreads" %in% control) if (!is.numeric(control$maxthreads)) stop("\"control$maxthreads\" must be numeric.")
  }
  
  
  # SET FLAGS
  if (!is.list(flags)) flags <- fit$specs$flags
  if (!"old_label" %in% names(flags)) {
    old_label <- FALSE
  } else {
    old_label <- fit$specs$flags$old_label
  }
  if (!"indices" %in% names(flags)) flags$indices <- fit$specs$flags$indices
  if (!"loglik" %in% names(flags)) flags$loglik <- fit$specs$flags$loglik
  flags$random_init <- FALSE
  
  
  # SET CONTROLS
  if (!is.list(control)) control <- fit$specs$control
  if (!"maxthreads" %in% names(control)) control$maxthreads <- fit$specs$control$maxthreads
  control$maxtreedepth1_3 <- fit$specs$control$maxtreedepth1_3
  control$maxtreedepth4 <- fit$specs$control$maxtreedepth4
  
  
  # PREPARE DATA
  keep_data_path <- FALSE
  data_frame <- NULL
  if (is.data.frame(data)) {
    temp_data <- to_drtmpt_data(data, model)
    data_frame <- temp_data$data
    if("transformation" %in% names(temp_data)) {
      transformation <- temp_data$transformation
    } else transformation <- list()
  } else if (is.character(data)) {
    temp_data <- to_drtmpt_data(read.table(file = data, header = TRUE), model)
    data_frame <- temp_data$data
    if("transformation" %in% names(temp_data)) {
      transformation <- temp_data$transformation
    } else transformation <- list()
    keep_data_path <- TRUE
  } else if (inherits(data, "drtmpt_data")) {
    data_frame <- data$data
    if("transformation" %in% names(data)) {
      transformation <- data$transformation
    } else transformation <- list()
  }
  rm(data); gc()
  
  data_elmnts <- c("subj", "group", "tree", "cat", "rt")
  if (!all(data_elmnts %in% names(data_frame))) stop("\"data\" must contain \"", data_elmnts[which(!(data_elmnts %in% names(data_frame)))[1]], "\".")
  if (!all(data_elmnts == names(data_frame))) {
    df <- data_frame
    data_frame <- df[,match(data_elmnts, names(data_frame))]
  }
  if (any(is.na(data_frame)) || min(data_frame) < 0) stop("All values in \"data\" need to be larger than or equal to zero and must not be NA.")
  if (max(data_frame[,1]+1) != length(unique(data_frame[,1]))) stop("\"max(data$subj+1)\" must equal \"length(unique(data$subj))\". You might want to use to_rtmpt_data().")
  if (max(data_frame[,2]+1) != length(unique(data_frame[,2]))) stop("\"max(data$group+1)\" must equal \"length(unique(data$group))\". You might want to use to_rtmpt_data().")
  if (max(data_frame[,3]+1) != length(unique(data_frame[,3]))) stop("\"max(data$tree+1)\" must equal \"length(unique(data$tree))\". You might want to use to_rtmpt_data().")
  if (max(data_frame[,4]+1) != length(unique(data_frame[,4]))) stop("\"max(data$cat+1)\" must equal \"length(unique(data$cat))\". You might want to use to_rtmpt_data().")
  
  
  # PREPARE PATHS
  # Inputs
  model_path <- gsub("\\\\", "/", tempfile(pattern = "model", tmpdir = tempdir(), fileext = ".info"))
  mdl_path <- gsub("\\\\", "/", tempfile(pattern = "model", tmpdir = tempdir(), fileext = ".txt"))
  data_path <- gsub("\\\\", "/", tempfile(pattern = "data", tmpdir = tempdir(), fileext = ".txt"))
  write.table(x = data_frame, file = data_path, sep = " ", row.names = FALSE, col.names = TRUE)
  # Outputs
  raus_path <- gsub("\\\\", "/", tempfile(pattern = "raus", tmpdir = tempdir(), fileext = ".out"))
  writeLines(text = fit$specs$raus, con = raus_path)
  ll_path <- gsub("\\\\", "/", tempfile(pattern = "loglik", tmpdir = tempdir(), fileext = ".out"))
  cont_path <- gsub("\\\\", "/", tempfile(pattern = "continue", tmpdir = tempdir(), fileext = ".out"))
  writeLines(text = fit$specs$continue, con = cont_path)
  means_path <- gsub("\\\\", "/", tempfile(pattern = "means", tmpdir = tempdir(), fileext = ".out"))
  tests_path <- gsub("\\\\", "/", tempfile(pattern = "tests", tmpdir = tempdir(), fileext = ".out"))
  rand_path <- gsub("\\\\", "/", tempfile(pattern = "random", tmpdir = tempdir(), fileext = ".out"))
  writeBin(fit$specs$rand, con = rand_path)
  
  
  # PREPARE INFOFILE
  model <- fit$specs$model
  infofile <- try(get_infofile(model, mdl_txt = mdl_path, mdl_info = model_path))
  if(is(infofile, "try-error")) stop("problem with S4 routines in makin infofile.\n")
  
  
  # Prior Parameter
  prior_params <- fit$specs$prior_params
  
  
  # PROCESS NAMES AND NUMBER
  proc_names <- names(model$params$threshold)
  nprocs <- length(proc_names)
  
  
  # KERN2FREE
  K2F <- integer(3*nprocs)
  cntk2f <- integer(3)
  consts <- rep(NaN, 3*nprocs)
  comp <- rep(FALSE, 3*nprocs)
  cntcomp <- integer(3)
  for (i in 1:3) {
    cntk2f[i] <- 0
    for (j in 1:nprocs) {
      tmp <- NULL
      if (i == 1) {
        tmp <- model$params$threshold[1, j]
      } else if (i == 2) {
        tmp <- model$params$driftrate[1, j]
      } else if (i == 3) {
        tmp <- model$params$startpoint[1, j]
      }
      
      ind <- (i-1)*nprocs+j
      if(is.na(tmp) | is.numeric(tmp)) {
        K2F[ind] <- cntk2f[i]
        cntk2f[i] <- cntk2f[i] + 1
        if (is.na(tmp)) {
          comp[(cntk2f[i]-1)*3 + i] <- TRUE
        }
        if (is.numeric(tmp)) {
          consts[(cntk2f[i]-1)*3 + i] <- tmp
        }
      } else if(tmp %in% proc_names) {
        K2F[ind] <- K2F[which(proc_names == as.character(tmp)) + (i>=2)*nprocs + (i==3)*nprocs]
      }
    }
  }
  
  
  
  # PATHS ARGUMENT
  CHAR1 <- as.character(c(Data = data_path,
                          Model = infofile,
                          Raus = raus_path,
                          LogLik = ll_path,
                          Continue = cont_path,
                          Means = means_path,
                          Tests = tests_path,
                          Random = rand_path,
                          tmpdir = paste0(tempdir(), "/")))
  
  
  # MCMC SETTINGS ARGUMENT
  n.chains <- fit$specs$n.chains
  INTEGER1 <- as.integer(c(IREP = Irep,
                           PHASE1 = fit$specs$n.phase1,
                           PHASE2 = fit$specs$n.phase2,
                           THIN = n.thin,
                           NOTHREADS = n.chains,
                           SAMPLE_SIZE = n.chains*n.iter,
                           MAXTHREADS = control$maxthreads,
                           nKERN = dim(model$responses)[1],
                           nPROCS = nprocs,
                           respno = length(unique(model$responses$MAP)),
                           cat2resp = model$responses$MAP))
  REAL1 <- c(Rhat_max = NaN)
  
  
  # FLAGS ARGUMENT
  BOOL1 <- as.integer(c(DIC_FL = flags$indices,
                        LOG_LIK_FL = (flags$loglik),
                        ML_INIT_FL = !flags$random_init))
  
  
  # PRIORS ARGUMENT
  INTEGER2 <- as.integer(c(degf = prior_params$delta_df))
  REAL2 <- c(PRIOR = prior_params$prec_epsilon,
             etat = prior_params$SIGMA_Corr_eta,
             taut = prior_params$SIGMA_SD_rho,
             etar = prior_params$GAMMA_Corr_eta,
             taur = prior_params$GAMMA_SD_rho,
             mu_prior = prior_params$delta_mu,
             rsd = prior_params$delta_scale,
             prioralpha = prior_params$Omega2_alpha,
             priorbeta = prior_params$Omega2_beta)
  
  
  # HMC OPTIONS ARGUMENT
  INTEGER3 <- as.integer(c(maxtreedepth1_3 = control$maxtreedepth1_3,
                           maxtreedepth4 = control$maxtreedepth4))
  
  
  # CONTINUATION ARGUMENT
  INTEGER4 <- as.integer(c(goon = TRUE,
                           ADDITION = n.chains*n.iter))
  
  
  # CONSTANTS AND EQUALIZATION
  INTEGER5 <- as.integer(c(kern2free = K2F,
                           comp = comp,
                           ifree = cntk2f))
  REAL3 <- c(consts = consts)
  
  
  # C++ FUNCTION DRTMPT CALL
  out <- .Call("drtmpt_fit", CHAR1, INTEGER1, REAL1, BOOL1,
               INTEGER2, REAL2,
               INTEGER3, INTEGER4,
               REAL3, INTEGER5)
  file.remove(data_path)

  
  # DATA INFORMATION
  data_info <- list(Nsubj = length(unique(data_frame[,1])),
                    Ngroups = length(unique(data_frame[,2])),
                    Nthreshold = sum(is.na(model$params$threshold[1,])),
                    Ndriftrate = sum(is.na(model$params$driftrate[1,])),
                    Nstartpoint = sum(is.na(model$params$startpoint[1,])),
                    Nresps = length(unique(model$responses$MAP)),
                    thresh_string = names(model$params$threshold[which(is.na(model$params$threshold[1,]))]),
                    drift_string = names(model$params$driftrate[which(is.na(model$params$driftrate[1,]))]),
                    start_string = names(model$params$startpoint[which(is.na(model$params$startpoint[1,]))]),
                    transformation = transformation)
  
  
  # PREPARE OUTPUT LIST
  drtmpt <- list()
  tmp_samples <- make_mcmc_list_d(file = out$pars_samples, infofile = infofile,
                                  Nchains = n.chains, Nsamples = n.iter,
                                  data_info = data_info, keep = old_label)
  drtmpt$samples <- mcmc.list(Map(function(o, a) mcmc(rbind(o, a)), fit$samples, tmp_samples))
  remove(tmp_samples)
  
  file.remove(raus_path)
  
  
  # DIAGNOSTICS
  drtmpt$diags <- get_diags_d(diag_file = tests_path, data_info = data_info, keep = old_label, DIC = flags$indices)
  drtmpt$diags$R_hat <- gelman.diag(drtmpt$samples, multivariate = FALSE)
  
  
  # SPECS
  infos <- readinfofile(infofile)
  drtmpt$specs <- list(model = model, n.chains = n.chains, n.iter = n.iter+fit$specs$n.iter, 
                       n.phase1 = fit$specs$n.phase1, n.phase2 = fit$specs$n.phase2,
                       n.thin = n.thin, n.groups = data_info$Ngroups, n.subj = data_info$Nsubj, Irep = Irep,
                       Rhat_max = NaN, prior_params = prior_params, infolist = infos, call = match.call())
  if(exists("transformation")) {
    drtmpt$specs$transformation <- transformation
  }
  
  
  # # COMPUTE WAIC AND LOO-IC
  if (all(out$loglik == 0)) {
    out$loglik <- NULL
  } else {
    if (flags$indices) {
      suppressWarnings( temp <- get_indices_x(out$loglik, n.chains, n.iter+fit$specs$n.iter, data_frame) )
      drtmpt$indices <- list(WAIC = temp$WAIC, LOO = temp$LOO)
      rm(temp)
    }
    if (flags$loglik) {
      drtmpt$loglik <- out$loglik
      if (file.exists(ll_path)) file.remove(ll_path)
    }
  }
  
  
  # CONTINUATION
  drtmpt$specs$continue <- readLines(con = cont_path)
  drtmpt$specs$raus <- fit$specs$raus
  drtmpt$specs$rand <- readBin(con = rand_path, what = "raw", n = 1000000)
  file.remove(cont_path)
  suppressWarnings(file.remove(raus_path))
  file.remove(rand_path)
  
  
  # SUMMARY
  drtmpt$summary <- writeSummaryDRTMPT(x = drtmpt, keep = old_label)
  
  
  # CLEAN-UP
  file.remove(means_path)
  file.remove(tests_path)
  file.remove(infofile)
  file.remove(mdl_path)
  suppressWarnings(file.remove(model_path))
  
  
  # OUTPUT
  class(drtmpt) <- "drtmpt_fit"
  gc()
  return(drtmpt)
  
}