
# function that gives variablenames/columnnames
labelnames <- function(data_info) {
  
  Gr <- data_info$Ngroups
  
  probs_string <- data_info$probs_string
  minus_string <- data_info$minus_string
  plus_string <- data_info$plus_string
  
  P <- length(probs_string)
  Pm <- length(minus_string)
  Pp <- length(plus_string)

  Ptot <- P+Pm+Pp
  Plam <- Ptot-P

  S <- data_info$Nsubj
  R <- data_info$Nresps
  
  
  # lable variables
  sig <- c("", "minus", "plus")
  Prime <- c("alpha_prime", "beta_prime")
  all_string <- c(probs_string, minus_string, plus_string)
  
  
  label <- c()
  index <- 0
  ## group mean probs
  for (gr in 1:Gr) {
    for (p in 1:P) label[index+p] <- paste0("mu_alpha_", probs_string[p], if (Gr>1) paste0("[", gr, "]"))
    index <- index + P
  }
  ## group mean rates
  for (gr in 1:Gr) {
    for (p in 1:Pm) label[index+p] <- paste0("exp_mu_beta_", minus_string[p], "_", sig[2], if (Gr>1) paste0("[", gr, "]"))
    index <- index + Pm
    for (p in 1:Pp) label[index+p] <- paste0("exp_mu_beta_", plus_string[p], "_", sig[3], if (Gr>1) paste0("[", gr, "]"))
    index <- index + Pp
  }
  # upper triangular of alpha' and beta' (SD and correlation)
  for (i in 1:( Ptot )) {
    for (j in i:( Ptot )) {
      label[index+j] <- paste0(ifelse(i<=P, Prime[1], Prime[2]), "_", all_string[i], 
                               ifelse(i<=P, sig[1], paste0("_", ifelse(i<=P+Pm, sig[2], sig[3]))), "__",
                               ifelse(j<=P, Prime[1], Prime[2]), "_", all_string[j],
                               ifelse(j<=P, sig[1], paste0("_", ifelse(j<=P+Pm, sig[2], sig[3]))))
    }
    index <- index + (Ptot-i)
  }
  index <- index + Ptot
  ## alpha'
  for (s in 1:S) label <- c(label, paste0(Prime[1], "_", probs_string, "[", s, "]"))
  index <- index + (S*P)
  ## beta'
  for (s in 1:S) label <- c(label, paste0(Prime[2], "_", all_string[(P+1):Ptot], "_", 
                                          c(rep(sig[2], Pm), rep(sig[3], Pp)), "[", s, "]"))
  index <- index + (S*Plam)
  ##  group mean encoding & response execution times
  for (gr in 1:Gr) {
    for (r in 1:R) label[index+r] <- paste0("mu_gamma", if(R>1) paste0("_R", r-1), if (Gr>1) paste0("[", gr, "]"))
    index <- index + R
  }
  ## omega^2
  label[index+1] <- paste("omega_square")
  index <- index + 1
  ## upper triangular of gamma' (SD and correlation)
  for (i in 1:R) {
    for (j in i:R) label[index+j] <- paste0("gamma_prime_", "R", i-1, "__", "gamma_prime_", "R", j-1)
    index <- index + (R-i)
  }
  index <- index + R
  ## gamma'
  for (s in 1:S) label <- c(label, paste0("gamma_prime", if(R>1) paste0("_R", (1:2)-1), "[", s, "]"))
  index <- index + (S*R)
  ## sigma square per participant
  for (s in 1:S) label[index+s] <- paste0("sigma_square[", s, "]")
  index <- index + S
  
  label <- c(label, "log_lik")
  index = index + 1
  
  return(label)

}


# data <- list(Nsubj = 47, Ngroups = 2, Nprobs = 3, Nminus = 3, Nplus = 4, Nresps = 1, epsilon = 1,
#              probs_string = c("A1", "C", "G"), minus_string = c("A2", "C", "G"),
#              plus_string = c("A1", "A2", "C", "G"))

make_mcmc_list <- function(file, infofile, Nchains, Nsamples, data_info) {
  
  # load necessary r-package and source
  # source("C:/Users/hartmann/Desktop/UNI_FREIBURG/U_Freiburg/CPP/DBDA2E-utilities.R")
  
  # read text file with chains
  if (is.character(file)) {
    temp <- as.vector(read.table(file=file,header=F,nrows=1))
    dt <- fread(file=file,skip=1)
  } else if (is.data.frame(file) || is.matrix(file)) {
    temp <- dim(file)
    dt <- as.data.table(file)
  }
  
  
  
  # specify parameters used in MCMC
  start <- 1
  end <- Nsamples
  npar <- temp[2]
  
  
  
  # generate MCMC-list for coda
  vec <- vector("list",Nchains)
  for (i in 1:Nchains) {
    vec[[i]] <- dt[(i-1)*end + start:end]
    vec[[i]] <- as.mcmc(vec[[i]],start=start,end=end,thin=1)
  }
  samples <- as.mcmc.list(vec,start=start,end=end,thin=1)
  rm(vec)
  
  
  # name chain columns
  # if(npar != length(varnames(samples))) stop("Number of columns in chains is wrong.")       ###---###
  # varstring <- system(paste("gawk 'END {print}' ", infofile), intern = TRUE)
  # var_char <- strsplit(varstring, "\\s+")[[1]]
  
  
  # label the samples
  varnames(samples) <- labelnames(data_info = data_info)
  
  return(samples)
  
}



