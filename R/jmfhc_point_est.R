#' Obtain the Point Estimation of the Joint Cure Model with Flexible Hazards Ratios Over Time
#' @description Fits a cure model accounting for longitudinal data and flexible patterns of hazard ratios over time
#' and returns the point estimation.
#' Tied failure times are dealt with Breslow approximation.
#' @param data a data.frame with no missing values contains observed survival time, event status, any continuous variable(s) and/or dummy variable(s) of categorical variable(s).
#' The categorical variables need to be converted into dummy variables before fitting the function.
#' @param event_time the variable name of observed survival time in the data.
#' @param event_status the variable name of event status in the data.
#' @param id the variable name of patient id in the data.
#' @param beta_variable the variable name(s) defined as long-term covariate(s) in the model, which cannot be NULL.
#' @param gamma_variable the variable name(s) defined as short-term covariate(s) in the model. By default gamma_variable = NULL.
#' If there is no defined short-term covariate, the joint model becomes a proportional hazards cure submodel.
#' @param fu_measure the variable name of repeatedly meausred outcome (e.g., a biomarker).
#' @param fu_time_variable the variable name of measurement times for the longitudinal measurements.
#' @param max.int maximum number of iterations. The default is 200.
#' @param no_cores the number of cores used during the estimation. The default is 9.
#' @param absconverge.par absolute difference \eqn{(current-previous)} for regression parameters as the convergence criteria. The default is 1e-3.
#' @param relconverge.par relative difference \eqn{[(current-previous)/previous]} for regression parameters as the convergence criteria. The default is 2e-3.
#' @param absconverge.F0t the average of absolute differences for \eqn{F_0(t)} as the convergence criteria. The default is 2e-3.
#' @param relconverge.F0t the average of relative differences for \eqn{F_0(t)} as the convergence criteria. The default is 5e-3.
#'
#' @return a list containing results of the fit. The following items coef, iter, dat_baseline, dat_long, re_cov, and am_random_effects are returned.
#'     \item{coef}{estimated regression parameters (beta and gamma)}
#'     \item{iter}{the number of iterations used to complete the point estimation}
#'     \item{dat_baseline}{the final data at baseline including each patient's mean random effects from adaptive Markov algorithm, estimated \eqn{F_0(t)}, and \eqn{f_0(t)}}
#'     \item{dat_long}{the final data with patients' all records including each patient's mean random effects from adaptive Markov algorithm, estimated \eqn{F_0(t)}, and \eqn{f_0(t)}}
#'     \item{re_cov}{estimated covariance matrix of random effects}
#'     \item{am_random_effects}{a list containing each patient's final adaptive Markov chain}
#' @export
#'
#' @examples data(jmfhc_dat)
#' result_coef <- jmfhc_point_est(data=jmfhc_dat, event_time="event.time", event_status="event", id="patient.id",
#'                                beta_variable="trt", gamma_variable="trt",
#'                                fu_measure="measure", fu_time_variable="mes.times")
#' result_coef <- jmfhc_point_est(data=jmfhc_dat, event_time="event.time", event_status="event", id="patient.id",
#'                                beta_variable="trt",
#'                                fu_measure="measure", fu_time_variable="mes.times")
#' result$coef
#' @import dplyr
#' @import lme4
#' @import survival
#' @importFrom zoo na.locf
#' @import MASS
#' @import mvtnorm
#' @import MHadaptive
#' @import parallel
#' @import foreach
#' @import doSNOW
#'
jmfhc_point_est <- function(data, event_time, event_status, id,
                            beta_variable, gamma_variable=NULL,
                            fu_measure, fu_time_variable,
                            max.int=200, no_cores=9,
                            absconverge.par=1e-3, relconverge.par=2e-3,
                            absconverge.F0t=2e-3, relconverge.F0t=5e-3){

  cat("Point estimation starts.","\n")

  colnames(data)[(colnames(data)==id)] <- "id"

  dat_base <- data %>%
    group_by(id) %>%
    slice_head(n = 1)
  n <- nrow(dat_base)

  length_lmm_var <- length(fu_time_variable)+1
  random_effects <- paste0("re",seq(length_lmm_var))

  Z_var <- c(beta_variable,random_effects)
  if (!is.null(gamma_variable)){
    X_var <- gamma_variable
  } else{
    X_var <- beta_variable[1]
  }

  #### Step 1: initial values ####
  # regression parameters in the longitudinal submodel
  lmem <- as.formula(paste(fu_measure,"~",fu_time_variable,"+(",fu_time_variable,"|id)"))
  lmem <- lmer(lmem,REML = TRUE,data=data)
  ss <- summary(lmem)
  prior_Sigma <- ss$varcor$id[seq(length_lmm_var),seq(length_lmm_var)]
  fixed <- matrix(coef(ss)[,1],nrow=length_lmm_var)
  sigma_error <- suppressWarnings(as.vector(attr(ss$varcor, "sc")))
  # predict the random effects
  dat_base[,random_effects] <-ranef(lmem)$id
  data <- merge(data,dat_base[,c(random_effects,"id")],by="id")

  # regression parameters in the cure submodel
  dat_base <- as.matrix(dat_base[order(as.matrix(dat_base[,event_time])),])
  t <- dat_base[dat_base[,event_status]==1,event_time]
  formula.exp <- paste0("Surv(",event_time,",",event_status,")~", paste(Z_var, collapse=" + "))
  dat_base <- as.data.frame(dat_base)
  initial <- coxph(as.formula(formula.exp), data = dat_base)
  H0 <- basehaz(initial, centered=FALSE) # considered ties using breslow estimator
  H0 <- H0[H0[, 2] %in% unique(t), ]
  #initial values of beta
  beta <- initial$coefficients
  #initial value of beta_0
  beta0 <- log(H0$hazard[H0$time==max(H0$time)])
  #initial values of F0(t)
  base_cdf <- sapply(1:nrow(H0),function(x) H0$hazard[x]/H0$hazard[length(H0$hazard)])
  #initial values of f0(t)
  base_f <- c(base_cdf[1],sapply(2:length(base_cdf), function(x) base_cdf[x]-base_cdf[x-1]))
  base <- data.frame(H0$time,base_cdf,base_f)
  names(base)[1] <- event_time
  dat_base[dat_base[,event_status]==1,c("base_cdf","base_f")] <- merge(dat_base[dat_base[,event_status]==1,],base,by=event_time)[,c("base_cdf","base_f")]
  dat_base$base_f[dat_base$event==0] <- 0
  # fill in F0(t) for censored patients
  cen <- which(dat_base$event==0)
  death <- which(dat_base$event==1)
  #check the beginning
  if (sum(cen<death[1])!=0){
    pos <- cen[cen<death[1]]
    dat_base$base_cdf[pos] <- 0
  }
  dat_base$base_cdf <- na.locf(dat_base$base_cdf)
  data <- merge(data,dat_base[,c("base_f","base_cdf","id")],by="id")
  # initial values of gamma
  gamma_loglik_init <- function(gamma,data,beta,beta0){

    data_temp <- as.matrix(data)
    Z_i <- data_temp[,Z_var]
    X_i <- data_temp[,X_var]
    betaZ <- Z_i%*%beta
    gammaX <- X_i%*%as.matrix(gamma,nrow=length(gamma))
    event_i <- data_temp[,"event"]

    failure_part <- event_i*(beta0+betaZ+gammaX+(exp(gammaX)-1)*log(data_temp[,"base_cdf"])+log(data_temp[,"base_f"]))
    failure_part[which(is.na(failure_part))] <- 0
    LL <- sum(failure_part-exp(beta0)*exp(betaZ)*(data_temp[,"base_cdf"]^exp(gammaX)))

    return(LL)
  }
  if (!is.null(gamma_variable)){
    if (length(gamma_variable)>1){
      gamma_est <- optim(par=rep(0,length(gamma_variable)),fn=gamma_loglik_init,data=dat_base,beta=beta,beta0=beta0,
                         control=list("fnscale"=-1), hessian=F, method = "Nelder-Mead")
    }else{
      gamma_est <- optim(par=0,fn=gamma_loglik_init,data=dat_base,beta=beta,beta0=beta0,
                         control=list("fnscale"=-1,maxit=10000), hessian=F, method = "Brent",
                         lower=-10,upper=10)
    }

    if (gamma_est$convergence==0){
      gamma<- gamma_est$par
    }
  } else {
    gamma <- rep(0,length(X_var))
  }
  dat_base <- dat_base[order(dat_base$id),]

  #### 1st iteration ####
  #### Step 2: generate random effects by using Adaptive Markov algorithm ####
  target <- function(pars,data,betas,gammas,sigma_prior,fixed,sigma_error){
    prior <- dmvnorm(pars, mean = rep(0,length(pars)), sigma = sigma_prior,log=T)
    h <- ifelse(data[1,event_status]==1,(betas%*%c(1,unlist(data[1,beta_variable]),pars)+gammas*data[1,X_var]+(exp(gammas*data[1,X_var])-1)*log(data$base_cdf[1])+log(data$base_f[1])),0)
    S <- -exp(betas%*%c(1,unlist(data[1,beta_variable]),pars))*(data$base_cdf[1]^(exp(gammas*data[1,X_var])))
    D_i <- matrix(c(rep(1,nrow(data)),data[,fu_time_variable]),byrow=F,ncol=length_lmm_var)
    f_biomarker <- dmvnorm(data[,fu_measure],mean=D_i%*%fixed+D_i%*%pars,sigma=diag(rep(sigma_error^2,nrow(data))),log=T)
    return(prior+h+S+f_biomarker)
  }
  cl <- makeCluster(no_cores)
  registerDoSNOW(cl)
  mcmc_r <- foreach(i=dat_base$id,.packages = c('MASS','mvtnorm','MHadaptive')) %dopar%
    {Metro_Hastings(li_func=target, pars=as.numeric(dat_base[dat_base$id==i,random_effects]),
                    iterations = 10000, burn_in=1000,
                    par_names=random_effects,
                    data=data[data$id==i,],
                    betas = c(beta0,beta),
                    gammas = gamma,
                    sigma_prior = prior_Sigma,
                    sigma_error = sigma_error,
                    fixed = fixed)
    }
  stopCluster(cl)

  # thinning by 5
  cl <- makeCluster(no_cores)
  registerDoSNOW(cl)
  thin_sample <- foreach(i=seq(n),.packages = c("MASS",'mvtnorm','MHadaptive')) %dopar%
    { mcmc_thin(mcmc_r[[i]]) }
  stopCluster(cl)

  # add random effects into the baseline data set
  sample <- lapply(seq(n),function(x) thin_sample[[x]]$trace)
  proposal_sigma <- lapply(seq(n),function(x) thin_sample[[x]]$prop_sigma)
  expect <- sapply(seq(n),function(x) apply(sample[[x]],2,mean))
  dat_base[,random_effects] <- t(expect)
  # add random effects into the long data set
  remove <- which(colnames(data)%in%random_effects)
  data <- data[,-remove]
  data <- merge(data,dat_base[,c("id",random_effects)],by="id")

  #### Step 3: Estimate longitudinal parameters ####
  # estimate fixed effects
  est_fixed <- function(data){ # long version of data
    D_full <- matrix(c(rep(1,nrow(data)),data[,fu_time_variable]), nrow=nrow(data),ncol=length_lmm_var)
    D_full_re <- sapply(dat_base$id,function(x)
      matrix(c(rep(1,nrow(data[data$id==x,])),data[data$id==x,fu_time_variable]),byrow=F,ncol=length_lmm_var)%*%t(data[data$id==x,random_effects][1,]))
    b_star <- matrix(data[,fu_measure]-unlist(D_full_re), nrow=nrow(data),ncol=1)
    fixed <- solve(t(D_full)%*%D_full)%*%t(D_full)%*%b_star
    return(fixed)
  }
  fixed_k1 <- est_fixed(data=data)

  # estimate sd of error term
  D_full_fixed <- sapply(dat_base$id,function(x)
    matrix(c(rep(1,nrow(data[data$id==x,])),data[data$id==x,fu_time_variable]),byrow=F,ncol=length_lmm_var)%*%fixed_k1)
  cl <- makeCluster(no_cores)
  registerDoSNOW(cl)
  quadratic <-  foreach(i=seq(nrow(sample[[1]])),.combine="c")%dopar%{
    each_D_full_alpha <- unlist(sapply(dat_base$id, function(x)
      matrix(c(rep(1,nrow(data[data$id==x,])),data[data$id==x,fu_time_variable]),byrow=F,ncol=length_lmm_var)%*%sample[[which(unique(data$id)==x)]][i,]))
    b_star_invQ_b_star <- data[,fu_measure]-each_D_full_alpha-unlist(D_full_fixed)
    t(b_star_invQ_b_star)%*%b_star_invQ_b_star
  }
  stopCluster(cl)
  check <- data %>% count(id)
  sigma_error_k1 <- sqrt((1/sum(check[,2]))*mean(quadratic))

  # estimate covariance matrix of random effects
  sum_mean_alpha_sqrd <- foreach(x=seq(dat_base$id),.combine="+")%do%{t(sample[[x]]/nrow(sample[[x]]))%*%sample[[x]]}
  prior_Sigma_k1 <- (1/n)*sum_mean_alpha_sqrd

  prior_Sigma_diag_k1 <- sqrt(diag(prior_Sigma_k1))
  prior_Sigma_rho_k1 <- c()
  for (i in seq(nrow(prior_Sigma_k1))){
    for (j in seq(nrow(prior_Sigma_k1))[-seq(i)]){
      prior_Sigma_rho_k1 <- c(prior_Sigma_rho_k1,
                              prior_Sigma_k1[i,j]/(prior_Sigma_diag_k1[i]*prior_Sigma_diag_k1[j]))
    }
  }

  #### Step 4: Estimate survival parameters except F0(t) ####
  # estimate beta and gamma
  beta_gamma_comp <- function(data,par,beta0){

    beta_gamma_loglik <- function(par){

      beta <- matrix(par[seq(Z_var)],nrow=1,ncol=length(Z_var))
      if (!is.null(gamma_variable)){
        gamma <- matrix(par[(length(Z_var)+1):length(par)],nrow=1,ncol=length(X_var))
      } else{
        gamma <- matrix(rep(0,length(X_var)),nrow=1,ncol=length(X_var))
      }
      re <- par[(length(beta_variable)+1):length(Z_var)]

      data$mu_exp_re <- foreach(i=seq(n),.combine = "c")%do%{
        exp_value <- exp(sample[[which(unique(data$id)==data$id[i])]]%*%re)
        return(mean(exp_value))
      }

      data_temp <- as.matrix(data)
      Z_i <- data_temp[,Z_var]
      X_i <- data_temp[,X_var]
      betaZ <-Z_i%*%t(beta)
      exp_betaZ <- exp(Z_i[,beta_variable]%*%matrix(beta[seq(length(beta_variable))],nrow=length(beta_variable)))*data_temp[,"mu_exp_re"]
      gammaX <- X_i%*%t(gamma)
      event_i <- data_temp[,"event"]

      failure_part <- event_i*(beta0+betaZ+gammaX+(exp(gammaX)-1)*log(data_temp[,"base_cdf"])+log(data_temp[,"base_f"]))
      failure_part[which(is.na(failure_part))] <- 0
      LL <- sum(failure_part-exp(beta0)*exp_betaZ*(data_temp[,"base_cdf"]^exp(gammaX)))

      return(LL)
    }

    beta_gamma_est <-optim(par=par,fn=beta_gamma_loglik,method = "Nelder-Mead",control=list("fnscale"=-1,maxit=10000), hessian=F)

    if (beta_gamma_est$convergence==0){
      beta_gamma_k <- beta_gamma_est$par
      beta <- matrix(beta_gamma_k[seq(Z_var)],nrow=1,ncol=length(Z_var))
      if (!is.null(gamma_variable)){
        gamma <- matrix(beta_gamma_k[(length(Z_var)+1):length(beta_gamma_k)],nrow=1,ncol=length(X_var))
      } else {
        gamma <- matrix(rep(0,length(X_var)),nrow=1,ncol=length(X_var))
      }
    }else{
      beta <- "beta does not converge"
      gamma <- "gamma does not converge"
    }

    return(list(beta,gamma))
  }
  beta_gamma_k1 <- beta_gamma_comp(data=dat_base,par=c(beta,gamma),beta0=beta0)
  beta_k1 <- beta_gamma_k1[[1]]
  gamma_k1 <- beta_gamma_k1[[2]]

  # add expectation of exp(random effects)
  re <- beta_k1[(length(beta_variable)+1):length(Z_var)]
  dat_base$mu_exp_re <- foreach(i=seq(n),.combine = "c")%do%{
    exp_value <- exp(sample[[which(unique(dat_base$id)==dat_base$id[i])]]%*%re)
    return(mean(exp_value))
  }

  # estimate beta0
  beta0_comp <- function(data,beta,gamma){
    data_temp <- as.matrix(data)
    X_i <- data_temp[,X_var]
    gammaX <- X_i%*%t(gamma)
    mu_exp_re <- data_temp[,"mu_exp_re"]
    exp_betaZ <- exp(data_temp[,beta_variable]%*%matrix(beta[seq(length(beta_variable))],nrow=length(beta_variable)))*mu_exp_re
    beta_0 <- log(sum(data_temp[,event_status])/(t(exp_betaZ)%*%(data_temp[,"base_cdf"]^(exp(gammaX)))))
    return(as.vector(beta_0))
  }
  beta0_k1 <- beta0_comp(data=dat_base, beta = beta_k1, gamma = gamma_k1)

  #### Step 5: estimate F0(t) in the cure submodel ####
  f_est <- function(data,gamma,beta,beta_0){

    data <- as.matrix(data)
    data <- data[order(data[,event_time]),]
    # failure times
    t <- data[data[,event_status]==1,event_time]
    # at risk
    exp_betaZ <- lapply(unique(t),function(x)
      exp(data[data[,event_time]>=x,beta_variable]%*%matrix(beta[seq(length(beta_variable))],nrow=length(beta_variable)))*data[data[,event_time]>=x,"mu_exp_re"])
    gammaX <- lapply(unique(t),function(x) data[data[,event_time]>=x,X_var]%*%gamma)
    base_F_i <- lapply(unique(t), function(x) data[data[,event_time]>=x,"base_cdf"])

    # at risk and fail in the future
    gammaX_failure <- lapply(unique(t),function(x) data[data[,event_time]>=x & data[,event_status]==1,X_var]%*%gamma)
    base_F_i_failure <- lapply(unique(t), function(x) data[data[,event_time]>=x & data[,event_status]==1,"base_cdf"])

    # number of tied events
    w_i <- sapply(unique(t),function(x) sum(data[,event_status]==1 &data[,event_time]==x))

    # denominator
    deno <- sapply(seq(unique(t)),function(x)
      sum(exp(beta_0)*exp_betaZ[[x]]*exp(gammaX[[x]])*((base_F_i[[x]])^(exp(gammaX[[x]])-1)))
      -sum((exp(gammaX_failure[[x]])-1)/base_F_i_failure[[x]]))

    #obtain the value of alpha
    alpha_est <- function(alpha){
      sum(w_i/(deno+alpha))-1
    }
    alpha <- "error"
    times <- 2
    while(alpha=="error" & times < 5){
      range_root <- lapply(seq(-10^times,10^times),function(x) c(x,(x+1)))
      alpha_result <- sapply(seq(range_root),function(i)tryCatch(uniroot(alpha_est,range_root[[i]],extendInt = "yes")$root,error=function(x) return("error")))
      alpha <- ifelse(sum(alpha_result=="error")<length(range_root),"root","error")
      if (alpha=="root"){
        alpha_select <- unique(round(as.numeric(alpha_result[alpha_result!="error"]),4))

        for (i in seq(alpha_select)){
          check <- as.numeric(alpha_select[i])
          if (abs(sum(w_i/(deno+check))-1) <1e-4 & sum(w_i/(deno+check)<0)==0){
            alpha_select[i] <- check
          } else{
            alpha_select[i] <- "error"
          }
        }
        if (sum(alpha_select!="error")>0){
          alpha <- as.numeric(alpha_select[alpha_select!="error"])
          if (length(alpha)>1){
            alpha <- alpha[1]
          }
        }else{alpha <- "error"}
      }
      times <- times +1
    }

    new_base_f <- w_i/(deno+alpha)
    new_base_cdf <- cumsum(new_base_f)
    event_base_f <- unlist(sapply(seq(length(unique(t))),function(x) rep(new_base_f[x],sum(data[data[,event_status]==1,event_time]==unique(t)[x]))))
    event_base_cdf <- unlist(sapply(seq(length(unique(t))),function(x) rep(new_base_cdf[x],sum(data[data[,event_status]==1,event_time]==unique(t)[x]))))

    remove <- which(colnames(data)%in%c("base_cdf","base_f"))
    data <- as.data.frame(data[,-remove])
    data$base_f <- NA
    data$base_cdf <- NA
    data$base_f[data$event==1] <- event_base_f
    data$base_f[data$event==0] <- 0
    data$base_cdf[data$event==1] <- event_base_cdf

    cen <- which(data$event==0)
    death <- which(data$event==1)
    #check the beginning
    if (sum(cen<death[1])!=0){
      pos <- cen[cen<death[1]]
      data$base_cdf[pos] <- 0
    }
    data$base_cdf <- na.locf(data$base_cdf)
    data <- data[order(data$id),]

    return(list(new_base_cdf,data,alpha))
  }
  base_f_0 <- f_est(data=dat_base,gamma=gamma_k1,beta=beta_k1,beta_0=beta0_k1)
  base_cdf_k1 <- base_f_0[[1]]
  dat_base <- base_f_0[[2]]
  remove <- which(colnames(data)%in%c("base_f","base_cdf"))
  data <- data[,-remove]
  data <- merge(data,dat_base[,c("id","base_f","base_cdf")],by="id")

  #### 2nd iteration ####
  #### Step 2 ####
  cl <- makeCluster(no_cores)
  registerDoSNOW(cl)
  mcmc_r <- foreach(i=dat_base$id,.packages = c("MASS",'mvtnorm','MHadaptive')) %dopar%
    {Metro_Hastings(li_func=target,pars=as.numeric(dat_base[dat_base$id==i,random_effects]),
                    iterations = 10000, burn_in=1000,prop_sigma=proposal_sigma[[which(dat_base$id==i)]],
                    par_names=random_effects,
                    data=data[data$id==i,],
                    betas=c(beta0_k1,beta_k1),
                    gammas = gamma_k1,
                    sigma_prior = prior_Sigma_k1,
                    sigma_error=sigma_error_k1,
                    fixed = fixed_k1)
    }
  stopCluster(cl)



  # thinning
  cl <- makeCluster(no_cores)
  registerDoSNOW(cl)
  thin_sample <- foreach(i=seq(n),.packages = c("MASS",'mvtnorm','MHadaptive')) %dopar%
    { mcmc_thin(mcmc_r[[i]]) }
  stopCluster(cl)
  # add random effects into the baseline data set
  sample <- lapply(seq(n),function(x) thin_sample[[x]]$trace)
  proposal_sigma <- lapply(seq(n),function(x) thin_sample[[x]]$prop_sigma)
  expect <- sapply(seq(n),function(x) apply(sample[[x]],2,mean))
  dat_base[,random_effects] <- t(expect)
  # add random effects into the long data set
  remove <- which(colnames(data)%in%random_effects)
  data <- data[,-remove]
  data <- merge(data,dat_base[,c("id",random_effects)],by="id")

  #### Step 3 ####
  # estimate fixed effects
  fixed_k2 <- est_fixed(data=data)
  # estimate sd of error term
  D_full_fixed <- sapply(dat_base$id,function(x)
    matrix(c(rep(1,nrow(data[data$id==x,])),data[data$id==x,fu_time_variable]),byrow=F,ncol=length_lmm_var)%*%fixed_k2)
  cl <- makeCluster(no_cores)
  registerDoSNOW(cl)
  quadratic <-  foreach(i=seq(nrow(sample[[1]])),.combine="c")%dopar%{
    each_D_full_alpha <- unlist(sapply(dat_base$id, function(x)
      matrix(c(rep(1,nrow(data[data$id==x,])),data[data$id==x,fu_time_variable]),byrow=F,ncol=length_lmm_var)%*%sample[[which(unique(data$id)==x)]][i,]))
    b_star_invQ_b_star <- data[,fu_measure]-each_D_full_alpha-unlist(D_full_fixed)
    t(b_star_invQ_b_star)%*%b_star_invQ_b_star
  }
  stopCluster(cl)
  check <- data %>% count(id)
  sigma_error_k2 <- sqrt((1/sum(check[,2]))*mean(quadratic))

  # estimate covariance matrix of random effects
  sum_mean_alpha_sqrd <- foreach(x=seq(dat_base$id),.combine="+")%do%{t(sample[[x]]/nrow(sample[[x]]))%*%sample[[x]]}
  prior_Sigma_k2 <- (1/n)*sum_mean_alpha_sqrd

  prior_Sigma_diag_k2 <- sqrt(diag(prior_Sigma_k2))
  prior_Sigma_rho_k2 <- c()
  for (i in seq(nrow(prior_Sigma_k2))){
    for (j in seq(nrow(prior_Sigma_k2))[-seq(i)]){
      prior_Sigma_rho_k2 <- c(prior_Sigma_rho_k2,
                              prior_Sigma_k2[i,j]/(prior_Sigma_diag_k2[i]*prior_Sigma_diag_k2[j]))
    }
  }

  #### Step 4 ####
  # estimate beta and gamma
  beta_gamma_k2 <- beta_gamma_comp(data=dat_base,par=c(beta_k1,gamma_k1),beta0=beta0_k1)
  beta_k2 <- beta_gamma_k2[[1]]
  gamma_k2 <- beta_gamma_k2[[2]]
  # add expectation of exp(random effects)
  re <- beta_k2[(length(beta_variable)+1):length(Z_var)]
  dat_base$mu_exp_re <- foreach(i=seq(n),.combine = "c")%do%{
    exp_value <- exp(sample[[which(unique(dat_base$id)==dat_base$id[i])]]%*%re)
    return(mean(exp_value))
  }
  # estimate beta0
  beta0_k2 <- beta0_comp(data=dat_base, beta = beta_k2, gamma = gamma_k2)

  #### Step 5 ####
  base_f_0 <- f_est(data=dat_base,gamma=gamma_k2,beta=beta_k2,beta_0=beta0_k2)
  base_cdf_k2 <- base_f_0[[1]]
  dat_base <- base_f_0[[2]]
  remove <- which(colnames(data)%in%c("base_f","base_cdf"))
  data <- data[,-remove]
  data <- merge(data,dat_base[,c("id","base_f","base_cdf")],by="id")

  #### Keep doing iterations until convergence ####
  iteration <- 2
  beta_0_diff <- min(c(abs((beta0_k2-beta0_k1)/beta0_k1) >= relconverge.par, abs(beta0_k2-beta0_k1) >= absconverge.par))
  beta_diff <- min(c(sum(abs((beta_k2-beta_k1)/beta_k1) >= relconverge.par),sum((abs(beta_k2-beta_k1) >= absconverge.par))))
  base_cdf_diff <- min(c(mean(abs((base_cdf_k2-base_cdf_k1)/base_cdf_k1)) >= relconverge.F0t, mean(abs((base_cdf_k2-base_cdf_k1))) >= absconverge.F0t))
  fixed_diff <- min(c(sum(abs((fixed_k2-fixed_k1)/fixed_k1) >= relconverge.par),sum(abs(fixed_k2-fixed_k1) >= absconverge.par)))
  sigma_error_diff <- min(c(sum(abs((sigma_error_k2-sigma_error_k1)/sigma_error_k1) >= relconverge.par),sum(abs(sigma_error_k2-sigma_error_k1) >= absconverge.par)))
  prior_sigma_diff <- min(c(sum(abs((prior_Sigma_diag_k2-prior_Sigma_diag_k1)/prior_Sigma_diag_k1) >= relconverge.par),sum(abs(prior_Sigma_diag_k2-prior_Sigma_diag_k1) >= absconverge.par)))
  rho_diff <- min(c(sum(abs((prior_Sigma_rho_k2-prior_Sigma_rho_k1)/prior_Sigma_rho_k1) >= relconverge.par),sum(abs(prior_Sigma_rho_k2-prior_Sigma_rho_k1) >= absconverge.par)))

  if(!is.null(gamma_variable)){
    gamma_diff <- min(c(sum(abs((gamma_k2-gamma_k1)/gamma_k1) >= relconverge.par),sum((abs(gamma_k2-gamma_k1) >= absconverge.par))))
    deter <- ( beta_0_diff!=0  | beta_diff !=0 | gamma_diff!=0 | base_cdf_diff!=0 | fixed_diff!=0 |sigma_error_diff!=0 |prior_sigma_diff!=0
               |rho_diff!=0 )
  }else{
    deter <- ( beta_0_diff!=0  | beta_diff !=0 | base_cdf_diff!=0 | fixed_diff!=0 |sigma_error_diff!=0 |prior_sigma_diff!=0 |rho_diff!=0 )
  }

  while(deter & iteration < max.int){

    iteration  <- iteration + 1
    beta0_k1 <- beta0_k2
    beta_k1 <- beta_k2
    gamma_k1 <- gamma_k2
    base_cdf_k1 <- base_cdf_k2
    fixed_k1 <- fixed_k2
    sigma_error_k1 <- sigma_error_k2
    prior_Sigma_k1 <- prior_Sigma_k2
    prior_Sigma_rho_k1 <- prior_Sigma_rho_k2

    #### Step 2 ####
    cl <- makeCluster(no_cores)
    registerDoSNOW(cl)
    mcmc_r <- foreach(i=dat_base$id,.packages = c("MASS",'mvtnorm','MHadaptive')) %dopar%
      {Metro_Hastings(li_func=target,pars=as.numeric(dat_base[dat_base$id==i,random_effects]),
                      iterations = 10000, burn_in=1000,prop_sigma=proposal_sigma[[which(dat_base$id==i)]],
                      par_names=random_effects,
                      data=data[data$id==i,],
                      betas=c(beta0_k1,beta_k1),
                      gammas = gamma_k1,
                      sigma_prior = prior_Sigma_k1,
                      sigma_error=sigma_error_k1,
                      fixed = fixed_k1)
      }
    stopCluster(cl)

    # thinning
    cl <- makeCluster(no_cores)
    registerDoSNOW(cl)
    thin_sample <- foreach(i=seq(n),.packages = c("MASS",'mvtnorm','MHadaptive')) %dopar%
      { mcmc_thin(mcmc_r[[i]]) }
    stopCluster(cl)
    # add random effects into the baseline data set
    sample <- lapply(seq(n),function(x) thin_sample[[x]]$trace)
    proposal_sigma <- lapply(seq(n),function(x) thin_sample[[x]]$prop_sigma)
    expect <- sapply(seq(n),function(x) apply(sample[[x]],2,mean))
    dat_base[,random_effects] <- t(expect)
    # add random effects into the long data set
    remove <- which(colnames(data)%in%random_effects)
    data <- data[,-remove]
    data <- merge(data,dat_base[,c("id",random_effects)],by="id")

    #### Step 3 ####
    # estimate fixed effects
    fixed_k2 <- est_fixed(data=data)
    # estimate sd of error term
    D_full_fixed <- sapply(dat_base$id,function(x)
      matrix(c(rep(1,nrow(data[data$id==x,])),data[data$id==x,fu_time_variable]),byrow=F,ncol=length_lmm_var)%*%fixed_k2)
    cl <- makeCluster(no_cores)
    registerDoSNOW(cl)
    quadratic <-  foreach(i=seq(nrow(sample[[1]])),.combine="c")%dopar%{
      each_D_full_alpha <- unlist(sapply(dat_base$id, function(x)
        matrix(c(rep(1,nrow(data[data$id==x,])),data[data$id==x,fu_time_variable]),byrow=F,ncol=length_lmm_var)%*%sample[[which(unique(data$id)==x)]][i,]))
      b_star_invQ_b_star <- data[,fu_measure]-each_D_full_alpha-unlist(D_full_fixed)
      t(b_star_invQ_b_star)%*%b_star_invQ_b_star
    }
    stopCluster(cl)
    check <- data %>% count(id)
    sigma_error_k2 <- sqrt((1/sum(check[,2]))*mean(quadratic))

    # estimate covariance matrix of random effects
    sum_mean_alpha_sqrd <- foreach(x=seq(dat_base$id),.combine="+")%do%{t(sample[[x]]/nrow(sample[[x]]))%*%sample[[x]]}
    prior_Sigma_k2 <- (1/n)*sum_mean_alpha_sqrd

    prior_Sigma_diag_k2 <- sqrt(diag(prior_Sigma_k2))
    prior_Sigma_rho_k2 <- c()
    for (i in seq(nrow(prior_Sigma_k2))){
      for (j in seq(nrow(prior_Sigma_k2))[-seq(i)]){
        prior_Sigma_rho_k2 <- c(prior_Sigma_rho_k2,
                                prior_Sigma_k2[i,j]/(prior_Sigma_diag_k2[i]*prior_Sigma_diag_k2[j]))
      }
    }

    #### Step 4 ####
    # estimate beta and gamma
    beta_gamma_k2 <- beta_gamma_comp(data=dat_base,par=c(beta_k1,gamma_k1),beta0=beta0_k1)
    beta_k2 <- beta_gamma_k2[[1]]
    gamma_k2 <- beta_gamma_k2[[2]]
    # add expectation of exp(random effects)
    re <- beta_k2[(length(beta_variable)+1):length(Z_var)]
    dat_base$mu_exp_re <- foreach(i=seq(n),.combine = "c")%do%{
      exp_value <- exp(sample[[which(unique(dat_base$id)==dat_base$id[i])]]%*%re)
      return(mean(exp_value))
    }
    # estimate beta0
    beta0_k2 <- beta0_comp(data=dat_base, beta = beta_k2, gamma = gamma_k2)

    #### Step 5 ####
    base_f_0 <- f_est(data=dat_base,gamma=gamma_k2,beta=beta_k2,beta_0=beta0_k2)
    base_cdf_k2 <- base_f_0[[1]]
    dat_base <- base_f_0[[2]]
    remove <- which(colnames(data)%in%c("base_f","base_cdf"))
    data <- data[,-remove]
    data <- merge(data,dat_base[,c("id","base_f","base_cdf")],by="id")

    ## Keep doing iterations until convergence
    beta_0_diff <- min(c(abs((beta0_k2-beta0_k1)/beta0_k1) >= relconverge.par, abs(beta0_k2-beta0_k1) >= absconverge.par))
    beta_diff <- min(c(sum(abs((beta_k2-beta_k1)/beta_k1) >= relconverge.par),sum((abs(beta_k2-beta_k1) >= absconverge.par))))
    base_cdf_diff <- min(c(mean(abs((base_cdf_k2-base_cdf_k1)/base_cdf_k1)) >= relconverge.F0t, mean(abs((base_cdf_k2-base_cdf_k1))) >= absconverge.F0t))
    fixed_diff <- min(c(sum(abs((fixed_k2-fixed_k1)/fixed_k1) >= relconverge.par),sum(abs(fixed_k2-fixed_k1) >= absconverge.par)))
    sigma_error_diff <- min(c(sum(abs((sigma_error_k2-sigma_error_k1)/sigma_error_k1) >= relconverge.par),sum(abs(sigma_error_k2-sigma_error_k1) >= absconverge.par)))
    prior_sigma_diff <- min(c(sum(abs((prior_Sigma_diag_k2-prior_Sigma_diag_k1)/prior_Sigma_diag_k1) >= relconverge.par),sum(abs(prior_Sigma_diag_k2-prior_Sigma_diag_k1) >= absconverge.par)))
    rho_diff <- min(c(sum(abs((prior_Sigma_rho_k2-prior_Sigma_rho_k1)/prior_Sigma_rho_k1) >= relconverge.par),sum(abs(prior_Sigma_rho_k2-prior_Sigma_rho_k1) >= absconverge.par)))

    if(!is.null(gamma_variable)){
      gamma_diff <- min(c(sum(abs((gamma_k2-gamma_k1)/gamma_k1) >= relconverge.par),sum((abs(gamma_k2-gamma_k1) >= absconverge.par))))
      deter <- ( beta_0_diff!=0  | beta_diff !=0 | gamma_diff!=0 | base_cdf_diff!=0 | fixed_diff!=0 |sigma_error_diff!=0 |prior_sigma_diff!=0
                 |rho_diff!=0 )
    }else{
      deter <- ( beta_0_diff!=0  | beta_diff !=0 | base_cdf_diff!=0 | fixed_diff!=0 |sigma_error_diff!=0 |prior_sigma_diff!=0 |rho_diff!=0 )
    }
  }

  rho_num <- unlist(sapply(seq(length_lmm_var-1),function(x) paste0(x,seq(length_lmm_var)[-seq(x)])))

  if (!is.null(gamma_variable)){
    result.coef <- matrix(c(beta0_k2,beta_k2,gamma_k2,fixed_k2,prior_Sigma_diag_k2,prior_Sigma_rho_k2, sigma_error_k2),nrow=1)
    colnames(result.coef) <- c("beta_intercept",paste0("beta_",Z_var),
                               paste0("gamma_",X_var),
                               paste0("fixed_",seq(length_lmm_var)),
                               paste0("re_sd_",seq(length_lmm_var)),
                               paste0("re_rho_",rho_num),
                               "error_sd")
  }else{
    result.coef <- matrix(c(beta0_k2,beta_k2,fixed_k2,prior_Sigma_diag_k2,prior_Sigma_rho_k2,sigma_error_k2),nrow=1)
    colnames(result.coef) <- c("beta_intercept",paste0("beta_",Z_var),
                               paste0("fixed_",seq(length_lmm_var)),
                               paste0("re_sd_",seq(length_lmm_var)),
                               paste0("re_rho_",rho_num),
                               "error_sd")
  }

  remove <- which(colnames(dat_base)=="mu_exp_re")
  dat_base <- dat_base[,-remove]
  colnames(dat_base)[(colnames(dat_base)=="id")] <- id
  colnames(data)[(colnames(data)=="id")] <- id

  cat("Point estimation is done.","\n")

  return(list(coef=result.coef,iter=iteration,dat_baseline=dat_base,dat_long = data,
              re_cov=prior_Sigma_k2,
              am_random_effects=sample))

}
