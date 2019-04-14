
# load in libraries
library(lme4) # library to calcualte GLM
library(geeM) # library to calculate GEE
library(jomo) # library to do MMI

# Function 1:
myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

# Function 2:
missing_icc <- function(icc){
  pi <- 3.142
  res <- pi^2/(3*(1/icc-1))
  return(res)
}

# Function 3:
missing_per <- function(data){
  con <- sum(data[data$arm == 0, ]$r / dim(data[data$arm == 0, ])[1])
  trt <- sum(data[data$arm == 1, ]$r / dim(data[data$arm == 1, ])[1])
  res = list(con = con, trt = trt)
  return(res)
}

# Function 4:
expit <- function(x){y <- exp(x)/(1+exp(x));return(y)}

# Function 5:
BCVexch <- function(y,X,beta,alpha,phi,id,w){
  
  require(MASS)
  
  # Creates two vectors that have the start and end points for each cluster
  BEGINEND <- function(n){
    last <- cumsum(n)
    first <- last-n+1
    return(cbind(first,last))
  }
  
  # Score function
  SCORE <- function(beta,alpha,phi,y,X,n,p){
    U <- rep(0,p)
    UUtran <- Ustar <- matrix(0,p,p)
    locx <- BEGINEND(n)
    
    for(i in 1:length(n)){
      X_c <- X[locx[i,1]:locx[i,2],,drop=FALSE]
      y_c <- y[locx[i,1]:locx[i,2]]
      w_c <- w[locx[i,1]:locx[i,2]]
      
      U_c <- rep(0,p)
      Ustar_c <- matrix(0,p,p)
      mu_c <- 1/(1+exp(c(-X_c%*%beta)))
      
      C <- X_c*(mu_c*(1-mu_c))
      A <- y_c-mu_c
      D <- diag(w_c,nrow=length(w_c))
      INVR <- diag(1/(1-alpha),n[i])-matrix(alpha/((1-alpha)*(1-alpha+n[i]*alpha)),n[i],n[i])
      INVB <- diag(1/sqrt(mu_c*(1-mu_c)),n[i]) %*% INVR %*% diag(1/sqrt(mu_c*(1-mu_c)),n[i]) %*% D/phi
      
      U_c <- t(C)%*%INVB%*%A
      UUtran_c <- tcrossprod(U_c)
      Ustar_c <- t(C)%*%INVB%*%C
      U <- U+U_c
      UUtran <- UUtran+UUtran_c
      Ustar <- Ustar+Ustar_c
    }
    return(list(U=U,UUtran=UUtran,Ustar=Ustar))
  }
  
  # Creates bias-corrected covariance matrix of beta
  p <- ncol(X)
  n <- as.numeric(table(id))
  SCORE_RES <- SCORE(beta,alpha,phi,y,X,n,p)
  U <- SCORE_RES$U
  UUtran <- SCORE_RES$UUtran
  Ustar <- SCORE_RES$Ustar
  
  # Naive or Model-based estimator
  naive <- ginv(Ustar)
  
  # BC0 or usual Sandwich estimator     
  robust <- naive%*%UUtran%*%t(naive)
  
  # Bias-corrected variance
  Ustar_c_array <- UUtran_c_array <- array(0,c(p,p,length(n)))
  UUtran <- UUbc <- UUbc2 <- UUbc3 <- Ustar <- matrix(0,p,p)
  
  locx <- BEGINEND(n)
  
  for(i in 1:length(n)){
    X_c <- X[locx[i,1]:locx[i,2],,drop=FALSE]
    y_c <- y[locx[i,1]:locx[i,2]]
    w_c <- w[locx[i,1]:locx[i,2]]
    mu_c <- 1/(1+exp(c(-X_c%*%beta)))
    
    U_i <- U_c <- rep(0,p)
    Ustar_c <- matrix(0,p,p)
    
    # commands for beta
    C <- X_c*(mu_c*(1-mu_c))
    A <- y_c-mu_c
    D <- diag(w_c,nrow=length(w_c))
    INVR <- diag(1/(1-alpha),n[i])-matrix(alpha/((1-alpha)*(1-alpha+n[i]*alpha)),n[i],n[i])
    INVB <- diag(1/sqrt(mu_c*(1-mu_c)),n[i]) %*% INVR %*% diag(1/sqrt(mu_c*(1-mu_c)),n[i]) %*% D/phi
    U_i <- t(C)%*%INVB%*%A
    
    L_i <- ginv(diag(1,nrow=n[i]) - C%*%naive%*%t(C)%*%INVB) %*% A
    U_c <- t(C)%*%INVB%*%L_i
    
    Ustar_c <- t(C)%*%INVB%*%D%*%C
    Ustar <- Ustar+Ustar_c
    UUtran_c <- tcrossprod(U_i)
    UUtran <- UUtran+UUtran_c
    UUbc_c <- tcrossprod(U_c)
    UUbc <- UUbc+UUbc_c
    UUbc_ic <- tcrossprod(U_c,U_i)
    UUbc2 <- UUbc2+UUbc_ic
    
    Ustar_c_array[,,i] <- Ustar_c
    UUtran_c_array[,,i] <- UUtran_c
  }
  
  # calculating adjustment factor for BC3
  for(i in 1:length(n)){      
    Hi <- diag(1/sqrt(1-pmin(0.75,c(diag(Ustar_c_array[,,i]%*%naive)))))
    UUbc3 <- UUbc3+Hi%*%UUtran_c_array[,,i]%*%Hi
  }
  
  # BC1 or Variance estimator due to Kauermann and Carroll (2001);
  varKC <- naive%*%(UUbc2+t(UUbc2))%*%t(naive)/2
  
  # BC2 or Variance estimator due to Mancl and DeRouen (2001);
  varMD <- naive%*%UUbc%*%t(naive)
  
  # BC3 or Variance estimator due to Fay and Graubard (2001);
  varFG <- naive%*%UUbc3%*%t(naive)
  
  ########################################
  # Output
  # naive: naive or model-based var
  # robust: robust sandwich var
  # varMD: bias-corrected sandwich var due to Mancl and DeRouen (2001)
  # varKC: bias-corrected sandwich var due to Kauermann and Carroll (2001)
  # varFG: bias-corrected sandwich var due to Fay and Graubard (2001)
  ########################################
  return(list(naive=naive,robust=robust,varMD=varMD,varKC=varKC,varFG=varFG))
}

datagen_interact <- function(k, M, 
                             icc, iccm,
                             mux=0, varx,
                             mud=0,slope_trt,slope_con,
                             intercept_trt, intercept_con){
  ########################################
  # Input:
  # k: a number of clusters in each intervention arm
  # M: the mean value of cluster size
  # mux: the mean value of covariate x, default as 0
  # varx: the variance of covariate x
  # icc: the icc of generated dataset
  # mud: the mean value of the parameter delta
  # iccm: the icc of missingness
  # intercept: the intercept of the missingness generation function
  ########################################
  
  # total cluster number
  K <- 2*k 
  # cluster sizes
  m <- rpois(K,M) 
  # total individual number
  N <- sum(m) 
  #  indicator of intervention arm; i = 1 treated; i = 0 non treated
  i <- rep(rep(c(0,1),each=k),times=m) 
  
  # cluster id
  cluster <- rep(1:K,times=m)
  
  # the variance of delta, based on dataset ICC
  vard <- missing_icc(icc)
  delta <- rep(rnorm(K,mud,sqrt(vard)),times=m) # the random effect in generating the outcome
  
  # different covariates in treatment arm and control arm
  N_trt <- sum(i); N_con = N - sum(i)
  
  x <- rnorm(N,mux,sqrt(varx))
  x_trt <- x[i == 1] # the covariate in treated 
  x_con <- x[i == 0] # the covariate  in controled
  
  p <- expit(1 + 1.36 * i + x + delta) # the probability of the outcome
  y <- rbinom(N,1,p) # the outcome, 1 or 0.
  
  ## generating the missingness
  alpha <- rep(rnorm(K,0,sqrt(missing_icc(iccm))),times=m)
  
  # missingness generation function
  mis_trt <- expit(intercept_trt + 1 + slope_trt * x_trt + alpha[i==1])
  mis_con <- expit(intercept_con + 0 + slope_con * x_con + alpha[i==0])
  
  mis <- c(mis_con, mis_trt)
  
  # the missingness indicator. 
  # r = 1 missing; r = 0 not missing 
  r <- rbinom(N,1,mis)
  
  # the whole dataset
  res <- data.frame(y = y, arm = i, x = x,
                    cluster = cluster, delta = delta, 
                    mis = mis, r = r)
  return(res)
}

mypool <- function(mean0,sd0,num=5,J=50){
  ########################################
  # input: 
  # mean0: a vector of the values of estimated beta parameter, with length equals to the imputation time.
  # sd0: a vector of the standard deviations of the estimated beta parameter, with length equals to the imputation time.
  # num: impuation time
  # J: the number of clusters in each intervention arm 
  ########################################
  
  # count the times of NA in the input vector.
  na_times <- sum(is.na(mean0)) 
  # the number of imputations without NA. 
  num_actual <- num-na_times 
  
  # the MI estimate of the beta parameter
  m <- mean(mean0,na.rm=TRUE) 
  
  # estimate of average wihtin-imputation variance 
  # i.e. based on the SE^2 of the beta parameter from each fitted model
  W <- mean(sd0^2,na.rm=TRUE) 
  
  # estimate of between-imputation variance 
  # i.e. empirical SD of the point estimates of beta parameter
  B <- var(mean0,na.rm=TRUE) 
  
  # estimate of total variance 
  # i.e. will need to take the sqrt of this to use for inference in making confidence intervals etc.
  v_hat <- W+(1+1/num_actual)*B 
  v_hat_sqrt<-sqrt(v_hat)
  
  # Confindence interval
  # Testing based on naive normal distribution assumption
  l <- m-1.96*v_hat_sqrt
  u <- m+1.96*v_hat_sqrt
  
  # Testing based on standard results from MI literature
  # i.e. df of t distribution for testing based on standard results from MI literature
  df_denom <- (1+1/num_actual)*B
  df_part <- 1+W/df_denom
  df_t <- (num_actual-1)*df_part^2 # adjusted df of t distribution
  l_t <- m-qt(0.975,df_t)*v_hat_sqrt
  u_t <- m+qt(0.975,df_t)*v_hat_sqrt
  
  # Testing based on results from MMI literature, Barnard and Rubin (1999), 
  # df of t distribution for testing based on results in re adjustment for MMI feature
  
  #df for full data
  df_com <- 2*J - 2 
  
  # calculate the adjusted df based on the literature.
  parenthesis <- 1+df_denom*(1/W) 
  df_obs <- df_com*((df_com+1)/(df_com+3))*(1/parenthesis) 
  df_adj_t <- 1/(1/df_t + 1/df_obs) 
  
  # Print the results
  #print("Standard t df and Barnard/Rubin adjusted t df");
  #print(c(df_t, df_adj_t))
  #print("97.5% quantiles from standard t df and Barnard/Rubin adjusted t df");
  #print(c(qt(0.975,df_t), qt(0.975,df_adj_t)))
  
  # The confidence interval calculated based on df of standard t distribution and adjusted df of t distribution 
  l_t <- m-qt(0.975,df_t)*v_hat_sqrt
  u_t <- m+qt(0.975,df_t)*v_hat_sqrt
  l_adj_t <- m-qt(0.975,df_adj_t)*v_hat_sqrt
  u_adj_t <- m+qt(0.975,df_adj_t)*v_hat_sqrt
  
  return(list(mean=m,std=v_hat_sqrt,
              df_t=df_t,
              df_adj_t=df_adj_t))
}


simulations_mmi_interact = function(k_set, icc_set, iccm_set, Times, M = 50,
                                    varx, Nimp = 15,
                                    mis_trt_per,mis_con_per,res_name, print_name){
  
  for(k in k_set){
    for(icc in icc_set){
      for(iccm in iccm_set){
        
        MISTRT = c(); MISCON = c()
        est_ind <- c();est_ex <- c() # Estimates for independent and exchangeable
        std_ind <- c();std_ex <- c() # SEs for independent and exchangeable
        ro_ind <- c();ro_ex <- c() # the robust SEs calculated by fan's function, which is same with SEs calculated by geeM
        warn_ind <- c();warn_ex <- c() # whether there is warning
        # the SEs from each small sample correction methods
        MD_ind <- c();MD_ex <- c();
        KC_ind <- c();KC_ex <- c();
        FG_ind <- c();FG_ex <- c()
        df_t <- c(); df_adj_t <- c()
        
        for(times in 1:Times){
          print(print_name)
          print(paste('k',k,'icc',icc,'iccm',iccm,'times',times))
          
          intersequence = seq(-3,0,0.1)
          # test 1
          intercept_con = 1
          mis_trt = c(); mis_con = c()
          for(ii in intersequence){
            set.seed(times)
            intercept_trt = ii
            temp <- datagen_interact(k, M, 
                                     icc, iccm,
                                     mux=0, varx,
                                     mud=0,slope_trt,slope_con,
                                     intercept_trt, intercept_con)
            mis_trt = c(mis_trt, missing_per(temp)$trt)
            mis_con = c(mis_con, missing_per(temp)$con)
          }
          intercept_trt = intersequence[which(abs(mis_trt - mis_trt_per) == min(abs(mis_trt - mis_trt_per)))]
          
          mis_trt = c(); mis_con = c()
          for(ii in intersequence){
            set.seed(times)
            intercept_con = ii
            temp <- datagen_interact(k, M, 
                                     icc, iccm,
                                     mux=0, varx,
                                     mud=0,slope_trt,slope_con,
                                     intercept_trt, intercept_con)
            mis_trt = c(mis_trt, missing_per(temp)$trt)
            mis_con = c(mis_con, missing_per(temp)$con)
          }
          
          intercept_con = intersequence[which(abs(mis_con - mis_con_per) == min(abs(mis_con - mis_con_per)))]
          intercept_trt
          intercept_con
          
          set.seed(times)
          d1 <- datagen_interact(k, M, 
                                 icc, iccm,
                                 mux=0, varx,
                                 mud=0,slope_trt,slope_con,
                                 intercept_trt, intercept_con)
          
          MISTRT = c(MISTRT, missing_per(d1)$trt)
          MISCON = c(MISCON, missing_per(d1)$con)
          
          d3 <- d1
          d3$y <- ifelse(d3$r==1,NA,d3$y)
          d3$missing <- d3$r
          # 3.4 MMI; just copy Hossain's code
          data.miss <- d3
          y.cat <- data.frame(outcome=data.miss$y)  # data frame for response variables with missing values
          y.numcat <- c(2)                                 # number of levels in outcome variable
          clus <- data.frame(clus=data.miss$cluster)          # data frame for clusters
          
          nobs <- dim(data.miss)[1]
          x <- data.frame(intercept=rep(1,nobs),covariate=data.miss$x,group=data.miss$arm)
          
          # 3.4.1 run to generate Nimp full datasets
          # "imp" represents full datasets 
          imp <- jomo1rancat(Y.cat = y.cat, 
                             Y.numcat = y.numcat, X = x,
                             clus = clus,nburn = 100, 
                             nbetween = 25, nimp = Nimp, output = 0)
          
          # set empty results and save them later
          mmi_est_ind <- c();mmi_std_ind <- c();mmi_warn_ind <- c()
          mmi_est_ex <- c();mmi_std_ex <- c();mmi_warn_ex <- c()
          mmi_MD_ex <- c();mmi_KC_ex <- c();mmi_FG_ex <- c()
          mmi_FG_ind <- c();mmi_KC_ind <- c(); mmi_MD_ind <- c()
          mmi_ro_ex <- c();  mmi_ro_ind <- c()
          
          # 3.4.2 Analyze each of the full dataset with GEE
          for(i in 1:Nimp){
            temp <- imp[imp$Imputation==i,]
            rownames(temp) <- NULL
            temp$outcome <- as.numeric(temp$outcome)-1
            
            mmi_ind <- myTryCatch(geem(formula = outcome~group,
                                       id = clus , data = temp,
                                       family =  binomial("logit"),
                                       corstr = "independence"))
            
            mmi_ex <- myTryCatch(geem(formula = outcome~group,
                                      id = clus , data = temp,
                                      family =  binomial("logit"),
                                      corstr = "exchangeable"))
            
            ## Save the results
            if(is.null(mmi_ind$value)==0){
              # if there is no error when calculating the MMI
              # calculate the alpha, beta, phi
              phi1 <- mmi_ind$value$phi
              alpha1 <- mmi_ind$value$alpha
              beta1 <- mmi_ind$value$beta
              t1 <- beta1[2]
              y1 <- temp$outcome
              X1 <- cbind(rep(1,length(temp$outcome)),temp$group)
              w1 <- rep(1,length(temp$outcome))
              id1 <- temp$clus
              # input value into the BCVexch function to calculate the 
              # small sample correction of standard deviations
              correction1 <- myTryCatch(BCVexch(y1,X1,beta1,alpha1,phi1,id1,w1))
              std1 <- summary(mmi_ind$value)$se.robust[2]
              
              if(is.null(correction1$error)==1){
                # if there is no error when calculating the correction values
                correction1 <- correction1$value
                mmi_est_ind <- c(mmi_est_ind,t1)
                mmi_std_ind <- c(mmi_std_ind,std1)
                mmi_ro_ind <- c(mmi_ro_ind,sqrt(diag(correction1$robust))[2])
                mmi_MD_ind <- c(mmi_MD_ind,sqrt(diag(correction1$varMD))[2])
                mmi_KC_ind <- c(mmi_KC_ind,sqrt(diag(correction1$varKC))[2])
                mmi_FG_ind <- c(mmi_FG_ind,sqrt(diag(correction1$varFG))[2])
                if(is.null(correction1$error)==0){mmi_warn_ind <- c(mmi_warn_ind,times)}
                if(is.null(correction1$error)==1){mmi_warn_ind <- c(mmi_warn_ind,0)}
              }else{
                # if there is an error when calculating the correction values
                mmi_est_ind <- c(mmi_est_ind,NA)
                mmi_std_ind <- c(mmi_std_ind,NA)
                mmi_ro_ind <- c(mmi_ro_ind,NA)
                mmi_MD_ind <- c(mmi_MD_ind,NA)
                mmi_KC_ind <- c(mmi_KC_ind,NA)
                mmi_FG_ind <- c(mmi_FG_ind,NA)
                mmi_warn_ind <- c(mmi_warn_ind,times)
              }
            }
            
            if(is.null(mmi_ind$value)==1){
              # if there is an error when calculating the MMI
              mmi_est_ind <- c(mmi_est_ind,NA)
              mmi_std_ind <- c(mmi_std_ind,NA)
              mmi_ro_ind <- c(mmi_ro_ind,NA)
              mmi_MD_ind <- c(mmi_MD_ind,NA)
              mmi_KC_ind <- c(mmi_KC_ind,NA)
              mmi_FG_ind <- c(mmi_FG_ind,NA)
              mmi_warn_ind <- c(mmi_warn_ind,times)
            }
            
            if( is.null(mmi_ex$value)==0){
              # if there is no error in calculating the MMI
              phi2 <- mmi_ex$value$phi
              alpha2 <- mmi_ex$value$alpha
              beta2 <- mmi_ex$value$beta
              t2 <- beta2[2]
              y2 <- temp$outcome
              X2 <- cbind(rep(1,length(temp$outcome)),temp$group)
              w2 <- rep(1,length(temp$outcome))
              id2 <- temp$clus
              correction2 <- myTryCatch(BCVexch(y2,X2,beta2,alpha2,phi2,id2,w2))
              std2 <- summary(mmi_ex$value)$se.robust[2]
              if(is.null(correction2$error)==1){
                # if there is no error in calculating the correction
                correction2 <- correction2$value
                mmi_est_ex <- c(mmi_est_ind,t2)
                mmi_std_ex <- c(mmi_std_ind,std2)
                mmi_ro_ex <- c(mmi_ro_ex,sqrt(diag(correction2$robust))[2])
                mmi_MD_ex <- c(mmi_MD_ex,sqrt(diag(correction2$varMD))[2])
                mmi_KC_ex <- c(mmi_KC_ex,sqrt(diag(correction2$varKC))[2])
                mmi_FG_ex <- c(mmi_FG_ex,sqrt(diag(correction2$varFG))[2])
                if(is.null(mmi_ex$value)==1){mmi_warn_ex <- c(mmi_warn_ex,times)}
                if(is.null(mmi_ex$value)==0){mmi_warn_ex <- c(mmi_warn_ex,0)}
              }else{
                #if there is an error in calculating the correction
                mmi_est_ex <- c(mmi_est_ex,NA)
                mmi_std_ex <- c(mmi_std_ex,NA)
                mmi_ro_ex <- c(mmi_ro_ex,NA)
                mmi_MD_ex <- c(mmi_MD_ex,NA)
                mmi_KC_ex <- c(mmi_KC_ex,NA)
                mmi_FG_ex <- c(mmi_FG_ex,NA)
                mmi_warn_ex <- c(mmi_warn_ex,times)
              }
            }
            if(is.null(mmi_ex$value)==1){
              # if there is an error in calculating the MMI
              mmi_est_ex <- c(mmi_est_ex,NA)
              mmi_ro_ex <- c(mmi_ro_ex,NA)
              mmi_std_ex <- c(mmi_std_ex,NA)
              mmi_MD_ex <- c(mmi_MD_ex,NA)
              mmi_KC_ex <- c(mmi_KC_ex,NA)
              mmi_FG_ex <- c(mmi_FG_ex,NA)
              mmi_warn_ex <- c(mmi_warn_ex,times)
            }
          }
          
          # 3.4.3 Pool the results together
          # pool with orginal SE
          temp1 <- mypool(mmi_est_ind,mmi_std_ind,num=Nimp, J = k)
          temp2 <- mypool(mmi_est_ex,mmi_std_ex,num=Nimp, J = k)
          
          temp11 <- mypool(mmi_est_ind,mmi_ro_ind,num=Nimp,J = k)
          temp22 <- mypool(mmi_est_ex,mmi_ro_ex,num=Nimp,J = k)
          
          # pool with corrected SE
          temp3 <- mypool(mmi_est_ind,mmi_MD_ind,num=Nimp,J = k)
          temp4 <- mypool(mmi_est_ex,mmi_MD_ex,num=Nimp,J = k)
          
          temp5 <- mypool(mmi_est_ind,mmi_KC_ind,num=Nimp,J = k)
          temp6 <- mypool(mmi_est_ex,mmi_KC_ex,num=Nimp,J = k)
          
          temp7 <- mypool(mmi_est_ind,mmi_FG_ind,num=Nimp,J = k)
          temp8 <- mypool(mmi_est_ex,mmi_FG_ex,num=Nimp,J = k)
          
          # save results
          est_ind <- c(est_ind,temp1$mean)
          std_ind <- c(std_ind,temp1$std)
          ro_ind <- c(ro_ind,temp11$std)
          ro_ex <- c(ro_ex,temp22$std)
          est_ex <- c(est_ex,temp2$mean)
          std_ex <- c(std_ex,temp2$std)
          
          MD_ind <- c(MD_ind,temp3$std)
          MD_ex <- c(MD_ex,temp4$std)
          KC_ind <- c(KC_ind,temp5$std)
          KC_ex <- c(KC_ex,temp6$std)
          FG_ind <- c(FG_ind,temp7$std)
          FG_ex <- c(FG_ex,temp8$std)
          
          df_t <- c(df_t,temp1$df_t)
          df_adj_t <- c(df_adj_t,temp1$df_adj_t)
          
          warn_ind <- c(warn_ind,sum(mmi_warn_ind))
          warn_ex <- c(warn_ex,sum(mmi_warn_ex))
        }
        
        result=list(est_ind=est_ind,est_ex=est_ex,
                    std_ind=std_ind,std_ex=std_ex,
                    ro_ind=ro_ind,ro_ex=ro_ex,
                    warn_ind=warn_ind,warn_ex=warn_ex,
                    MD_ind=MD_ind,MD_ex=MD_ex,
                    KC_ind=KC_ind,KC_ex=KC_ex,
                    FG_ind=FG_ind,FG_ex=FG_ex,
                    df_t=df_t,
                    df_adj_t=df_adj_t,
                    mistrt = MISTRT,
                    miscon = MISCON)
        
        # print(result)
        names <- paste(res_name,k,icc,iccm,'.RData',sep='')
        save(result,file=names)
      }
    }
  }
}

          