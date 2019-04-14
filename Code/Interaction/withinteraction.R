## with interaction

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


simulations_ipw_interact = function(k_set, icc_set, iccm_set, Times, M = 50,
                                    varx, 
                                    mis_trt_per,mis_con_per,res_name, print_name){
  
  
 
   for(k in k_set){
    for(icc in icc_set){
      for(iccm in iccm_set){
        
        
        no_interaction = c()
        
        MISTRT = c(); MISCON = c()
        
        # 4.1 set empty results and save for later
        # the results for true values.
        robust_true_ind <- c();MD_true_ind <- c();KC_true_ind <- c();FG_true_ind <- c()
        est_true_ind <- c();std_true_ind <- c()
        robust_true_ex <- c();MD_true_ex <- c();KC_true_ex <- c();FG_true_ex <- c()
        est_true_ex <- c();std_true_ex <- c()
        
        # the results for CRA-GEE
        robust_ucra_ind <- c();MD_ucra_ind <- c();KC_ucra_ind <- c();FG_ucra_ind <- c()
        est_ucra_ind <- c();std_ucra_ind <- c()
        robust_ucra_ex <- c();MD_ucra_ex <- c();KC_ucra_ex <- c();FG_ucra_ex <- c()
        est_ucra_ex <- c();std_ucra_ex <- c()
        
        # the results for A-CRA-GEE
        robust_cra_ind <- c();MD_cra_ind <- c();KC_cra_ind <- c();FG_cra_ind <- c()
        est_cra_ind <- c();std_cra_ind <- c()
        robust_cra_ex <- c();MD_cra_ex <- c();KC_cra_ex <- c();FG_cra_ex <- c()
        est_cra_ex <- c();std_cra_ex <- c()
        
        # the results for W-GEE
        robust_ipw_ex <- c();MD_ipw_ex <- c();KC_ipw_ex <- c();FG_ipw_ex <- c()
        est_ipw_ex <- c();std_ipw_ex <- c()
        robust_ipw_clu_ind <- c();MD_ipw_clu_ind <- c();KC_ipw_clu_ind <- c();FG_ipw_clu_ind <- c()
        est_ipw_clu_ind <- c();std_ipw_clu_ind <- c()
        
        # the results for CW-GEE
        robust_ipw_clu_ex <- c();MD_ipw_clu_ex <- c();KC_ipw_clu_ex <- c();FG_ipw_clu_ex <- c()
        est_ipw_clu_ex <- c();std_ipw_clu_ex <- c()
        robust_ipw_ind <- c();MD_ipw_ind <- c();KC_ipw_ind <- c();FG_ipw_ind <- c()
        est_ipw_ind <- c();std_ipw_ind <- c()
        
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
          
          # 4.4 calculate the weights for IPW
          # w1: the weights without consideration of cluster effects
          # w2: the weights with consideration of cluster effects
          
          w1 <- glm(missing ~x + arm + x*arm, data = d3,
                    family = binomial(link='logit'))
          w2 <- glmer(missing ~x + arm + x*arm +(1|cluster) , data = d3,
                      family = binomial(link='logit'))
          
          if(summary(w1)$coefficients[4,1] * summary(w2)$coefficients[4,1] == 0){
            no_interaction = c(no_interaction, times)
          }
          w1 <- predict(w1,type="response")  # get the weights value from the glm
          w2 <- predict(w2,type="response")  # get the weights value form the glmer
          
          d3$weight <- 1/w1
          d3$weight2 <- 1/w2
          
          d2 <- na.omit(d3)
          
          # 4.5 Calculation
          ### 4.5.1 True effect
          trues_ind <- myTryCatch(geem(formula = y~arm,id = cluster, data = d1,
                                       family =  binomial("logit"),
                                       corstr = "independence"))
          trues_ex <- myTryCatch(geem(formula = y~arm,id = cluster, data = d1,
                                      family =  binomial("logit"),
                                      corstr = "exchangeable"))
          
          ### 4.5.2 Unadjusted CRA
          ucra_ind <- myTryCatch(geem(formula = y~arm,id = cluster, data = d2,
                                      family =  binomial("logit"),
                                      corstr = "independence"))
          
          ucra_ex <- myTryCatch(geem(formula = y~arm,id = cluster, data = d2,
                                     family =  binomial("logit"),
                                     corstr = "exchangeable"))
          
          ### 4.5.3 Adjusted CRA
          cra_ind <- myTryCatch(geem(formula = y~x+arm,id = cluster, data = d2,
                                     family =  binomial("logit"),
                                     corstr = "independence"))
          cra_ex <- myTryCatch(geem(formula = y~x+arm,id = cluster, data = d2,
                                    family =  binomial("logit"),
                                    corstr = "exchangeable"))
          
          ### 4.5.4 IPW without cluster effects 
          ipw_ind <- myTryCatch(geem(formula = y~arm,id = cluster, data = d2,
                                     family =  binomial("logit"),
                                     weights = d2$weight,
                                     corstr = "independence"))
          ipw_ex <- myTryCatch(geem(formula = y~arm,id = cluster, data = d2,
                                    family =  binomial("logit"),
                                    weights = d2$weight,
                                    corstr = "exchangeable"))
          
          ### 4.5.5 IPW with cluster effects
          ipw_clu_ind <- myTryCatch(geem(formula = y~arm,id = cluster, data = d2,
                                         family =  binomial("logit"),
                                         weights = d2$weight2,
                                         corstr = "independence"))
          ipw_clu_ex <- myTryCatch(geem(formula=y~arm,id=cluster, data = d2,
                                        family =  binomial("logit"),
                                        weights = d2$weight2,
                                        corstr = "exchangeable"))
          
          # 4.6 Save the results
          ### true independent
          if(is.null(trues_ind$value)==0){
            if(trues_ind$value$converged==1){
              # if there is no error and the gee results are converged:
              # calculate the alpha, phi, beta value to calculate the corrected sd
              phi1 <- trues_ind$value$phi
              alpha1 <- trues_ind$value$alpha
              beta1 <- trues_ind$value$beta
              est_true_ind <- c(est_true_ind,beta1[2])
              std_true_ind <- c(std_true_ind,summary(trues_ind$value)$se.robust[2])
              w1 <- rep(1,dim(d1)[1])
              y1 <- d1$y
              X1 <- cbind(rep(1,dim(d1)[1]),d1$arm)
              id1 <- d1$cluster
              correction1 <- BCVexch(y1,X1,beta1,alpha1,phi1,id1,w1)
              robust_true_ind <- c(robust_true_ind,sqrt(diag(correction1$robust))[2])
              MD_true_ind <- c(MD_true_ind,sqrt(diag(correction1$varMD))[2])
              KC_true_ind <- c(KC_true_ind,sqrt(diag(correction1$varKC))[2])
              FG_true_ind <- c(FG_true_ind,sqrt(diag(correction1$varFG))[2])
            }else{
              # if the results are not converged:
              est_true_ind <- c(est_true_ind,NA)
              std_true_ind <- c(std_true_ind,NA)
              robust_true_ind <- c(robust_true_ind,NA)
              MD_true_ind <- c(MD_true_ind,NA)
              KC_true_ind <- c(KC_true_ind,NA)
              FG_true_ind <- c(FG_true_ind,NA)
            }
          }else{
            # if there is an error:
            est_true_ind <- c(est_true_ind,NA)
            std_true_ind <- c(std_true_ind,NA)
            robust_true_ind <- c(robust_true_ind,NA)
            MD_true_ind <- c(MD_true_ind,NA)
            KC_true_ind <- c(KC_true_ind,NA)
            FG_true_ind <- c(FG_true_ind,NA)
          }
          
          robust_true_ind;MD_true_ind;KC_true_ind;FG_true_ind
          est_true_ind;std_true_ind
          
          ## true exchangeable
          if(is.null(trues_ex$value)==0){
            if( trues_ex$value$converged==1){
              phi2 <- trues_ex$value$phi
              alpha2 <- trues_ex$value$alpha
              beta2 <- trues_ex$value$beta
              est_true_ex <- c(est_true_ex,beta2[2])
              std_true_ex <- c(std_true_ex,summary(trues_ex$value)$se.robust[2])
              w2 <- rep(1,dim(d1)[1])
              y2 <- d1$y
              X2 <- cbind(rep(1,dim(d1)[1]),d1$arm)
              id2 <- d1$cluster
              correction2 <- BCVexch(y2,X2,beta2,alpha2,phi2,id2,w2)
              robust_true_ex <- c(robust_true_ex,sqrt(diag(correction2$robust))[2])
              MD_true_ex <- c(MD_true_ex,sqrt(diag(correction2$varMD))[2])
              KC_true_ex <- c(KC_true_ex,sqrt(diag(correction2$varKC))[2])
              FG_true_ex <- c(FG_true_ex,sqrt(diag(correction2$varFG))[2])
            }else{
              est_true_ex <- c(est_true_ex,NA)
              std_true_ex <- c(std_true_ex,NA)
              robust_true_ex  <- c(robust_true_ex,NA)
              MD_true_ex <- c(MD_true_ex,NA)
              KC_true_ex <- c(KC_true_ex,NA)
              FG_true_ex <- c(FG_true_ex,NA)
            }
          }else{
            est_true_ex <- c(est_true_ex,NA)
            std_true_ex <- c(std_true_ex,NA)
            robust_true_ex  <- c(robust_true_ex,NA)
            MD_true_ex <- c(MD_true_ex,NA)
            KC_true_ex <- c(KC_true_ex,NA)
            FG_true_ex <- c(FG_true_ex,NA)
          }
          
          robust_true_ex;MD_true_ex;KC_true_ex;FG_true_ex
          est_true_ex;std_true_ex
          
          ### UCRA independent
          if(is.null(ucra_ind$value)==0){
            if( ucra_ind$value$converged==1){
              phi3 <- ucra_ind$value$phi
              alpha3 <- ucra_ind$value$alpha
              beta3 <- ucra_ind$value$beta
              est_ucra_ind <- c(est_ucra_ind,beta3[2])
              std_ucra_ind <- c(std_ucra_ind,summary(ucra_ind$value)$se.robust[2])
              w3 <- rep(1,dim(d2)[1])
              y3 <- d2$y
              X3 <- cbind(rep(1,dim(d2)[1]),d2$arm)
              id3 <- d2$cluster
              correction3 <- BCVexch(y3,X3,beta3,alpha3,phi3,id3,w3)
              robust_ucra_ind <- c(robust_ucra_ind,sqrt(diag(correction3$robust))[2])
              MD_ucra_ind <- c(MD_ucra_ind,sqrt(diag(correction3$varMD))[2])
              KC_ucra_ind <- c(KC_ucra_ind,sqrt(diag(correction3$varKC))[2])
              FG_ucra_ind <- c(FG_ucra_ind,sqrt(diag(correction3$varFG))[2])
            }else{
              est_ucra_ind <- c(est_ucra_ind,NA)
              std_ucra_ind <- c(std_ucra_ind,NA)
              robust_ucra_ind <- c(robust_ucra_ind,NA)
              MD_ucra_ind <- c(MD_ucra_ind,NA)
              KC_ucra_ind <- c(KC_ucra_ind,NA)
              FG_ucra_ind <- c(FG_ucra_ind,NA)
            }
          }else{
            est_ucra_ind <- c(est_ucra_ind,NA)
            std_ucra_ind <- c(std_ucra_ind,NA)
            robust_ucra_ind <- c(robust_ucra_ind,NA)
            MD_ucra_ind <- c(MD_ucra_ind,NA)
            KC_ucra_ind <- c(KC_ucra_ind,NA)
            FG_ucra_ind <- c(FG_ucra_ind,NA)
          }
          
          robust_ucra_ind;MD_ucra_ind;KC_ucra_ind;FG_ucra_ind
          est_ucra_ind;std_ucra_ind
          
          ## ucra exchangeable
          if(is.null(ucra_ex$value)==0 ){
            if(ucra_ex$value$converged==1){
              phi4 <- ucra_ex$value$phi
              alpha4 <- ucra_ex$value$alpha
              beta4 <- ucra_ex$value$beta
              est_ucra_ex <- c(est_ucra_ex,beta4[2])
              std_ucra_ex <- c(std_ucra_ex,summary(ucra_ex$value)$se.robust[2])
              w4 <- rep(1,dim(d2)[1])
              y4 <- d2$y
              X4 <- cbind(rep(1,dim(d2)[1]),d2$arm)
              id4 <- d2$cluster
              correction4 <- BCVexch(y4,X4,beta4,alpha4,phi4,id4,w4)
              robust_ucra_ex <- c(robust_ucra_ex,sqrt(diag(correction4$robust))[2])
              MD_ucra_ex <- c(MD_ucra_ex,sqrt(diag(correction4$varMD))[2])
              KC_ucra_ex <- c(KC_ucra_ex,sqrt(diag(correction4$varKC))[2])
              FG_ucra_ex <- c(FG_ucra_ex,sqrt(diag(correction4$varFG))[2])
              
            }else{
              est_ucra_ex <- c(est_ucra_ex,NA)
              std_ucra_ex <- c(std_ucra_ex,NA)
              robust_ucra_ex  <- c(robust_ucra_ex,NA)
              MD_ucra_ex <- c(MD_ucra_ex,NA)
              KC_ucra_ex <- c(KC_ucra_ex,NA)
              FG_ucra_ex <- c(FG_ucra_ex,NA)
            }
          }else{
            est_ucra_ex <- c(est_ucra_ex,NA)
            std_ucra_ex <- c(std_ucra_ex,NA)
            robust_ucra_ex  <- c(robust_ucra_ex,NA)
            MD_ucra_ex <- c(MD_ucra_ex,NA)
            KC_ucra_ex <- c(KC_ucra_ex,NA)
            FG_ucra_ex <- c(FG_ucra_ex,NA)
          }
          
          robust_ucra_ex;MD_ucra_ex;KC_ucra_ex;FG_ucra_ex
          est_ucra_ex;std_ucra_ex
          
          ## cra independent
          if(is.null(cra_ind$value)==0){
            if(cra_ind$value$converged==1){
              phi5 <- cra_ind$value$phi
              alpha5 <- cra_ind$value$alpha
              beta5 <- cra_ind$value$beta
              est_cra_ind <- c(est_cra_ind,beta5[3])
              std_cra_ind <- c(std_cra_ind,summary(cra_ind$value)$se.robust[2])
              w5 <- rep(1,dim(d2)[1])
              y5 <- d2$y
              X5 <- cbind(rep(1,dim(d2)[1]),d2$x,d2$arm)
              id5 <- d2$cluster
              correction5 <- BCVexch(y5,X5,beta5,alpha5,phi5,id5,w5)
              robust_cra_ind <- c(robust_cra_ind,sqrt(diag(correction5$robust))[2])
              MD_cra_ind <- c(MD_cra_ind,sqrt(diag(correction5$varMD))[2])
              KC_cra_ind <- c(KC_cra_ind,sqrt(diag(correction5$varKC))[2])
              FG_cra_ind <- c(FG_cra_ind,sqrt(diag(correction5$varFG))[2])
            }else{
              est_cra_ind <- c(est_cra_ind,NA)
              std_cra_ind <- c(std_cra_ind,NA)
              robust_cra_ind <- c(robust_cra_ind,NA)
              MD_cra_ind <- c(MD_cra_ind,NA)
              KC_cra_ind <- c(KC_cra_ind,NA)
              FG_cra_ind <- c(FG_cra_ind,NA)
            }
          }else{
            est_cra_ind <- c(est_cra_ind,NA)
            std_cra_ind <- c(std_cra_ind,NA)
            robust_cra_ind <- c(robust_cra_ind,NA)
            MD_cra_ind <- c(MD_cra_ind,NA)
            KC_cra_ind <- c(KC_cra_ind,NA)
            FG_cra_ind <- c(FG_cra_ind,NA)
          }
          
          robust_cra_ind;MD_cra_ind;KC_cra_ind;FG_cra_ind
          est_cra_ind;std_cra_ind
          
          ## cra exchangeable
          if(is.null(cra_ex$value)==0){
            if(cra_ex$value$converged==1){
              phi6 <- cra_ex$value$phi
              alpha6 <- cra_ex$value$alpha
              beta6 <- cra_ex$value$beta
              est_cra_ex <- c(est_cra_ex,beta6[3])
              std_cra_ex <- c(std_cra_ex,summary(cra_ex$value)$se.robust[2])
              w6 <- rep(1,dim(d1)[1])
              y6 <- d2$y
              X6 <- cbind(rep(1,dim(d2)[1]),d2$x,d2$arm)
              id6 <- d2$cluster
              correction6 <- BCVexch(y6,X6,beta6,alpha6,phi6,id6,w6)
              robust_cra_ex <- c(robust_cra_ex,sqrt(diag(correction6$robust))[2])
              MD_cra_ex <- c(MD_cra_ex,sqrt(diag(correction6$varMD))[2])
              KC_cra_ex <- c(KC_cra_ex,sqrt(diag(correction6$varKC))[2])
              FG_cra_ex <- c(FG_cra_ex,sqrt(diag(correction6$varFG))[2])
            }else{
              est_cra_ex <- c(est_cra_ex,NA)
              std_cra_ex <- c(std_cra_ex,NA)
              robust_cra_ex  <- c(robust_cra_ex,NA)
              MD_cra_ex <- c(MD_cra_ex,NA)
              KC_cra_ex <- c(KC_cra_ex,NA)
              FG_cra_ex <- c(FG_cra_ex,NA)
            }
          }else{
            est_cra_ex <- c(est_cra_ex,NA)
            std_cra_ex <- c(std_cra_ex,NA)
            robust_cra_ex  <- c(robust_cra_ex,NA)
            MD_cra_ex <- c(MD_cra_ex,NA)
            KC_cra_ex <- c(KC_cra_ex,NA)
            FG_cra_ex <- c(FG_cra_ex,NA)
          } 
          
          robust_cra_ex;MD_cra_ex;KC_cra_ex;FG_cra_ex
          est_cra_ex;std_cra_ex
          
          ### IPW independent
          if(is.null(ipw_ind$value)==0){
            if(ipw_ind$value$converged==1){
              phi7 <- ipw_ind$value$phi
              alpha7 <- ipw_ind$value$alpha
              beta7 <- ipw_ind$value$beta
              est_ipw_ind <- c(est_ipw_ind,beta7[2])
              std_ipw_ind <- c( std_ipw_ind,summary(ipw_ind$value)$se.robust[2])
              w7 <- d2$weight
              y7 <- d2$y
              X7 <- cbind(rep(1,dim(d2)[1]),d2$arm)
              id7 <- d2$cluster
              correction7 <- BCVexch(y7,X7,beta7,alpha7,phi7,id7,w7)
              robust_ipw_ind <- c(robust_ipw_ind,sqrt(diag(correction7$robust))[2])
              MD_ipw_ind <- c(MD_ipw_ind,sqrt(diag(correction7$varMD))[2])
              KC_ipw_ind <- c(KC_ipw_ind,sqrt(diag(correction7$varKC))[2])
              FG_ipw_ind <- c(FG_ipw_ind,sqrt(diag(correction7$varFG))[2])
            }else{
              est_ipw_ind <- c(est_ipw_ind,NA)
              std_ipw_ind <- c(std_ipw_ind,NA)
              robust_ipw_ind <- c(robust_ipw_ind,NA)
              MD_ipw_ind <- c(MD_ipw_ind,NA)
              KC_ipw_ind <- c(KC_ipw_ind,NA)
              FG_ipw_ind <- c(FG_ipw_ind,NA)
            }
          }else{
            est_ipw_ind <- c(est_ipw_ind,NA)
            std_ipw_ind <- c(std_ipw_ind,NA)
            robust_ipw_ind <- c(robust_ipw_ind,NA)
            MD_ipw_ind <- c(MD_ipw_ind,NA)
            KC_ipw_ind <- c(KC_ipw_ind,NA)
            FG_ipw_ind <- c(FG_ipw_ind,NA)
          } 
          
          robust_ipw_ind;MD_ipw_ind;KC_ipw_ind;FG_ipw_ind
          est_ipw_ind;std_ipw_ind
          
          ## ipw exchangeable
          if(is.null(ipw_ex$value)==0){
            if(ipw_ex$value$converged==1){
              phi8 <- ipw_ex$value$phi
              alpha8 <- ipw_ex$value$alpha
              beta8 <- ipw_ex$value$beta
              est_ipw_ex <- c(est_ipw_ex,beta8[2])
              std_ipw_ex <- c(std_ipw_ex,summary(ipw_ex$value)$se.robust[2])
              w8 <- d2$weight
              y8 <- d2$y
              X8 <- cbind(rep(1,dim(d2)[1]),d2$arm)
              id8 <- d2$cluster
              correction8 <- BCVexch(y8,X8,beta8,alpha8,phi8,id8,w8)
              robust_ipw_ex <- c(robust_ipw_ex,sqrt(diag(correction8$robust))[2])
              MD_ipw_ex <- c(MD_ipw_ex,sqrt(diag(correction8$varMD))[2])
              KC_ipw_ex <- c(KC_ipw_ex,sqrt(diag(correction8$varKC))[2])
              FG_ipw_ex <- c(FG_ipw_ex,sqrt(diag(correction8$varFG))[2])
            }else{
              est_ipw_ex <- c(est_ipw_ex,NA)
              std_ipw_ex <- c(std_ipw_ex,NA)
              robust_ipw_ex  <- c(robust_ipw_ex,NA)
              MD_ipw_ex <- c(MD_ipw_ex,NA)
              KC_ipw_ex <- c(KC_ipw_ex,NA)
              FG_ipw_ex <- c(FG_ipw_ex,NA)
            }
          }else{
            est_ipw_ex <- c(est_ipw_ex,NA)
            std_ipw_ex <- c(std_ipw_ex,NA)
            robust_ipw_ex  <- c(robust_ipw_ex,NA)
            MD_ipw_ex <- c(MD_ipw_ex,NA)
            KC_ipw_ex <- c(KC_ipw_ex,NA)
            FG_ipw_ex <- c(FG_ipw_ex,NA)
          } 
          
          robust_ipw_ex;MD_ipw_ex;KC_ipw_ex;FG_ipw_ex
          est_ipw_ex;std_ipw_ex
          
          ### IPW_CLU independent
          if(is.null(ipw_clu_ind$value)==0){
            if(ipw_clu_ind$value$converged==1){
              phi9 <- ipw_clu_ind$value$phi
              alpha9 <- ipw_clu_ind$value$alpha
              beta9 <- ipw_clu_ind$value$beta
              est_ipw_clu_ind <- c(est_ipw_clu_ind,beta9[2])
              std_ipw_clu_ind <- c(std_ipw_clu_ind,summary(ipw_clu_ind$value)$se.robust[2])
              w9 <- d2$weight2
              y9 <- d2$y
              X9 <- cbind(rep(1,dim(d2)[1]),d2$arm)
              id9 <- d2$cluster
              correction9 <- BCVexch(y9,X9,beta9,alpha9,phi9,id9,w9)
              robust_ipw_clu_ind <- c(robust_ipw_clu_ind,sqrt(diag(correction9$robust))[2])
              MD_ipw_clu_ind <- c(MD_ipw_clu_ind,sqrt(diag(correction9$varMD))[2])
              KC_ipw_clu_ind <- c(KC_ipw_clu_ind,sqrt(diag(correction9$varKC))[2])
              FG_ipw_clu_ind <- c(FG_ipw_clu_ind,sqrt(diag(correction9$varFG))[2])
            }else{
              est_ipw_clu_ind <- c(est_ipw_clu_ind,NA)
              std_ipw_clu_ind <- c(std_ipw_clu_ind,NA)
              robust_ipw_clu_ind <- c(robust_ipw_clu_ind,NA)
              MD_ipw_clu_ind <- c(MD_ipw_clu_ind,NA)
              KC_ipw_clu_ind <- c(KC_ipw_clu_ind,NA)
              FG_ipw_clu_ind <- c(FG_ipw_clu_ind,NA)
            }
          }else{
            est_ipw_clu_ind <- c(est_ipw_clu_ind,NA)
            std_ipw_clu_ind <- c(std_ipw_clu_ind,NA)
            robust_ipw_clu_ind <- c(robust_ipw_clu_ind,NA)
            MD_ipw_clu_ind <- c(MD_ipw_clu_ind,NA)
            KC_ipw_clu_ind <- c(KC_ipw_clu_ind,NA)
            FG_ipw_clu_ind <- c(FG_ipw_clu_ind,NA)
          }
          
          robust_ipw_clu_ind;MD_ipw_clu_ind;KC_ipw_clu_ind;FG_ipw_clu_ind
          est_ipw_clu_ind;std_ipw_clu_ind
          
          ## ipw_clu exchangeable
          
          if(is.null(ipw_clu_ex$value)==0){
            if( ipw_clu_ex$value$converged==1){
              phi10 <- ipw_clu_ex$value$phi
              alpha10 <- ipw_clu_ex$value$alpha
              beta10 <- ipw_clu_ex$value$beta
              est_ipw_clu_ex <- c(est_ipw_clu_ex,beta10[2])
              std_ipw_clu_ex <- c(std_ipw_clu_ex,summary(ipw_clu_ex$value)$se.robust[2])
              w10 <- d2$weight2
              y10 <- d2$y
              X10 <- cbind(rep(1,dim(d2)[1]),d2$arm)
              id10 <- d2$cluster
              correction10 <- BCVexch(y10,X10,beta10,alpha10,phi10,id10,w10)
              robust_ipw_clu_ex <- c(robust_ipw_clu_ex,sqrt(diag(correction10$robust))[2])
              MD_ipw_clu_ex <- c(MD_ipw_clu_ex,sqrt(diag(correction10$varMD))[2])
              KC_ipw_clu_ex <- c(KC_ipw_clu_ex,sqrt(diag(correction10$varKC))[2])
              FG_ipw_clu_ex <- c(FG_ipw_clu_ex,sqrt(diag(correction10$varFG))[2])
            }else{
              est_ipw_clu_ex <- c(est_ipw_clu_ex,NA)
              std_ipw_clu_ex <- c(std_ipw_clu_ex,NA)
              robust_ipw_clu_ex  <- c(robust_ipw_clu_ex,NA)
              MD_ipw_clu_ex <- c(MD_ipw_clu_ex,NA)
              KC_ipw_clu_ex <- c(KC_ipw_clu_ex,NA)
              FG_ipw_clu_ex <- c(FG_ipw_clu_ex,NA)
            }
          }else{
            est_ipw_clu_ex <- c(est_ipw_clu_ex,NA)
            std_ipw_clu_ex <- c(std_ipw_clu_ex,NA)
            robust_ipw_clu_ex  <- c(robust_ipw_clu_ex,NA)
            MD_ipw_clu_ex <- c(MD_ipw_clu_ex,NA)
            KC_ipw_clu_ex <- c(KC_ipw_clu_ex,NA)
            FG_ipw_clu_ex <- c(FG_ipw_clu_ex,NA)
          } 
          
          robust_ipw_clu_ex;MD_ipw_clu_ex;KC_ipw_clu_ex;FG_ipw_clu_ex
          est_ipw_clu_ex;std_ipw_clu_ex
        }
        
        # save the results as data frame
        est_ind <- data.frame(est_true_ind=est_true_ind,est_ucra_ind=est_ucra_ind,
                              est_cra_ind=est_cra_ind,est_ipw_ind=est_ipw_ind,
                              est_ipw_clu_ind=est_ipw_clu_ind)
        est_ex <- data.frame(est_true_ex=est_true_ex,est_ucra_ex=est_ucra_ex,
                             est_cra_ex=est_cra_ex,est_ipw_ex=est_ipw_ex,
                             est_ipw_clu_ex=est_ipw_clu_ex)
        
        std_ind <- data.frame(std_true_ind=std_true_ind,std_ucra_ind=std_ucra_ind,
                              std_cra_ind=std_cra_ind,std_ipw_ind=std_ipw_ind,
                              std_ipw_clu_ind=std_ipw_clu_ind)
        std_ex <- data.frame(std_true_ex=std_true_ex,std_ucra_ex=std_ucra_ex,
                             std_cra_ex=std_cra_ex,std_ipw_ex=std_ipw_ex,
                             std_ipw_clu_ex=std_ipw_clu_ex)
        
        robust_ind <- data.frame(robust_true_ind=robust_true_ind,robust_ucra_ind=robust_ucra_ind,
                                 robust_cra_ind=robust_cra_ind,robust_ipw_ind=robust_ipw_ind,
                                 robust_ipw_clu_ind=robust_ipw_clu_ind)
        robust_ex <- data.frame(robust_true_ex=robust_true_ex,robust_ucra_ex=robust_ucra_ex,
                                robust_cra_ex=robust_cra_ex,robust_ipw_ex=robust_ipw_ex,
                                robust_ipw_clu_ex=robust_ipw_clu_ex)
        
        MD_ind <- data.frame(MD_true_ind=MD_true_ind,MD_ucra_ind=MD_ucra_ind,
                             MD_cra_ind=MD_cra_ind,MD_ipw_ind=MD_ipw_ind,
                             MD_ipw_clu_ind=MD_ipw_clu_ind)
        MD_ex <- data.frame(MD_true_ex=MD_true_ex,MD_ucra_ex=MD_ucra_ex,
                            MD_cra_ex=MD_cra_ex,MD_ipw_ex=MD_ipw_ex,
                            MD_ipw_clu_ex=MD_ipw_clu_ex)
        
        
        KC_ind <- data.frame(KC_true_ind=KC_true_ind,KC_ucra_ind=KC_ucra_ind,
                             KC_cra_ind=KC_cra_ind,KC_ipw_ind=KC_ipw_ind,
                             KC_ipw_clu_ind=KC_ipw_clu_ind)
        KC_ex <- data.frame(KC_true_ex=KC_true_ex,KC_ucra_ex=KC_ucra_ex,
                            KC_cra_ex=KC_cra_ex,KC_ipw_ex=KC_ipw_ex,
                            KC_ipw_clu_ex=KC_ipw_clu_ex)
        
        FG_ind <- data.frame(FG_true_ind=FG_true_ind,FG_ucra_ind=FG_ucra_ind,
                             FG_cra_ind=FG_cra_ind,FG_ipw_ind=FG_ipw_ind,
                             FG_ipw_clu_ind=FG_ipw_clu_ind)
        FG_ex <- data.frame(FG_true_ex=FG_true_ex,FG_ucra_ex=FG_ucra_ex,
                            FG_cra_ex=FG_cra_ex,FG_ipw_ex=FG_ipw_ex,
                            FG_ipw_clu_ex=FG_ipw_clu_ex)
        
        result <- list(est_ind=est_ind,est_ex=est_ex,std_ind=std_ind,std_ex=std_ex,
                       robust_ind=robust_ind,robust_ex=robust_ex,
                       MD_ind=MD_ind,KC_ind=KC_ind,FG_ind=FG_ind,
                       MD_ex=MD_ex,KC_ex=KC_ex,FG_ex=FG_ex,
                       mistrt = MISTRT,
                       miscon = MISCON,
                       no_interaction = no_interaction)
        names <- paste(res_name,k,icc,iccm,'.RData',sep='')
        save(result,file=names)
      }
    }
  }
}
