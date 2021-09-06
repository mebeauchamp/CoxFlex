###################################################################################################
##                                                                                               ##
## R program for flexible Cox model including non-linear (NL) and/or time-dependent (TD) effects ##
##                                                                                               ##
###################################################################################################

# Original code developed by Willy Wynant. Subsequent corrections and modifications made by Yishu Wang.
# Latest modifications and corrections by Marie-Eve Beauchamp (since September 2018).

# Date of last major update: February 26, 2018 
# Last update: March 25, 2020

# For questions or comments about the code contact Marie-Eve Beauchamp (marie-eve.beauchamp at rimuhc.ca).

# Reference:
# Wynant W, Abrahamowicz M. Impact of the model building strategy on the inference 
# about time-dependent and non-linear covariate effects in survival analysis. 
# Stat Med 2014;33(19):3318-3337.


### Main functions:

  ## CoxFlex(data, Type, variables, TD, NL, m, p, knots) ##
    # The main function, which estimates the flexible Cox model (by calling the function last_prog) and 
    # evaluate the significance of each TD and LN effect requested.  

  ## last_prog(data, Type, variables, TD, NL, m, p, knots=-999) ##
    # Function estimating the flexible Cox model with TD and/or NL effects.

    # Arguments:
    #   data: data frame containing the data. The first column must be the variable identifying the subjects.
    #   Type: vector of variables in data indicating the start and stop of time intervals, for time-varying 
    #         data (counting process data), and the event (1=event, 0=censored). E.g., Type=c("Start","Stop","Event"). 
    #         If time-invariant data: Type=c("Time","Event").     
    #   variables: a vector containing names of independent variables in the model.
    #   TD: a vector indicating for each independent variable if the TD effect is modeled (1=yes, 0=no).
    #   NL: a vector indicating for each independent variable if the NL effect is modeled (1=yes, 0=no). Can be 1  
    #       only for continuous variables.
    #   m: number of interior knots (the same for all TD and NL effects).
    #   p: degree of the splines (the same for all TD and NL effects).
    #   knots: position of interior knots.


  ## backward_selection2(data, Type, variables, continuous, TD, NL, m=1, p=2, alpha_back=0.05, knots=-999) ## 
    # Backward selection of TD and NL effects and variables.

    # Arguments:
    #   data: data frame containing the data. The first column must be the variable identifying the subjects.
    #   Type: vector of variables in data indicating the start and stop of time intervals, for time-varying 
    #         data (counting process data), and the event (1=event, 0=censored). E.g., Type=c("Start","Stop","Event"). 
    #         If time-invariant data: Type=c("Time","Event").      
    #   variables: a vector containing names of independent variables in the model.
    #   continuous: a vector indicating whether each variable is continuous (1=yes, 0=no).
    #   TD: a vector indicating if TD effect of corresponding variable is forced (1=yes, 0=no).
    #   NL: a vector indicating if NL effect of corresponding variable is forced  (1=yes, 0=no). Can be 1 only for 
    #       continuous variables.           
    #   m: number of interior knots (the same for all TD and NL effects). Default m=1.
    #   p: degree of the splines (the same for all TD and NL effects). Default p=2.
    #   alpha_back: the alpha value used to select effects and variables. Default alpha_back=0.05.
    #   knots: position of interior knots. -999 indicates that knots are automatically allocated.


### Plotting functions:

  ## plot.FlexSurv(model.FlexSurv, variable, TD, NL, TimePoint=-999, ref.value.NL=0,...) ## 
    # Plot (or add a line for) the estimated TD or NL effect, or their product, for one variable.

  ## lines.FlexSurv(model.FlexSurv, variable, TD, NL, TimePoint=-999, ref.value.NL=0,...) ## 
    # Add a line for the estimated TD or NL effect, or their product, for one variable to an existing plot.

    # Arguments:
    #   model.FlexSurv: object of result from CoxFlex or last_prog functions.
    #   variable: name of variable for which the effect(s) is plotted. 
    #   TD: binary indicator of whether the estimated TD effect is plotted. 
    #   NL: binary indicator of whether the estimated NL effect is plotted. 
    #   TimePoint: when plotting TD and NL effects of a variable requested in on function call, indicate the 
    #              time point at which the NL effect is plotted. Default TimePoint=-999 means both effects 
    #              are plotted sseperately.
    #   ref.value.NL: log(HR) for NL effect plotted with respect to this reference value. Default ref.value.NL=0.
    #   ...: other options possible, e.g. type, col, xlab, ylab, main, etc.


### Internal functions:

  ## DvlpMatrix(data, listeT, ncol, TypeNEW) ##
    # Rearrange the data to recreate the risk sets. 
    # It deletes rows not in any risk sets to improve efficiency. 
    # And creates row for each event times so that splines TD effect are correctly represented in each risk set.

    # Arguments:
    #  data: the data
    #  listeT: list of unique event times in data
    #  ncol: number of columns in data
    #  TypeNEW: position of variables in argument Type in data


  ## spli(x, j, p, knots) ##
    # Calulate the spline value for x, for one spline in the spline basis

    # Arguments:
    #   x: value
    #   j: number of the spline in the basis
    #   p: degree of splines
    #   knots: position of interior knots


library(survival)
library(splines)


CoxFlex <- function(data, Type, variables, TD, NL, m, p, knots){
  
  # Estimating the flexible Cox model
  res_principal <- last_prog(data, Type, variables, TD, NL, m, p, knots)
  
  TDtestP <- rep(-999, length(variables))
  TDtestN <- rep(-999, length(variables))
  TDtest <- rep(-999, length(variables))
  NLtestP <- rep(-999, length(variables))
  NLtestN <- rep(-999, length(variables))
  NLtest <- rep(-999, length(variables))
  
  # Evaluating the significance of each TD and NL effect requested
  for (kt in 1:length(variables)){
    if (TD[kt]==1){
      TDnew <- TD
      TDnew[kt] <- 0
      resSEC <- last_prog(data, Type, variables, TDnew, NL, m, p, knots)
      TDtestP[kt] <- resSEC$Partial_Log_Likelihood
      TDtestN[kt] <- resSEC$Number_of_parameters
      TDtest[kt] <- 1-pchisq(-2*(TDtestP[kt]-res_principal$Partial_Log_Likelihood), res_principal$Number_of_parameters-TDtestN[kt])
    }
    if (NL[kt]==1){
      NLnew <- NL
      NLnew[kt] <- 0
      resSEC <- last_prog(data, Type, variables, TD, NLnew, m, p, knots)
      NLtestP[kt] <- resSEC$Partial_Log_Likelihood
      NLtestN[kt] <- resSEC$Number_of_parameters
      NLtest[kt] <- 1-pchisq(-2*(NLtestP[kt]-res_principal$Partial_Log_Likelihood), res_principal$Number_of_parameters-NLtestN[kt])
    }
  }
  
  # Preparing output of results (MEB: EXPLORE Rescox to improve output)
  Rescox <- matrix(ncol=6, nrow=sum(2*(NL+TD==2)+1*(NL+TD!=2)))
  colnames(Rescox) <- c("", "coef", "exp(coef)", "se(coef)", "z", "p")
  rownames(Rescox) <- rep(c(""), sum(2*(NL+TD==2)+1*(NL+TD!=2)))
  
  indexrescox <- 0
  for (yu in 1:length(variables)){
    indexrescox <- indexrescox+1
    if (TD[yu]==0 & NL[yu]==0){
      Rescox[indexrescox,1] <- variables[yu]
      Rescox[indexrescox,2] <- round(as.numeric(res_principal$coefficients[yu]), 3)
      Rescox[indexrescox,3] <- round(as.numeric(exp(res_principal$coefficients[yu])), 3)
      Rescox[indexrescox,4] <- round(as.numeric(res_principal$Standard_Error[yu]), 3)
      Rescox[indexrescox,5] <- round(as.numeric(res_principal$coefficients[yu]/res_principal$Standard_Error[yu]), 3)
      Rescox[indexrescox,6] <- round(as.numeric(1-pchisq((res_principal$coefficients[yu]/res_principal$Standard_Error[yu])^2,df=1)), 3)
    }
    if (TD[yu]==0 & NL[yu]==1){
      Rescox[indexrescox,1] <- paste("NL(",variables[yu],")", sep="")
      Rescox[indexrescox,2] <- "---"
      Rescox[indexrescox,3] <- "splines"
      Rescox[indexrescox,4] <- "---"
      Rescox[indexrescox,5] <- "---"
      Rescox[indexrescox,6] <- round(NLtest[yu], 3)
    }
    if (TD[yu]==1 & NL[yu]==0){
      Rescox[indexrescox,1] <- paste("TD(",variables[yu],")", sep="")
      Rescox[indexrescox,2] <- "---"
      Rescox[indexrescox,3] <- "splines"
      Rescox[indexrescox,4] <- "---"
      Rescox[indexrescox,5] <- "---"
      Rescox[indexrescox,6] <- round(TDtest[yu], 3)
    }
    if (TD[yu]==1 & NL[yu]==1){
      Rescox[indexrescox,1] <- paste("NL(",variables[yu],")", sep="")
      Rescox[indexrescox,2] <- "---"
      Rescox[indexrescox,3] <- "splines"
      Rescox[indexrescox,4] <- "---"
      Rescox[indexrescox,5] <- "---"
      Rescox[indexrescox,6] <- round(NLtest[yu], 3)
      indexrescox<-indexrescox+1
      Rescox[indexrescox,1] <- paste("TD(",variables[yu],")", sep="")
      Rescox[indexrescox,2] <- "---"
      Rescox[indexrescox,3] <- "splines"
      Rescox[indexrescox,4] <- "---"
      Rescox[indexrescox,5] <- "---"
      Rescox[indexrescox,6] <- round(TDtest[yu], 3)
    }
  }
  
  list(Partial_Log_Likelihood = res_principal$Partial_Log_Likelihood, 
       Number_of_parameters = res_principal$Number_of_parameters, 
       Number_events=res_principal$Number_events, 
       Number_knots = res_principal$Number_knots, 
       Degree_of_splines = res_principal$Degree_of_splines, 
       knots_covariates = res_principal$knots_covariates, 
       knots_time = res_principal$knots_time, 
       coefficients = res_principal$coefficients, 
       Standard_Error=res_principal$Standard_Error,
       coefficients_splines_NL = res_principal$coefficients_splines_NL,
       coefficients_splines_TD = res_principal$coefficients_splines_TD,
       variables=res_principal$variables,
       # Modif by MEB: correction of the wrong variance in the next line
       #coef=suppressWarnings(as.numeric(Rescox[,2])), var=suppressWarnings((as.numeric(Rescox[,2]))^2), pvalue=suppressWarnings(as.numeric(Rescox[,6])))
       coef=suppressWarnings(as.numeric(Rescox[,2])), var=suppressWarnings((as.numeric(Rescox[,4]))^2), pvalue=suppressWarnings(as.numeric(Rescox[,6])))
}



last_prog <- function(data, Type, variables, TD, NL, m, p, knots=-999){
 
  if (length(Type)==2){
    Type2 <- Type
    data$StartV0 <- rep(0, dim(data)[1])
    Type <- c("StartV0", Type2[1], Type2[2])
  } 
 
  i1 <- sum((NL + TD) == 0)          # number of non-nl non-td variables   
  i2 <- sum(((NL == 1) & (TD == 0))) # number of NL non-td variables
  i3 <- sum(((NL == 0) & (TD == 1))) # number of non-nl TD variables 
  i4 <- sum((NL + TD) == 2)          # number of NL and TD variables 
  #nonpara <- TD + NL 

  V <- length(variables)             # number ofcovariates
  
  variablesNEW <- match(variables, names(data)) # positions in data of items in argument variables 
  TypeNEW <- match(Type, names(data))           # positions in data of items in argument Type

  ## Store interior and exterior knots in (2*(p+1)+m) columns, for variables (rows 1:V) and time (row V+1)
  
  listeprobaquantile <- seq(1, m)/(m + 1) 
  knotsNEW <- matrix(nrow = V + 1, ncol = p + 1 + m + p + 1) 

  # If no user-defined knots
  if (is.matrix(knots)==FALSE){  
    if (is.numeric(knots)==TRUE & knots==-999) { # use default knots 
      
      # Knots for variables:
        # p+1 exterior knots at min(variables[i])
        # interior knots at quantiles
        # p+1 exterior knots equally spaced between max(variable) and max(variable)+p
      for (i in 1:V) {
        knotsNEW[i, ] <- c(rep(min(data[, variables[i]]), p+1), 
                           quantile(data[, variables[i]], probs=listeprobaquantile),  
                           seq(max(data[, variables[i]]), max(data[, variables[i]])+p, 1)) 
      }
      
      # Knots for time:
        # p+1 exterior knots at 0. 
	  # MEB's note: This assumes f-up always start at 0. Need to make it more general by using min(start). Need to change listeT below to start at min(start) instead of 0.
        # interior knots at quantiles of event times
        # p+1 exterior knots equally spaced between max(time) and max(time)+p
      knotsNEW[V+1, ] <- c(rep(0, p + 1), 
                          quantile(data[data[, TypeNEW[3]]==1, TypeNEW[2]], probs=listeprobaquantile), 
                          seq(max(data[data[, TypeNEW[3]]==1, TypeNEW[2]]), max(data[data[, TypeNEW[3]]==1, TypeNEW[2]])+p, 1)) 
    }                                                                          
  }

  # If (some) user-defined interior knots are specified   
  if (is.matrix(knots) == TRUE){
    
    if (dim(knots)[1] != (length(variables)+1) | dim(knots)[2] != m) 
      stop("Error Message: variable knots should be a matrix of dimension (length(variables)+1)*m")

    # Knots for variables     
    for (i in 1:V) {
      if (is.na(knots[i,1])==TRUE){ # if user-defined knots for that variable is NA, use default 
        knotsNEW[i, ] <- c(rep(min(data[, variables[i]]), p + 1),   
                          quantile(data[, variables[i]], probs=listeprobaquantile), 
                          seq(max(data[, variables[i]]), max(data[, variables[i]])+p, 1))

      } else { # use user-defined interior knots 
        knotsNEW[i, ] <- c(rep(min(data[, variables[i]]), p + 1),   
                           knots[i,], 
                           seq(max(data[, variables[i]]), max(data[, variables[i]])+p, 1))
      }
    }

    # Knots for time      
    if (is.na(knots[(V+1),1]) == TRUE){ # if user-defined knots for time is NA, use default  
      knotsNEW[(V+1), ] <- c(rep(0, p + 1), 
                            quantile(data[data[, TypeNEW[3]]==1, TypeNEW[2]], probs=listeprobaquantile), 
                            seq(max(data[data[, TypeNEW[3]]==1, TypeNEW[2]]), max(data[data[, TypeNEW[3]]==1, TypeNEW[2]])+p, 1))

    } else { # use user-defined interior knots
      knotsNEW[(V+1), ] <- c(rep(0, p + 1), 
                             knots[(V+1),], 
                             seq(max(data[data[, TypeNEW[3]]==1, TypeNEW[2]]), max(data[data[, TypeNEW[3]]==1, TypeNEW[2]])+p, 1))
    }
  }
 
  ## Prepare data (QWR) for estimation 
  
  data <- as.matrix(data)
  listeT <- c(0, sort(unique(data[data[, TypeNEW[3]]==1, TypeNEW[2]]))) # list of unique event times 
  ncol <- dim(data)[2] 
  
  X <- split(data, data[, 1]) 
  matX <- sapply(X, DvlpMatrix, listeT=listeT, ncol=ncol, TypeNEW=TypeNEW) 
  QWR <- do.call(rbind, matX)  # Creates data for risk sets

  # Add splines for NL effects in QWR    
  nbNL <- sum(NL)
  if (nbNL != 0){  
    for (i in 1:nbNL){ 
      QWR <- cbind(QWR, splineDesign(knotsNEW[seq(1,V,1)[NL==1][i],], x=QWR[, variablesNEW[NL==1][i]], ord=p+1)[,-1])
    }
  }

  # Add splines for TD effects in QWR (splines for time)   
  QWR <- cbind(QWR, splineDesign(knotsNEW[V+1,], x=QWR[, TypeNEW[2]], ord=p+1))

  ## Construct equation (tt) of Cox model (part to be estimated at all steps of all iterations)  
  tt <- paste("modX<-coxph(Surv(QWR[,", TypeNEW[1], "],QWR[,", 
              TypeNEW[2], "],QWR[,", TypeNEW[3], "])~", sep="")

  Nn <- 0 # number of parameters added to model up to now
  
  # Add parameters for variables without NL and TD effects 
  nbpara <- sum((NL == 0 & TD == 0)) 
  if (nbpara != 0){ 
    for (k in 1:nbpara){
      tt <- paste(tt, "QWR[,", variablesNEW[NL==0 & TD==0][k], "]+", sep="")
      Nn <- Nn + 1
    }
  }
  
  # Add parameters for variables with NL effect only 
  nbonlyNL <- sum((NL == 1 & TD == 0)) 
  if (nbonlyNL != 0){
    onlyNL <- match(variablesNEW[NL == 1 & TD == 0], variablesNEW[NL==1]) 
    for (k in 1:nbonlyNL){
      covp <- paste("QWR[,", (dim(data)[2] + (onlyNL[k]-1)*(m+p) + 1):(dim(data)[2] + onlyNL[k]*(m+p)), "]", sep="")
      tt <- paste(tt, paste(c(covp,""), collapse="+"))
      Nn <- Nn + m+p
    }
  }
  
  # Add parameters for variables with TD effect only
  nbonlyTD <- sum((NL == 0 & TD == 1)) 
  if (nbonlyTD != 0){
    for (k in 1:nbonlyTD){
      flag <- dim(QWR)[2] + 1
      QWR <- cbind(QWR, QWR[, variablesNEW[NL==0 & TD==1][k]] * 
                     QWR[, (dim(data)[2] + (m+p)*nbNL + 1):(dim(data)[2] + (m+p)*(nbNL+1) + 1)]) 
      covp <- paste("QWR[,", flag:dim(QWR)[2], "]", sep="")
      tt <- paste(tt, paste(c(covp,""), collapse="+"))  
      Nn <- Nn + m+p+1
    }
  }
  
  tt2 <- tt

  nbNLTD <- sum((NL == 1 & TD == 1)) 
  NOTonlyNL <- match(variablesNEW[NL == 1 & TD == 1], variablesNEW[NL==1])
  
  ### If some variables have up to 2 effects to be estimated ###
  
  if (nbNLTD != 0){ 
 
    ## Iteration 1, step 1 (estimate NL effects) ##
    
    for (k in 1:nbNLTD){
      covp <- paste("QWR[,", (dim(data)[2] + (NOTonlyNL[k]-1)*(m+p)+1):(dim(data)[2] + NOTonlyNL[k]*(m+p)), "]", sep="")
      tt2 <- paste(tt2, paste(c(covp,""), collapse="+")) 
    }
    
    vrais <- c()
    modX <- coxph(eval(parse(text=substr(tt2, 13, nchar(tt2) - 1))), method="efron")
    vrais <- c(vrais, modX$loglik[2])
    
    ## Iteration 1, step 2 (estimate TD effects) ##
    
    mod <- modX
    tt2 <- tt
    VV <- matrix(ncol=nbNLTD*(m+p+1), nrow=dim(QWR)[1])
    
    for (k in 1:nbNLTD){
      VV[,((1+m+p)*(k-1)+1):((1+m+p)*k)] <- QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))] * 
                                            as.vector(QWR[,(dim(data)[2]+(NOTonlyNL[k]-1)*(m+p)+1):(dim(data)[2]+NOTonlyNL[k]*(m+p))] %*% 
                                                      mod$coef[(Nn+(k-1)*(m+p)+1):(Nn+k*(m+p))])
      covp <- paste( "VV[,", ((m+p+1)*(k-1)+1):((k)*(m+p+1)), "]", sep="")
      tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
    }
    
    modX <- coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))), method="efron")
    vrais <- c(vrais, modX$loglik[2])

    diff <- 1
    
    ## Subsequent iterations until convergence ##
    
    while(diff > 0.00001){ 

      ## Step 1 (estimate NL effects) ##
      
      mod <- modX
      tt2 <- tt
      VV <- matrix(ncol=nbNLTD*(m+p),nrow=dim(QWR)[1])
      
      for (k in 1:nbNLTD){
        VV[,((m+p)*(k-1)+1):((m+p)*(k))] <- QWR[,(dim(data)[2]+(NOTonlyNL[k]-1)*(m+p)+1):(dim(data)[2]+NOTonlyNL[k]*(m+p))] * 
                                            as.vector(QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))] %*% 
                                                      mod$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+(k)*(m+p+1))])
        covp <- paste( "VV[,", ((m+p)*(k-1)+1):((k)*(m+p)), "]", sep="")
        tt2 <- paste(tt2, paste(c(covp,""), collapse="+"))
      }
      
      modX <- coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))), method="efron")
      vrais <- c(vrais, modX$loglik[2])
      
      ## Step 2 (estimate TD effects) ##
      
      mod <- modX
      tt2 <- tt
      VV <- matrix(ncol=nbNLTD*(m+p+1), nrow=dim(QWR)[1])

      for (k in 1:nbNLTD){
        VV[,((1+m+p)*(k-1)+1):((1+m+p)*(k))] <- QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))] * 
                                                as.vector(QWR[,(dim(data)[2]+(NOTonlyNL[k]-1)*(m+p)+1):(dim(data)[2]+NOTonlyNL[k]*(m+p))] %*% 
                                                          mod$coef[(Nn+(k-1)*(m+p)+1):(Nn+(k)*(m+p))])
        covp <- paste( "VV[,", ((m+p+1)*(k-1)+1):((k)*(m+p+1)), "]", sep = "")
        tt2 <- paste(tt2, paste(c(covp,""), collapse= "+"))
      }
      
      modX <- coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))), method="efron")
      vrais <- c(vrais, modX$loglik[2])
      diff <- abs(vrais[length(vrais)] - vrais[length(vrais)-2])
    }
    
  } else {
    
    ### ELSE IF (All variables have 0 or 1 effect)
    
    modX <- coxph(eval(parse(text = substr(tt2, 13, nchar(tt2) - 1))), method="efron")
    vrais <- c(modX$loglik[2])
  }
  rm(QWR, X, matX)
  gc() # Added by MEB
  
  ### Store parameters estimated ###
  
  MAT <- matrix(ncol=V, nrow=1+m+p+m+p+1) 
  if (i1 != 0){ # for non-nl non-td variables
    for (j in 1:i1){
      MAT[1, j] <- modX$coef[j]
    }
  }
  if (i2 != 0){ # for NL non-td variables
    for (j in 1:i2){
      MAT[2:(m+p+1), i1+j] <- modX$coef[(i1 + (j-1)*(m+p) + 1):(i1 + (j-1)*(m+p) + m+p)]
    }
  }
  if (i3 != 0){ # for non-nl TD variables
    for (j in 1:i3){
      MAT[(m+p+2):(2*m+2*p+2), i1+i2+j] <- modX$coef[(i1 + i2*(m+p) + 
                                                    (j-1)*(m+p+1) + 1):(i1 + i2*(m+p) + (j-1)*(m+p+1) + m+p+1)]
    }
  }
  if (i4 != 0){ # for NL and TD variables
    for (j in 1:i4){
      MAT[2:(m+p+1), i1+i2+i3+j] <- mod$coef[(i1 + i2*(m+p) + i3*(m+p+1) + 
                                            (j-1)*(m+p) + 1):(i1 + i2*(m+p) + i3*(m+p+1) + (j-1)*(m+p) + m+p)]
      MAT[(m+p+2):(2*(m+p+1)), i1+i2+i3+j] <- modX$coef[(i1 + i2*(m+p) + i3*(m+p + 1) + 
                                                        (j-1)*(m+p+1) + 1):(i1 + i2*(m+p) + i3*(m+p+1) + (j-1)*(m+p+1) + m+p+1)]
    }
  }
  
  ### Store SE estimated ###
  
  MATse <- matrix(ncol=V, nrow=1+m+p+m+p+1) 
  if (i1 != 0) {
    for (j in 1:i1) {
      MATse[1, j] <- sqrt(diag(modX$var)[j])
    }
  }
  if (i2 != 0){
    for (j in 1:i2){
      MATse[2:(m+p+1), i1+j] <- sqrt(diag(modX$var)[(i1 + (j-1)*(m+p) + 1):(i1 + (j-1)*(m+p) + m+p)])
    }
  }
  if (i3 != 0){
    for (j in 1:i3){
      MATse[(m+p+2):(2*m+2*p+2), i1+i2+j] <- sqrt(diag(modX$var)[(i1 + i2*(m+p) + 
                                                                (j-1)*(m+p+1) + 1):(i1 + i2*(m+p) + (j-1)*(m+p+1) + m+p+1)])
    }
  }
  if (i4 != 0){
    for (j in 1:i4){
      MATse[2:(m+p+1), i1+i2+i3+j] <- sqrt(diag(mod$var)[(i1 + i2*(m+p) + i3*(m+p+1) + 
                                                        (j-1)*(m+p) + 1):(i1 + i2*(m+p) + i3*(m+p+1) + (j-1)*(m+p) + m+p)])
      MATse[(m+p+2):(2*(m+p+1)), i1+i2+i3+j] <- sqrt(diag(modX$var)[(i1 + i2*(m+p) + i3*(m+p+1) + (j-1)*(m+p+1) + 1):(i1 + i2*(m+p) + 
                                                                    i3*(m+p+1) + (j-1)*(m+p+1) + m + p + 1)])
    }
  }
  
  ### Prepare output ###
  
  var_order <- c(variablesNEW[NL==0 & TD==0], variablesNEW[NL==1 & TD==0], 
                 variablesNEW[NL==0 & TD==1], variablesNEW[NL==1 & TD==1])
  
  coefficients <- MAT[1,match(variablesNEW, var_order)] 
  names(coefficients) <- variables
  
  se_coef <- MATse[1, match(variablesNEW, var_order)] 
  
  coefficients_splines_NL <- as.matrix(MAT[2:(m+p+1), match(variablesNEW, var_order)])
  coefficients_splines_NL <- rbind(rep(0,V), coefficients_splines_NL)
  coefficients_splines_NL[1, (NL==0)] <- NA
  colnames(coefficients_splines_NL) <- variables
  
  coefficients_splines_TD <- as.matrix(MAT[(m+p+2):(2*(m+p+1)), match(variablesNEW, var_order)])
  colnames(coefficients_splines_TD) <- variables
  
  knots_covariates<-knotsNEW[1:V,]
  if (V>1) {rownames(knots_covariates) <- variables}
  if (V>1) {
    knots_covariates[(NL==0),] <- rep(NA, p+1+m+p+1)
  } else { 
    if (NL==0) {
      knots_covariates<-rep(NA, p+1+m+p+1)
    } 
  }
  knots_time <- knotsNEW[V+1,]
  
  nEvents <- sum(data[,TypeNEW[3]]==1) 
  
  rm(data, modX)
  gc()
  
  list(Partial_Log_Likelihood = vrais[length(vrais)], 
       Number_of_parameters = nbpara + nbonlyNL*(m+p) + nbonlyTD*(1+m+p) + nbNLTD*(m+p+m+p+1), 
       Number_events = nEvents, 
       Number_knots = m, 
       Degree_of_splines = p, 
       knots_covariates = knots_covariates, 
       knots_time = knots_time, 
       coefficients = coefficients, 
       Standard_Error = se_coef,
       coefficients_splines_NL = coefficients_splines_NL,
       coefficients_splines_TD = coefficients_splines_TD,
       variables = variables)
}



backward_selection2 <- function(data, Type, variables, continuous, TD, NL, m=1, p=2, alpha_back=0.05, knots=-999){

  V <- length(variables) # number of variables
  
  covB<-variables[order(1-continuous)] 
  Ftd<-TD[order(1-continuous)]
  Fnl<-NL[order(1-continuous)]
  
  TDf<-rep(1,V) 
  NLf<-c(rep(1,sum(continuous)),rep(0,sum(1-continuous)))
  nNLf<-sum(continuous) 
  nNLfNOT=V-nNLf 

  m1nl<-matrix(nrow=nNLf,ncol=nNLf,1)
  diag(m1nl)<-0
  m2nl<-matrix(nrow=nNLf,ncol=nNLfNOT,0)
  m3nl<-cbind(m1nl,m2nl)
  m4nl<-matrix(nrow=V,ncol=V,c(rep(1,nNLf),rep(0,nNLfNOT)),byrow=T)
  mNL<-rbind(NLf,m3nl,m4nl)
  
  ## TD matrix
  m1td<-matrix(nrow=nNLf,ncol=V,1)
  m2td<-matrix(nrow=nNLfNOT,ncol=nNLf,1)
  m3td<-matrix(nrow=nNLfNOT,ncol=nNLfNOT,1)
  diag(m3td)<-0
  m4td<-cbind(m2td,m3td)
  m5td<-matrix(nrow=nNLf,ncol=nNLf,1)
  diag(m5td)<-0
  m6td<-matrix(nrow=nNLf,ncol=nNLfNOT,1)
  m7td<-cbind(m5td,m6td)
  mTD<-rbind(TDf,m1td,m4td,m7td)

  for (i in 1:V){
    if(Ftd[i]==1) mTD[,i]<-1
    else if(Ftd[i]==-1) mTD[,i]<-0
    
    if(Fnl[i]==1) mNL[,i]<-1
    else if(Fnl[i]==-1) mNL[,i]<-0
  }
  
  PLback<-rep(0,V+nNLf+1) 
  DFback<-rep(-999,V+nNLf+1) 
  
  for (i in 1:(V+nNLf+1)){
    res<-CoxFlex(data,Type,variables=covB,TD=mTD[i,],NL=mNL[i,],m,p,knots=-999)
    PLback[i]<-res$Partial_Log_Likelihood 
    DFback[i]<-res$Number_of_parameters
  }

  Mind<-matrix(nrow=4,ncol=V,c(rep(1,2*V),rep(1,nNLf),rep(0,nNLfNOT),rep(2,nNLf),rep(1,nNLfNOT)),byrow=TRUE)
  for (i in 1:V){
    if(Ftd[i]==-1) Mind[2,i]<-0 # 2nd index TD 
    if(Fnl[i]==-1) Mind[3,i]<-0 # 3rd index NL 
  }
  Mind[4,]<-Mind[2,]+Mind[3,]

  pB<-rep(0,(V+nNLf)) # stores the test results of likelihood ratio tests for each model 
  
  for (j in 1:(V+nNLf)){
    if(DFback[1]-DFback[j+1]!=0)
      pB[j]<-1-pchisq(-2*(PLback[j+1]-PLback[1]),DFback[1]-DFback[j+1])
  }

  MpB<-matrix(nrow=5,ncol=V,0)
  
  MpB[4,1:nNLf]<-pB[1:nNLf] 
  if(sum(1-continuous)!=0){
  MpB[2,(nNLf+1):V]<-pB[(nNLf+1):V] 
  }
  MpB[5,1:nNLf]<-pB[(V+1):(V+nNLf)] 

  a<-which.max(MpB) 
  nc<-floor(a/5)+((a-5*floor(a/5))!=0) 
  nr<-a-5*floor(a/5)+5*((a-5*floor(a/5))==0)

  indexrow<-1                                       

  MatResPvalue<-matrix(nrow=2*length(variables)+sum(continuous),ncol=1)   
  MatResVariables<-matrix(nrow=2*length(variables)+sum(continuous),ncol=1) 
  
  if(MpB[nr,nc]<alpha_back){ 
    mbase<-CoxFlex(data,Type,variables=covB,TD=mTD[1,],NL=mNL[1,],m=1,p=2,knots=-999)
  }
  
  while(MpB[nr,nc]>=alpha_back) {     
    
    MatResPvalue[indexrow,1]<-MpB[nr,nc]       
    MatResVariables[indexrow,1]<-paste(nr,nc,sep="") 
    indexrow<-indexrow+1                              

    if (nr==1) Mind[1,nc]<-0 
    if (nr==2 | nr==5) Mind[2,nc]<-0 
    if (nr==3 | nr==4) Mind[3,nc]<-0 
    Mind[4,]<-Mind[2,]+Mind[3,] 
    
    mbase<-CoxFlex(data,Type,variables=covB[Mind[1,]==1],TD=Mind[2,][Mind[1,]==1],NL=Mind[3,][Mind[1,]==1],m=1,p=2,knots=-999)
    
    MpB<-matrix(nrow=5,ncol=V,rep(0,V*5))

    for (k in 1:V){ 
      
      if (Mind[2,k]==1 & Mind[3,k]==0) { 
        MindNew<-Mind
        if(Ftd[k]!=1) {
          MindNew[2,k]<-0 
          mnew<-CoxFlex(data,Type,variables=covB[MindNew[1,]==1],TD=MindNew[2,][MindNew[1,]==1],NL=MindNew[3,][MindNew[1,]==1],m=1,p=2,knots=-999)
          MpB[2,k]<-1-pchisq(-2*(mnew$Partial_Log_Likelihood-mbase$Partial_Log_Likelihood),mbase$Number_of_parameters-mnew$Number_of_parameters)
         
        }
        else MpB[2,k]<-0
      }
      
      if (Mind[2,k]==1 & Mind[3,k]==1) {
        MindNew<-Mind
        if(Ftd[k]!=1){
          MindNew[2,k]<-0 
          mnew<-CoxFlex(data,Type,variables=covB[MindNew[1,]==1],TD=MindNew[2,][MindNew[1,]==1],NL=MindNew[3,][MindNew[1,]==1],m=1,p=2,knots=-999)
          MpB[5,k]<-1-pchisq(-2*(mnew$Partial_Log_Likelihood-mbase$Partial_Log_Likelihood),mbase$Number_of_parameters-mnew$Number_of_parameters)
          
        }
        else MpB[5,k]<-0
      }
      
      if (Mind[3,k]==1 & Mind[2,k]==0) { #if NL + non-TD
        MindNew<-Mind
        if(Fnl[k]!=1){
          MindNew[3,k]<-0 
          mnew<-CoxFlex(data,Type,variables=covB[MindNew[1,]==1],TD=MindNew[2,][MindNew[1,]==1],NL=MindNew[3,][MindNew[1,]==1],m=1,p=2,knots=-999)
          MpB[3,k]<-1-pchisq(-2*(mnew$Partial_Log_Likelihood-mbase$Partial_Log_Likelihood),mbase$Number_of_parameters-mnew$Number_of_parameters)
         
        }
        else MpB[3,k]<-0
      }
      
      if (Mind[3,k]==1 & Mind[2,k]==1) { 
        MindNew<-Mind
        if(Fnl[k]!=1){
          MindNew[3,k]<-0 
          mnew<-CoxFlex(data,Type,variables=covB[MindNew[1,]==1],TD=MindNew[2,][MindNew[1,]==1],NL=MindNew[3,][MindNew[1,]==1],m=1,p=2,knots=-999)
          MpB[4,k]<-1-pchisq(-2*(mnew$Partial_Log_Likelihood-mbase$Partial_Log_Likelihood),mbase$Number_of_parameters-mnew$Number_of_parameters)
          
        }
        else MpB[4,k]<-0
      }
      
      if (Mind[3,k]==0 & Mind[2,k]==0 & Mind[1,k]==1) { 
        MindNew<-Mind
        if(Fnl[k]!=1 & Fnl[k]!=-1 & Ftd[k]!=1 & Ftd[k]!=-1){
          MindNew[1,k]<-0 # remove this variable 
          mnew<-CoxFlex(data,Type,variables=covB[MindNew[1,]==1],TD=MindNew[2,][MindNew[1,]==1],NL=MindNew[3,][MindNew[1,]==1],m=1,p=2,knots=-999)
          MpB[1,k]<-1-pchisq(-2*(mnew$Partial_Log_Likelihood-mbase$Partial_Log_Likelihood),mbase$Number_of_parameters-mnew$Number_of_parameters)
         
        }
        else MpB[1,k]<-0
      }
    } 
    
    a<-which.max(MpB)
    nc<-floor(a/5)+((a-5*floor(a/5))!=0)
    nr<-a-5*floor(a/5)+5*((a-5*floor(a/5))==0)
  }
  
  list(final_model=mbase)
}



plot.FlexSurv <- function(model.FlexSurv, variable, TD, NL, TimePoint=-999, ref.value.NL=0,...){
  
  variableNEW <- match(variable, model.FlexSurv$variables)
  
  # NOTE: line added by MEB because was causing error messages when only one variable was modeled as NL 
  #       (model.FlexSurv$knots_covariates is treated as a matrix below at 2 places)
  if (is.vector(model.FlexSurv$knots_covariates)){
    model.FlexSurv$knots_covariates <- matrix(model.FlexSurv$knots_covariates, nrow=1)
  }
  
  tdestim1 <- function(x){
    tdfct <- 0
    for (k in 1:(model.FlexSurv$Number_knots+model.FlexSurv$Degree_of_splines+1)){
      tdfct <- tdfct + model.FlexSurv$coefficients_splines_TD[k, variableNEW] * 
               spli(x, k, model.FlexSurv$Degree_of_splines, model.FlexSurv$knots_time)
    }
    return(tdfct)
  }
  
  nlestim1 <- function(y){
    nlfct <- 0
    for (k in 1:(model.FlexSurv$Number_knots+model.FlexSurv$Degree_of_splines+1)){
      nlfct <- nlfct + model.FlexSurv$coefficients_splines_NL[k, variableNEW] * 
               spli(y, k, model.FlexSurv$Degree_of_splines, model.FlexSurv$knots_covariates[variableNEW,])
    }
    return(nlfct)
  }
  
  if(TD==1){
    # Modif by MEB: seq by 1 instead of 0.01 because this is precise enough to have a smooth line in most analyses,
    # e.g. 100 units of follow-up imply 100 points. Seq by 0.01 was time consuming and create very large graphs
    # with normal follow-up, e.g. 365 days implied 36,500 points and 10 years or daily units implied 365,000 points.
    #axist<-seq(0,model.FlexSurv$knots_time[model.FlexSurv$Degree_of_splines+2+model.FlexSurv$Number_knots],0.01)
    axist <- seq(0, model.FlexSurv$knots_time[model.FlexSurv$Degree_of_splines + 2 + model.FlexSurv$Number_knots])
  }
  
  if(NL==1){
    axisx <- seq(model.FlexSurv$knots_covariates[variableNEW, 1], 
                 model.FlexSurv$knots_covariates[variableNEW, model.FlexSurv$Degree_of_splines + 2 + model.FlexSurv$Number_knots], 
                 0.01)
  }
  
  if (TD==1 & NL==0){
    plot(axist, tdestim1(axist),...)
  }
  
  if (TD==0 & NL==1){
    plot(axisx, nlestim1(axisx) - nlestim1(ref.value.NL),...)
  }  
  
  if (TD==1 & NL==1 & TimePoint==-999){
    op <- par(mfrow=c(2,1))
    plot(axist, tdestim1(axist),...)
    plot(axisx, nlestim1(axisx) - nlestim1(ref.value.NL),...)
    par(op)
  }
  
  if (TD==1 & NL==1 & TimePoint!=-999){
    plot(axisx, (nlestim1(axisx) - nlestim1(ref.value.NL)) * tdestim1(TimePoint),...)
  }
}



lines.FlexSurv <- function(model.FlexSurv, variable, TD, NL, TimePoint=-999, ref.value.NL=0,...){
  
  variableNEW <- match(variable, model.FlexSurv$variables)

  # NOTE: line added by MEB because was causing error messages when only one variable was modeled as NL 
  #       (model.FlexSurv$knots_covariates treated as a matrix below at 2 places)
  if (is.vector(model.FlexSurv$knots_covariates)){
    model.FlexSurv$knots_covariates <- matrix(model.FlexSurv$knots_covariates, nrow=1)
  }
  
  tdestim1 <- function(x){
    tdfct <- 0
    for (k in 1:(model.FlexSurv$Number_knots+model.FlexSurv$Degree_of_splines+1)){
      tdfct <- tdfct + model.FlexSurv$coefficients_splines_TD[k, variableNEW] * 
                spli(x, k, model.FlexSurv$Degree_of_splines, model.FlexSurv$knots_time)
    }
    return(tdfct)
  }
  
  nlestim1 <- function(y){
    nlfct <- 0
    for (k in 1:(model.FlexSurv$Number_knots+model.FlexSurv$Degree_of_splines+1)){
      nlfct <- nlfct + model.FlexSurv$coefficients_splines_NL[k, variableNEW] * 
                spli(y, k, model.FlexSurv$Degree_of_splines, model.FlexSurv$knots_covariates[variableNEW,])
    }
    return(nlfct)
  }
  
  if(TD==1){
    # Modif by MEB: seq by 1 instead of 0.01 because this is precise enough to have a smooth line in most analyses,
      # e.g. 100 units of follow-up imply 100 points. Seq by 0.01 was time consuming and create very large graphs
      # with normal follow-up, e.g. 365 days implied 36,500 points and 10 years or daily units implied 365,000 points.
  #axist <- seq(0, model.FlexSurv$knots_time[model.FlexSurv$Degree_of_splines+2+model.FlexSurv$Number_knots],0.01)
  axist <- seq(0, model.FlexSurv$knots_time[model.FlexSurv$Degree_of_splines + 2 + model.FlexSurv$Number_knots])
  }
  if(NL==1){
  axisx <- seq(model.FlexSurv$knots_covariates[variableNEW,1], 
               model.FlexSurv$knots_covariates[variableNEW, model.FlexSurv$Degree_of_splines+2+model.FlexSurv$Number_knots],
               0.01)
  }
  
  if (TD==1 & NL==0){
    lines(axist, tdestim1(axist),...)
  }
  
  if (TD==0 & NL==1){
    lines(axisx, nlestim1(axisx) - nlestim1(ref.value.NL),...)
  }
  
  if (TD==1 & NL==1 & TimePoint!=-999){
    lines(axisx, (nlestim1(axisx) - nlestim1(ref.value.NL)) * tdestim1(TimePoint),...)
  }
  
}



DvlpMatrix <- function(data, listeT, ncol, TypeNEW){
  data <- matrix(data, ncol = ncol) # data is a list of full data on each individual, first transform this list to a matrix 
  if (max(data[, TypeNEW[2]]) < min(listeT[listeT != 0])) {  
    XX <- data                      # if max(stop time) is less then min(event time), then data remain the same
  }
  else {                            # if max(stop time) > min(event time)
    aindex <- rep(0, (sum(listeT <= max(data[, TypeNEW[2]])) - 1))
    for (i in 1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1)) { # for the length of unique event times 
      for (j in 1:(dim(data)[1])) { # for each row of an individual, if start time < ith ordered event time<=stop time of row j, stores
                                    # the largest row satisfy this condition in aindex[i]
        if (as.numeric(data[j, TypeNEW[1]]) < as.numeric(listeT[1 + i]) & as.numeric(data[j, TypeNEW[2]]) >= as.numeric(listeT[1 + i]))
          aindex[i] <- j
      }
    }   
    # XX first column repeat ID, column for start day assigns c(0, unique event time-last one), column for stop day assigns unique event time 
    XX <- matrix(nrow = sum(listeT <= max(data[, TypeNEW[2]])) - 1, ncol = ncol)
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), 1] <- rep(data[1, 1], sum(listeT <= max(data[, TypeNEW[2]])) - 1)
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), TypeNEW[1]] <- listeT[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1)]
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), TypeNEW[2]] <- listeT[2:(sum(listeT <= max(data[, TypeNEW[2]])))]
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), TypeNEW[3]] <- c(rep(0,(sum(listeT <= max(data[, TypeNEW[2]])) - 2)), data[dim(data)[1],TypeNEW[3]])
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), -c(1, TypeNEW[1], TypeNEW[2], TypeNEW[3])] <- as.matrix(data[aindex, 
                                                                                                      -c(1, TypeNEW[1], TypeNEW[2], TypeNEW[3])])
    
  }
  X <- XX
  
  list(X)
}



spli <- function(x, j, p, knots){
  if (p == 0) {
    b <- ifelse(x >= knots[j] & x < knots[j + 1], 1, 0)
    return(b)
  }
  else {
    a1 <- ifelse(rep(knots[j] != knots[j + p], length(x)), (x - knots[j])/(knots[j + p] - knots[j]), 0)
    a2 <- ifelse(rep(knots[j + p + 1] != knots[j + 1], length(x)), (knots[j + p + 1] - x)/(knots[j + p + 1] - knots[j + 1]), 0)
    return(a1 * spli(x, j, p - 1, knots) + a2 * spli(x, j + 1, p - 1, knots))
  }
}

