#' DataCreation
#'
#' This function is used to synthesise data for record linkage. It creates 2 data sources of specific sizes, with a common set of records of a specific size, with a certain amount of Partially Identifying Variables (PIVs).
#' For each PIV, we specify the number of unique values, the desired proportion of mistakes and missing values. They can be stable or evolving over time (e.g. representing the postal code).
#' For the unstable PIVs, we can specify the parameter(s) to be used in the survival exponential model that generates changes over time between the records referring to the same entities.
#' When a PIV is unstable, it is later harder to estimate its parameters (probability of mistake vs. probability of change across time). Therefore we may want to enforce estimability in the synthetic data, which we do by
#' enforcing half of the links to have a near-zero time gaps.
#'
#' There are more details to understand the method in our paper, or on the experiments repository of our paper, or in the vignettes.
#'
#' @param PIVs_config A list (of size number of PIVs) where element names are the PIVs and element values are lists with elements: stable (boolean for whether the PIV is stable), conditionalHazard (boolean for whether there are external covariates available to model instability, only required if stable is FALSE), pSameH.cov.A and pSameH.covB (vectors with strings corresponding to the names of the covariates to use to model instability from file A and file B, only required if stable is FALSE, empty vectors may be provided if conditionalHazard is FALSE)
#' @param Nval A vector (of size number of PIVs) with the number fo unique values per PIVs (in the order of the PIVs defined in PIVs_config)
#' @param NRecords A vector (of size 2) with the number of records to be generated in file A and in file B
#' @param Nlinks An integer with the number of links (record referring to the same entities) to be generated
#' @param PmistakesA A vector (of size number of PIVs) with the proportion of mistakes to be generated per PIVs in file A (in the order of the PIVs defined in PIVs_config)
#' @param PmistakesB A vector (of size number of PIVs) with the proportion of mistakes to be generated per PIVs in file B (in the order of the PIVs defined in PIVs_config)
#' @param PmissingA A vector (of size number of PIVs) with the proportion of missing to be generated per PIVs in file A (in the order of the PIVs defined in PIVs_config)
#' @param PmissingB A vector (of size number of PIVs) with the proportion of missing to be generated per PIVs in file A (in the order of the PIVs defined in PIVs_config)
#' @param moving_params A list (of size number of PIVs) where element names are the PIVs and element values are vectors (of size: 1 + number of covariates to use from A + number of covariates to use from B) with the log hazards coefficient (1st one: log baseline hazard, then: the coefficients for conditional hazard covariates from A, then: the coefficients for conditional hazard covariates from B)
#' @param enforceEstimability A boolean value for whether half of the links should have near-0 time gaps (useful for modeling instability and avoiding estimability issues as discussed in the paper)
#'
#' @return A list with generated
#' - dataframe A (encoded: the categorical values of the PIVs are matched to sets of natural numbers),
#' - dataframe B (encoded),
#' - vector of Nvalues (Nval),
#' - vector of TimeDifference (for the links, when thewre is instability),
#' - matrix proba_same_H (number of links, number fo PIVs) with the proba that true values coincide (e.g. 1 - proba of moving)
#' @export
#'
#' @examples
#' PIVs_config = list( V1 = list(stable = TRUE),
#'                     V2 = list(stable = TRUE),
#'                     V3 = list(stable = TRUE),
#'                     V4 = list(stable = TRUE),
#'                     V5 = list( stable = FALSE,
#'                                conditionalHazard = FALSE,
#'                                pSameH.cov.A = c(),
#'                                pSameH.cov.B = c()) )
#' Nval = c(6, 7, 8, 9, 15)
#' NRecords = c(500, 800)
#' Nlinks = 300
#' PmistakesA = c(0.02, 0.02, 0.02, 0.02, 0.02)
#' PmistakesB = c(0.02, 0.02, 0.02, 0.02, 0.02)
#' PmissingA = c(0.007, 0.007, 0.007, 0.007, 0.007)
#' PmissingB = c(0.007, 0.007, 0.007, 0.007, 0.007)
#' moving_params = list(V1=c(),V2=c(),V3=c(),V4=c(),V5=c(0.28))
#' enforceEstimability = TRUE
#' DataCreation( PIVs_config,
#'               Nval,
#'               NRecords,
#'               Nlinks,
#'               PmistakesA,
#'               PmistakesB,
#'               PmissingA,
#'               PmissingB,
#'               moving_params,
#'               enforceEstimability)
DataCreation = function(PIVs_config, Nval, NRecords, Nlinks, PmistakesA, PmistakesB, PmissingA, PmissingB, moving_params, enforceEstimability){

  PIVs = names(PIVs_config)
  PIVs_stable = sapply(PIVs_config, function(x) x$stable)
  PIVs_conditionalHazard_variables = sapply(PIVs_config, function(x) if(!x$stable){x$conditionalHazard}else{FALSE})

  for(i in 1:2)
  {
    dataSet=c()
    for(u in 1:length(Nval))
    {
      # Simulate different probabilities for the true values (xp = exp(0 *(0:(Nval[u]-1))) for uniform distribution of PIVs)
      xp            = exp(0.23 *(0:(Nval[u]-1)))
      probx         = xp/sum(xp)
      dataSet       = cbind(dataSet, sample(1:Nval[u], NRecords[i], replace=TRUE, prob=probx))
    }
    dataSet         = as.data.frame(dataSet)
    names(dataSet)  = PIVs
    assign( paste("dataSet", i, sep=""), dataSet )
  }

  # Add overlapping units
  dataSet2[1:Nlinks,] = dataSet1[1:Nlinks,]

  # Add typos
  for(x in 1:ncol(dataSet1))
  {
    biased = as.logical(stats::rbinom(nrow(dataSet1),1,PmistakesA[x]))
    if(sum(biased)>0)
      dataSet1[,x][biased] = sapply(dataSet1[,x][biased], function(i) sample((1:Nval[x])[-c(i)], 1))
  }

  for(x in 1:ncol(dataSet2))
  {
    biased = as.logical(stats::rbinom(nrow(dataSet2),1,PmistakesB[x]))
    if(sum(biased)>0)
      dataSet2[,x][biased] = sapply(dataSet2[,x][biased], function(i) sample((1:Nval[x])[-c(i)], 1))
  }

  # Add missings
  for(x in 1:ncol(dataSet1))
  {
    biased = as.logical(stats::rbinom(nrow(dataSet1),1,PmissingA[x]))
    if(sum(biased)>0)
      dataSet1[,x][biased] = NA
  }

  for(x in 1:ncol(dataSet2))
  {
    biased = as.logical(stats::rbinom(nrow(dataSet2),1,PmissingB[x]))
    if(sum(biased)>0)
      dataSet2[,x][biased] = NA
  }

  instability = any(!PIVs_stable)
  dataSet1$change                   = FALSE
  dataSet2$change                   = FALSE
  if(instability)
  {
    # Add registration time
    dataSet1$date                     = stats::runif(nrow(dataSet1), 0, 3)
    dataSet2$date                     = stats::runif(nrow(dataSet2), 3, 6)

    if(enforceEstimability){
      nullTimeDiff                      = as.integer(Nlinks/2)
      dataSet1[1:nullTimeDiff, "date"]  = stats::runif(nullTimeDiff, 0.00, 0.01)
      dataSet2[1:nullTimeDiff, "date"]  = stats::runif(nullTimeDiff, 0.00, 0.01)
    }

    TimeDifference                    = abs( dataSet2[1:Nlinks, "date"] - dataSet1[1:Nlinks, "date"] )

    intercept                         = rep(1, Nlinks)
    proba_same_H                      = base::matrix(1,Nlinks,length(PIVs))
    for(k in 1:length(PIVs)){
      if(!PIVs_stable[[k]]){
        # unstable piv
        if(PIVs_conditionalHazard_variables[[k]]){
          # create covariates
          for(c in 1:length(PIVs_config[[k]]$pSameH.cov.A)){
            covariate = PIVs_config[[k]]$pSameH.cov.A[[c]]
            dataSet1[, covariate] = stats::rnorm(nrow(dataSet1), 1, 1)
          }
          for(c in 1:length(PIVs_config[[k]]$pSameH.cov.B)){
            covariate = PIVs_config[[k]]$pSameH.cov.B[[c]]
            dataSet2[, covariate] = stats::rnorm(nrow(dataSet2), 2, 1)
          }
          cov                         = cbind( intercept,
                                               dataSet1[1:Nlinks, PIVs_config[[k]]$pSameH.cov.A, drop=FALSE],
                                               dataSet2[1:Nlinks, PIVs_config[[k]]$pSameH.cov.B, drop=FALSE] )
        }else{
          cov                         = cbind( intercept )
        }

        # generate instability
        for(i in 1:Nlinks)
        {
          is_not_moving = stats::rbinom(1, 1, proba_same_H[i,k])
          if(!is_not_moving)
          {
            if(!is.na(dataSet1[i,x])){
              dataSet2[i,x]         = sample((1:Nval[x])[-c(dataSet1[i,x])], 1)
              dataSet2[i,"change"]  = TRUE
            }
          }
        }
      }
    }
  }

  A = dataSet1
  B = dataSet2

  # Recode the PIVs
  levels_PIVs = lapply(PIVs, function(x) levels(factor(as.character(c(A[,x], B[,x])))))

  for(i in 1:length(PIVs))
  {
    A[,PIVs[i]] = as.numeric(factor(as.character(A[,PIVs[i]]), levels=levels_PIVs[[i]]))
    B[,PIVs[i]] = as.numeric(factor(as.character(B[,PIVs[i]]), levels=levels_PIVs[[i]]))
  }

  Nvalues       = sapply(levels_PIVs, length)

  A$localID     = 1:nrow(A)
  B$localID     = 1:nrow(B)

  A$source      = "A"
  B$source      = "B"

  list(A=A, B=B, Nvalues=Nvalues, TimeDifference=TimeDifference, proba_same_H=proba_same_H)
}

#' createDataAlpha
#'
#' @param nCoefUnstable An integer value with the number of covariates (including the intercept): number of cov from A + number of cov from B + 1.
#' @param stable A boolean value indicating whether the Partially Identifying Variable (PIV) concerned is stable.
#'
#' @return An empty data frame (if stable=FALSE) with nCoefUnstable + 2 columns for book keeping of the elements necessary to update the parameter for instability. There are
#' nCoefUnstable + 2 of those elements: number of cov from A + number of cov from B + intercept + boolean vector indicating where the true values of the records (for the concerned PIV) are equal, vector of time gaps between records
#' @export
#'
#' @examples
#' PIVs_config = list( V1 = list(stable = TRUE),
#'                     V2 = list(stable = TRUE),
#'                     V3 = list(stable = TRUE),
#'                     V4 = list(stable = TRUE),
#'                     V5 = list( stable = FALSE,
#'                                conditionalHazard = FALSE,
#'                                pSameH.cov.A = c(),
#'                                pSameH.cov.B = c()) )
#'
#' PIVs_stable = sapply(PIVs_config, function(x) x$stable)
#'
#' nCoefUnstable = c(0,0,0,0,1)
#'
#' Valpha = mapply( createDataAlpha,
#'                  nCoefUnstable = nCoefUnstable,
#'                  stable = PIVs_stable,
#'                  SIMPLIFY=FALSE)
createDataAlpha <- function(nCoefUnstable, stable){
  nCoef = nCoefUnstable + 2 # cov A + cov B + 1 + Hequal + times
  if (!stable){ return ( data.frame(base::matrix(nrow = 0, ncol = nCoef)) ) } }

#' logPossibleConfig
#'
#' This function helps calculating the number of possible designs for Delta given by nB!/(nB-nLinks)! Needed to compute the log likelihood of the linkage matrix.
#'
#' @param Brecords Number of records in data source B (the largest).
#' @param sumD Number of linked records (at a specific time point of the algorithm).
#'
#' @return The sum of logs of the vector going from nB - nLinks + 1 to nB.
#' @export
#'
#' @examples
#' sumColD = c(0,0,0,0,1,0,0,0,1,0,1,0,1,0,1)
#' links = base::data.frame(list(idxA=c(5,9,11,12,13), idxB=c(5,9,11,13,15)))
#' possconfig = logPossibleConfig(length(sumColD),nrow(links))
logPossibleConfig = function(Brecords,sumD)
{
  return = 0
  if(sumD>0)
    return = sum( log(Brecords:(Brecords-sumD+1))  )
  return
}

#' loglik
#'
#' Log(likelihood) of the linkage matrix.
#'
#' @param LLL A (sparse) matrix with contributions to the complete likelihood of the linked records.
#' @param LLA A vector with contributions to the complete likelihood of the non linked records from A.
#' @param LLB A vector with contributions to the complete likelihood of the non linked records from B.
#' @param links A matrix of 2 columns with indices of the linked records.
#' @param sumRowD A boolean vector indicating, for each row of the linkage matrix, i.e. for each record in the smallest file A, whether the record has a link in B or not.
#' @param sumColD A boolean vector indicating, for each column of the linkage matrix, i.e. for each record in the largest file B, whether the record has a link in A or not.
#' @param gamma The proportion of linked records as a fraction of the smallest file.
#'
#' @return The Log(likelihood) of the linkage matrix.
#' @export
#'
#' @examples
#' LLL = Matrix::Matrix(0, nrow=13, ncol=15, sparse=TRUE)
#' LLA = c(0.001,0.001,0.001,0.001,1.43,0.02,0.007,0.001,2.1,0.0003,1.67,1.5,10)
#' LLB = c(0.001,0.001,0.001,0.001,1.22,0.008,0.01,0.04,3.9,0.0002,1.99,0,2.4,0.009,12)
#' linksR = as.matrix(base::data.frame(list(idxA=c(5,9,11,12,13), idxB=c(5,9,11,13,15))))
#' LLL[linksR] = 0.67
#' sumRowD = c(0,0,0,0,1,0,0,0,1,0,1,1,1)
#' sumColD = c(0,0,0,0,1,0,0,0,1,0,1,0,1,0,1)
#' gamma = 0.5
#' LL0 = loglik(LLL=LLL, LLA=LLA, LLB=LLB, links=linksR, sumRowD=sumRowD, sumColD=sumColD,
#'              gamma=gamma)
loglik = function(LLL, LLA, LLB, links, sumRowD, sumColD, gamma)
{
  # Sanity check for using logPossibleConfig(.) below
  if(length(sumColD) - nrow(links) + 1 <= 0){
    warning("\nThe number of records in B has to be >= to the number of linked records.\nNumber of records in B: ", length(sumColD), "\nNumber of linked records: ", nrow(links), "\nThe problem may come from the fact that file A is bigger than file B, which should not happen.\nNumber of records in A: ", length(sumRowD), "\n")
  }
  logPossD = sum(log(gamma) * sumRowD + log(1-gamma) * (1-sumRowD)) - logPossibleConfig(length(sumColD),nrow(links))
  logPossD + sum(LLA[sumRowD==0]) + sum(LLB[sumColD==0]) + sum(LLL[links])
}

#' simulateH
#'
#' @param data A list with elements:
#' - A: the smallest data source (encoded: the categorical values of the PIVs have to be mapped to sets of natural numbers and missing values are encoded as 0).
#' - B: the largest data source (encoded).
#' - Nvalues: A vector (of size number of PIVs) with the number fo unique values per PIVs (in the order of the PIVs defined in PIVs_config).
#' - PIVs_config: A list (of size number of PIVs) where element names are the PIVs and element values are lists with elements: stable (boolean for whether the PIV is stable), conditionalHazard (boolean for whether there are external covariates available to model instability, only required if stable is FALSE), pSameH.cov.A and pSameH.covB (vectors with strings corresponding to the names of the covariates to use to model instability from file A and file B, only required if stable is FALSE, empty vectors may be provided if conditionalHazard is FALSE).
#' - controlOnMistakes: A vector (of size number of PIVs) of booleans indicating potential bounds on the mistakes probabilities for each PIV. For each PIV, if TRUE there will be control on mistake and the mistake probability will not go above 10%. If FALSE there is no bound on the probability of mistake. WATCH OUT, if you suspect that a variable is unstable but you do not have data to model its dynamics the boolean value should be set to FALSE to allow the parameter for mistake to adapt for the instability. However if you model this instability, the boolean value should be set to TRUE to help the algorithm differenciate the mistakes from the changes over time.
#' - sameMistakes: A boolean value for whether there should be one parameter for the mistakes in A and B or whether each source should have its own parameter. Setting sameMistakes=TRUE is recommended in case of small data sources; the estimation with 2 parameters in that case will fail to capture the mistakes correctly while 1 parameter will be more adapted.
#' - phiMistakesAFixed A vector (of size number of PIVs) of booleans indicating whether the parameters for mistakes should be fixed in case of instability. It should be FALSE, except for unstable PIVs for which it may be set to TRUE in order to avoid estimability problems between the parameter for mistake and the parameter for changes across time.
#' - phiMistakesBFixed A vector (of size number of PIVs) of booleans indicating whether the parameters for mistakes should be fixed in case of instability. It should be FALSE, except for unstable PIVs for which it may be set to TRUE in order to avoid estimability problems between the parameter for mistake and the parameter for changes across time.
#' - phiForMistakesA A vector (of size number of PIVs) of NA or fixed values for the parameters for mistakes. It should be NA, except for unstable PIVs for which one wants to fix the parameter to avoid estimability problem (as indicated with the boolean values in phiMistakesAFixed). In that case it should be set the the expected value for the probability of mistake. If you have no idea: you can put it to 0, the algorithm is quite robust to wrongly fixed parameters.
#' - phiForMistakesB A vector (of size number of PIVs) of NA or fixed values for the parameters for mistakes. It should be NA, except for unstable PIVs for which one wants to fix the parameter to avoid estimability problem (as indicated with the boolean values in phiMistakesBFixed). In that case it should be set the the expected value for the probability of mistake. If you have no idea: you can put it to 0, the algorithm is quite robust to wrongly fixed parameters.
#' @param links A matrix of 2 columns with indices of the linked records.
#' @param survivalpSameH a matrix of size (nrow=number of linked records, ncol=number of PIVs), filled for each PIV in column, with 1 if the PIV is stable and with the probability for true values of the records to coincide as calculate by survival function if the PIV is unstable.
#' @param sumRowD A boolean vector indicating, for each row of the linkage matrix, i.e. for each record in the smallest file A, whether the record has a link in B or not.
#' @param sumColD A boolean vector indicating, for each column of the linkage matrix, i.e. for each record in the largest file B, whether the record has a link in A or not.
#' @param eta The distribution weights for the PIVs.
#' @param phi The proportion of mistakes and missing for the PIVs.
#'
#' @return A list with 2 matrices of the shapes of both data sources, representing the true values of the PIVs underlying the registered values present in the data sources.
#' @export
#'
#' @examples
#' PIVs_config = list( V1 = list(stable = TRUE),
#'                     V2 = list(stable = TRUE),
#'                     V3 = list(stable = TRUE),
#'                     V4 = list(stable = TRUE),
#'                     V5 = list( stable = FALSE,
#'                                conditionalHazard = FALSE,
#'                                pSameH.cov.A = c(),
#'                                pSameH.cov.B = c()) )
#' PIVs = names(PIVs_config)
#' PIVs_stable = sapply(PIVs_config, function(x) x$stable)
#' Nval = c(6, 7, 8, 9, 15)
#' NRecords = c(13, 15)
#' Nlinks = 6
#' PmistakesA = c(0.02, 0.02, 0.02, 0.02, 0.02)
#' PmistakesB = c(0.02, 0.02, 0.02, 0.02, 0.02)
#' PmissingA = c(0.007, 0.007, 0.007, 0.007, 0.007)
#' PmissingB = c(0.007, 0.007, 0.007, 0.007, 0.007)
#' moving_params = list(V1=c(),V2=c(),V3=c(),V4=c(),V5=c(0.28))
#' enforceEstimability = TRUE
#' DATA = DataCreation( PIVs_config,
#'                      Nval,
#'                      NRecords,
#'                      Nlinks,
#'                      PmistakesA,
#'                      PmistakesB,
#'                      PmissingA,
#'                      PmissingB,
#'                      moving_params,
#'                      enforceEstimability)
#' A                    = DATA$A
#' B                    = DATA$B
#' Nvalues              = DATA$Nvalues
#'
#' encodedA = A
#' encodedB = B
#'
#' encodedA[,PIVs][ is.na(encodedA[,PIVs]) ] = 0
#' encodedB[,PIVs][ is.na(encodedB[,PIVs]) ] = 0
#'
#' dataForStEM = list( A                    = encodedA,
#'                     B                    = encodedB,
#'                     Nvalues              = Nvalues,
#'                     PIVs_config          = PIVs_config,
#'                     controlOnMistakes    = c(TRUE, TRUE, TRUE, TRUE, TRUE),
#'                     sameMistakes         = TRUE,
#'                     phiMistakesAFixed    = TRUE,
#'                     phiMistakesBFixed    = TRUE,
#'                     phiForMistakesA      = c(NA, NA, NA, NA, 0),
#'                     phiForMistakesB      = c(NA, NA, NA, NA, 0)
#' )
#'
#' initDeltaMap()
#' linksR = base::matrix(0,0,2)
#' linksCpp = linksR
#' sumRowD = rep(0, nrow(dataForStEM$A))
#' sumColD = rep(0, nrow(dataForStEM$B))
#' nlinkrec = 0
#' survivalpSameH = base::matrix(1, nrow(linksR), length(dataForStEM$Nvalues))
#' gamma = 0.5
#' eta = lapply(dataForStEM$Nvalues, function(x) rep(1/x,x))
#' phi = lapply(dataForStEM$Nvalues, function(x)  c(0.9,0.9,0.1,0.1))
#' nCoefUnstable = lapply( seq_along(PIVs_stable),
#'                         function(idx)
#'                         if(PIVs_stable[idx]){ 0 }
#'                         else{
#'                         ncol(dataForStEM$A[, dataForStEM$PIVs_config[[idx]]$pSameH.cov.A,
#'                              drop=FALSE]) +
#'                         ncol(dataForStEM$B[, dataForStEM$PIVs_config[[idx]]$pSameH.cov.B,
#'                              drop=FALSE]) +
#'                         1 } )
#' alpha = lapply( seq_along(PIVs_stable),
#'                 function(idx) if(PIVs_stable[idx]){ c(-Inf) }else{
#'                   rep(log(0.05), nCoefUnstable[[idx]]) })
#' newTruePivs = simulateH(data=dataForStEM, links=linksCpp, survivalpSameH=survivalpSameH,
#'                         sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
#' truepivsA = newTruePivs$truepivsA
#' truepivsB = newTruePivs$truepivsB
simulateH = function(data, links, survivalpSameH, sumRowD, sumColD, eta, phi)
{
  PIVs = names(data$PIVs_config)
  PIVs_stable = sapply(data$PIVs_config, function(x) x$stable)
  nonlinkedA = sumRowD==0
  nonlinkedB = sumColD==0
  truePIVs = sampleH(nA=dim(data$A[,PIVs]), nB=dim(data$B[,PIVs]), links=links, survivalpSameH=as.matrix(survivalpSameH), pivs_stable=PIVs_stable, pivsA=data$A[,PIVs], pivsB=data$B[,PIVs], nvalues=data$Nvalues, nonlinkedA=nonlinkedA, nonlinkedB=nonlinkedB, eta=eta, phi=phi)
  list(truepivsA=truePIVs$truepivsA, truepivsB=truePIVs$truepivsB)
}

#' simulateD
#'
#' @param data A list with elements:
#' - A: the smallest data source (encoded: the categorical values of the PIVs have to be mapped to sets of natural numbers and missing values are encoded as 0).
#' - B: the largest data source (encoded).
#' - Nvalues: A vector (of size number of PIVs) with the number fo unique values per PIVs (in the order of the PIVs defined in PIVs_config).
#' - PIVs_config: A list (of size number of PIVs) where element names are the PIVs and element values are lists with elements: stable (boolean for whether the PIV is stable), conditionalHazard (boolean for whether there are external covariates available to model instability, only required if stable is FALSE), pSameH.cov.A and pSameH.covB (vectors with strings corresponding to the names of the covariates to use to model instability from file A and file B, only required if stable is FALSE, empty vectors may be provided if conditionalHazard is FALSE).
#' - controlOnMistakes: A vector (of size number of PIVs) of booleans indicating potential bounds on the mistakes probabilities for each PIV. For each PIV, if TRUE there will be control on mistake and the mistake probability will not go above 10%. If FALSE there is no bound on the probability of mistake. WATCH OUT, if you suspect that a variable is unstable but you do not have data to model its dynamics the boolean value should be set to FALSE to allow the parameter for mistake to adapt for the instability. However if you model this instability, the boolean value should be set to TRUE to help the algorithm differenciate the mistakes from the changes over time.
#' - sameMistakes: A boolean value for whether there should be one parameter for the mistakes in A and B or whether each source should have its own parameter. Setting sameMistakes=TRUE is recommended in case of small data sources; the estimation with 2 parameters in that case will fail to capture the mistakes correctly while 1 parameter will be more adapted.
#' - phiMistakesAFixed A vector (of size number of PIVs) of booleans indicating whether the parameters for mistakes should be fixed in case of instability. It should be FALSE, except for unstable PIVs for which it may be set to TRUE in order to avoid estimability problems between the parameter for mistake and the parameter for changes across time.
#' - phiMistakesBFixed A vector (of size number of PIVs) of booleans indicating whether the parameters for mistakes should be fixed in case of instability. It should be FALSE, except for unstable PIVs for which it may be set to TRUE in order to avoid estimability problems between the parameter for mistake and the parameter for changes across time.
#' - phiForMistakesA A vector (of size number of PIVs) of NA or fixed values for the parameters for mistakes. It should be NA, except for unstable PIVs for which one wants to fix the parameter to avoid estimability problem (as indicated with the boolean values in phiMistakesAFixed). In that case it should be set the the expected value for the probability of mistake. If you have no idea: you can put it to 0, the algorithm is quite robust to wrongly fixed parameters.
#' - phiForMistakesB A vector (of size number of PIVs) of NA or fixed values for the parameters for mistakes. It should be NA, except for unstable PIVs for which one wants to fix the parameter to avoid estimability problem (as indicated with the boolean values in phiMistakesBFixed). In that case it should be set the the expected value for the probability of mistake. If you have no idea: you can put it to 0, the algorithm is quite robust to wrongly fixed parameters.
#' @param linksR A matrix of 2 columns with indices of the linked records.
#' @param sumRowD A boolean vector indicating, for each row of the linkage matrix, i.e. for each record in the smallest file A, whether the record has a link in B or not.
#' @param sumColD A boolean vector indicating, for each column of the linkage matrix, i.e. for each record in the largest file B, whether the record has a link in A or not.
#' @param truepivsA A matrix of the shape of data source A, representing the true values of the PIVs underlying the registered values present in A.
#' @param truepivsB A matrix of the shape of data source B, representing the true values of the PIVs underlying the registered values present in B.
#' @param gamma The proportion of linked records as a fraction of the smallest file.
#' @param eta The distribution weights for the PIVs.
#' @param alpha The parameter involved in the survival model for the probability of true values to coincide (parameter for instability).
#' @param phi The proportion of mistakes and missing for the PIVs.
#'
#' @return A list with:
#' - new set of links
#' - new sumRowD
#' - new sumColD
#' - new value of the complete log likelihood
#' - new number fo linked records
#' @export
#'
#' @examples
#' PIVs_config = list( V1 = list(stable = TRUE),
#'                     V2 = list(stable = TRUE),
#'                     V3 = list(stable = TRUE),
#'                     V4 = list(stable = TRUE),
#'                     V5 = list( stable = FALSE,
#'                                conditionalHazard = FALSE,
#'                                pSameH.cov.A = c(),
#'                                pSameH.cov.B = c()) )
#' PIVs = names(PIVs_config)
#' PIVs_stable = sapply(PIVs_config, function(x) x$stable)
#' Nval = c(6, 7, 8, 9, 15)
#' NRecords = c(13, 15)
#' Nlinks = 6
#' PmistakesA = c(0.02, 0.02, 0.02, 0.02, 0.02)
#' PmistakesB = c(0.02, 0.02, 0.02, 0.02, 0.02)
#' PmissingA = c(0.007, 0.007, 0.007, 0.007, 0.007)
#' PmissingB = c(0.007, 0.007, 0.007, 0.007, 0.007)
#' moving_params = list(V1=c(),V2=c(),V3=c(),V4=c(),V5=c(0.28))
#' enforceEstimability = TRUE
#' DATA = DataCreation( PIVs_config,
#'                      Nval,
#'                      NRecords,
#'                      Nlinks,
#'                      PmistakesA,
#'                      PmistakesB,
#'                      PmissingA,
#'                      PmissingB,
#'                      moving_params,
#'                      enforceEstimability)
#' A                    = DATA$A
#' B                    = DATA$B
#' Nvalues              = DATA$Nvalues
#'
#' encodedA = A
#' encodedB = B
#'
#' encodedA[,PIVs][ is.na(encodedA[,PIVs]) ] = 0
#' encodedB[,PIVs][ is.na(encodedB[,PIVs]) ] = 0
#'
#' dataForStEM = list( A                    = encodedA,
#'                     B                    = encodedB,
#'                     Nvalues              = Nvalues,
#'                     PIVs_config          = PIVs_config,
#'                     controlOnMistakes    = c(TRUE, TRUE, TRUE, TRUE, TRUE),
#'                     sameMistakes         = TRUE,
#'                     phiMistakesAFixed    = TRUE,
#'                     phiMistakesBFixed    = TRUE,
#'                     phiForMistakesA      = c(NA, NA, NA, NA, 0),
#'                     phiForMistakesB      = c(NA, NA, NA, NA, 0)
#' )
#'
#' initDeltaMap()
#' linksR = base::matrix(0,0,2)
#' linksCpp = linksR
#' sumRowD = rep(0, nrow(dataForStEM$A))
#' sumColD = rep(0, nrow(dataForStEM$B))
#' nlinkrec = 0
#' survivalpSameH = base::matrix(1, nrow(linksR), length(dataForStEM$Nvalues))
#' gamma = 0.5
#' eta = lapply(dataForStEM$Nvalues, function(x) rep(1/x,x))
#' phi = lapply(dataForStEM$Nvalues, function(x)  c(0.9,0.9,0.1,0.1))
#' nCoefUnstable = lapply( seq_along(PIVs_stable),
#'                         function(idx)
#'                         if(PIVs_stable[idx]){ 0 }
#'                         else{
#'                         ncol(dataForStEM$A[, dataForStEM$PIVs_config[[idx]]$pSameH.cov.A,
#'                              drop=FALSE]) +
#'                         ncol(dataForStEM$B[, dataForStEM$PIVs_config[[idx]]$pSameH.cov.B,
#'                              drop=FALSE]) +
#'                         1 } )
#' alpha = lapply( seq_along(PIVs_stable),
#'                 function(idx) if(PIVs_stable[idx]){ c(-Inf) }else{
#'                   rep(log(0.05), nCoefUnstable[[idx]]) })
#' newTruePivs = simulateH(data=dataForStEM, links=linksCpp, survivalpSameH=survivalpSameH,
#'                         sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
#' truepivsA = newTruePivs$truepivsA
#' truepivsB = newTruePivs$truepivsB
#' Dsample = simulateD(data=dataForStEM, linksR=linksR, sumRowD=sumRowD, sumColD=sumColD,
#'                     truepivsA=truepivsA, truepivsB=truepivsB,
#'                     gamma=gamma, eta=eta, alpha=alpha, phi=phi)
#' linksCpp = Dsample$links
#' linksR = linksCpp + 1
simulateD = function(data, linksR, sumRowD, sumColD, truepivsA, truepivsB, gamma, eta, alpha, phi)
{
  PIVs = names(data$PIVs_config)
  PIVs_stable = sapply(data$PIVs_config, function(x) x$stable)
  #Determine which observation pairs can be matches (or not) based on the true values of the PIVS
  UA = sspaste2(as.matrix(truepivsA[,PIVs_stable]))
  UB = sspaste2(as.matrix(truepivsB[,PIVs_stable]))
  tmpUA = UA
  tmpUB = UB
  valuesU = unique(c(UA,UB))
  UA = as.numeric(factor(UA,levels=valuesU))
  UB = as.numeric(factor(UB,levels=valuesU))
  tmpA = F2(UA, length(valuesU))
  tmpB = F2(UB, length(valuesU))
  select = F33(tmpA, tmpB, length(tmpA))
  pLink = rep(gamma, nrow(data$A[,PIVs]))
  # What should we add/substract from the loglikelihood if an observation is a not linked
  LLA = rep(0,nrow(data$A[,PIVs]))
  for(k in 1:length(data$Nvalues))
  {
    logpTrue = log(eta[[k]])[truepivsA[,k]]
    pMissingA = phi[[k]][3]
    pTypoA = (1-pMissingA) * (1-phi[[k]][1]) / (data$Nvalues[k]-1)
    pAgreeA = (1-pMissingA) * phi[[k]][1]
    # Contribution to the likelihood
    contr = rep(pAgreeA, nrow(data$A[,PIVs]))
    contr[data$A[,PIVs][,k] != truepivsA[,k]] = pTypoA
    contr[data$A[,PIVs][,k] == 0] = pMissingA
    LLA = LLA + logpTrue + log(contr)
  }
  # What should we add/substract from the loglikelihood if an observation is a not linked
  LLB = rep(0,nrow(data$B[,PIVs]))
  for(k in 1:length(data$Nvalues))
  {
    logpTrue = log(eta[[k]])[truepivsB[,k]]
    pMissingB = phi[[k]][4]
    pTypoB = (1-pMissingB) * (1-phi[[k]][2]) / (data$Nvalues[k]-1)
    pAgreeB = (1-pMissingB) * phi[[k]][2]
    # Contribution to the likelihood
    contr = rep(pAgreeB, nrow(data$B[,PIVs]))
    contr[data$B[,PIVs][,k] != truepivsB[,k]] = pTypoB
    contr[data$B[,PIVs][,k] == 0] = pMissingB
    LLB = LLB + logpTrue + log(contr)
  }
  # What do we add/substract from the loglikelihood if a pair is linked
  LLL = Matrix::Matrix(0, nrow=nrow(data$A[,PIVs]), ncol=nrow(data$B[,PIVs]), sparse=TRUE)
  for(k in 1:length(data$Nvalues))
  {
    HA = truepivsA[ select[,1],k ]
    HB = truepivsB[ select[,2],k ]
    logpTrue = log(eta[[k]])[HA]
    pMissingA = phi[[k]][3]
    pTypoA = (1-pMissingA) * (1-phi[[k]][1]) / (data$Nvalues[k]-1)
    pAgreeA = (1-pMissingA) * phi[[k]][1]
    pMissingB = phi[[k]][4]
    pTypoB = (1-pMissingB) * (1-phi[[k]][2]) / (data$Nvalues[k]-1)
    pAgreeB = (1-pMissingB) * phi[[k]][2]
    # Contribution to the likelihood of linked observation from A
    helpA = rep(pAgreeA, length(HA))
    helpA[data$A[,PIVs][select[,1],k] != HA] = pTypoA
    helpA[data$A[,PIVs][select[,1],k] == 0] = pMissingA
    # Contribution to the likelihood of linked observation from B
    helpB = rep(pAgreeB, length(HB))
    helpB[data$B[,PIVs][select[,2],k] != HB] = pTypoB
    helpB[data$B[,PIVs][select[,2],k] == 0] = pMissingB
    LLL[select] = LLL[select] + logpTrue + log(helpA) + log(helpB)
    # Add unstable part if unstable
    if(!PIVs_stable[k])
    {
      times = abs(data$B[data$B$source!="synthetic",][select[,2], "date"] - data$A[data$A$source!="synthetic",][select[,1], "date"])
      intercept = rep(1, nrow(select))
      cov_k = cbind( intercept,
                     data$A[select[,1], data$PIVs_config[[k]]$pSameH.cov.A, drop=FALSE],
                     data$B[select[,2], data$PIVs_config[[k]]$pSameH.cov.B, drop=FALSE] )
      pSameH = SurvivalUnstable(cov_k, alpha[[k]], times)
      helpH = pSameH^(HA==HB) * ((1-pSameH)/(data$Nvalues[k]-1))^(HA!=HB)
      LLL[select] = LLL[select] + log(helpH)
    }
  }
  if(sum(is.na(LLL[select]))>0){
    warning("\nSomething went wrong filling LLL\n")
  }
  # Complete data likelihood
  LL0 = loglik(LLL=LLL, LLA=LLA, LLB=LLB, links=linksR, sumRowD=sumRowD, sumColD=sumColD, gamma=gamma)
  # Single run through D
  Dsample = sampleD(S=as.matrix(select),
                    LLA=LLA,
                    LLB=LLB,
                    LLL=LLL[select],
                    gamma=pLink,
                    loglik=LL0,
                    nlinkrec=as.integer(nrow(linksR)),
                    sumRowD=sumRowD>0,
                    sumColD=sumColD>0)
  linksR = Dsample$links+1
  # Sanity check: does it give the same likelihood?
  if (round(Dsample$loglik, digits = 3) != round(loglik( LLL=LLL, LLA=LLA, LLB=LLB, links=linksR, sumRowD=Dsample$sumRowD, sumColD=Dsample$sumColD, gamma=pLink), digits = 3))
  {
    warning("\nSanity check failed.\nLog likelihood associated with new Delta:", Dsample$loglik, "\nLog likelihood computed on the new Delta:", loglik( LLL=LLL, LLA=LLA, LLB=LLB, links=linksR, sumRowD=Dsample$sumRowD, sumColD=Dsample$sumColD, gamma=pLink), "\n")
  }
  Dsample
}

#' The log likelihood of the survival function with exponential model (-)
#'
#' Log(likelihood) of the survival function with exponential model (as proposed in our paper), representing the probability that true values of a pair of records referring to the same entity coincide. See ?FlexRL::SurvivalUnstable.
#' This function is only used if the PIV is unstable and evolve over time. If so the true values of a linked pair of records may not coincide.
#' If you want to use a different survival function to model instability, you can change the function 'SurvivalUnstable' as well as this function 'loglikSurvival'.
#'
#' In our Stochastic Expectation Maximisation (StEM) algorithm (see ?FlexRL::StEM) we minimise - log(likelihood), which is equivalent to maximise log(likelihood). Therefore this function actually returns (and should return if you create your own) the opposite (-) of the log(likelihood) associated with the survival function defining the probabilities that true values coincide.
#'
#' @param alphas A vector of size 1+cov in A+cov in B with coefficients of the hazard (baseline hazard and conditional hazard)
#' @param X A matrix with number of linked records rows and 1+cov in A+cov in B columns (first column: intercept, following columns: covariates from A and then from B to model instability) (used for optimisation: X concatenate the X obtained in each iteration of the Gibbs sampler)
#' @param times A vector of size number of linked records with the time gaps between the record from each sources (used for optimisation: times concatenate the times vectors obtained in each iteration of the Gibbs sampler)
#' @param Hequal A vector of size number of linked records with boolean values indicating wether the values in A and in B coincide (used for optimisation: times concatenate the times vectors obtained in each iteration of the Gibbs sampler)
#'
#' @return The value of the opposite (-) of the log(likelihood) associated with the survival function defining the probabilities that true values coincide (as defined in the paper) (the algorithm minimises -log(likelihood) i.e. maximises the log(likelihood)).
#' @export
#'
#' @examples
#' nCoefUnstable = 1
#' alphaInit = rep(-0.05, nCoefUnstable)
#' Valpha = base::data.frame(list(cov=c(2,2.1,3.4,2.5,2.9),
#'                                times=c(0.001,0.2,1.3,1.5,2),
#'                                Hequal=c(TRUE, TRUE, TRUE, FALSE, FALSE)))
#' X = Valpha[,1:nCoefUnstable]
#' times = Valpha$times
#' Hequal = Valpha$Hequal
#' optim = stats::nlminb(alphaInit, loglikSurvival, control=list(trace=FALSE),
#'                       X=X, times=times, Hequal=Hequal)
#' alpha = optim$par
loglikSurvival = function(alphas, X, times, Hequal)
{
  loglikelihood = sum( - exp( as.matrix(X) %*% alphas ) * times + (!Hequal) * log( exp( exp( as.matrix(X) %*% alphas ) * times ) - 1 ) )
  - loglikelihood
}

#' The survival function with exponential model
#'
#' The survival function with exponential model (as proposed in our paper) is representing the probability that true values of a pair of records referring to the same entity coincide. For a linked pair of record (i,j) from sources A, B respectively, P(HAik = HBjk | t_ij, ...) = exp( - exp(X.ALPHA) t_ij ).
#' This function is only used if the PIV is unstable and evolve over time. If so the true values of a linked pair of records may not coincide.
#' If you want to use a different survival function to model instability, you can change this function 'SurvivalUnstable' as well as the associated log(likelihood) function 'loglikSurvival'.
#' Also see ?FlexRL::loglikSurvival.
#'
#' The simplest model (without covariates) just writes P(HAik = HBjk | t_ij, ...) = exp( - exp(alpha) . t_ij )
#' The more complex model (with covariates) writes P(HAik = HBjk | t_ij, ...) = exp( - exp(X.ALPHA) . t_ij ) and uses a matrix X (nrow=nbr of linked records, ncol=1 + nbr of cov from A + nbr of cov from B) where the first column is filled with 1 (intercept) and the subsequent columns are the covariates values from source A and/or from source B to be used. The ALPHA in this case is a vector of parameters, the first one being associated with the intercept is the same one than for the simplest model, the subsequent ones are associated with the covariates from A and/or from B.
#'
#' @param Xlinksk A matrix with number of linked records rows and 1+cov in A+cov in B columns, with a first column filled with 1 (intercept), and following columns filled with the values of the covariates useful for modelling instability for the linked records
#' @param alphask A vector of size 1+cov in A+cov in B, with as first element the baseline hazard and following elements being the coefficient of the conditional hazard associated with the covariates given in X
#' @param times A vector of size number of linked records with the time gaps between the record from each sources
#'
#' @return A vector (for an unstable PIV) of size number of linked records with the probabilities that true values coincide (e.g. 1 - proba to move if the PIV is postal code) defined according to the survival function with exponential model proposed in the paper
#' @export
#'
#' @examples
#' nCoefUnstable = 1
#' intercept = rep(1,5)
#' cov_k = cbind( intercept )
#' times = c(0.001,0.2,1.3,1.5,2)
#' survivalpSameH = SurvivalUnstable(cov_k, log(0.28), times)
SurvivalUnstable = function(Xlinksk, alphask, times)
{
  # returns the proba that true values of PIV indexed by k coincide
  # P(HAik = HBjk | t_ij) = exp( - exp(X.ALPHA) t_ij )
  # first column in Xlinks is the intercept
  # following columns in Xlinks are the covariates to include
  # number of rows in Xlinks is the number of linked records
  exp( - exp( as.matrix(Xlinksk) %*% alphask ) * times )
}

#' Stochastic Expectation Maximisation (StEM) for Record Linkage
#'
#' @param data A list with elements:
#' - A: the smallest data source (encoded: the categorical values of the Partially Identifying Variables (PIVs) have to be mapped to sets of natural numbers and missing values are encoded as 0).
#' - B: the largest data source (encoded).
#' - Nvalues: A vector (of size number of PIVs) with the number fo unique values per PIVs (in the order of the PIVs defined in PIVs_config).
#' - PIVs_config: A list (of size number of PIVs) where element names are the PIVs and element values are lists with elements: stable (boolean for whether the PIV is stable), conditionalHazard (boolean for whether there are external covariates available to model instability, only required if stable is FALSE), pSameH.cov.A and pSameH.covB (vectors with strings corresponding to the names of the covariates to use to model instability from file A and file B, only required if stable is FALSE, empty vectors may be provided if conditionalHazard is FALSE).
#' - controlOnMistakes: A vector (of size number of PIVs) of booleans indicating potential bounds on the mistakes probabilities for each PIV. For each PIV, if TRUE there will be control on mistake and the mistake probability will not go above 10%. If FALSE there is no bound on the probability of mistake. WATCH OUT, if you suspect that a variable is unstable but you do not have data to model its dynamics the boolean value should be set to FALSE to allow the parameter for mistake to adapt for the instability. However if you model this instability, the boolean value should be set to TRUE to help the algorithm differenciate the mistakes from the changes over time.
#' - sameMistakes: A boolean value for whether there should be one parameter for the mistakes in A and B or whether each source should have its own parameter. Setting sameMistakes=TRUE is recommended in case of small data sources; the estimation with 2 parameters in that case will fail to capture the mistakes correctly while 1 parameter will be more adapted.
#' - phiMistakesAFixed A vector (of size number of PIVs) of booleans indicating whether the parameters for mistakes should be fixed in case of instability. It should be FALSE, except for unstable PIVs for which it may be set to TRUE in order to avoid estimability problems between the parameter for mistake and the parameter for changes across time.
#' - phiMistakesBFixed A vector (of size number of PIVs) of booleans indicating whether the parameters for mistakes should be fixed in case of instability. It should be FALSE, except for unstable PIVs for which it may be set to TRUE in order to avoid estimability problems between the parameter for mistake and the parameter for changes across time.
#' - phiForMistakesA A vector (of size number of PIVs) of NA or fixed values for the parameters for mistakes. It should be NA, except for unstable PIVs for which one wants to fix the parameter to avoid estimability problem (as indicated with the boolean values in phiMistakesAFixed). In that case it should be set the the expected value for the probability of mistake. If you have no idea: you can put it to 0, the algorithm is quite robust to wrongly fixed parameters.
#' - phiForMistakesB A vector (of size number of PIVs) of NA or fixed values for the parameters for mistakes. It should be NA, except for unstable PIVs for which one wants to fix the parameter to avoid estimability problem (as indicated with the boolean values in phiMistakesBFixed). In that case it should be set the the expected value for the probability of mistake. If you have no idea: you can put it to 0, the algorithm is quite robust to wrongly fixed parameters.
#' @param StEMIter An integer with the total number of iterations of the Stochastic
#' EM algorithm (including the period to discard as burn-in)
#' @param StEMBurnin An integer with the number of iterations to discard as burn-in
#' @param GibbsIter An integer with the total number of iterations of the Gibbs sampler
#' (done in each iteration of the StEM) (including the period to discard as burn-in)
#' @param GibbsBurnin An integer with the number of iterations to discard as burn-in
#' @param musicOn A boolean value, if TRUE the algorithm will play music at the end of
#' the algorithm, useful if you have to wait for the record linkage to run and to act as
#' an alarm when record linkage is done
#' @param newDirectory A NULL value or: A string with the name of (or path to) the directory
#' (which should already exist) where to save the environment variables at the end of each
#' iteration (useful when record linkage is very long, to not loose everything and not restart
#' from scratch in case your computer shut downs before record linkage is finished)
#' @param saveInfoIter A boolean value to indicate whether you want the environment variables
#' to be saved at the end of each iteration (useful when record linkage is very long, to not
#' loose everything and not restart from scratch in case your computer shut downs before
#' record linkage is finished)
#'
#' @return A list with:
#' - Delta, the summary of a sparse matrix, i.e. a data frame with 3 columns: the indices from the first data source A, the indices from the second data source B, the non-zero probability that the records associated with this pair of indices are linked (i.e. the
#' posterior probabilities to be linked). One has to select the pairs where this probability > 0.5 to get a valid set of linked records, (this threshold on the linkage probability is necessary to ensure the one-to-one assignment constraint of record linkage stating that one record in one file can at most be linked to one record in the other file).
#' - gamma, a vector with the chain of the parameter gamma representing
#' the proportion of linked records as a fraction of the smallest file,
#' - eta, a vector with the
#' chain of the parameter eta representing the distribution of the PIVs,
#' - alpha, a vector with
#' the chain of the parameter alpha representing the hazard coefficient of the model for instability,
#' - phi, a vector with the chain of the parameter phi representing the registration errors parameters).
#'
#' There are more details to understand the method in our paper, or on the experiments repository of our paper, or in the vignettes.
#' @export
#'
#' @examples
#' \donttest{
#' PIVs_config = list( V1 = list(stable = TRUE),
#'                     V2 = list(stable = TRUE),
#'                     V3 = list(stable = TRUE),
#'                     V4 = list(stable = TRUE),
#'                     V5 = list( stable = FALSE,
#'                                conditionalHazard = FALSE,
#'                                pSameH.cov.A = c(),
#'                                pSameH.cov.B = c()) )
#' PIVs = names(PIVs_config)
#' PIVs_stable = sapply(PIVs_config, function(x) x$stable)
#' Nval = c(6, 7, 8, 9, 15)
#' NRecords = c(500, 800)
#' Nlinks = 300
#' PmistakesA = c(0.02, 0.02, 0.02, 0.02, 0.02)
#' PmistakesB = c(0.02, 0.02, 0.02, 0.02, 0.02)
#' PmissingA = c(0.007, 0.007, 0.007, 0.007, 0.007)
#' PmissingB = c(0.007, 0.007, 0.007, 0.007, 0.007)
#' moving_params = list(V1=c(),V2=c(),V3=c(),V4=c(),V5=c(0.28))
#' enforceEstimability = TRUE
#' DATA = DataCreation( PIVs_config,
#'                      Nval,
#'                      NRecords,
#'                      Nlinks,
#'                      PmistakesA,
#'                      PmistakesB,
#'                      PmissingA,
#'                      PmissingB,
#'                      moving_params,
#'                      enforceEstimability)
#' A                    = DATA$A
#' B                    = DATA$B
#' Nvalues              = DATA$Nvalues
#'
#' encodedA = A
#' encodedB = B
#'
#' encodedA[,PIVs][ is.na(encodedA[,PIVs]) ] = 0
#' encodedB[,PIVs][ is.na(encodedB[,PIVs]) ] = 0
#'
#' data = list( A = encodedA,
#'              B = encodedB,
#'              Nvalues = Nvalues,
#'              PIVs_config = PIVs_config,
#'              controlOnMistakes = c(TRUE,TRUE,FALSE,FALSE,FALSE),
#'              sameMistakes = TRUE,
#'              phiMistakesAFixed = FALSE,
#'              phiMistakesBFixed = FALSE,
#'              phiForMistakesA = c(NA,NA,NA,NA,NA),
#'              phiForMistakesB = c(NA,NA,NA,NA,NA))
#'  fit = stEM( data = data,
#'              StEMIter = 50,
#'              StEMBurnin = 30,
#'              GibbsIter = 50,
#'              GibbsBurnin = 30,
#'              musicOn = TRUE,
#'              newDirectory = NULL,
#'              saveInfoIter = FALSE )
#' }
stEM = function(data, StEMIter, StEMBurnin, GibbsIter, GibbsBurnin, musicOn=TRUE, newDirectory=NULL, saveInfoIter=FALSE)
{

  PIVs = names(data$PIVs_config)
  PIVs_stable = sapply(data$PIVs_config, function(x) x$stable)
  PIVs_conditionalHazard_variables = sapply(data$PIVs_config, function(x) if(!x$stable){x$conditionalHazard}else{FALSE})

  # data$A, data$B, data$PIVs_stable, PIVs are used globally without necessarily being passed as parameters

  # ANY CONTROL ON MISTAKES? i.e. BOUND ON PHI
  if(any(data$controlOnMistakes)){
    controlMistakes = which(data$controlOnMistakes)
  }else{
    controlMistakes = c()
  }

  # ANY UNSTABLE PIV?
  instability = any(!PIVs_stable)
  if(instability){
    unstablePIVs = which(!PIVs_stable)
    if(any(PIVs_conditionalHazard_variables)){
      for(k in 1:length(unstablePIVs)){
        unstablePIV = unstablePIVs[k]
        if( length(data$PIVs_config[[unstablePIV]]$pSameH.cov.A)>0 )
          testit::assert( "Some variables from A to include in the survival model for unstable PIVs do not exist.", data$PIVs_config[[unstablePIV]]$pSameH.cov.A %in% colnames(data$A) )
        if( length(data$PIVs_config[[unstablePIV]]$pSameH.cov.B)>0 )
          testit::assert( "Some variables from B to include in the survival model for unstable PIVs do not exist.", data$PIVs_config[[unstablePIV]]$pSameH.cov.B %in% colnames(data$B) )
      }
    }
  }

  nGibbsIter = GibbsIter - GibbsBurnin
  testit::assert( "Number of iterations for StEM should be positive.", StEMIter - StEMBurnin > 0 )
  testit::assert( "Number of iterations for StEM Gibbs sampler should be positive.", nGibbsIter > 0 )

  # Parameters for PIVs

  # Parameter for the probability for a pair to be linked
  gamma = 0.5

  # Parameters for the distributions of the true values
  eta = lapply(data$Nvalues, function(x) rep(1/x,x))

  # Parameters for the survival model describing pSameH for unstable PIVs (over time and potentially more covariates)
  # Number of coefficients: baseline + covariates from A + covariates from B
  nCoefUnstable = lapply(seq_along(PIVs_stable), function(idx) if(PIVs_stable[idx]){ 0 }else{ ncol(data$A[, data$PIVs_config[[idx]]$pSameH.cov.A, drop=FALSE]) + ncol(data$B[, data$PIVs_config[[idx]]$pSameH.cov.B, drop=FALSE]) + 1 } )
  alpha = lapply(seq_along(PIVs_stable), function(idx) if(PIVs_stable[idx]){ c(-Inf) }else{ rep(log(0.05), nCoefUnstable[[idx]]) })

  # Parameters for the registration errors (agreement in A, agreement in B, missing in A, missing in B)
  phi = lapply(data$Nvalues, function(x)  c(0.9,0.9,0.1,0.1))

  NmissingA = lapply(seq_along(data$Nvalues), function(k) sum(data$A[,PIVs][,k]==0))
  NmissingB = lapply(seq_along(data$Nvalues), function(k) sum(data$B[,PIVs][,k]==0))

  gamma.iter = array(NA, c(StEMIter, length(gamma)))
  eta.iter = lapply(data$Nvalues, function(x) array(NA, c(StEMIter, x)))
  phi.iter = lapply(data$Nvalues, function(x) array(NA, c(StEMIter, 4)))
  alpha.iter = lapply(nCoefUnstable, function(x) array(NA, c(StEMIter, x)))

  time.iter=c()
  pb = progress::progress_bar$new(format = "Running StEM algorithm [:bar] :percent in :elapsed",       total = StEMIter, clear = FALSE, width= 60)

  # Stochastic EM
  for(iter in 1:StEMIter)
  {
    tijdM = Sys.time()

    initDeltaMap()
    linksR = base::matrix(0,0,2)
    linksCpp = linksR
    sumRowD = rep(0, nrow(data$A))
    sumColD = rep(0, nrow(data$B))
    nlinkrec = 0
    survivalpSameH = base::matrix(1, nrow(linksR), length(data$Nvalues))

    # if burnin value is 0, the algorithm will explore the necessary burnin for the number of linked records to stagnate
    countBurnin = 10
    Burnin_total = 0

    # Burn-in period Gibbs sampler
    if(GibbsBurnin == 0){
      pb_burnin_explore = progress::progress_bar$new(format = "(:spin)     Burn-in period: :whatburnin     Number of linked records: :whatlinkrec", total = 300, clear = FALSE, width= 100)
      while(countBurnin != 0)
      {
        Burnin_total = Burnin_total + 1
        newTruePivs = simulateH(data=data, links=linksCpp, survivalpSameH=survivalpSameH, sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
        truepivsA = newTruePivs$truepivsA
        truepivsB = newTruePivs$truepivsB
        Dsample = simulateD(data=data, linksR=linksR, sumRowD=sumRowD, sumColD=sumColD, truepivsA=truepivsA, truepivsB=truepivsB, gamma=gamma, eta=eta, alpha=alpha, phi=phi)

        linksCpp = Dsample$links
        linksR = linksCpp + 1
        sumRowD = Dsample$sumRowD
        sumColD = Dsample$sumColD
        loglikelihood = Dsample$loglik
        new_nlinkred = Dsample$nlinkrec

        if( abs(new_nlinkred - nlinkrec) < 10 ){
          countBurnin = countBurnin - 1
          nlinkrec = new_nlinkred
        }else{
          countBurnin = 15
          nlinkrec = new_nlinkred
        }

        survivalpSameH = base::matrix(1, nrow(linksR), length(data$Nvalues))
        if(instability){
          if(nrow(linksR)>0){
            times = abs(data$B[data$B$source!="synthetic",][linksR[,2], "date"] - data$A[data$A$source!="synthetic",][linksR[,1], "date"])
            intercept = rep(1, nrow(linksR))
            for(k in 1:length(PIVs)){
              if(k %in% unstablePIVs){
                cov_k = cbind( intercept,
                               data$A[linksR[,1], data$PIVs_config[[k]]$pSameH.cov.A, drop=FALSE],
                               data$B[linksR[,2], data$PIVs_config[[k]]$pSameH.cov.B, drop=FALSE] )
                survivalpSameH[,k] = SurvivalUnstable(cov_k, alpha[[k]], times)
              }
            }
          }
        }

        pb_burnin_explore$tick(tokens = list(whatburnin = Burnin_total, whatlinkrec = nlinkrec))
        if(Burnin_total>=300){
          countBurnin = 0
          warning("\nThe burn-in period exceeded 300 iterations, this is the default maximum burn-in in that case.\nIf you want to increase the burn-in period fix it with the parameter 'GibbsBurnin' of the StEM function.\n")
          pb_burnin_explore$terminate()
          GibbsBurnin = 300
          stop()
        }
      }
    }else{
      for(j in 1:GibbsBurnin)
      {
        newTruePivs = simulateH(data=data, links=linksCpp, survivalpSameH=survivalpSameH, sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
        truepivsA = newTruePivs$truepivsA
        truepivsB = newTruePivs$truepivsB
        Dsample = simulateD(data=data, linksR=linksR, sumRowD=sumRowD, sumColD=sumColD, truepivsA=truepivsA, truepivsB=truepivsB, gamma=gamma, eta=eta, alpha=alpha, phi=phi)

        linksCpp = Dsample$links
        linksR = linksCpp + 1
        sumRowD = Dsample$sumRowD
        sumColD = Dsample$sumColD
        loglikelihood = Dsample$loglik
        nlinkred = Dsample$nlinkrec

        survivalpSameH = base::matrix(1, nrow(linksR), length(data$Nvalues))
        if(instability){
          if(nrow(linksR)>0){
            times = abs(data$B[data$B$source!="synthetic",][linksR[,2], "date"] - data$A[data$A$source!="synthetic",][linksR[,1], "date"])
            intercept = rep(1, nrow(linksR))
            for(k in 1:length(PIVs)){
              if(k %in% unstablePIVs){
                cov_k = cbind( intercept,
                               data$A[linksR[,1], data$PIVs_config[[k]]$pSameH.cov.A, drop=FALSE],
                               data$B[linksR[,2], data$PIVs_config[[k]]$pSameH.cov.B, drop=FALSE] )
                survivalpSameH[,k] = SurvivalUnstable(cov_k, alpha[[k]], times)
              }
            }
          }
        }
      }
    }

    # Administration for the M-step
    Vgamma = c()
    Veta = lapply(data$Nvalues, function(x) c())
    Valpha = mapply(createDataAlpha, nCoefUnstable = nCoefUnstable, stable = PIVs_stable, SIMPLIFY=FALSE)
    Vphi = lapply(data$Nvalues, function(x) c())

    for(j in 1:nGibbsIter)
    {
      newTruePivs = simulateH(data=data, links=linksCpp, survivalpSameH=survivalpSameH, sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
      truepivsA = newTruePivs$truepivsA
      truepivsB = newTruePivs$truepivsB
      Dsample = simulateD(data=data, linksR=linksR, sumRowD=sumRowD, sumColD=sumColD, truepivsA=truepivsA, truepivsB=truepivsB, gamma=gamma, eta=eta, alpha=alpha, phi=phi)

      linksCpp = Dsample$links
      linksR = linksCpp + 1
      sumRowD = Dsample$sumRowD
      sumColD = Dsample$sumColD
      loglikelihood = Dsample$loglik
      nlinkrec = Dsample$nlinkrec

      survivalpSameH = base::matrix(1, nrow(linksR), length(data$Nvalues))
      if(instability){
        if(nrow(linksR)>0){
          times = abs(data$B[data$B$source!="synthetic",][linksR[,2], "date"] - data$A[data$A$source!="synthetic",][linksR[,1], "date"])
          intercept = rep(1, nrow(linksR))
          for(k in 1:length(PIVs)){
            if(k %in% unstablePIVs){
              cov_k = cbind( intercept,
                             data$A[linksR[,1], data$PIVs_config[[k]]$pSameH.cov.A, drop=FALSE],
                             data$B[linksR[,2], data$PIVs_config[[k]]$pSameH.cov.B, drop=FALSE] )
              survivalpSameH[,k] = SurvivalUnstable(cov_k, alpha[[k]], times)
              Hequal = truepivsA[linksR[,1],k] == truepivsB[linksR[,2],k]
              Vtmp_k = cbind( cov_k,
                              times,
                              Hequal )
              Valpha[[k]] = rbind(Valpha[[k]], Vtmp_k)
            }
          }
        }
      }

      # Update gamma
      Vgamma[j] = nrow(linksR)
      if(nrow(linksR)==0){
        warning("\nloglikelihood = ", loglikelihood, "\nSo far, gamma = ", gamma, "\nIn the last iteration no link has been made.\nThis is probably due to a difference in the support of some PIV between file A and file B.\n")
      }

      # Update eta
      for(k in 1:length(data$Nvalues))
      {
        facpivsA = factor(truepivsA[,k], levels=1:data$Nvalues[k])
        facpivsB = factor(truepivsB[,k], levels=1:data$Nvalues[k])
        Veta[[k]] = rbind(Veta[[k]], table(facpivsA[sumRowD==0]) + table(facpivsB[sumColD==0]) + table(facpivsA[sumRowD==1]))
      }

      # Update phi: count agreements / missings
      for(k in 1:length(data$Nvalues)){
        agreements_in_A = sum(truepivsA[,k]==data$A[,PIVs][,k])
        agreements_in_B = sum(truepivsB[,k]==data$B[,PIVs][,k])
        Vphi[[k]] = rbind( Vphi[[k]],
                           c( agreements_in_A,
                              agreements_in_B,
                              NmissingA[[k]],
                              NmissingB[[k]]) )
      }
    }

    # Calculate new parameters gamma/eta/phi/alpha

    # New gamma
    gamma = sum(Vgamma) / (nGibbsIter*nrow(truepivsA))

    # New eta
    for(k in 1:length(data$Nvalues)){
      eta[[k]] = colSums(Veta[[k]])/sum(Veta[[k]])
      if(sum(eta[[k]]==0)>0){
        warning("\nAn error may appear due to some rare values in PIV", k, "\nSo far, eta ", k, ":\n", eta[[k]], "\nCount values of PIV", k, "in A: \n", table(data$A[,PIVs][,k]), "\nand in B: \n", table(data$B[,PIVs][,k]), "\n")
      }
    }

    # New alpha
    if(instability){
      if(nrow(linksR)>0){
        for(k in 1:length(PIVs)){
          if(!PIVs_stable[k]){
            alphaInit = rep(-0.05, nCoefUnstable[[k]])
            X = Valpha[[k]][,1:nCoefUnstable[[k]]]
            times = Valpha[[k]]$times
            Hequal = Valpha[[k]]$Hequal
            optim = stats::nlminb(alphaInit, loglikSurvival, control=list(trace=FALSE), X=X, times=times, Hequal=Hequal)
            alpha[[k]] = optim$par
          }
        }
      }
    }

    # New phi
    for(k in 1:length(data$Nvalues))
    {
      NtotalA = nGibbsIter * nrow(data$A)
      NtotalB = nGibbsIter * nrow(data$B)
      Ntotal  = nGibbsIter * (nrow(data$A) + nrow(data$B))
      if(k %in% controlMistakes){
        phi[[k]][1] = max( 0.9, sum(Vphi[[k]][,1]) / (NtotalA - sum(Vphi[[k]][,3])) )
        phi[[k]][2] = max( 0.9, sum(Vphi[[k]][,2]) / (NtotalB - sum(Vphi[[k]][,4])) )
      }else{
        phi[[k]][1] = sum(Vphi[[k]][,1]) / (NtotalA - sum(Vphi[[k]][,3]))
        phi[[k]][2] = sum(Vphi[[k]][,2]) / (NtotalB - sum(Vphi[[k]][,4]))
      }
      if(data$sameMistakes){
        phi[[k]][1] = (sum(Vphi[[k]][,1]) + sum(Vphi[[k]][,2])) / (Ntotal - sum(Vphi[[k]][,3]) - sum(Vphi[[k]][,4]))
        phi[[k]][2] = (sum(Vphi[[k]][,1]) + sum(Vphi[[k]][,2])) / (Ntotal - sum(Vphi[[k]][,3]) - sum(Vphi[[k]][,4]))
      }
      if(data$phiMistakesAFixed | data$phiMistakesBFixed){
        for(i in 1:length(unstablePIVs)){
          if(data$phiMistakesAFixed){
            phi[[ unstablePIVs[i] ]][1] = 1-data$phiForMistakesA[[ unstablePIVs[i] ]]
          }
          if(data$phiMistakesBFixed){
            phi[[ unstablePIVs[i] ]][2] = 1-data$phiForMistakesB[[ unstablePIVs[i] ]]
          }
        }
      }
      phi[[k]][3] = sum(Vphi[[k]][,3]) / (nGibbsIter*nrow(data$A))
      phi[[k]][4] = sum(Vphi[[k]][,4]) / (nGibbsIter*nrow(data$B))
    }

    # Administration

    gamma.iter[iter,] = gamma

    for(k in 1:length(data$Nvalues))
    {
      eta.iter[[k]][iter,] = eta[[k]]
      alpha.iter[[k]][iter,] = alpha[[k]]
      phi.iter[[k]][iter,] = phi[[k]]
    }

    pb$tick()

    # save current environment, global variables and local variables
    if(!is.null(newDirectory) & saveInfoIter)
    {
      save.image(file=file.path(newDirectory, 'myEnvironment.RData'))
      save(iter, gamma, eta, alpha, phi, gamma.iter, eta.iter, phi.iter, alpha.iter, file=file.path(newDirectory, 'myEnvironmentLocal.RData'))
    }
  }

  pb$terminate()

  Delta = Matrix::Matrix(0, nrow=nrow(data$A), ncol=nrow(data$B), sparse=TRUE)

  gamma_avg = apply(gamma.iter, 2, function(x) mean(x[StEMBurnin:StEMIter , drop=FALSE]))
  eta_avg = lapply(eta.iter, function(x) apply(x[StEMBurnin:StEMIter, , drop=FALSE], 2, mean))
  alpha_avg = lapply(alpha.iter, function(x) apply(x[StEMBurnin:StEMIter, , drop=FALSE], 2, mean))
  phi_avg = lapply(phi.iter, function(x) apply(x[StEMBurnin:StEMIter, , drop=FALSE], 2, mean))

  pbfinal = progress::progress_bar$new(format = "Drawing Delta          [:bar] :percent in :elapsed",       total = 1000, clear = FALSE, width= 60)

  for(m in 1:1000)
  {
    newTruePivs = simulateH(data=data, links=linksCpp, survivalpSameH=survivalpSameH, sumRowD=sumRowD, sumColD=sumColD, eta=eta_avg, phi=phi_avg)
    truepivsA = newTruePivs$truepivsA
    truepivsB = newTruePivs$truepivsB
    Dsample = simulateD(data=data, linksR=linksR, sumRowD=sumRowD, sumColD=sumColD, truepivsA=truepivsA, truepivsB=truepivsB, gamma=gamma_avg, eta=eta_avg, alpha=alpha_avg, phi=phi_avg)

    linksCpp = Dsample$links
    linksR = linksCpp+1
    sumRowD = Dsample$sumRowD
    sumColD = Dsample$sumColD
    loglikelihood = Dsample$loglik
    nlinkrec = Dsample$nlinkrec

    survivalpSameH = base::matrix(1, nrow(linksR), length(data$Nvalues))
    if(instability){
      if(nrow(linksR)>0){
        times = abs(data$B[data$B$source!="synthetic",][linksR[,2], "date"] - data$A[data$A$source!="synthetic",][linksR[,1], "date"])
        intercept = rep(1, nrow(linksR))
        for(k in 1:length(PIVs)){
          if(k %in% unstablePIVs){
            cov_k = cbind( intercept,
                           data$A[linksR[,1], data$PIVs_config[[k]]$pSameH.cov.A, drop=FALSE],
                           data$B[linksR[,2], data$PIVs_config[[k]]$pSameH.cov.B, drop=FALSE] )
            survivalpSameH[,k] = SurvivalUnstable(cov_k, alpha[[k]], times)
          }
        }
      }
    }

    for(l in 1:nrow(linksR)){
      i = linksR[l,1]
      j = linksR[l,2]
      Delta[i,j] = Delta[i,j] + 1
    }
    pbfinal$tick()
  }
  Delta = Delta / 1000

  pbfinal$terminate()

  if(musicOn){
    utils::browseURL('https://www.youtube.com/watch?v=NTa6Xbzfq1U')
  }

  list(Delta=as.data.frame(Matrix::summary(Delta)), gamma=gamma.iter, eta=eta.iter, alpha=alpha.iter, phi=phi.iter)
}

#' launchNaive
#'
#' @param PIVs A vector of size the number of Partially Identifying Variables (PIVs) with their names (as columns names in the data sources).
#' @param encodedA One data source (encoded: the categorical values of the PIVs have to be mapped to sets of natural numbers and missing values are encoded as 0).
#' @param encodedB The other data source (encoded).
#'
#' @return
#' The linkage set, a dataframe of 2 columns with indices from the first data source (A) and the second one (B) for which all the PIVs (when non-missing) matches in their values. When a PIV is missing, this method will match the record to all the records in the other file for which all other values of the PIVs match.
#' Therefore, this 'naive' (or 'simplistic') method does not enforce the one-to-one assignment constraint of record linkage (one record in one file can at most be linked to one record in the other file).
#' This method should only be used to judge the difficulty of the record linkage task: it gives information about the amount of duplicates between files and the discriminative power of all the PIVs together as a way to link the records.
#' @export
#'
#' @examples
#' PIVs_config = list( V1 = list(stable = TRUE),
#'                     V2 = list(stable = TRUE),
#'                     V3 = list(stable = TRUE),
#'                     V4 = list(stable = TRUE),
#'                     V5 = list( stable = FALSE,
#'                                conditionalHazard = FALSE,
#'                                pSameH.cov.A = c(),
#'                                pSameH.cov.B = c()) )
#' PIVs = names(PIVs_config)
#' PIVs_stable = sapply(PIVs_config, function(x) x$stable)
#' Nval = c(6, 7, 8, 9, 15)
#' NRecords = c(500, 800)
#' Nlinks = 300
#' PmistakesA = c(0.02, 0.02, 0.02, 0.02, 0.02)
#' PmistakesB = c(0.02, 0.02, 0.02, 0.02, 0.02)
#' PmissingA = c(0.007, 0.007, 0.007, 0.007, 0.007)
#' PmissingB = c(0.007, 0.007, 0.007, 0.007, 0.007)
#' moving_params = list(V1=c(),V2=c(),V3=c(),V4=c(),V5=c(0.28))
#' enforceEstimability = TRUE
#' DATA = DataCreation( PIVs_config,
#'                      Nval,
#'                      NRecords,
#'                      Nlinks,
#'                      PmistakesA,
#'                      PmistakesB,
#'                      PmissingA,
#'                      PmissingB,
#'                      moving_params,
#'                      enforceEstimability)
#' A                    = DATA$A
#' B                    = DATA$B
#' Nvalues              = DATA$Nvalues
#' TimeDifference       = DATA$TimeDifference
#' proba_same_H         = DATA$proba_same_H
#'
#' TrueDelta = base::data.frame( matrix(0, nrow=0, ncol=2) )
#' for (i in 1:Nlinks)
#' {
#'   TrueDelta = rbind(TrueDelta, cbind(rownames(A[i,]),rownames(B[i,])))
#' }
#' true_pairs = do.call(paste, c(TrueDelta, list(sep="_")))
#'
#' encodedA = A
#' encodedB = B
#'
#' encodedA[,PIVs][ is.na(encodedA[,PIVs]) ] = 0
#' encodedB[,PIVs][ is.na(encodedB[,PIVs]) ] = 0
#'
#' DeltaResult = launchNaive(PIVs, encodedA, encodedB)
#'
#' results = base::data.frame( Results=matrix(NA, nrow=6, ncol=1) )
#' rownames(results) = c("tp","fp","fn","f1score","fdr","sens.")
#' if(nrow(DeltaResult)>1){
#'   linked_pairs    = do.call(paste, c(DeltaResult[,c("idxA","idxB")], list(sep="_")))
#'   truepositive    = length( intersect(linked_pairs, true_pairs) )
#'   falsepositive   = length( setdiff(linked_pairs, true_pairs) )
#'   falsenegative   = length( setdiff(true_pairs, linked_pairs) )
#'   precision       = truepositive / (truepositive + falsepositive)
#'   fdr             = 1 - precision
#'   sensitivity     = truepositive / (truepositive + falsenegative)
#'   f1score         = 2 * (precision * sensitivity) / (precision + sensitivity)
#'   results[,"Naive"] = c(truepositive,falsepositive,falsenegative,f1score,fdr,sensitivity)
#' }
#' results
launchNaive <- function(PIVs, encodedA, encodedB){
  testit::assert( "NA in A should be encoded as 0.", sum(is.na(encodedA))==0 )
  testit::assert( "NA in B should be encoded as 0.", sum(is.na(encodedB))==0 )
  rownames(encodedA) = 1:nrow(encodedA)
  rownames(encodedB) = 1:nrow(encodedB)
  DeltaNaiveLinked = data.frame( base::matrix(0, nrow=0, ncol=2) )
  colnames(DeltaNaiveLinked) = c("idxA", "idxB")
  # A NOT MISSING THAT MATCH WITH B
  isNotMissingA = apply(encodedA[,PIVs]!=0, 1, all)
  A_PIVs_notMissing = encodedA[isNotMissingA,PIVs]
  A_PIVs_notMissing_ID = rownames(encodedA[isNotMissingA,])
  UA = sspaste2(as.matrix(A_PIVs_notMissing))
  UB = sspaste2(as.matrix(encodedB[,PIVs]))
  valuesU = unique(c(UA,UB))
  UA = as.numeric(factor(UA,levels=valuesU))
  UB = as.numeric(factor(UB,levels=valuesU))
  tmpA = F2(UA, length(valuesU))
  tmpB = F2(UB, length(valuesU))
  select = F33(tmpA, tmpB, length(tmpA))
  if(nrow(select>0)){
    for (l in 1:nrow(select))
    {
      idxA = as.integer(select[l,1])
      idxA = A_PIVs_notMissing_ID[idxA]
      idxB = as.integer(select[l,2])
      idxB = rownames(encodedB[idxB,])
      DeltaNaiveLinked = rbind(DeltaNaiveLinked, cbind(idxA=idxA, idxB=idxB))
    }
  }
  # A MISSING THAT MATCH WITH B
  for(k in 1:length(PIVs)){
    isMissingA_k = encodedA[,k]==0
    A_PIVs_k_Missing = encodedA[isMissingA_k,PIVs]
    A_PIVs_k_Missing_ID = rownames(encodedA[isMissingA_k,])
    UA = sspaste2(as.matrix(A_PIVs_k_Missing[,-k]))
    UB = sspaste2(as.matrix(encodedB[,PIVs][,-k]))
    valuesU = unique(c(UA,UB))
    UA = as.numeric(factor(UA,levels=valuesU))
    UB = as.numeric(factor(UB,levels=valuesU))
    tmpA = F2(UA, length(valuesU))
    tmpB = F2(UB, length(valuesU))
    select = F33(tmpA, tmpB, length(tmpA))
    if(nrow(select>0)){
      for (l in 1:nrow(select))
      {
        idxA = as.integer(select[l,1])
        idxA = A_PIVs_k_Missing_ID[idxA]
        idxB = as.integer(select[l,2])
        idxB = rownames(encodedB[idxB,])
        DeltaNaiveLinked = rbind(DeltaNaiveLinked, cbind(idxA=idxA, idxB=idxB))
      }
    }
  }
  # B NOT MISSING THAT MATCH WITH A
  isNotMissingB = apply(encodedB[,PIVs]!=0, 1, all)
  B_PIVs_notMissing = encodedB[isNotMissingB,PIVs]
  B_PIVs_notMissing_ID = rownames(encodedB[isNotMissingB,])
  UB = sspaste2(as.matrix(B_PIVs_notMissing))
  UA = sspaste2(as.matrix(encodedA[,PIVs]))
  valuesU = unique(c(UA,UB))
  UA = as.numeric(factor(UA,levels=valuesU))
  UB = as.numeric(factor(UB,levels=valuesU))
  tmpA = F2(UA, length(valuesU))
  tmpB = F2(UB, length(valuesU))
  select = F33(tmpA, tmpB, length(tmpA))
  if(nrow(select>0)){
    for (l in 1:nrow(select))
    {
      idxA = as.integer(select[l,1])
      idxA = rownames(encodedA[idxA,])
      idxB = as.integer(select[l,2])
      idxB = B_PIVs_notMissing_ID[idxB]
      DeltaNaiveLinked = rbind(DeltaNaiveLinked, cbind(idxA=idxA, idxB=idxB))
    }
  }
  # B MISSING THAT MATCH WITH A
  for(k in 1:length(PIVs)){
    isMissingB_k = encodedB[,k]==0
    B_PIVs_k_Missing = encodedB[isMissingB_k,PIVs]
    B_PIVs_k_Missing_ID = rownames(encodedB[isMissingB_k,])
    UB = sspaste2(as.matrix(B_PIVs_k_Missing[,-k]))
    UA = sspaste2(as.matrix(encodedA[,PIVs][,-k]))
    valuesU = unique(c(UA,UB))
    UA = as.numeric(factor(UA,levels=valuesU))
    UB = as.numeric(factor(UB,levels=valuesU))
    tmpA = F2(UA, length(valuesU))
    tmpB = F2(UB, length(valuesU))
    select = F33(tmpA, tmpB, length(tmpA))
    if(nrow(select>0)){
      for (l in 1:nrow(select))
      {
        idxA = as.integer(select[l,1])
        idxA = rownames(encodedA[idxA,])
        idxB = as.integer(select[l,2])
        idxB = B_PIVs_k_Missing_ID[idxB]
        DeltaNaiveLinked = rbind(DeltaNaiveLinked, cbind(idxA=idxA, idxB=idxB))
      }
    }
  }
  DeltaNaiveLinked = DeltaNaiveLinked[!duplicated(DeltaNaiveLinked), ]
  DeltaNaiveLinked
}
