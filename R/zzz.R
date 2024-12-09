#' FlexRL
#'
#' A Flexible Model For Record Linkage
#'
#' The example below aims to link 2 synthetic data sources, with 5 PIVs (4 stable ones and 1 unstable, which evolve over time).
#' PIVs may be considered stable if there is not enough information to model their dynamics. Since we synthesise the data we can model it.
#' We may not give a bound on some PIVs in real settings, since there may be a lot of disagreements among the links for those variables (in situations where we would have liked to model their dynamics otherwise but we do not have enough information for this).
#' Here we bound all the mistakes parameters since we know that the mistakes probabilities are inferior to 10%.
#' For a small example we prefer having one parameter for the probabilities of mistakes over the 2 data sources.
#' We do need to fix the mistake parameter of the 5th PIV to avoid estimability problems here (since it is an 'unstable' variable).
#' We know the true linkage structure in this example so we can compute performances of the method at the end.
#'
#' There are more details to understand the method in our paper, or on the experiments repository of our paper, or in the vignettes, or in the documentation of the main algorithm ?FlexRL::StEM.
#'
#' @author Kayané Robach
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib FlexRL, .registration=TRUE
#' @name FlexRL
#'
#' @param data is a list gathering information on the data to be linked
#' - sources 'A' and 'B',
#' - a vector 'Nvalues' gathering the number of unique values in each PIV,
#' - 'PIVs_config' the list of PIVs to use for Record Linkage and details on how each should be handled by the algorithm,
#' - potential bounds on the mistakes probabilities for each PIV: 'controlOnMistakes',
#' - 'sameMistakes' whether there should be one parameter for the mistakes in A and B or whether each source should have its own (in case of small data sources it is recommended to set sameMistakes=TRUE)
#' - whether the parameters for mistakes should be fixed in case of instability 'phiMistakesAFixed' and 'phiMistakesBFixed',
#' - as well as the values they should be fixed to 'phiForMistakesA' and 'phiForMistakesB'
#' @param StEMIter The total number of iterations of the Stochastic Expectation Maximisation (StEM) algorithm (including the period to discard as burn-in)
#' @param StEMBurnin The number of iterations to discard as burn-in
#' @param GibbsIter The total number of iterations of the Gibbs sampler
#' (run in each iteration of the StEM) (including the period to discard as burn-in)
#' @param GibbsBurnin The number of iterations to discard as burn-in
#'
#' @return The Stochastic Expectation Maximisation (StEM) function returns w list with:
#' - Delta, the (sparse) matrix with the pairs of records linked and their posterior probabilities to be linked (select the pairs where the proba>0.5 to get a valid set of linked records),
#' - as well as the model parameters chains:
#'    - gamma,
#'    - eta,
#'    - alpha,
#'    - phi.
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
#' TimeDifference       = DATA$TimeDifference
#' proba_same_H         = DATA$proba_same_H
#'
#' # the first 1:Nlinks records of each files created are links
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
#' data = list( A                    = encodedA,
#'              B                    = encodedB,
#'              Nvalues              = Nvalues,
#'              PIVs_config          = PIVs_config,
#'              controlOnMistakes    = c(TRUE, TRUE, TRUE, TRUE, TRUE),
#'              sameMistakes         = TRUE,
#'              phiMistakesAFixed    = TRUE,
#'              phiMistakesBFixed    = TRUE,
#'              phiForMistakesA      = c(NA, NA, NA, NA, 0),
#'              phiForMistakesB      = c(NA, NA, NA, NA, 0)
#'            )
#'
#' fit = FlexRL::stEM(  data                 = data,
#'                      StEMIter             = 50,
#'                      StEMBurnin           = 30,
#'                      GibbsIter            = 50,
#'                      GibbsBurnin          = 30,
#'                      musicOn              = TRUE,
#'                      newDirectory         = NULL,
#'                      saveInfoIter         = FALSE
#'                   )
#'
#' DeltaResult = fit$Delta
#' colnames(DeltaResult) = c("idxA","idxB","probaLink")
#' DeltaResult = DeltaResult[DeltaResult$probaLink>0.5,]
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
#'   results[,"FlexRL"] = c(truepositive,falsepositive,falsenegative,f1score,fdr,sensitivity)
#' }
#' }
#'
NULL

.onLoad <- function(...) {
  base::packageStartupMessage("If you are happy with FlexRL, please cite us! Also, if you are unhappy, please cite us anyway.\nHERE ADD\nbibtex format in CITATION.", appendLF = TRUE)
}
