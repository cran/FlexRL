
# FlexRL

<!-- badges: start -->
<!-- badges: end -->

FlexRL is a package for Flexible Record Linkage, to find the common set of records among 2 data sources. The area of applications of Record Linkage is broad, it can be used to link data from any sources where an overlap in the populations is expected, like healthcare monitoring studies at 2 different time points, registries of casualties in conflict zones collected by distinct organisations, ... 

FlexRL models registration errors (missing values and mistakes in the data) and handles dynamic **P**artially **I**dentifying **V**ariable**s** that evolve over time (e.g. postal code can change between the 2 data collections) in order to identify a final set of linked records (and their posterior probabilities to be linked). 

The algorithm can take time to run on large data sets but has a low memory footprint and can easily run on standard computers.

This package implements the **St**ochastic **EM** approach to Record Linkage described in '[A flexible model for Record Linkage](https://arxiv.org/abs/2407.06835)'. The main article and the supplementary material are available on arxiv.

Please [open an issue](https://github.com/robachowyk/FlexRL/issues) to report any bug, to make a request, or to ask for help :-)

## Installation

You can install FlexRL from CRAN with:
``` r
install.packages("FlexRL")
library(FlexRL)
```

Or you can install the development version of FlexRL from its [GitHub](https://github.com/robachowyk/FlexRL) with one of the following:

``` r
pak::pak("robachowyk/FlexRL")
```
``` r
remotes::install_github("robachowyk/FlexRL")
```
``` r
devtools::install_github("robachowyk/FlexRL")
```

FlexRL relies on Rcpp; when imported from Github, it may require gfortran and gcc.

## How to use `FlexRL`

Here is a basic example which shows how to solve a common record linkage task:

``` r
library(FlexRL)

# load real data subsets from the vignettes
df2016 = read.csv("FlexRL/vignettes/exA.csv", row.names = 1)
df2020 = read.csv("FlexRL/vignettes/exB.csv", row.names = 1)

# use 5 PIVs birth year, sex, marital status, educational level, regional code; 
# we do not have enough information to model instability
# all PIVs are considered stable
PIVs_config = list(
  ANASCI     = list(stable = TRUE),
  SESSO      = list(stable = TRUE),
  STACIV     = list(stable = TRUE),
  STUDIO     = list(stable = TRUE),
  IREG       = list(stable = TRUE)
)
PIVs = names(PIVs_config)
PIVs_stable = sapply(PIVs_config, function(x) x$stable)

# we put bounds on the probability of mistakes: 
# there should realistically be less than 10% of mistakes in the PIVs
# however for some PIVs which are probably dynamic (though we cannot model it here) it is good to
# not bound the mistakes parameter, which will adapt to the dynamics not taken into account
boundMistakes = c(TRUE, TRUE, FALSE, FALSE, FALSE)

# filter the data to their common support
for(i in 1:length(PIVs)){
  intersect_support_piv = intersect( unique(df2016[,PIVs[i]]), unique(df2020[,PIVs[i]]) )
  df2016 = df2016[df2016[,PIVs[i]] %in% c(NA,intersect_support_piv),]
  df2020 = df2020[df2020[,PIVs[i]] %in% c(NA,intersect_support_piv),]
}

# reset the rownames
rownames(df2016) = 1:nrow(df2016)
rownames(df2020) = 1:nrow(df2020)

# we know the true linkage structure so we know the true links
links = intersect(df2016$ID, df2020$ID)
Nlinks = length(links)

TrueDelta = data.frame( matrix(0, nrow=0, ncol=2) )
for (i in 1:Nlinks)
{
  id = links[i]
  id16 = which(df2016$ID == id)
  id20 = which(df2020$ID == id)
  TrueDelta = rbind(TrueDelta, cbind(rownames(df2016[id16,]),rownames(df2020[id20,])))
}
true_pairs = do.call(paste, c(TrueDelta, list(sep="_")))

# we need a source column
df2016$source = "df2016"
df2020$source = "df2020"

# the first dataset (namely source A) has to be the smallest one 
# i.e. the second dataset (namely source B) has to be the largest one
if(nrow(df2020)>nrow(df2016)){
  encodedA = df2016
  encodedB = df2020
  cat("df2020 is the largest file, denoted encodedB")
}else{
  encodedB = df2016
  encodedA = df2020
  cat("df2016 is the largest file, denoted encodedB")
}

# encode the PIVs
levels_PIVs = lapply(PIVs, function(x) levels(factor(as.character(c(encodedA[,x], encodedB[,x])))))

for(i in 1:length(PIVs))
{
  encodedA[,PIVs[i]] = as.numeric(factor(as.character(encodedA[,PIVs[i]]), levels=levels_PIVs[[i]]))
  encodedB[,PIVs[i]] = as.numeric(factor(as.character(encodedB[,PIVs[i]]), levels=levels_PIVs[[i]]))
}
nvalues = sapply(levels_PIVs, length)
names(nvalues) = PIVs

encodedA[,PIVs][ is.na(encodedA[,PIVs]) ] = 0
encodedB[,PIVs][ is.na(encodedB[,PIVs]) ] = 0

data = list( A                    = encodedA,
             B                    = encodedB,
             Nvalues              = nvalues,
             PIVs_config          = PIVs_config,
             controlOnMistakes    = boundMistakes,
             sameMistakes         = TRUE,
             phiMistakesAFixed    = FALSE,
             phiMistakesBFixed    = FALSE,
             phiForMistakesA      = c(NA, NA, NA, NA, NA),
             phiForMistakesB      = c(NA, NA, NA, NA, NA)
)

# launch FlexRL algorithm
fit = stEM(  data                 = data,
             StEMIter             = 50,
             StEMBurnin           = 30,
             GibbsIter            = 50,
             GibbsBurnin          = 30,
             musicOn              = TRUE,
             newDirectory         = NULL,
             saveInfoIter         = FALSE
)

# collect the output and build the final set of linked records with posterior probability > 0.5 
# to ensure that the one-to-one assignment constraint is fulfilled
DeltaResult = fit$Delta
colnames(DeltaResult) = c("idx20","idx16","probaLink")
DeltaResult = DeltaResult[DeltaResult$probaLink>0.5,]
DeltaResult

# compute the results by comparing the true links to the linked records
# watch out for the order of the source in the outcome, source A here corresponds to df2020 
# and source B corresponds to df2016
# in the set of true links we ordered pairs with first df2016 and then df2020
results = data.frame( Results=matrix(NA, nrow=6, ncol=1) )
rownames(results) = c("tp","fp","fn","f1score","fdr","sens.")
if(nrow(DeltaResult)>1){
  linked_pairs    = do.call(paste, c(DeltaResult[,c("idx16","idx20")], list(sep="_")))
  truepositive    = length( intersect(linked_pairs, true_pairs) )
  falsepositive   = length( setdiff(linked_pairs, true_pairs) )
  falsenegative   = length( setdiff(true_pairs, linked_pairs) )
  precision       = truepositive / (truepositive + falsepositive)
  fdr             = 1 - precision
  sensitivity     = truepositive / (truepositive + falsenegative)
  f1score         = 2 * (precision * sensitivity) / (precision + sensitivity)
  results[,"FlexRL"] = c(truepositive,falsepositive,falsenegative,f1score,fdr,sensitivity)
}
results

# one can use the simplistic method to assess the difficulty of the task 
DeltaResult = launchNaive(PIVs, encodedA, encodedB)

if(nrow(DeltaResult)>1){
  linked_pairs    = do.call(paste, c(DeltaResult[,c("idxB","idxA)], list(sep="_")))
  truepositive    = length( intersect(linked_pairs, true_pairs) )
  falsepositive   = length( setdiff(linked_pairs, true_pairs) )
  falsenegative   = length( setdiff(true_pairs, linked_pairs) )
  precision       = truepositive / (truepositive + falsepositive)
  fdr             = 1 - precision
  sensitivity     = truepositive / (truepositive + falsenegative)
  f1score         = 2 * (precision * sensitivity) / (precision + sensitivity)
  results[,"Naive"] = c(truepositive,falsepositive,falsenegative,f1score,fdr,sensitivity)
}
results
```

Here is another example (using synthetic data to illustrate more elaborated options) which shows you how to solve a common record linkage problem:


``` r
library(FlexRL)

PIVs_config = list( V1 = list(stable = TRUE),
                    V2 = list(stable = TRUE),
                    V3 = list(stable = TRUE),
                    V4 = list(stable = TRUE),
                    V5 = list( stable = FALSE,
                               conditionalHazard = FALSE,
                               pSameH.cov.A = c(),
                               pSameH.cov.B = c()) )
PIVs = names(PIVs_config)
PIVs_stable = sapply(PIVs_config, function(x) x$stable)
Nval = c(6, 7, 8, 9, 15)
NRecords = c(500, 800)
Nlinks = 300
PmistakesA = c(0.02, 0.02, 0.02, 0.02, 0.02)
PmistakesB = c(0.02, 0.02, 0.02, 0.02, 0.02)
PmissingA = c(0.007, 0.007, 0.007, 0.007, 0.007)
PmissingB = c(0.007, 0.007, 0.007, 0.007, 0.007)
moving_params = list(V1=c(),V2=c(),V3=c(),V4=c(),V5=c(0.28))
enforceEstimability = TRUE
DATA = DataCreation( PIVs_config,
                     Nval,
                     NRecords,
                     Nlinks,
                     PmistakesA,
                     PmistakesB,
                     PmissingA,
                     PmissingB,
                     moving_params,
                     enforceEstimability)
A                    = DATA$A
B                    = DATA$B
Nvalues              = DATA$Nvalues
TimeDifference       = DATA$TimeDifference
proba_same_H         = DATA$proba_same_H

# we generate data with an unstable PIV (representing that people move for instance)
proba_same_H_5 = proba_same_H[,5]
plot( sort(proba_same_H_5, decreasing=TRUE), ylim=c(0,1) )
plot( TimeDifference, proba_same_H_5, ylim=c(0,1) )

# it is realistic to bound the parameter for mistakes to not exceeds 10% 
# (also for the unstable PIV since we model its dynamics)
boundMistakes = c(TRUE, TRUE, TRUE, TRUE, TRUE)

# the first 1:Nlinks records of each files created are links
TrueDelta = data.frame( matrix(0, nrow=0, ncol=2) )
for (i in 1:Nlinks)
{
  TrueDelta = rbind(TrueDelta, cbind(rownames(A[i,]),rownames(B[i,])))
}
true_pairs = do.call(paste, c(TrueDelta, list(sep="_")))

encodedA = A
encodedB = B

encodedA[,PIVs][ is.na(encodedA[,PIVs]) ] = 0
encodedB[,PIVs][ is.na(encodedB[,PIVs]) ] = 0

data = list( A                    = encodedA,
             B                    = encodedB,
             Nvalues              = Nvalues,
             PIVs_config          = PIVs_config,
             controlOnMistakes    = boundMistakes,
             sameMistakes         = TRUE,
             phiMistakesAFixed    = TRUE,
             phiMistakesBFixed    = TRUE,
             phiForMistakesA      = c(NA, NA, NA, NA, 0),
             phiForMistakesB      = c(NA, NA, NA, NA, 0)
           )

fit = stEM(  data                 = data,
             StEMIter             = 100,
             StEMBurnin           = 70,
             GibbsIter            = 200,
             GibbsBurnin          = 100,
             musicOn              = TRUE,
             newDirectory         = NULL,
             saveInfoIter         = FALSE
          )

DeltaResult = fit$Delta
colnames(DeltaResult) = c("idxA","idxB","probaLink")
DeltaResult = DeltaResult[DeltaResult$probaLink>0.5,]

results = data.frame( Results=matrix(NA, nrow=6, ncol=1) )
rownames(results) = c("tp","fp","fn","f1score","fdr","sens.")
if(nrow(DeltaResult)>1){
  linked_pairs    = do.call(paste, c(DeltaResult[,c("idxA","idxB")], list(sep="_")))
  truepositive    = length( intersect(linked_pairs, true_pairs) )
  falsepositive   = length( setdiff(linked_pairs, true_pairs) )
  falsenegative   = length( setdiff(true_pairs, linked_pairs) )
  precision       = truepositive / (truepositive + falsepositive)
  fdr             = 1 - precision
  sensitivity     = truepositive / (truepositive + falsenegative)
  f1score         = 2 * (precision * sensitivity) / (precision + sensitivity)
  results[,"FlexRL"] = c(truepositive,falsepositive,falsenegative,f1score,fdr,sensitivity)
}
results
```

More documentation is accessible on CRAN.

More examples are available in the vignette and in the repository [FlexRL-experiments](https://github.com/robachowyk/FlexRL-experiments).
 
For support requests, contact _k dot c dot robach at amsterdamumc dot nl_.
