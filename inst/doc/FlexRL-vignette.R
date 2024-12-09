## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(FlexRL)

## -----------------------------------------------------------------------------
df2016 = read.csv("exA.csv", row.names = 1)
df2020 = read.csv("exB.csv", row.names = 1)

## -----------------------------------------------------------------------------
head(df2016)

## -----------------------------------------------------------------------------
head(df2020)

## -----------------------------------------------------------------------------
PIVs_config = list( 
  ANASCI     = list(stable = TRUE),
  SESSO      = list(stable = TRUE),
  STACIV     = list(stable = TRUE),
  STUDIO     = list(stable = TRUE), 
  IREG       = list(stable = TRUE) 
  )
PIVs = names(PIVs_config)
PIVs_stable = sapply(PIVs_config, function(x) x$stable)

## -----------------------------------------------------------------------------
for(i in 1:length(PIVs)){
  intersect_support_piv = intersect( unique(df2016[,PIVs[i]]), unique(df2020[,PIVs[i]]) )
  df2016 = df2016[df2016[,PIVs[i]] %in% c(NA,intersect_support_piv),]
  df2020 = df2020[df2020[,PIVs[i]] %in% c(NA,intersect_support_piv),]
}

rownames(df2016) = 1:nrow(df2016)
rownames(df2020) = 1:nrow(df2020)

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
nrow(df2016)

## -----------------------------------------------------------------------------
nrow(df2020)

## -----------------------------------------------------------------------------
Nlinks

## -----------------------------------------------------------------------------
colSums(df2016[TrueDelta[,1],PIVs] == df2020[TrueDelta[,2],PIVs]) / Nlinks

## -----------------------------------------------------------------------------
colSums(is.na(df2016[TrueDelta[,1],PIVs])) / Nlinks

## -----------------------------------------------------------------------------
colSums(is.na(df2020[TrueDelta[,2],PIVs])) / Nlinks

## -----------------------------------------------------------------------------
df2016$source = "df2016"
df2020$source = "df2020"

if(nrow(df2020)>nrow(df2016)){
  encodedA = df2016
  encodedB = df2020
  cat("df2020 is the largest file, denoted encodedB")
}else{
  encodedB = df2016
  encodedA = df2020
  cat("df2016 is the largest file, denoted encodedB")
}

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

## -----------------------------------------------------------------------------
nvalues

## -----------------------------------------------------------------------------
data = list( A                    = encodedA,
             B                    = encodedB, 
             Nvalues              = nvalues,
             PIVs_config          = PIVs_config,
             controlOnMistakes    = c(TRUE, TRUE, FALSE, FALSE, FALSE),
             sameMistakes         = TRUE,
             phiMistakesAFixed    = FALSE,
             phiMistakesBFixed    = FALSE,
             phiForMistakesA      = c(NA, NA, NA, NA, NA),
             phiForMistakesB      = c(NA, NA, NA, NA, NA)
             )

fit = FlexRL::stEM(  data                 = data,
             StEMIter             = 50,
             StEMBurnin           = 30,
             GibbsIter            = 100, 
             GibbsBurnin          = 70,
             musicOn              = TRUE,
             newDirectory         = NULL,
             saveInfoIter         = FALSE
             )

## -----------------------------------------------------------------------------
DeltaResult = fit$Delta
colnames(DeltaResult) = c("idx20","idx16","probaLink")
DeltaResult = DeltaResult[DeltaResult$probaLink>0.5,]
DeltaResult

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
colSums(encodedA[DeltaResult[,"idx20"],PIVs] == encodedB[DeltaResult[,"idx16"],PIVs]) / nrow(DeltaResult)

