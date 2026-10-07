
# ============================================================================
# FlexRL — flexible probabilistic record linkage with data quality diagnostics 
# for downstream inference on linked data.
# ============================================================================


#' Simulate two linked data sources for record linkage benchmarking
#'
#' Creates two synthetic data sources of given sizes sharing a given number
#' of common entities ("links"), each described by a set of Partially
#' Identifying Variables (PIVs). For every PIV, choose the number of
#' possible values, the proportion of mistakes and of missing values, and
#' whether it is stable over time, flexible (may change but change is not
#' modelled), or structured (expected to change over time, with a survival
#' model for the hazard of change). For structured PIVs, `enforce_estimability`
#' forces half of the linked pairs to have a near-zero time gap, which helps
#' separate "mistake" from "change over time" when fitting the model.
#'
#' @param PIVs_config Named list, one entry per PIV, each a list with:
#'   `dynamics` (`"stable"`, `"flexible"`, or `"structured"`),
#'   `bound_mistakes` (length-2 numeric/NA, upper bound on the mistake
#'   probability in file 1 / file 2), `fix_mistakes` (length-2 numeric/NA,
#'   mistake probability fixed to this value in file 1 / file 2), and, only
#'   for `dynamics = "structured"`, `cond_hazard_cov` (a list with `cov1` and
#'   `cov2`, the names of covariates in file 1 / file 2 used to model the
#'   hazard of change).
#' @param n_values Integer vector, number of unique values per PIV (same order
#'   as `PIVs_config`).
#' @param n_records Integer vector of length 2, number of records to generate
#'   in file 1 and file 2 (file 2 must be the larger of the two).
#' @param n_links Integer, number of records shared between the two files.
#' @param p_mistake Named list (one entry per PIV) of length-2 numeric vectors,
#'   proportion of mistakes to introduce in file 1 / file 2.
#' @param p_missing Named list (one entry per PIV) of length-2 numeric vectors,
#'   proportion of missing values to introduce in file 1 / file 2.
#' @param cond_hazard_params Named list (one entry per PIV) of numeric vectors
#'   for the survival model generating the changes (see [survival_model()]); 
#'   only used for `dynamics = "structured"` PIVs.
#' @param enforce_estimability Logical; if `TRUE`, half of the linked pairs are
#'   given a near-zero time gap to help estimate the instability parameters.
#' @param model_dynamics Object from [survival_model()] used to generate the
#'   changes of the structured PIVs (default exponential).
#'
#' @return A list with: `data1`, `data2` the two simulated (encoded)
#'   data frames, `n_values` number of unique values per PIV, `time_difference`
#'   time gap between linked records (`NA` if no PIV is structured),
#'   `proba_same_H` matrix (n links x n PIVs) of probabilities that the true
#'   values coincide, `true_pairs` data frame with the true `1`/`2` indices
#'   of the linked records
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V2 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V3 = list(dynamics = "flexible",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(NA,NA)),
#'                      V4 = list(dynamics = "structured",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(0.03,0.03),
#'                                cond_hazard_cov = list(cov1=c("Xe", "Xf"),
#'                                                     cov2=c())) )
#' n_values  <- c( 5, 6, 7, 12 )
#' p_mistake <- list( V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                   V3 = c(0.05, 0.05), V4 = c(0.02, 0.02) )
#' p_missing <- list( V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                   V3 = c(0.005, 0.005), V4 = c(0.005, 0.005) )
#' cond_hazard_params <- list( V1 = c(), V2 = c(), 
#'                             V3 = c(), V4 = log(c(0.7, 0.6, 0.5)) )
#' gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
#'                            p_mistake, p_missing, cond_hazard_params, TRUE )
#' str(gen_data, max.level = 1)
simulate_data <- function(PIVs_config, n_values, n_records, n_links, p_mistake, p_missing,
                          cond_hazard_params, enforce_estimability,
                          model_dynamics = survival_model("exponential")) {

  if (n_records[2] < n_records[1]) {
    stop("`n_records[2]` (file B) must be >= `n_records[1]` (file A).", call. = FALSE)
  }
  if (n_links > min(n_records)) {
    stop("`n_links` cannot exceed the size of the smallest file.", call. = FALSE)
  }

  PIVs <- names(PIVs_config)
  PIVs_stable <- sapply(PIVs_config, function(x) x$dynamics != "structured")
  PIVs_conditionalHazardVariables <- sapply(PIVs_config, function(x) {
    if (x$dynamics == "structured") x$cond_hazard_cov else FALSE
  })
  modelDynaPIVs <- which(!PIVs_stable)

  Pmistake1 <- sapply(p_mistake, function(x) x[1])
  Pmistake2 <- sapply(p_mistake, function(x) x[2])
  Pmissing1 <- sapply(p_missing, function(x) x[1])
  Pmissing2 <- sapply(p_missing, function(x) x[2])

  # Simulate the true PIV values for both files
  draw_pivs <- function(n) {
    out <- lapply(n_values, function(nv) {
      xp <- exp(0.23 * (0:(nv - 1)))
      sample(seq_len(nv), n, replace = TRUE, prob = xp / sum(xp))
    })
    stats::setNames(as.data.frame(out), PIVs)
  }
  data1 <- draw_pivs(n_records[1])
  data2 <- draw_pivs(n_records[2])

  # First `n_links` records of file 2 are copies of the first `n_links` of file 1
  data2[1:n_links, ] <- data1[1:n_links, ]

  # Introduce mistakes
  for (x in seq_len(ncol(data1))) {
    biased <- as.logical(stats::rbinom(nrow(data1), 1, Pmistake1[x]))
    if (any(biased)) {
      data1[, x][biased] <- sapply(data1[, x][biased], function(i) sample((1:n_values[x])[-c(i)], 1))
    }
  }
  for (x in seq_len(ncol(data2))) {
    biased <- as.logical(stats::rbinom(nrow(data2), 1, Pmistake2[x]))
    if (any(biased)) {
      data2[, x][biased] <- sapply(data2[, x][biased], function(i) sample((1:n_values[x])[-c(i)], 1))
    }
  }

  # Introduce missing values
  for (x in seq_len(ncol(data1))) {
    biased <- as.logical(stats::rbinom(nrow(data1), 1, Pmissing1[x]))
    if (any(biased)) data1[, x][biased] <- NA
  }
  for (x in seq_len(ncol(data2))) {
    biased <- as.logical(stats::rbinom(nrow(data2), 1, Pmissing2[x]))
    if (any(biased)) data2[, x][biased] <- NA
  }

  data1$change <- FALSE
  data2$change <- FALSE

  if (any(!PIVs_stable)) {

    # Registration times (file 2 always registered after file 1)
    data1$date <- stats::runif(nrow(data1), 0, 3)
    data2$date <- stats::runif(nrow(data2), 3, 6)

    if (enforce_estimability) {
      null_time_diff <- as.integer(n_links / 2)
      data1[1:null_time_diff, "date"] <- stats::runif(null_time_diff, 0.00, 0.01)
      data2[1:null_time_diff, "date"] <- stats::runif(null_time_diff, 0.00, 0.01)
    }

    time_difference <- abs(data2[1:n_links, "date"] - data1[1:n_links, "date"])
    intercept <- rep(1, n_links)
    proba_same_H <- matrix(1, n_links, length(PIVs))

    for (k in modelDynaPIVs) {

      hasCov <- !is.logical(PIVs_conditionalHazardVariables[[k]]) &&
        any(!sapply(PIVs_conditionalHazardVariables[[k]], is.null))

      if (hasCov) {
        for (covariate in PIVs_config[[k]]$cond_hazard_cov$cov1) {
          data1[, covariate] <- stats::rnorm(nrow(data1), 1, 1)
        }
        for (covariate in PIVs_config[[k]]$cond_hazard_cov$cov2) {
          data2[, covariate] <- stats::rnorm(nrow(data2), 2, 1)
        }
        cov <- cbind(
          intercept,
          data1[1:n_links, PIVs_config[[k]]$cond_hazard_cov$cov1, drop = FALSE],
          data2[1:n_links, PIVs_config[[k]]$cond_hazard_cov$cov2, drop = FALSE]
        )
      } else {
        cov <- cbind(intercept)
      }

      n_par_k <- model_dynamics$n_par(ncol(as.matrix(cov)))
      if (length(cond_hazard_params[[k]]) != n_par_k) {
        stop(sprintf(
          "`cond_hazard_params[['%s']]` must have length %d for the %s survival model (see ?survival_model).",
          PIVs[k], n_par_k, model_dynamics$type
        ), call. = FALSE)
      }

      proba_same_H[, k] <- model_dynamics$S(as.matrix(cov), cond_hazard_params[[k]], time_difference)

      # Generate instability: for each linked pair, decide (via the survival
      # probability above) whether the true value changed between the two
      # registrations
      for (i in seq_len(n_links)) {
        is_not_changing <- stats::rbinom(1, 1, proba_same_H[i, k])
        if (!is_not_changing && !is.na(data1[i, x])) {
          data2[i, x] <- sample((1:n_values[x])[-c(data1[i, x])], 1)
          data2[i, "change"] <- TRUE
        }
      }
    }
  } else {
    time_difference <- NA
    proba_same_H <- NA
  }

  # Recode the PIVs
  levels_PIVs <- lapply(PIVs, function(x) levels(factor(as.character(c(data1[, x], data2[, x])))))
  for (i in seq_along(PIVs)) {
    data1[, PIVs[i]] <- as.numeric(factor(as.character(data1[, PIVs[i]]), levels = levels_PIVs[[i]]))
    data2[, PIVs[i]] <- as.numeric(factor(as.character(data2[, PIVs[i]]), levels = levels_PIVs[[i]]))
  }
  n_values_recoded <- sapply(levels_PIVs, length)

  data1$local_id <- seq_len(nrow(data1))
  data2$local_id <- seq_len(nrow(data2))
  if(nrow(data1) > n_links){
    data1$entity_id <- c(seq_len(n_links), seq.int(n_links + 1, nrow(data1)))
  }else{
    data1$entity_id <- seq_len(n_links)
  }
  if(nrow(data2) > n_links){
    data2$entity_id <- c(seq_len(n_links), seq.int(nrow(data1) + 1, nrow(data1) + nrow(data2) - n_links))
  }else{
    data2$entity_id <- seq_len(n_links)
  }
  data1$source <- "1"
  data2$source <- "2"

  true_Delta <- data.frame(`1` = seq_len(n_links), `2` = seq_len(n_links), check.names = FALSE)

  list(
    data1 = data1, data2 = data2, n_values = n_values_recoded,
    time_difference = time_difference, proba_same_H = proba_same_H, true_pairs = true_Delta
  )
}

#' Book-keeping data frame for parameters of PIVs dynamics
#'
#' Internal helper used by [StEM()] to accumulate, across Gibbs iterations,
#' the covariates, true-value agreement indicator, and time gaps needed to
#' re-estimate the survival (hazard) parameters of a dynamic structured PIV.
#'
#' @param n_coef_unstable Integer, number of hazard coefficients for this PIV
#'   (1 for the baseline hazard, plus one per covariate from file A and file B).
#' @param stable Logical, whether this PIV is stable
#'   (`dynamics != "structured"`).
#'
#' @return An empty data frame with `n_coef_unstable + 2` columns (the
#'   covariates/intercept, plus `Hequal` and `times`) if `stable` is `FALSE`;
#'   `NULL` if the PIV is stable (nothing to accumulate).
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V2 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V3 = list(dynamics = "flexible",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(NA,NA)),
#'                      V4 = list(dynamics = "structured",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(0.03,0.03),
#'                                cond_hazard_cov = list(cov1=c("Xe", "Xf"),
#'                                                     cov2=c())) )
#' PIVs_stable <- sapply(PIVs_config, function(x) x$dynamics != "structured")
#' n_coef_unstable = c(0,0,0,3)
#' Valpha <- mapply(create_data_alpha, n_coef_unstable = n_coef_unstable,
#'                  stable = PIVs_stable, SIMPLIFY = FALSE)
create_data_alpha <- function(n_coef_unstable, stable) {
  if (stable) return(invisible(NULL))
  n_coef <- n_coef_unstable + 2  # covariates/intercept + Hequal + times
  data.frame(matrix(nrow = 0, ncol = n_coef))
}

#' Log number of possible linkage configurations
#'
#' Computes `log(nB! / (nB - sumD)!)`, i.e. the log number of ways to choose
#' `sumD` ordered links among `n_records_B` records of the larger file; used in
#' [log_lik()] to normalise the likelihood of the linkage matrix.
#'
#' @param n_records_B Integer, number of records in the larger data source (B).
#' @param sumD Integer, number of currently linked records.
#'
#' @return Numeric, `sum(log((n_records_B - sumD + 1):n_records_B))`, or `0` if 
#'   `sumD == 0`.
#' @export
#'
#' @examples
#' log_possible_config(n_records_B = 15, sumD = 5)
log_possible_config <- function(n_records_B, sumD) {
  if (sumD > 0) sum(log(n_records_B:(n_records_B - sumD + 1))) else 0
}

#' Log-likelihood of the linkage matrix
#'
#' @param LLL (Sparse) matrix of log-likelihood contributions for linked
#'   records.
#' @param LLA Numeric vector of log-likelihood contributions for non-linked
#'   records from A.
#' @param LLB Numeric vector of log-likelihood contributions for non-linked
#'   records from B.
#' @param links 2-column matrix of indices (A, B) for the currently linked
#'   records.
#' @param sumRowD Logical vector, one entry per record in A: does it form a
#'   link?
#' @param sumColD Logical vector, one entry per record in B: does it form a
#'   link?
#' @param gamma Numeric, proportion of linked records as a fraction of the
#'   smaller file.
#'
#' @return Numeric, the log-likelihood of the linkage matrix.
#' @export
#'
#' @examples
#' LLL <- Matrix::Matrix(0, nrow = 13, ncol = 15, sparse = TRUE)
#' LLA <- stats::runif(13, 0, 2)
#' LLB <- stats::runif(15, 0, 2)
#' links <- as.matrix(data.frame(idxA = c(5, 9, 11, 12, 13),
#'                               idxB = c(5, 9, 11, 13, 15)))
#' LLL[links] <- 0.67
#' sumRowD <- (seq_len(13) %in% links[, 1])
#' sumColD <- (seq_len(15) %in% links[, 2])
#' gamma <- 0.5
#' log_lik(LLL, LLA, LLB, links, sumRowD, sumColD, gamma)
log_lik <- function(LLL, LLA, LLB, links, sumRowD, sumColD, gamma) {
  if (length(sumColD) - nrow(links) + 1 <= 0) {
    warning(
      "File B has fewer records (", length(sumColD), ") than linked records (",
      nrow(links), "); file A should never be larger than file B.",
      call. = FALSE
    )
  }
  logPossD <- sum(log(gamma) * sumRowD + log(1 - gamma) * (1 - sumRowD)) -
    log_possible_config(length(sumColD), nrow(links))
  logPossD + sum(LLA[sumRowD == 0]) + sum(LLB[sumColD == 0]) + sum(LLL[links])
}

#' Simulate the true PIV values underlying the registered records
#'
#' Draw the latent true values of each PIV given the currently registered 
#' (possibly mistaken or missing) values, the current linkage status, 
#' and the current parameters.
#'
#' @param data List with `encodedA`, `encodedB` (the two encoded data
#'   sources, missing values as `0`), `n_values`, `PIVs_config` and
#'   `same_mistakes`; see [StEM()].
#' @param links 2-column matrix of (A, B) indices for the currently linked
#'   records.
#' @param survivalpSameH Matrix (n links x n PIVs); `1` for stable PIVs, and the
#'   survival probability that the true value is unchanged for unstable PIVs.
#' @param sumRowD Logical vector, one entry per record in A: does it form a
#'   link?
#' @param sumColD Logical vector, one entry per record in B: does it form a
#'   link?
#' @param eta List (one per PIV) of the distribution of true values.
#' @param phi List (one per PIV) of length-4 vectors: agreement prob. in A,
#'   agreement prob. in B, missing prob. in A, missing prob. in B.
#'
#' @return List with `truepivsA` and `truepivsB`, matrices (same shape as the
#'   input data) of simulated true PIV values.
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V2 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V3 = list(dynamics = "flexible",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(NA,NA)),
#'                      V4 = list(dynamics = "structured",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(0.03,0.03),
#'                                cond_hazard_cov = list(cov1=c("Xe", "Xf"),
#'                                                     cov2=c())) )
#' n_values  <- c( 5, 6, 7, 12 )
#' p_mistake <- list( V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                   V3 = c(0.05, 0.05), V4 = c(0.02, 0.02) )
#' p_missing <- list( V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                   V3 = c(0.005, 0.005), V4 = c(0.005, 0.005) )
#' cond_hazard_params <- list( V1 = c(), V2 = c(), 
#'                             V3 = c(), V4 = log(c(0.7, 0.6, 0.5)) )
#' gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
#'                            p_mistake, p_missing, cond_hazard_params, TRUE )
#' data_StEM <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
#'                            PIVs_config, TRUE, "entity_id", TRUE )
#' PIVs_stable <- sapply(data_StEM$PIVs_config, function(x)
#'                         x$dynamics != "structured")
#' FlexRL:::initDeltaMap()
#' linksR = base::matrix(0,0,2)
#' linksCpp = linksR
#' sumRowD = rep(0, nrow(data_StEM$encodedA))
#' sumColD = rep(0, nrow(data_StEM$encodedB))
#' nlinkrec = 0
#' survivalpSameH = base::matrix(1, nrow(linksR), length(data_StEM$n_values))
#' gamma = 0.5
#' eta = lapply(data_StEM$n_values, function(x) rep(1/x,x))
#' phi = lapply(data_StEM$n_values, function(x)  c(0.9,0.9,0.1,0.1))
#' n_coef_unstable = lapply( seq_along(PIVs_stable), function(idx)
#'  if(PIVs_stable[idx]){ 0 }else{
#'    ncol(data_StEM$encodedA[, data_StEM$PIVs_config[[idx]]$cond_hazard_cov$covA,
#'                              drop=FALSE]) +
#'    ncol(data_StEM$encodedB[, data_StEM$PIVs_config[[idx]]$cond_hazard_cov$covB,
#'                              drop=FALSE]) + 1 } )
#' alpha = lapply( seq_along(PIVs_stable),
#'                 function(idx) if(PIVs_stable[idx]){ c(-Inf) }
#'                 else{ rep(log(0.05), n_coef_unstable[[idx]]) }
#'               )
#' newTruePivs = simulateH(data=data_StEM, links=linksCpp,
#'                         survivalpSameH=survivalpSameH,
#'                         sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
#' truepivsA = newTruePivs$truepivsA
#' truepivsB = newTruePivs$truepivsB
simulateH <- function(data, links, survivalpSameH, sumRowD, sumColD, eta, phi) {
  PIVs <- names(data$PIVs_config)
  PIVs_stable <- sapply(data$PIVs_config, function(x) x$dynamics != "structured")
  truePIVs <- sampleH(
    dim(data$encodedA[, PIVs, drop=FALSE]), dim(data$encodedB[, PIVs, drop=FALSE]),
    links, as.matrix(survivalpSameH), PIVs_stable,
    data$encodedA[, PIVs, drop=FALSE], data$encodedB[, PIVs, drop=FALSE],
    data$n_values, sumRowD == 0, sumColD == 0,
    eta, phi
  )
  list(truepivsA = truePIVs$truepivsA, truepivsB = truePIVs$truepivsB)
}

#' Simulate the linkage matrix D
#'
#' Given the current draw of true PIV values, samples the linkage matrix D from 
#' its conditional distribution, and returns the updated log-likelihood.
#'
#' @param data List with `encodedA`, `encodedB` (the two encoded data
#'   sources, missing values as `0`), `n_values`, `PIVs_config` and
#'   `same_mistakes`; see [StEM()].
#' @param linksR 2-column matrix (1-indexed) of the currently linked (A, B)
#'   indices.
#' @param sumRowD Logical vector, one entry per record in A: does it form a
#'   link?
#' @param sumColD Logical vector, one entry per record in B: does it form a
#'   link?
#' @param truepivsA Matrix of true PIV values, as returned by [simulateH()].
#' @param truepivsB Matrix of true PIV values, as returned by [simulateH()].
#' @param gamma Numeric, proportion of linked records as a fraction of the
#'   smaller file.
#' @param eta List (one per PIV) of the distribution of true values.
#' @param alpha List (one per PIV) of survival-model parameters
#'   (see [survival_model()]).
#' @param phi List (one per PIV) of registration-error parameters.
#' @param model_dynamics Object from [survival_model()] giving the survival
#'   function of the structured PIVs.
#'
#' @return A list (`Dsample`, Rcpp `sampleD()` output) with: `links` updated set
#'   of links, `sumRowD` updated sumRowD, `sumColD` updated sumColD, `loglik`
#'   updated value of the complete log likelihood, `nlinkrec` updated number of
#'   linked records
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V2 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V3 = list(dynamics = "flexible",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(NA,NA)),
#'                      V4 = list(dynamics = "structured",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(0.03,0.03),
#'                                cond_hazard_cov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())) )
#' n_values  <- c( 5, 6, 7, 12 )
#' p_mistake <- list( V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                    V3 = c(0.05, 0.05), V4 = c(0.02, 0.02) )
#' p_missing <- list( V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                    V3 = c(0.005, 0.005), V4 = c(0.005, 0.005) )
#' cond_hazard_params <- list( V1 = c(), V2 = c(), 
#'                             V3 = c(), V4 = log(c(0.7, 0.6, 0.5)) )
#' gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
#'                            p_mistake, p_missing, cond_hazard_params, TRUE )
#' data_StEM <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
#'                            PIVs_config, TRUE, "entity_id", TRUE )
#' PIVs_stable <- sapply(data_StEM$PIVs_config, function(x)
#'                         x$dynamics != "structured")
#' FlexRL:::initDeltaMap()
#' linksR = base::matrix(0,0,2)
#' linksCpp = linksR
#' sumRowD = rep(0, nrow(data_StEM$encodedA))
#' sumColD = rep(0, nrow(data_StEM$encodedB))
#' nlinkrec = 0
#' survivalpSameH = base::matrix(1, nrow(linksR), length(data_StEM$n_values))
#' gamma = 0.5
#' eta = lapply(data_StEM$n_values, function(x) rep(1/x,x))
#' phi = lapply(data_StEM$n_values, function(x)  c(0.9,0.9,0.1,0.1))
#' n_coef_unstable = lapply( seq_along(PIVs_stable), function(idx)
#'  if(PIVs_stable[idx]){ 0 }else{
#'    ncol(data_StEM$encodedA[, data_StEM$PIVs_config[[idx]]$cond_hazard_cov$covA,
#'                              drop=FALSE]) +
#'    ncol(data_StEM$encodedB[, data_StEM$PIVs_config[[idx]]$cond_hazard_cov$covB,
#'                              drop=FALSE]) + 1 } )
#' alpha = lapply( seq_along(PIVs_stable),
#'                 function(idx) if(PIVs_stable[idx]){ c(-Inf) }
#'                 else{ rep(log(0.05), n_coef_unstable[[idx]]) }
#'               )
#' newTruePivs = simulateH(data=data_StEM, links=linksCpp,
#'                         survivalpSameH=survivalpSameH,
#'                         sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
#' truepivsA = newTruePivs$truepivsA
#' truepivsB = newTruePivs$truepivsB
#' Dsample = simulateD(data=data_StEM, linksR=linksR, sumRowD=sumRowD,
#'                    sumColD=sumColD, truepivsA=truepivsA, truepivsB=truepivsB,
#'                    gamma=gamma, eta=eta, alpha=alpha, phi=phi)
#' linksCpp = Dsample$links
#' linksR = linksCpp + 1
simulateD <- function(data, linksR, sumRowD, sumColD, truepivsA, truepivsB, gamma, eta, alpha, phi,
                      model_dynamics = survival_model("exponential")) {
  PIVs <- names(data$PIVs_config)
  PIVs_stable <- sapply(data$PIVs_config, function(x) x$dynamics != "structured")

  # Restrict to pairs that are compatible on the stable PIVs (candidate links)
  UA <- pasteIntoPattern(as.matrix(truepivsA[, PIVs_stable]))
  UB <- pasteIntoPattern(as.matrix(truepivsB[, PIVs_stable]))
  valuesU <- unique(c(UA, UB))
  UA <- as.numeric(factor(UA, levels = valuesU))
  UB <- as.numeric(factor(UB, levels = valuesU))
  tmpA <- indexPatterns(UA, length(valuesU))
  tmpB <- indexPatterns(UB, length(valuesU))
  select <- pairPatterns(tmpA, tmpB, length(tmpA))

  pLink <- rep(gamma, nrow(data$encodedA[, PIVs, drop=FALSE]))

  # Contribution to the log-likelihood if a record from A is NOT linked
  LLA <- rep(0, nrow(data$encodedA[, PIVs, drop=FALSE]))
  for (k in seq_along(data$n_values)) {
    logpTrue <- log(eta[[k]])[truepivsA[, k]]
    pMissingA <- phi[[k]][3]
    pTypoA <- (1 - pMissingA) * (1 - phi[[k]][1]) / (data$n_values[k] - 1)
    pAgreeA <- (1 - pMissingA) * phi[[k]][1]
    contr <- rep(pAgreeA, nrow(data$encodedA[, PIVs, drop=FALSE]))
    contr[data$encodedA[, PIVs, drop=FALSE][, k] != truepivsA[, k]] <- pTypoA
    contr[data$encodedA[, PIVs, drop=FALSE][, k] == 0] <- pMissingA
    LLA <- LLA + logpTrue + log(contr)
  }

  # Contribution to the log-likelihood if a record from B is NOT linked
  LLB <- rep(0, nrow(data$encodedB[, PIVs, drop=FALSE]))
  for (k in seq_along(data$n_values)) {
    logpTrue <- log(eta[[k]])[truepivsB[, k]]
    pMissingB <- phi[[k]][4]
    pTypoB <- (1 - pMissingB) * (1 - phi[[k]][2]) / (data$n_values[k] - 1)
    pAgreeB <- (1 - pMissingB) * phi[[k]][2]
    contr <- rep(pAgreeB, nrow(data$encodedB[, PIVs, drop=FALSE]))
    contr[data$encodedB[, PIVs, drop=FALSE][, k] != truepivsB[, k]] <- pTypoB
    contr[data$encodedB[, PIVs, drop=FALSE][, k] == 0] <- pMissingB
    LLB <- LLB + logpTrue + log(contr)
  }

  # Contribution to the log-likelihood if a candidate pair IS linked
  LLL <- Matrix::Matrix(0, nrow = nrow(data$encodedA[, PIVs, drop=FALSE]), ncol = nrow(data$encodedB[, PIVs, drop=FALSE]), sparse = TRUE)
  for (k in seq_along(data$n_values)) {
    HA <- truepivsA[select[, 1], k]
    HB <- truepivsB[select[, 2], k]
    logpTrue <- log(eta[[k]])[HA]
    pMissingA <- phi[[k]][3]
    pTypoA <- (1 - pMissingA) * (1 - phi[[k]][1]) / (data$n_values[k] - 1)
    pAgreeA <- (1 - pMissingA) * phi[[k]][1]
    pMissingB <- phi[[k]][4]
    pTypoB <- (1 - pMissingB) * (1 - phi[[k]][2]) / (data$n_values[k] - 1)
    pAgreeB <- (1 - pMissingB) * phi[[k]][2]
    # Contribution to the likelihood of linked observation from A
    helpA <- rep(pAgreeA, length(HA))
    helpA[data$encodedA[, PIVs, drop=FALSE][select[, 1], k, drop=FALSE] != HA] <- pTypoA
    helpA[data$encodedA[, PIVs, drop=FALSE][select[, 1], k, drop=FALSE] == 0] <- pMissingA
    # Contribution to the likelihood of linked observation from B
    helpB <- rep(pAgreeB, length(HB))
    helpB[data$encodedB[, PIVs, drop=FALSE][select[, 2], k, drop=FALSE] != HB] <- pTypoB
    helpB[data$encodedB[, PIVs, drop=FALSE][select[, 2], k, drop=FALSE] == 0] <- pMissingB

    LLL[select] <- LLL[select] + logpTrue + log(helpA) + log(helpB)

    # Add unstable part if unstable
    if (!PIVs_stable[k]) {
      times <- abs(data$encodedB[select[, 2], "date"] - data$encodedA[select[, 1], "date"])
      intercept <-  rep(1, nrow(select))
      cov_k <- cbind(
        intercept,
        data$encodedA[select[, 1], data$PIVs_config[[k]]$cond_hazard_cov$covA, drop = FALSE],
        data$encodedB[select[, 2], data$PIVs_config[[k]]$cond_hazard_cov$covB, drop = FALSE]
      )
      pSameH <- model_dynamics$S(cov_k, alpha[[k]], times)
      helpH <- pSameH^(HA == HB) * ((1 - pSameH) / (data$n_values[k] - 1))^(HA != HB)
      LLL[select] <- LLL[select] + log(helpH)
    }
  }

  if (anyNA(LLL[select])) {
    warning("Some entries of the linkage log-likelihood matrix are NA.", call. = FALSE)
  }

  # Complete data likelihood
  LL0 <- log_lik(LLL = LLL, LLA = LLA, LLB = LLB, links = linksR, sumRowD = sumRowD, sumColD = sumColD, gamma = gamma)

  Dsample <- sampleD(as.matrix(select), LLA, LLB, LLL[select], pLink, LL0, as.integer(nrow(linksR)), sumRowD > 0, sumColD > 0)
  linksR <- Dsample$links + 1

  # Sanity check: recomputed log-likelihood should match the sampler's own value
  ll_check <- log_lik(LLL = LLL, LLA = LLA, LLB = LLB, links = linksR, sumRowD = Dsample$sumRowD, sumColD = Dsample$sumColD, gamma = pLink)
  if (round(Dsample$loglik, 3) != round(ll_check, 3)) {
    stop(
      "Log-likelihood sanity check failed (sampler: ", round(Dsample$loglik, 3),
      ", recomputed: ", round(ll_check, 3), ").",
      call. = FALSE
    )
  }
  Dsample
}

#' Survival models for the dynamics of a PIV
#'
#' A structured PIV may change between the two registrations. The probability
#' that the true value of a linked pair is unchanged after a time gap `t` is
#' modelled by a survival function `S(t | X, alpha)`, where `X` holds an
#' intercept and the covariates given in `cond_hazard_cov` and `alpha` the
#' parameters estimated by [StEM()]. This function builds the model object
#' used by [StEM()], [simulateD()] and [simulate_data()].
#'
#' Available models (`h` is the hazard, `lambda = exp(X alpha_cov)` the
#' proportional-hazards term):
#' \describe{
#'   \item{`"exponential"`}{`S(t) = exp(-lambda t)`; `alpha` = coefficients of
#'     `X` (default, the model of the methodology paper).}
#'   \item{`"weibull"`}{`S(t) = exp(-(lambda t)^k)` with shape
#'     `k = exp(alpha[1])`; `alpha` = log-shape, then coefficients of `X`.}
#'   \item{`"gompertz"`}{`h(t) = lambda exp(g t)`, so
#'     `S(t) = exp(-lambda (exp(g t) - 1) / g)` with `g = alpha[1]`;
#'     `alpha` = `g`, then coefficients of `X`.}
#'   \item{`"piecewise"`}{piecewise-constant baseline hazard on the intervals
#'     defined by `cuts`, times `lambda`; `alpha` = one log-hazard per
#'     interval, then coefficients of the covariates (no intercept).}
#'   \item{`"custom"`}{`S`, `n_par` and `init` supplied by the user.}
#' }
#' For every model the negative log-likelihood used in the M-step is
#' `-sum(Hequal log S + (1 - Hequal) log(1 - S))`, where `Hequal` indicates
#' that the true values of the linked pair agree. Parameters are estimated
#' with [stats::nlminb()] (numerical gradient).
#'
#' @param type One of `"exponential"`, `"weibull"`, `"gompertz"`,
#'   `"piecewise"`, `"custom"`.
#' @param cuts Numeric vector of cut points for `type = "piecewise"`
#'   (e.g. `c(1, 3)` gives three intervals `[0,1)`, `[1,3)`, `[3, Inf)`).
#' @param S For `type = "custom"`: function `(X, alpha, times)` returning the
#'   survival probability of each linked pair; `X` is a matrix with an
#'   intercept column first, then the covariates.
#' @param n_par For `type = "custom"`: function of `ncol(X)` returning the
#'   length of `alpha`.
#' @param init For `type = "custom"`: function of `ncol(X)` returning the
#'   starting values of `alpha`.
#'
#' @return An object of class `"survival_model"`: a list with `type`, `S`,
#'   `negloglik`, `n_par`, `init` and `par_names`.
#' @export
#'
#' @examples
#' X <- cbind(intercept = rep(1, 5))
#' times <- c(0.001, 0.2, 1.3, 1.5, 2)
#' Hequal <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
#'
#' expo <- survival_model("exponential")
#' expo$S(X, alpha = log(0.3), times)
#' stats::nlminb(expo$init(ncol(X)), expo$negloglik, X = X, times = times,
#'               Hequal = Hequal)$par
#'
#' weib <- survival_model("weibull")
#' weib$S(X, alpha = c(0.5, -1), times)
#' stats::nlminb(weib$init(ncol(X)), weib$negloglik, X = X, times = times,
#'               Hequal = Hequal)$par
#'
#' gomp <- survival_model("gompertz")
#' gomp$S(X, alpha = c(0.5, -1), times)
#' stats::nlminb(gomp$init(ncol(X)), gomp$negloglik, X = X, times = times,
#'               Hequal = Hequal)$par
#'               
#' pw <- survival_model("piecewise", cuts = 1)
#' pw$S(X, alpha = c(-2, -1), times)
#' stats::nlminb(pw$init(ncol(X)), pw$negloglik, X = X, times = times,
#'               Hequal = Hequal)$par
#'
#' # a custom model: log-logistic with scale exp(X alpha[-1]) 
#' # and shape exp(alpha[1])
#' loglogistic <- survival_model("custom",
#'   S = function(X, alpha, times) {
#'     lambda <- exp(as.matrix(X) %*% alpha[-1])
#'     1 / ( 1 + (lambda * times)^exp(alpha[1]) )
#'   },
#'   n_par = function(n_cov) n_cov + 1,
#'   init = function(n_cov) c(stats::rnorm(1, 0, 0.1), stats::runif(n_cov, log(0.01), log(1))))
#' loglogistic$S(X, alpha = c(0, -2), times)
#' stats::nlminb(loglogistic$init(ncol(X)), loglogistic$negloglik, X = X, times = times,
#'               Hequal = Hequal)$par
survival_model <- function(type = c("exponential", "weibull", "gompertz", "piecewise", "custom"),
                           cuts = NULL, S = NULL, n_par = NULL, init = NULL) {
  type <- match.arg(type)
  ph <- function(X, beta) as.vector(exp(as.matrix(X) %*% beta))
  
  # Random starting values on a sensible scale: log baseline hazard for a rate
  # between 0.01 and 1 per time unit, covariate coefficients near 0 (no effect).
  log_h0 <- function(n = 1) stats::runif(n, log(0.01), log(1))
  beta0  <- function(n)     stats::rnorm(n, 0, 0.1)
  
  model <- switch(type,
                  exponential = list(
                    S = function(X, alpha, times) exp(-ph(X, alpha) * times),
                    n_par = function(n_cov) n_cov,
                    init = function(n_cov) c(log_h0(), beta0(n_cov - 1)),
                    par_names = function(cov_names) cov_names
                  ),
                  weibull = list(
                    S = function(X, alpha, times) exp(-(ph(X, alpha[-1]) * times)^exp(alpha[1])),
                    n_par = function(n_cov) n_cov + 1,
                    init = function(n_cov) c(beta0(1), log_h0(), beta0(n_cov - 1)),
                    par_names = function(cov_names) c("log_shape", cov_names)
                  ),
                  gompertz = list(
                    S = function(X, alpha, times) {
                      g <- alpha[1]
                      if (abs(g) < 1e-8) return(exp(-ph(X, alpha[-1]) * times))
                      exp(-ph(X, alpha[-1]) * (exp(g * times) - 1) / g)
                    },
                    n_par = function(n_cov) n_cov + 1,
                    init = function(n_cov) c(beta0(1), log_h0(), beta0(n_cov - 1)),
                    par_names = function(cov_names) c("g", cov_names)
                  ),
                  piecewise = {
                    if (is.null(cuts) || any(diff(c(0, cuts)) <= 0)) {
                      stop("`cuts` must be an increasing vector of positive cut points for `type = \"piecewise\"`.", call. = FALSE)
                    }
                    edges <- c(0, cuts, Inf)
                    n_int <- length(cuts) + 1
                    list(
                      S = function(X, alpha, times) {
                        h <- exp(alpha[seq_len(n_int)])
                        # cumulative baseline hazard: time spent in each interval times its hazard
                        H0 <- vapply(times, function(t) sum(h * pmax(0, pmin(t, edges[-1]) - edges[-length(edges)])), numeric(1))
                        beta <- alpha[-seq_len(n_int)]
                        lam <- if (length(beta)) ph(as.matrix(X)[, -1, drop = FALSE], beta) else 1
                        exp(-lam * H0)
                      },
                      n_par = function(n_cov) n_int + n_cov - 1,
                      init = function(n_cov) c(log_h0(n_int), beta0(n_cov - 1)),
                      par_names = function(cov_names) c(paste0("log_h", seq_len(n_int)), cov_names[-1])
                    )
                  },
                  custom = {
                    if (!is.function(S) || !is.function(n_par) || !is.function(init)) {
                      stop("`S`, `n_par` and `init` must be functions for `type = \"custom\"`.", call. = FALSE)
                    }
                    list(S = S, n_par = n_par, init = init,
                         par_names = function(cov_names) paste0("alpha", seq_len(n_par(length(cov_names)))))
                  }
  )
  model$type <- type
  model$negloglik <- function(alpha, X, times, Hequal) {
    s <- pmin(pmax(model$S(X, alpha, times), 1e-12), 1 - 1e-12)
    -sum(Hequal * log(s) + (!Hequal) * log(1 - s))
  }
  structure(model, class = "survival_model")
}

#' Stochastic Expectation-Maximisation for record linkage
#'
#' Fits the FlexRL model with a Stochastic EM algorithm: each iteration runs
#' a Gibbs sampler (alternating between simulating the true PIV values and
#' the linkage matrix D) and then updates the model parameters (`gamma`,
#' `eta`, `alpha`, `phi`) from the post-burn-in Gibbs draws. See the
#' methodology paper https://10.1093/jrsssc/qlaf016 for details.
#'
#' @param data List, typically the output of [prepare_data()], with:
#'   `encodedA` (smaller source, PIVs encoded to natural numbers,
#'    `0` = missing), `encodedB` (larger source, encoded), `n_values`,
#'    `same_mistakes` (logical: same mistake parameter shared by A and B?),
#'    `PIVs_config` (named list, one entry per PIV, each a list with: `dynamics`
#'    (`"stable"`, `"flexible"`, or `"structured"`), `bound_mistakes` (length-2
#'    numeric/NA, upper bound on the mistake probability in file 1 / file 2),
#'    `fix_mistakes` (length-2 numeric/NA, mistake probability fixed to this
#'    value in file 1 / file 2), and, only for `dynamics = "structured"`,
#'    `cond_hazard_cov` (a list with `cov1` and `cov2`, the names of covariates
#'    in file 1 / file 2 used to model the hazard of change).
#' @param StEM_iter Integer, total number of StEM iterations (including burn-in).
#' @param StEM_burnin Integer, number of StEM iterations discarded as burn-in.
#' @param gibbs_iter Integer, total number of Gibbs iterations per StEM step
#'   (including burn-in).
#' @param gibbs_burnin Integer, number of Gibbs iterations discarded as burn-in
#'   (`0` lets the algorithm auto-detect burn-in from the stabilisation of the
#'   linked count).
#' @param music_on Logical; if `TRUE`, opens a short tune in the browser when the
#'   algorithm finishes.
#' @param new_directory Path to an existing directory to save progress after
#'   each iteration, or `NULL` to disable.
#' @param save_info_iter Logical; save the environment at the end of each
#'   iteration (only used if `new_directory` is not `NULL`).
#' @param gamma0 Optional starting values for `gamma`; at random if `NULL`.
#' @param phiA0 Optional starting values for `phi`; at random if `NULL`.
#' @param phiB0 Optional starting values for `phi`; at random if `NULL`.
#' @param n_post_sample Integer, number of posterior draws used to estimate the
#'   final linkage probabilities `Delta`. Default is set to 1000, we recommend
#'   not lowering it.
#' @param model_dynamics Object from [survival_model()]: the survival model
#'   for the change over time of the structured PIVs (default exponential).
#'
#' @return A list with: `Delta` sparse-matrix summary (`i`, `j`, `x`) of
#'   posterior linkage probabilities; a pair is a valid link candidate once
#'   `x > 0.5` (one-to-one constraint), `gamma`, `eta`, `alpha`, `phi` the StEM
#'   chains for each parameter.
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V2 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V3 = list(dynamics = "flexible",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(NA,NA)),
#'                      V4 = list(dynamics = "structured",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(0.03,0.03),
#'                                cond_hazard_cov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())) )
#' n_values  <- c( 5, 6, 7, 12 )
#' p_mistake <- list( V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                    V3 = c(0.05, 0.05), V4 = c(0.02, 0.02) )
#' p_missing <- list( V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                    V3 = c(0.005, 0.005), V4 = c(0.005, 0.005) )
#'                    
#' cond_hazard_params1 <- list( V1 = c(), V2 = c(), 
#'                              V3 = c(), V4 = log(c(0.7, 0.6, 0.5)) )
#' gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
#'                            p_mistake, p_missing, cond_hazard_params1, TRUE )
#' prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
#'                            PIVs_config, TRUE, "entity_id", TRUE )
#' fit1 <- StEM( data = prep_data, StEM_iter = 5, StEM_burnin = 2, 
#'               gibbs_iter = 5, gibbs_burnin = 2, n_post_sample = 5,
#'               model_dynamics = survival_model("exponential") )
#' apply(fit1$alpha$V4[2:5,], 2, mean)
#' cond_hazard_params1$V4
#' head(fit1$Delta[fit1$Delta$x > 0.5, ])
#' 
#' cond_hazard_params2 <- list( V1 = c(), V2 = c(), 
#'                              V3 = c(), V4 = log(c(0.3, 0.7, 0.6, 0.5)) )
#' gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
#'                            p_mistake, p_missing, cond_hazard_params2, TRUE,
#'                            model_dynamics = survival_model("weibull") )
#' prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
#'                            PIVs_config, TRUE, "entity_id", TRUE )
#' fit2 <- StEM( data = prep_data, StEM_iter = 5, StEM_burnin = 2, 
#'               gibbs_iter = 5, gibbs_burnin = 2, n_post_sample = 5,
#'               model_dynamics = survival_model("weibull"))
#' apply(fit2$alpha$V4[2:5,], 2, mean)
#' cond_hazard_params2$V4
#' head(fit2$Delta[fit2$Delta$x > 0.5, ])
#' 
#' cond_hazard_params3 <- list( V1 = c(), V2 = c(), 
#'                              V3 = c(), V4 = log(c(0.3, 0.7, 0.6, 0.5)) )
#' gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
#'                            p_mistake, p_missing, cond_hazard_params3, TRUE,
#'                            model_dynamics = survival_model("gompertz") )
#' prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
#'                            PIVs_config, TRUE, "entity_id", TRUE )
#' fit3 <- StEM( data = prep_data, StEM_iter = 5, StEM_burnin = 2, 
#'               gibbs_iter = 5, gibbs_burnin = 2, n_post_sample = 5,
#'               model_dynamics = survival_model("gompertz"))
#' apply(fit3$alpha$V4[2:5,], 2, mean)
#' cond_hazard_params3$V4
#' head(fit3$Delta[fit3$Delta$x > 0.5, ])
#' 
#' cond_hazard_params4 <- list( V1 = c(), V2 = c(), 
#'                              V3 = c(), V4 = log(c(0.6, 0.5, 0.4, 0.3, 0.2, 0.1)) )
#' gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
#'                            p_mistake, p_missing, cond_hazard_params4, TRUE,
#'                            model_dynamics = survival_model("piecewise", cuts = c(1,2,3)) )
#' prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
#'                            PIVs_config, TRUE, "entity_id", TRUE )
#' fit4 <- StEM( data = prep_data, StEM_iter = 2, StEM_burnin = 1, 
#'               gibbs_iter = 2, gibbs_burnin = 1, n_post_sample = 2,
#'               model_dynamics = survival_model("piecewise", cuts = c(1,2,3)))
#' apply(fit4$alpha$V4, 2, mean)
#' cond_hazard_params4$V4
#' head(fit4$Delta[fit4$Delta$x > 0.5, ])
StEM <- function(data, StEM_iter = 30, StEM_burnin = 15, gibbs_iter = 20, gibbs_burnin = 10,
                 music_on = FALSE, new_directory = NULL, save_info_iter = FALSE,
                 gamma0 = NULL, phiA0 = NULL, phiB0 = NULL, n_post_sample = 1000,
                 model_dynamics = survival_model("exponential")) {

  message("FlexRL")

  required <- c("encodedA", "encodedB", "PIVs_config", "n_values", "same_mistakes")
  if (!all(required %in% names(data))) {
    stop("`data` must contain: ", paste(required, collapse = ", "), ".", call. = FALSE)
  }

  nGibbsIter <- gibbs_iter - gibbs_burnin
  if (StEM_iter - StEM_burnin <= 0) stop("`StEM_iter` must be greater than `StEM_burnin`.", call. = FALSE)
  if (nGibbsIter <= 0) stop("`gibbs_iter` must be greater than `gibbs_burnin`.", call. = FALSE)

  encodedA <- data$encodedA
  encodedB <- data$encodedB

  PIVs <- names(data$PIVs_config)
  PIVs_stable <- sapply(data$PIVs_config, function(x) x$dynamics != "structured")
  PIVs_boundMistakes <- sapply(data$PIVs_config, function(x) x$bound_mistakes)
  PIVs_fixMistakes <- sapply(data$PIVs_config, function(x) x$fix_mistakes)
  PIVs_conditionalHazardVariables <- sapply(data$PIVs_config, function(x) {
    if (x$dynamics == "structured") x$cond_hazard_cov else FALSE
  })

  N_PIVs <- length(PIVs)
  modelDynaPIVs <- which(!PIVs_stable)

  if (length(modelDynaPIVs) > 0) {
    if (!"date" %in% colnames(encodedA) || !"date" %in% colnames(encodedB)) {
      stop("Modelling unstable PIVs requires a `date` column in both file A and file B.", call. = FALSE)
    } else {
      tryCatch(
        { times_test = abs( mean(encodedB$date,na.rm=TRUE) - mean(encodedA$date,na.rm=TRUE) )
        }, error = function(msg){
          warning("Error in the `date` columns in encodedA or encodedB, format does not allow to compute absolute time gaps.\n")
        })
    }
    for (k in modelDynaPIVs) {
      covA <- PIVs_conditionalHazardVariables[[k]]$covA
      covB <- PIVs_conditionalHazardVariables[[k]]$covB
      if (length(covA) > 0 && !all(covA %in% colnames(encodedA))) {
        stop("Some `cond_hazard_cov$covA` variables for PIV '", PIVs[k], "' are missing from file A.", call. = FALSE)
      }
      if (length(covB) > 0 && !all(covB %in% colnames(encodedB))) {
        stop("Some `cond_hazard_cov$covB` variables for PIV '", PIVs[k], "' are missing from file B.", call. = FALSE)
      }
    }
  }

  # =Initial parameter values
  gamma <- if (is.null(gamma0)) stats::runif(1, 0.2, 0.8) else gamma0
  eta <- lapply(data$n_values, function(x) rep(1 / x, x))

  # n_coef_unstable = ncol(X) (intercept + covariates)
  n_coef_unstable <- lapply(seq_along(PIVs_stable), function(idx) {
    if (PIVs_stable[idx]) 0 else {
      ncol(encodedA[, PIVs_conditionalHazardVariables[[idx]]$covA, drop = FALSE]) +
        ncol(encodedB[, PIVs_conditionalHazardVariables[[idx]]$covB, drop = FALSE]) + 1
    }
  })
  if (!inherits(model_dynamics, "survival_model")) {
    stop("`model_dynamics` must be built with `survival_model()`.", call. = FALSE)
  }
  # n_alpha = length(alpha)
  n_alpha <- lapply(n_coef_unstable, function(n) if (n == 0) 1 else model_dynamics$n_par(n))
  alpha <- lapply(seq_along(PIVs_stable), function(idx) {
    if (PIVs_stable[idx]) -Inf else model_dynamics$init(n_coef_unstable[[idx]])
  })

  # (agreement in A, agreement in B, missing in A, missing in B)
  if (!is.null(phiA0) && !is.null(phiB0)) {
    phi <- lapply(data$n_values, function(x) c(phiA0, phiB0, 0.1, 0.1))
  } else if (!is.null(phiA0)) {
    phi <- lapply(data$n_values, function(x) c(phiA0, stats::runif(1, 0.8, 0.97), 0.1, 0.1))
  } else if (!is.null(phiB0)) {
    phi <- lapply(data$n_values, function(x) c(stats::runif(1, 0.8, 0.97), phiB0, 0.1, 0.1))
  } else {
    phi <- lapply(data$n_values, function(x) c(stats::runif(1, 0.8, 0.97), stats::runif(1, 0.8, 0.97), 0.1, 0.1))
  }

  NmissingA <- lapply(seq_along(data$n_values), function(k) sum(encodedA[, PIVs, drop=FALSE][, k] == 0))
  NmissingB <- lapply(seq_along(data$n_values), function(k) sum(encodedB[, PIVs, drop=FALSE][, k] == 0))

  gamma.iter <- array(NA, c(StEM_iter, length(gamma)))
  eta.iter   <- lapply(data$n_values, function(x) array(NA, c(StEM_iter, x)))
  phi.iter   <- lapply(data$n_values, function(x) array(NA, c(StEM_iter, 4)))
  alpha.iter <- lapply(seq_along(n_alpha), function(k) {
    m <- array(NA, c(StEM_iter, n_alpha[[k]]))
    if (!PIVs_stable[k]) {
      colnames(m) <- model_dynamics$par_names(c("intercept",
                                                PIVs_conditionalHazardVariables[[k]]$covA,
                                                PIVs_conditionalHazardVariables[[k]]$covB))
    }
    m
  })
  names(alpha.iter) <- PIVs

  boundHitCountA <- stats::setNames(rep(0L, N_PIVs), PIVs)
  boundHitCountB <- stats::setNames(rep(0L, N_PIVs), PIVs)
  boundHitThreshold <- 0.25

  pb_id <- cli::cli_progress_bar(
    format = "Running StEM algorithm {cli::pb_bar} {cli::pb_percent} | iter {cli::pb_current}/{cli::pb_total} [{cli::pb_elapsed}]",
    total = StEM_iter, clear = FALSE
  )

  # Runs one Gibbs iteration and refreshes `survivalpSameH` for unstable PIVs
  gibbs_step <- function(linksCpp, linksR, sumRowD, sumColD, survivalpSameH, gamma, eta, alpha, phi) {
    newTruePivs <- simulateH(data = data, links = linksCpp, survivalpSameH = survivalpSameH,
                             sumRowD = sumRowD, sumColD = sumColD, eta = eta, phi = phi)
    Dsample <- simulateD(data = data, linksR = linksR, sumRowD = sumRowD, sumColD = sumColD,
                         truepivsA = newTruePivs$truepivsA, truepivsB = newTruePivs$truepivsB,
                         gamma = gamma, eta = eta, alpha = alpha, phi = phi,
                         model_dynamics = model_dynamics)
    linksCpp <- Dsample$links
    linksR <- linksCpp + 1
    survivalpSameH <- matrix(1, nrow(linksR), N_PIVs)
    if (length(modelDynaPIVs) > 0 && nrow(linksR) > 0) {
      times <- abs(encodedB[linksR[, 2], "date"] - encodedA[linksR[, 1], "date"])
      intercept <- rep(1, nrow(linksR))
      for (k in modelDynaPIVs) {
        cov_k <- cbind(intercept,
                       encodedA[linksR[, 1], PIVs_conditionalHazardVariables[[k]]$covA, drop = FALSE],
                       encodedB[linksR[, 2], PIVs_conditionalHazardVariables[[k]]$covB, drop = FALSE])
        survivalpSameH[, k] <- model_dynamics$S(cov_k, alpha[[k]], times)
      }
    }
    list(linksCpp = linksCpp, linksR = linksR, sumRowD = Dsample$sumRowD, sumColD = Dsample$sumColD,
         survivalpSameH = survivalpSameH, truepivsA = newTruePivs$truepivsA, truepivsB = newTruePivs$truepivsB,
         loglik = Dsample$loglik, nlinkrec = Dsample$nlinkrec)
  }
  run_auto_burnin <- function(linksCpp, linksR, sumRowD, sumColD, survivalpSameH, nlinkrec, gamma, eta, alpha, phi) {
    burnin_id <- cli::cli_progress_bar(
      format = "  Burn-in (auto) {cli::pb_spin} iter {Burnin_total} | linked records: {nlinkrec}",
      total = 500, current = FALSE
    )
    countBurnin <- 10
    Burnin_total <- 0
    while (countBurnin != 0) {
      Burnin_total <- Burnin_total + 1
      step <- gibbs_step(linksCpp, linksR, sumRowD, sumColD, survivalpSameH, gamma, eta, alpha, phi)
      linksCpp <- step$linksCpp
      linksR <- step$linksR
      sumRowD <- step$sumRowD
      sumColD <- step$sumColD
      survivalpSameH <- step$survivalpSameH
      countBurnin <- if (abs(step$nlinkrec - nlinkrec) < 10) countBurnin - 1 else 15
      nlinkrec <- step$nlinkrec
      cli::cli_progress_update(id = burnin_id)
      if (Burnin_total >= 500) {
        cli::cli_progress_done(id = burnin_id)
        stop("Auto burn-in exceeded 500 iterations without stabilising; set `gibbs_burnin` explicitly.", call. = FALSE)
      }
    }
    cli::cli_progress_done(id = burnin_id)
    list(linksCpp = linksCpp, linksR = linksR, sumRowD = sumRowD, sumColD = sumColD,
         survivalpSameH = survivalpSameH, nlinkrec = nlinkrec)
  }

  for (iter in seq_len(StEM_iter)) {

    initDeltaMap()
    linksR <- matrix(0, 0, 2)
    linksCpp <- linksR
    sumRowD <- rep(0, nrow(encodedA))
    sumColD <- rep(0, nrow(encodedB))
    nlinkrec <- 0
    survivalpSameH <- matrix(1, nrow(linksR), N_PIVs)

    # Burn-in
    if (gibbs_burnin == 0) {
      # Auto burn-in: keep sampling until the number of linked records stops
      # moving by more than 10 for 10 consecutive iterations, capped at 500.
      burnin <- run_auto_burnin(linksCpp, linksR, sumRowD, sumColD, survivalpSameH, nlinkrec, gamma, eta, alpha, phi)
      linksCpp <- burnin$linksCpp
      linksR <- burnin$linksR
      sumRowD <- burnin$sumRowD
      sumColD <- burnin$sumColD
      survivalpSameH <- burnin$survivalpSameH
      nlinkrec <- burnin$nlinkrec
    } else {
      for (j in seq_len(gibbs_burnin)) {
        step <- gibbs_step(linksCpp, linksR, sumRowD, sumColD, survivalpSameH, gamma, eta, alpha, phi)
        linksCpp <- step$linksCpp
        linksR <- step$linksR
        sumRowD <- step$sumRowD
        sumColD <- step$sumColD
        survivalpSameH <- step$survivalpSameH
      }
    }

    # Post burn-in: accumulate sufficient statistics for the M-step
    Vgamma <- c()
    Veta <- lapply(data$n_values, function(x) c())
    Valpha <- mapply(create_data_alpha, n_coef_unstable = n_coef_unstable, stable = PIVs_stable, SIMPLIFY = FALSE)
    Vphi <- lapply(data$n_values, function(x) c())

    for (j in seq_len(nGibbsIter)) {
      step <- gibbs_step(linksCpp, linksR, sumRowD, sumColD, survivalpSameH, gamma, eta, alpha, phi)
      linksCpp <- step$linksCpp
      linksR <- step$linksR
      sumRowD <- step$sumRowD
      sumColD <- step$sumColD
      survivalpSameH <- step$survivalpSameH
      truepivsA <- step$truepivsA
      truepivsB <- step$truepivsB

      if (length(modelDynaPIVs) > 0 && nrow(linksR) > 0) {
        times <- abs(encodedB[linksR[, 2], "date"] - encodedA[linksR[, 1], "date"])
        intercept <- rep(1, nrow(linksR))
        for (k in modelDynaPIVs) {
          cov_k <- cbind(intercept,
                         encodedA[linksR[, 1], PIVs_conditionalHazardVariables[[k]]$covA, drop = FALSE],
                         encodedB[linksR[, 2], PIVs_conditionalHazardVariables[[k]]$covB, drop = FALSE])
          Hequal <- truepivsA[linksR[, 1], k] == truepivsB[linksR[, 2], k]
          Valpha[[k]] <- rbind(Valpha[[k]], cbind(cov_k, times = times, Hequal = Hequal))
        }
      }

      # Update gamma
      Vgamma[j] <- nrow(linksR)
      if (nrow(linksR) == 0) {
        warning("No link was made in a Gibbs iteration (log-lik = ", step$loglik,
                ", gamma so far = ", gamma, "); check for a support mismatch between A and B.", call. = FALSE)
      }

      # Update eta and phi
      for (k in seq_len(N_PIVs)) {
        facpivsA <- factor(truepivsA[, k], levels = 1:data$n_values[k])
        facpivsB <- factor(truepivsB[, k], levels = 1:data$n_values[k])
        Veta[[k]] <- rbind(Veta[[k]], table(facpivsA[sumRowD == 0]) + table(facpivsB[sumColD == 0]) + table(facpivsA[sumRowD == 1]))
        Vphi[[k]] <- rbind(Vphi[[k]], c(
          sum(truepivsA[, k] == encodedA[, PIVs, drop=FALSE][, k, drop=FALSE]),
          sum(truepivsB[, k] == encodedB[, PIVs, drop=FALSE][, k, drop=FALSE]),
          NmissingA[[k]], NmissingB[[k]]
        ))
      }
    }

    # M-step
    gamma <- sum(Vgamma) / (nGibbsIter * nrow(truepivsA))
    if (gamma == 1) { gamma <- 0.99; warning("`gamma` hit 1 and was reset to 0.99.", call. = FALSE) }
    if (gamma == 0) { gamma <- 0.01; warning("`gamma` hit 0 and was reset to 0.01.", call. = FALSE) }

    for (k in seq_len(N_PIVs)) {
      eta[[k]] <- colSums(Veta[[k]]) / sum(Veta[[k]])
      if (any(eta[[k]] == 0)) {
        warning("PIV '", PIVs[k], "' has one or more values with zero estimated frequency; check for rare categories.", call. = FALSE)
      }
    }

    if (length(modelDynaPIVs) > 0 && nrow(linksR) > 0) {
      for (k in modelDynaPIVs) {
        X <- as.matrix(Valpha[[k]][, 1:n_coef_unstable[[k]], drop = FALSE])
        opt <- stats::nlminb(model_dynamics$init(n_coef_unstable[[k]]), model_dynamics$negloglik,
                             X = X, times = Valpha[[k]]$times, Hequal = Valpha[[k]]$Hequal)
        alpha[[k]] <- opt$par
      }
    }

    for (k in seq_len(N_PIVs)) {
      NtotalA <- nGibbsIter * nrow(encodedA)
      NtotalB <- nGibbsIter * nrow(encodedB)
      Ntotal  <- NtotalA + NtotalB

      if (data$same_mistakes) {
        phi[[k]][1] <- (sum(Vphi[[k]][, 1]) + sum(Vphi[[k]][, 2])) / (Ntotal - sum(Vphi[[k]][, 3]) - sum(Vphi[[k]][, 4]))
        phi[[k]][2] <-(sum(Vphi[[k]][, 1]) + sum(Vphi[[k]][, 2])) / (Ntotal - sum(Vphi[[k]][, 3]) - sum(Vphi[[k]][, 4]))
      } else {
        phi[[k]][1] <- sum(Vphi[[k]][, 1]) / (NtotalA - sum(Vphi[[k]][, 3]))
        phi[[k]][2] <- sum(Vphi[[k]][, 2]) / (NtotalB - sum(Vphi[[k]][, 4]))
      }

      if (!is.na(PIVs_fixMistakes[,k][1])) phi[[k]][1] <- 1 - PIVs_fixMistakes[,k][1]
      if (!is.na(PIVs_fixMistakes[,k][2])) phi[[k]][2] <- 1 - PIVs_fixMistakes[,k][2]

      if (!is.na(PIVs_boundMistakes[,k][1])) {
        if (phi[[k]][1] < 1 - PIVs_boundMistakes[,k][1]) {
          phi[[k]][1] <- 1 - PIVs_boundMistakes[,k][1]
          boundHitCountA[k] <- boundHitCountA[k] + 1L
        }
      }
      if (!is.na(PIVs_boundMistakes[,k][2])) {
        if (phi[[k]][2] < 1 - PIVs_boundMistakes[,k][2]) {
          phi[[k]][2] <- 1 - PIVs_boundMistakes[,k][2]
          boundHitCountB[k] <- boundHitCountB[k] + 1L
        }
      }

      phi[[k]][3] <- sum(Vphi[[k]][, 3]) / NtotalA
      phi[[k]][4] <- sum(Vphi[[k]][, 4]) / NtotalB
    }

    gamma.iter[iter, ] <- gamma
    for (k in seq_len(N_PIVs)) {
      eta.iter[[k]][iter, ] <- eta[[k]]
      alpha.iter[[k]][iter, ] <- alpha[[k]]
      phi.iter[[k]][iter, ] <- phi[[k]]
    }

    cli::cli_progress_update(id = pb_id)

    if (!is.null(new_directory) && save_info_iter) {
      save.image(file=file.path(new_directory, "myEnvironment.RData"))
      save(iter, gamma, eta, alpha, phi, gamma.iter, eta.iter, phi.iter, alpha.iter,
           file = file.path(new_directory, "myEnvironmentLocal.RData"))
    }
  }

  freqA <- boundHitCountA / StEM_iter
  freqB <- boundHitCountB / StEM_iter
  warning_msg_bound_phi <- c(
    sprintf("PIV '%s', file A hit its mistake bound on %.0f%% of StEM iterations.", PIVs[freqA > boundHitThreshold], 100 * freqA[freqA > boundHitThreshold]),
    sprintf("PIV '%s', file B hit its mistake bound on %.0f%% of StEM iterations.", PIVs[freqB > boundHitThreshold], 100 * freqB[freqB > boundHitThreshold])
  )
  if (length(warning_msg_bound_phi) > 0) {
    warning(paste(warning_msg_bound_phi, collapse = "\n  "), call. = FALSE)
  }

  cli::cli_progress_done(id = pb_id)

  # Posterior draws of the linkage matrix
  gamma_avg <- apply(gamma.iter, 2, function(x) mean(x[StEM_burnin:StEM_iter], drop = FALSE))
  eta_avg   <- lapply(eta.iter,     function(x) apply(x[StEM_burnin:StEM_iter, , drop = FALSE], 2, mean))
  alpha_avg <- lapply(alpha.iter,   function(x) apply(x[StEM_burnin:StEM_iter, , drop = FALSE], 2, mean))
  phi_avg   <- lapply(phi.iter,     function(x) apply(x[StEM_burnin:StEM_iter, , drop = FALSE], 2, mean))

  Delta <- Matrix::Matrix(0, nrow = nrow(encodedA), ncol = nrow(encodedB), sparse = TRUE)
  pbfinal_id <- cli::cli_progress_bar(
    format = "Drawing Delta          {cli::pb_bar} {cli::pb_percent} [{cli::pb_elapsed}]",
    total = n_post_sample, clear = FALSE
  )
  for (m in seq_len(n_post_sample)) {
    step <- gibbs_step(linksCpp, linksR, sumRowD, sumColD, survivalpSameH, gamma_avg, eta_avg, alpha_avg, phi_avg)
    # Note: eta/alpha/phi are frozen at their posterior mean for these draws
    linksCpp <- step$linksCpp
    linksR <- step$linksR
    sumRowD <- step$sumRowD
    sumColD <- step$sumColD
    survivalpSameH <- step$survivalpSameH

    if (nrow(linksR) > 0) {
      for (l in seq_len(nrow(linksR))) Delta[linksR[l, 1], linksR[l, 2]] <- Delta[linksR[l, 1], linksR[l, 2]] + 1
    }
    cli::cli_progress_update(id = pbfinal_id)
  }
  Delta <- Delta / n_post_sample
  cli::cli_progress_done(id = pbfinal_id)

  if (music_on) utils::browseURL("https://www.youtube.com/watch?v=NTa6Xbzfq1U")

  list(
    Delta = as.data.frame(Matrix::summary(Delta)),
    gamma = gamma.iter, eta = eta.iter, alpha = alpha.iter, phi = phi.iter,
    model_dynamics = model_dynamics
  )
}

#' Naive (deterministic exact-match) record linkage
#'
#' Links records that agree exactly on every non-missing PIV. Does not
#' enforce the one-to-one assignment constraint, so should be used only to
#' gauge the difficulty of the linkage task (amount of duplication, and
#' discriminative power of the PIVs together).
#'
#' @param PIVs Character vector, names of the PIVs (columns present in both
#'   files).
#' @param encodedA The data source (PIVs encoded to natural numbers).
#' @param encodedB The data source (PIVs encoded to natural numbers).
#' @param na_match Logical; if `TRUE`, a missing PIV value is treated as
#'   matching any value in the other file (default `TRUE`).
#' @param na_is_zero Logical; if `TRUE`, missing values are already coded as
#'   `0` (FlexRL's convention); if `FALSE`, `NA` will be recoded to `0`
#'   (default `TRUE`).
#'
#' @return Data frame with columns `idxA`, `idxB`: pairs of record indices
#'   that match.
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V2 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V3 = list(dynamics = "flexible",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(NA,NA)),
#'                      V4 = list(dynamics = "structured",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(0.03,0.03),
#'                                cond_hazard_cov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())) )
#' n_values  <- c( 5, 6, 7, 12 )
#' p_mistake <- list( V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                    V3 = c(0.05, 0.05), V4 = c(0.02, 0.02) )
#' p_missing <- list( V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                    V3 = c(0.005, 0.005), V4 = c(0.005, 0.005) )
#' cond_hazard_params <- list( V1 = c(), V2 = c(), 
#'                              V3 = c(), V4 = log(c(0.7, 0.6, 0.5)) )
#' gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
#'                            p_mistake, p_missing, cond_hazard_params, TRUE )
#' prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
#'                            PIVs_config, TRUE, "entity_id", TRUE )
#' naive_linkage( names(prep_data$PIVs_config), 
#'                prep_data$encodedA, prep_data$encodedB )
naive_linkage <- function(PIVs, encodedA, encodedB, na_match = TRUE, na_is_zero = TRUE) {

  if (!na_is_zero) {
    encodedA[PIVs][is.na(encodedA[PIVs])] <- 0
    encodedB[PIVs][is.na(encodedB[PIVs])] <- 0
  }
  if (!isTRUE(na_match) && !isFALSE(na_match)) {
    stop("`na_match` must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(na_match) && length(PIVs)==1){
    na_match = FALSE
    warning("`na_match` is set to FALSE with one PIV.", call. = FALSE)
  }

  rownames(encodedA) <- seq_len(nrow(encodedA))
  rownames(encodedB) <- seq_len(nrow(encodedB))
  DeltaNaiveLinked <- data.frame(idxA = character(0), idxB = character(0))

  match_block <- function(idA, idB) {
    valuesU <- unique(c(idA$U, idB$U))
    a <- as.numeric(factor(idA$U, levels = valuesU))
    b <- as.numeric(factor(idB$U, levels = valuesU))
    tmpA <- indexPatterns(a, length(valuesU))
    tmpB <- indexPatterns(b, length(valuesU))
    select <- pairPatterns(tmpA, tmpB, length(tmpA))
    if (nrow(select) == 0) return(NULL)
    data.frame(idxA = idA$ID[as.integer(select[, 1])], idxB = idB$ID[as.integer(select[, 2])])
  }

  # A (non-missing) vs. B (exact match on all PIVs), and vice versa
  isNotMissingA <- apply(encodedA[, PIVs, drop=FALSE] != 0, 1, all)
  isNotMissingB <- apply(encodedB[, PIVs, drop=FALSE] != 0, 1, all)
  DeltaNaiveLinked <- rbind(
    DeltaNaiveLinked,
    match_block(
      list(U = pasteIntoPattern(as.matrix(encodedA[isNotMissingA, PIVs, drop=FALSE])), ID = rownames(encodedA)[isNotMissingA]),
      list(U = pasteIntoPattern(as.matrix(encodedB[, PIVs, drop=FALSE])), ID = rownames(encodedB))
    ),
    match_block(
      list(U = pasteIntoPattern(as.matrix(encodedA[, PIVs, drop=FALSE])), ID = rownames(encodedA)),
      list(U = pasteIntoPattern(as.matrix(encodedB[isNotMissingB, PIVs, drop=FALSE])), ID = rownames(encodedB)[isNotMissingB])
    )
  )

  # Records with exactly one missing PIV: match on the remaining PIVs
  # (only run when na_match = TRUE; missing values otherwise never count
  # as matches), (there is no match when several values are missing)
  if (isTRUE(na_match)) {
    for (k in seq_along(PIVs)) {
      isMissingA_k <- encodedA[, PIVs, drop=FALSE][, k] == 0
      isMissingB_k <- encodedB[, PIVs, drop=FALSE][, k] == 0
      DeltaNaiveLinked <- rbind(
        DeltaNaiveLinked,
        match_block(
          list(U = pasteIntoPattern(as.matrix(encodedA[isMissingA_k, PIVs, drop=FALSE][, -k])), ID = rownames(encodedA)[isMissingA_k]),
          list(U = pasteIntoPattern(as.matrix(encodedB[, PIVs, drop=FALSE][, -k])), ID = rownames(encodedB))
        ),
        match_block(
          list(U = pasteIntoPattern(as.matrix(encodedA[, PIVs, drop=FALSE][, -k])), ID = rownames(encodedA)),
          list(U = pasteIntoPattern(as.matrix(encodedB[isMissingB_k, PIVs, drop=FALSE][, -k])), ID = rownames(encodedB)[isMissingB_k])
        )
      )
    }
  }

  DeltaNaiveLinked[!duplicated(DeltaNaiveLinked), ]
}

# ============================================================================
# Wrappers around external record linkage packages.
#
# Every link_with_* wrapper takes two data sets (`dataA`, `dataB`, e.g. as
# returned by synthesise(), or the original encodedA/encodedB) and a list
# `arguments` forwarded to the underlying package, runs that package, and
# returns list(idxA, idxB, linked_scores). The scores are
# `NULL` for packages that do not return a per-pair score. This uniform
# output is what compute_RL_FDP_score() and compute_augmRL_FDP_synth() use.
# FDP estimation methods: https://doi.org/10.1002/sim.70292.
# ============================================================================

#' Run FlexRL through the uniform wrapper interface
#'
#' Links `dataA` and `dataB` with the package and returns the linked pairs in the
#' common format used by [compute_RL_FDP_score()] and
#' [compute_augmRL_FDP_synth()]. `dataA` and `dataB` can be the outputs of
#' [synthesise()] or the original encoded data sources.
#'
#' More details about FlexRL on
#' https://cran.r-project.org/web/packages/FlexRL/index.html
#'
#' @param dataA Data source to link: the encoded source from
#'   [prepare_data()], or the augmented sources from [synthesise()].
#' @param dataB Data source to link: the encoded source from
#'   [prepare_data()], or the augmented sources from [synthesise()].
#' @param arguments List of extra arguments forwarded to [StEM()] (e.g.
#'   `data`, `StEM_iter`, `StEM_burnin`, `gibbs_iter`, `gibbs_burnin`).
#' @param ... Ignored; kept for compatibility with [RL_diagnostics()].
#' 
#' @return Named list with `idxA`, `idxB` and `LinkScore` (row indices 
#'   in `dataA`, `dataB` and linkage scores of potential linked records).
#' @export
#'
link_with_FlexRL <- function(dataA, dataB, arguments, ...) {
  PIVs <- names(arguments$data$PIVs_config)
  dataA[PIVs][is.na(dataA[PIVs])] <- 0 # FlexRL sentinel for missing values
  dataB[PIVs][is.na(dataB[PIVs])] <- 0
  arguments$data[["encodedA"]] <- dataA
  arguments$data[["encodedB"]] <- dataB
  fit <- do.call(StEM, arguments[intersect(names(arguments), names(formals(StEM)))])
  list(idxA = fit$Delta$i, idxB = fit$Delta$j, LinkScore = fit$Delta$x)
}

#' Run fedmatch through the uniform wrapper interface
#'
#' Links `dataA` and `dataB` with the package and returns the linked pairs in the
#' common format used by [compute_RL_FDP_score()] and
#' [compute_augmRL_FDP_synth()]. `dataA` and `dataB` can be the outputs of
#' [synthesise()] or the original encoded data sources.
#'
#' More details about fedmatch on
#' https://cran.r-project.org/web/packages/fedmatch/index.html
#'
#' @param dataA Data source to link: the encoded source from
#'   [prepare_data()], or the augmented sources from [synthesise()].
#' @param dataB Data source to link: the encoded source from
#'   [prepare_data()], or the augmented sources from [synthesise()].
#' @param arguments List of extra arguments forwarded to
#'   `fedmatch::merge_plus()` (e.g. `by`, `match_type`, `unique_key_1`,
#'   `unique_key_2`, `multivar_settings`).
#' @param ... Ignored; kept for compatibility with [RL_diagnostics()].
#' 
#' @return Named list with `idxA`, `idxB` and `LinkScore` (row indices 
#'   in `dataA`, `dataB` and linkage scores of potential linked records).
#' @export
#'
link_with_fedmatch <- function(dataA, dataB, arguments, ...) {
  dataA[, arguments$unique_key_1] <- rownames(dataA)
  dataB[, arguments$unique_key_2] <- rownames(dataB)
  arguments$data1 <- dataA
  arguments$data2 <- dataB
  results <- do.call(fedmatch::merge_plus, arguments[intersect(names(arguments), names(formals(fedmatch::merge_plus)))])
  list(idxA = results$matches$unique_key_A, idxB = results$matches$unique_key_B,
       LinkScore = results$matches$multivar_score)
}

#' Run reclin2 through the uniform wrapper interface
#'
#' Links `dataA` and `dataB` with the package and returns the linked pairs in the
#' common format used by [compute_RL_FDP_score()] and
#' [compute_augmRL_FDP_synth()]. `dataA` and `dataB` can be the outputs of
#' [synthesise()] or the original encoded data sources.
#'
#' More details about reclin2 on
#' https://cran.r-project.org/web/packages/reclin2/index.html
#'
#' @param dataA Data source to link: the encoded source from
#'   [prepare_data()], or the augmented sources from [synthesise()].
#' @param dataB Data source to link: the encoded source from
#'   [prepare_data()], or the augmented sources from [synthesise()].
#' @param arguments List of extra arguments forwarded across the `reclin2`
#'   pipeline (`pair()`, `compare_pairs()`, `problink_em()`, `predict()`,
#'   `select_threshold()`): e.g. `on`, `formula`, `type`, `add`, `variable`,
#'   `score`, `threshold`.
#' @param ... Ignored; kept for compatibility with [RL_diagnostics()].
#' 
#' @return Named list with `idxA`, `idxB` and `LinkScore` (row indices 
#'   in `dataA`, `dataB` and linkage scores of potential linked records
#'   for the `select_threshold()` argument).
#' @export
#'
link_with_reclin2 <- function(dataA, dataB, arguments, ...) {
  arguments$x <- dataA
  arguments$y <- dataB
  pairs <- do.call(reclin2::pair, arguments[intersect(names(arguments), names(formals(reclin2::pair)))])

  arguments$pairs <- pairs
  pairs <- do.call(reclin2::compare_pairs, arguments[intersect(names(arguments), names(formals(reclin2::compare_pairs)))])

  arguments$data <- pairs
  model <- do.call(reclin2::problink_em, arguments[intersect(names(arguments), names(formals(reclin2::problink_em)))])

  arguments$object <- model
  arguments$pairs <- pairs
  predict_formals <- union(names(formals(stats::predict)), c("object", "pairs", "newdata", "type", "binary", "add", "comparators", "inplace", "new_name"))
  pairs <- do.call(stats::predict, arguments[intersect(names(arguments), predict_formals)])

  arguments$pairs <- pairs
  selected <- do.call(reclin2::select_threshold, arguments[intersect(names(arguments), names(formals(reclin2::select_threshold)))])
  selected <- as.data.frame(selected)
  selected <- selected[selected[, arguments$variable] == TRUE, ]

  list(idxA = selected$.x, idxB = selected$.y, LinkScore = selected[, arguments$type])
}

#' Run BRL through the uniform wrapper interface
#'
#' Links `dataA` and `dataB` with the package and returns the linked pairs in the
#' common format used by [compute_RL_FDP_score()] and
#' [compute_augmRL_FDP_synth()]. `dataA` and `dataB` can be the outputs of
#' [synthesise()] or the original encoded data sources.
#'
#' More details about BRL on
#' https://cran.r-project.org/web/packages/BRL/index.html
#'
#' @param dataA Data source to link: the encoded source from
#'   [prepare_data()], or the augmented sources from [synthesise()].
#' @param dataB Data source to link: the encoded source from
#'   [prepare_data()], or the augmented sources from [synthesise()].
#' @param arguments List of extra arguments forwarded to `BRL::compareRecords()`
#'   and `BRL::bipartiteGibbs()` (e.g. `flds`, `types`, `nIter`).
#' @param ... Ignored; kept for compatibility with [RL_diagnostics()].
#' 
#' @return Named list with `idxA`, `idxB` and `LinkScore` (row indices 
#'   in `dataA`, `dataB` and linkage scores of potential linked records).
#' @export
#'
link_with_BRL <- function(dataA, dataB, arguments, ...) {

  arguments$df1 <- dataB
  arguments$df2 <- dataA
  comp_data <- do.call(BRL::compareRecords, arguments[intersect(names(arguments), names(formals(BRL::compareRecords)))])

  arguments$cd <- comp_data
  chain <- do.call(BRL::bipartiteGibbs, arguments[intersect(names(arguments), names(formals(BRL::bipartiteGibbs)))])

  if (!"nIter" %in% names(arguments)) arguments$nIter <- 1000
  n1 <- nrow(dataB)
  n2 <- nrow(dataA)
  z_chain <- chain$Z[, setdiff(seq_len(arguments$nIter), seq_len(0.10 * arguments$nIter)), drop = FALSE]
  z_chain[z_chain > n1 + 1] <- n1 + 1

  table_labels <- apply(z_chain, 1, tabulate, nbins = n1 + 1) / ncol(z_chain)
  prob_no_link <- table_labels[n1 + 1, ]
  max_prob_option <- apply(table_labels, 2, which.max)
  prob_max_prob_option <- apply(table_labels, 2, max)
  max_prob_option_is_link <- max_prob_option <= n1
  
  # Bayes-optimal decision rule under a symmetric loss
  # (see BRL package documentation)
  lFM1 <- 1
  lFM2 <- 2
  lFNM <- 1
  # threshold is computed as: lFM1 / (lFM1 + lFNM) + (lFM2 - lFM1 - lFNM) * (1 - prob_no_link - prob_max_prob_option) / (lFM1 + lFNM)
  # but we return all pairs (i.e. threshold = 0)
  isLink <- max_prob_option_is_link & (prob_max_prob_option > 0)
  z_hat <- (n1+1):(n1+n2)
  z_hat[isLink] <- max_prob_option[isLink]

  list(idxA = which(z_hat <= n1), idxB = z_hat[z_hat <= n1], LinkScore = prob_max_prob_option[isLink])
}

#' Run fastLink through the uniform wrapper interface
#'
#' Links `dataA` and `dataB` with the package and returns the linked pairs in the
#' common format used by [compute_RL_FDP_score()] and
#' [compute_augmRL_FDP_synth()]. `dataA` and `dataB` can be the outputs of
#' [synthesise()] or the original encoded data sources.
#'
#' More details about fastLink on
#' https://cran.r-project.org/web/packages/fastLink/index.html
#'
#' @param dataA Data source to link: the encoded source from
#'   [prepare_data()], or the augmented sources from [synthesise()].
#' @param dataB Data source to link: the encoded source from
#'   [prepare_data()], or the augmented sources from [synthesise()].
#' @param arguments List of extra arguments forwarded to `fastLink::fastLink()`
#'   (e.g. `varnames`, `threshold.match`, `tol.em`, `return.all`).
#' @param ... Ignored; kept for compatibility with [RL_diagnostics()].
#' 
#' @return Named list with `idxA`, `idxB` and `LinkScore` (row indices 
#'   in `dataA`, `dataB` and linkage scores of potential linked records).
#' @export
#'
link_with_fastLink <- function(dataA, dataB, arguments, ...) {
  arguments$dfA <- dataA
  arguments$dfB <- dataB
  out <- do.call(fastLink::fastLink, arguments[intersect(names(arguments), names(formals(fastLink::fastLink)))])
  list(idxA = out$matches$inds.a, idxB = out$matches$inds.b, LinkScore = out$posterior)
}

#' Run multilink through the uniform wrapper interface
#'
#' Links `dataA` and `dataB` with the package and returns the linked pairs in the
#' common format used by [compute_RL_FDP_score()] and
#' [compute_augmRL_FDP_synth()]. `dataA` and `dataB` can be the outputs of
#' [synthesise()] or the original encoded data sources.
#'
#' More details about multilink on
#' https://cran.r-project.org/web/packages/multilink/index.html
#'
#' @param dataA Data source to link: the encoded source from
#'   [prepare_data()], or the augmented sources from [synthesise()].
#' @param dataB Data source to link: the encoded source from
#'   [prepare_data()], or the augmented sources from [synthesise()].
#' @param arguments List of extra arguments forwarded across the `multilink`
#'   pipeline (`create_comparison_data()`, `reduce_comparison_data()`,
#'   `specify_prior()`, `gibbs_sampler()`, `find_bayes_estimate()`,
#'   `relabel_bayes_estimate()`): e.g. `records`, `types`, `breaks`, 
#'   `duplicates`, `n_iter`.
#' @param ... Ignored; kept for compatibility with [RL_diagnostics()].
#' 
#' @return Named list with `idxA`, `idxB` and `LinkScore = NULL` (row indices 
#'   of default potential linked records in `dataA`, `dataB`, this package 
#'   does not return scores).
#' @export
#'
link_with_multilink <- function(dataA, dataB, arguments, ...) {
  PIVs <- names(arguments$breaks)
  arguments$records <- as.data.frame(lapply(rbind(dataA[, PIVs, drop=FALSE], dataB[, PIVs, drop=FALSE]), as.character))
  arguments$file_sizes <- c(nrow(dataA), nrow(dataB))
  comparison_list <- do.call(multilink::create_comparison_data, arguments[intersect(names(arguments), names(formals(multilink::create_comparison_data)))])

  arguments$comparison_list <- comparison_list
  arguments$pairs_to_keep <- rep(TRUE, nrow(dataA) * nrow(dataB))
  reduced <- do.call(multilink::reduce_comparison_data, arguments[intersect(names(arguments), names(formals(multilink::reduce_comparison_data)))])

  if (!"dup_count_prior_family" %in% names(arguments)) arguments$dup_count_prior_family <- rep(NA, reduced$K)
  if (!"dup_count_prior_pars" %in% names(arguments)) arguments$dup_count_prior_pars <- rep(NA, reduced$K)
  prior <- do.call(multilink::specify_prior, arguments[intersect(names(arguments), names(formals(multilink::specify_prior)))])

  arguments$comparison_list <- reduced
  arguments$prior_list <- prior
  results <- do.call(multilink::gibbs_sampler, arguments[intersect(names(arguments), names(formals(multilink::gibbs_sampler)))])

  arguments$burn_in <- if ("n_iter" %in% names(arguments)) as.integer(arguments$n_iter / 5) else 400
  arguments$partitions <- results$partitions
  estimate <- do.call(multilink::find_bayes_estimate, arguments[intersect(names(arguments), names(formals(multilink::find_bayes_estimate)))])

  arguments$reduced_comparison_list <- reduced
  arguments$bayes_estimate <- estimate
  relabel <- do.call(multilink::relabel_bayes_estimate, arguments[intersect(names(arguments), names(formals(multilink::relabel_bayes_estimate)))])

  linked_pairs <- data.frame(idxA = integer(0), idxB = integer(0))
  for (i in seq_along(relabel$link_id)) {
    l <- relabel$link_id[i]
    if (l > 0 && sum(relabel$link_id == l) > 1) {
      idx <- which(relabel$link_id == l)
      linked_pairs <- rbind(linked_pairs, data.frame(idxA = idx[1], idxB = idx[2] - nrow(dataA)))
    }
  }
  linked_pairs <- unique(linked_pairs)
  list(idxA = linked_pairs$idxA, idxB = linked_pairs$idxB, LinkScore = NULL)
}

#' Run diyar through the uniform wrapper interface
#'
#' Links `dataA` and `dataB` with the package and returns the linked pairs in the
#' common format used by [compute_RL_FDP_score()] and
#' [compute_augmRL_FDP_synth()]. `dataA` and `dataB` can be the outputs of
#' [synthesise()] or the original encoded data sources.
#'
#' More details about diyar on
#' https://cran.r-project.org/web/packages/diyar/index.html
#'
#' @param dataA Data source to link: the encoded source from
#'   [prepare_data()], or the augmented sources from [synthesise()].
#' @param dataB Data source to link: the encoded source from
#'   [prepare_data()], or the augmented sources from [synthesise()].
#' @param arguments List of extra arguments forwarded to
#'   `diyar::prob_score_range()` and `diyar::links_wf_probabilistic()` (e.g.
#'   `attribute`, `probabilistic`, `return_weights`).
#' @param ... Ignored; kept for compatibility with [RL_diagnostics()].
#' 
#' @return Named list with `idxA`, `idxB` and `LinkScore = NULL` (row indices 
#'   of default potential linked records in `dataA`, `dataB`, this package 
#'   does not return scores).
#' @export
#'
link_with_diyar <- function(dataA, dataB, arguments, ...) {
  all_records <- rbind(dataA, dataB)
  PIVs <- setdiff(names(all_records),c("local_id", "source"))
  arguments$attribute <- as.list(all_records[, PIVs, drop=FALSE])
  link_scores <- do.call(diyar::prob_score_range, arguments[intersect(names(arguments), names(formals(diyar::prob_score_range)))])

  arguments$score_threshold <- link_scores$mid_scorce
  arguments$data_source <- all_records$source
  id_linkage <- do.call(diyar::links_wf_probabilistic, arguments[intersect(names(arguments), names(formals(diyar::links_wf_probabilistic)))])

  linked_pairs <- data.frame(idxA = integer(0), idxB = integer(0))
  for (clusterid in unique(id_linkage$pid)) {
    idx <- which(id_linkage$pid == clusterid)
    if (length(idx) == 2 && idx[1] <= nrow(dataA) && idx[2] > nrow(dataA)) {
      linked_pairs <- rbind(linked_pairs, data.frame(idxA = idx[1], idxB = idx[2] - nrow(dataA)))
    }
  }
  linked_pairs <- unique(linked_pairs)
  list(idxA = linked_pairs$idxA, idxB = linked_pairs$idxB, LinkScore = NULL)
}

#' Agreement rate of pairs on shared variables
#'
#' For a set of pairs, computes how often the two records agree on each variable 
#' in `common_vars`.
#'
#' @param data1 Data frame containing `common_vars`.
#' @param data2 Data frame containing `common_vars`.
#' @param common_vars Character vector, names of the variables to compare
#'   (must exist in both `data1` and `data2`).
#' @param pairs Data frame/matrix/list with 2 columns of indices
#'   (into `data1`, `data2`) for the pairs to evaluate.
#' @param na.rm Logical; if `TRUE` (default), pairs with a missing value on a
#'   variable are excluded from that variable's agreement rate.
#' @param na_match Logical, required if `na.rm = FALSE`: should a missing value
#'   be treated as agreeing (`TRUE`) or disagreeing (`FALSE`) with any value?
#'
#' @return List with `agreements` (named numeric vector, one entry per
#'   variable in `common_vars`) and, if available, `true_agreements`.
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V2 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V3 = list(dynamics = "flexible",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(NA,NA)),
#'                      V4 = list(dynamics = "structured",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(0.03,0.03),
#'                                cond_hazard_cov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())) )
#' n_values  <- c( 5, 6, 7, 12 )
#' p_mistake <- list( V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                    V3 = c(0.05, 0.05), V4 = c(0.02, 0.02) )
#' p_missing <- list( V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                    V3 = c(0.005, 0.005), V4 = c(0.005, 0.005) )
#' cond_hazard_params <- list( V1 = c(), V2 = c(), 
#'                              V3 = c(), V4 = log(c(0.7, 0.6, 0.5)) )
#' gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
#'                            p_mistake, p_missing, cond_hazard_params, TRUE )
#' prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
#'                            PIVs_config, TRUE, "entity_id", TRUE )
#' fit <- StEM( data = prep_data, StEM_iter = 10, StEM_burnin = 5, 
#'              gibbs_iter = 10, gibbs_burnin = 5, n_post_sample = 10 )                         
#' linked_pairs <- fit$Delta[fit$Delta$x > 0.5, ]
#' RL_agreement( prep_data$encodedA, prep_data$encodedB,
#'               names(PIVs_config), linked_pairs )
#' RL_agreement( prep_data$encodedA, prep_data$encodedB,
#'               names(PIVs_config), prep_data$true_pairs )
RL_agreement <- function(data1, data2, common_vars, pairs,
                         na.rm = TRUE, na_match = NULL) {
  stopifnot(is.data.frame(data1), is.data.frame(data2))
  if (!is.character(common_vars) || length(common_vars) == 0) {
    stop("`common_vars` must be a non-empty character vector.", call. = FALSE)
  }
  missing1 <- setdiff(common_vars, names(data1))
  missing2 <- setdiff(common_vars, names(data2))
  if (length(missing1) || length(missing2)) {
    stop("`common_vars` missing from data1: [", paste(missing1, collapse = ", "),
         "]; from data2: [", paste(missing2, collapse = ", "), "]", call. = FALSE)
  }
  pairs <- as.data.frame(pairs)
  if (ncol(pairs) < 2) stop("`pairs` must have two columns.", call. = FALSE)
  idx1 <- pairs[[1]]
  idx2 <- pairs[[2]]
  agreements <- vapply(common_vars, function(v) {
    is_match <- data1[idx1, v] == data2[idx2, v]
    if (isTRUE(na.rm)) {
      is_match <- is_match[!is.na(is_match)]
    } else {
      is_match[is.na(is_match)] <- isTRUE(na_match)
    }
    if (length(is_match) == 0) NA_real_ else mean(is_match)
  }, numeric(1))
  list(agreements = agreements)
}

# Cramer's V for a two-way contingency table
.cramer_v <- function(tab) {
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  if (min(dim(tab)) < 2) return(0)
  chi2 <- suppressWarnings(stats::chisq.test(tab, correct = FALSE)$statistic)
  unname(sqrt(chi2 / (sum(tab) * (min(dim(tab)) - 1))))
}

#' Prepare two data sources for [StEM()]
#'
#' Wraps the data-preparation steps needed before calling [StEM()]: tags each
#' source with a `source` column, labels the larger source as `B`, adjusting 
#' `cond_hazard_cov`/`bound_mistakes`/`fix_mistakes` accordingly if A and 
#' B are swapped), drops records whose PIV values fall outside the support 
#' shared by both files, warns if two PIVs are strongly associated 
#' (Cramer's V > 0.3), which may degrade record linkage performance, encodes 
#' every PIV to natural numbers using levels pooled across both sources, and 
#' encodes missing values to `0`.
#'
#' @param data1 Data frame, the raw data source (whichever has more rows
#' becomes `B`).
#' @param data2 Data frame, the raw data source (whichever has more rows
#' becomes `B`).
#' @param label1 Character, label recorded in the `source` column
#'   for `data1`.
#' @param label2 Character, label recorded in the `source` column
#'   for `data2`.
#' @param PIVs_config Named list describing each PIV — see [simulate_data()].
#' @param same_mistakes Logical, will A and B share one mistake-probability
#'   parameter per PIV.
#' @param uniq_id Optional column name (present in both files) with the true
#'   entity identifier, used to build `true_pairs` for evaluation; `NULL` if
#'   unavailable.
#' @param restrict_support_intersection Logical; if `TRUE` (default), records
#'   with an out-of-common-support PIV value are dropped (and a warning issued);
#'   if `FALSE`, only the warning is issued.
#'
#' @return A list ready to use as the `data` argument of [StEM()]: `encodedA`,
#'   `encodedB`, `n_values`, `PIVs_config`, `same_mistakes`, and `true_pairs`
#'   (`NULL` if `uniq_id` was not supplied).
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V2 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V3 = list(dynamics = "flexible",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(NA,NA)),
#'                      V4 = list(dynamics = "structured",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(0.03,0.03),
#'                                cond_hazard_cov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())) )
#' n_values  <- c( 5, 6, 7, 12 )
#' p_mistake <- list( V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                    V3 = c(0.05, 0.05), V4 = c(0.02, 0.02) )
#' p_missing <- list( V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                    V3 = c(0.005, 0.005), V4 = c(0.005, 0.005) )
#' cond_hazard_params <- list( V1 = c(), V2 = c(), 
#'                              V3 = c(), V4 = log(c(0.7, 0.6, 0.5)) )
#' gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
#'                            p_mistake, p_missing, cond_hazard_params, TRUE )
#' prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
#'                            PIVs_config, TRUE, "entity_id", TRUE )
#' str(prep_data, max.level = 1)
prepare_data <- function(data1, data2, label1, label2, PIVs_config,
                         same_mistakes = TRUE, uniq_id = NULL,
                         restrict_support_intersection = TRUE) {

  rownames(data1) <- seq_len(nrow(data1))
  rownames(data2) <- seq_len(nrow(data2))

  PIVs <- names(PIVs_config)
  n_pivs <- length(PIVs)

  # PIVs must be categorical, restricted to the common support
  for (i in seq_len(n_pivs)) {
    unique1 <- unique(data1[, PIVs[i]])
    unique2 <- unique(data2[, PIVs[i]])
    common <- intersect(unique1, unique2)
    if (length(common) > 100) {
      stop("PIV '", PIVs[i], "' has more than 100 unique values in the common support; PIVs must be categorical.", call. = FALSE)
    }
    missing1 <- setdiff(unique1, common)
    missing1 <- missing1[!is.na(missing1)]
    missing2 <- setdiff(unique2, common)
    missing2 <- missing2[!is.na(missing2)]
    if (length(missing1) > 0 || length(missing2) > 0) {
      warning(sprintf(
        "PIV '%s': values out of the common support%s%s.",
        PIVs[i],
        if (length(missing1) > 0) sprintf(" in data1 [%s]", paste(missing1, collapse = ", ")) else "",
        if (length(missing2) > 0) sprintf(" in data2 [%s]", paste(missing2, collapse = ", ")) else ""
      ), call. = FALSE)
      if (restrict_support_intersection) {
        data1 <- data1[data1[, PIVs[i]] %in% c(NA, common), ]
        data2 <- data2[data2[, PIVs[i]] %in% c(NA, common), ]
      }
    }
  }
  rownames(data1) <- seq_len(nrow(data1))
  rownames(data2) <- seq_len(nrow(data2))
  
  # Ground truth (optional)
  true_Delta <- NULL
  if (!is.null(uniq_id)) {
    true_links <- intersect(data1[[uniq_id]], data2[[uniq_id]])
    true_Delta <- data.frame(matrix(0, nrow = 0, ncol = 2))
    for (id in true_links) {
      id1 <- which(data1[[uniq_id]] == id)
      id2 <- which(data2[[uniq_id]] == id)
      true_Delta <- rbind(true_Delta, cbind(rownames(data1[id1, ]), rownames(data2[id2, ])))
    }
    colnames(true_Delta) <- c(label1, label2)
  }

  # Warn about strongly associated PIVs
  # (FlexRL assumes conditional independence)
  for (i in seq_len(n_pivs)) for (j in seq_len(i - 1)) {
    v_A <- .cramer_v(table(data1[, PIVs[i]], data1[, PIVs[j]]))
    v_B <- .cramer_v(table(data2[, PIVs[i]], data2[, PIVs[j]]))
    if (v_A > 0.3 || v_B > 0.3) {
      warning(sprintf(
        "PIVs '%s' and '%s' are associated (Cramer's V = %.2f in data1, %.2f in data2); FlexRL assumes conditional independence between PIVs, consider merging or dropping one.",
        PIVs[i], PIVs[j], v_A, v_B
      ), call. = FALSE)
    }
  }

  # Validate PIVs_config 
  stopifnot(is.data.frame(data1), is.data.frame(data2))
  if (!is.list(PIVs_config) || length(PIVs_config) == 0 || is.null(names(PIVs_config))) {
    stop("`PIVs_config` must be a non-empty named list, one entry per PIV.", call. = FALSE)
  }
  allowed <- c("dynamics", "bound_mistakes", "fix_mistakes", "cond_hazard_cov")
  for (k in seq_len(n_pivs)) {
    if (!all(names(PIVs_config[[k]]) %in% allowed)) {
      stop("`PIVs_config[['", PIVs[k], "']]` may only contain: ", paste(allowed, collapse = ", "), ".", call. = FALSE)
    }
  }
  missing1 <- setdiff(PIVs, names(data1))
  missing2 <- setdiff(PIVs, names(data2))
  if (length(missing1) || length(missing2)) {
    stop("PIVs missing from data1: [", paste(missing1, collapse = ", "),
         "]; from data2: [", paste(missing2, collapse = ", "), "]", call. = FALSE)
  }

  bad_dynamics <- !vapply(PIVs_config, function(x) is.character(x$dynamics) && !is.na(x$dynamics), logical(1))
  if (any(bad_dynamics)) {
    stop("`dynamics` must be one of 'stable', 'flexible', 'structured'; problem for: ",
         paste(PIVs[bad_dynamics], collapse = ", "), ".", call. = FALSE)
  }
  is_numeric_or_na <- function(v) is.vector(v) && length(v) == 2 && all(is.na(v) | grepl("^-?[0-9.]+$", v))
  bad_fix <- !vapply(PIVs_config, function(x) is_numeric_or_na(x$fix_mistakes), logical(1))
  if (any(bad_fix)) stop("`fix_mistakes` must be a length-2 numeric/NA vector; problem for: ", paste(PIVs[bad_fix], collapse = ", "), ".", call. = FALSE)
  bad_bound <- !vapply(PIVs_config, function(x) is_numeric_or_na(x$bound_mistakes), logical(1))
  if (any(bad_bound)) stop("`bound_mistakes` must be a length-2 numeric/NA vector; problem for: ", paste(PIVs[bad_bound], collapse = ", "), ".", call. = FALSE)

  PIVs_stable <- vapply(PIVs_config, function(x) x$dynamics != "structured", logical(1))
  modelDynaPIVs <- PIVs[!PIVs_stable]

  if (length(modelDynaPIVs) > 0) {
    cov_issues <- character(0)
    for (p in modelDynaPIVs) {
      cov1 <- PIVs_config[[p]]$cond_hazard_cov$cov1
      cov2 <- PIVs_config[[p]]$cond_hazard_cov$cov2
      if (length(cov1) > 0) {
        if (!all(cov1 %in% names(data1))) cov_issues <- c(cov_issues, sprintf("PIV '%s': cov1 not found in data1: %s", p, paste(setdiff(cov1, names(data1)), collapse = ", ")))
        if (!"date" %in% names(data1)) cov_issues <- c(cov_issues, sprintf("PIV '%s' is structured but data1 has no `date` column.", p))
      }
      if (length(cov2) > 0) {
        if (!all(cov2 %in% names(data2))) cov_issues <- c(cov_issues, sprintf("PIV '%s': cov2 not found in data2: %s", p, paste(setdiff(cov2, names(data2)), collapse = ", ")))
        if (!"date" %in% names(data2)) cov_issues <- c(cov_issues, sprintf("PIV '%s' is structured but data2 has no `date` column.", p))
      }
    }
    if (length(cov_issues) > 0) stop("Invalid `cond_hazard_cov` configuration:\n  ", paste(cov_issues, collapse = "\n  "), call. = FALSE)
  }

  if (!is.logical(same_mistakes)) stop("`same_mistakes` must be logical.", call. = FALSE)
  if (same_mistakes) {
    for (k in seq_len(n_pivs)) {
      bm <- PIVs_config[[k]]$bound_mistakes
      fm <- PIVs_config[[k]]$fix_mistakes
      if (!identical(bm[1], bm[2]) || !identical(fm[1], fm[2])) {
        stop("`same_mistakes = TRUE` requires identical `bound_mistakes`/`fix_mistakes` for file 1 and file 2 (PIV '", PIVs[k], "').", call. = FALSE)
      }
    }
  }

  # For each PIV, `bound_mistakes`/`fix_mistakes` should match its `dynamics`:
  # stable     -> recommend `bound_mistakes`, discourage `fix_mistakes`
  # flexible   -> discourage both (nothing to bound/fix: change is not modelled)
  # structured -> discourage `bound_mistakes`, recommend `fix_mistakes`
  #               (to avoid confounding mistakes with genuine change over time)
  cov_issues <- c()
  for (k in seq_len(n_pivs)) {
    dyn <- PIVs_config[[k]]$dynamics
    bm  <- PIVs_config[[k]]$bound_mistakes
    fm  <- PIVs_config[[k]]$fix_mistakes

    if (dyn == "stable") {
      if (!all(is.numeric(bm))) {
        cov_issues <- c(cov_issues, sprintf("We recommend bounding the mistakes with `bound_mistakes` when `dynamics` is `stable` (PIV: %s).", PIVs[k]))
      }
      if (!all(is.na(fm))) {
        cov_issues <- c(cov_issues, sprintf("We do not recommend fixing the mistakes with `fix_mistakes` when `dynamics` is `stable` (PIV: %s).", PIVs[k]))
      }
    } else if (dyn == "flexible") {
      if (!all(is.na(bm))) {
        cov_issues <- c(cov_issues, sprintf("We do not recommend bounding the mistakes with `bound_mistakes` when `dynamics` is `flexible` (PIV: %s).", PIVs[k]))
      }
      if (!all(is.na(fm))) {
        cov_issues <- c(cov_issues, sprintf("We do not recommend fixing the mistakes with `fix_mistakes` when `dynamics` is `flexible` (PIV: %s).", PIVs[k]))
      }
    } else if (dyn == "structured") {
      if (!all(is.na(bm))) {
        cov_issues <- c(cov_issues, sprintf("We do not recommend bounding the mistakes with `bound_mistakes` when `dynamics` is `structured` (PIV: %s).", PIVs[k]))
      }
      if (!all(is.numeric(fm))) {
        cov_issues <- c(cov_issues, sprintf("We recommend fixing the mistakes with `fix_mistakes` when `dynamics` is `structured` (PIV: %s).", PIVs[k]))
      }
    }
  }
  if (length(cov_issues) > 0) warning(paste(cov_issues, collapse = "\n  "), call. = FALSE)

  # Configuration entries outside the 4 recognised fields are ignored
  useless_config <- vapply(PIVs_config, function(x) {
    any(!names(x) %in% allowed)
  }, logical(1))
  if (any(useless_config)) {
    warning(
      "Configuration fields outside of `dynamics`, `bound_mistakes`, `fix_mistakes`, `cond_hazard_cov` are ignored (PIV: ",
      paste(PIVs[useless_config], collapse = ", "), ").",
      call. = FALSE
    )
  }

  # `cond_hazard_cov` only has an effect when `dynamics = "structured"`; warn if
  # it was supplied for a stable/flexible PIV, since it will not be used.
  useless_given_cov <- vapply(PIVs_config, function(x) {
    non_empty <- any(vapply(x[["cond_hazard_cov"]], length, integer(1)) > 0)
    is_stable_or_flexible <- x$dynamics %in% c("stable", "flexible")
    non_empty && is_stable_or_flexible
  }, logical(1))
  if (any(useless_given_cov)) {
    warning(
      sprintf(
        "`cond_hazard_cov` was supplied but `dynamics` is not `structured`, so no dynamics will be modelled (PIV: %s).",
        paste(PIVs[useless_given_cov], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # source column + local_id column + assign the larger file as B
  if (!"source" %in% names(data1)) data1$source <- label1
  if (!"source" %in% names(data2)) data2$source <- label2
  if (!"local_id" %in% names(data1)) data1$local_id <- seq_len(nrow(data1))
  if (!"local_id" %in% names(data2)) data2$local_id <- seq_len(nrow(data2))

  swap <- nrow(data1) > nrow(data2)
  encodedA <- if (swap) data2 else data1
  encodedB <- if (swap) data1 else data2
  message(sprintf("'%s' is the larger source, saved as B; '%s' saved as A.", if (swap) label1 else label2, if (swap) label2 else label1))
  if (swap) {
    for (p in modelDynaPIVs) {
      names(PIVs_config[[p]]$cond_hazard_cov) <- gsub("^cov1$", "covB", gsub("^cov2$", "covA", names(PIVs_config[[p]]$cond_hazard_cov)))
    }
    for (k in seq_len(n_pivs)) {
      PIVs_config[[k]]$bound_mistakes <- rev(PIVs_config[[k]]$bound_mistakes)
      PIVs_config[[k]]$fix_mistakes <- rev(PIVs_config[[k]]$fix_mistakes)
    }
    if (!is.null(true_Delta)) {
      true_Delta <- true_Delta[c(label2, label1)]
    }
  } else {
    for (p in modelDynaPIVs) {
      names(PIVs_config[[p]]$cond_hazard_cov) <- gsub("^cov1$", "covA", gsub("^cov2$", "covB", names(PIVs_config[[p]]$cond_hazard_cov)))
    }
  }

  # Encode PIVs to natural numbers, pooling levels across A and B
  levels_PIVs <- stats::setNames(lapply(PIVs, function(x) levels(factor(as.character(c(encodedA[[x]], encodedB[[x]]))))), PIVs)
  for (x in PIVs) {
    encodedA[[x]] <- as.numeric(factor(as.character(encodedA[[x]]), levels = levels_PIVs[[x]]))
    encodedB[[x]] <- as.numeric(factor(as.character(encodedB[[x]]), levels = levels_PIVs[[x]]))
  }
  n_values <- stats::setNames(vapply(levels_PIVs, length, integer(1)), PIVs)
  encodedA[PIVs][is.na(encodedA[PIVs])] <- 0  # FlexRL sentinel for missing values
  encodedB[PIVs][is.na(encodedB[PIVs])] <- 0

  list(
    encodedA = encodedA, encodedB = encodedB, n_values = n_values,
    PIVs_config = PIVs_config, same_mistakes = same_mistakes, true_pairs = true_Delta
  )
}

#' Compare distributions of shared variables across data sets
#' (histograms/barplots)
#'
#' Overlays, for each variable in `common_vars`, the empirical distribution
#' in every data set of `data_list` (e.g. a baseline file vs. the linked set).
#'
#' @param data_list Named list of data frames, each containing `common_vars`.
#' @param common_vars Character vector, variables to compare.
#' @param threshold Numeric, the linkage decision rule that defined the linked
#'   data; only used to label the linked data in the legend.
#' @param colours Optional vector of colours, one per element of `data_list`.
#'
#' @return `NULL`, invisibly; called for its plotting side effect.
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V2 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V3 = list(dynamics = "flexible",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(NA,NA)),
#'                      V4 = list(dynamics = "structured",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(0.03,0.03),
#'                                cond_hazard_cov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())) )
#' n_values  <- c( 5, 6, 7, 12 )
#' p_mistake <- list( V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                    V3 = c(0.05, 0.05), V4 = c(0.02, 0.02) )
#' p_missing <- list( V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                    V3 = c(0.005, 0.005), V4 = c(0.005, 0.005) )
#' cond_hazard_params <- list( V1 = c(), V2 = c(), 
#'                              V3 = c(), V4 = log(c(0.7, 0.6, 0.5)) )
#' gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
#'                            p_mistake, p_missing, cond_hazard_params, TRUE )
#' prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
#'                            PIVs_config, TRUE, "entity_id", TRUE )
#' fit <- StEM( data = prep_data, StEM_iter = 10, StEM_burnin = 5,
#'              gibbs_iter = 10, gibbs_burnin = 5, n_post_sample = 10 )
#' threshold_strict <- stats::quantile(fit$Delta$x, 0.75)
#' data_list = list(data_baseline = prep_data$encodedA, 
#' data_select = prep_data$encodedA[fit$Delta[fit$Delta$x > threshold_strict, "i"],])
#' common_vars = names(PIVs_config)
#' plot_distributions(data_list, common_vars)
#' common_vars = c("Xe", "Xf")
#' plot_distributions(data_list, common_vars)
plot_distributions <- function(data_list, common_vars, threshold, colours = NULL) {
  if (!is.character(common_vars) || length(common_vars) == 0) {
    stop("`common_vars` must be a non-empty character vector.", call. = FALSE)
  }
  for (data in data_list) {
    stopifnot(is.data.frame(data))
    missing <- setdiff(common_vars, names(data))
    if (length(missing)) stop("`common_vars` not found in a data set: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  N <- length(data_list)
  if (is.null(colours)) colours <- grDevices::rgb(stats::runif(N), stats::runif(N), stats::runif(N), alpha = 0.6)

  for (v in common_vars) {
    if (length(unique(data_list[[1]][, v])) > 20 && is.numeric(data_list[[1]][, v])) {
      b <- diff(range(unique(unlist(data_list[[1]][, v])), na.rm = TRUE)) + 2
      h <- graphics::hist(as.numeric(data_list[[1]][, v]), ylim = c(0, 1), breaks = b, col = colours[1],
                          main = sprintf("Distribution of %s", v), prob = TRUE, xlab = v, ylab = "Density")
      for (d in 2:N) {
        graphics::hist(as.numeric(data_list[[d]][, v]), ylim = c(0, 1), col = colours[d], breaks = h$breaks, prob = TRUE, add = TRUE)
      }
    } else {
      levels_v <- sort(unique(unlist(lapply(data_list, function(d) unique(d[[v]])))))
      dens <- prop.table(table(factor(data_list[[1]][[v]], levels_v)))
      graphics::barplot(dens, ylim = c(0, 1), main = sprintf("Distribution of %s", v), col = colours[1], xlab = v, ylab = "Density")
      for (d in 2:N) {
        graphics::barplot(prop.table(table(factor(data_list[[d]][[v]], levels_v))), ylim = c(0, 1), col = colours[d], add = TRUE)
      }
    }
    names_data_list <- c()
    for (n in names(data_list)){
      if (length(grep("linked", n))>0){
        new_name <- sprintf("%s (%s)", n, threshold)
        names_data_list <- c(names_data_list, new_name)
      } else {
        names_data_list <- c(names_data_list, n)
      }
    }
    graphics::legend("topright", paste(names_data_list, sapply(data_list, nrow), sep = ": obs. "), col = colours, lwd = 10)
  }
  invisible(NULL)
}

#' Maximum Mean Discrepancy between two sets of variables
#'
#' Joint (multivariate) measure of distributional discrepancy between `x`
#' and `y`, using a Gaussian (RBF) kernel with bandwidth set by the median
#' pairwise distance (biased estimator, Gretton et al. 2012).
#'
#' @param x,y Numeric matrices with the same number of columns.
#'
#' @return Numeric, the (biased) MMD estimate.
#' @export
#'
#' @examples
#' x <- matrix(rnorm(50), ncol = 5)
#' y <- matrix(rnorm(30) + 0.5, ncol = 5)
#' mmd(x, y)
mmd <- function(x, y) {
  stopifnot(ncol(x) == ncol(y))
  n <- nrow(x)
  m <- nrow(y)
  dists <- as.matrix(stats::dist(rbind(x, y)))
  sigma <- stats::median(dists) / 2
  k <- exp(-(dists^2) / (2 * sigma^2)) + diag(1e-5, n + m)

  k_x <- k[1:n, 1:n]
  k_y <- k[(n + 1):(n + m), (n + 1):(n + m)]
  k_xy <- k[1:n, (n + 1):(n + m)]

  sum(k_x) / (n * (n - 1)) + sum(k_y) / (m * (m - 1)) - 2 * sum(k_xy) / (n * m)
}


#' Intersection-over-union of two histogram supports
#'
#' @param h1 `histogram` object (as returned by [graphics::hist()]).
#' @param h2 `histogram` object (as returned by [graphics::hist()]).
#'
#' @return Numeric, IoU of the ranges over which `h1` and `h2`
#'   have non-zero counts.
#' @export
#'
#' @examples
#' h1 <- graphics::hist(rnorm(200), plot = FALSE)
#' h2 <- graphics::hist(rnorm(200) + 1, plot = FALSE)
#' support_iou(h1, h2)
support_iou <- function(h1, h2) {
  r1 <- range(h1$breaks[h1$counts > 0])
  r2 <- range(h2$breaks[h2$counts > 0])
  intersection <- max(0, min(r1[2], r2[2]) - max(r1[1], r2[1]))
  union <- max(r1[2], r2[2]) - min(r1[1], r2[1])
  intersection / union
}

#' Proportion of a data set with a given variable at a given level
#'
#' @param df Data frame.
#' @param var Character, column name in `df`.
#' @param level Value to match against `df[[var]]`.
#'
#' @return Numeric, proportion of rows where `df[[var]] == level` (na.rm=TRUE).
#' @export
#'
#' @examples
#' df <- data.frame(colour = sample(c("orange", "purple"), 100, replace = TRUE))
#' prop_level(df, "colour", "purple")
prop_level <- function(df, var, level) {
  mean(df[, var] == level, na.rm=TRUE)
}

#' Augment file B with synthetic records
#'
#' Fits a generative model on `encodedB[, PIVs]` and draws `n_synth`
#' new synthetic records from it, appended to `encodedB` with `source =
#' "synthetic"`; used by [compute_augmRL_FDP_synth()] to estimate the false
#' discovery proportion of a record-linkage method without ground truth.
#'
#' @param method One of `"arf"` ([arf::adversarial_rf()]), `"synthpop"`
#'   ([synthpop::syn()]), or `"mice"` ([mice::mice()]).
#' @param encodedA The (already prepared / encoded) data source.
#' @param encodedB The (already prepared / encoded) data source.
#' @param PIVs Character vector, names of the PIVs to synthesise.
#' @param n_synth Integer, number of synthetic records to generate.
#' @param restrict_support_intersection Logical; drop synthetic records whose
#'   PIV values fall outside the support shared with `encodedA` (default `TRUE`).
#'
#' @return List with `dataA` (unchanged `encodedA`, support-restricted) and
#'   `dataB` (`encodedB` plus the synthetic records).
#' @export
#'
synthesise <- function(method, encodedA, encodedB, PIVs, n_synth, restrict_support_intersection = TRUE) {

  if (!method %in% c("arf", "synthpop", "mice")) {
    stop("`method` must be one of 'arf', 'synthpop', 'mice'.", call. = FALSE)
  }

  if (method == "arf") {
    arf_model <- arf::adversarial_rf(encodedB[, PIVs, drop=FALSE])
    psi <- arf::forde(arf_model, encodedB[, PIVs, drop=FALSE])
    syntheticNewB <- arf::forge(psi, n_synth)
  } else {
    # synthpop / mice work better with continuous coding for
    # very-high-cardinality PIVs
    for (p in PIVs) {
      if (length(unique(encodedB[, p])) >= 60) encodedB[, p] <- as.numeric(encodedB[, p])
      if (length(unique(encodedA[, p])) >= 60) encodedA[, p] <- as.numeric(encodedA[, p])
    }
    if (method == "synthpop") {
      syntheticNewB <- synthpop::syn(encodedB[, PIVs, drop=FALSE], k = n_synth)$syn
    } else {
      empty <- matrix(NA, nrow = n_synth, ncol = length(PIVs), dimnames = list(NULL, PIVs))
      imputed <- mice::complete(mice::mice(rbind(encodedB[, PIVs, drop=FALSE], empty), m = 1))
      syntheticNewB <- imputed[(nrow(encodedB) + 1):nrow(imputed), ]
    }
  }
  rownames(syntheticNewB) <- seq_len(nrow(syntheticNewB))

  extraCols <- setdiff(names(encodedB), c(PIVs, "local_id", "source"))
  for (col in extraCols) syntheticNewB[, col] <- NA

  syntheticNewB$local_id <- nrow(encodedB) + seq_len(nrow(syntheticNewB))
  syntheticNewB$source <- "synthetic"
  encodedNewB <- rbind(encodedB, syntheticNewB)

  for (p in PIVs) {
    encodedNewB[, p] <- as.integer(encodedNewB[, p])
    encodedA[, p] <- as.integer(encodedA[, p])
  }

  for (p in PIVs) {
    common <- intersect(unique(encodedA[, p]), unique(encodedNewB[, p]))
    missing1 <- setdiff(unique(encodedA[, p]), common)
    missing2 <- setdiff(unique(encodedNewB[, p]), common)
    if (length(missing1) > 0 || length(missing2) > 0) {
      warning(sprintf("PIV '%s': synthetic records out of the common support were %s.",
                      p, if (restrict_support_intersection) "removed" else "kept (support not restricted)"), call. = FALSE)
      if (restrict_support_intersection) {
        encodedA <- encodedA[encodedA[, p] %in% c(NA, common), ]
        encodedNewB <- encodedNewB[encodedNewB[, p] %in% c(NA, common), ]
      }
    }
  }
  rownames(encodedA) <- seq_len(nrow(encodedA))
  rownames(encodedNewB) <- seq_len(nrow(encodedNewB))

  list(dataA = encodedA, dataB = encodedNewB)
}

#' Standardised mean difference between a selected set and a baseline
#' 
#' When missing values are encoded as 0m the smd will report information on the 
#' missingness in the linked sample vs. the source.
#'
#' @param data_select Data frame to compare (e.g. the linked
#'   set vs. the original file).
#' @param data_baseline Data frame to compare (e.g. the linked
#'   set vs. the original file).
#' @param var Character, column name to compare.
#' @param continuous Logical; if `TRUE`, compares means of `var` directly; if
#'   `FALSE`, compares the proportion at each observed level of `var`
#'   (default `TRUE`).
#'
#' @return Named list, one SMD per variable (`continuous = TRUE`) or per level
#'   (`continuous = FALSE`).
#' @export
#'
#' @examples
#' base <- data.frame(age = rnorm(200, 40, 10), sex = sample(c("M", "F"),
#'                     200, replace = TRUE))
#' select <- data.frame(age = rnorm(80, 43, 10), sex = sample(c("M", "F"), 80,
#'                     replace = TRUE, prob = c(0.6, 0.4)))
#' smd(select, base, "age")
#' smd(select, base, "sex", continuous = FALSE)
smd <- function(data_select, data_baseline, var, continuous = TRUE) {
  res <- list()
  if (continuous) {
    res[[var]] <- (mean(data_select[, var], na.rm = TRUE) - mean(data_baseline[, var], na.rm = TRUE)) / stats::sd(data_baseline[, var], na.rm = TRUE)
  } else {
    for (val in sort(unique(data_baseline[, var]))) {
      res[[paste(var, val, sep = "_")]] <-
        (mean(data_select[, var] == val, na.rm = TRUE) - mean(data_baseline[, var] == val, na.rm = TRUE)) /
        stats::sd(data_baseline[, var] == val, na.rm = TRUE)
    }
  }
  res
}

#' Plot the distribution of linkage scores
#'
#' @param n_pairs Integer, total number of candidate pairs considered
#'   (`nrow(A) * nrow(B)`).
#' @param LinkScore Numeric vector, linkage scores of the pairs above
#'   0 (e.g. `Delta$x`).
#'
#' @return `NULL`, invisibly; called for its plotting side effect.
#' @export
#'
#' @examples
#' fit_Delta_x <- c(0,0,0,0,0,0,0,0,0,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,
#' 0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.5,0.6,0.7,0.7,0.7,0.7,0.7,0.7,0.8,0.8)
#' plot_linkage_scores(1000, fit_Delta_x)
plot_linkage_scores <- function(n_pairs, LinkScore) {
  breaks <- seq(0, 1, by = 0.05)
  h <- graphics::hist(LinkScore, plot = FALSE, breaks = breaks)
  h$counts[1] <- h$counts[1] + (n_pairs - length(LinkScore))
  h$counts <- pmax(log10(h$counts), 0)
  graphics::plot(h, xlim = c(0, 1), main = "Linkage score distribution",
                 xlab = "Posterior linkage score", ylab = "Log10 frequency")
  invisible(NULL)
}

#' Model specific, linkage score based, false discovery proportion at 
#' a given threshold
#'
#' @param LinkScore Numeric vector, linkage scores of the candidate pairs.
#' @param threshold Numeric, score threshold above which a pair is declared
#'   linked.
#'
#' @return A list with `FDP_score` Numeric, `1 - mean(score | score > threshold)`, 
#'   the estimated FDP at `threshold` and `n_linked` Integer, number of linked
#'   records.
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V2 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V3 = list(dynamics = "flexible",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(NA,NA)),
#'                      V4 = list(dynamics = "structured",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(0.03,0.03),
#'                                cond_hazard_cov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())) )
#' n_values  <- c( 5, 6, 7, 12 )
#' p_mistake <- list( V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                    V3 = c(0.05, 0.05), V4 = c(0.02, 0.02) )
#' p_missing <- list( V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                    V3 = c(0.005, 0.005), V4 = c(0.005, 0.005) )
#' cond_hazard_params <- list( V1 = c(), V2 = c(), 
#'                              V3 = c(), V4 = log(c(0.7, 0.6, 0.5)) )
#' gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
#'                            p_mistake, p_missing, cond_hazard_params, TRUE )
#' prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
#'                            PIVs_config, TRUE, "entity_id", TRUE )
#' fit <- StEM( data = prep_data, StEM_iter = 10, StEM_burnin = 5,
#'              gibbs_iter = 10, gibbs_burnin = 5, n_post_sample = 10 )
#' FDP_score(fit$Delta$x, 0.5)
#' FDP_score(fit$Delta$x, 0.75)
FDP_score <- function(LinkScore, threshold) {
  linked <- LinkScore > threshold
  list(
     FDP_score = 1 - sum(LinkScore[linked]) / sum(linked),
     n_linked = sum(linked)
   )
}

#' Model agnostic, based on synthetic data, false discovery proportion at 
#' a given threshold
#'
#' @param idxA Integer vector, indices in A of linked records (for a 
#'   previously set threshold). 
#' @param idxB Integer vector, indices in B of linked records (for a 
#'   previously set threshold). 
#' @param LinkScore Numeric vector, linkage scores of the candidate 
#'   pairs. If `NULL`, consider all given pairs as linked.
#' @param threshold Numeric, score threshold above which a pair is declared
#'   linked.
#' @param n_records_A Integer, number of records in A.
#' @param n_records_B Integer, number of records in B.
#' @param n_synth Integer, number of synthetic records to generate per
#'   iteration.
#'
#' @return A list with `FDP_synth` Numeric, proportion of synthetic falsely
#'   linked records, `n_linked_real` Integer, number of real linked records and
#'   `n_linked_all` Integer, total number of linked records (real and synthetic).
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V2 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V3 = list(dynamics = "flexible",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(NA,NA)),
#'                      V4 = list(dynamics = "structured",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(0.03,0.03),
#'                                cond_hazard_cov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())) )
#' n_values  <- c( 5, 6, 7, 12 )
#' p_mistake <- list( V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                    V3 = c(0.05, 0.05), V4 = c(0.02, 0.02) )
#' p_missing <- list( V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                    V3 = c(0.005, 0.005), V4 = c(0.005, 0.005) )
#' cond_hazard_params <- list( V1 = c(), V2 = c(), 
#'                              V3 = c(), V4 = log(c(0.7, 0.6, 0.5)) )
#' gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
#'                            p_mistake, p_missing, cond_hazard_params, TRUE )
#' prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
#'                            PIVs_config, TRUE, "entity_id", TRUE )
#' PIVs <- names(PIVs_config)
#' n_synth <-as.integer(0.10 * nrow(prep_data$encodedB))
#' new_data <- synthesise("arf", prep_data$encodedA[, c(PIVs,"local_id","source",
#'                       "date",prep_data$PIVs_config$V4$cond_hazard_cov$covA)],
#'                       prep_data$encodedB[, c(PIVs,"local_id","source","date",
#'                       prep_data$PIVs_config$V4$cond_hazard_cov$covB)], PIVs, 
#'                       n_synth, TRUE)
#' new_data$dataA[PIVs][is.na(new_data$dataA[PIVs])] <- 0
#' new_data$dataB[PIVs][is.na(new_data$dataB[PIVs])] <- 0
#' # cannot model dynamics for synthetic data, set dates to 0
#' new_data$dataA$date[is.na(new_data$dataA$date)] <- 0
#' new_data$dataB$date[is.na(new_data$dataB$date)] <- 0
#' arguments <- list(data = prep_data, StEM_iter = 5, StEM_burnin = 2, 
#'                   gibbs_iter = 5, gibbs_burnin = 2, n_post_sample = 5)
#' fit_flexrl <- link_with_FlexRL(new_data$dataA, new_data$dataB, arguments)  
#' FDP_synth(fit_flexrl$idxA, fit_flexrl$idxB, fit_flexrl$LinkScore, 0.5,
#'           nrow(prep_data$encodedA), nrow(prep_data$encodedB), 
#'           n_synth)
FDP_synth <- function(idxA, idxB, LinkScore, threshold,
                      n_records_A, n_records_B, n_synth) {
  if (!is.null(LinkScore)){
    keep <- LinkScore > threshold
    idxA <- idxA[keep]
    idxB <- idxB[keep]
  }
  
  real <- idxB <= n_records_B & idxA <= n_records_A
  synthfp <- sum(!real)
  n_linked <- length(idxA)
  list(
    FDP_synth = (synthfp * (n_records_B / n_synth)) / (n_linked - synthfp),
    n_linked_real = sum(real),
    n_linked_all = n_linked
  )
}

#' Estimate the false discovery proportion of a record-linkage method via
#' model specific linkage scores
#'
#' @param encodedA Encoded data source as returned by [prepare_data()] 
#'   (`encodedA` is the smaller one).
#' @param encodedB Encoded data source as returned by [prepare_data()] 
#'   (`encodedB` is the larger one).
#' @param PIVs Character vector, names of the PIVs.
#' @param maxIter4CV Integer, max number of retries per iteration if no valid
#'   FDP estimate is obtained.
#' @param n_repeats Integer, number of augmentation iterations to average over.
#' @param RL_method One of `"multilink"`, `"fastLink"`, `"BRL"`,
#'   `"reclin2"`, `"diyar"`, `"fedmatch"`, `"FlexRL"`.
#' @param ... Extra arguments forwarded to the chosen link_with_* wrapper
#'   (i.e. to the underlying record-linkage package).
#'
#' @return List with `FDP_score_estimator`, `Linked_pairs`:
#'   data frames (`n_repeats` rows x 50 thresholds, 0.50 to 0.99).
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V2 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V3 = list(dynamics = "flexible",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(NA,NA)),
#'                      V4 = list(dynamics = "structured",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(0.03,0.03),
#'                                cond_hazard_cov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())) )
#' n_values  <- c( 5, 6, 7, 12 )
#' p_mistake <- list( V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                    V3 = c(0.05, 0.05), V4 = c(0.02, 0.02) )
#' p_missing <- list( V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                    V3 = c(0.005, 0.005), V4 = c(0.005, 0.005) )
#' cond_hazard_params <- list( V1 = c(), V2 = c(), 
#'                              V3 = c(), V4 = log(c(0.7, 0.6, 0.5)) )
#' gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
#'                            p_mistake, p_missing, cond_hazard_params, TRUE )
#' prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
#'                            PIVs_config, TRUE, "entity_id", TRUE )
#' PIVs <- names(PIVs_config)                          
#' compute_RL_FDP_score( prep_data$encodedA, prep_data$encodedB, PIVs, 1, 2,
#'                       "BRL", flds = PIVs, types = rep("bi",length(PIVs)) )
#' compute_RL_FDP_score( prep_data$encodedA, prep_data$encodedB, PIVs, 1, 2,
#'                       "FlexRL", data = prep_data, StEM_iter = 5, 
#'                       StEM_burnin = 2, gibbs_iter = 5, gibbs_burnin = 2,
#'                       n_post_sample = 5 )
compute_RL_FDP_score <- function(encodedA, encodedB, PIVs, maxIter4CV = 10, n_repeats = 10,
                                 RL_method, ...) {
  n_records_A <- nrow(encodedA)
  n_records_B <- nrow(encodedB)
  if (n_records_A > n_records_B) stop("`encodedA` must be smaller than `encodedB`.", call. = FALSE)
  
  supported <- c("multilink", "fastLink", "BRL", "reclin2", "diyar", "fedmatch", "FlexRL")
  if (!RL_method %in% supported) stop("`RL_method` must be one of: ", paste(supported, collapse = ", "), ".", call. = FALSE)
  
  thresholds <- seq(0.5, 0.99, by = 0.01)
  th_names <- sprintf("%.2f", thresholds)
  FDP_score_res <- stats::setNames(data.frame(matrix(NA, n_repeats, 50)), th_names)
  n_linked_score_res <- stats::setNames(data.frame(matrix(NA, n_repeats, 50)), th_names)
  arguments <- list(...)
  
  run_RL <- function(encodedA, encodedB) {
    switch(RL_method,
           FlexRL     = link_with_FlexRL(encodedA, encodedB, arguments),
           fedmatch   = link_with_fedmatch(encodedA, encodedB, arguments),
           reclin2    = link_with_reclin2(encodedA, encodedB, arguments),
           BRL        = link_with_BRL(encodedA, encodedB, arguments),
           fastLink   = link_with_fastLink(encodedA, encodedB, arguments),
           multilink  = link_with_multilink(encodedA, encodedB, arguments),
           diyar      = link_with_diyar(encodedA, encodedB, arguments)
    )
  }
  
  for (i in seq_len(n_repeats)) {
    
    any_valid_estimate <- FALSE
    n_tries <- 0
    
    while (n_tries < maxIter4CV && !any_valid_estimate) {
      res <- run_RL(encodedA, encodedB)
      res_idxA <- res$idxA
      res_idxB <- res$idxB
      res_LinkScore <- res$LinkScore
      
      if (!is.null(res_LinkScore)) {
        for (j in seq_along(thresholds)) {
          keep <- res_LinkScore > thresholds[j]
          n_linked <- sum(keep, na.rm = TRUE)
          if (n_linked > 0) {
            est_score <- FDP_score(res_LinkScore, thresholds[j])
            FDP_score_res[i, j] <- est_score$FDP_score
            n_linked_score_res[i, j] <- est_score$n_linked
          } else {
            if (j == 1) warning(sprintf("Nothing linked at iteration %s at threshold 0.50.", i), call. = FALSE)
            FDP_score_res[i, j:50] <- 0
            n_linked_score_res[i, j:50] <- 0
            break
          }
        }
      }
      
      n_tries <- n_tries + 1
      any_valid_estimate <- any(!is.na(FDP_score_res[i, ]) & FDP_score_res[i, ] <= 1 & n_linked_score_res[i, ] > 0)
    }
    
    if (n_tries == maxIter4CV && !any_valid_estimate) {
      FDP_score_res[FDP_score_res>1] <- NA
      warning(sprintf(
        "No valid FDP estimate after %s attempts at iteration %s. Increase `maxIter4CV`, or the estimator may be unreliable for this method/data.",
        maxIter4CV, i
      ), call. = FALSE)
      break
    }
  }
  
  to_show <- data.frame(
    `FDP model score estimator `     = round(colMeans(FDP_score_res, na.rm = TRUE), 2),
    `Linked pairs (RL)         `     = round(colMeans(n_linked_score_res, na.rm = TRUE)),
    check.names = FALSE
  )
  out <- t(to_show[c(1,6,11,16,21,26,31,36,41,46),])
  out <- rbind("Decision rule threshold   " = as.numeric(colnames(out)), out)
  colnames(out) <- rep("", ncol(out))
  
  message(sprintf("%s results (average over %s iterations) (original RL task) `compute_RL_FDP_score`", RL_method, n_repeats))
  utils::write.table(round(out,2), quote = FALSE, col.names = FALSE, sep = "\t")
  
  list(
    FDP_score_estimator = FDP_score_res,
    Linked_pairs        = n_linked_score_res
  )
}

#' Estimate the false discovery proportion of a record-linkage method via
#' synthetic augmentation
#'
#' Repeatedly augments file B with synthetic records ([synthesise()]), runs
#' the chosen record-linkage method link_with_*, and compares the 
#' synthetic-vs-real proportion among linked pairs to estimate the false 
#' discovery proportion, for a range of score thresholds (0.50 to 0.99). 
#' See https://doi.org/10.1002/sim.70292 for the method.
#'
#' @param synth_method Passed to [synthesise()]: `"arf"`, `"synthpop"`, or
#'   `"mice"`.
#' @param encodedA,encodedB The two encoded data sources as returned by
#'   [prepare_data()] (`encodedA` is the smaller one).
#' @param PIVs Character vector, names of the PIVs.
#' @param n_synth Integer, number of synthetic records to generate per
#'   iteration (default: 10% of `nrow(encodedB)`).
#' @param restrict_support_intersection Passed to [synthesise()].
#' @param maxIter4CV Integer, max number of retries per iteration if no valid
#'   FDP estimate is obtained.
#' @param n_repeats Integer, number of augmentation iterations to average over.
#' @param RL_method One of `"multilink"`, `"fastLink"`, `"BRL"`,
#'   `"reclin2"`, `"diyar"`, `"fedmatch"`, `"FlexRL"`.
#' @param ... Extra arguments forwarded to the chosen link_with_* wrapper
#'   (i.e. to the underlying record-linkage package).
#'
#' @return List with `FDP_score_estimator`, `FDP_synth_estimator`,
#'   `Linked_pairs_augm`, `Linked_pairs`: data frames (`n_repeats` rows x 50
#'   thresholds, 0.50 to 0.99).
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V2 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V3 = list(dynamics = "flexible",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(NA,NA)),
#'                      V4 = list(dynamics = "structured",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(0.03,0.03),
#'                                cond_hazard_cov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())) )
#' n_values  <- c( 5, 6, 7, 12 )
#' p_mistake <- list( V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                    V3 = c(0.05, 0.05), V4 = c(0.02, 0.02) )
#' p_missing <- list( V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                    V3 = c(0.005, 0.005), V4 = c(0.005, 0.005) )
#' cond_hazard_params <- list( V1 = c(), V2 = c(), 
#'                              V3 = c(), V4 = log(c(0.7, 0.6, 0.5)) )
#' gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
#'                            p_mistake, p_missing, cond_hazard_params, TRUE )
#' prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
#'                            PIVs_config, TRUE, "entity_id", TRUE )
#' PIVs <- names(PIVs_config)                         
#' compute_augmRL_FDP_synth( "arf", prep_data$encodedA, prep_data$encodedB, PIVs, 
#'                           NULL, TRUE, 1, 2, "BRL", flds = PIVs, 
#'                           types = rep("bi",length(PIVs)) )
#' compute_augmRL_FDP_synth( "arf", prep_data$encodedA, prep_data$encodedB, PIVs, 
#'                           NULL, TRUE, 1, 2, "FlexRL", data = prep_data, 
#'                           StEM_iter = 5, StEM_burnin = 2, 
#'                           gibbs_iter = 5, gibbs_burnin = 2, n_post_sample = 5 )
compute_augmRL_FDP_synth <- function(synth_method, encodedA, encodedB, PIVs, n_synth = NULL,
                                    restrict_support_intersection = TRUE, maxIter4CV = 10, n_repeats = 10,
                                    RL_method, ...) {

  n_records_A <- nrow(encodedA)
  n_records_B <- nrow(encodedB)
  if (n_records_A > n_records_B) stop("`encodedA` must be smaller than `encodedB`.", call. = FALSE)

  supported <- c("multilink", "fastLink", "BRL", "reclin2", "diyar", "fedmatch", "FlexRL")
  if (!RL_method %in% supported) stop("`RL_method` must be one of: ", paste(supported, collapse = ", "), ".", call. = FALSE)

  if (is.null(n_synth)) n_synth <- as.integer(0.10 * n_records_B)

  thresholds <- seq(0.5, 0.99, by = 0.01)
  th_names <- sprintf("%.2f", thresholds)
  FDP_synth_res <- stats::setNames(data.frame(matrix(NA, n_repeats, 50)), th_names)
  FDP_score_res <- stats::setNames(data.frame(matrix(NA, n_repeats, 50)), th_names)
  n_linked_synth_res <- stats::setNames(data.frame(matrix(NA, n_repeats, 50)), th_names)
  n_linked_score_res <- stats::setNames(data.frame(matrix(NA, n_repeats, 50)), th_names)
  arguments <- list(...)

  run_RL <- function(dataA, dataB) {
    switch(RL_method,
           FlexRL     = link_with_FlexRL(dataA, dataB, arguments),
           fedmatch   = link_with_fedmatch(dataA, dataB, arguments),
           reclin2    = link_with_reclin2(dataA, dataB, arguments),
           BRL        = link_with_BRL(dataA, dataB, arguments),
           fastLink   = link_with_fastLink(dataA, dataB, arguments),
           multilink  = link_with_multilink(dataA, dataB, arguments),
           diyar      = link_with_diyar(dataA, dataB, arguments)
    )
  }

  for (i in seq_len(n_repeats)) {

    if (RL_method == "FlexRL") {
      PIVs_stable <- sapply(arguments$data$PIVs_config, function(x) x$dynamics != "structured")
      if (any(!PIVs_stable)) {
        covariates_to_add_A <- c()
        covariates_to_add_B <- c()
        for (k in seq_len(length(arguments$data$PIVs_config))){
          if (!PIVs_stable[k]){
            covariates_to_add_A <- c(covariates_to_add_A, arguments$data$PIVs_config[[k]]$cond_hazard_cov$covA)
            covariates_to_add_B <- c(covariates_to_add_B, arguments$data$PIVs_config[[k]]$cond_hazard_cov$covB)
          }
        }
        new_data <- synthesise(synth_method, 
                              encodedA[, c(PIVs, "local_id", "source", "date", covariates_to_add_A)],
                              encodedB[, c(PIVs, "local_id", "source", "date", covariates_to_add_B)], 
                              PIVs, n_synth, 
                              restrict_support_intersection)
        # cannot model dynamics for synthetic data, set dates to 0
        new_data$dataA$date[is.na(new_data$dataA$date)] <- 0
        new_data$dataB$date[is.na(new_data$dataB$date)] <- 0
      } else {
        new_data <- synthesise(synth_method, encodedA[, c(PIVs, "local_id", "source")],
                              encodedB[, c(PIVs, "local_id", "source")], PIVs, n_synth, 
                              restrict_support_intersection)
      }
    } else {
      new_data <- synthesise(synth_method, encodedA[, c(PIVs, "local_id", "source")],
                            encodedB[, c(PIVs, "local_id", "source")], PIVs, n_synth, 
                            restrict_support_intersection)
    }

    any_valid_estimate <- FALSE
    n_tries <- 0

    while (n_tries < maxIter4CV && !any_valid_estimate) {
      res <- run_RL(new_data$dataA, new_data$dataB)
      res_idxA <- res$idxA
      res_idxB <- res$idxB
      res_LinkScore <- res$LinkScore

      if (is.null(res_LinkScore)) {
        n_linked <- length(res_idxA)
        if (n_linked > 0) {
          est <- FDP_synth(res_idxA, res_idxB, res_LinkScore, NULL, n_records_A, n_records_B, n_synth)
          FDP_synth_res[i, 1] <- est$FDP_synth
          n_linked_synth_res[i, 1] <- est$n_linked_real
        } else {
          warning(sprintf("Nothing linked at iteration %s with default parameters.", i), call. = FALSE)
          FDP_synth_res[i, 1] <- 0
          n_linked_synth_res[i, 1] <- 0
        }
      } else {
        for (j in seq_along(thresholds)) {
          n_linked <- sum(res_LinkScore > thresholds[j], na.rm = TRUE)
          if (n_linked > 0) {
            est_synth <- FDP_synth(res_idxA, res_idxB, res_LinkScore, thresholds[j], n_records_A, n_records_B, n_synth)
            est_score <- FDP_score(res_LinkScore, thresholds[j])
            FDP_synth_res[i, j] <- est_synth$FDP_synth
            FDP_score_res[i, j] <- est_score$FDP_score
            n_linked_synth_res[i, j] <- est_synth$n_linked_real
            n_linked_score_res[i, j] <- est_score$n_linked
          } else {
            if (j == 1) warning(sprintf("Nothing linked at iteration %s at threshold 0.50.", i), call. = FALSE)
            FDP_synth_res[i, j:50] <- 0
            FDP_score_res[i, j:50] <- 0
            n_linked_synth_res[i, j:50] <- 0
            n_linked_score_res[i, j:50] <- 0
            break
          }
        }
      }

      n_tries <- n_tries + 1
      any_valid_estimate <- any(!is.na(FDP_synth_res[i, ]) & FDP_synth_res[i, ] <= 1 & n_linked_synth_res[i, ] > 0)
    }

    if (n_tries == maxIter4CV && !any_valid_estimate) {
      FDP_synth_res[FDP_synth_res>1] <- NA
      FDP_score_res[FDP_score_res>1] <- NA
      warning(sprintf(
        "No valid FDP estimate after %s attempts at iteration %s. Increase `maxIter4CV`, or the estimator may be unreliable for this method/data.",
        maxIter4CV, i
      ), call. = FALSE)
      break
    }
  }

  to_show <- data.frame(
    `FDP model score estimator `     = round(colMeans(FDP_score_res, na.rm = TRUE), 2),
    `FDP synth data estimator  `     = round(colMeans(FDP_synth_res, na.rm = TRUE), 2),
    `Linked pairs (augm. RL)   `     = round(colMeans(n_linked_synth_res, na.rm = TRUE)),
    `Linked pairs (RL)         `     = round(colMeans(n_linked_score_res, na.rm = TRUE)),
    check.names = FALSE
  )
  out <- t(to_show[c(1,6,11,16,21,26,31,36,41,46),])
  out <- rbind("Decision rule threshold   " = as.numeric(colnames(out)), out)
  colnames(out) <- rep("", ncol(out))

  message(sprintf("%s results (average over %s iterations) (augmented RL task) `compute_augmRL_FDP_synth`", RL_method, n_repeats))
  utils::write.table(round(out,2), quote = FALSE, col.names = FALSE, sep = "\t")

  list(
    FDP_score_estimator = FDP_score_res,
    FDP_synth_estimator = FDP_synth_res,
    Linked_pairs_augm   = n_linked_synth_res,
    Linked_pairs        = n_linked_score_res
  )
}

# ============================================================================
# S3 class "RL_diagnostics" for post-linkage diagnostics, built on 
# RL_agreement(), smd(), support_iou(), mmd(), FDP_score(), FDP_synth(), 
# plot_linkage_scores(), plot_distributions(), plot_StEM_convergence(),
# plot.FDP_curves(), plot.discrepancy_curves().
# 
# diag <- RL_diagnostics(fit_flexrl, prep_data$encodedA, prep_data$encodedB,
#                        PIVs, PIVs_type, PIVs, prep_data$true_pairs, TRUE, 
#                        RL_method = "FlexRL", data = prep_data,
#                        StEM_iter = 10, StEM_burnin = 5, 
#                        gibbs_iter = 10, gibbs_burnin = 5,
#                        maxIter4CV = 3, n_repeats = 5)
# 
# diag                              # print(), diagnostics summary
# plot(diag,"scores")               # linkage score histogram
# plot(diag,"distributions")        # linked subset vs. data source
# plot(diag,"convergence")          # StEM trace plots
# plot(diag,"FDP")                  # FDP estimates
# plot(diag,"discrepancy")          # measures of linked data discrepancy
# ============================================================================

#' Monte Carlo convergence (trace) plots for a fitted StEM model
#'
#' Trace-plots the raw StEM chains of `gamma`, `eta`, `alpha` and `phi`
#' across iterations, to visually judge whether the chains have stabilised
#' and pick an adequate `StEM_burnin` for [StEM()]. One plot per parameter:
#' `gamma` (proportion linked), `eta` (PIVs distribution), `alpha` (hazard 
#' coefficients, for dynamic PIVs only), and `phi` (agreement/missing rates).
#'
#' @param fit List as returned by [StEM()], containing the raw chains
#'   `gamma`, `eta`, `alpha`, `phi` (StEM_iter rows each).
#'
#' @return `NULL`, invisibly; called for its plotting side effect.
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V2 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V3 = list(dynamics = "flexible",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(NA,NA)),
#'                      V4 = list(dynamics = "structured",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(0.03,0.03),
#'                                cond_hazard_cov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())) )
#' n_values  <- c( 5, 6, 7, 12 )
#' p_mistake <- list( V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                    V3 = c(0.05, 0.05), V4 = c(0.02, 0.02) )
#' p_missing <- list( V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                    V3 = c(0.005, 0.005), V4 = c(0.005, 0.005) )
#' cond_hazard_params <- list( V1 = c(), V2 = c(), 
#'                              V3 = c(), V4 = log(c(0.7, 0.6, 0.5)) )
#' gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
#'                            p_mistake, p_missing, cond_hazard_params, TRUE )
#' prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
#'                            PIVs_config, TRUE, "entity_id", TRUE )
#' fit <- StEM( data = prep_data, StEM_iter = 10, StEM_burnin = 3,
#'              gibbs_iter = 10, gibbs_burnin = 3, n_post_sample = 10 )
#' plot_StEM_convergence(fit)
plot_StEM_convergence <- function(fit) {
  required <- c("gamma", "eta", "alpha", "phi")
  if (!all(required %in% names(fit))) {
    stop("`fit` must contain the StEM chains: ", paste(required, collapse = ", "), ".", call. = FALSE)
  }
  PIVs <- names(fit$eta)

  op <- graphics::par(mfrow = c(1, 1))
  on.exit(graphics::par(op))
  # gamma: proportion of linked records (a fraction of the smaller file)
  graphics::matplot(fit$gamma, type = "l", lty = 1, col = "steelblue",
                    main = "Convergence: gamma (proportion linked)",
                    xlab = "StEM iteration", ylab = expression(gamma))

  # eta: distribution of true values, one panel per PIV
  for (k in seq_along(fit$eta)) {
    graphics::matplot(fit$eta[[k]], type = "l", lty = 1,
                      main = sprintf("Convergence: eta['%s']", PIVs[k]),
                      xlab = "StEM iteration", ylab = expression(eta),
                      ylim=c(0,0.5))
  }

  # alpha: hazard coefficients (structured PIVs only)
  for (k in seq_along(fit$alpha)) {
    ak <- fit$alpha[[k]]
    if (is.null(ncol(ak)) || ncol(ak) == 0 || all(!is.finite(ak))) next
    graphics::matplot(ak, type = "l", lty = 1, col = seq_len(ncol(ak)),
                      main = sprintf("Convergence: alpha['%s']", PIVs[k]),
                      xlab = "StEM iteration", ylab = expression(alpha))
    if (!is.null(colnames(ak))) {
      graphics::legend("bottomright", colnames(ak), col = seq_len(ncol(ak)), lty = 1, cex = 0.7, bty = "n")
    }
  }

  # phi: registration-error parameters (agreement A, agreement B, missing A, missing B)
  for (k in seq_along(fit$phi)) {
    graphics::matplot(fit$phi[[k]], type = "l", lty = 1,
                      col = c("steelblue", "firebrick", "darkgreen", "orange"),
                      main = sprintf("Convergence: phi['%s']", PIVs[k]),
                      xlab = "StEM iteration", ylab = expression(phi),
                      ylim = c(0,1),
                      lwd = rep_len(c(2, 1), ncol(fit$phi[[k]]))    )
    graphics::legend("right", c("agree A", "agree B", "missing A", "missing B"),
                     col = c("steelblue", "firebrick", "darkgreen", "orange"), 
                     lty = 1, cex = 0.7, bty = "n",
                     lwd = rep_len(c(2, 1), ncol(fit$phi[[k]])))
  }

  invisible(NULL)
}

# Divergence measures between the linked subset (for the given threshold) and 
# the source files (agreement among PIVs of linked profiles, smd, iou, mmd). 
.discrepancy_measures <- function(fit, encodedA, encodedB, compare_vars, vars_type_cont, threshold) {
  
  if (is.null(fit$LinkScore)) {
    n_linked <- length(fit$idxA)
    linked <- rep(TRUE, n_linked)
    linkedA <- encodedA[fit$idxA, , drop = FALSE]
    linkedB <- encodedB[fit$idxB, , drop = FALSE]
  } else {
    linked <- fit$LinkScore > threshold
    n_linked <- sum(linked)
    linkedA <- encodedA[fit$idxA[linked], , drop = FALSE]
    linkedB <- encodedB[fit$idxB[linked], , drop = FALSE]
  }
  
  varsA  <- intersect(compare_vars, names(encodedA))
  varsB  <- intersect(compare_vars, names(encodedB))
  varsAB <- intersect(varsA, varsB)
  
  agreement <- if (n_linked > 0 && length(varsAB) > 0) {
    RL_agreement(encodedA, encodedB, varsAB, data.frame(fit$idxA[linked], fit$idxB[linked]))$agreements
  } else {
    stats::setNames(rep(NA_real_, length(varsAB)), varsAB)
  }
  
  smd_by_var <- function(linkedX, encodedX, vars) {
    stats::setNames(lapply(vars, function(v) {
      if (n_linked == 0) return(NA_real_)
      res <- smd(linkedX, encodedX, v, continuous = vars_type_cont[[v]])
      if (vars_type_cont[[v]]) res[[1]] else res
    }), vars)
  }
  iou_by_var <- function(linkedX, encodedX, vars) {
    stats::setNames(vapply(vars, function(v) {
      if (n_linked < 2) return(NA_real_)
      support_iou(graphics::hist(as.numeric(linkedX[[v]]), plot = FALSE),
                  graphics::hist(as.numeric(encodedX[[v]]), plot = FALSE))
    }, numeric(1)), vars)
  }
  mmd_all <- function(linkedX, encodedX, vars) {
    if (n_linked < 2 || length(vars) == 0) return(NA_real_)
    mmd(stats::na.omit(as.matrix(linkedX[, vars, drop = FALSE])), stats::na.omit(as.matrix(encodedX[, vars, drop = FALSE])))
  }
  
  list(n_linked = n_linked, linkedA = linkedA, linkedB = linkedB, agreement = agreement,
       smdA = smd_by_var(linkedA, encodedA, varsA), smdB = smd_by_var(linkedB, encodedB, varsB),
       iouA = iou_by_var(linkedA, encodedA, varsA), iouB = iou_by_var(linkedB, encodedB, varsB),
       mmdA = mmd_all(linkedA, encodedA, varsA), mmdB = mmd_all(linkedB, encodedB, varsB))
}

#' Plot post-linkage diagnostics over thresholds
#'
#' @param x A `discrepancy_curves` object.
#' @param ... Ignored.
#'
#' @return `x`, invisibly; called for its plotting side effect.
#' @export
plot.discrepancy_curves <- function(x, ...) {
  xi <- x$xi
  cols <- function(prefix) {
    v <- as.data.frame(x[grep(paste0("^", prefix, "\\."), names(x))])
    names(v) <- sub(paste0("^", prefix, "\\.([^.]*\\.)?"), "", names(v))
    v
  }
  lines_by_col <- function(vals, ylab, main, ylim = NULL, where = "bottomleft") {
    if (ncol(vals) == 0 || all(is.na(vals))) return(invisible(NULL))
    if (is.null(ylim)) ylim <- range(vals, na.rm = TRUE)
    if (any(duplicated(as.list(vals)) )){
      perturbation <- duplicated(as.list(vals))*c(stats::runif(length(duplicated(as.list(vals))),0,0.005))
      perturbation <- matrix(perturbation, nrow = nrow(vals), ncol = length(perturbation), byrow = TRUE)
    }else{
      perturbation <- 0
    }
    graphics::matplot(xi, vals+perturbation, type = "l", lty = 1, col = seq_len(ncol(vals)),
                      xlab = "decision rule threshold", ylab = ylab, main = main, ylim = ylim)
    graphics::legend(where, colnames(vals), col = seq_len(ncol(vals)), lty = 1,
                     bty = "n", cex = 0.7, ncol = max(1, ncol(vals) %/% 8))
  }
  two_files <- function(a, b, ylab, main, ylim = NULL, where = "topleft") {
    if (all(is.na(a)) && all(is.na(b))) return(invisible(NULL))
    if (is.null(ylim)) ylim <- range(c(a, b), na.rm = TRUE)
    graphics::plot(xi, a, type = "l", xlab = "decision rule threshold", ylab = ylab, main = main, ylim = ylim)
    graphics::lines(xi, b, lty = 2)
    graphics::legend(where, c("linked A vs. A", "linked B vs. B"), lty = c(1, 2), bty = "n", cex = 0.8)
  }
  op <- graphics::par(mfrow = c(1, 1))
  on.exit(graphics::par(op))
  two_files(x$mmdA, x$mmdB, "MMD", "Multivariate discrepancy")
  agreements <- cbind(cols("agreementAB"), cols("agreementLinks"))
  names(agreements) <- c(names(cols("agreementAB")), sprintf("%s true pairs", names(cols("agreementLinks"))))
  lines_by_col(agreements, "agreement rate", "linked A vs. linked B", c(0, 1.03))
  op <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(op))
  lines_by_col(cols("iouA"), "support IoU", "linked A vs. A", c(0, 1.03))
  lines_by_col(cols("iouB"), "support IoU", "linked B vs. B", c(0, 1.03))
  op <- graphics::par(mfrow = c(1, 2))
  on.exit(graphics::par(op))
  allsmd <- c(as.matrix(cols("smdA")), as.matrix(cols("smdB")))
  ylim <- if (any(is.finite(allsmd))) range(allsmd[is.finite(allsmd)]) else NULL
  lines_by_col(cols("smdA"), "SMD", "linked A vs. A", ylim, "topleft")
  lines_by_col(cols("smdB"), "SMD", "linked B vs. B", ylim, "topleft")
  invisible(x)
}

#' Post-linkage diagnostics
#'
#' Gathers, from a fitted linkage, the diagnostics needed to judge the linked 
#' data before using it for inference: false discovery proportion estimate, 
#' agreement among PIVs of linked pairs, comparison of the linked subset with 
#' each source file (standardised mean differences, support overlap, maximum 
#' mean discrepancy). 
#'
#' @param fit List with either `Delta` (a data frame with columns `i`, `j`,
#'   `x`, as returned by [StEM()]) or the three elements `idxA`, `idxB`,
#'   `LinkScore` (e.g. the output of a link_with_* wrapper).
#' @param encodedA Data source used for linkage (same encoding and row order as 
#'   passed to the linkage method).
#' @param encodedB Data source used for linkage (same encoding and row order as 
#'   passed to the linkage method).
#' @param compare_vars Character vector, variables (PIVs or others) to
#'   compare between the linked subset and each source file.
#' @param vars_type_cont Named list or logical vector, one entry per
#'   `compare_vars`: `TRUE` if the variable is treated as continuous,
#'   `FALSE` if categorical (one SMD per level).
#' @param PIVs Character vector, partially identifying variables to
#'   compare between the linked records from A and from B.
#' @param true_pairs Optional data frame with 2 columns of true (A, B)
#'   indices, when known, to report the realised FDP and sensitivity.
#' @param FDP_estimation Logical; if `TRUE`, run [compute_RL_FDP_score()] and
#'   [compute_augmRL_FDP_synth()] to estimate the FDP over thresholds 0.50 to
#'   0.95. `RL_method` and other arguments must then be given in `...`.
#' @param ... Arguments passed to [compute_RL_FDP_score()] and
#'   [compute_augmRL_FDP_synth()] when `FDP_estimation = TRUE`: `RL_method`,
#'   `n_repeats`, `maxIter4CV`, optionally `synth_method` (default `"arf"`),
#'   `n_synth`, `PIVs` (default `compare_vars`), and the arguments of
#'   the chosen link_with_* wrapper.
#'
#' @details All measures are computed for every threshold from 0.50 to 0.99
#'   (by 0.01) on the linkage scores; for a method without scores (or with a
#'   single score value) they are computed once, for all returned pairs.
#'   `print(x, threshold = )` summarises diagnostics at one threshold, `plot()`
#'   shows the curves.
#'
#' @return An object of class `"RL_diagnostics"`, a list with:
#'   \item{discrepancy_measures}{data frame (class `"discrepancy_curves"`),
#'      one row per threshold `xi`: number of linked pairs, agreement rate per
#'      variable (see [RL_agreement()]), standardised mean differences (one
#'      per continuous variable or per level, see [smd()]), support overlap
#'      per variable (see [support_iou()]) and multivariate discrepancy (see
#'      [mmd()]), each for the linked subset of A vs. A and of B vs. B}
#'   \item{FDP_measures}{FDP estimates by threshold (class `"FDP_curves"`) if
#'      `FDP_estimation = TRUE`, else `NULL`}
#'   \item{true_agreement}{only if `true_pairs` is supplied: agreement rate of
#'      the true pairs on each of `compare_vars`, to compare with the agreement
#'      of the linked pairs}     
#'   \item{true_performance}{only if `true_pairs` is supplied: data frame with
#'      the realised FDP and sensitivity by threshold}
#'   \item{idxA, idxB, LinkScore, RL_method, n_pairs, compare_vars, A, B}{the
#'      linkage and data needed by the `print()` and `plot()` methods}
#'   \item{gamma, eta, alpha, phi}{the StEM chains from `fit`, if present}
#' @export
#'
#' @examples
#' PIVs_config <- list( V1 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V2 = list(dynamics = "stable",
#'                                bound_mistakes = c(0.10,0.10),
#'                                fix_mistakes = c(NA,NA)),
#'                      V3 = list(dynamics = "flexible",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(NA,NA)),
#'                      V4 = list(dynamics = "structured",
#'                                bound_mistakes = c(NA,NA),
#'                                fix_mistakes = c(0.03,0.03),
#'                                cond_hazard_cov = list(cov1=c("Xe", "Xf"),
#'                                                       cov2=c())) )
#' n_values  <- c( 5, 6, 7, 12 )
#' p_mistake <- list( V1 = c(0.02, 0.02), V2 = c(0.02, 0.02),
#'                    V3 = c(0.05, 0.05), V4 = c(0.02, 0.02) )
#' p_missing <- list( V1 = c(0.005, 0.005), V2 = c(0.005, 0.005),
#'                    V3 = c(0.005, 0.005), V4 = c(0.005, 0.005) )
#' cond_hazard_params <- list( V1 = c(), V2 = c(), 
#'                              V3 = c(), V4 = log(c(0.7, 0.6, 0.5)) )
#' gen_data <- simulate_data( PIVs_config, n_values, c(150, 200), 100, 
#'                            p_mistake, p_missing, cond_hazard_params, TRUE )
#' prep_data <- prepare_data( gen_data$data1, gen_data$data2, "1", "2",
#'                            PIVs_config, TRUE, "entity_id", TRUE )
#' PIVs <- names(PIVs_config)  
#' PIVs_type <- list(V1=FALSE, V2=FALSE, V3=FALSE, V4=TRUE)
#' fit_flexrl <- StEM( data = prep_data, StEM_iter = 5, StEM_burnin = 2,
#'                     gibbs_iter = 5, gibbs_burnin = 2, n_post_sample = 5 )
#' diag_flexrl <- RL_diagnostics(fit_flexrl, prep_data$encodedA, prep_data$encodedB,
#'                               PIVs, PIVs_type, PIVs, true_pairs = prep_data$true_pairs,
#'                               FDP_estimation = TRUE, RL_method = "FlexRL", 
#'                               data = prep_data,
#'                               StEM_iter = 5, StEM_burnin = 2, 
#'                               gibbs_iter = 5, gibbs_burnin = 2,
#'                               n_post_sample = 5,
#'                               maxIter4CV = 1, n_repeats = 1)
#' diag_flexrl # print(diag_flexrl)
#' print(diag_flexrl, threshold = 0.75)
#' plot(diag_flexrl, "scores")
#' plot(diag_flexrl, "distributions", threshold = 0.75)
#' plot(diag_flexrl, "convergence")
#' plot(diag_flexrl, "FDP")
#' plot(diag_flexrl, "discrepancy")
#' 
#' fit_brl <- link_with_BRL( prep_data$encodedA, prep_data$encodedB, 
#'                           list( flds = PIVs, 
#'                                 types = rep("bi",length(PIVs)) ) )
#' diag_brl <- RL_diagnostics(fit_brl, prep_data$encodedA, prep_data$encodedB,
#'                            PIVs, PIVs_type, PIVs, true_pairs = prep_data$true_pairs,
#'                            FDP_estimation = TRUE, RL_method = "BRL", 
#'                            flds = PIVs, types = rep("bi",length(PIVs)),
#'                            maxIter4CV = 1, n_repeats = 1)
#' diag_brl # print(diag_brl)
#' print(diag_brl, threshold = 0.75)
#' plot(diag_brl, "scores")
#' plot(diag_brl, "distributions", threshold = 0.75)
#' plot(diag_brl, "FDP")
#' plot(diag_brl, "discrepancy")
RL_diagnostics <- function(fit, encodedA, encodedB, compare_vars, vars_type_cont,
                           PIVs, true_pairs = NULL, FDP_estimation = TRUE, ...) {

  arguments <- list(...)
  
  varsA  <- intersect(compare_vars, names(encodedA))
  varsB  <- intersect(compare_vars, names(encodedB))
  varsAB <- intersect(varsA, varsB)
  if (length(setdiff(compare_vars, union(varsA, varsB))) > 0) {
    stop("Some `compare_vars` are in neither `encodedA` nor `encodedB`.", call. = FALSE)
  }
  if (!all(compare_vars %in% names(vars_type_cont))) {
    stop("All `compare_vars` should be specified in `vars_type_cont`.", call. = FALSE)
  }
  
  if ("Delta" %in% names(fit) && all(c("i", "j", "x") %in% names(fit$Delta))) {
      fit$idxA <- fit$Delta$i
      fit$idxB <- fit$Delta$j
      fit$LinkScore <- fit$Delta$x
  }
  if (!all(c("idxA", "idxB", "LinkScore") %in% names(fit))) {
    stop("`fit` must contain `idxA`, `idxB`, `LinkScore` (or `Delta` with `i`, `j`, `x`).", call. = FALSE)
  }
  if (length(fit$idxA) == 0) {
    if (is.null(arguments$RL_method)){
      stop("The linkage contains no linked pair: nothing to diagnose.", call. = FALSE)
    } else {
      stop(sprintf("%s linkage contains no linked pair: nothing to diagnose.", arguments$RL_method), call. = FALSE)
    } 
  }

  if (!is.null(fit$LinkScore) && length(unique(fit$LinkScore))>1) {
    th <- seq(0.5, 0.99, by = 0.01)
  } else {
    th <- 0.5
    if (!is.null(fit$LinkScore)) {
      single_score <- fit$LinkScore[1]
      if (is.na(single_score)) {
        warning("No valid linked pair, `LinkScore` is NA.", call. = FALSE)
      } else if (single_score < th) {
        warning("No valid linked pair, `LinkScore` < 0.5.", call. = FALSE)
      }
    }
  }

  true_performance <- NULL
  if (!is.null(true_pairs)) {
    true_key <- paste(true_pairs[[1]], true_pairs[[2]], sep = "_")
    true_performance <- do.call(rbind, lapply(th, function(t) {
      if (!is.null(fit$LinkScore)) {
        keep <- fit$LinkScore > t
      } else {
        keep <- rep(TRUE, length(fit$idxA))
      }
      linked_key <- paste(fit$idxA[keep], fit$idxB[keep], sep = "_")
      tp <- length(intersect(linked_key, true_key))
      fp <- length(setdiff(linked_key, true_key))
      fn <- length(setdiff(true_key, linked_key))
      data.frame(xi = t,
                 FDP = if (tp + fp > 0) fp / (tp + fp) else NA_real_,
                 sensitivity = if (tp + fn > 0) tp / (tp + fn) else NA_real_)
    }))
  }
  
  true_agreement <- NULL
  if (!is.null(true_pairs) && length(varsAB) > 0) {
    true_agreement <- RL_agreement(encodedA, encodedB, varsAB, true_pairs)$agreements
  }

  discrepancy_measures <- do.call(rbind, lapply(th, function(t) {
    m <- .discrepancy_measures(fit, encodedA, encodedB, compare_vars, vars_type_cont, t)
    if (m$n_linked > 1) {
      t(data.frame(unlist(list(xi = t,
                               n_linked = m$n_linked,
                               agreementAB = unlist(m$agreement),
                               agreementLinks = true_agreement,
                               smdA = unlist(m$smdA),
                               smdB = unlist(m$smdB),
                               iouA = unlist(m$iouA),
                               iouB = unlist(m$iouB),
                               mmdA = m$mmdA,
                               mmdB = m$mmdB))))
    }}))
  if (is.null(discrepancy_measures)) {
    stop("Fewer than two linked pairs at every threshold: no diagnostics can be computed.", call. = FALSE)
  }
  discrepancy_measures <- as.data.frame(discrepancy_measures)
  rownames(discrepancy_measures) <- seq_len(nrow(discrepancy_measures))
  
  FDP_measures <- NULL
  if (isTRUE(FDP_estimation)){

    if (!"RL_method" %in% names(arguments)){
      stop("`RL_method` must be given when FDP_estimation is `TRUE`, see [compute_augmRL_FDP_synth()].", call. = FALSE)
    }

    arguments$synth_method <- if ("synth_method" %in% names(arguments)) arguments$synth_method else "arf"
    arguments$encodedA <- encodedA
    arguments$encodedB <- encodedB
    arguments$PIVs <- if ("PIVs" %in% names(arguments)) arguments$PIVs else PIVs

    cat("\nFDP estimation on augmented RL task")
    FDP_synth_res <- do.call(compute_augmRL_FDP_synth, arguments)
    cat("\nFDP estimation on original RL task")
    FDP_score_res <- do.call(compute_RL_FDP_score, arguments)

    FDP_measures <- .FDP_measures(FDP_score_res, FDP_synth_res, true_performance)
  }

  structure(
    list(
      compare_vars = compare_vars, LinkScore = fit$LinkScore,
      idxA = fit$idxA, idxB = fit$idxB, RL_method = arguments$RL_method,
      n_pairs = nrow(encodedA) * nrow(encodedB), true_performance = true_performance,
      true_agreement = true_agreement,
      A = encodedA[, varsA, drop = FALSE],
      B = encodedB[, varsB, drop = FALSE],
      gamma = fit$gamma, eta = fit$eta, alpha = fit$alpha, phi = fit$phi,
      FDP_measures = FDP_measures,
      discrepancy_measures = structure(discrepancy_measures, class = c("discrepancy_curves", "data.frame"))
    ),
    class = "RL_diagnostics"
  )
}

#' Print post-linkage diagnostics
#'
#' @param x An `RL_diagnostics` object, see [RL_diagnostics()].
#' @param threshold Numeric, the linkage decision rule at which the measures
#'   are displayed (default `0.5`).
#' @param ... Ignored.
#'
#' @return `x`, invisibly; called for its printing side effect.
#' @export
print.RL_diagnostics <- function(x, threshold = 0.5, ...) {
  cat("<RL_diagnostics>\n\n")
  
  if (is.null(x$LinkScore)){
    if (is.null(x$RL_method)){
      cat("  The linkage returned no LinkScore.\n")
    } else {
      cat(sprintf("  %s linkage returned no LinkScore.\n", x$RL_method))
    }
  }
  if (!is.null(x$LinkScore) && length(unique(x$LinkScore))==1) {
    if (is.null(x$RL_method)){
      cat(sprintf("  The linkage returned one unique LinkScore: %s.\n", round(x$LinkScore[1],3)))
    }
    cat(sprintf("  %s linkage returned one unique LinkScore: %s.\n", x$RL_method, round(x$LinkScore[1],3)))
  }
  
  dm <- x$discrepancy_measures
  if (min(abs(dm$xi - threshold)) > 1e-8) {
    cat(sprintf("  No pair linked at threshold %s.\n", threshold))
  }
  at_dm <- dm[which.min(abs(dm$xi - threshold)), , drop = FALSE]
  threshold <- at_dm$xi
  
  cat(sprintf("  Linked pairs  (score > %.2f):  %d\n\n", threshold, at_dm$n_linked))
  if (!is.null(x$true_performance)) {
    at <- x$true_performance[which.min(abs(x$true_performance$xi - threshold)), , drop = FALSE]
    cat(sprintf("  FDP           (true pairs):    %.3f\n", at$FDP))
    cat(sprintf("  Sensitivity   (true pairs):    %.3f\n\n", at$sensitivity))
  }
  
  if (!is.null(x$FDP_measures)) {
    valid <- !apply(apply(x$FDP_measures$curve, 2, is.na),1,any)
    first <- x$FDP_measures$curve[valid, ][1, ]
    idx <- which.min(abs(x$FDP_measures$curve$xi - threshold))
    middle <- x$FDP_measures$curve[idx, ]
    last <- x$FDP_measures$curve[valid, ][sum(valid), ]
    cat(sprintf("  FDP score estimate (RL task):            min threshold %.2f: FDP ~ %.3f\n%-*s     threshold %.2f: FDP ~ %.3f\n%-*s max threshold %.2f: FDP ~ %.3f\n",     first$xi, first$RL_FDP_score,     42, "", threshold, middle$RL_FDP_score    , 42, "", last$xi, last$RL_FDP_score     ))
    cat(sprintf("  FDP score estimate (augmented RL task):  min threshold %.2f: FDP ~ %.3f\n%-*s     threshold %.2f: FDP ~ %.3f\n%-*s max threshold %.2f: FDP ~ %.3f\n",     first$xi, first$augmRL_FDP_score, 42, "", threshold, middle$augmRL_FDP_score, 42, "", last$xi, last$augmRL_FDP_score ))
    cat(sprintf("  FDP synth estimate (augmented RL task):  min threshold %.2f: FDP ~ %.3f\n%-*s     threshold %.2f: FDP ~ %.3f\n%-*s max threshold %.2f: FDP ~ %.3f\n",     first$xi, first$augmRL_FDP_synth, 42, "", threshold, middle$augmRL_FDP_synth, 42, "", last$xi, last$augmRL_FDP_synth ))
    cat(sprintf("  The score-based estimator is valid when the linkage model is well calibrated to the data.\n  The synthetic-data estimator is valid when the augmented task is equivalent to the\n  original one, which requires links and non-links to have similar distributions; similar\n  score-based FDP estimates on the original and augmented tasks support this prerequisite.\n\n" ))
  }
  
  cat(sprintf("  Multivariate MMD  (linked A vs. A):  %.4f\n", at_dm$mmdA))
  cat(sprintf("  Multivariate MMD  (linked B vs. B):  %.4f\n\n", at_dm$mmdB))
  
  fmt <- function(z) formatC(z, format = "f", digits = 4)
  show_lines <- function(left, right, side) {
    cat(sprintf("  %s  %s:  %s\n", formatC(left, width = -max(nchar(left))), side, right), sep = "")
  }
  show_block <- function(prefix, label, side) {
    max_levels = 3
    v <- unlist(at_dm[grep(paste0("^", prefix, "\\."), names(at_dm))])
    if (length(v) == 0) return (invisible(NULL))
    nm <- sub(paste0("^", prefix, "\\."), "", names(v))
    var <- sub("\\..*$", "", nm)
    left <- right <- character(0)
    for (vr in unique(var)) {
      vals <- v[var == vr]
      if (length(vals) == 1 && !grepl("\\.", nm[var == vr])) {
        left  <- c(left, sprintf("%s %s", label, vr))
        right <- c(right, fmt(vals))
        next
      }
      lev <- sub(paste0("^.*\\.", vr, "_"), "", nm[var == vr])                           
      lev_txt <- lev
      val_txt <- fmt(vals)
      if (length(vals) > max_levels) {
        keep    <- sort(order(abs(vals), decreasing = TRUE)[seq_len(max_levels)])
        lev_txt <- c(lev[keep], "...")
        val_txt <- c(fmt(vals[keep]), "...")
      }
      left <- c(left, sprintf("%s %s values: {%s}", label, vr, paste(lev_txt, collapse = ", ")))
      right <- c(right, paste(val_txt, collapse = ", "))
    }
    show_lines(left, right, side)
  }
  
  show_block("iouA", "IoU support", "(linked A vs. A)")
  show_block("iouB", "IoU support", "(linked B vs. B)")
  cat("\n")
  show_block("agreementAB", "Agreement", "(linked A vs. linked B)")
  if (!is.null(x$true_agreement)) {
    show_block("agreementLinks", "Agreement", "(true pairs)") 
  }
  cat("\n")
  show_block("smdA", "SMD", "(linked A vs. A)")
  show_block("smdB", "SMD", "(linked B vs. B)")
  invisible(x)
}

#' Plot post-linkage diagnostics
#'
#' @param x An `RL_diagnostics` object, see [RL_diagnostics()].
#' @param type One of `"scores"` (linkage score histogram),
#'   `"distributions"` (linked subset vs. data sources, per variable, at
#'   `threshold`), `"convergence"` (StEM trace plots, see
#'   [plot_StEM_convergence()]), `"FDP"` (FDP estimates by threshold, only if
#'   `FDP_estimation = TRUE` was used) or `"discrepancy"` (MMD, SMD, IoU and
#'   agreement by threshold, see [plot.discrepancy_curves()]).
#' @param threshold Numeric, the linkage decision rule defining the linked
#'   subset for `type = "distributions"` (ignored for methods without scores).
#' @param ... Passed on to the underlying plotting helper.
#'
#' @return `x`, invisibly; called for its plotting side effect.
#' @export
plot.RL_diagnostics <- function(x, type, threshold = NULL, ...) {
  types <- c("scores", 
             "distributions", 
             "convergence", 
             "FDP", 
             "discrepancy")
  if (!type %in% types) {
    stop("`type` must be one of: ", paste(types, collapse = ", "), ".", call. = FALSE)
  }
  if (type == "scores") {
    if (is.null(x$LinkScore)) {
      stop("[plot(`RL_diagnostics`, `scores`)] is only available for `RL_method` that return `LinkScore`.", call. = FALSE)
    }
    op <- graphics::par(mfrow = c(1, 1))
    on.exit(graphics::par(op))
    plot_linkage_scores(x$n_pairs, x$LinkScore[x$LinkScore > 0])
  } else if (type == "distributions") {
    if (is.null(x$LinkScore)) {
      linked <- rep(TRUE, length(x$idxA))
      threshold <- NA
    } else {
      if (is.null(threshold)) {
        stop("`threshold` is needed to define the linked subset for `type = \"distributions\"`.", call. = FALSE)
      }
      linked <- x$LinkScore > threshold
      if (sum(linked) == 0) {
        warning(sprintf("No pair linked at threshold %s, using 0.5.", threshold), call. = FALSE)
        threshold <- 0.5
        linked <- x$LinkScore > threshold
      }
    }
    linkedA <- x$A[x$idxA[linked], , drop = FALSE]
    linkedB <- x$B[x$idxB[linked], , drop = FALSE]
    varsA <- colnames(x$A)
    if (length(varsA) > 0) {
      op <- graphics::par(mfrow = grDevices::n2mfrow(length(varsA)))
      on.exit(graphics::par(op))
      plot_distributions(list(A = x$A, linkedA = linkedA), varsA, threshold, ...)
    }
    varsB <- colnames(x$B)
    if (length(varsB) > 0) {
      op <- graphics::par(mfrow = grDevices::n2mfrow(length(varsB)))
      on.exit(graphics::par(op))
      plot_distributions(list(B = x$B, linkedB = linkedB), varsB, threshold, ...)
    }
  } else if (type == "convergence") {
    if (is.null(x$gamma)) {
      stop("No StEM chains stored on this object (the `fit` passed to `RL_diagnostics()` had no `gamma`/`eta`/`alpha`/`phi`).", call. = FALSE)
    }
    op <- graphics::par(mfrow = c(1,1))
    on.exit(graphics::par(op))
    plot_StEM_convergence(list(gamma = x$gamma, eta = x$eta, alpha = x$alpha, phi = x$phi), ...)
  } else if (type == "FDP") {
    if (!is.null(x$FDP_measures)) {
      if (length(x$discrepancy_measures$xi)<=1) {
        stop("[plot(`RL_diagnostics`, `FDP`)] is only available for `RL_method` that return `LinkScore`. Point estimates are available with `print`.", call. = FALSE)
      }
      op <- graphics::par(mfrow = c(1,1))
      on.exit(graphics::par(op))
      plot(x$FDP_measures, ...)
    } else {
      warning("No FDP measures to plot: use `FDP_estimation = TRUE` in `RL_diagnostics()`.", call. = FALSE)
    }
  } else if (type == "discrepancy") {
    if (!is.null(x$discrepancy_measures)) {
      if (length(x$discrepancy_measures$xi)<=1) {
        stop("[plot(`RL_diagnostics`, `discrepancy`)] is only available for `RL_method` that return `LinkScore`. Point estimates are available with `print`.", call. = FALSE)
      }
      op <- graphics::par(mfrow = c(1,1))
      on.exit(graphics::par(op))
      plot(x$discrepancy_measures, ...)
    } else {
      warning("No discrepancy measures to plot.", call. = FALSE)
    }
  }
  invisible(x)
}

# FDP estimation measures for given linked subsets (corresponding to different
# thresholds). FDP based on model scores on the original record linkage task, 
# FDP based on model scores and FDP based on synthetic data on the augmented
# record linkage task,
.FDP_measures <- function(FDP_score_res, FDP_synth_res, true_performance) {
  grid <- seq(0.5, 0.99, by = 0.01)
  # FDP score-based estimate on the record linkage task
  curve <- data.frame(xi = grid,
                      RL_n_linked_score = colMeans(FDP_score_res$Linked_pairs, na.rm = TRUE),
                      RL_FDP_score      = colMeans(FDP_score_res$FDP_score_estimator, na.rm = TRUE))
  curve[!is.na(curve$RL_n_linked_score) & curve$RL_n_linked_score == 0, "RL_FDP_score"] <- NA
  # FDP score-based estimate on the augmented task
  curve$augmRL_n_linked_score <- colMeans(FDP_synth_res$Linked_pairs, na.rm = TRUE)
  curve$augmRL_FDP_score      <- colMeans(FDP_synth_res$FDP_score_estimator, na.rm = TRUE)
  curve[!is.na(curve$augmRL_n_linked_score) & curve$augmRL_n_linked_score == 0, "augmRL_FDP_score"] <- NA
  # FDP synthetic-data estimate on the augmented task
  curve$augmRL_n_linked_synth <- colMeans(FDP_synth_res$Linked_pairs_augm, na.rm = TRUE)
  curve$augmRL_FDP_synth      <- colMeans(FDP_synth_res$FDP_synth_estimator, na.rm = TRUE)
  curve[!is.na(curve$augmRL_n_linked_synth) & curve$augmRL_n_linked_synth == 0, "augmRL_FDP_synth"] <- NA
  structure(list(curve = curve, true_performance = true_performance), class = "FDP_curves")
}

#' Plot FDP estimates against the linkage decision threshold
#'
#' @param x An object of class `"FDP_curves"`.
#' @param ... Ignored.
#'
#' @return `x`, invisibly; called for its plotting side effect.
#' @export
plot.FDP_curves <- function(x, ...) {
  op <- graphics::par(mfrow = c(1, 1), mar = c(4, 4, 2, 4))
  on.exit(graphics::par(op))
  cv <- x$curve
  ymax <- max(0.5, cv$RL_FDP_score, cv$augmRL_FDP_score, cv$augmRL_FDP_synth,
              if (!is.null(x$true_performance)) x$true_performance$FDP, na.rm = TRUE)
  graphics::plot(cv$xi, cv$RL_FDP_score, type = "l", lwd = 2, ylim = c(0, ymax),
                 xlab = "decision rule threshold", ylab = "estimated FDP",
                 main = "FDP versus linkage decision rule")
  graphics::lines(cv$xi, cv$augmRL_FDP_synth, lwd = 2, lty = 2)
  graphics::lines(cv$xi, cv$augmRL_FDP_score, lwd = 2, lty = 3)
  has_truth <- !is.null(x$true_performance)
  if (has_truth) graphics::lines(x$true_performance$xi, x$true_performance$FDP, col = "grey50")
  graphics::par(new = TRUE)
  graphics::plot(cv$xi, cv$RL_n_linked_score, type = "l", col = "firebrick", axes = FALSE,
                 xlab = "", ylab = "",
                 ylim = c(0, max(cv$RL_n_linked_score, cv$augmRL_n_linked_synth, cv$augmRL_n_linked_score, na.rm = TRUE)))
  graphics::lines(cv$xi, cv$augmRL_n_linked_synth, col = "red")
  graphics::lines(cv$xi, cv$augmRL_n_linked_score, col = "orange")
  graphics::axis(4, col = "firebrick", col.axis = "firebrick")
  graphics::mtext("number of links", side = 4, line = 2.5, col = "firebrick")
  lab <- c("FDP score (RL)", "FDP synth (augm. RL)", "FDP score (augm. RL)",
           "links (RL)", "real links (augm. RL)", "links (augm. RL)")
  lty <- c(1, 2, 3, 1, 1, 1)
  col <- c("black", "black", "black", "firebrick", "red", "orange")
  if (has_truth) { lab <- c(lab, "true FDP"); lty <- c(lty, 1); col <- c(col, "grey50") }
  graphics::legend("topright", lab, lty = lty, col = col, bty = "n", cex = 0.8)
  invisible(x)
}

