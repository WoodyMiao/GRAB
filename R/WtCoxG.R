#' WtCoxG method in GRAB package
#'
#' WtCoxG is an accurate, powerful, and computationally efficient Cox-based approach to perform genome-wide time-to-event data analyses in study cohorts with case ascertainment.
#'
#' @details
#' Additional arguments in \code{GRAB.NullModel()}:
#' \itemize{
#'   \item \code{RefAfFile}: A character string specifying a reference allele frequency file, which is a csv file (with a header) and includes columns of \code{CHROM}, \code{POS}, \code{ID}, \code{REF}, \code{ALT}, \code{AF_ref}, and \code{AN_ref}.
#'   \item \code{OutputFile}: A character string specifying the output file name.
#'   \item \code{SampleIDColumn}: A character string specifying the column name in the input data that contains sample IDs.
#'   \item \code{SurvTimeColumn}: A character string specifying the column name in the input data that contains survival time information.
#'   \item \code{IndicatorColumn}: A character string specifying the column name in the input data that indicates case-control status (should be 0 for controls and 1 for cases).
#' }
#'
#' Additional arguments in list \code{control} in \code{GRAB.NullModel()}:
#' \itemize{
#'   \item \code{RefPrevalence}: A numeric value specifying the population-level disease prevalence used for weighting in the analysis.
#'   \item \code{SNPnum}: Minimum number of SNPs. Default is 1e4.
#' }
#'
#' Additional arguments in list \code{control} in \code{GRAB.Marker()}:
#' \itemize{
#'   \item \code{cutoff}: A numeric value specifying the batch effect p-value cutoff for method selection of an assocition test. Default is 0.1.
#' }
#' 
#' @examples
#' # Step0&1: fit a null model and estimate parameters according to batch effect p-values
#' PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData <- data.table::fread(PhenoFile, header = T)
#' SparseGRMFile <- system.file("SparseGRM", "SparseGRM.txt", package = "GRAB")
#' 
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' RefAfFile <- system.file("extdata", "simuRefMAF.txt", package = "GRAB")
#' RefPrevalence <- 0.1 # population-level disease prevalence
#' 
#' OutputDir <- system.file("results", package = "GRAB")
#' OutputStep1 <- paste0(OutputDir, "/wtCoxG_step1_out.txt")
#' OutputStep2 <- paste0(OutputDir, "/wtCoxG_step2_out.txt")
#'
#' obj.WtCoxG <- GRAB.NullModel(
#'   formula = survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER,
#'   data = PhenoData,
#'   subjData = PhenoData$IID,
#'   method = "WtCoxG",
#'   traitType = "time-to-event",
#'   GenoFile = GenoFile,
#'   SparseGRMFile = SparseGRMFile,
#'   control = list(AlleleOrder = "ref-first", AllMarkers = T, RefPrevalence = RefPrevalence),
#'   RefAfFile = RefAfFile,
#'   OutputFile = OutputStep1,
#'   SampleIDColumn = "IID",
#'   SurvTimeColumn = "SurvTime",
#'   IndicatorColumn = "SurvEvent"
#' )
#'
#' resultStep1 <- data.table::fread(OutputStep1)
#' 
#' # Step2: conduct association testing
#' GRAB.Marker(
#'   objNull = obj.WtCoxG,
#'   GenoFile = GenoFile,
#'   OutputFile = OutputStep2,
#'   control = list(AlleleOrder = "ref-first", AllMarkers = T, cutoff = 0.1, nMarkersEachChunk = 5000)
#' )
#' 
#' resultStep2 <- data.table::fread(OutputStep2)
#' @export
GRAB.WtCoxG <- function() {
  print("Check ?GRAB.WtCoxG for more details about 'WtCoxG' method.")
}


# check the control list in null model fitting for SPACox method
checkControl.NullModel.WtCoxG <- function(control, traitType, ...) {
  if (!traitType %in% c("time-to-event", "binary", "Residual")) {
    stop("For 'WtCoxG' method, only traitType of 'time-to-event', 'binary' or 'Residual' is supported.")
  }

  default.control <- list(RefPrevalence = 0)

  control <- updateControl(control, default.control)

  if (control$RefPrevalence <= 0 | control$RefPrevalence >= 0.5) {
    stop("control$Refprevalence is required and should be between (0, 0.5).")
  }

  return(control)
}


getWeight.WtCoxG <- function(Indicator, RefPrevalence) {
  # BWJ (2023-08-08): Check NA later.
  # if(any(!unique(Indicator) %in% c(0, 1, NA)))
  #   stop("The value of Indicator should be 0, 1, or NA")

  if (any(!unique(Indicator) %in% c(0, 1))) {
    stop("The value of Indicator should be 0 or 1.")
  }

  sumOnes <- sum(Indicator)
  sumZeros <- sum(1 - Indicator)
  ratio <- sumOnes / sumZeros

  weight <- ifelse(Indicator == 1, 1,
    (1 - RefPrevalence) / RefPrevalence * ratio
  )
  return(weight)
}


# fit null model using WtCoxG method
fitNullModel.WtCoxG <- function(response, designMat, subjData,
                                control = list(OutlierRatio = 1.5), ...) {
  if (!class(response) %in% c("Surv", "Residual")) { # later support binary trait and Residual
    stop("For WtCoxG, the response variable should be of class 'Surv' or 'Residual'.")
  }

  if (class(response) == "Surv") {
    formula <- response ~ designMat

    Indicator <- as.matrix(response)[, "status"]
    RefPrevalence <- control$RefPrevalence

    weight <- getWeight.WtCoxG(Indicator, RefPrevalence)

    obj.coxph <- survival::coxph(formula, x = T, weight = weight, robust = T, ...)

    y <- obj.coxph$y
    yVec <- y[, ncol(y)] # status

    mresid <- obj.coxph$residuals
    # Cova = obj.coxph$x
    Cova <- designMat
  } else {
    stop("We only support 'time-to-event' trait for WtCoxG by 2023-08-08.")
  }

  # if(class(response) == "Residual")
  # {
  #   yVec = mresid = response
  #   Cova = designMat
  # }

  # var.resid = var(mresid, na.rm = T)
  ## outliers or not depending on the residuals (0.25%-1.5IQR, 0.75%+1.5IQR)
  q25 <- quantile(mresid, 0.25, na.rm = T)
  q75 <- quantile(mresid, 0.75, na.rm = T)
  IQR <- q75 - q25

  # r.outlier = 1.5    # put this to the control argument later
  r.outlier <- ifelse(is.null(control$OutlierRatio), 1.5, control$OutlierRatio) # put this to the control argument later
  cutoff <- c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
  posOutlier <- which(mresid < cutoff[1] | mresid > cutoff[2])
  cat("The r.outlier is:", r.outlier, "\n")
  while (length(posOutlier) == 0) {
    r.outlier <- r.outlier * 0.8
    cutoff <- c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
    posOutlier <- which(mresid < cutoff[1] | mresid > cutoff[2])
    cat("The curent outlier ratio is:", r.outlier, "\n")
    cat("The number of outlier is:", length(posOutlier), "\n")
  }

  # The original code is from SPAmix in which multiple residuals were analysis simultaneously
  # posValue = which(!is.na(mresid))
  # posNonOutlier = setdiff(posValue, posOutlier)
  posValue <- 1:length(mresid)
  posNonOutlier <- setdiff(posValue, posOutlier)

  cat("The outlier of residuals will be passed to SPA analysis.\n")
  cat("Cutoffs to define residuals:\t", signif(cutoff, 2), "\n")
  cat("Totally, ", length(posOutlier), "/", length(posValue), " are defined as outliers.\n")

  if (length(posOutlier) == 0) {
    stop("No outlier is observed. SPA is not required in this case.")
  }

  # "-1" is to convert R style (index starting from 1) to C++ style (index starting from 0)
  outLierList <- list(
    posOutlier = posOutlier - 1,
    posNonOutlier = posNonOutlier - 1,
    resid = mresid,
    resid2 = mresid^2,
    residOutlier = mresid[posOutlier],
    residNonOutlier = mresid[posNonOutlier],
    resid2NonOutlier = mresid[posNonOutlier]^2
  )

  re <- list(
    mresid = mresid, Cova = Cova, yVec = yVec, weight = weight,
    RefPrevalence = RefPrevalence,
    N = length(mresid),
    outLierList = outLierList
  )

  class(re) <- "WtCoxG_NULL_Model"
  return(re)
}


#' Quality control to check batch effect between study cohort and reference population.
#'
#' This function performs quality control to test for the batch effect between a study cohort and a reference population. And fit a weighted null model.
#'
#' @param objNull a \code{WtCoxG_NULL_Model} object, which is the output of \code{\link{GRAB.NullModel}}.
#' @param data a data.frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data.frame), containing the variables in formula. Neither a matrix nor an array will be accepted.
#' @param GenoFile A character string of the genotype file. See Details section for more details.
#' @param GenoFileIndex Additional index file(s) corresponding to GenoFile. See Details section for more details.
#' @param SparseGRMFile a path to file of output to be passed to \code{\link{GRAB.NullModel}}.
#' @param RefAFfile A character string of the reference file. The reference file must be a \code{txt} file (header required) including at least 7 columns: \code{CHROM}, \code{POS}, \code{ID}, \code{REF}, \code{ALT}, \code{AF_ref}, \code{AN_ref}.
#' @param OutputFile A character string of the output file name. The output file will be a \code{txt} file.
#' @param IndicatorColumn A character string of the column name in \code{data} that indicates the case-control status. The value should be 0 for controls and 1 for cases.
#' @param SurvTimeColumn A character string of the column name in \code{data} that indicates the survival time.
#' @param SampleIDColumn A character string of the column name in \code{data} that indicates the sample ID.
#' @return A dataframe of marker info and reference MAF.
#'
#' @export
TestforBatchEffect <- function(
    objNull,
    data,
    GenoFile, # a character of file names of genotype files
    GenoFileIndex = NULL, # additional index file(s) corresponding to GenoFile
    SparseGRMFile = NULL, # sparse genotype relatedness matrix
    RefAfFile, # header should include c("CHROM", "POS", "ID", "REF", "ALT", "AF_ref","AN_ref")
    OutputFile,
    IndicatorColumn,
    SurvTimeColumn,
    SampleIDColumn) {

  if (!is.null(SparseGRMFile)) {
    sparseGRM <- data.table::fread(SparseGRMFile)
  } else {
    sparseGRM <- NULL
  }

  control <- objNull$control
  RefPrevalence <- control$RefPrevalence # refernce population prevalence, the proportion of indicator == 1.

  if (is.null(control$SNPnum)) {
    SNPnum <- 1e4 # default least number of SNPs is 1e4
  } else {
    SNPnum <- control$SNPnum
  }

  posCol <- which(colnames(data) == IndicatorColumn)
  colnames(data)[posCol] <- "Indicator"

  posCol <- which(colnames(data) == SurvTimeColumn)
  colnames(data)[posCol] <- "SurvTime"

  posCol <- which(colnames(data) == SampleIDColumn)
  colnames(data)[posCol] <- "SampleID"


  # step1: quality control--------------------------------------------------------
  ## reference genoInfo----------------------------------------------------------
  refGenoInfo <- data.table::fread(RefAfFile) %>% as_tibble()

  # check if there are 7 columns in RefAfFile
  for (colname in c("CHROM", "POS", "ID", "REF", "ALT", "AF_ref", "AN_ref")) {
    if (!colname %in% colnames(refGenoInfo)) {
      stop(paste0(colname, " is missing in RefAfFile!"))
    }
  }

  ## merge sample genoInfo and ref genoInfo--------------------------------------
  GenoInfo.ctrl <- GRAB.getGenoInfo(
    GenoFile = GenoFile,
    GenoFileIndex = GenoFileIndex,
    SampleIDs = with(data, SampleID[Indicator == 0]), # MAF in cases
    control = control
  ) %>%
    rename(mu0 = altFreq, mr0 = missingRate) %>%
    select(mu0, mr0)

  if (nrow(GenoInfo.ctrl) < SNPnum) {
    stop("The number of genetic variants < ", SNPnum)
  }

  GenoInfo <- GRAB.getGenoInfo(
    GenoFile = GenoFile,
    GenoFileIndex = GenoFileIndex,
    SampleIDs = with(data, SampleID[Indicator == 1]), # MAF in controls
    control = control
  ) %>%
    rename(mu1 = altFreq, mr1 = missingRate) %>%
    cbind(., GenoInfo.ctrl) %>%
    as_tibble() %>%
    mutate(RA = paste0(pmin(REF, ALT), pmax(REF, ALT))) %>%
    mutate(index = 1:n())

  mergeGenoInfo <- refGenoInfo %>%
    mutate(RA = paste0(pmin(REF, ALT), pmax(REF, ALT))) %>%
    merge(., GenoInfo, by = c("CHROM", "POS", "RA"), all.y = T, sort = F) %>%
    rename(REF = REF.y, ALT = ALT.y, ID = ID.y) %>%
    mutate(AF_ref = ifelse(REF == REF.x, AF_ref, 1 - AF_ref)) %>%
    select(-REF.x, -ALT.x, -ID.x, -RA) %>%
    arrange(index) %>%
    mutate(
      n1 = sum(data$Indicator) * (1 - mr1),
      n0 = sum(1 - data$Indicator) * (1 - mr0),
      mu.int = 0.5 * mu1 + 0.5 * mu0,
      mu.int = ifelse(mu.int > 0.5, 1 - mu.int, mu.int)
    )

  #### calculate batch effect p-value for each genetic variant------------------------------------
  w1 <- objNull$weight / (2 * sum(objNull$weight))
  names(w1) <- data$SampleID
  meanR <- mean(objNull$mresid)
  R_tilde <- objNull$mresid - meanR
  names(R_tilde) <- data$SampleID

  if (!is.null(sparseGRM)) {
    sparseGRM <- sparseGRM %>% mutate(
      cov = Value * w1[as.character(ID1)] * w1[as.character(ID2)],
      cov_R = Value * R_tilde[as.character(ID1)] * R_tilde[as.character(ID2)]
    )
    var.ratio.w0 <- (sum(sparseGRM$cov) + 1 / (2 * mergeGenoInfo$AN_ref)) / (sum(w1^2) + 1 / (2 * mergeGenoInfo$AN_ref))
    var.ratio.int <- sum(sparseGRM$cov_R) / sum(R_tilde^2)
  } else {
    var.ratio.w0 <- var.ratio.int <- 1
  }

  mergeGenoInfo <- mergeGenoInfo %>%
    mutate(
      var.ratio.w0 = var.ratio.w0,
      var.ratio.int = var.ratio.int
    )

  pvalue_bat <- lapply(1:nrow(mergeGenoInfo), function(ind) {
    p.test <- Batcheffect.Test(
      n0 = mergeGenoInfo$n0[ind],
      n1 = mergeGenoInfo$n1[ind],
      n.ext = mergeGenoInfo$AN_ref[ind] / 2,
      maf0 = mergeGenoInfo$mu0[ind],
      maf1 = mergeGenoInfo$mu1[ind],
      maf.ext = mergeGenoInfo$AF_ref[ind],
      pop.prev = RefPrevalence,
      var.ratio = mergeGenoInfo$var.ratio.w0[ind]
    )
  }) %>% unlist()

  mergeGenoInfo <- mergeGenoInfo %>% mutate(pvalue_bat)
  rm(pvalue_bat)

  #### estimate unknown parameters according to batch effect p-values---------------------------------
  cat("Estimate TPR and sigma2--------------\n")
  maf.group <- c(seq(0, 0.4, 0.05), max(mergeGenoInfo$mu.int))
  mergeGenoInfo <- lapply(1:(length(maf.group) - 1), function(i) {
    cat(i, "\n")

    ## assume that genotypes with MAF in [ maf.group[i] , maf.group[i+1]] have the same mixture distribution
    mergeGenoInfo_1 <- mergeGenoInfo %>% filter(mu.int > maf.group[i] & mu.int <= maf.group[i + 1])

    ## using batcheffect p-values with MAF in [maf.group[i]-0.1 , maf.group[i+1]+0.1] to estimate parameters
    mergeGenoInfo_2 <- mergeGenoInfo %>%
      filter(mu.int >= max(maf.group[i] - 0.1, 0) & mu.int < min(1, maf.group[i + 1] + 0.1))

    mu <- (maf.group[i] + maf.group[i + 1]) / 2

    n.ext <- mean(na.omit(mergeGenoInfo_1$AN_ref)[1]) / 2
    var_mu_ext <- mu * (1 - mu) / (2 * n.ext)

    var_Sbat <- ifelse(is.null(sparseGRM), sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext,
      na.omit(mergeGenoInfo$var.ratio.w0)[1] * (sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext)
    )

    obj <- fun.est.param(
      vec_p_bat = mergeGenoInfo_2$pvalue_bat,
      vec_var_Sbat = var_Sbat
    )
    TPR <- obj[1]
    sigma2 <- obj[2]

    w.ext <- optim(
      par = 0.5, method = "L-BFGS-B", lower = 0, upper = 1,
      fn = fun.optimalWeight,
      pop.prev = RefPrevalence,
      y = data$Indicator,
      R = objNull$mresid,
      w = objNull$weight,
      mu = mu,
      N = nrow(data),
      n.ext = n.ext,
      sigma2 = obj$sigma2,
      TPR = obj$TPR
    )$par[1]

    if (is.null(sparseGRM)) {
      var.ratio.ext <- 1
    } else {
      R_tilde_w <- objNull$mresid - mean(objNull$mresid) * w.ext
      names(R_tilde_w) <- data$SampleID
      sparseGRM <- sparseGRM %>% mutate(cov_Rext = Value * R_tilde_w[as.character(ID1)] * R_tilde_w[as.character(ID2)])
      var.ratio.ext <- (sum(sparseGRM$cov_Rext) + w.ext^2 * sum(objNull$mresid)^2 / n.ext) / (sum(R_tilde_w^2) + w.ext^2 * sum(objNull$mresid)^2 / n.ext)
    }

    mergeGenoInfo_1 <- mergeGenoInfo_1 %>% cbind(., TPR, sigma2, w.ext, var.ratio.ext)
  }) %>%
    do.call("rbind", .) %>%
    as_tibble() %>%
    arrange(index) %>%
    select(-index)

  #### output-----------------------------------------------------------
  data.table::fwrite(
    mergeGenoInfo,
    OutputFile,
    row.names = F, col.names = T, quote = F, sep = "\t"
  )
  return(mergeGenoInfo)
}


## function to test for batch effect--------------------------------------------
#' Test for batch effect
#'
#' This function test for the allele frequency difference between the study cohort and the external datasets.
#'
#' @param n0 A numeric. The sample size of cases in the study cohort
#' @param n1 A numeric. The sample size of controls in the study cohort
#' @param n.ext A numeric. The sample size of external datasets
#' @param maf0 A numeric. The MAF of the cases.
#' @param maf1 A numeric. The MAF of the controls
#' @param maf.ext A numeric. The MAF of the external datasets.
#' @param pop.prev A numeric. The population prevalence of the disease.
#' @param var.ratio A numeric. The variance ratio calculated by sparseGRM.
#' @return A numeric of batch effect p-value
Batcheffect.Test <- function(n0, # number of controls
                             n1, # number of cases
                             n.ext, # number of external dataset
                             maf0, # estimated MAF in controls
                             maf1, # estimated MAF in cases
                             maf.ext, # estimated MAF in external dataset
                             pop.prev,
                             var.ratio = 1) {
  er <- n1 / (n1 + n0)
  w0 <- (1 - pop.prev) / pop.prev / ((1 - er) / er)
  w1 <- 1

  weight.maf <- sum(maf0 * w0 * n0 + maf1 * w1 * n1) / sum(w0 * n0 + w1 * n1) ## weighted mean of genotypes
  est.maf <- sum(maf0 * w0 * n0 + maf1 * w1 * n1 + maf.ext * n.ext * w0) / sum(n1 * w1 + n0 * w0 + n.ext * w0) ## MAF estimates

  v <- ((n1 * w1^2 + n0 * w0^2) / (2 * (n1 * w1 + n0 * w0)^2) + 1 / (2 * n.ext)) * est.maf * (1 - est.maf) ## variance of test statistics
  z <- (weight.maf - maf.ext) / sqrt(v) ## standardized statistics
  z.adj <- z / sqrt(var.ratio) ## adjusted statistics by variance ratio
  p <- 2 * pnorm(-abs(z.adj), lower.tail = TRUE)

  return(p)
}

# estimate TPR and sigma2----------------------------------------------------
fun.est.param <- function(vec_p_bat,
                          vec_var_Sbat,
                          vec_cutoff = seq(0.01, 0.4, 0.1)) {
  ########  step1: the proportion of p_bat>cutoff
  vec_p_deno <- lapply(vec_cutoff, function(p_cut) {
    p_deno <- mean(na.omit(vec_p_bat > p_cut))
  }) %>% unlist()

  ######## optimization function
  opti_fun <- function(var_Sbat, vec_p_deno, par) {
    diff <- lapply(1:length(vec_cutoff), function(j) {
      p_cut <- vec_cutoff[j]
      lb <- -qnorm(1 - p_cut / 2) * sqrt(var_Sbat)
      ub <- qnorm(1 - p_cut / 2) * sqrt(var_Sbat)
      p_deno <- vec_p_deno[j]

      c <- pnorm(ub, 0, sqrt(var_Sbat + par[2]), log.p = T)
      d <- pnorm(lb, 0, sqrt(var_Sbat + par[2]), log.p = T)

      pro.cut <- par[1] * (exp(d) * (exp(c - d) - 1)) + (1 - par[1]) * (1 - p_cut)
      t <- ((p_deno - pro.cut) / p_deno)^2
    }) %>% do.call("sum", .)

    return(diff)
  }

  ####### estimate TPR and sigma2 for each SNP
  var.diff <- lapply(1:length(vec_var_Sbat), function(i) {
    if (i %% 100 == 0) cat(i, "\n")

    obj <- optim(
      par = c(0.01, 0.01),
      # ,method = "SANN"
      # , lower = 0, upper = 1
      fn = opti_fun,
      vec_p_deno = vec_p_deno,
      var_Sbat = vec_var_Sbat[i]
    )
    TPR <- min(1, max(0, obj$par[1]))
    sigma2 <- min(1, max(0, obj$par[2]))
    return(cbind(TPR, sigma2))
  }) %>%
    do.call("rbind", .) %>%
    as_tibble()

  return(var.diff)
}

# optimal weight for mu_ext -----------------------------------------------
fun.optimalWeight <- function(par, pop.prev, R, y, mu1, w, mu, N, n.ext, sigma2, TPR) {
  b <- par[1]

  p.fun <- function(b, pop.prev, R, y, mu1, mu, w, N, n.ext, sigma2, TPR) {
    meanR <- mean(R)
    sumR <- sum(R)

    mu0 <- mu
    mu.pop <- mu1 * pop.prev + mu0 * (1 - pop.prev)

    mu.i <- ifelse(y == 1, 2 * mu1, 2 * mu0)

    S <- sum((R - (1 - b) * meanR) * mu.i) - sumR * 2 * b * mu.pop

    w1 <- w / (2 * sum(w))
    mu <- mean(mu.i) / 2

    var_mu_ext <- mu * (1 - mu) / (2 * n.ext)
    var_Sbat <- sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext

    p_cut <- 0.1
    lb <- -qnorm(1 - p_cut / 2) * sqrt(var_Sbat)
    ub <- qnorm(1 - p_cut / 2) * sqrt(var_Sbat)
    c <- pnorm(ub, 0, sqrt(var_Sbat + sigma2), log.p = T)
    d <- pnorm(lb, 0, sqrt(var_Sbat + sigma2), log.p = T)
    p_deno <- TPR * (exp(d) * (exp(c - d) - 1)) + (1 - TPR) * (1 - p_cut)

    ## sigma2=0
    var.int <- sum((R - (1 - b) * meanR)^2) * 2 * mu * (1 - mu)
    var_S <- var.int + 4 * b^2 * sumR^2 * var_mu_ext
    cov_Sbat_S <- sum(w1 * (R - (1 - b) * meanR)) * 2 * mu * (1 - mu) + 2 * b * sumR * var_mu_ext
    VAR <- matrix(c(var_S, cov_Sbat_S, cov_Sbat_S, var_Sbat), nrow = 2)
    p0 <- max(0, mvtnorm::pmvnorm(lower = c(-Inf, lb), upper = c(-abs(S), ub), mean = c(0, 0), sigma = VAR))

    ## sigma2!=0
    var_S1 <- var.int + 4 * b^2 * sumR^2 * (var_mu_ext + sigma2)
    cov_Sbat_S1 <- sum(w1 * (R - (1 - b) * meanR)) * 2 * mu * (1 - mu) + 2 * b * sumR * (var_mu_ext + sigma2)
    var_Sbat1 <- var_Sbat + sigma2
    VAR1 <- matrix(c(var_S1, cov_Sbat_S1, cov_Sbat_S1, var_Sbat1), nrow = 2)
    p1 <- max(0, mvtnorm::pmvnorm(lower = c(-Inf, lb), upper = c(-abs(S), ub), mean = c(0, 0), sigma = VAR1))

    p.con <- 2 * (TPR * p1 + (1 - TPR) * p0) / p_deno
    # diff = -log10(p.con)+log10(5e-8)
    diff <- -log10(p.con / 5e-8)

    return(diff)
  }

  mu1 <- uniroot(p.fun,
    lower = mu, upper = 1,
    b = b, pop.prev = pop.prev, mu = mu,
    R = R, y = y, w = w, N = N, n.ext = n.ext, sigma2 = sigma2, TPR = TPR
  )$root

  return(mu1)
}


mainMarker.WtCoxG.inR <- function(genoType, markerInfo, control, objNull) {
  mergeGenoInfo <- objNull$mergeGenoInfo
  RefPrevalence <- control$RefPrevalence
  cutoff <- if (!is.null(control$cutoff)) control$cutoff else 0.1

  R <- objNull$mresid
  w <- objNull$weight
  control <- checkControl.ReadGeno(control)
  G <- getGenoInCPP(genoType, markerInfo, length(w), control$ImputeMethod)

  ## GWAS analysis ----------------------------------------------------------------------
  GWAS <- lapply(1:ncol(G), function(i) {
    if (i %% 1000 == 0) cat("Complete ", i, "/", ncol(G), "\n")

    g <- G[, i]
    mu.ext <- mergeGenoInfo$AF_ref[i]
    n.ext <- mergeGenoInfo$AN_ref[i] / 2
    TPR <- mergeGenoInfo$TPR[i]
    sigma2 <- mergeGenoInfo$sigma2[i]
    p_bat <- mergeGenoInfo$pvalue_bat[i]
    w.ext <- mergeGenoInfo$w.ext[i]
    var.ratio.w0 <- mergeGenoInfo$var.ratio.w0[i]
    var.ratio.int <- mergeGenoInfo$var.ratio.int[i]
    var.ratio0 <- mergeGenoInfo$var.ratio.ext[i]

    WtCoxG.ext <- WtCoxG.test(
      g = g,
      R = R,
      w = w,
      TPR = TPR,
      sigma2 = sigma2,
      b = w.ext,
      var.ratio.w0 = var.ratio.w0,
      var.ratio.w1 = var.ratio.w0,
      var.ratio0 = var.ratio0,
      var.ratio1 = var.ratio0,
      mu.ext = mu.ext,
      n.ext = n.ext,
      p_bat = p_bat,
      p_cut = cutoff
    )

    WtCoxG.noext <- WtCoxG.test(
      g = g,
      R = R,
      w = w,
      var.ratio.int = var.ratio.int,
      p_bat = p_bat,
      p_cut = cutoff
    )

    return(cbind(WtCoxG.ext, WtCoxG.noext))
  }) %>%
    do.call("rbind", .) %>%
    cbind(., mergeGenoInfo)

  return(GWAS)
}


# SPA method----------------------------------------------------------------------
# calculate p-value --------------------------------------------------
WtCoxG.test <- function(
    g,
    R,
    w,
    p_bat,
    TPR = NA,
    sigma2 = NA,
    b = 0,
    var.ratio.int = 1, ### variance ratio of S.int
    var.ratio.w0 = 1, ### variance ratio of S.bat when bathceffect = 0
    var.ratio.w1 = 1, ### variance ratio of S.bat when bathceffect = 1
    var.ratio0 = 1, ### variance ratio of Score when bathceffect = 0
    var.ratio1 = 1, ### variance ratio of Score when bathceffect = 1
    mu.ext = NA,
    n.ext = NA,
    p_cut = 0.1 ### batcheffect p value cut off
    ) {
  ## imputation missing SNP
  missing.rate <- mean(is.na(g))
  pos.na <- which(is.na(g))
  if (missing.rate != 0) {
    g[pos.na] <- mean(na.omit(g))
  }


  ### if external MAF is unavailable, then S.int = sum(R*(g-mu.int))
  if (is.na(mu.ext)) {
    mu.int <- mean(g) / 2
    p.con <- SPA_G.one.SNP_homo(g = g, R = R, mu.ext = NA, n.ext = 0, sigma2 = 0, var.ratio = var.ratio.int)[1]
    p_deno <- NA
    return(p.con)
  }

  if (p_bat < p_cut | is.na(p_bat) | sum(g) < 10 | sum(2 - g) < 10) {
    p.con <- NA
    return(p.con)
  } else {
    meanR <- mean(R)
    sumR <- sum(R)
    mu.int <- mean(g) / 2
    N <- length(g)


    mu <- (1 - b) * mu.int + b * mu.ext
    S <- sum(R * (g - 2 * mu))
    # S=sum((R-(1-b)*meanR)*g)-sumR*2*b*mu.ext

    w1 <- w / (2 * sum(w))

    var_mu_ext <- mu * (1 - mu) / (2 * n.ext)
    var_Sbat <- sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext


    lb <- -qnorm(1 - p_cut / 2) * sqrt(var_Sbat) * sqrt(var.ratio.w0)
    ub <- qnorm(1 - p_cut / 2) * sqrt(var_Sbat) * sqrt(var.ratio.w0)
    c <- pnorm(ub / sqrt(var.ratio.w1), 0, sqrt(var_Sbat + sigma2), log.p = T)
    d <- pnorm(lb / sqrt(var.ratio.w1), 0, sqrt(var_Sbat + sigma2), log.p = T)
    p_deno <- TPR * (exp(d) * (exp(c - d) - 1)) + (1 - TPR) * (1 - p_cut)

    ## sigma2=0

    p_spa_s0 <- SPA_G.one.SNP_homo(
      g = g, R = R, b = b, mu.ext = mu.ext, n.ext = n.ext, sigma2 = 0,
      var.ratio = var.ratio0
    )[1]
    var_S <- S^2 / var.ratio0 / qchisq(p_spa_s0, 1, ncp = 0, lower.tail = F)


    var.int <- sum((R - (1 - b) * meanR)^2) * 2 * mu * (1 - mu)
    # var_S = var.int +  4*b^2 *sumR^2* var_mu_ext
    cov_Sbat_S <- sum(w1 * (R - (1 - b) * meanR)) * 2 * mu * (1 - mu) + 2 * b * sumR * var_mu_ext
    cov_Sbat_S <- cov_Sbat_S * sqrt(var_S / (var.int + 4 * b^2 * sumR^2 * var_mu_ext))
    VAR <- matrix(c(var_S, cov_Sbat_S, cov_Sbat_S, var_Sbat), nrow = 2)
    p0 <- max(0, min(1, mvtnorm::pmvnorm(
      lower = c(-Inf, lb / sqrt(var.ratio.w0)),
      upper = c(-abs(S / sqrt(var.ratio0)), ub / sqrt(var.ratio.w0)), mean = c(0, 0), sigma = VAR
    )))

    ## sigma2!=0
    p_spa_s1 <- SPA_G.one.SNP_homo(
      g = g, R = R, b = b, mu.ext = mu.ext, n.ext = n.ext, sigma2 = sigma2,
      var.ratio = var.ratio1
    )[1]

    var_S1 <- S^2 / var.ratio1 / qchisq(p_spa_s1, 1, ncp = 0, lower.tail = F)

    # var_S1 = var.int +  4*b^2 *sumR^2* (var_mu_ext+sigma2)
    cov_Sbat_S1 <- sum(w1 * (R - (1 - b) * meanR)) * 2 * mu * (1 - mu) + 2 * b * sumR * (var_mu_ext + sigma2)
    cov_Sbat_S1 <- cov_Sbat_S1 * sqrt(var_S1 / (var.int + 4 * b^2 * sumR^2 * (var_mu_ext + sigma2)))
    var_Sbat1 <- var_Sbat + sigma2
    VAR1 <- matrix(c(var_S1, cov_Sbat_S1, cov_Sbat_S1, var_Sbat1), nrow = 2)
    p1 <- max(0, min(1, mvtnorm::pmvnorm(
      lower = c(-Inf, lb / sqrt(var.ratio.w1)),
      upper = c(-abs(S / sqrt(var.ratio1)), ub / sqrt(var.ratio.w1)), mean = c(0, 0), sigma = VAR1
    )))

    p.con <- 2 * (TPR * p1 + (1 - TPR) * p0) / p_deno

    return(p.con)
  }
}

#############################################################################
# hidden function----------------------------------------------------------------------
# MGF function --------------------------------------------------
M_G0 <- function(t, MAF) {
  re <- (1 - MAF + MAF * exp(t))^2
  return(re)
}

# The first derivative of the MGF of G (genotype)
M_G1 <- function(t, MAF) {
  re <- 2 * (MAF * exp(t)) * (1 - MAF + MAF * exp(t))
  return(re)
}

# The second derivative of the MGF of G (genotype)
M_G2 <- function(t, MAF) {
  re <- 2 * (MAF * exp(t))^2 + 2 * (MAF * exp(t)) * (1 - MAF + MAF * exp(t))
  return(re)
}

# The CGF of G (genotype)
K_G0 <- function(t, MAF) {
  re <- log(M_G0(t, MAF))
  return(re)
}

K_G1 <- function(t, MAF) {
  re <- M_G1(t, MAF) / M_G0(t, MAF)
  return(re)
}

# The second derivative of the CGF of G (genotype)
K_G2 <- function(t, MAF) {
  re <- M_G0(t, MAF) / M_G0(t, MAF) * M_G2(t, MAF) / M_G0(t, MAF) - (M_G1(t, MAF) / M_G0(t, MAF))^2
  return(re)
}

# The CGF of score test statistic
H_org <- function(t, R, MAF, n.ext, N.all, sumR, var_mu_ext, g.var.est, meanR, b) {
  n.t <- length(t)
  out <- rep(0, n.t)
  R_hat <- sumR / N.all

  mu.adj <- -2 * b * sumR * MAF
  var.adj <- 4 * b^2 * sumR^2 * var_mu_ext
  for (i in 1:n.t) {
    t1 <- t[i]
    out[i] <- sum(K_G0(t1 * (R - (1 - b) * meanR), MAF)) + mu.adj * t1 + var.adj / 2 * t1^2
  }

  return(out)
}

# The first derivative of the CGF of score test statistic
H1_adj <- function(t, R, s, MAF, n.ext, N.all, sumR, var_mu_ext, g.var.est, meanR, b) {
  n.t <- length(t)
  out <- rep(0, n.t)
  R_hat <- sumR / N.all

  mu.adj <- -2 * b * sumR * MAF
  var.adj <- 4 * b^2 * sumR^2 * var_mu_ext

  for (i in 1:n.t) {
    t1 <- t[i]
    out[i] <- sum((R - (1 - b) * meanR) * K_G1(t1 * ((R - (1 - b) * meanR)), MAF)) + mu.adj + var.adj * t1 - s
  }
  return(out)
}

# The second derivative of the CGF of score test statistic
H2 <- function(t, R, MAF, n.ext, N.all, sumR, var_mu_ext, g.var.est, meanR, b) {
  n.t <- length(t)
  out <- rep(0, n.t)
  R_hat <- sumR / N.all

  mu.adj <- -n.ext * R_hat * 2 * MAF
  var.adj <- n.ext * R_hat^2 * 2 * MAF * (1 - MAF)


  for (i in 1:n.t) {
    t1 <- t[i]
    out[i] <- sum((R - (1 - b) * meanR)^2 * K_G2(t1 * (R - (1 - b) * meanR), MAF)) + var.adj
  }
  return(out)
}

GetProb_SPA_G <- function(MAF, R, s, n.ext, N.all, sumR, var_mu_ext, g.var.est, meanR, b, lower.tail) {
  out <- uniroot(H1_adj, c(-1, 1),
    extendInt = "yes",
    R = R, s = s, MAF = MAF, n.ext = n.ext, N.all = N.all, sumR = sumR,
    var_mu_ext = var_mu_ext, g.var.est = g.var.est, meanR = meanR, b = b
  )
  zeta <- out$root

  k1 <- H_org(zeta,
    R = R, MAF = MAF, n.ext = n.ext, N.all = N.all, sumR = sumR,
    var_mu_ext = var_mu_ext, g.var.est = g.var.est, meanR = meanR, b = b
  )
  k2 <- H2(zeta,
    R = R, MAF = MAF, n.ext = n.ext, N.all = N.all, sumR = sumR,
    var_mu_ext = var_mu_ext, g.var.est = g.var.est, meanR = meanR, b = b
  )

  temp1 <- zeta * s - k1

  w <- sign(zeta) * (2 * temp1)^{
    1 / 2
  }
  v <- zeta * (k2)^{
    1 / 2
  }

  pval <- pnorm(w + 1 / w * log(v / w), lower.tail = lower.tail)
  # pval.norm = pnorm(q2, lower.tail = lower.tail)
  re <- pval
  return(re)
}


# saddlepoint approximation (SPA) to calicrate the p value in the case that batcheffect = 0 or 1
SPA_G.one.SNP_homo <- function(
    g, ### genotype vector
    R, ### null model residual vector
    mu.ext = NA, ### external MAF
    n.ext = NA, ### external sample size
    b = 0, ### weight of external MAF
    sigma2 = NA, ###
    var.ratio = 1,
    Cutoff = 2,
    impute.method = "fixed",
    missing.cutoff = 0.15,
    min.mac = 10, # update on 2022-08-16 : replace 0.0001 by 0.000001
    G.model = "Add") {
  ## calculate MAF and update genotype vector

  ## imputation missing SNP
  missing.rate <- mean(is.na(g))
  pos.na <- which(is.na(g))
  if (missing.rate != 0) {
    g[pos.na] <- mean(na.omit(g))
  }


  if (is.na(mu.ext)) {
    mu.ext <- 0
    n.ext <- 0
  }


  if (sum(g) < min.mac | sum(2 - g) < min.mac | missing.rate > missing.cutoff) {
    MAF <- mean(na.omit(g)) / 2

    return(c(NA, NA))
  }


  ######################################################################

  ## Score statistic
  N <- length(g)
  mu.int <- mean(g) / 2
  MAF <- (1 - b) * mu.int + b * mu.ext
  sumR <- sum(R)
  N.all <- N + n.ext
  S <- sum(R * (g - 2 * MAF))
  S <- S / var.ratio

  ## estimated variance without adjusting for covariates
  g.var.est <- 2 * MAF * (1 - MAF)
  var_mu_ext <- ifelse(n.ext == 0, 0, MAF * (1 - MAF) / (2 * n.ext) + sigma2)


  #  S.var = sum(R^2 * g.var.est + sum(R)^2 * g.var.est/n.ext)
  meanR <- mean(R)
  S.var <- sum((R - (1 - b) * meanR)^2) * g.var.est + 4 * b^2 * sumR^2 * var_mu_ext

  z <- S / sqrt(S.var)

  if (abs(z) < Cutoff) {
    pval.norm <- pnorm(abs(z), lower.tail = FALSE) * 2
    return(c(pval.norm, pval.norm)) # update on 2022-10-05 : MAF.est.negative.num
  } else {
    pval1 <- GetProb_SPA_G(MAF,
      R = R, abs(S), n.ext = n.ext, N.all = N.all,
      var_mu_ext = var_mu_ext, g.var.est = g.var.est, meanR = meanR,
      sumR = sumR, b = b, lower.tail = FALSE
    ) # EmpSPA-G p value
    pval2 <- GetProb_SPA_G(MAF,
      R = R, -abs(S), n.ext = n.ext, N.all = N.all,
      var_mu_ext = var_mu_ext, g.var.est = g.var.est, meanR = meanR,
      sumR = sumR, b = b, lower.tail = TRUE
    ) # EmpSPA-G p value

    pval3 <- pnorm(abs(z), lower.tail = FALSE) # Normal
    pval4 <- pnorm(-abs(z), lower.tail = TRUE) # Normal

    pval.spa.G <- pval1 + pval2
    pval.norm <- pval3 + pval4

    # if(abs(z) < Cutoff){
    #   pval.spa.G = pval.norm
    # }

    return(c(pval.spa.G, pval.norm))
  }
}


##############################################################################
# GRAB.marker(obj.WtCoxG) is calculated by R now and will be rewriten by Rcpp.
# The functions below are from WtSPAG for Rcpp and is needed for adaption now.

checkControl.Marker.WtCoxG <- function(control) {
  return(control)
  
  default.control <- list(SPA_Cutoff = 2)
  control <- updateControl(control, default.control)
  return(control)
}


setMarker.WtCoxG = function(objNull, control) {
  return()

  QCFile = control$QCFile
  QCCutoff = control$QCCutoff
  
  QCData = data.table::fread(QCFile)
  QCData = QCData %>% select(genoIndex, AF_ref, AN_ref, pvalue_bat)
  
  mresid = objNull$mresid
  N = objNull$N
  weight = objNull$weight
  
  SPA_Cutoff = control$SPA_Cutoff
  
  obj.setMarker = list(QCData = QCData,
                       QCCutoff = QCCutoff)
  
  # The following function is in Main.cpp
  setWtSPAGobjInCPP(mresid,
                    # weight,
                    N,
                    SPA_Cutoff,
                    objNull$outLierList)
  
  return(obj.setMarker)
}


mainMarker.WtCoxG = function(genoType, genoIndex, outputColumns, obj.setMarker) {
  return()

  QCData = obj.setMarker$QCData
  QCCutoff = obj.setMarker$QCCutoff
  
  genoIndexSet = genoIndex
  QCData = QCData %>% filter(genoIndex %in% genoIndexSet)
  
  # QCData %>% head() %>% print()
  genoIndex = QCData$genoIndex
  AF_ref = QCData$AF_ref
  AN_ref = QCData$AN_ref
  pvalue_bat = QCData$pvalue_bat
  
  # The below function is in Main.cpp
  updateQCInCPP(AF_ref, AN_ref, pvalue_bat, QCCutoff)
  
  # The below function is in Main.cpp
  OutList = mainMarkerInCPP("WtSPAG", genoType, genoIndex)
  obj.mainMarker = data.frame(Marker = OutList$markerVec,           # marker IDs
                              Info = OutList$infoVec,               # marker information: CHR:POS:REF:ALT
                              AltFreq = OutList$altFreqVec,         # alternative allele frequencies
                              AltCounts = OutList$altCountsVec,     # alternative allele counts
                              MissingRate = OutList$missingRateVec, # alternative allele counts
                              Pvalue = OutList$pvalVec)             # marker-level p-values
  
  # optionalColumns = c("beta", "seBeta", "zScore", "PvalueNorm", "AltFreqInGroup", "AltCountsInGroup", "nSamplesInGroup")
  optionalColumns = c("zScore", "PvalueNorm", "AltFreqInGroup", "AltCountsInGroup", "nSamplesInGroup")
  additionalColumns = intersect(optionalColumns, outputColumns)
  
  if(length(additionalColumns) > 0)
    obj.mainMarker = cbind.data.frame(obj.mainMarker, 
                                      as.data.frame(OutList[additionalColumns]))
  
  return(obj.mainMarker)
}