# This file allows you to source in the many functions we have implemented 
# in lab.

# Sourcing is functionally very similar to installing a package with a few less
# perks.

# Required Packages ####
library(psych)
library(car)
library(tidyverse)
library(lme4)
###############################################################################
# ggplotPredict.####

ggplotPredict = function (Model, Data = NULL, Label = NULL, Type = "response") 
{
  if (is.null(Data) & class(Model)[1] == "lm") {
    return(fitted.values(Model))
  }
  else {
    if (is.null(Label)) {
      PredictName = "Predicted"
      CILoName = "CILo"
      CIHiName = "CIHi"
      SEName = "SE"
    }
    else {
      PredictName = paste0("Predicted", Label)
      CILoName = paste0("CILo", Label)
      CIHiName = paste0("CIHi", Label)
      SEName = paste0("SE", Label)
    }
    Predictions = matrix(data = NA, nrow = nrow(Data), ncol = 4, 
                         dimnames = list(1:nrow(Data), c(PredictName, CILoName, 
                                                         CIHiName, SEName)))
    if (class(Model)[1] == "lm") {
      CILevel = 1 - 2 * pt(c(1), df = Model$df.residual, 
                           lower.tail = FALSE)
      Predictions[, 1:3] = predict(Model, newdata = Data, 
                                   interval = "confidence", level = CILevel)
      Predictions[, 4] = Predictions[, 1] - Predictions[, 
                                                        2]
      Predictions = as.data.frame(Predictions)
    }
    if (class(Model)[1] == "glm") {
      tmpPred = predict(Model, newdata = Data, type = "link", 
                        se.fit = TRUE)
      upr <- tmpPred$fit + tmpPred$se.fit
      lwr <- tmpPred$fit - tmpPred$se.fit
      fit <- tmpPred$fit
      if (Type == "response") {
        fit <- Model$family$linkinv(fit)
        upr <- Model$family$linkinv(upr)
        lwr <- Model$family$linkinv(lwr)
      }
      Predictions[, 1] = fit
      Predictions[, 2] = lwr
      Predictions[, 3] = upr
      Predictions[, 4] = Predictions[, 1] - Predictions[, 
                                                        2]
      Predictions = as.data.frame(Predictions)
    }
    if ((class(Model)[1] == "lmerMod") || (class(Model)[1] == 
                                           "glmerMod")) {
      Predictions[, c(1, 4)] = predictSE(Model, Data, 
                                         se.fit = TRUE, type = Type, level = 0, print.matrix = TRUE)
      Predictions[, 2] = Predictions[, 1] - Predictions[, 
                                                        4]
      Predictions[, 3] = Predictions[, 1] + Predictions[, 
                                                        4]
    }
    if (any(names(Data) == PredictName) || any(names(Data) == 
                                               CILoName) || any(names(Data) == CIHiName) || any(names(Data) == 
                                                                                                SEName)) {
      warning("Variable names (Predicted, CILo, CIHi, SE with Label PostFix) used in Data.  These variables removed before merging in predicted values")
      Data[, c(PredictName, CILoName, CIHiName, SEName)] = list(NULL)
    }
    Data = data.frame(Predictions, Data)
    return(Data)
  }
}
###############################################################################
# corr_table ####
corr_table <- function(vars, use_method) {
  # Run the corr.test function
  corr_results <- corr.test(d[, vars], use = use_method)
  # Make a dataframe called d_corrs that's empty other than row.names that are our vars
  d_corrs <- data.frame(row.names = vars)
  # Our first 'for' loop!! This loops through the code between the two {} once 
  # for every item in 'vars,' calling each one 'v' in turn. 
  # Once it runs everything between the {} for that item, then it loops back up 
  # to do it all again for the next item, until it runs out of items. 
  for (v in vars) {
    print(v)
    for (pair_v in vars) { # a nested for loop
      # Obtain the r value from the correlation table where row = v, column = pair_v
      # and round it to 3 digits
      r <- round(corr_results$r[v, pair_v], digits = 3) 
      # Obtain the p value from the correlation table where row = v, column = pair_v
      # and round it to 3 digits
      p <- round(corr_results$p[v, pair_v], digits = 3)
      # paste() puts it all together
      v_r_p <- paste("r = ", r, ", (p = ", p, ")", sep = "")
      # Now, in our d_corrs dataframe where row = v, column = pair_v, put our 
      # pasted together result
      d_corrs[v, pair_v] <- v_r_p
    }
  }
  return(d_corrs)
}
###############################################################################
# center ####
center <- function(var) {
  (var - mean(var, na.rm = T))
} # mean center a variable 
###############################################################################
# modelCaseAnalysis ####
modelCaseAnalysis = function (Model, Type = "RESIDUALS", Term = NULL, ID = row.names(Model$model)) 
{
  switch(toupper(Type), UNIVARIATE = {
    d = Model$model
    {
      Vars = names(d)
    }
    for (varname in Vars) {
      par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2, ask = TRUE)
      if (is.factor(d[[varname]])) {
        plot(d[varname], xlab = varname, ylab = "Frequency")
      } else {
        hist(d[[varname]], xlab = varname, main = "Red: Mean +- 3SD; Green: Median +- 2.2IQR")
        points(d[[varname]], rep(0, length(d[[varname]])), 
               pch = "|", col = "blue")
        abline(v = c(-3, 0, 3) * sd(d[[varname]]) + 
                 mean(d[[varname]]), col = "red", lty = c(1, 
                                                          2, 1))
        abline(v = c(-2.2, 0, 2.2) * IQR(d[[varname]]) + 
                 median(d[[varname]]), col = "green", lty = c(1, 
                                                              2, 1))
        if (!is.null(ID)) {
          Indices = identify(d[[varname]], rep(0, length(d[[varname]])), 
                             labels = ID)
          Values = d[[varname]][Indices]
        }
      }
    }
  }, HATVALUES = {
    par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2)
    TheTitle = paste("Model: ", Model$call[2], "\n", "Small sample cut (green) = 3 * mean(Hat)\nLarge sample cut: 2 * mean(Hat)", 
                     sep = "")
    hist(hatvalues(Model), xlab = "Hat Values", main = TheTitle)
    abline(v = c(2, 3) * mean(hatvalues(Model)), col = c("red", 
                                                         "green"))
    points(hatvalues(Model), rep(0, length(hatvalues(Model))), 
           pch = "|", col = "blue")
    if (!is.null(ID)) {
      Indices = identify(hatvalues(Model), rep(0, length(hatvalues(Model))), 
                         labels = ID)
      Values = hatvalues(Model)[Indices]
    }
  }, RESIDUALS = {
    par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2)
    N = length(rstudent(Model))
    k = length(coef(Model)) - 1
    TCut <- qt(p = 0.025/N, df = N - k - 2, lower.tail = FALSE)
    TheTitle = paste("Model: ", Model$call[2], "\n", "Bonferroni corrected p < .05 cut-off in red", 
                     sep = "")
    hist(rstudent(Model), xlab = "Studentized Residuals", 
         main = TheTitle)
    abline(v = c(-1, 1) * TCut, col = "red")
    points(rstudent(Model), rep(0, length(rstudent(Model))), 
           pch = "|", col = "blue")
    if (!is.null(ID)) {
      Indices = identify(rstudent(Model), rep(0, length(rstudent(Model))), 
                         labels = ID)
      Values = rstudent(Model)[Indices]
    }
  }, COOKSD = {
    par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2)
    N = length(cooks.distance(Model))
    k = length(coef(Model)) - 1
    TheTitle = paste("Model: ", Model$call[2], "\n", "4/(N-P) cut-off (red)\nqf(.5,P,N-P) cut-off (green)", 
                     sep = "")
    hist(cooks.distance(Model), xlab = "Cooks d", main = TheTitle)
    abline(v = c((4/(N - k - 1)), qf(0.5, k + 1, N - k - 
                                       1)), col = c("red", "green"))
    points(cooks.distance(Model), rep(0, length(cooks.distance(Model))), 
           pch = "|", col = "blue")
    if (!is.null(ID)) {
      Indices = identify(cooks.distance(Model), rep(0, 
                                                    length(cooks.distance(Model))), labels = ID)
      Values = cooks.distance(Model)[Indices]
    }
  }, DFBETAS = {
    if (is.null(Term)) {
      {
        Vars = dimnames(dfbetas(Model))[[2]]
      }
    } else {
      if (!(Term %in% dimnames(dfbetas(Model))[[2]])) {
        stop("Term specified for DFBETAS not valid")
      } else Vars = Term
    }
    for (varname in Vars) {
      par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2, ask = TRUE)
      TheTitle = paste("Model: ", Model$call[2], "\n", 
                       "B= ", coef(Model)[varname], sep = "")
      hist(dfbetas(Model)[, varname], xlab = paste("DFBETAS:", 
                                                   varname), main = TheTitle)
      points(dfbetas(Model)[, varname], rep(0, length(dfbetas(Model)[, 
                                                                     varname])), pch = "|", col = "blue")
      abline(v = c(-2, 2), col = "red")
      if (!is.null(ID)) {
        Indices = identify(dfbetas(Model)[, varname], 
                           rep(0, length(dfbetas(Model)[, varname])), 
                           labels = ID)
      }
    }
    par(ask = FALSE)
    if (is.null(Term)) {
      avPlots(Model, intercept = TRUE, id.method = "identify", 
              id.n = nrow(dfbetas(Model)), labels = ID)
    } else {
      avPlots(Model, terms = Term, id.method = "identify", 
              id.n = nrow(dfbetas(Model)), labels = ID)
      Values = dfbetas(Model)[Indices, Term]
    }
    if (is.null(Term)) {
      Indices = NULL
      Values = NULL
    }
  }, INFLUENCEPLOT = {
    par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2)
    TheTitle = paste("Influence Bubble plot", "\nModel: ", 
                     Model$call[2], sep = "")
    plot(hatvalues(Model), rstudent(Model), type = "n", 
         xlab = "Hat Values", ylab = "Studentized Residuals", 
         main = TheTitle)
    cooksize = 10 * sqrt(cooks.distance(Model))/max(cooks.distance(Model))
    points(hatvalues(Model), rstudent(Model), cex = cooksize)
    N = length(rstudent(Model))
    k = length(coef(Model)) - 1
    TCut <- qt(p = 0.025/N, df = N - k - 2, lower.tail = FALSE)
    abline(h = c(-1, 0, 1) * TCut, col = "red", lty = c(1, 
                                                        2, 1))
    abline(v = c(1, 2, 3) * mean(hatvalues(Model)), col = "red", 
           lty = c(2, 1, 1))
    if (!is.null(ID)) {
      Indices = identify(hatvalues(Model), rstudent(Model), 
                         labels = ID)
      Values = c(hatvalues(Model)[Indices], rstudent(Model)[Indices])
    }
  }, COVRATIO = {
    N = length(covratio(Model))
    k = length(coef(Model)) - 1
    par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2)
    TheTitle = paste("Model: ", Model$call[2], "\n", "abs((3*P)/N)-1 cut-off in red", 
                     sep = "")
    hist(covratio(Model), xlab = "CovRatio", main = TheTitle)
    abline(v = abs((3 * (k + 1)/N) - 1), col = "red")
    points(covratio(Model), rep(0, length(covratio(Model))), 
           pch = "|", col = "blue")
    if (!is.null(ID)) {
      Indices = identify(covratio(Model), rep(0, length(covratio(Model))), 
                         labels = ID)
      Values = covratio(Model)[Indices]
    }
  }, {
    print("Valid options for type: hatvalues, residuals, cooksd, dfbetas, influenceplot, covratio, univariate")
  })
  Rownames = row.names(Model$model)[Indices]
  Cases = list(Rownames = Rownames, Values = Values)
  return(Cases)
} # Outlier analyses
###############################################################################
# standardize ####
standardize <- function(var) {
  ((var - mean(var, na.rm = T))/(sd(var, na.rm = T)))
} # create z scores
###############################################################################
# figStripChart ####
figStripChart = function (x, side = 1, sshift = 0.3, adjoffset = 1, strip.col = "gray", 
                          strip.pch = 15, strip.cex = 0.2) 
{
  if (is.null(strip.col)) 
    strip.col = getOption("FigPars")$strip.col
  if (is.null(strip.pch)) 
    strip.pch = getOption("FigPars")$strip.pch
  if (is.null(strip.cex)) 
    strip.cex = getOption("FigPars")$strip.cex
  if (side == 1) {
    flip <- 1
    yaxis <- FALSE
    parside <- 3
  }
  else if (side == 2) {
    flip <- 1
    yaxis <- TRUE
    parside <- 1
  }
  else if (side == 3) {
    flip <- -1
    yaxis <- FALSE
    parside <- 4
  }
  else if (side == 4) {
    flip <- -1
    yaxis <- TRUE
    parside <- 2
  }
  base <- par("usr")[parside]
  plotwidth <- diff(par("usr")[1:2])
  plotheight <- diff(par("usr")[3:4])
  shift <- par("pin")[1] * 0.003 * flip
  gap <- par("pin")[1] * 0.003
  meanshift <- par("cin")[1] * 0.5 * flip
  stripshift <- par("cin")[1] * sshift * flip
  if (yaxis) {
    shift <- shift/par("pin")[1] * plotwidth
    meanshift <- meanshift/par("pin")[1] * plotwidth
    stripshift <- stripshift/par("pin")[1] * plotwidth
    gap <- gap/par("pin")[2] * plotheight
  }
  else {
    shift <- shift/par("pin")[2] * plotheight
    meanshift <- meanshift/par("pin")[2] * plotheight
    stripshift <- stripshift/par("pin")[2] * plotheight
    gap <- gap/par("pin")[1] * plotwidth
  }
  if (yaxis) 
    offset = flip * par("cin")[2]/par("cin")[1] * adjoffset
  else offset = flip * adjoffset
  oldxpd <- par(xpd = TRUE)
  on.exit(par(oldxpd))
  stripchart(x, method = "stack", vertical = yaxis, offset = offset, 
             pch = strip.pch, cex = strip.cex, add = TRUE, at = base + 
               shift + stripshift, col = strip.col)
} # This function is needed for varPlot to function
################################################################################
# varPlot ####
varPlot = function (TheVar, VarName = "", IDs = NULL, AddPoints = "Strip", 
          AddDensity = TRUE, Detail = 2) 
{
  print(paste("Descriptive statistics: ", VarName, sep = ""))
  print(describe(TheVar, Detail))
  par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2)
  HistData = hist(TheVar, main = "", xlab = VarName)
  switch(toupper(AddPoints), STRIP = {
    figStripChart(TheVar, strip.col = "red", , strip.cex = 0.3)
  }, RUG = {
    rug(TheVar, col = "red")
  })
  if (AddDensity) {
    DensityData = density(TheVar)
    lines(DensityData$x, DensityData$y * (max(HistData$counts)/max(DensityData$y)), 
          col = "blue")
  }
  if (!is.null(IDs)) {
    Indices = identify(TheVar, rep(0, length(TheVar)), labels = IDs)
    return(Cases = list(Indices = Indices, Rownames = IDs[Indices], 
                        Values = TheVar[Indices]))
  }
}
################################################################################
# modelAssumptions ####
modelAssumptions = function (Model, Type = "NORMAL", ID = row.names(Model$model), 
          one.page = TRUE) 
{
  switch(toupper(Type), NORMAL = {
    if (one.page) {
      dev.new(width = 14, height = 7)
      par(mfrow = c(1, 2))
    } else {
      dev.new(width = 7, height = 7, record = TRUE)
    }
    par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2)
    qqPlot(Model, labels = FALSE, sim = TRUE, main = "Quantile-Comparison Plot to Assess Normality", 
           xlab = "t Quantiles", ylab = "Studentized Residuals")
    plot(density(rstudent(Model)), main = "Density Plot to Assess Normality of Residuals", 
         xlab = "Studentized Residual")
    zx <- seq(-4, 4, length.out = 100)
    lines(zx, dnorm(zx, mean = 0, sd = sd(rstudent(Model))), 
          lty = 2, col = "blue")
    cat("Descriptive Statistics for Studentized Residuals\n")
    describe(rstudent(Model))
  }, CONSTANT = {
    if (one.page) {
      dev.new(width = 14, height = 7)
      par(mfrow = c(1, 2))
    } else {
      dev.new(width = 7, height = 7, record = TRUE)
    }
    par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2)
    plot(rstudent(Model) ~ fitted.values(Model), main = "Studentized Residuals vs. Fitted Values", 
         xlab = "Fitted Values", ylab = "Studentized Residuals")
    abline(h = 0, lty = 2, col = "blue")
    print(spreadLevelPlot(Model))
    cat("\n\n")
    print(ncvTest(Model))
  }, LINEAR = {
    dev.new(width = 7, height = 7, record = TRUE)
    par(cex.lab = 1.5, cex.axis = 1.2, lwd = 2)
    crPlots(Model, ask = TRUE)
  }, {
    print("Valid options for type: normal, constant, linear")
  })
  print(gvlma(Model))
}


################################################################################
# modelBoxCox ####
modelBoxCox = function (Model, Lambdas = seq(-2, 2, by = 0.1)) 
{
  LR <- boxCox(Model, lambda = Lambdas)
  Lambda1Index <- sum(LR$x < 1)
  Chi1 <- mean(c(LR$y[Lambda1Index], LR$y[Lambda1Index + 1]))
  ChiLambda <- LR$y[which.max(LR$y)]
  ChiDiff <- 2 * (ChiLambda - Chi1)
  print(paste("Best Lambda=", round(LR$x[which.max(LR$y)], 
                                    2)))
  print(paste("Chi-square (df=1)=", round(ChiDiff, 2)))
  print(paste("p-value=", round(pchisq(ChiDiff, df = 1, lower.tail = FALSE), 
                                5)))
}

################################################################################
# quiet ####
# Allows you to run modelCaseAnalysis in r markdown without it forcing your to
# click on a point. You MUST still run modelCaseAnalysis in a separate r file if
# you wish to identify problematic points.
quiet = function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

################################################################################
# read_merge_many ####
read_merge_many <- function(files, merge_by) {
  for (f in files) {
    print(f)
    df_f <- read_csv(paste("", f, sep = ""))
    if (exists('df') && is.data.frame(get('df'))) {
      df <- merge(df, df_f, by = merge_by)
    }
    else {
      df <- df_f
    }
  }
  return(df)
}

################################################################################
# power_analysis ####
#' Determine the minimum detectable effect size
#'
#' Given desired power, degrees of freedom, and alpha, this function returns the
#' minimum detectable effect size as an f-squared statistic.
#'
#' @param desired_power Defaults to .8, but this can be changed to test for the
#' minimum detectable effect size for other levels of power.
#' @param num_df Numerator degrees of freedom for the test.
#' @param den_df Denominator degrees of freedom for the test.
#' @param alpha Alpha level for the test.
#' @return Returns an f-squared statistic representing the minimum detectable effect
#' size with the desired power, degrees of freedom, and alpha.
#' @export
mdes <- function(desired_power = .8, num_df, den_df, alpha = .05) {
  
  # overshoot, nearest tenth
  power <- 0
  fsq <- 0
  while (power < desired_power) {
    fsq <- fsq + .01
    power <- pwr::pwr.f2.test(u = num_df, v = den_df, f2 = fsq, sig.level = alpha)$power
  }
  
  # overshoot, nearest thousandth
  fsq <- fsq - .01
  power <- 0
  while (power < desired_power) {
    fsq <- fsq + .0001
    power <- pwr::pwr.f2.test(u = num_df, v = den_df, f2 = fsq, sig.level = alpha)$power
  }
  
  # overshoot, nearest ten thousandth
  fsq <- fsq - .0001
  power <- 0
  while (power < desired_power) {
    fsq <- fsq + .00001
    power <- pwr::pwr.f2.test(u = num_df, v = den_df, f2 = fsq, sig.level = alpha)$power
  }
  
  return(fsq)
  
}


#' Calculate power and for GLM tests
#'
#' Wrapper to calculate power for tests of parameter estimates or full model in
#' GLM based on Cohen's tables and using pwr.f2.test in pwr package. Inspired by
#' modelPower in John Curtin's lmSupport package. Allows the use
#' of partial eta squared or delta R2 rather than just f2 as effect size. If
#' you provide power, it returns N, if you provide N, it returns power. If you
#' provide N and power, it returns the minimum detectable effect size given the
#' specified N and power. You must alwasy specify an effect size as either
#' f2, partial eta2, or delta R2 with model R2. You must also specify the
#' number of parameters in the compact (pc) and  augmented (pa) for the model
#' comparison that will test the effect.
#'
#' @param pc Number of parameters in the compact model; i.e., intercept + all
#' parameters excluding the effect of interest; This is the numerator df of the
#' F test for the effect.
#' @param pa Number of parameters in the augmented model; i.e., the intercept
#' and all parameters including the effect of interest.
#' @param n Sample size.
#' @param alpha Alpha for statistical test.
#' @param power Power for statistical test.
#' @param f2 Effect size.
#' @param peta2 = Partial eta2 effect size.
#' @param dr2 Delta r2 effect; if provided must also specify r2.
#' @param r2 Model r2, only needed if using delta r2.
#' @return Returns a list with n, power, possibly minimum detectable effect size.
#' @examples
#' # return the minimum detectable effect size with 200 participants, power of
#' # .8 (the default), pa of 5, and pc of 4:
#' power_analysis(pa = 5, pc = 4, n = 200)
#'
#' # return the number of participants needed for 70% power given
#' # pa of 3, pc of 2, and peta2 of .01
#' power_analysis(pa = 3, pc = 2, peta2 = .01, power = .7)
#'
#' # return the power of a study with peta2 of .02 and 50 participants, with
#' # pa of 5 and pc of 4.
#' power_analysis(pa = 5, pc = 4, peta2 = .02, n = 50)
#' @export
power_analysis <- function(pc = NULL, pa = NULL, n = NULL, alpha = 0.05, power = NULL,
                           f2 = NULL, peta2 = NULL, dr2 = NULL, r2 = NULL) {
  if (is.null(pa) | is.null(pc)) {
    stop("Must provide pa and pc")
  }
  
  u <- pa - pc
  mdes_peta2 <- NULL
  mdes_f2 <- NULL
  
  nEffs <- 0
  if (!is.null(f2)) {
    nEffs <- nEffs + 1
  }
  if (!is.null(peta2)) {
    f2 <- peta2 / (1 - peta2)
    nEffs <- nEffs + 1
  }
  if (!is.null(dr2) & !is.null(r2)) {
    f2 <- dr2 / (1 - r2)
    nEffs <- nEffs + 1
  }
  if (nEffs > 1) {
    stop("Must not specify more than one of the following: f2, peta2, or both dr2 and r2")
  }
  
  if (!is.null(n)) {
    v <- n - pa
  }
  else {
    v <- NULL
  }
  
  if (!is.null(pa) & !is.null(pc) & !is.null(n)) {
    mdes_power <- power
    if(is.null(mdes_power)) mdes_power <- .8
    mdes_f2 <- mdes(desired_power = mdes_power, u, v, alpha)
    mdes_peta2 <- mdes_f2 / (1 + mdes_f2)
  }
  
  if (nEffs != 0) {
    results <- pwr::pwr.f2.test(
      u = u, v = v, f2 = f2, sig.level = alpha,
      power = power
    )
    
    output <- list(
      n = pa + results$v,
      peta2 = results$f2 / (1 + results$f2),
      f2 = results$f2,
      power = results$power,
      mdes_peta2 = mdes_peta2,
      mdes_f2 = mdes_f2
    )
  }
  
  if (nEffs == 0) {
    print("2")
    output <- list(
      n = n,
      peta2 = NULL,
      f2 = NULL,
      power = power,
      mdes_peta2 = mdes_peta2,
      mdes_f2 = mdes_f2
    )
  }
  
  return(output)
}

################################################################################
# indirect.mlm ####
indirect.mlm <- function( data, indices, y, x, mediator, group.id, covariates=NULL, 
                          uncentered.x=T, 
                          random.a=T, random.b=T, random.c=T, 
                          between.x=F, between.m=F ) {
  # Create resampled dataset
  resampled.data <- data[array(indices),]
  
  # Centre variables
  # Predictor
  if ( uncentered.x == T ) {
    x.by.id <- as.data.frame( tapply(resampled.data[,x], resampled.data[,group.id], mean, na.rm=T ) )
    names( x.by.id ) <- "x.by.id"
    x.by.id[group.id] <- row.names(x.by.id)
    resampled.data <- merge( resampled.data, x.by.id, by=group.id )
    resampled.data[paste("c.", x, sep="")] <- resampled.data[x] - resampled.data["x.by.id"]
  }
  
  # Mediator
  mediator.by.id <- as.data.frame( tapply(resampled.data[,mediator], resampled.data[,group.id], mean, na.rm=T ) )
  names(mediator.by.id ) <- "mediator.by.id"
  mediator.by.id[group.id] <- row.names( mediator.by.id )
  resampled.data <- merge( resampled.data, mediator.by.id, by=group.id )
  resampled.data[paste("c.", mediator, sep="")] <- resampled.data[mediator] - resampled.data["mediator.by.id"]
  
  # Create between-group variables
  if (between.x == T) {
    if ( uncentered.x == T ) {
      resampled.data[paste("c.average.", x, sep="")] <- resampled.data["x.by.id"] - mean( unlist(resampled.data["x.by.id"]), na.rm=T )
    }
    if (uncentered.x == F ) {
      between.x = F
      cat("Cannot calculate between-group effect for predictor without the uncentered, raw values of the predictor. Please provide the raw values of the predictor if you want to calculate the between-group effect for the predictor. Proceeding without calculating a between-group effect for slope a (i.e., as if 'between.x=F').")
    }
  }
  if (between.m == T) {
    resampled.data[paste("c.average.", mediator, sep="")] <- resampled.data["mediator.by.id"] - mean( unlist(resampled.data["mediator.by.id"]), na.rm=T )
  }
  
  # Create fixed effects formulas
  a.path.formula <- as.formula(paste(mediator, " ~ (1", ifelse(random.a, paste( " + ", x, sep=""), ""), "|", group.id, ") + ", ifelse(!is.null(covariates), paste(paste(covariates, collapse=" + "), " + ", sep=""), ""), x, ifelse( between.x==T, paste(" + c.average.", x, sep=""), ""), sep=""))
  b.path.formula <- as.formula(paste(y, " ~ (1", ifelse(random.c&&!random.b, paste( " + ", x, sep=""), ifelse(random.c&&random.b, paste(" + ", x, " + c.", mediator, sep=""), ifelse(!random.c&&random.b, paste(" + ", "c.", mediator, sep=""), ""))), "|", group.id, ") + ", ifelse(!is.null(covariates), paste(paste(covariates, collapse=" + "), " + ", sep=""), ""), x, ifelse( between.x==T, paste(" + c.average.", x, sep=""), ""), " + ", paste("c.", mediator, sep=""), ifelse( between.m==T, paste(" + c.average.", mediator, sep=""), ""), sep=""))
  c.path.formula <- as.formula(paste(y, " ~ (1", ifelse(random.c, paste( " + ", x, sep=""), ""), "|", group.id, ") + ", ifelse(!is.null(covariates), paste(paste(covariates, collapse=" + "), " + ", sep=""), ""), x, ifelse( between.x==T, paste(" + c.average.", x, sep=""), ""), sep=""))
  
  # Run models with resampled dataset 
  model.a <- do.call( "lmer", list( formula=a.path.formula, data=resampled.data, na.action="na.exclude", REML=F) )
  model.b <- do.call( "lmer", list( formula=b.path.formula, data=resampled.data, na.action="na.exclude", REML=F) )
  model.c <- do.call( "lmer", list( formula=c.path.formula, data=resampled.data, na.action="na.exclude", REML=F) )
  
  # Save Fixed Slopes
  fixed.slope.within.a <- coef(summary(model.a))[x,1]
  fixed.slope.within.b <- coef(summary(model.b))[paste("c.", mediator, sep=""),1]
  fixed.slope.between.a <- NA
  fixed.slope.between.b <- NA
  if (between.x == T ) { fixed.slope.between.a <- coef(summary(model.a))[paste("c.average.", x, sep=""),1] }
  if (between.m == T ) { fixed.slope.between.b <- coef(summary(model.b))[paste("c.average.", mediator, sep=""),1] }
  fixed.slope.c <- coef(summary(model.c))[x,1]
  fixed.slope.c.prime <- coef(summary(model.b))[x,1]
  
  # Save Random Slopes
  random.slope.a <- NA
  random.slope.b <- NA
  random.slope.c.prime <- NA
  if (random.a==T) { random.slope.a <- as.numeric(unlist(coef(model.a)[group.id][[1]][x])) }
  if (random.b==T) { random.slope.b <- as.numeric(unlist(coef(model.b)[group.id][[1]][paste("c.", mediator, sep="")])) }
  if (random.c==T) { random.slope.c.prime <- as.numeric(unlist(coef(model.b)[group.id][[1]][x])) }
  
  #"Population Covariance" of random slopes
  random.covariance <- NA
  if (random.a==T&random.b==T) {random.covariance <- cov( random.slope.a, random.slope.b )}
  
  #Within-Group Indirect Effects
  within.indirect.effect <- mean( random.slope.a * random.slope.b, na.rm=T )
  within.indirect.effect.biased <- fixed.slope.within.a*fixed.slope.within.b
  within.indirect.bias <- abs( within.indirect.effect.biased - within.indirect.effect )
  
  #Between-Group Indirect Effects
  between.indirect.effect <- NA
  between.indirect.effect <- ifelse( between.x==T&between.m==T, fixed.slope.between.a*fixed.slope.between.b,
                                     ifelse( between.x==F&between.m==T, mean( random.slope.a * fixed.slope.between.b, na.rm=T ),
                                             ifelse( between.x==T&between.m==F, mean( fixed.slope.between.a * random.slope.b, na.rm=T ),
                                                     NA) ) )
  between.indirect.effect.biased <- NA
  between.indirect.effect.biased <- ifelse( between.x==T&between.m==T, fixed.slope.between.a*fixed.slope.between.b,
                                            ifelse( between.x==F&between.m==T, mean( fixed.slope.within.a * fixed.slope.between.b, na.rm=T ),
                                                    ifelse( between.x==T&between.m==F, mean( fixed.slope.between.a * fixed.slope.within.b, na.rm=T ),
                                                            NA) ) )
  between.indirect.bias <- abs( between.indirect.effect - between.indirect.effect.biased )
  
  # Total Effects
  total.effect <- mean( random.slope.a * random.slope.b + random.slope.c.prime, na.rm=T )
  if (random.c==F) {total.effect <- mean( random.slope.a * random.slope.b + fixed.slope.c.prime, na.rm=T )}
  total.effect.bias <- abs( fixed.slope.c - total.effect )
  
  return( c( fixed.slope.c, fixed.slope.c.prime, fixed.slope.within.a, fixed.slope.between.a, fixed.slope.within.b, fixed.slope.between.b, 
             random.covariance, within.indirect.effect, within.indirect.effect.biased, within.indirect.bias, between.indirect.effect, 
             between.indirect.effect.biased, between.indirect.bias, total.effect, total.effect.bias) )
}



################################################################################
# indirect.mlm.summary ####
indirect.mlm.summary <- function( boot.object ) {
  with( boot.object, {
    if ((is.null(call$random.a)||(call$random.a=="T")) & (is.null(call$random.b)||call$random.b=="T")) {
      cat( "#### Population Covariance ####\n" )
      cat( paste("Covariance of Random Slopes a and b: ", round(t0[ 7], 3), " [", round(boot.ci( boot.object, index= 7, type="perc" )$perc[4], 3), ", " , round(boot.ci( boot.object, index= 7, type="perc" )$perc[5], 3), "]", "\n", sep="") )
      cat( "\n\n" )
      cat( "#### Indirect Effects ####\n" )
      cat( "# Within-subject Effects\n" )
      cat( paste("Unbiased Estimate of Within-subjects Indirect Effect: ", round(t0[ 8], 3), " [", round(boot.ci( boot.object, index= 8, type="perc" )$perc[4], 3), ", " , round(boot.ci( boot.object, index= 8, type="perc" )$perc[5], 3), "]", "\n", sep="") )
      cat( paste("Biased Estimate of Within-subjects Indirect Effect: ", round(t0[ 9], 3), " [", round(boot.ci( boot.object, index= 9, type="perc" )$perc[4], 3), ", " , round(boot.ci( boot.object, index= 9, type="perc" )$perc[5], 3), "]", "\n", sep="") )
      cat( paste("Bias in Within-subjects Indirect Effect: ", round(t0[ 10], 3), " [", round(boot.ci( boot.object, index= 10, type="perc" )$perc[4], 3), ", " , round(boot.ci( boot.object, index= 10, type="perc" )$perc[5], 3), "]", "\n", sep="") )
      cat( "\n" )
      if ( (!is.null(call$between.x)&&call$between.x=="T") || (!is.null(call$between.m)&&call$between.m=="T") ) {
        cat( "# Between-subject Effects\n" )
        cat( paste("Unbiased Estimate of Between-subjects Indirect Effect: ", round(t0[ 11], 3), " [", round(boot.ci( boot.object, index= 11, type="perc" )$perc[4], 3), ", " , round(boot.ci( boot.object, index= 11, type="perc" )$perc[5], 3), "]", "\n", sep="") )
        cat( paste("Biased Estimate of Between-subjects Indirect Effect: ", round(t0[ 12], 3), " [", round(boot.ci( boot.object, index= 12, type="perc" )$perc[4], 3), ", " , round(boot.ci( boot.object, index= 12, type="perc" )$perc[5], 3), "]", "\n", sep="") )
        cat( paste("Bias in Between-subjects Indirect Effect: ", round(t0[ 13], 3), " [", round(boot.ci( boot.object, index= 13, type="perc" )$perc[4], 3), ", " , round(boot.ci( boot.object, index= 13, type="perc" )$perc[5], 3), "]", "\n", sep="") )
        cat( "\n" )
      }
    }
    else {
      cat( "#### Indirect Effects ####\n" )
      cat( paste("Within-subjects Indirect Effect: ", round(t0[ 9], 3), " [", round(boot.ci( boot.object, index= 9, type="perc" )$perc[4], 3), ", " , round(boot.ci( boot.object, index= 9, type="perc" )$perc[5], 3), "]", "\n", sep="") )
      cat( "\n" )
      if ( (!is.null(call$between.x)&&call$between.x=="T") | (!is.null(call$between.m)&&call$between.m=="T") ) {
        cat( "# Between-subject Effects\n" )
        cat( paste("Between-subjects Indirect Effect: ", round(t0[ 12], 3), " [", round(boot.ci( boot.object, index= 12, type="perc" )$perc[4], 3), ", " , round(boot.ci( boot.object, index= 12, type="perc" )$perc[5], 3), "]", "\n", sep="") )
        cat("\n")
      }
    }
    cat("\n")
    cat( "#### Total Effect ####\n" )
    if (((!is.null(call$random.a)&&call$random.a==T)||is.null(call$random.a)) && ((!is.null(call$random.b)&&call$random.b==T)||is.null(call$random.b))) {
      cat( paste("Unbiased Estimate of Total Effect: ", round(t0[ 14], 3), " [", round(boot.ci( boot.object, index= 14, type="perc" )$perc[4], 3), ", " , round(boot.ci( boot.object, index= 14, type="perc" )$perc[5], 3), "]", "\n", sep="") )
      cat( paste("Biased Total Effect of X on Y (c path): ", round(t0[ 1], 3), " [", round(boot.ci( boot.object, index= 1, type="perc" )$perc[4], 3), ", " , round(boot.ci( boot.object, index= 1, type="perc" )$perc[5], 3), "]", "\n", sep="") )
      cat( paste("Bias in Total Effect: ", round(t0[ 15], 3), " [", round(boot.ci( boot.object, index= 15, type="perc" )$perc[4], 3), ", " , round(boot.ci( boot.object, index= 15, type="perc" )$perc[5], 3), "]", "\n", sep="") )
    }
    else {
      cat( paste("Total Effect of X on Y (c path): ", round(t0[ 1], 3), " [", round(boot.ci( boot.object, index= 1, type="perc" )$perc[4], 3), ", " , round(boot.ci( boot.object, index= 1, type="perc" )$perc[5], 3), "]", "\n", sep="") )
    }
    cat("\n")
    cat("\n")
    cat( "#### Direct Effects ####\n" )
    cat( paste("Direct Effect of Predictor on Dependent Variable (c' path): ", round(t0[ 2], 3), " [", round(boot.ci( boot.object, index= 2, type="perc" )$perc[4], 3), ", " , round(boot.ci( boot.object, index= 2, type="perc" )$perc[5], 3), "]", "\n", sep="") )
    cat( paste("Within-subjects Effect of Predictor on Mediator (a path for group-mean centered predictor): ", round(t0[ 3], 3), " [", round(boot.ci( boot.object, index= 3, type="perc" )$perc[4], 3), ", " , round(boot.ci( boot.object, index= 3, type="perc" )$perc[5], 3), "]", "\n", sep="") )
    if ( (!is.null(call$between.x)&&call$between.x=="T") ) {cat( paste("Between-subjects Effect of Predictor on Mediator (a path for grand-mean centered average predictor): ", round(t0[ 4], 3), " [", round(boot.ci( boot.object, index= 4, type="perc" )$perc[4], 3), ", " , round(boot.ci( boot.object, index= 4, type="perc" )$perc[5], 3), "]", "\n", sep="") ) }
    cat( paste("Within-subjects Effect of Mediator on Dependent Variable (b path for group-mean centered mediator): ", round(t0[ 5], 3), " [", round(boot.ci( boot.object, index= 5, type="perc" )$perc[4], 3), ", " , round(boot.ci( boot.object, index= 5, type="perc" )$perc[5], 3), "]", "\n", sep="") )
    if (!is.null(call$between.m)&&call$between.m=="T") {cat( paste("Between-subjects Effect of Mediator on Dependent Variable (b path for grand-mean centered average mediator): ", round(t0[ 6], 3), " [", round(boot.ci( boot.object, index= 6, type="perc" )$perc[4], 3), ", " , round(boot.ci( boot.object, index= 6, type="perc" )$perc[5], 3), "]", "\n", sep="") ) }
  })
}
################################################################################