context("FunctionsPoDParameters")

# Data
data(diseased)
data(nondiseased)
data(PoDParams)
tolerance <- 0.001

# GenerateNondiseased
titersMockup <- 1:10
nondiseasedGen1 <- GenerateNondiseased(blindNondiseasedTiters = titersMockup, nondiseasedCount = length(titersMockup) + 0.01)
nondiseasedGen2 <- GenerateNondiseased(blindNondiseasedTiters = titersMockup, nondiseasedCount = length(titersMockup))
nondiseasedGen3 <- GenerateNondiseased(blindNondiseasedTiters = titersMockup, nondiseasedCount = length(titersMockup) - 0.01)

test_that("GenerateNondiseased",
          {expect_equal(nondiseasedGen1,titersMockup)
           expect_equal(nondiseasedGen2,titersMockup)
           expect_length(nondiseasedGen3, length(titersMockup) - 1)
           expect_true(any(!is.na( match( nondiseasedGen3, titersMockup))))})

# PoDMLE
PoDMLEOk <- PoDMLE(nondiseased$titers, diseased$titers)
PoDMLEOklowTiterPercent <- PoDMLE(nondiseased$titers, diseased$titers, lowTiterPercent = 0)
PoDMLEOkinitialSlope <- PoDMLE(nondiseased$titers, diseased$titers, initialSlope = 35)
PoDMLENull <- PoDMLE(0, nondiseased$titers)
PoDMLEOk1 <- PoDMLE(nondiseased$titers, diseased$titers, adjustTiters = TRUE, adjustFrom = 7, adjustTo = log2(5))

test_that("PoDMLE",
          {expect_named(PoDMLEOk, c("pmax", "et50", "slope"), ignore.order = TRUE)
           expect_false(PoDMLEOk$et50 == PoDMLEOk1$et50)
           expect_false(PoDMLEOk$pmax == PoDMLEOk1$pmax)
           expect_false(PoDMLEOk$slope == PoDMLEOk1$slope)
           expect_is(PoDMLEOk, "list")
           expect_is(PoDMLEOklowTiterPercent, "list")
           expect_is(PoDMLEOkinitialSlope, "list")
           expect_null(PoDMLENull)
           })

# cppMLE
data(diseased)
data(nondiseased)
PoDParamsList <- list("pmax" = 0.05, "et50" = 5, "slope" = 7)
cppMLEList <- cppMLE(PoDParamsList, nondiseased$titers, diseased$titers)
cppMLEOk <- cppMLE(PoDParams, nondiseased$titers, diseased$titers, adjustTiters = FALSE)

PoDParamsError1 <- data.frame("Pmax" = 0.05, "et50" = 5, "slope" = 7)
PoDParamsError2 <- data.frame("pmax" = 0.05, 5, "slope" = 7)
PoDParamsError3 <- data.frame("pmax" = 0.05, "et50" = 5, " " = 7)

cppMLEOk1 <- cppMLE(PoDParams, nondiseased$titers, diseased$titers, adjustTiters = TRUE, adjustFrom = log2(10), adjustTo = log2(5))

test_that("cppMLE",
          {expect_false(cppMLEOk == cppMLEOk1)
           expect_equal(cppMLEOk, cppMLEList)
           expect_error(cppMLE(PoDParamsError1, nondiseased$titers, diseased$titers))
           expect_error(cppMLE(PoDParamsError2, nondiseased$titers, diseased$titers))
           expect_error(cppMLE(PoDParamsError3, nondiseased$titers, diseased$titers))
          })

# MLE
data(diseased)
data(nondiseased)
MLEOk <- MLE(PoDParams, nondiseased$titers, diseased$titers)
MLEOk1 <- MLE(PoDParams, nondiseased$titers, diseased$titers, adjustTiters = TRUE, adjustFrom = log2(10), adjustTo = log2(5))

PoDParamsList <- list("pmax" = 0.05, "et50" = 5, "slope" = 7)
MLEList <- MLE(PoDParamsList, nondiseased$titers, diseased$titers)

PoDParamsError1 <- data.frame("Pmax" = 0.05, "et50" = 5, "slope" = 7)
PoDParamsError2 <- data.frame("pmax" = 0.05, 5, "slope" = 7)
PoDParamsError3 <- data.frame("pmax" = 0.05, "et50" = 5, " " = 7)

test_that("MLE",
          {expect_false(MLEOk == MLEOk1)
           expect_equal(MLEOk, MLEList)
           expect_error(MLE(PoDParamsError1, nondiseased$titers, diseased$titers))
           expect_error(MLE(PoDParamsError2, nondiseased$titers, diseased$titers))
           expect_error(MLE(PoDParamsError3, nondiseased$titers, diseased$titers))
          })

test_that("cppMLE - MLE",
          {expect_equal(MLEOk, cppMLEOk)
           expect_equal(MLEOk1, cppMLEOk1)
           expect_equal(MLEList, cppMLEList)
          })


# PoDParamEstimation
data(diseased)
data(nondiseased)
repeatCount <- 2

# Full method
set.seed(2)
PoDParamEst <- PoDParamEstimation(diseased$titers,
                                  nondiseased$titers,
                                  nondiseased$N,
                                  repeatCount = repeatCount,
                                  adjustTiters = FALSE)

test_that("PoDParamEstimation names - FULL method",
          {expect_false(any(is.na(match(names(PoDParamEst), c("results", "resultsPriorReset", "failCount")))))
           #results
           expect_is(PoDParamEst$results,"data.frame")
           expect_named(PoDParamEst$results, c("pmax", "et50", "slope"), ignore.order = TRUE)
           expect_gt(max(PoDParamEst$results$pmax) - min(PoDParamEst$results$pmax),0)
           expect_gt(max(PoDParamEst$results$et50) - min(PoDParamEst$results$et50),0)
           expect_gt(max(PoDParamEst$results$slope) - min(PoDParamEst$results$slope),0)
           #resultsPriorReset
           expect_is(PoDParamEst$resultsPriorReset,"data.frame") 
           expect_named(PoDParamEst$resultsPriorReset, c("pmax", "et50", "slope"), ignore.order = TRUE)
           expect_equal(max(PoDParamEst$resultsPriorReset$pmax) - min(PoDParamEst$resultsPriorReset$pmax), 0)
           expect_equal(max(PoDParamEst$resultsPriorReset$et50) - min(PoDParamEst$resultsPriorReset$et50), 0)
           expect_equal(max(PoDParamEst$resultsPriorReset$slope) - min(PoDParamEst$resultsPriorReset$slope), 0)
           #failcount
           expect_lte(PoDParamEst$failCount, repeatCount)
          })

Rversion_short <- paste(R.version$major, R.version$minor, sep = ".")

if (Rversion_short >= "3.6.0") {
  
  # R version 3.6.3
  test_that("PoDParamEstimation values - FULL method - tests might fail if method changes",
            {expect_equal(PoDParamEst$results$pmax, c(0.05241564, 0.04800023), tolerance = tolerance)
             expect_equal(PoDParamEst$results$et50, c(5.397124, 5.684373), tolerance = tolerance)
             expect_equal(PoDParamEst$results$slope, c(11.37913, 20.19611), tolerance = tolerance)
             #resultsPriorReset
             expect_equal(PoDParamEst$resultsPriorReset$pmax, c(0.03429965, 0.03429965), tolerance = tolerance)
             expect_equal(PoDParamEst$resultsPriorReset$et50, c(6.051441, 6.051441), tolerance = tolerance)
             expect_equal(PoDParamEst$resultsPriorReset$slope, c(28.49322, 28.49322), tolerance = tolerance)
             #failcount
             expect_equal(PoDParamEst$failCount, 0)
            })
} else {
  # R version <3.6.0 (different seeding)
  test_that("PoDParamEstimation values - FULL method - tests might fail if method changes",
            {expect_equal(PoDParamEst$results$pmax, c(0.03714693, 0.03607800), tolerance = tolerance)
             expect_equal(PoDParamEst$results$et50, c(6.175141, 6.280964), tolerance = tolerance)
             expect_equal(PoDParamEst$results$slope, c(21.05802, 28.63538), tolerance = tolerance)
             #resultsPriorReset
             expect_equal(PoDParamEst$resultsPriorReset$pmax, c(0.03429965, 0.03429965), tolerance = tolerance)
             expect_equal(PoDParamEst$resultsPriorReset$et50, c(6.051441, 6.051441), tolerance = tolerance)
             expect_equal(PoDParamEst$resultsPriorReset$slope, c(28.49322, 28.49322), tolerance = tolerance)
             #failcount
             expect_equal(PoDParamEst$failCount, 0)
            })
}

# diseased and unifected titers switched
set.seed(1)
PoDParamEst1 <- PoDParamEstimation(nondiseased$titers,
                                   diseased$titers,
                                   nondiseased$N,
                                   repeatCount = repeatCount)

test_that("PoDParamEstimation names - diseased and unifected titers switched",
          {expect_equal(nrow(PoDParamEst1$results),0)
           expect_equal(nrow(PoDParamEst1$resultsPriorReset),0)
           expect_equal(PoDParamEst1$failCount, repeatCount)
          })

# Method = Ratio - immunogenicity ratio
set.seed(1)
NondiseasedImmunogenicitySubset <- ImmunogenicitySubset(diseased, nondiseased, method = list(name = "Ratio", value = 4))

# Number of Nondiseased patient identified in the Clinical Trial
nondiseasedGenerationCount <- nondiseased$N

PoDParamEstImuno <- PoDParamEstimation(diseased$titers,
                                       NondiseasedImmunogenicitySubset$titers,
                                       nondiseased$N,
                                       repeatCount = repeatCount)

test_that("PoDParamEstimation general - RATIO method",
          {#results
           expect_gt(max(PoDParamEstImuno$results$pmax) - min(PoDParamEstImuno$results$pmax),0)
           expect_gt(max(PoDParamEstImuno$results$et50) - min(PoDParamEstImuno$results$et50),0)
           expect_gt(max(PoDParamEstImuno$results$slope) - min(PoDParamEstImuno$results$slope),0)
           #resultsPriorReset
           expect_gt(max(PoDParamEstImuno$resultsPriorReset$pmax) - min(PoDParamEstImuno$resultsPriorReset$pmax),0)
           expect_gt(max(PoDParamEstImuno$resultsPriorReset$et50) - min(PoDParamEstImuno$resultsPriorReset$et50),0)
           expect_gt(max(PoDParamEstImuno$resultsPriorReset$slope) - min(PoDParamEstImuno$resultsPriorReset$slope),0)
           #failcount
           expect_lte(PoDParamEst$failCount, repeatCount)
          })

if (Rversion_short >= "3.6.0") {
  # R version >= 3.6.3 
  test_that("PoDParamEstimation values - RATIO method - tests might fail if method changes",
            {expect_equal(PoDParamEstImuno$results$pmax, c(0.05668154, 0.09605827), tolerance = tolerance)
             expect_equal(PoDParamEstImuno$results$et50, c( 5.009985, 3.996407), tolerance = tolerance)
             expect_equal(PoDParamEstImuno$results$slope, c(7.760187, 5.196349), tolerance = tolerance)
             #resultsPriorReset
             expect_equal(PoDParamEstImuno$resultsPriorReset$pmax, c(0.05454088, 0.07329043), tolerance = tolerance)
             expect_equal(PoDParamEstImuno$resultsPriorReset$et50, c(4.612283, 4.010578), tolerance = tolerance)
             expect_equal(PoDParamEstImuno$resultsPriorReset$slope, c(6.099464, 5.060878), tolerance = tolerance)
             #failcount
             expect_equal(PoDParamEstImuno$failCount, 0)
            })
} else {
  # R version < 3.6.0 (different seeding)
  test_that("PoDParamEstimation values - RATIO method - tests might fail if method changes",
            {expect_equal(PoDParamEstImuno$results$pmax, c(0.08196423, 0.09766574), tolerance = tolerance)
             expect_equal(PoDParamEstImuno$results$et50, c(3.419386, 2.974790), tolerance = tolerance)
             expect_equal(PoDParamEstImuno$results$slope, c(4.540833, 4.859419), tolerance = tolerance)
             #resultsPriorReset
             expect_equal(PoDParamEstImuno$resultsPriorReset$pmax, c(0.07597830, 0.08971953), tolerance = tolerance)
             expect_equal(PoDParamEstImuno$resultsPriorReset$et50, c(3.890352, 3.555434), tolerance = tolerance)
             expect_equal(PoDParamEstImuno$resultsPriorReset$slope, c(4.930218, 4.499795), tolerance = tolerance)
             #failcount
             expect_equal(PoDParamEstImuno$failCount, 0)
            })
}


# Warning test
test_that("PoDParamEstimation warning",
          {expect_warning(PoDParamEstimation(diseased$titers,
                                             nondiseased$titers,
                                             500,
                                             repeatCount = 1))
          })

# adjustTiters = TRUE
set.seed(1)
PoDParamEstFALSE <- PoDParamEstimation(diseased$titers,
                                       nondiseased$titers,
                                       nondiseased$N,
                                       repeatCount = repeatCount,
                                       adjustTiters = FALSE,
                                       adjustFrom = log2(100),
                                       adjustTo = log2(5))


set.seed(1)
PoDParamEstTRUE <- PoDParamEstimation(diseased$titers,
                                      nondiseased$titers,
                                      nondiseased$N,
                                      repeatCount = repeatCount,
                                      adjustTiters = TRUE,
                                      adjustFrom = log2(100),
                                      adjustTo = log2(5))

set.seed(1)
PoDParamEstTRUE50 <- PoDParamEstimation(diseased$titers,
                                        nondiseased$titers,
                                        nondiseased$N,
                                        repeatCount = repeatCount,
                                        adjustTiters = TRUE,
                                        adjustFrom = log2(50),
                                        adjustTo = log2(5))

test_that("PoDParamEstimation adjustTiters = TRUE, FALSE",
          {expect_false(any(colSums(PoDParamEstFALSE$results - PoDParamEstTRUE$results) == 0))
           expect_false(any(colSums(PoDParamEstFALSE$resultsPriorReset - PoDParamEstTRUE$resultsPriorReset) == 0))
           expect_false(any(colSums(PoDParamEstTRUE50$results - PoDParamEstTRUE$results) == 0))
           expect_false(any(colSums(PoDParamEstTRUE50$resultsPriorReset - PoDParamEstTRUE$resultsPriorReset) == 0))
          })

# PoDParamsCI
df <- data.frame("pmax" = seq(from = 0, to = 1, by = 0.01),
                 "et50" = seq(from = 5, to = 6, by = 0.01),
                 "slope" = 0:100)

PoDParamsCIres <- PoDParamsCI(df)

test_that("PoDParamsCI",
          {expect_equal(PoDParamsCIres$PmaxCILow, 0.025)
           expect_equal(PoDParamsCIres$PmaxCIHigh, 0.975)
           expect_equal(PoDParamsCIres$Et50CILow, 5.025)
           expect_equal(PoDParamsCIres$Et50CIHigh, 5.975)
           expect_equal(PoDParamsCIres$SlopeCILow, 2.5)
           expect_equal(PoDParamsCIres$SlopeCIHigh, 97.5)})

# PoDParamsCICoverage
PoDParamsCICoverageRes <- PoDParamsCICoverage(df)

test_that("PoDParamsCICoverage",
          {expect_equal(PoDParamsCICoverageRes$PmaxCILow95, 0.025)
           expect_equal(PoDParamsCICoverageRes$PmaxCIHigh95, 0.975)
           expect_equal(PoDParamsCICoverageRes$PmaxCILow90, 0.05)
           expect_equal(PoDParamsCICoverageRes$PmaxCIHigh90, 0.95)
           expect_equal(PoDParamsCICoverageRes$PmaxCILow80, 0.1)
           expect_equal(PoDParamsCICoverageRes$PmaxCIHigh80, 0.9)
            
           expect_equal(PoDParamsCICoverageRes$Et50CILow95, 5.025)
           expect_equal(PoDParamsCICoverageRes$Et50CIHigh95, 5.975)
           expect_equal(PoDParamsCICoverageRes$Et50CILow90, 5.05)
           expect_equal(PoDParamsCICoverageRes$Et50CIHigh90, 5.95)
           expect_equal(PoDParamsCICoverageRes$Et50CILow80, 5.1)
           expect_equal(PoDParamsCICoverageRes$Et50CIHigh80, 5.9)
            
           expect_equal(PoDParamsCICoverageRes$SlopeCILow95, 2.5)
           expect_equal(PoDParamsCICoverageRes$SlopeCIHigh95, 97.5)
           expect_equal(PoDParamsCICoverageRes$SlopeCILow90, 5)
           expect_equal(PoDParamsCICoverageRes$SlopeCIHigh90, 95)
           expect_equal(PoDParamsCICoverageRes$SlopeCILow80, 10)
           expect_equal(PoDParamsCICoverageRes$SlopeCIHigh80, 90)})

# PoDParamPointEstimation
data(estimatedParameters)

PoDParamPointEstimationRes <- PoDParamPointEstimation(estimatedParameters$resultsPriorReset)

test_that("PoDParamPointEstimation",
          {expect_true(all.equal(PoDParamPointEstimationRes, c( 0.0766969, 5.928621, 4.200404), check.attributes = FALSE, tolerance = tolerance))
           expect_named(PoDParamPointEstimationRes, c("pmax", "et50", "slope"), ignore.order = TRUE)
          })

# fitPoD
data(PoDParams)

TitersInput <- seq(from = 0, to = 20, by = 0.01)

CurveTitersMedian <- PoD(TitersInput, 0.05, 6, 7)

fitPoDRes <- fitPoD(PoDParams, TitersInput, CurveTitersMedian)

test_that("fitPoDRes", 
          {expect_equal(fitPoDRes, -0.005119542)
          })

# PoDCI
PoDCIpmax <- PoDCI(df$pmax)

test_that("PoDCI",
          {expect_true(all.equal(PoDCIpmax, c(0.025, 0.5, 0.975), check.attributes = FALSE, tolerance = tolerance))})

# PoDCurvePlot
PoDCurvePlotRes <- PoDCurvePlot(-5:10, estimatedParameters)

test_that("PoDCurvePlot",
          {expect_is(PoDCurvePlotRes, "ggplot")
           })

# efficacySquaredError
data(vaccinated)
data(control)

TrueEfficacy <- 0.5

effSE <- efficacySquaredError(PoDParams, TrueEfficacy, titerFun = list(vaccinated$popFun(), control$popFun()))
effSEAdjustTiters <- efficacySquaredError(PoDParams, TrueEfficacy, titerFun = list(vaccinated$popFun(), control$popFun()), adjustTiters = TRUE, adjustFrom = log2(10), adjustTo = log2(5))
effSEAdjustTitersFALSE <- efficacySquaredError(PoDParams, TrueEfficacy, titerFun = list(vaccinated$popFun(), control$popFun()), adjustTiters = FALSE, adjustFrom = log2(10), adjustTo = log2(5))
effSESwitch <- efficacySquaredError(PoDParams, TrueEfficacy, titerFun = list(control$popFun(),vaccinated$popFun()))

test_that("efficacySquaredError",
          {expect_true(all.equal(effSE, effSEAdjustTitersFALSE))
            expect_false(isTRUE(all.equal(effSEAdjustTiters, effSEAdjustTitersFALSE)))
            expect_gt(effSESwitch, effSE)
            })

# PoDEfficacySquaredError
PoDEfficacySquaredErrorRes <- PoDEfficacySquaredError(TrueEfficacy, vaccinated, control)

PoDEfficacySquaredErrorResAdjustTiters <- PoDEfficacySquaredError(TrueEfficacy, vaccinated, control, adjustTiters = TRUE, adjustFrom = log2(10), adjustTo = log2(5))
PoDEfficacySquaredErrorResAdjustTiters1 <- PoDEfficacySquaredError(TrueEfficacy, vaccinated, control, adjustTiters = TRUE, adjustFrom = log2(10), adjustTo = 0)
PoDEfficacySquaredErrorResAdjustTitersFALSE <- PoDEfficacySquaredError(TrueEfficacy, vaccinated, control, adjustTiters = FALSE, adjustFrom = log2(10), adjustTo = log2(5))

test_that("PoDEfficacySquaredError",
          {expect_true(all.equal(PoDEfficacySquaredErrorRes, c(5.511748, 6.113515), check.attributes = FALSE, tolerance = tolerance))
           expect_true(all.equal(PoDEfficacySquaredError(TrueEfficacy, control, vaccinated), c(25.51815, 22.88361), check.attributes = FALSE, tolerance = tolerance))
           expect_warning(PoDEfficacySquaredError(TrueEfficacy, control, vaccinated))
           expect_true(all.equal(PoDEfficacySquaredErrorRes, PoDEfficacySquaredErrorResAdjustTitersFALSE))
           expect_false(isTRUE(all.equal(PoDEfficacySquaredErrorResAdjustTitersFALSE, PoDEfficacySquaredErrorResAdjustTiters)))
           expect_false(isTRUE(all.equal(PoDEfficacySquaredErrorResAdjustTiters1, PoDEfficacySquaredErrorResAdjustTiters)))
          })

# PmaxEstimation
IncidenceRate <- 0.2
params <- list("et50" = 4, "slope" = 6)
paramsError <- list("Et50" = 4, "slope" = 6)
paramsError2 <- list("Et50" = 4, "slope" = 6)

PmaxEstimationRes <- PmaxEstimation(IncidenceRate, params, control)

test_that("PmaxEstimation",
  {expect_equal(PmaxEstimationRes$pmax, 0.5650941, tolerance = tolerance)
   expect_error(PmaxEstimation(IncidenceRate, paramsError, control))
   expect_error(PmaxEstimation(IncidenceRate, paramsError2, control))})

