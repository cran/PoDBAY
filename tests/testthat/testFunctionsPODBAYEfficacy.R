context("FunctionsPoDBAYEfficacy")

# Data
data(vaccinated)
data(control)
data(estimatedParameters)
data(PoDParams)
tolerance <- 0.001

# PoDBAYEfficacy
set.seed(1)
PoDBAYEfficacyRes <- PoDBAYEfficacy(estimatedParameters$results,
                                    vaccinated,
                                    control)

estimatedParametersError <- estimatedParameters
names(estimatedParametersError$results) <- c("Pmax", "slope", "et50")

test_that("PoDBAYEfficacy",
          {expect_equal(PoDBAYEfficacyRes, c(0.6337973, 0.5618014, 0.6759340, 0.6733055, 0.7199953, 0.7312204, 0.5784680, 0.7561478, 0.6694463, 0.6722989), tolerance = tolerance)
           expect_warning(PoDBAYEfficacy(estimatedParameters$results, control, vaccinated))
           expect_error(PoDBAYEfficacy(estimatedParametersError$results, vaccinated, control))})

# EfficacyCI
EfficacyCIRes <- EfficacyCI(0:100)
EfficacyCIRes80 <- EfficacyCI(0:100, ci = 0.8)

test_that("EfficacyCI",
          {expect_true(all.equal(EfficacyCIRes, c(50, 50, 2.5, 97.5), check.attributes = FALSE))
           expect_true(all.equal(EfficacyCIRes80, c(50, 50, 10, 90), check.attributes = FALSE))})

# EfficacyCICoverage
EfficacyCICoverageRes <- EfficacyCICoverage(0:100)

test_that("EfficacyCICoverage",
          {expect_true(all.equal(EfficacyCICoverageRes, c(50.0, 50.0, 2.5, 97.5, 5.0, 95.0, 10.0, 90.0), check.attributes = FALSE))
            })

# JitterMean
set.seed(1)

JitterMeanControl <- JitterMean(control)
JitterMeanVaccinated <- JitterMean(vaccinated)

test_that("JitterMean",
          {expect_equal(JitterMeanControl, 4.96038, tolerance = tolerance)
           expect_equal(JitterMeanVaccinated, 7.011615, tolerance = tolerance)})

# efficacyComputation
means <- list(vaccinated = vaccinated$mean, control = control$mean)
standardDeviations <- list(vaccinated = vaccinated$stdDev, control = control$stdDev)

meansError1 <- list(vacciated = vaccinated$mean, control = control$mean)
meansError2 <- list(vaccinated = vaccinated$mean, Control = control$mean)

standardDeviationsError1 <- list(Vaccinated = vaccinated$stdDev, control = control$stdDev)
standardDeviationsError2 <- list(vaccinated = vaccinated$stdDev, cntrol = control$stdDev)

PoDParamsError <- c(pmax = 0.05, Et50 = 5, slope = 7)
  
efficacyComputationRes <- efficacyComputation(PoDParams, means, standardDeviations)
efficacyComputationResAdjustTiters <- efficacyComputation(PoDParams, means, standardDeviations, adjustTiters = TRUE, adjustFrom = log2(10), adjustTo = log2(5))
efficacyComputationResAdjustTitersFALSE <- efficacyComputation(PoDParams, means, standardDeviations, adjustTiters = FALSE, adjustFrom = log2(10), adjustTo = log2(5))

test_that("efficacyComputation",
          {expect_equal(efficacyComputationRes, 0.5821759, tolerance = tolerance)
           expect_false(isTRUE(all.equal(efficacyComputationRes,efficacyComputationResAdjustTiters)))
           expect_equal(efficacyComputationRes,efficacyComputationResAdjustTitersFALSE)
           expect_error(efficacyComputation(PoDParamsError, means, standardDeviations))
           expect_error(efficacyComputation(PoDParams, meansError1, standardDeviations))
           expect_error(efficacyComputation(PoDParams, means, standardDeviationsError1))
           expect_error(efficacyComputation(PoDParams, meansError2, standardDeviations))
           expect_error(efficacyComputation(PoDParams, means, standardDeviationsError2))})

# ExpectedPoD
funPoD <- function(x) PoD(x, pmax = PoDParams$pmax, et50 = PoDParams$et50, slope = PoDParams$slope)
fun1 <- function(x) dnorm(x, mean = 6, sd = 2)
fun2 <- function(x) dnorm(x, mean = 5.99, sd = 2)
fun3 <- function(x) dnorm(x, mean = 6.01, sd = 2)

ExpectedPoDFun1 <- ExpectedPoD(funPoD, fun1)
ExpectedPoDFun2 <- ExpectedPoD(funPoD, fun2)
ExpectedPoDFun3 <- ExpectedPoD(funPoD, fun3)

test_that("ExpectedPoD",
          {expect_lt(ExpectedPoDFun1, ExpectedPoDFun2)
           expect_lt(ExpectedPoDFun3, ExpectedPoDFun1)})
