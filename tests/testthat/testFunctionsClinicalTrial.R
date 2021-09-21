context("FunctionsClinicalTrial")

# Data
data(vaccinated)
data(control)
tolerance <- 0.001

# ClinicalTrial
set.seed(1)
CT <- ClinicalTrial(vaccinated, control)

test_that("ClinicalTrial", 
          {expect_equal(CT$efficacy, 0.6071429, tolerance = tolerance)
           expect_equal(CT$confidenceInterval$lowerBound, 0.2152594, tolerance = tolerance)
           expect_equal(CT$confidenceInterval$upperBound, 0.8033277, tolerance = tolerance)
           expect_equal(length(vaccinated$diseaseStatus), floor(vaccinated$N))
           expect_equal(length(control$diseaseStatus), floor(control$N))
          })

# ClinicalTrialCoverage
set.seed(1)
CTCoverage <- ClinicalTrialCoverage(vaccinated, control)

test_that("ClinicalTrial", 
          {expect_equal(CTCoverage$efficacy, 0.6071429, tolerance = tolerance)
           expect_equal(CTCoverage$confidenceInterval95$lowerBound, 0.2152594, tolerance = tolerance)
           expect_equal(CTCoverage$confidenceInterval95$upperBound, 0.8033277, tolerance = tolerance)
           expect_equal(CTCoverage$confidenceInterval90$lowerBound, 0.2978740, tolerance = tolerance)
           expect_equal(CTCoverage$confidenceInterval90$upperBound, 0.7801866, tolerance = tolerance)
           expect_equal(CTCoverage$confidenceInterval80$lowerBound, 0.3823884, tolerance = tolerance)
           expect_equal(CTCoverage$confidenceInterval80$upperBound, 0.7501071, tolerance = tolerance)
           expect_equal(length(vaccinated$diseaseStatus), floor(vaccinated$N))
           expect_equal(length(control$diseaseStatus), floor(control$N))
          })

# WaldCI
ConfidenceIntervals95 <- waldCI(vaccinated, control)
ConfidenceIntervals70 <- waldCI(vaccinated, control, 0.70)

test_that("WaldCI",
          {expect_equal(ConfidenceIntervals95$lowerBound, 0.2152594, tolerance = tolerance)
           expect_equal(ConfidenceIntervals95$upperBound, 0.8033277, tolerance = tolerance)
           expect_equal(ConfidenceIntervals70$lowerBound, 0.4335844, tolerance = tolerance)
           expect_equal(ConfidenceIntervals70$upperBound, 0.7275203, tolerance = tolerance)})
