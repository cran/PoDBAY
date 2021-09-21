context("FunctionsImmunogenicitySample")

# Data
data(vaccinated)
data(control)
tolerance <- 0.001

# Assign disease status to the vaccinated and control population
set.seed(1)
CT <- ClinicalTrial(vaccinated, 
                    control)

# ExtractDiseased
diseased <- ExtractDiseased(vaccinated, 
                            control)

test_that("ExtractDiseased",
          {expect_is(diseased, "Population")
           expect_equal(diseased$N, 39)
           expect_equal(diseased$mean, 4.045505, tolerance = tolerance)
           expect_equal(diseased$stdDev, 1.837925, tolerance = tolerance)
           expect_length(diseased$titers, floor(diseased$N))
           expect_true(all(diseased$diseaseStatus))})


# ExtractNondiseased
nondiseased <- ExtractNondiseased(vaccinated, 
                                  control)

test_that("ExtractNondiseased",
          {expect_is(nondiseased, "Population")
           expect_equal(nondiseased$N, 1961)
           expect_equal(nondiseased$mean, 6.010405, tolerance = tolerance)
           expect_equal(nondiseased$stdDev, 2.29733, tolerance = tolerance)
           expect_length(nondiseased$titers, floor(nondiseased$N))
           expect_false(all(nondiseased$diseaseStatus))})

# ImmunogenicitySubset methodFULL
methodFULL = list(name = "Full", value = NA)

ImmunogenicitySample <- ImmunogenicitySubset(diseased = diseased,
                                             nondiseased = nondiseased,
                                             method = methodFULL)

test_that("ImmunogenicitySubset",
          {expect_is(ImmunogenicitySample, "Population")
           expect_equal(ImmunogenicitySample$N, (nondiseased$N + diseased$N))
           expect_equal(ImmunogenicitySample$mean, mean(c(nondiseased$titers, diseased$titers)), tolerance = tolerance)
           expect_equal(ImmunogenicitySample$stdDev, sd(c(nondiseased$titers, diseased$titers)), tolerance = tolerance)
           expect_equal(mean(ImmunogenicitySample$titers), ImmunogenicitySample$mean)})

# BlindSampling methodFULL
ImmunogenicitySample <- BlindSampling(diseased,
                                      nondiseased,
                                      method = methodFULL)

test_that("BlindSampling methodFULL",
          {expect_is(ImmunogenicitySample$ImmunogenicityNondiseased, "Population")
           expect_equal(nondiseased$N, ImmunogenicitySample$ImmunogenicityNondiseased$N)
           expect_equal(nondiseased$mean, ImmunogenicitySample$ImmunogenicityNondiseased$mean, tolerance = tolerance)
           expect_equal(nondiseased$stdDev, ImmunogenicitySample$ImmunogenicityNondiseased$stdDev, tolerance = tolerance)
           expect_equal(mean(nondiseased$titers), mean(ImmunogenicitySample$ImmunogenicityNondiseased$titers))
            
           expect_is(ImmunogenicitySample$ImmunogenicityVaccinated, "Population")
           expect_equal(ImmunogenicitySample$ImmunogenicityVaccinated$N, vaccinated$N)
           expect_equal(ImmunogenicitySample$ImmunogenicityVaccinated$mean, mean(vaccinated$titers), tolerance = tolerance)
           expect_equal(ImmunogenicitySample$ImmunogenicityVaccinated$stdDev, sd(vaccinated$titers), tolerance = tolerance)
           expect_equal(mean(ImmunogenicitySample$ImmunogenicityVaccinated$titers), mean(vaccinated$titers))
            
           expect_is(ImmunogenicitySample$ImmunogenicityControl, "Population")
           expect_equal(ImmunogenicitySample$ImmunogenicityControl$N, control$N)
           expect_equal(ImmunogenicitySample$ImmunogenicityControl$mean, mean(control$titers), tolerance = tolerance)
           expect_equal(ImmunogenicitySample$ImmunogenicityControl$stdDev, sd(control$titers), tolerance = tolerance)
           expect_equal(mean(ImmunogenicitySample$ImmunogenicityControl$titers), mean(control$titers))})


# ImmunogenicitySubset methodFIXED
value <- (vaccinated$N + control$N)/10
methodFIXED = list(name = "Fixed", value = value)

nondiseasedSample <- ImmunogenicitySubset(diseased = diseased,
                                          nondiseased = nondiseased,
                                          method = methodFIXED)

test_that("ImmunogenicitySubset",
          {expect_is(nondiseasedSample, "Population")
           expect_equal(nondiseasedSample$N, value)
           expect_is(nondiseasedSample$mean, "numeric")
           expect_is(nondiseasedSample$stdDev, "numeric")
           expect_length(nondiseasedSample$titers, value)})

# BlindSampling methodFIXED
ImmunogenicitySample <- BlindSampling(diseased,
                                      nondiseased,
                                      method = methodFIXED)

test_that("BlindSampling methodFIXED",
          {expect_is(ImmunogenicitySample$ImmunogenicityNondiseased, "Population")
           expect_lte(ImmunogenicitySample$ImmunogenicityNondiseased$N, value)
           expect_is(ImmunogenicitySample$ImmunogenicityNondiseased$mean, "numeric")
           expect_is(ImmunogenicitySample$ImmunogenicityNondiseased$stdDev, "numeric")
           expect_lte(length(ImmunogenicitySample$ImmunogenicityNondiseased$titers), value)
            
           expect_is(ImmunogenicitySample$ImmunogenicityVaccinated, "Population")
           expect_is(ImmunogenicitySample$ImmunogenicityVaccinated$mean, "numeric")
           expect_is(ImmunogenicitySample$ImmunogenicityVaccinated$stdDev, "numeric")
            
           expect_is(ImmunogenicitySample$ImmunogenicityControl, "Population")
           expect_is(ImmunogenicitySample$ImmunogenicityControl$mean, "numeric")
           expect_is(ImmunogenicitySample$ImmunogenicityControl$stdDev, "numeric")
            
           expect_equal(ImmunogenicitySample$ImmunogenicityControl$N + ImmunogenicitySample$ImmunogenicityVaccinated$N, value)
           expect_length(c(ImmunogenicitySample$ImmunogenicityVaccinated$titers, ImmunogenicitySample$ImmunogenicityControl$titers) , value)})

# ImmunogenicitySubset methodRatio
value <- 4
methodRATIO = list(name = "Ratio", value = value)

nondiseasedSample <- ImmunogenicitySubset(diseased = diseased,
                                          nondiseased = nondiseased,
                                          method = methodRATIO)

test_that("ImmunogenicitySubset methodRATIO",
          {expect_is(nondiseasedSample, "Population")
           expect_equal(nondiseasedSample$N, value * diseased$N)
           expect_is(nondiseasedSample$mean, "numeric")
           expect_is(nondiseasedSample$stdDev, "numeric")
           expect_length(nondiseasedSample$titers, value * diseased$N)})

# BlindSampling methodRATIO
ImmunogenicitySample <- BlindSampling(diseased,
                                      nondiseased,
                                      method = methodRATIO)

test_that("BlindSampling methodRATIO",
          {expect_is(ImmunogenicitySample$ImmunogenicityNondiseased, "Population")
           expect_equal(ImmunogenicitySample$ImmunogenicityNondiseased$N, value * diseased$N)
           expect_is(ImmunogenicitySample$ImmunogenicityNondiseased$mean, "numeric")
           expect_is(ImmunogenicitySample$ImmunogenicityNondiseased$stdDev, "numeric")
           expect_length(ImmunogenicitySample$ImmunogenicityNondiseased$titers, value * diseased$N)
            
           expect_is(ImmunogenicitySample$ImmunogenicityVaccinated, "Population")
           expect_is(ImmunogenicitySample$ImmunogenicityVaccinated$mean, "numeric")
           expect_is(ImmunogenicitySample$ImmunogenicityVaccinated$stdDev, "numeric")
            
           expect_is(ImmunogenicitySample$ImmunogenicityControl, "Population")
           expect_is(ImmunogenicitySample$ImmunogenicityControl$mean, "numeric")
           expect_is(ImmunogenicitySample$ImmunogenicityControl$stdDev, "numeric")
            
           expect_equal(ImmunogenicitySample$ImmunogenicityControl$N + ImmunogenicitySample$ImmunogenicityVaccinated$N, value * diseased$N)
           expect_length(c(ImmunogenicitySample$ImmunogenicityVaccinated$titers, ImmunogenicitySample$ImmunogenicityControl$titers) , value * diseased$N)})

# Errors expectations
methodFullError <- list(name = "Ful", value = NA)
methodFixedError <-  list(name = "FIXED", value = 300)
methodRatioError <-  list(name = " ", value = 4)

methodFixedError1 <-  list(name = "Fixed", value = NA)
methodFixedError2 <-  list(name = "Fixed", value = -10)
methodFixedError3 <-  list(name = "Fixed", value = 1000000)

methodRatioError1 <-  list(name = "Ratio", value = NA)
methodRatioError2 <-  list(name = "Ratio", value = 0)
methodRatioError3 <-  list(name = "Ratio", value = 1000)

test_that("Method errors",
          {expect_error(BlindSampling(diseased, nondiseased, method = methodFullError))
           expect_error(BlindSampling(diseased, nondiseased, method = methodFixedError))
           expect_error(BlindSampling(diseased, nondiseased, method = methodRatioError))
            
           expect_error(BlindSampling(diseased, nondiseased, method = methodFixedError1))
           expect_error(BlindSampling(diseased, nondiseased, method = methodFixedError2))
           expect_error(BlindSampling(diseased, nondiseased, method = methodFixedError3))
            
           expect_error(BlindSampling(diseased, nondiseased, method = methodRatioError1))
           expect_error(BlindSampling(diseased, nondiseased, method = methodRatioError2))
           expect_warning(BlindSampling(diseased, nondiseased, method = methodRatioError3))})
