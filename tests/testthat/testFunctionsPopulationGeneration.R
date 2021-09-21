context("FunctionsPopulationGeneration")

# Data
N <- 10.8
mean <- 5
stdDev <- 2

# generatePopulation
population <- generatePopulation(N = N,
                                 mean = mean,
                                 stdDev = stdDev)

population$assignPoD(rep(0.3,N))
population$diseaseStatus <- c(rep(TRUE, round(N / 3)), rep(FALSE, N - round(N / 3)))

test_that("generatePopulation", 
          {expect_is(population, "Population")
           expect_is(population$popFun(), "function")
           expect_equal(population$N, N)
           expect_equal(population$mean, mean)
           expect_equal(population$stdDev, stdDev)
           expect_equal(length(population$titers), floor(N))
           expect_equal(length(population$PoDs), floor(N))
           expect_equal(length(population$diseaseStatus), floor(N))
           expect_equal(population$getDiseasedCount(), round(N / 3))
           expect_equal(population$getNondiseasedCount(), floor(N - round(N / 3)))
           expect_equal(length(population$getDiseasedTiters()), round(N / 3))
           expect_equal(length(population$getNondiseasedTiters()), floor(N - round(N / 3)))
           })


# PoD
PoDs <- PoD(-10:10, 0.05, 5, 7, FALSE, log2(10), log2(5))
cppPoDs <- cppPoD(-10:10, 0.05, 5, 7, FALSE, log2(10), log2(5))

PoDsTRUE <- PoD(-10:10, 0.05, 5, 7, TRUE, log2(10), log2(5))
cppPoDsTRUE <- cppPoD(-10:10, 0.05, 5, 7, TRUE, log2(10), log2(5))

test_that("PoD",
          {expect_true(all(PoDs >= 0 ))
           expect_true(all(PoDsTRUE >= 0 ))
           expect_equal(PoDs, cppPoDs)
           expect_equal(PoDsTRUE, cppPoDsTRUE)})
