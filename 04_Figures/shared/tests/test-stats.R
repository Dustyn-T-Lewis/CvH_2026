library(testthat)
setwd(here::here())
source("04_Figures/shared/stats.R")

test_that("fisher_z_ci matches the interval stats::cor.test computes", {
  x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  y <- c(2, 1, 4, 3, 6, 5, 8, 7, 10, 9)
  ct <- cor.test(x, y)
  expect_equal(
    unname(fisher_z_ci(ct$estimate, length(x))),
    as.numeric(ct$conf.int),
    tolerance = 1e-8
  )
})

test_that("fisher_z_ci returns lo and hi bracketing r", {
  ci <- fisher_z_ci(0.5, 28)
  expect_named(ci, c("lo", "hi"))
  expect_equal(unname(ci), c(0.1560284, 0.7358185), tolerance = 1e-6)
  expect_lt(ci[["lo"]], 0.5)
  expect_gt(ci[["hi"]], 0.5)
})

test_that("fisher_z_ci is antisymmetric in the sign of r", {
  pos <- fisher_z_ci(0.5, 28)
  neg <- fisher_z_ci(-0.5, 28)
  expect_equal(neg[["lo"]], -pos[["hi"]], tolerance = 1e-12)
  expect_equal(neg[["hi"]], -pos[["lo"]], tolerance = 1e-12)
  expect_equal(unname(fisher_z_ci(0, 28)), c(-0.3730769, 0.3730769), tolerance = 1e-6)
})

test_that("fisher_z_ci widens with a higher level and narrows with a larger n", {
  expect_gt(
    diff(fisher_z_ci(0.5, 28, level = 0.99)),
    diff(fisher_z_ci(0.5, 28, level = 0.95))
  )
  expect_lt(diff(fisher_z_ci(0.5, 200)), diff(fisher_z_ci(0.5, 28)))
})

test_that("fisher_z_ci degenerates to the full range at n = 3", {
  expect_equal(unname(fisher_z_ci(0.5, 3)), c(-1, 1))
})

test_that("classify_proteins_f4 assigns each of the four classes", {
  cls <- classify_proteins_f4(
    pi_CvH = c(0.01, 0.01, 0.90, 0.90),
    pi_TR  = c(0.01, 0.90, 0.01, 0.90)
  )
  expect_equal(
    as.character(cls),
    c("Sig Both", "Sig Cancer only", "Sig Training only", "NS")
  )
})

test_that("classify_proteins_f4 keeps the documented level order", {
  expect_equal(
    levels(classify_proteins_f4(0.5, 0.5)),
    c("Sig Both", "Sig Cancer only", "Sig Training only", "NS")
  )
})

test_that("classify_proteins_f4 treats the threshold as strict and honours the argument", {
  expect_equal(as.character(classify_proteins_f4(0.05, 0.05)), "NS")
  expect_equal(as.character(classify_proteins_f4(0.049, 0.9)), "Sig Cancer only")
  expect_equal(
    as.character(classify_proteins_f4(0.2, 0.2, threshold = 0.5)),
    "Sig Both"
  )
})

test_that("classify_proteins_f4 sends non-estimable Pi values to NS, not to a class", {
  expect_equal(as.character(classify_proteins_f4(NA_real_, NA_real_)), "NS")
  expect_equal(as.character(classify_proteins_f4(NA_real_, 0.01)), "Sig Training only")
  expect_equal(as.character(classify_proteins_f4(0.01, NA_real_)), "Sig Cancer only")
})
