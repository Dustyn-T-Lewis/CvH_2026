# Real-DB tests for the HPA-derived blood-contaminant filter.
# Uses the actual 00_input/HPA_annotations.tsv (no mocks).

library(testthat)
source(here::here("R", "blood_filter.R"))

hpa_path <- here::here("00_input", "HPA_annotations.tsv")
ref <- classify_blood_reference(hpa_path)

test_that("classification has the expected shape and vocabulary", {
  expect_true(all(c("gene", "uniprot", "secretome", "ery", "myo",
                    "reason", "verdict") %in% names(ref)))
  expect_setequal(unique(ref$verdict), c("keep", "remove"))
  expect_gt(nrow(ref), 10000)
})

test_that("muscle marker proteins are KEPT", {
  kept <- ref$verdict[match(c("MB", "CKM", "CA3", "ACTA1", "MYH2"), ref$gene)]
  expect_equal(unname(kept), rep("keep", 5))
})

test_that("muscle ECM / cytoskeleton are KEPT (not swept by a noisy arm)", {
  kept <- ref$verdict[match(c("COL1A1", "VIM", "FLNA"), ref$gene)]
  expect_equal(unname(kept), rep("keep", 3))
})

test_that("plasma, immunoglobulin and hemoglobin proteins are REMOVED", {
  expect_equal(ref$verdict[ref$gene == "ALB"],   "remove")
  expect_equal(ref$verdict[ref$gene == "C3"],    "remove")
  expect_equal(ref$verdict[ref$gene == "TF"],    "remove")
  expect_equal(ref$verdict[ref$gene == "PPBP"],  "remove")   # secreted platelet factor
  expect_equal(ref$verdict[ref$gene == "HBB"],   "remove")
  expect_equal(ref$verdict[ref$gene == "HBA1"],  "remove")
  expect_equal(ref$verdict[ref$gene == "IGLC2"], "remove")
})

test_that("reasons map to the correct arm", {
  expect_match(ref$reason[ref$gene == "HBB"],   "erythrocyte")
  expect_match(ref$reason[ref$gene == "IGLC2"], "immunoglobulin")
  expect_match(ref$reason[ref$gene == "C3"],    "secreted-to-blood")
})

test_that("A2M is rescued by the Myonuclei keep-override", {
  expect_equal(ref$verdict[ref$gene == "A2M"], "keep")
  expect_match(ref$reason[ref$gene == "A2M"], "rescued")
})

test_that("blood_contaminant_genes returns a unique character vector", {
  genes <- blood_contaminant_genes(hpa_path)
  expect_type(genes, "character")
  expect_false(any(duplicated(genes)))
  expect_true("HBB" %in% genes)
  expect_false("MB" %in% genes)
})
