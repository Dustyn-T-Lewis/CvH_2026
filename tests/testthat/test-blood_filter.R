# Real-DB tests for the HPA-derived blood-contaminant filter.
# Uses the actual 00_input/HPA_blood_annotations.tsv (no mocks).

library(testthat)
source(here::here("R", "blood_filter.R"))

hpa_path <- here::here("00_input", "HPA_blood_annotations.tsv")
ref <- classify_blood_reference(hpa_path)

test_that("classification has the expected shape and vocabulary", {
  expect_true(all(c("gene", "uniprot", "secretome", "ih_myocyte",
                    "ms_blood", "reason", "verdict") %in% names(ref)))
  expect_setequal(unique(ref$verdict), c("keep", "remove"))
  expect_gt(nrow(ref), 4000)
})

test_that("muscle marker proteins are KEPT", {
  kept <- ref$verdict[match(c("MB", "CKM", "CA3"), ref$gene)]
  expect_equal(unname(kept), rep("keep", 3))
})

test_that("blood/plasma proteins are REMOVED", {
  removed <- ref$verdict[match(c("HBB", "HBA1", "HBD", "TF", "C3", "APOA1"), ref$gene)]
  expect_equal(unname(removed), rep("remove", 6))
})

test_that("hemoglobins are flagged as erythrocyte, plasma proteins as secreted", {
  expect_match(ref$reason[ref$gene == "HBB"], "erythrocyte")
  expect_match(ref$reason[ref$gene == "C3"],  "secreted-to-blood")
})

test_that("muscle ECM and cytoskeleton are NOT removed", {
  # the prototype's earlier ratio rule wrongly removed these; the secretome rule keeps them
  ecm <- ref$verdict[match(c("COL1A1", "COL6A3", "VIM", "FLNA"), ref$gene)]
  expect_true(all(ecm == "keep", na.rm = TRUE))
})

test_that("blood_contaminant_genes returns a non-empty unique character vector", {
  genes <- blood_contaminant_genes(hpa_path)
  expect_type(genes, "character")
  expect_false(any(duplicated(genes)))
  expect_true("HBB" %in% genes)
  expect_false("MB" %in% genes)
})
