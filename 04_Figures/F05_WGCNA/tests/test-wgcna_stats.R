library(testthat)
source(here::here("04_Figures/shared/wgcna_stats.R"))

test_that("run_module_fgsea returns one row per module x contrast with bounded NES", {
  set.seed(1)
  genes <- paste0("G", 1:200)
  rank_wide <- data.frame(
    gene = genes,
    t_A = c(rnorm(20, 3), rnorm(180)),
    t_B = rnorm(200)
  )
  module_genes <- list(m1 = genes[1:20], m2 = genes[50:80])
  out <- run_module_fgsea(rank_wide, module_genes, c("A", "B"))
  expect_setequal(unique(out$contrast), c("A", "B"))
  expect_equal(nrow(out), 4L)
  expect_true(all(is.finite(out$NES)))
  expect_true(all(out$padj >= 0 & out$padj <= 1))
  a_m1 <- out$NES[out$contrast == "A" & out$module == "m1"]
  expect_gt(a_m1, 0)
})

test_that("module_trait_cor recovers a planted positive correlation", {
  set.seed(2)
  eig <- matrix(rnorm(30 * 2), 30, 2, dimnames = list(NULL, c("blue", "brown")))
  traits <- cbind(age = eig[, "blue"] + rnorm(30, sd = 0.1), grip = rnorm(30))
  out <- module_trait_cor(eig, traits)
  expect_named(out, c("module", "trait", "r", "p", "n", "padj"))
  expect_equal(unique(out$n), 30L)
  r_blue_age <- out$r[out$module == "blue" & out$trait == "age"]
  expect_gt(r_blue_age, 0.9)
  p_blue_age <- out$p[out$module == "blue" & out$trait == "age"]
  expect_lt(p_blue_age, 0.001)
  expect_true(all(out$padj >= out$p - 1e-9))
})

test_that("module_trait_lmm sign matches the planted slope and respects pairing", {
  set.seed(3)
  subj <- rep(paste0("S", 1:15), each = 2)
  sid <- paste0(subj, "_", rep(c("T1", "T2"), 15))
  trait_val <- rnorm(30)
  eig_val <- 0.8 * trait_val + rnorm(30, sd = 0.3)
  eig_long <- data.frame(sample_id = sid, module = "blue", eigengene = eig_val)
  trait_long <- data.frame(sample_id = sid, trait = "age", value = trait_val)
  block <- setNames(subj, sid)
  out <- module_trait_lmm(eig_long, trait_long, block)
  expect_gt(out$beta[1], 0)
  expect_lt(out$p[1], 0.05)
})
