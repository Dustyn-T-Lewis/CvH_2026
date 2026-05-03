# methods/11_msimpute_v2.R
# msImpute v2 — low-rank approximation via barycenter
# Hediyeh-zadeh et al. 2023
# Pure MAR (low-rank)

impute_msImpute_v2 <- function(mat, meta, is_mnar, ...) {
  requireNamespace("msImpute", quietly = TRUE)
  result <- msImpute::msImpute(mat, method = "v2", group = factor(meta$Group_Time))
  dimnames(result) <- dimnames(mat)
  result
}
