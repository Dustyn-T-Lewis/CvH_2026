# methods/12_msimpute_v2_mnar.R
# msImpute v2-mnar — hybrid: v2 for MAR, left-shift for MNAR
# Hediyeh-zadeh et al. 2023
# Hybrid; classifier-dependent

impute_msImpute_v2_mnar <- function(mat, meta, is_mnar, ...) {
  requireNamespace("msImpute", quietly = TRUE)
  # v2 handles the MAR component; MNAR proteins get left-shifted
  result <- msImpute::msImpute(mat, method = "v2-mnar", group = factor(meta$Group_Time))
  dimnames(result) <- dimnames(mat)
  result
}
