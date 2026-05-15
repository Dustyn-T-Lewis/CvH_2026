# Sample-ID normalization for F02.
# Sam uses zero-padded 3-digit subject IDs (CR006_T1); our pipeline uses
# CR6_T1. Always canonicalize to the short form before joining.

normalize_sample_id <- function(x) {
  x <- sub("^CR0*([0-9]+)_T([12])$", "CR\\1_T\\2", x)
  x <- sub("^PPS0*([0-9]+)$", "PPS\\1", x)
  x
}

# Quick self-check (only runs if sourced as main).
if (sys.nframe() == 0) {
  stopifnot(
    normalize_sample_id("CR006_T1") == "CR6_T1",
    normalize_sample_id("CR6_T1")   == "CR6_T1",
    normalize_sample_id("CR10_T2")  == "CR10_T2",
    normalize_sample_id("CR017_T2") == "CR17_T2",
    normalize_sample_id("PPS02")    == "PPS2",
    normalize_sample_id("PPS2")     == "PPS2",
    normalize_sample_id("PPS44")    == "PPS44"
  )
  message("_shared_keys.R: self-check passed")
}
