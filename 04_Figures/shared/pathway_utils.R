# Unified Pathway Enrichment Utility
# Post-hoc Jaccard deduplication (Reimand et al. 2019, Nat Protocols)
# Sources: MSigDB Hallmark (H), Canonical Pathways (C2:CP), GO:BP (C5:GO:BP)

#' Greedy Jaccard deduplication (single-pass, no stratification)
#'
#' Sorts results by padj ascending. For each term, checks Jaccard overlap with
#' all previously kept terms. If Jaccard > cutoff with any kept term, the term
#' is dropped (the more-significant one was already kept).
#'
#' @param results tibble with at least columns: pathway, padj
#' @param pathways named list of gene sets (character vectors)
#' @param jaccard_cutoff numeric, drop if Jaccard > this (default 0.5)
#' @return filtered tibble with redundant terms removed
deduplicate_enrichment_flat <- function(results, pathways, jaccard_cutoff = 0.5) {
  if (nrow(results) == 0) return(results)

  results <- results[order(results$padj), ]
  kept_names <- character(0)
  kept_sets  <- list()
  keep_mask  <- logical(nrow(results))

  for (i in seq_len(nrow(results))) {
    pw_name <- results$pathway[i]
    pw_genes <- pathways[[pw_name]]
    if (is.null(pw_genes)) { keep_mask[i] <- TRUE; next }

    is_redundant <- FALSE
    for (j in seq_along(kept_sets)) {
      inter <- length(intersect(pw_genes, kept_sets[[j]]))
      union <- length(union(pw_genes, kept_sets[[j]]))
      if (union > 0 && (inter / union) > jaccard_cutoff) {
        is_redundant <- TRUE
        break
      }
    }

    if (!is_redundant) {
      keep_mask[i] <- TRUE
      kept_names <- c(kept_names, pw_name)
      kept_sets[[length(kept_sets) + 1]] <- pw_genes
    }
  }

  results[keep_mask, ]
}

#' Database-stratified Jaccard deduplication of enrichment results
#'
#' Two-pass strategy: (1) dedup within each database to remove internal
#' redundancy (Reactome sub-pathways, nested GO terms), then (2) interleave
#' databases by rank so each gets fair representation. Without stratification,
#' databases with many granular sub-pathways (e.g. Reactome) dominate the
#' top results by winning every Jaccard comparison.
#'
#' Falls back to flat (unstratified) dedup if no 'database' column exists.
#'
#' @param results tibble with columns: pathway, padj, and optionally database
#' @param pathways named list of gene sets (character vectors)
#' @param jaccard_cutoff numeric, drop if Jaccard > this (default 0.5)
#' @return filtered tibble with redundant terms removed
deduplicate_enrichment <- function(results, pathways, jaccard_cutoff = 0.5) {
  if (nrow(results) == 0) return(results)

  # Fall back to flat dedup if no database column
  if (!"database" %in% names(results)) {
    return(deduplicate_enrichment_flat(results, pathways, jaccard_cutoff))
  }

  # Pass 1: dedup within each database independently
  dbs <- unique(results$database)
  within_dedup <- list()
  for (db in dbs) {
    db_rows <- results[results$database == db, ]
    within_dedup[[db]] <- deduplicate_enrichment_flat(db_rows, pathways,
                                                      jaccard_cutoff)
  }

  # Combine survivors and sort by padj (no cross-database dedup).
  # Cross-database overlap reflects complementary annotation of the same biology
  # from different curation perspectives, not true redundancy.
  survivors <- do.call(rbind, within_dedup)
  survivors[order(survivors$padj), ]
}


#' Build unified pathway collection from MSigDB
#'
#' Combines Hallmark (H), KEGG Medicus, Reactome, and GO:BP.
#' Optionally includes GO Slim gene sets. Filters disease/cancer terms from
#' curated pathways. Applies size filters.
#'
#' @param species character, default "Homo sapiens"
#' @param min_size integer, minimum gene set size (default 10)
#' @param max_size integer, maximum gene set size (default 500)
#' @param include_goslim logical, include GO Slim gene sets (default TRUE)
#' @return named list of character vectors (pathway name -> gene symbols)
build_pathway_collection <- function(species = "Homo sapiens",
                                     min_size = 10, max_size = 500,
                                     include_goslim = TRUE,
                                     exclude_variants = FALSE) {
  requireNamespace("msigdbr", quietly = TRUE)

  hallmark <- msigdbr::msigdbr(species = species, collection = "H")
  kegg     <- msigdbr::msigdbr(species = species, collection = "C2",
                                subcollection = "CP:KEGG_MEDICUS")
  reactome <- msigdbr::msigdbr(species = species, collection = "C2",
                                subcollection = "CP:REACTOME")
  gobp     <- msigdbr::msigdbr(species = species, collection = "C5",
                                subcollection = "GO:BP")

  disease_pat <- paste0("DISEASE|CANCER|TUMOR|CARCINOMA|LEUKEMIA|LYMPHOMA|",
                        "MELANOMA|GLIOMA|HEPATITIS|HIV|INFECTION|VIRAL|",
                        "BACTERIAL|PARASIT")
  kegg     <- kegg[!grepl(disease_pat, kegg$gs_name, ignore.case = TRUE), ]
  reactome <- reactome[!grepl(disease_pat, reactome$gs_name, ignore.case = TRUE), ]

  if (exclude_variants) {
    kegg <- kegg[!grepl("_VARIANT_", kegg$gs_name), ]
  }

  cols <- c("gs_name", "gene_symbol")
  sets_list <- list(hallmark[, cols], kegg[, cols], reactome[, cols], gobp[, cols])
  dbs <- c("H", "KEGG", "Reactome", "GO:BP")
  all_sets <- do.call(rbind, sets_list)

  pw_list <- split(all_sets$gene_symbol, all_sets$gs_name)
  pw_list <- lapply(pw_list, unique)

  if (include_goslim) {
    goslim_sets <- build_goslim_gene_sets(
      species = species, min_size = min_size, max_size = max_size
    )
    pw_list <- c(pw_list, goslim_sets)
    dbs <- c(dbs, "GO Slim")
  }

  sizes <- vapply(pw_list, length, integer(1))
  pw_list <- pw_list[sizes >= min_size & sizes <= max_size]

  message(sprintf("Pathway collection: %d sets (%s), size %d-%d",
                  length(pw_list), paste(dbs, collapse = " + "),
                  min_size, max_size))
  pw_list
}


#' Build GO Slim gene sets from GO.db hierarchy
build_goslim_gene_sets <- function(species = "Homo sapiens",
                                   min_size = 10, max_size = 500) {
  requireNamespace("GO.db", quietly = TRUE)
  requireNamespace("org.Hs.eg.db", quietly = TRUE)
  requireNamespace("AnnotationDbi", quietly = TRUE)

  bp_slim <- c(
    "GO:0000278", "GO:0000910", "GO:0002181", "GO:0002376", "GO:0003012",
    "GO:0003013", "GO:0003014", "GO:0003016", "GO:0005975", "GO:0006091",
    "GO:0006260", "GO:0006281", "GO:0006310", "GO:0006325", "GO:0006351",
    "GO:0006355", "GO:0006399", "GO:0006457", "GO:0006520", "GO:0006629",
    "GO:0006766", "GO:0006886", "GO:0006913", "GO:0006914", "GO:0006954",
    "GO:0007005", "GO:0007010", "GO:0007018", "GO:0007031", "GO:0007059",
    "GO:0007126", "GO:0007155", "GO:0007163", "GO:0007586", "GO:0009100",
    "GO:0012501", "GO:0016071", "GO:0016192", "GO:0023052", "GO:0030154",
    "GO:0030163", "GO:0030198", "GO:0032200", "GO:0034330", "GO:0042060",
    "GO:0042180", "GO:0042254", "GO:0044782", "GO:0048856", "GO:0048870",
    "GO:0050877", "GO:0051604", "GO:0055085", "GO:0055086", "GO:0061024",
    "GO:0065003", "GO:0071941", "GO:0072659", "GO:0098542", "GO:0098754",
    "GO:0140014", "GO:1901135"
  )

  offspring <- as.list(GO.db::GOBPOFFSPRING)

  suppressMessages({
    go_genes <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "GO"),
      keytype = "GO",
      columns = c("SYMBOL", "ONTOLOGY")
    )
  })
  go_bp_genes <- go_genes[!is.na(go_genes$ONTOLOGY) & go_genes$ONTOLOGY == "BP", ]
  go_to_symbols <- split(go_bp_genes$SYMBOL, go_bp_genes$GO)

  goslim_sets <- list()
  slim_names <- vapply(bp_slim, function(id) {
    tryCatch(AnnotationDbi::Term(GO.db::GOTERM[[id]]),
             error = function(e) NA_character_)
  }, character(1))

  for (i in seq_along(bp_slim)) {
    go_id <- bp_slim[i]
    go_term <- slim_names[i]
    if (is.na(go_term)) next

    all_terms <- go_id
    desc <- offspring[[go_id]]
    if (!is.null(desc)) all_terms <- c(all_terms, desc)

    genes <- unique(unlist(go_to_symbols[intersect(all_terms, names(go_to_symbols))],
                           use.names = FALSE))
    genes <- genes[!is.na(genes)]

    if (length(genes) >= min_size && length(genes) <= max_size) {
      set_name <- paste0("GOSLIM_", toupper(gsub(" ", "_", go_term)))
      goslim_sets[[set_name]] <- genes
    }
  }

  message(sprintf("GO Slim: %d/%d terms passed size filter (%d-%d)",
                  length(goslim_sets), length(bp_slim), min_size, max_size))
  goslim_sets
}


#' Run fGSEA on unified pathway collection with post-hoc deduplication
run_fgsea_deduplicated <- function(ranks, pathways, jaccard_cutoff = 0.5,
                                   nperm = 10000, min_size = 15,
                                   max_size = 500) {
  requireNamespace("fgsea", quietly = TRUE)

  res <- fgsea::fgseaMultilevel(
    pathways    = pathways,
    stats       = ranks,
    minSize     = min_size,
    maxSize     = max_size,
    nPermSimple = nperm,
    eps         = 0
  )
  res <- as.data.frame(res)
  res$database <- classify_database(res$pathway)
  res <- tibble::as_tibble(res)

  keep_cols <- c("pathway", "padj", "NES", "size", "leadingEdge",
                 "database", "pval", "ES", "log2err")
  res <- res[, intersect(keep_cols, names(res))]

  sig   <- res[!is.na(res$padj) & res$padj < 0.05, ]
  nonsig <- res[is.na(res$padj) | res$padj >= 0.05, ]

  sig_dedup <- deduplicate_enrichment(sig, pathways, jaccard_cutoff)

  n_removed <- nrow(sig) - nrow(sig_dedup)
  pct <- if (nrow(sig) > 0) round(100 * n_removed / nrow(sig), 1) else 0
  message(sprintf("fGSEA dedup: %d sig -> %d kept (removed %d, %.1f%%)",
                  nrow(sig), nrow(sig_dedup), n_removed, pct))

  rbind(sig_dedup, nonsig)
}

#' Run over-representation analysis with post-hoc deduplication
run_ora_deduplicated <- function(genes, universe, pathways,
                                 jaccard_cutoff = 0.5,
                                 min_size = 10, max_size = 500,
                                 padj_cutoff = 0.05) {
  requireNamespace("fgsea", quietly = TRUE)

  genes <- intersect(genes, universe)

  pw_by_db <- split(names(pathways), classify_database(names(pathways)))
  db_results <- list()
  for (db in names(pw_by_db)) {
    db_pw <- pathways[pw_by_db[[db]]]
    if (length(db_pw) < 2) next
    db_res <- fgsea::fora(
      pathways = db_pw,
      genes    = genes,
      universe = universe,
      minSize  = min_size,
      maxSize  = max_size
    )
    db_res <- as.data.frame(db_res)
    db_res$database <- db
    db_results[[db]] <- db_res
  }
  res <- do.call(rbind, db_results)

  N <- length(universe)
  K <- length(genes)
  res$odds_ratio <- vapply(seq_len(nrow(res)), function(i) {
    a <- res$overlap[i]
    b <- K - a
    c <- res$size[i] - a
    d <- N - K - c
    if (b == 0 || c == 0) Inf else (a * d) / (b * c)
  }, numeric(1))

  res <- tibble::as_tibble(res)

  sig <- res[!is.na(res$padj) & res$padj < padj_cutoff, ]
  sig_dedup <- deduplicate_enrichment(sig, pathways, jaccard_cutoff)

  n_removed <- nrow(sig) - nrow(sig_dedup)
  pct <- if (nrow(sig) > 0) round(100 * n_removed / nrow(sig), 1) else 0
  message(sprintf("ORA dedup: %d sig -> %d kept (removed %d, %.1f%%)",
                  nrow(sig), nrow(sig_dedup), n_removed, pct))

  sig_dedup
}


#' Classify pathway name to database source
classify_database <- function(pathway_names) {
  dplyr::case_when(
    grepl("^HALLMARK_",       pathway_names) ~ "Hallmark",
    grepl("^REACTOME_",       pathway_names) ~ "Reactome",
    grepl("^KEGG_MEDICUS_",   pathway_names) ~ "KEGG",
    grepl("^KEGG_",           pathway_names) ~ "KEGG",
    grepl("^GOSLIM_",         pathway_names) ~ "GO Slim",
    grepl("^GOBP_",           pathway_names) ~ "GO:BP",
    TRUE ~ "Other"
  )
}


# MSigDB pathway ID -> 15 consolidated categories (keyword rules)
CONSOLIDATED_PATHWAY_ORDER <- c(
  "Muscle & Contractile", "Cytoskeleton & Motility", "ECM & Adhesion",
  "Lipid Metabolism", "Carbohydrate & Energy Metabolism",
  "Amino Acid & Cofactor Metabolism",
  "Mitochondria & Energy", "Protein Homeostasis",
  "Transport", "Translation & Ribosome", "Transcription & Chromatin",
  "Immune & Inflammation", "DNA & Cell Cycle", "Circulatory System",
  "Development", "Other"
)

CONSOLIDATED_COLORS <- c(
  "Muscle & Contractile"              = "#E57373",
  "Cytoskeleton & Motility"           = "#FFB74D",
  "ECM & Adhesion"                    = "#FFF176",
  "Lipid Metabolism"                  = "#AED581",
  "Carbohydrate & Energy Metabolism"  = "#81C784",
  "Amino Acid & Cofactor Metabolism"  = "#66BB6A",
  "Mitochondria & Energy"             = "#4DB6AC",
  "Protein Homeostasis"               = "#4FC3F7",
  "Transport"                         = "#7986CB",
  "Translation & Ribosome"            = "#BA68C8",
  "Transcription & Chromatin"         = "#AB47BC",
  "Immune & Inflammation"             = "#A1887F",
  "DNA & Cell Cycle"                  = "#90A4AE",
  "Circulatory System"                = "#CE93D8",
  "Development"                       = "#B0BEC5",
  "Other"                             = "#D0D0D0"
)

