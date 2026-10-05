#!/usr/bin/env Rscript
"""
PredictAP — rankcal function
=============================
Compute ECDF-based percentile ranks of PrediXcan predicted gene expression
against a reference population (e.g. NFE).

Standalone R script (requires tidyverse, data.table, tidyselect)

Usage
-----
source("PredictAP_rankcal.R")

# Load reference data (x1, x2, ..., x49, result_table)
load_reference_data(tissue_results_dir = "dosage/results/",
                    result_table_path = "result_table.csv")

# Run rankcal
data <- fread("example_data.csv")
result <- rankcal(data, tissue = "x1")
head(result)
fwrite(result, "rank_result.csv")
"""

library(tidyverse)
library(data.table)
library(tidyselect)

# ── Global cache for reference data ────────────────────────────────────────
.predictap_env <- new.env()


# ── Core functions ─────────────────────────────────────────────────────────

#' ECDF-based percentile rank
#'
#' For each value in query, compute the fraction of reference values
#' that are less than or equal to it (empirical CDF at query).
#'
#' @param reference  Reference population values (numeric vector)
#' @param query      Query subject values (numeric vector)
#'
#' @return Numeric vector of percentile ranks in [0, 1]
#' @export
ecdf_fn <- function(reference, query) {
  ref <- sort(as.numeric(reference))
  q   <- as.numeric(query)
  rank <- sapply(q, function(x) sum(ref <= x) / length(ref))
  return(rank)
}


#' Load reference data for rankcal
#'
#' Loads x1...x49 tissue reference data and result_table from disk.
#' Cached in package environment to avoid repeated I/O.
#'
#' @param tissue_results_dir  Directory containing x1.txt ... x49.txt
#' @param result_table_path   Path to result_table.csv
#'
#' @return NULL (invisibly). Data loaded into package environment.
#' @export
load_reference_data <- function(tissue_results_dir, result_table_path) {
  message("[1] Loading result_table...")
  result_table <- fread(result_table_path)
  assign("result_table", result_table, envir = .predictap_env)

  message("[2] Loading tissue reference data (x1...x49)...")
  tissue_files <- list.files(tissue_results_dir, pattern = "^x\\d+\\.txt$", full.names = TRUE)

  if (length(tissue_files) == 0) {
    stop(sprintf("No tissue files (x1.txt ... x49.txt) found in %s", tissue_results_dir))
  }

  tissue_files <- sort(tissue_files)

  for (i in seq_along(tissue_files)) {
    tissue_data <- fread(tissue_files[i])
    tissue_name <- paste0("x", i)
    assign(tissue_name, tissue_data, envir = .predictap_env)
    if (i %% 10 == 0) message(sprintf("   Loaded %d / %d tissues", i, length(tissue_files)))
  }

  message(sprintf("   Loaded %d / %d tissues", length(tissue_files), length(tissue_files)))
  message(sprintf("   result_table: %d rows × %d cols", nrow(result_table), ncol(result_table)))
  invisible(NULL)
}


#' Compute ECDF percentile ranks for predicted gene expression
#'
#' Compares each subject's PrediXcan-predicted expression against the
#' reference population distribution for the specified tissue.
#'
#' @param data      PrediXcan output table (FID | IID | ENSG... | ENSG... | ...)
#' @param tissue    Tissue index, e.g. "x1" (Adipose_Subcutaneous) through "x49" (Whole_Blood)
#'
#' @return data.table with columns:
#'   FID, IID, gene, expression, percentile.rank,
#'   N.snps, mean, sd, ks.test.pval, European.Asian.mean.difference
#'
#' @export
rankcal <- function(data, tissue) {
  data <- as.data.table(data)

  # Get reference data from cache
  ref_data <- get(tissue, envir = .predictap_env)
  result_table <- get("result_table", envir = .predictap_env)

  nsubject <- nrow(data)
  nc <- ncol(ref_data)

  # ── Compute percentile rank for each gene ──────────────────────────────
  rankm <- matrix(NA, nsubject, nc, byrow = FALSE)

  for (i in 3:nc) {
    rankm[, i] <- ecdf_fn(ref_data[[i]], data[[i]])
  }

  rankm[, 1] <- data[[1]]  # FID
  rankm[, 2] <- data[[2]]  # IID
  rankm <- data.table(rankm)
  colnames(rankm) <- colnames(data)

  # ── Reshape to long format ─────────────────────────────────────────────
  l_data <- pivot_longer(data,
    cols = starts_with("ENSG"),
    names_to = "gene",
    values_to = "expression"
  )
  l_rank <- pivot_longer(rankm,
    cols = starts_with("ENSG"),
    names_to = "gene",
    values_to = "percentile.rank"
  )
  l_mer <- merge(l_data, l_rank, by = c("FID", "IID", "gene"))

  # ── Merge with result_table (tissue metadata) ──────────────────────────
  res_por <- result_table %>%
    select(gene, all_of(tissue), starts_with(paste0(tissue, "."))) %>%
    filter(!is.na(!!sym(tissue))) %>%
    merge(l_mer, by = "gene")

  # ── Final result with renamed columns ───────────────────────────────────
  result <- res_por %>%
    select(FID, IID, gene, expression, percentile.rank, starts_with(tissue)) %>%
    arrange(FID, IID)

  colnames(result)[c(6, 7, 8, 9, 10)] <-
    c("N.snps", "mean", "sd", "ks.test.pval", "European.Asian.mean.difference")

  return(result)
}


# ── Example usage ──────────────────────────────────────────────────────────
if (interactive() || !is.null(commandArgs(trailingOnly = FALSE))) {
  # Run example if sourced interactively

  cat("\n=== PredictAP Example Usage ===\n")
  cat("1. Load reference data:\n")
  cat('   load_reference_data(tissue_results_dir = "dosage/results/",\n')
  cat('                       result_table_path = "result_table.csv")\n\n')

  cat("2. Run rankcal:\n")
  cat('   data <- fread("example_data.csv")\n')
  cat('   result <- rankcal(data, tissue = "x1")\n')
  cat('   head(result)\n\n')

  cat("3. Save results:\n")
  cat('   fwrite(result, "rank_result.csv")\n\n')
}
