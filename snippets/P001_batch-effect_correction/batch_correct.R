#!/usr/bin/env Rscript
# =============================================================================
# batch_correct.R -- batch-effect removal + before/after PCA
# =============================================================================
# Methods (--method):
#   combat      sva::ComBat           log-scale / continuous data (empirical Bayes)
#   combat_seq  sva::ComBat_seq       RAW integer counts, returns counts
#   sva         sva::sva              surrogate variables; Batch label NOT used, needs --covariates
#   limma       limma::removeBatchEffect   log-scale data (for visualisation / clustering)
#   ruv         RUVSeq::RUVg          needs negative-control features (spike-ins e.g. ERCC,
#                                     or empirical control genes); Batch label NOT used
#
# Inputs
#   --matrix   csv / tsv / xlsx / parquet. Default orientation: features in rows,
#              column 1 = feature ID, sample columns named by the design's SampleID.
#              Any extra non-sample columns (Sequence, Length, Annotation, ...) are
#              treated as annotation and re-attached, untouched, to the output.
#   --design   same formats. Column 1 = SampleID, column 2 = Batch,
#              remaining columns = optional covariates (Group, Sex, Age, ...).
#
# Outputs (in --outdir)
#   <prefix>_<method>_corrected.<fmt>   corrected matrix (same orientation as input)
#   <prefix>_PCA_before_after.pdf/png   PC1 vs PC2, coloured by Batch
#   <prefix>_pca_scores.tsv             PC scores, both stages
#   <prefix>_pca_summary.tsv            % variance and Batch R^2 per PC, both stages
#   <prefix>_run_info.txt               arguments + sessionInfo()
#   (+ SVs / RUV factors / dropped feature list when applicable)
#
# Example (note the quoted shell variables -- always "${var}"):
#   MAT="counts.parquet"; DES="design.tsv"; OUT="results/bc"
#   Rscript batch_correct.R -m "${MAT}" -d "${DES}" -o "${OUT}" \
#       --method combat_seq --covariates Group --pca-shape Group
#   Rscript batch_correct.R -m "${MAT}" -d "${DES}" -o "${OUT}" \
#       --method ruv --control-regex '^ERCC-' --ruv-k 2
#
# Packages:
#   install.packages(c("optparse", "data.table", "ggplot2", "readxl", "writexl", "arrow"))
#   BiocManager::install(c("sva", "limma", "RUVSeq"))
#   optional: ggrepel (sample labels on the PCA)
# =============================================================================

# ---- locate and source data_io.R (same folder as this script) ---------------
script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", grep("^--file=", args, value = TRUE))   # REGEX, anchored
  if (length(f) == 1L) return(dirname(normalizePath(f)))
  for (i in rev(seq_len(sys.nframe()))) {                           # when source()d
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  getwd()
}
.data_io <- file.path(script_dir(), "data_io.R")
if (!file.exists(.data_io)) stop("data_io.R must sit next to batch_correct.R: ", .data_io)
source(.data_io)

METHODS <- c("combat", "combat_seq", "sva", "limma", "ruv")
METHOD_LABEL <- c(combat = "ComBat", combat_seq = "ComBat-seq", sva = "SVA",
                  limma = "limma removeBatchEffect", ruv = "RUVg")
# colour-blind-safe (Okabe-Ito based)
PALETTE <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#999999", "#F0E442")

normalize_method <- function(x) {
  key <- gsub("[-. ]", "_", tolower(trimws(x)))
  alias <- c(combatseq = "combat_seq", removebatcheffect = "limma", ruvseq = "ruv", ruvg = "ruv")
  if (key %in% names(alias)) key <- alias[[key]]
  if (!key %in% METHODS) {
    stop(sprintf("Unknown --method '%s'. Choose one of: %s", x, paste(METHODS, collapse = ", ")), call. = FALSE)
  }
  key
}

# =============================================================================
# 1. Input validation / preparation
# =============================================================================

check_ids <- function(ids, what) {
  if (anyNA(ids) || any(!nzchar(trimws(ids)))) stop(sprintf("Empty or missing values in %s.", what), call. = FALSE)
  if (anyDuplicated(ids)) {
    stop(sprintf("Duplicated values in %s (e.g. %s).", what,
                 paste(head(unique(ids[duplicated(ids)]), 5), collapse = ", ")), call. = FALSE)
  }
}

#' Design: col 1 = SampleID, col 2 = Batch (by POSITION), the rest = covariates.
validate_design <- function(des) {
  if (ncol(des) < 2L) stop("Design needs at least 2 columns: SampleID, Batch [, covariates...].", call. = FALSE)
  if (!identical(names(des)[1:2], c("SampleID", "Batch"))) {
    warning(sprintf("Design columns 1-2 are '%s','%s'; treating them as SampleID and Batch.",
                    names(des)[1], names(des)[2]), call. = FALSE)
    names(des)[1:2] <- c("SampleID", "Batch")
  }
  des$SampleID <- id_chr(des$SampleID)
  check_ids(des$SampleID, "design SampleID")
  b <- id_chr(des$Batch)
  if (anyNA(b) || any(!nzchar(trimws(b)))) stop("Missing Batch value(s) in the design.", call. = FALSE)
  des$Batch <- factor(b, levels = unique(b))              # numeric codes are treated as labels
  if (nlevels(des$Batch) < 2L) stop("Need at least 2 batches in the design.", call. = FALSE)
  des
}

#' Align the matrix with the design -> numeric features x samples matrix.
prepare_inputs <- function(mat_df, des, orientation = c("features_rows", "samples_rows")) {
  orientation <- match.arg(orientation)
  id_name <- names(mat_df)[1]
  ids <- id_chr(mat_df[[1]])
  annot <- NULL

  if (orientation == "features_rows") {
    check_ids(ids, sprintf("feature ID column '%s'", id_name))
    if (anyDuplicated(names(mat_df))) {
      stop(sprintf("Duplicated column names in the matrix (e.g. %s).",
                   paste(head(unique(names(mat_df)[duplicated(names(mat_df))]), 5), collapse = ", ")), call. = FALSE)
    }
    miss <- setdiff(des$SampleID, names(mat_df))
    if (length(miss)) {
      stop(sprintf("%d design SampleID(s) are not columns of the matrix (e.g. %s). Check spelling or --orientation.",
                   length(miss), paste(head(miss, 5), collapse = ", ")), call. = FALSE)
    }
    non_num <- des$SampleID[!vapply(mat_df[des$SampleID], is.numeric, logical(1))]
    if (length(non_num)) {
      stop(sprintf("Non-numeric sample column(s): %s. Check for text cells or decimal commas.",
                   paste(head(non_num, 5), collapse = ", ")), call. = FALSE)
    }
    X <- as.matrix(mat_df[, des$SampleID, drop = FALSE])
    storage.mode(X) <- "double"
    rownames(X) <- ids
    ann_cols <- setdiff(names(mat_df)[-1], des$SampleID)
    if (length(ann_cols)) {
      annot <- mat_df[, ann_cols, drop = FALSE]
      msg("Carrying %d non-sample column(s) through as annotation: %s%s", length(ann_cols),
          paste(head(ann_cols, 8), collapse = ", "), if (length(ann_cols) > 8L) ", ..." else "")
    }
  } else {                                                # samples in rows
    check_ids(ids, sprintf("sample ID column '%s'", id_name))
    miss <- setdiff(des$SampleID, ids)
    if (length(miss)) {
      stop(sprintf("%d design SampleID(s) are not rows of the matrix (e.g. %s).",
                   length(miss), paste(head(miss, 5), collapse = ", ")), call. = FALSE)
    }
    sub <- mat_df[match(des$SampleID, ids), -1, drop = FALSE]
    if (anyDuplicated(names(sub))) stop("Duplicated feature names in the matrix header.", call. = FALSE)
    non_num <- names(sub)[!vapply(sub, is.numeric, logical(1))]
    if (length(non_num)) {
      stop(sprintf("Non-numeric feature column(s): %s.", paste(head(non_num, 5), collapse = ", ")), call. = FALSE)
    }
    X <- t(as.matrix(sub))
    storage.mode(X) <- "double"
    colnames(X) <- des$SampleID
  }
  list(X = X, annot = annot, id_name = id_name, orientation = orientation, all_ids = rownames(X))
}

#' counts / linear / log. "auto" only recognises integer counts -- anything else must be declared.
resolve_scale <- function(X, scale) {
  if (any(!is.finite(X))) stop("Matrix contains NA/NaN/Inf. Impute or filter before batch correction.", call. = FALSE)
  if (!identical(scale, "auto")) return(scale)
  if (all(X >= 0) && all(abs(X - round(X)) < 1e-8)) return("counts")
  stop("Input is not integer-valued, so I cannot tell if it is log-scale or linear-scale.\n",
       "  --input-scale log     already log-transformed (log2 intensities, logCPM, VST, ...)\n",
       "  --input-scale linear  non-negative, not logged (CPM, TPM, peak areas, ...) -> log2(x+1) is applied",
       call. = FALSE)
}

#' Model matrix of the biological covariates to PROTECT (intercept always column 1).
build_model <- function(des, covariates, factor_covariates = character(0)) {
  if (length(covariates) == 0L) return(stats::model.matrix(~ 1, data = des))
  miss <- setdiff(covariates, names(des))
  if (length(miss)) stop(sprintf("Covariate(s) not in the design: %s", paste(miss, collapse = ", ")), call. = FALSE)
  if (any(covariates %in% c("SampleID", "Batch"))) stop("Batch/SampleID cannot be used as covariates.", call. = FALSE)
  cov_df <- des[, covariates, drop = FALSE]
  for (nm in names(cov_df)) {
    v <- cov_df[[nm]]
    if (anyNA(v)) stop(sprintf("Covariate '%s' has missing values.", nm), call. = FALSE)
    if (is.character(v) || is.logical(v) || nm %in% factor_covariates) v <- factor(v)
    if ((is.factor(v) && nlevels(v) < 2L) || (is.numeric(v) && length(unique(v)) < 2L)) {
      stop(sprintf("Covariate '%s' is constant; remove it from --covariates.", nm), call. = FALSE)
    }
    cov_df[[nm]] <- v
  }
  stats::model.matrix(~ ., data = cov_df)
}

#' A covariate perfectly nested in Batch cannot be separated from it.
check_confounding <- function(des, mod) {
  full <- cbind(stats::model.matrix(~ Batch, data = des),
                mod[, setdiff(colnames(mod), "(Intercept)"), drop = FALSE])
  if (qr(full)$rank < ncol(full)) {
    stop("Batch is confounded with the chosen covariate(s) (design matrix is rank-deficient).\n",
         "  Batch effect and biology cannot be separated; drop the covariate or fix the study design.",
         call. = FALSE)
  }
  invisible(TRUE)
}

within_batch_zero_var <- function(X, batch) {
  bad <- rep(FALSE, nrow(X))
  for (b in levels(batch)) {
    idx <- which(batch == b)
    if (length(idx) < 2L) next
    bad <- bad | (apply(X[, idx, drop = FALSE], 1, stats::var) == 0)
  }
  bad
}

find_controls <- function(feature_ids, ids_file = NULL, regex = NULL) {
  ctrl <- character(0)
  if (!is.null(ids_file)) {
    ids <- trimws(readLines(ids_file, warn = FALSE))
    ids <- ids[nzchar(ids)]
    n_miss <- sum(!ids %in% feature_ids)
    if (n_miss) warning(sprintf("%d control ID(s) in %s are not in the matrix.", n_miss, ids_file), call. = FALSE)
    ctrl <- union(ctrl, intersect(ids, feature_ids))
  }
  if (!is.null(regex)) ctrl <- union(ctrl, feature_ids[grepl(regex, feature_ids)])   # REGEX, not glob
  idx <- which(feature_ids %in% ctrl)
  if (length(idx) == 0L) {
    stop("method 'ruv' needs negative-control features: use --control-regex '^ERCC-' (a REGEX, not a glob)\n",
         "  and/or --control-ids file.txt (one feature ID per line). No control feature matched.", call. = FALSE)
  }
  if (length(idx) < 10L) warning(sprintf("Only %d control features; RUVg works best with many (>= 10-20).", length(idx)), call. = FALSE)
  idx
}

# =============================================================================
# 2. Correction (pure function: matrix in, matrix out)
# =============================================================================

#' Remove SV effects but keep the intercept/covariate effects in `mod`
#' (same regression as the sva vignette's cleaningP()).
clean_with_svs <- function(X, mod, sv) {
  design <- cbind(mod, sv)
  fit <- stats::lm.fit(design, t(X))
  sv_coef <- fit$coefficients[-seq_len(ncol(mod)), , drop = FALSE]
  X - t(sv %*% sv_coef)
}

correct_batch <- function(X, des, method, input_scale,
                          covariates = character(0), factor_covariates = character(0),
                          nonparametric = FALSE, mean_only = FALSE, ref_batch = NULL,
                          n_sv = NULL, ruv_k = 1L, control_idx = NULL, seed = 1L) {
  method <- match.arg(method, METHODS)
  stopifnot(identical(colnames(X), des$SampleID))
  batch <- des$Batch
  bsize <- table(batch)
  mod <- build_model(des, covariates, factor_covariates)
  check_confounding(des, mod)
  log2p <- function(M) log2(M + 1)
  out <- list(X = NULL, scale = NULL, dropped = character(0), W = NULL, sv = NULL)

  # ---- count-scale methods ---------------------------------------------------
  if (method == "combat_seq") {
    need_pkg("sva", "BiocManager::install(\"sva\")")
    if (input_scale != "counts") {
      stop("ComBat-seq models raw integer counts; the input is not (--input-scale ", input_scale, ").", call. = FALSE)
    }
    if (any(bsize < 2L)) stop("ComBat-seq needs >= 2 samples in every batch.", call. = FALSE)
    cov_mat <- if (ncol(mod) > 1L) mod[, -1, drop = FALSE] else NULL
    out$X <- sva::ComBat_seq(counts = X, batch = batch, group = NULL, covar_mod = cov_mat)
    dimnames(out$X) <- dimnames(X)
    out$scale <- "counts"
    return(out)
  }

  if (method == "ruv") {
    need_pkg("RUVSeq", "BiocManager::install(\"RUVSeq\")")
    if (length(control_idx) == 0L) stop("method 'ruv' needs control features.", call. = FALSE)
    is_log <- input_scale != "counts"
    Xin <- if (input_scale == "linear") log2p(X) else X
    res <- RUVSeq::RUVg(Xin, cIdx = control_idx, k = ruv_k, isLog = is_log)
    out$X <- res$normalizedCounts
    dimnames(out$X) <- dimnames(X)
    out$W <- res$W
    rownames(out$W) <- colnames(X)
    out$scale <- if (is_log) "log" else "counts"
    return(out)
  }

  # ---- log-scale methods -----------------------------------------------------
  Xl <- if (input_scale == "log") X else log2p(X)

  if (method == "limma") {
    need_pkg("limma", "BiocManager::install(\"limma\")")
    out$X <- limma::removeBatchEffect(Xl, batch = batch, design = mod)
  }

  if (method == "combat") {
    need_pkg("sva", "BiocManager::install(\"sva\")")
    if (!is.null(ref_batch) && !ref_batch %in% levels(batch)) {
      stop(sprintf("--ref-batch '%s' is not one of the batches: %s", ref_batch,
                   paste(levels(batch), collapse = ", ")), call. = FALSE)
    }
    if (any(bsize < 2L) && !mean_only) {
      warning("A batch has a single sample: switching to mean.only = TRUE.", call. = FALSE)
      mean_only <- TRUE
    }
    bad <- within_batch_zero_var(Xl, batch)       # ComBat stops on zero-variance rows
    if (any(bad)) {
      msg("ComBat: dropping %d feature(s) with zero variance inside a batch.", sum(bad))
      out$dropped <- rownames(Xl)[bad]
      Xl <- Xl[!bad, , drop = FALSE]
    }
    out$X <- sva::ComBat(dat = Xl, batch = batch, mod = mod, par.prior = !nonparametric,
                         mean.only = mean_only, ref.batch = ref_batch)
  }

  if (method == "sva") {
    need_pkg("sva", "BiocManager::install(\"sva\")")
    if (ncol(mod) < 2L) {
      stop("method 'sva' needs at least one biological covariate to protect (--covariates Group).\n",
           "  Without one, the surrogate variables would absorb the real signal as well.", call. = FALSE)
    }
    bad <- apply(Xl, 1, stats::var) == 0
    if (any(bad)) {
      msg("SVA: dropping %d constant feature(s).", sum(bad))
      out$dropped <- rownames(Xl)[bad]
      Xl <- Xl[!bad, , drop = FALSE]
    }
    mod0 <- stats::model.matrix(~ 1, data = des)
    n <- n_sv %||% sva::num.sv(Xl, mod, method = "be", seed = seed)
    msg("SVA: using %d surrogate variable(s).", as.integer(n))
    if (n < 1L) {
      warning("0 surrogate variables estimated; data returned unchanged.", call. = FALSE)
      out$X <- Xl
    } else {
      sv <- sva::sva(Xl, mod, mod0, n.sv = n)$sv
      dimnames(sv) <- list(colnames(Xl), paste0("SV", seq_len(ncol(sv))))
      out$X <- clean_with_svs(Xl, mod, sv)
      out$sv <- sv
    }
  }

  out$scale <- "log"
  out
}

# =============================================================================
# 3. PCA + batch metrics + plot
# =============================================================================

to_pca_scale <- function(M, scale) if (identical(scale, "log")) M else log2(M + 1)

run_pca <- function(M, top_n = 0L, scale_unit = FALSE, n_pc = 10L) {
  v <- apply(M, 1, stats::var)
  ok <- is.finite(v) & v > 0
  if (sum(ok) < 3L) stop("Fewer than 3 variable features available for PCA.", call. = FALSE)
  M <- M[ok, , drop = FALSE]
  v <- v[ok]
  if (top_n > 0L && nrow(M) > top_n) M <- M[order(v, decreasing = TRUE)[seq_len(top_n)], , drop = FALSE]
  pc <- stats::prcomp(t(M), center = TRUE, scale. = scale_unit)
  k <- min(n_pc, ncol(pc$x))
  scores <- data.frame(SampleID = rownames(pc$x), pc$x[, seq_len(k), drop = FALSE],
                       check.names = FALSE, stringsAsFactors = FALSE)
  rownames(scores) <- NULL
  list(scores = scores, pve = pc$sdev^2 / sum(pc$sdev^2), k = k, n_features = nrow(M))
}

#' R^2 of PC_i ~ Batch: how much of each PC is explained by batch (0 = none, 1 = all).
batch_r2 <- function(scores, des, k) {
  b <- des$Batch[match(scores$SampleID, des$SampleID)]
  vapply(seq_len(k), function(i) {
    y <- scores[[paste0("PC", i)]]
    summary(stats::lm(y ~ b))$r.squared
  }, numeric(1))
}

plot_pca_before_after <- function(pcs, des, shape_by = NULL, labels = FALSE, ellipse = FALSE, title = NULL) {
  need_pkg("ggplot2")
  df <- do.call(rbind, lapply(names(pcs), function(st) {
    p <- pcs[[st]]
    data.frame(Stage = sprintf("%s\nPC1 %.1f%% | PC2 %.1f%%", st, 100 * p$pve[1], 100 * p$pve[2]),
               SampleID = p$scores$SampleID, PC1 = p$scores$PC1, PC2 = p$scores$PC2,
               stringsAsFactors = FALSE)
  }))
  df$Stage <- factor(df$Stage, levels = unique(df$Stage))
  df$Batch <- des$Batch[match(df$SampleID, des$SampleID)]
  if (!is.null(shape_by)) df$Shape <- factor(des[[shape_by]][match(df$SampleID, des$SampleID)])

  g <- ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, colour = Batch)) +
    ggplot2::geom_point(size = 2.6, alpha = 0.9) +
    ggplot2::facet_wrap(~ Stage, scales = "free") +
    ggplot2::labs(x = "PC1", y = "PC2", colour = "Batch", title = title) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey95"))

  nb <- nlevels(df$Batch)
  if (nb <= length(PALETTE)) g <- g + ggplot2::scale_colour_manual(values = PALETTE[seq_len(nb)])
  if (!is.null(shape_by)) g <- g + ggplot2::aes(shape = Shape) + ggplot2::labs(shape = shape_by)
  if (ellipse) g <- g + ggplot2::stat_ellipse(ggplot2::aes(group = Batch), linetype = 2, na.rm = TRUE)
  if (labels) {
    lab <- if (requireNamespace("ggrepel", quietly = TRUE)) {
      ggrepel::geom_text_repel(ggplot2::aes(label = SampleID), size = 2.6, show.legend = FALSE, max.overlaps = Inf)
    } else {
      ggplot2::geom_text(ggplot2::aes(label = SampleID), size = 2.6, vjust = -0.8,
                         check_overlap = TRUE, show.legend = FALSE)
    }
    g <- g + lab
  }
  g
}

# =============================================================================
# 4. Output assembly
# =============================================================================

assemble_output <- function(X, inp) {
  if (inp$orientation == "features_rows") {
    out <- data.frame(rownames(X), stringsAsFactors = FALSE)
    names(out) <- inp$id_name
    if (!is.null(inp$annot)) {
      ann <- inp$annot[match(rownames(X), inp$all_ids), , drop = FALSE]
      rownames(ann) <- NULL
      out <- data.frame(out, ann, check.names = FALSE, stringsAsFactors = FALSE)
    }
    out <- data.frame(out, as.data.frame(X), check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    out <- data.frame(colnames(X), t(X), check.names = FALSE, stringsAsFactors = FALSE)
    names(out)[1] <- inp$id_name
  }
  rownames(out) <- NULL
  out
}

# =============================================================================
# 5. Pipeline + CLI
# =============================================================================

run_batch_correction <- function(opt) {
  method <- normalize_method(opt$method)
  set.seed(opt$seed)
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  prefix <- opt$prefix %||% path_stem(opt$matrix)
  covariates <- split_csv(opt$covariates)
  stem <- file.path(opt$outdir, prefix)

  msg("Reading design: %s", opt$design)
  des <- validate_design(read_table_any(opt$design, sheet = parse_sheet(opt$design_sheet)))
  if (!is.null(opt$pca_shape) && !opt$pca_shape %in% names(des)) {
    stop(sprintf("--pca-shape '%s' is not a design column.", opt$pca_shape), call. = FALSE)
  }
  msg("Reading matrix: %s", opt$matrix)
  inp <- prepare_inputs(read_table_any(opt$matrix, sheet = parse_sheet(opt$matrix_sheet)), des, opt$orientation)

  bt <- table(des$Batch)
  msg("%d features x %d samples; batches: %s", nrow(inp$X), ncol(inp$X),
      paste(sprintf("%s=%d", names(bt), as.vector(bt)), collapse = ", "))

  scale <- resolve_scale(inp$X, opt$input_scale)
  msg("Input scale: %s | method: %s", scale, METHOD_LABEL[[method]])
  if (method %in% c("combat_seq") && scale == "counts" && any(abs(inp$X - round(inp$X)) > 1e-8)) {
    stop("ComBat-seq needs integer counts.", call. = FALSE)
  }

  ctrl_idx <- NULL
  if (method == "ruv") {
    ctrl_idx <- find_controls(rownames(inp$X), opt$control_ids, opt$control_regex)
    msg("RUVg: %d control feature(s), k = %d", length(ctrl_idx), opt$ruv_k)
  }

  res <- correct_batch(inp$X, des, method, scale,
                       covariates = covariates, factor_covariates = split_csv(opt$factor_covariates),
                       nonparametric = opt$nonparametric, mean_only = opt$mean_only,
                       ref_batch = opt$ref_batch, n_sv = opt$n_sv, ruv_k = opt$ruv_k,
                       control_idx = ctrl_idx, seed = opt$seed)
  msg("Correction done: %d features kept, output scale = %s", nrow(res$X), res$scale)

  # ---- PCA before / after on the SAME features, both on the log2 scale -------
  X_before <- inp$X[rownames(res$X), , drop = FALSE]
  pb <- run_pca(to_pca_scale(X_before, scale), opt$pca_top, opt$pca_scale)
  pa <- run_pca(to_pca_scale(res$X, res$scale), opt$pca_top, opt$pca_scale)
  k <- min(pb$k, pa$k)

  summ <- function(stage, p) {
    data.frame(Stage = stage, PC = paste0("PC", seq_len(k)), var_explained = p$pve[seq_len(k)],
               batch_R2 = batch_r2(p$scores, des, k), stringsAsFactors = FALSE)
  }
  s_before <- summ("before", pb)
  s_after <- summ("after", pa)
  wr2 <- function(s) sum(s$var_explained * s$batch_R2, na.rm = TRUE) / sum(s$var_explained[!is.na(s$batch_R2)])
  msg("Batch R2 over the top %d PCs (variance-weighted): before %.3f -> after %.3f", k, wr2(s_before), wr2(s_after))

  as_long <- function(stage, p) {
    s <- p$scores[, c("SampleID", paste0("PC", seq_len(k))), drop = FALSE]
    data.frame(Stage = stage, Batch = as.character(des$Batch[match(s$SampleID, des$SampleID)]), s,
               check.names = FALSE, stringsAsFactors = FALSE)
  }
  write_table_any(rbind(as_long("before", pb), as_long("after", pa)), paste0(stem, "_pca_scores.tsv"))
  write_table_any(rbind(s_before, s_after), paste0(stem, "_pca_summary.tsv"))

  # ---- plot -------------------------------------------------------------------
  stages <- stats::setNames(list(pb, pa), c("Before correction", sprintf("After correction (%s)", METHOD_LABEL[[method]])))
  g <- plot_pca_before_after(stages, des, shape_by = opt$pca_shape, labels = opt$pca_labels,
                             ellipse = opt$ellipse,
                             title = sprintf("Batch effect before vs after %s", METHOD_LABEL[[method]]))
  for (ext in split_csv(opt$plot_format)) {
    f <- paste0(stem, "_PCA_before_after.", ext)
    tryCatch(suppressWarnings(ggplot2::ggsave(f, g, width = 10.5, height = 4.9, dpi = 300, bg = "white")),
             error = function(e) warning(sprintf("Could not write %s: %s", f, conditionMessage(e)), call. = FALSE))
  }

  # ---- corrected matrix -------------------------------------------------------
  X_out <- res$X
  out_scale <- res$scale
  if (isTRUE(opt$backtransform) && res$scale == "log" && scale != "log") {
    X_out <- pmax(2^X_out - 1, 0)
    out_scale <- "linear (back-transformed 2^x - 1, floored at 0)"
  }
  out_file <- sprintf("%s_%s_corrected.%s", stem, method, opt$out_format)
  write_table_any(assemble_output(X_out, inp), out_file)

  if (length(res$dropped)) writeLines(res$dropped, paste0(stem, "_dropped_features.txt"))
  if (!is.null(res$sv)) write_table_any(data.frame(SampleID = rownames(res$sv), res$sv, check.names = FALSE), paste0(stem, "_sva_surrogate_variables.tsv"))
  if (!is.null(res$W)) write_table_any(data.frame(SampleID = rownames(res$W), res$W, check.names = FALSE), paste0(stem, "_ruv_factors.tsv"))

  args_txt <- vapply(names(opt), function(n) sprintf("  --%s = %s", n, paste(format(opt[[n]]), collapse = ",")), character(1))
  writeLines(c("# batch_correct.R run info", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "",
               sprintf("method: %s | input scale: %s | output scale: %s", method, scale, out_scale),
               sprintf("batch R2 (variance-weighted, top %d PCs): before %.4f, after %.4f", k, wr2(s_before), wr2(s_after)),
               "", "## arguments", args_txt, "", "## session", utils::capture.output(print(utils::sessionInfo()))),
             paste0(stem, "_run_info.txt"))

  msg("Done. Outputs in %s", opt$outdir)
  invisible(list(corrected = X_out, pca_before = pb, pca_after = pa, plot = g))
}

parse_cli <- function(argv) {
  need_pkg("optparse")
  O <- function(flags, dest, type = "character", default = NULL, help = "", flag = FALSE) {
    if (flag) {
      optparse::make_option(flags, dest = dest, action = "store_true", default = FALSE, help = help)
    } else {
      optparse::make_option(flags, dest = dest, type = type, default = default, help = help)
    }
  }
  opts <- list(
    O(c("-m", "--matrix"), "matrix", help = "Data matrix (csv/tsv/xlsx/parquet). Default: features in rows, column 1 = feature ID."),
    O(c("-d", "--design"), "design", help = "Design table: col 1 = SampleID, col 2 = Batch, other columns = covariates."),
    O(c("-o", "--outdir"), "outdir", default = "batch_correction_out", help = "Output directory [default %default]"),
    O("--prefix", "prefix", help = "Output file prefix [default: matrix file stem]"),
    O("--method", "method", default = "combat", help = "combat | combat_seq | sva | limma | ruv [default %default]"),
    O("--orientation", "orientation", default = "features_rows", help = "features_rows | samples_rows [default %default]"),
    O("--input-scale", "input_scale", default = "auto", help = "auto (integer counts only) | counts | linear | log [default %default]"),
    O("--covariates", "covariates", help = "Design columns whose biology to PRESERVE, comma-separated (e.g. Group,Sex)."),
    O("--factor-covariates", "factor_covariates", help = "Numeric-coded covariates to force to categorical (e.g. 1/2/3 coded groups)."),
    O("--nonparametric", "nonparametric", flag = TRUE, help = "ComBat: non-parametric prior (slow)."),
    O("--mean-only", "mean_only", flag = TRUE, help = "ComBat: adjust batch means only (use with tiny batches)."),
    O("--ref-batch", "ref_batch", help = "ComBat: batch to use as reference (kept unchanged)."),
    O("--n-sv", "n_sv", type = "integer", help = "SVA: number of surrogate variables [default: estimated]"),
    O("--ruv-k", "ruv_k", type = "integer", default = 1L, help = "RUVg: number of unwanted-variation factors [default %default]"),
    O("--control-ids", "control_ids", help = "RUVg: text file with one control feature ID per line."),
    O("--control-regex", "control_regex", help = "RUVg: REGEX (not a glob) selecting control features, e.g. '^ERCC-'."),
    O("--backtransform", "backtransform", flag = TRUE, help = "Write 2^x - 1 instead of log2(x+1) values (counts/linear input only)."),
    O("--out-format", "out_format", default = "tsv", help = "csv | tsv | xlsx | parquet | rds [default %default]"),
    O("--pca-top", "pca_top", type = "integer", default = 0L, help = "Use only the N most variable features for PCA; 0 = all [default %default]"),
    O("--pca-scale", "pca_scale", flag = TRUE, help = "Unit-variance scaling of features before PCA."),
    O("--pca-shape", "pca_shape", help = "Design column mapped to point shape (e.g. Group)."),
    O("--pca-labels", "pca_labels", flag = TRUE, help = "Label points with SampleID."),
    O("--ellipse", "ellipse", flag = TRUE, help = "Draw a per-batch ellipse (needs >= 4 samples per batch)."),
    O("--plot-format", "plot_format", default = "pdf,png", help = "Comma-separated plot formats [default %default]"),
    O("--matrix-sheet", "matrix_sheet", default = "1", help = "xlsx matrix: sheet position or name [default %default]"),
    O("--design-sheet", "design_sheet", default = "1", help = "xlsx design: sheet position or name [default %default]"),
    O("--seed", "seed", type = "integer", default = 1L, help = "Random seed [default %default]")
  )
  parser <- optparse::OptionParser(
    usage = "Rscript %prog -m MATRIX -d DESIGN -o OUTDIR --method combat|combat_seq|sva|limma|ruv [options]",
    option_list = opts,
    description = "Batch-effect removal with before/after PCA. Requires data_io.R in the same folder.")
  opt <- optparse::parse_args(parser, args = argv)
  if (is.null(opt$matrix) || is.null(opt$design)) {
    optparse::print_help(parser)
    stop("--matrix and --design are required.", call. = FALSE)
  }
  opt
}

batch_correct_cli <- function(argv = commandArgs(trailingOnly = TRUE)) {
  tryCatch(run_batch_correction(parse_cli(argv)), error = function(e) {
    message("ERROR: ", conditionMessage(e))
    quit(status = 1L)
  })
}

# run the CLI only when executed as `Rscript batch_correct.R ...`, never when source()d
if (sys.nframe() == 0L && !interactive()) batch_correct_cli()