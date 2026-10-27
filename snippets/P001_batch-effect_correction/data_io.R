#!/usr/bin/env Rscript
# =============================================================================
# data_io.R -- read / write / convert tabular data
#              (csv, tsv, xlsx, parquet, rds)  ->  clean data.frame
# =============================================================================
# Use it two ways:
#   1. As a module : source("data_io.R"); df <- read_table_any("counts.xlsx")
#   2. As a CLI    : Rscript data_io.R -i counts.xlsx --out-format parquet
#
# Recommended hand-off format between R and Python: PARQUET
#   * typed columns (int / double / string / bool survive the round trip,
#     whereas csv/tsv make both sides re-guess types and xlsx may turn
#     gene names like SEPT2 into dates)
#   * columnar + compressed: ~5-10x smaller and much faster than csv/xlsx
#   * pandas.read_parquet() and arrow::read_parquet() read the same file
#   Use tsv when a human has to open it; xlsx only as an *input* format.
#   Caveat: parquet has no row names -> keep IDs in a normal column
#   (this module never relies on row names).
#
# Packages (only what you actually use is needed):
#   install.packages(c("optparse", "data.table", "readxl", "writexl"))
#   install.packages("arrow")        # or the lighter: install.packages("nanoparquet")
# =============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

msg <- function(fmt, ...) {
  txt <- if (...length() > 0L) sprintf(fmt, ...) else fmt
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), txt))
}

need_pkg <- function(pkg, hint = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required for this step but is not installed.\n  %s",
                 pkg, hint %||% sprintf("install.packages(\"%s\")", pkg)),
         call. = FALSE)
  }
  invisible(TRUE)
}

FORMATS <- c("csv", "tsv", "xlsx", "parquet", "rds")
.EXT_MAP <- c(csv = "csv", tsv = "tsv", tab = "tsv", txt = "tsv",
              xlsx = "xlsx", xlsm = "xlsx", xls = "xlsx",
              parquet = "parquet", pq = "parquet", rds = "rds")

# ---- small helpers ----------------------------------------------------------

# ".csv.gz" -> "csv" ; the sub() pattern is a REGEX (anchored with $), not a glob
detect_format <- function(path) {
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", path, ignore.case = TRUE)))
  fmt <- unname(.EXT_MAP[ext])
  if (length(fmt) != 1L || is.na(fmt)) {
    stop(sprintf("Cannot infer the format of '%s' from its extension. Known extensions: %s",
                 path, paste(names(.EXT_MAP), collapse = ", ")), call. = FALSE)
  }
  fmt
}

path_stem <- function(path) {
  tools::file_path_sans_ext(sub("\\.gz$", "", basename(path), ignore.case = TRUE))
}

# "2" -> 2L (sheet position) ; "Counts" -> "Counts" (sheet name). REGEX: all digits
parse_sheet <- function(x) {
  if (is.null(x)) return(1L)
  if (is.character(x) && grepl("^[0-9]+$", x)) as.integer(x) else x
}

split_csv <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

# numeric IDs (e.g. 100000) must not become "1e+05"
id_chr <- function(x) {
  if (is.numeric(x)) format(x, scientific = FALSE, trim = TRUE, digits = 15) else as.character(x)
}

parquet_backend <- function() {
  if (requireNamespace("arrow", quietly = TRUE) &&
      isTRUE(tryCatch(arrow::arrow_with_parquet(), error = function(e) TRUE))) return("arrow")
  if (requireNamespace("nanoparquet", quietly = TRUE)) return("nanoparquet")
  stop("Parquet support needs the 'arrow' or 'nanoparquet' package.\n",
       "  install.packages(\"arrow\")   # or the lighter: install.packages(\"nanoparquet\")",
       call. = FALSE)
}

# Make the column names safe: no empty names, no duplicates, no pandas index artefact
tidy_df <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  nm <- trimws(names(df))
  if ("__index_level_0__" %in% nm) {
    msg("Parquet carries a pandas index column; renamed '__index_level_0__' to 'index'.")
    nm[nm == "__index_level_0__"] <- "index"
  }
  empty <- is.na(nm) | nm == ""
  if (any(empty)) {
    warning(sprintf("%d column(s) had no name; named V<position>.", sum(empty)), call. = FALSE)
    nm[empty] <- paste0("V", which(empty))
  }
  if (anyDuplicated(nm)) {
    warning(sprintf("Duplicated column names (%s); made unique with a numeric suffix.",
                    paste(head(unique(nm[duplicated(nm)]), 5), collapse = ", ")), call. = FALSE)
    nm <- make.unique(nm, sep = "_")
  }
  names(df) <- nm
  df
}

# ---- reader -----------------------------------------------------------------

#' Read csv / tsv / xlsx / parquet / rds into a plain data.frame.
#' Column names are kept exactly as in the file (no make.names mangling).
read_table_any <- function(path, format = "auto", sheet = 1L, sep = NULL,
                           guess_max = 10000L,
                           na_strings = c("NA", "NaN", "N/A", "NULL")) {
  if (!file.exists(path)) stop(sprintf("File not found: %s", path), call. = FALSE)
  format <- if (identical(format, "auto")) detect_format(path) else match.arg(format, FORMATS)

  df <- switch(
    format,
    csv = , tsv = {
      need_pkg("data.table")
      # integer64 = "double": huge integer columns stay ordinary numerics
      data.table::fread(path, sep = sep %||% (if (format == "csv") "," else "\t"),
                        data.table = FALSE, check.names = FALSE,
                        na.strings = na_strings, integer64 = "double")
    },
    xlsx = {
      need_pkg("readxl")
      # readxl guesses column types from the first `guess_max` rows
      readxl::read_excel(path, sheet = sheet, col_names = TRUE, guess_max = guess_max,
                         .name_repair = "minimal", na = c("", na_strings))
    },
    parquet = {
      if (parquet_backend() == "arrow") arrow::read_parquet(path, as_data_frame = TRUE)
      else nanoparquet::read_parquet(path)
    },
    rds = {
      obj <- readRDS(path)
      if (is.data.frame(obj)) obj else as.data.frame(obj)
    }
  )
  tidy_df(df)
}

# ---- writer -----------------------------------------------------------------

#' Write a data.frame as csv / tsv / xlsx / parquet / rds.
#' row_names_col: if given, the row names are written as a first column of that name.
write_table_any <- function(df, path, format = "auto", sheet = "Sheet1",
                            compression = "snappy", row_names_col = NULL) {
  format <- if (identical(format, "auto")) detect_format(path) else match.arg(format, FORMATS)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!is.null(row_names_col)) {
    df <- data.frame(rownames(df), df, check.names = FALSE, stringsAsFactors = FALSE)
    names(df)[1] <- row_names_col
    rownames(df) <- NULL
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  switch(
    format,
    csv = , tsv = {
      need_pkg("data.table")
      # NB: do not pass na = "NA" -- in data.table 1.14 that makes fwrite quote EVERY string.
      # Missing values are therefore written as empty cells (fread reads them back as NA).
      data.table::fwrite(df, path, sep = if (format == "csv") "," else "\t")
    },
    xlsx = {
      need_pkg("writexl")
      if (nrow(df) + 1L > 1048576L || ncol(df) > 16384L) {
        stop("Table exceeds Excel's limits (1,048,576 rows x 16,384 columns); use parquet or tsv.",
             call. = FALSE)
      }
      writexl::write_xlsx(stats::setNames(list(df), sheet), path)
    },
    parquet = {
      if (parquet_backend() == "arrow") arrow::write_parquet(df, path, compression = compression)
      else nanoparquet::write_parquet(df, path, compression = compression)
    },
    rds = saveRDS(df, path)
  )
  invisible(path)
}

#' Read in one format, write in another.
convert_table <- function(input, output, in_format = "auto", out_format = "auto", ...) {
  df <- read_table_any(input, format = in_format, ...)
  write_table_any(df, output, format = out_format)
  invisible(df)
}

# ---- data.frame -> numeric matrix ------------------------------------------

#' Turn "ID column + numeric columns" into a numeric matrix with row names = IDs.
df_to_matrix <- function(df, id_col = 1L, data_cols = NULL) {
  id_idx <- if (is.character(id_col)) match(id_col, names(df)) else id_col
  if (is.na(id_idx)) stop(sprintf("ID column '%s' not found.", id_col), call. = FALSE)
  ids <- id_chr(df[[id_idx]])
  if (anyNA(ids) || any(!nzchar(trimws(ids)))) stop("Empty or missing values in the ID column.", call. = FALSE)
  if (anyDuplicated(ids)) {
    stop(sprintf("Duplicated IDs (e.g. %s).",
                 paste(head(unique(ids[duplicated(ids)]), 5), collapse = ", ")), call. = FALSE)
  }
  if (is.null(data_cols)) data_cols <- setdiff(seq_along(df), id_idx)
  sub <- df[, data_cols, drop = FALSE]
  bad <- names(sub)[!vapply(sub, is.numeric, logical(1))]
  if (length(bad)) {
    stop(sprintf("Non-numeric data column(s): %s%s. Check for text cells, decimal commas or stray header rows.",
                 paste(head(bad, 5), collapse = ", "), if (length(bad) > 5L) ", ..." else ""),
         call. = FALSE)
  }
  m <- as.matrix(sub)
  storage.mode(m) <- "double"
  rownames(m) <- ids
  m
}

# ---- CLI --------------------------------------------------------------------

convert_cli <- function(argv = commandArgs(trailingOnly = TRUE)) {
  need_pkg("optparse")
  opts <- list(
    optparse::make_option(c("-i", "--input"), dest = "input", type = "character",
      help = "Input file (.csv .tsv .txt .xlsx .parquet .rds; .gz allowed for csv/tsv)."),
    optparse::make_option(c("-o", "--output"), dest = "output", type = "character",
      help = "Output file. Default: <input stem>.<out-format> next to the input."),
    optparse::make_option("--in-format", dest = "in_format", type = "character", default = "auto",
      help = "csv | tsv | xlsx | parquet | rds | auto (from extension) [default %default]"),
    optparse::make_option("--out-format", dest = "out_format", type = "character", default = "auto",
      help = "csv | tsv | xlsx | parquet | rds | auto (from --output) [default %default]"),
    optparse::make_option("--sheet", dest = "sheet", type = "character", default = "1",
      help = "xlsx input: sheet position or name [default %default]"),
    optparse::make_option("--sep", dest = "sep", type = "character",
      help = "Override the field separator for csv/tsv input (e.g. ';')."),
    optparse::make_option("--out-sheet", dest = "out_sheet", type = "character", default = "Sheet1",
      help = "xlsx output: sheet name [default %default]"),
    optparse::make_option("--compression", dest = "compression", type = "character", default = "snappy",
      help = "parquet output: snappy | zstd | gzip | uncompressed [default %default]")
  )
  parser <- optparse::OptionParser(
    usage = "Rscript %prog -i INPUT (-o OUTPUT | --out-format FMT) [options]",
    option_list = opts,
    description = "Convert csv / tsv / xlsx / parquet / rds tables into one another (R and Python friendly).")
  opt <- optparse::parse_args(parser, args = argv)

  if (is.null(opt$input)) {
    optparse::print_help(parser)
    stop("--input is required.", call. = FALSE)
  }
  if (is.null(opt$output)) {
    if (identical(opt$out_format, "auto")) stop("Give --output or --out-format.", call. = FALSE)
    opt$output <- file.path(dirname(opt$input), paste0(path_stem(opt$input), ".", opt$out_format))
  }
  if (identical(normalizePath(opt$output, mustWork = FALSE), normalizePath(opt$input, mustWork = FALSE))) {
    stop("Refusing to overwrite the input file; choose a different --output.", call. = FALSE)
  }

  df <- read_table_any(opt$input, format = opt$in_format, sheet = parse_sheet(opt$sheet), sep = opt$sep)
  msg("Read %d rows x %d columns from %s", nrow(df), ncol(df), opt$input)
  cls <- vapply(df, function(x) class(x)[1], character(1))
  msg("First columns: %s", paste(sprintf("%s<%s>", head(names(df), 8), head(cls, 8)), collapse = ", "))

  write_table_any(df, opt$output, format = opt$out_format, sheet = opt$out_sheet,
                  compression = opt$compression)
  msg("Wrote %s", opt$output)

  fmt_out <- if (identical(opt$out_format, "auto")) detect_format(opt$output) else opt$out_format
  hint <- switch(fmt_out,
    parquet = sprintf("Python: pandas.read_parquet(\"%s\")   R: arrow::read_parquet(\"%s\")", opt$output, opt$output),
    csv     = sprintf("Python: pandas.read_csv(\"%s\")", opt$output),
    tsv     = sprintf("Python: pandas.read_csv(\"%s\", sep=\"\\t\")", opt$output),
    xlsx    = sprintf("Python: pandas.read_excel(\"%s\")  (needs openpyxl)", opt$output),
    rds     = "R only: readRDS()")
  msg("%s", hint)
  invisible(df)
}

# run the CLI only when executed as `Rscript data_io.R ...`, never when source()d
if (sys.nframe() == 0L && !interactive()) {
  tryCatch(convert_cli(), error = function(e) {
    message("ERROR: ", conditionMessage(e))
    quit(status = 1L)
  })
}