#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(dplyr)
})

# ------------------------------ #
# simple arg parser: --key=value #
# ------------------------------ #
parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a) == 0) return(list())
  kv <- strsplit(a, "=", fixed = TRUE)
  # handle args without '=' gracefully
  keys   <- vapply(kv, function(z) if (length(z) >= 1) z[1] else NA_character_, character(1))
  values <- vapply(kv, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
  args <- setNames(values, keys)
  # defaults
  if (is.na(args[["--tolerance"]])) args[["--tolerance"]] <- "50"
  if (is.na(args[["--outdir"]]))    args[["--outdir"]]    <- "."
  args
}

args <- parse_args()
annotsv_path <- args[["--annotsv"]]           # required
common_path  <- args[["--common"]]            # required (5 or 6 columns allowed)
outdir       <- args[["--outdir"]] %||% "."
tolerance    <- as.integer(args[["--tolerance"]] %||% "50")

`%||%` <- function(a,b) if (is.null(a) || is.na(a)) b else a

if (any(is.na(c(annotsv_path, common_path)))) {
  stop("Usage: find_overlapping_common_SV.R --annotsv=sample.annotsv.tsv --common=common.bed [--outdir=dir] [--tolerance=N]")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# prefix from sample filename
sample_prefix <- {
  nm <- basename(annotsv_path)
  sub("\\.[^.]*$", "", nm)
}

# ------------------------------ #
# helpers                        #
# ------------------------------ #
canon_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x <- toupper(x)
  x[x == "X"]  <- "23"
  x[x == "Y"]  <- "24"
  x[x %in% c("M","MT")] <- "25"
  suppressWarnings(as.integer(x))
}

to_chr_prefix <- function(i) {
  v <- ifelse(i == 23, "X",
       ifelse(i == 24, "Y",
       ifelse(i == 25, "MT", as.character(i))))
  paste0("chr", v)
}

norm_sv <- function(x) substr(ifelse(toupper(x) == "BND", "TRA", toupper(x)), 1, 3)

# ------------------------------ #
# load COMMON SV file (5 or 6 c) #
# accepts tabs OR spaces         #
# ------------------------------ #
load_common <- function(path) {
  # first pass: try header=TRUE & auto-sep
  df <- tryCatch(
    suppressWarnings(fread(path, sep = "auto", header = TRUE, data.table = FALSE, fill = TRUE)),
    error = function(e) NULL
  )

  # helper to standardize names when present
  std_names <- function(nms) {
    nms <- tolower(gsub("^#+", "", nms))
    nms <- gsub("\\s+", "", nms)
    nms
  }

  if (!is.null(df) && ncol(df) >= 5) {
    nms <- std_names(colnames(df))
    # 6-col header
    if (all(c("chr1","pos1","chr2","end","pos2","svtype") %in% nms)) {
      names(df) <- nms
      pos2 <- suppressWarnings(as.integer(df$pos2))
      end  <- suppressWarnings(as.integer(df$end))
      pos2_fixed <- ifelse(is.na(pos2) | pos2 == 0L, end, pos2)
      out <- df %>%
        transmute(chr1 = chr1, pos1 = pos1, chr2 = chr2, pos2 = pos2_fixed, svtype = svtype)
      return(out)
    }
    # 5-col header
    if (all(c("chr1","pos1","chr2","pos2","svtype") %in% nms)) {
      names(df) <- nms
      out <- df %>%
        transmute(chr1 = chr1, pos1 = pos1, chr2 = chr2, pos2 = pos2, svtype = svtype)
      return(out)
    }
    # Some tools use "end" instead of pos2 in 5-col
    if (all(c("chr1","pos1","chr2","end","svtype") %in% nms)) {
      names(df) <- nms
      out <- df %>%
        transmute(chr1 = chr1, pos1 = pos1, chr2 = chr2, pos2 = end, svtype = svtype)
      return(out)
    }
    # If header exists but names are non-standard, fall through to no-header read below
  }

  # second pass: treat as no header, auto-sep
  df2 <- suppressWarnings(fread(path, sep = "auto", header = FALSE, data.table = FALSE, fill = TRUE))
  if (is.null(df2) || ncol(df2) < 5) {
    stop("common SV file must have either 5 columns: chr1 pos1 chr2 pos2 svtype OR 6 columns: chr1 pos1 chr2 end pos2 svtype")
  }

  # If 6+ columns, take first 6 as 6-col schema
  if (ncol(df2) >= 6) {
    colnames(df2)[1:6] <- c("chr1","pos1","chr2","end","pos2","svtype")
    pos2 <- suppressWarnings(as.integer(df2$pos2))
    end  <- suppressWarnings(as.integer(df2$end))
    pos2_fixed <- ifelse(is.na(pos2) | pos2 == 0L, end, pos2)
    out <- df2 %>%
      transmute(chr1 = .data$chr1, pos1 = .data$pos1, chr2 = .data$chr2, pos2 = pos2_fixed, svtype = .data$svtype)
    return(out)
  }

  # Else exactly 5 columns → map directly
  colnames(df2)[1:5] <- c("chr1","pos1","chr2","pos2","svtype")
  df2[,1:4] <- lapply(df2[,1:4, drop=FALSE], as.character)
  out <- df2 %>%
    transmute(chr1 = .data$chr1, pos1 = .data$pos1, chr2 = .data$chr2, pos2 = .data$pos2, svtype = .data$svtype)
  return(out)
}

common <- load_common(common_path) %>%
  mutate(
    svtype = norm_sv(svtype),
    chr1   = canon_chr(chr1),
    chr2   = canon_chr(chr2),
    pos1   = suppressWarnings(as.integer(pos1)),
    pos2   = suppressWarnings(as.integer(pos2))
  ) %>%
  filter(complete.cases(chr1, pos1, chr2, pos2, svtype)) %>%
  as.data.table()

# create windows
common[, `:=`(
  pos1_low  = pos1 - tolerance,
  pos1_high = pos1 + tolerance,
  pos2_low  = pos2 - tolerance,
  pos2_high = pos2 + tolerance
)]

# ------------------------------ #
# load AnnotSV TSV               #
# needs: SV_chrom, SV_start, INFO, SV_type, Annotation_mode
# ------------------------------ #
ann <- fread(annotsv_path, sep = "\t", header = TRUE, data.table = TRUE)
required_cols <- c("SV_chrom","SV_start","INFO","SV_type","Annotation_mode")
missing <- setdiff(required_cols, names(ann))
if (length(missing)) {
  stop("AnnotSV file is missing required column(s): ", paste(missing, collapse = ", "))
}

# keep only full annotations & extract CHR2/END from INFO
ann <- ann[Annotation_mode == "full"]

ann[, chr2_raw := {
  m <- str_match(INFO, "(?:^|;)CHR2=([^;]+)")
  ifelse(is.na(m[,2]), NA_character_, m[,2])
}]
ann[, end_raw := {
  m <- str_match(INFO, "(?:^|;)END=([0-9]+)")
  ifelse(is.na(m[,2]), NA_character_, m[,2])
}]

ann[, `:=`(
  chr1   = canon_chr(SV_chrom),
  pos1   = suppressWarnings(as.integer(SV_start)),
  chr2   = canon_chr(ifelse(is.na(chr2_raw), NA, ifelse(str_detect(chr2_raw, "^chr"), chr2_raw, paste0("chr", chr2_raw)))),
  pos2   = suppressWarnings(as.integer(end_raw)),
  svtype = norm_sv(SV_type)
)]

sample_dt <- ann[, .(chr1, pos1, chr2, pos2, svtype)][complete.cases(chr1,pos1,chr2,pos2,svtype)]
sample_dt[, rid := .I]

# ------------------------------ #
# non-equi join for overlaps     #
# chr1, chr2, svtype + position #
# ------------------------------ #
setkey(common, chr1, chr2, svtype, pos1_low, pos1_high, pos2_low, pos2_high)
matches <- sample_dt[common,
  on = .(chr1, chr2, svtype,
         pos1 >= pos1_low, pos1 <= pos1_high,
         pos2 >= pos2_low, pos2 <= pos2_high),
  nomatch = 0L,
  allow.cartesian = TRUE
]

overlap_idx <- unique(matches$rid)
overlap_dt  <- sample_dt[rid %in% overlap_idx]
tumor_only  <- sample_dt[!rid %in% overlap_idx]

# ------------------------------ #
# write outputs (CHROM, START)   #
# ------------------------------ #
overlap_bed <- if (nrow(overlap_dt)) data.table(CHROM = to_chr_prefix(overlap_dt$chr1), START = overlap_dt$pos1) else data.table(CHROM=character(), START=integer())
tumoronly_bed <- if (nrow(tumor_only)) data.table(CHROM = to_chr_prefix(tumor_only$chr1), START = tumor_only$pos1) else data.table(CHROM=character(), START=integer())

overlap_path   <- file.path(outdir, paste0(sample_prefix, ".overlapping_SV_calls.bed"))
tumoronly_path <- file.path(outdir, paste0(sample_prefix, ".tumor_only_SV_calls.bed"))

fwrite(overlap_bed,   overlap_path,   sep = "\t", col.names = FALSE)
fwrite(tumoronly_bed, tumoronly_path, sep = "\t", col.names = FALSE)

message("Wrote: ", overlap_path)
message("Wrote: ", tumoronly_path)
