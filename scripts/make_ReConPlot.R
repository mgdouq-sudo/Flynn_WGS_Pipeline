# install.packages("BiocManager"); BiocManager::install("VariantAnnotation")
suppressPackageStartupMessages({
  library(VariantAnnotation)
})

parse_survivor_vcf <- function(vcf_path) {
  vcf <- readVcf(vcf_path)
  gr  <- rowRanges(vcf)
  inf <- info(vcf)
  
  # Safely fetch fields that might not exist
  get_info <- function(x, name) if (name %in% colnames(x)) x[, name] else rep(NA, nrow(x))
  
  svtype   <- as.character(get_info(inf, "SVTYPE"))
  end_pos  <- as.integer(get_info(inf, "END"))
  strandsI <- as.character(get_info(inf, "STRANDS"))  # many callers keep this, SURVIVOR often preserves it
  
  chrom1 <- as.character(seqnames(gr))
  pos1   <- as.integer(start(gr))
  
  # Initialize outputs
  chrom2 <- rep(".", length(chrom1))
  pos2   <- rep(NA_integer_, length(chrom1))
  strands <- rep(NA_character_, length(chrom1))
  
  # 1) Insertions: no second coordinate in your desired format
  is_ins <- svtype == "INS"
  chrom2[is_ins]  <- "."
  pos2[is_ins]    <- NA_integer_
  strands[is_ins] <- "INS"
  
  # 2) Single-breakends (some tools encode as SVTYPE=SGL or BND without mate; SURVIVOR may tag SBE)
  is_sbe <- svtype %in% c("SGL", "SBE")
  chrom2[is_sbe]  <- "."
  pos2[is_sbe]    <- NA_integer_
  strands[is_sbe] <- "SBE"
  
  # 3) Simple intra-chromosomal events with END (DEL/DUP/INV/CNV/etc.)
  has_end <- !is.na(end_pos) & !is_ins & !is_sbe
  chrom2[has_end] <- chrom1[has_end]
  pos2[has_end]   <- end_pos[has_end]
  
  # If STRANDS present, use it; otherwise provide a sensible fallback by SVTYPE
  fallback_strand <- function(t) {
    if (t == "DEL") return("+-")
    if (t == "DUP") return("-+")
    if (t == "INV") return("++")   # convention varies; INFO/STRANDS is preferred
    if (t == "CNV") return("+-")   # generic fallback
    return(NA_character_)
  }
  need_strand <- is.na(strands) & has_end
  strands[need_strand] <- ifelse(!is.na(strandsI[need_strand]), strandsI[need_strand],
                                 vapply(svtype[need_strand], fallback_strand, character(1)))
  
  # 4) Inter-chromosomal BNDs: parse ALT to get mate chr:pos if END missing
  is_bnd <- svtype == "BND"
  if (any(is_bnd)) {
    alt <- as.character(unlist(alt(vcf)[is_bnd]))
    # ALT encodes mate like one of: ]chr:pos]SEQ, SEQ]chr:pos], [chr:pos[SEQ, SEQ[chr:pos[
    # Extract chr:pos robustly:
    m <- regexpr("(\\[|\\])([^:\\[\\]]+):(\\d+)(\\[|\\])", alt, perl = TRUE)
    mate_chr <- rep(NA_character_, sum(is_bnd))
    mate_pos <- rep(NA_integer_,  sum(is_bnd))
    has_mate <- m != -1
    if (any(has_mate)) {
      caps <- regmatches(alt, m)
      # Pull the chr:pos inside the brackets
      mate_chr[has_mate] <- sub(".*[\\[\\]]([^:\\[\\]]+):(\\d+)[\\[\\]].*", "\\1", caps[has_mate], perl = TRUE)
      mate_pos[has_mate] <- as.integer(sub(".*[\\[\\]]([^:\\[\\]]+):(\\d+)[\\[\\]].*", "\\2", caps[has_mate], perl = TRUE))
    }
    
    idx_bnd <- which(is_bnd)
    chrom2[idx_bnd[has_mate]] <- mate_chr[has_mate]
    pos2[idx_bnd[has_mate]]   <- mate_pos[has_mate]
    
    # Prefer INFO/STRANDS if present
    sI_bnd <- strandsI[is_bnd]
    strands[is_bnd & !is.na(sI_bnd)] <- sI_bnd[!is.na(sI_bnd)]
    
    # If STRANDS absent, leave as NA (deriving from ALT bracket orientation is caller-specific and brittle).
  }
  
  # Finalize dataframe with exact column names/order
  out <- data.frame(
    chr1    = chrom1,
    pos1    = pos1,
    chr2    = chrom2,
    pos2    = pos2,
    strands = strands,
    stringsAsFactors = FALSE
  )
  
  # Make dots for missing pos2 like in your example
  out$pos2 <- ifelse(is.na(out$pos2), ".", as.character(out$pos2))
  
  out
}

# ---- Example ----
# df <- parse_survivor_vcf("merged.survivor.vcf.gz")
# head(df)
