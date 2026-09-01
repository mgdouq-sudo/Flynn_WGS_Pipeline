# shatterseek.R
library(ShatterSeek)

# Example: load SVs and CNVs
svs <- read.delim("svs.bedpe")
cnvs <- read.delim("cnvs.bed")

# Run ShatterSeek on a given chromosome or sample
results <- shatterseek(
  SV = svs,
  CNV = cnvs,
  genome = "hg38"
)

# Save results
write.table(results, "shatterseek_results.tsv", sep = "\t", row.names = FALSE)
