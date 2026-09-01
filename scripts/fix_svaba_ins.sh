#!/usr/bin/env bash
set -euo pipefail

vcf="$1"
out="$2"

{
  # print header exactly as-is
  grep '^#' "$vcf"

  # process body
  grep -v '^#' "$vcf" | \
  awk '
    {
      # only if this row contains INSERTION anywhere
      if ($0 ~ /INSERTION=/) {

        # column 8 is INFO
        # ONLY if INFO ends with BND
        if ($8 ~ /BND$/) {
          sub(/BND$/, "INS", $8)
        }
      }

      print
    }
  '
} > "$out"
