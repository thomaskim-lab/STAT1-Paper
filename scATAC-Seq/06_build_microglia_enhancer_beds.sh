#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# 06_build_microglia_enhancer_beds.sh
# Builds permissive and strict microglial enhancer BED files from Ciernia-derived
# H3K27ac and H3K4me1 BEDs. Requires bedtools.
#
# Definitions:
#   permissive enhancer set:
#     merged H3K27ac peaks excluding genome-wide promoters
#
#   strict enhancer set:
#     merged H3K27ac peaks overlapping H3K4me1, excluding genome-wide promoters
#
# These are external microglial enhancer annotations only.
# They are not evidence of direct STAT1/IRF1 binding.
#
# GenomeDK / Hypothalamus version.
# =============================================================================

PROJECT=/faststorage/project/Hypothalamus
DATA=${PROJECT}/data
SCATAC_DIR=${DATA}/scATAC/merged_microglia_STAT1Project
FIG2_DIR=${SCATAC_DIR}/Fig2_scATAC

CIERNIA_DIR=${FIG2_DIR}/data/ciernia
CIERNIA_BED_DIR=${CIERNIA_DIR}/bed
OUT_DIR=${CIERNIA_DIR}/processed

PROMOTERS_2KB=${FIG2_DIR}/bed/all_promoters_2kb_merged.bed

LOG_DIR=${FIG2_DIR}/logs

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

echo "============================================================"
echo "06_build_microglia_enhancer_beds.sh"
echo "Started: $(date)"
echo "Running on: $(hostname)"
echo "============================================================"

echo "PROJECT: ${PROJECT}"
echo "FIG2_DIR: ${FIG2_DIR}"
echo "CIERNIA_BED_DIR: ${CIERNIA_BED_DIR}"
echo "OUT_DIR: ${OUT_DIR}"
echo "PROMOTERS_2KB: ${PROMOTERS_2KB}"

# =============================================================================
# Checks
# =============================================================================

if ! command -v bedtools >/dev/null 2>&1; then
  echo "ERROR: bedtools is not installed or not in PATH." >&2
  echo "Activate my_r_env first:" >&2
  echo "  source /home/tkim/miniforge3/etc/profile.d/conda.sh" >&2
  echo "  conda activate my_r_env" >&2
  exit 1
fi

echo "bedtools: $(command -v bedtools)"
bedtools --version

if [ ! -d "${CIERNIA_BED_DIR}" ]; then
  echo "ERROR: Ciernia BED directory not found:" >&2
  echo "  ${CIERNIA_BED_DIR}" >&2
  echo "Run 05_download_ciernia_tracks.sh first." >&2
  exit 1
fi

if [ ! -f "${PROMOTERS_2KB}" ]; then
  echo "ERROR: promoter BED not found:" >&2
  echo "  ${PROMOTERS_2KB}" >&2
  echo "Run 01_export_peakunion_and_promoters.R first." >&2
  exit 1
fi

if [ ! -s "${PROMOTERS_2KB}" ]; then
  echo "ERROR: promoter BED exists but is empty:" >&2
  echo "  ${PROMOTERS_2KB}" >&2
  exit 1
fi

# =============================================================================
# Helper paths
# =============================================================================

H3K27AC_FILE_LIST=${OUT_DIR}/h3k27ac_bed_files.txt
H3K4ME1_FILE_LIST=${OUT_DIR}/h3k4me1_bed_files.txt

H3K27AC_CLEAN=${OUT_DIR}/microglia_H3K27ac_clean.sorted.bed
H3K4ME1_CLEAN=${OUT_DIR}/microglia_H3K4me1_clean.sorted.bed

H3K27AC_MERGED=${OUT_DIR}/microglia_H3K27ac_merged.bed
H3K4ME1_MERGED=${OUT_DIR}/microglia_H3K4me1_merged.bed

PERMISSIVE_ENHANCERS=${OUT_DIR}/microglia_enhancers_permissive_H3K27ac_noPromoter.bed
STRICT_OVERLAP=${OUT_DIR}/microglia_H3K27ac_H3K4me1_overlap.bed
STRICT_ENHANCERS=${OUT_DIR}/microglia_enhancers_strict_H3K27ac_H3K4me1_noPromoter.bed

SUMMARY=${OUT_DIR}/microglia_enhancer_bed_summary.tsv

# =============================================================================
# Locate H3K27ac and H3K4me1 BED files
# =============================================================================

echo "Finding H3K27ac BED files..."

find "${CIERNIA_BED_DIR}" -type f \( \
    -iname '*H3K27ac*.bed' -o \
    -iname '*H3k27ac*.bed' -o \
    -iname '*K27ac*.bed' \
  \) | sort > "${H3K27AC_FILE_LIST}"

if [ ! -s "${H3K27AC_FILE_LIST}" ]; then
  echo "ERROR: no H3K27ac BED files found in:" >&2
  echo "  ${CIERNIA_BED_DIR}" >&2
  echo "Available BED files:" >&2
  find "${CIERNIA_BED_DIR}" -type f -iname '*.bed' | sort >&2
  exit 1
fi

echo "H3K27ac BED files:"
cat "${H3K27AC_FILE_LIST}"

echo "Finding H3K4me1 BED files..."

find "${CIERNIA_BED_DIR}" -type f \( \
    -iname '*H3K4me1*.bed' -o \
    -iname '*H3k4me1*.bed' -o \
    -iname '*K4me1*.bed' \
  \) | sort > "${H3K4ME1_FILE_LIST}"

if [ -s "${H3K4ME1_FILE_LIST}" ]; then
  echo "H3K4me1 BED files:"
  cat "${H3K4ME1_FILE_LIST}"
else
  echo "WARNING: no H3K4me1 BED files found."
  echo "Strict enhancer set will not be created."
fi

# =============================================================================
# Clean, sort, and merge H3K27ac peaks
# =============================================================================

echo "Cleaning and sorting H3K27ac BED intervals..."

xargs -r cat < "${H3K27AC_FILE_LIST}" \
  | awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X|Y)$/ && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $3 > $2 {print $1,$2,$3}' \
  | sort -k1,1 -k2,2n \
  > "${H3K27AC_CLEAN}"

if [ ! -s "${H3K27AC_CLEAN}" ]; then
  echo "ERROR: cleaned H3K27ac BED is empty:" >&2
  echo "  ${H3K27AC_CLEAN}" >&2
  exit 1
fi

echo "Merging H3K27ac intervals..."

bedtools merge \
  -i "${H3K27AC_CLEAN}" \
  > "${H3K27AC_MERGED}"

if [ ! -s "${H3K27AC_MERGED}" ]; then
  echo "ERROR: merged H3K27ac BED is empty:" >&2
  echo "  ${H3K27AC_MERGED}" >&2
  exit 1
fi

# =============================================================================
# Build permissive enhancer set: H3K27ac minus promoters
# =============================================================================

echo "Building permissive enhancer set: H3K27ac minus promoters..."

bedtools subtract \
  -a "${H3K27AC_MERGED}" \
  -b "${PROMOTERS_2KB}" \
  | awk 'BEGIN{OFS="\t"} $3 > $2 {print $1,$2,$3}' \
  | sort -k1,1 -k2,2n \
  | bedtools merge \
  > "${PERMISSIVE_ENHANCERS}"

if [ ! -s "${PERMISSIVE_ENHANCERS}" ]; then
  echo "ERROR: permissive enhancer BED is empty:" >&2
  echo "  ${PERMISSIVE_ENHANCERS}" >&2
  exit 1
fi

# Validate no promoter overlap remains.
if bedtools intersect -u -a "${PERMISSIVE_ENHANCERS}" -b "${PROMOTERS_2KB}" | head -1 | grep -q .; then
  echo "ERROR: permissive enhancer BED still overlaps promoters." >&2
  exit 1
fi

echo "Permissive enhancer BED:"
ls -lh "${PERMISSIVE_ENHANCERS}"

# =============================================================================
# Build strict enhancer set if H3K4me1 is available
# =============================================================================

STRICT_CREATED=0

if [ -s "${H3K4ME1_FILE_LIST}" ]; then
  echo "Cleaning and sorting H3K4me1 BED intervals..."

  xargs -r cat < "${H3K4ME1_FILE_LIST}" \
    | awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X|Y)$/ && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $3 > $2 {print $1,$2,$3}' \
    | sort -k1,1 -k2,2n \
    > "${H3K4ME1_CLEAN}"

  if [ ! -s "${H3K4ME1_CLEAN}" ]; then
    echo "WARNING: cleaned H3K4me1 BED is empty. Strict enhancer set will not be created." >&2
  else
    echo "Merging H3K4me1 intervals..."

    bedtools merge \
      -i "${H3K4ME1_CLEAN}" \
      > "${H3K4ME1_MERGED}"

    if [ ! -s "${H3K4ME1_MERGED}" ]; then
      echo "WARNING: merged H3K4me1 BED is empty. Strict enhancer set will not be created." >&2
    else
      echo "Building strict overlap: H3K27ac overlapping H3K4me1..."

      bedtools intersect \
        -a "${H3K27AC_MERGED}" \
        -b "${H3K4ME1_MERGED}" \
        -u \
        | sort -k1,1 -k2,2n \
        | bedtools merge \
        > "${STRICT_OVERLAP}"

      if [ ! -s "${STRICT_OVERLAP}" ]; then
        echo "WARNING: no H3K27ac/H3K4me1 overlap found. Strict enhancer set will not be created." >&2
      else
        echo "Subtracting promoters from strict enhancer overlap..."

        bedtools subtract \
          -a "${STRICT_OVERLAP}" \
          -b "${PROMOTERS_2KB}" \
          | awk 'BEGIN{OFS="\t"} $3 > $2 {print $1,$2,$3}' \
          | sort -k1,1 -k2,2n \
          | bedtools merge \
          > "${STRICT_ENHANCERS}"

        if [ -s "${STRICT_ENHANCERS}" ]; then
          if bedtools intersect -u -a "${STRICT_ENHANCERS}" -b "${PROMOTERS_2KB}" | head -1 | grep -q .; then
            echo "ERROR: strict enhancer BED still overlaps promoters." >&2
            exit 1
          fi

          STRICT_CREATED=1

          echo "Strict enhancer BED:"
          ls -lh "${STRICT_ENHANCERS}"
        else
          echo "WARNING: strict enhancer BED is empty after promoter subtraction." >&2
        fi
      fi
    fi
  fi
fi

# =============================================================================
# Summary
# =============================================================================

echo -e "file\tlabel\tn_intervals\ttotal_bp" > "${SUMMARY}"

summarise_bed() {
  local bed_file="$1"
  local label="$2"

  if [ -s "${bed_file}" ]; then
    local n
    local bp
    n=$(wc -l < "${bed_file}")
    bp=$(awk '{sum += ($3 - $2)} END {print sum + 0}' "${bed_file}")
    echo -e "${bed_file}\t${label}\t${n}\t${bp}" >> "${SUMMARY}"
  else
    echo -e "${bed_file}\t${label}\t0\t0" >> "${SUMMARY}"
  fi
}

summarise_bed "${H3K27AC_CLEAN}" "H3K27ac_clean"
summarise_bed "${H3K27AC_MERGED}" "H3K27ac_merged"
summarise_bed "${PERMISSIVE_ENHANCERS}" "permissive_H3K27ac_noPromoter"

if [ -s "${H3K4ME1_CLEAN:-}" ]; then
  summarise_bed "${H3K4ME1_CLEAN}" "H3K4me1_clean"
fi

if [ -s "${H3K4ME1_MERGED:-}" ]; then
  summarise_bed "${H3K4ME1_MERGED}" "H3K4me1_merged"
fi

if [ -s "${STRICT_OVERLAP:-}" ]; then
  summarise_bed "${STRICT_OVERLAP}" "H3K27ac_H3K4me1_overlap"
fi

if [ "${STRICT_CREATED}" -eq 1 ]; then
  summarise_bed "${STRICT_ENHANCERS}" "strict_H3K27ac_H3K4me1_noPromoter"
fi

echo "Enhancer BED summary:"
cat "${SUMMARY}"

echo "Processed enhancer BEDs:"
ls -lh "${OUT_DIR}"

echo "============================================================"
echo "06_build_microglia_enhancer_beds.sh complete"
echo "Finished: $(date)"
echo "============================================================"