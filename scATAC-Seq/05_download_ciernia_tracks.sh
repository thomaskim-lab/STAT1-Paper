#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# 05_download_ciernia_tracks.sh
# Downloads Ciernia/MicrogliaChIPseq track hub file, extracts H3K27ac/H3K4me1
# bigBed URLs, downloads .bb files, and converts them to BED.
#
# GenomeDK / Hypothalamus version.
# =============================================================================

PROJECT=/faststorage/project/Hypothalamus
DATA=${PROJECT}/data
SCATAC_DIR=${DATA}/scATAC/merged_microglia_STAT1Project
FIG2_DIR=${SCATAC_DIR}/Fig2_scATAC

TOOLS_DIR=${PROJECT}/tools
UCSC_DIR=${TOOLS_DIR}/ucsc

CIERNIA_DIR=${FIG2_DIR}/data/ciernia
HUB_DIR=${CIERNIA_DIR}/hub
BB_DIR=${CIERNIA_DIR}/bb
BED_DIR=${CIERNIA_DIR}/bed
LOG_DIR=${FIG2_DIR}/logs

BIGBEDTOBED=${UCSC_DIR}/bigBedToBed

mkdir -p "${HUB_DIR}" "${BB_DIR}" "${BED_DIR}" "${LOG_DIR}"

HUB_URL="https://raw.githubusercontent.com/ciernialab/MicrogliaChIPseq/master/TrackHubmm10_MGEpi_Master.txt"
HUB_FILE="${HUB_DIR}/TrackHubmm10_MGEpi_Master.txt"

ALL_URLS="${HUB_DIR}/all_bigbed_urls.txt"
ENHANCER_URLS="${HUB_DIR}/enhancer_bigbed_urls.txt"
DOWNLOAD_MANIFEST="${HUB_DIR}/download_manifest.tsv"
CONVERSION_MANIFEST="${HUB_DIR}/conversion_manifest.tsv"

echo "============================================================"
echo "05_download_ciernia_tracks.sh"
echo "Started: $(date)"
echo "Running on: $(hostname)"
echo "============================================================"

echo "PROJECT: ${PROJECT}"
echo "FIG2_DIR: ${FIG2_DIR}"
echo "CIERNIA_DIR: ${CIERNIA_DIR}"
echo "BIGBEDTOBED: ${BIGBEDTOBED}"

# =============================================================================
# Checks
# =============================================================================

if ! command -v curl >/dev/null 2>&1; then
  echo "ERROR: curl not found in PATH." >&2
  exit 1
fi

if [ ! -x "${BIGBEDTOBED}" ]; then
  echo "ERROR: bigBedToBed not found or not executable:" >&2
  echo "  ${BIGBEDTOBED}" >&2
  echo "Install/check it first, e.g.:" >&2
  echo "  ls -lh ${BIGBEDTOBED}" >&2
  exit 1
fi

echo "curl: $(command -v curl)"
echo "bigBedToBed: ${BIGBEDTOBED}"

# =============================================================================
# Download track hub
# =============================================================================

echo "Downloading Ciernia track hub:"
echo "  ${HUB_URL}"

curl -L --fail --retry 3 --retry-delay 5 \
  "${HUB_URL}" \
  -o "${HUB_FILE}"

if [ ! -s "${HUB_FILE}" ]; then
  echo "ERROR: Hub file was not downloaded or is empty:" >&2
  echo "  ${HUB_FILE}" >&2
  exit 1
fi

echo "Downloaded hub file:"
ls -lh "${HUB_FILE}"

# =============================================================================
# Extract bigBed URLs
# =============================================================================

# The hub file can be irregularly spaced, so split on whitespace.
# Keep http/https URLs ending in .bb.
tr '[:space:]' '\n' < "${HUB_FILE}" \
  | grep -E '^https?://.*\.bb$' \
  | sort -u \
  > "${ALL_URLS}" || true

if [ ! -s "${ALL_URLS}" ]; then
  echo "ERROR: No .bb URLs found in hub file." >&2
  echo "Inspect hub file manually:" >&2
  echo "  ${HUB_FILE}" >&2
  exit 1
fi

echo "All bigBed URLs found: $(wc -l < "${ALL_URLS}")"
cat "${ALL_URLS}"

# Keep enhancer-related histone marks.
grep -Ei 'H3K27ac|H3K4me1|H3k27ac|H3k4me1|K27ac|K4me1' \
  "${ALL_URLS}" \
  | sort -u \
  > "${ENHANCER_URLS}" || true

if [ ! -s "${ENHANCER_URLS}" ]; then
  echo "ERROR: No H3K27ac/H3K4me1 bigBed URLs found." >&2
  echo "All URLs are in:" >&2
  echo "  ${ALL_URLS}" >&2
  exit 1
fi

echo "Enhancer-related bigBed URLs found: $(wc -l < "${ENHANCER_URLS}")"
cat "${ENHANCER_URLS}"

# =============================================================================
# Download .bb files
# =============================================================================

echo -e "url\tbb_file\tstatus\tbytes" > "${DOWNLOAD_MANIFEST}"

while read -r url; do
  [ -z "${url}" ] && continue

  fname=$(basename "${url}")
  out_bb="${BB_DIR}/${fname}"

  echo "Downloading:"
  echo "  URL: ${url}"
  echo "  OUT: ${out_bb}"

  if [ -s "${out_bb}" ]; then
    echo "  Existing non-empty file found; skipping download."
  else
    curl -L --fail --retry 3 --retry-delay 5 \
      "${url}" \
      -o "${out_bb}"
  fi

  if [ ! -s "${out_bb}" ]; then
    echo "ERROR: Download failed or file is empty:" >&2
    echo "  ${out_bb}" >&2
    echo -e "${url}\t${out_bb}\tFAILED\t0" >> "${DOWNLOAD_MANIFEST}"
    exit 1
  fi

  bytes=$(stat -c%s "${out_bb}")
  echo -e "${url}\t${out_bb}\tOK\t${bytes}" >> "${DOWNLOAD_MANIFEST}"

done < "${ENHANCER_URLS}"

echo "Download manifest:"
cat "${DOWNLOAD_MANIFEST}"

# =============================================================================
# Convert .bb to .bed
# =============================================================================

echo -e "bb_file\tbed_file\tstatus\tbed_lines" > "${CONVERSION_MANIFEST}"

shopt -s nullglob

bb_files=("${BB_DIR}"/*.bb)

if [ "${#bb_files[@]}" -eq 0 ]; then
  echo "ERROR: No .bb files found in:" >&2
  echo "  ${BB_DIR}" >&2
  exit 1
fi

for bb in "${bb_files[@]}"; do
  base=$(basename "${bb}" .bb)
  out_bed="${BED_DIR}/${base}.bed"

  echo "Converting:"
  echo "  BB:  ${bb}"
  echo "  BED: ${out_bed}"

  if [ -s "${out_bed}" ]; then
    echo "  Existing non-empty BED found; overwriting to ensure consistency."
    rm -f "${out_bed}"
  fi

  "${BIGBEDTOBED}" "${bb}" "${out_bed}"

  if [ ! -s "${out_bed}" ]; then
    echo "ERROR: Conversion failed or output BED is empty:" >&2
    echo "  ${out_bed}" >&2
    echo -e "${bb}\t${out_bed}\tFAILED\t0" >> "${CONVERSION_MANIFEST}"
    exit 1
  fi

  bed_lines=$(wc -l < "${out_bed}")
  echo -e "${bb}\t${out_bed}\tOK\t${bed_lines}" >> "${CONVERSION_MANIFEST}"

done

echo "Conversion manifest:"
cat "${CONVERSION_MANIFEST}"

echo "Converted BED files:"
ls -lh "${BED_DIR}"

echo "============================================================"
echo "05_download_ciernia_tracks.sh complete"
echo "Finished: $(date)"
echo "============================================================"