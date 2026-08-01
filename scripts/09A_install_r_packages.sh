#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_DIR="${1:?Usage: 09A_install_r_packages.sh <project_dir>}"
SOURCE_DIR="${PROJECT_DIR}/environment/package_sources/stage9A_bioc318"
LIB_DIR="${PROJECT_DIR}/environment/Rlib_stage9A"
PREPROCESSCORE_SOURCE="${PROJECT_DIR}/environment/sources/preprocessCore_1.64.0.tar.gz"
mkdir -p "${LIB_DIR}"

required_sources=(
  ff_4.5.3.tar.gz
  affxparser_1.74.0.tar.gz
  oligoClasses_1.64.0.tar.gz
  oligo_1.66.0.tar.gz
  pd.hta.2.0_3.12.2.tar.gz
  hta20transcriptcluster.db_8.8.0.tar.gz
)
for source in "${required_sources[@]}"; do
  test -s "${SOURCE_DIR}/${source}"
done
test -s "${PREPROCESSCORE_SOURCE}"
printf '%s  %s\n' \
  '2116c6363074b59becdaf7a1e88caf91' \
  "${PREPROCESSCORE_SOURCE}" | md5sum -c -
(
  cd "${SOURCE_DIR}"
  md5sum -c <<'CHECKSUMS'
1b03d5fbdda174d9d106db44767b8dc5  ff_4.5.3.tar.gz
2b6b9373d749a0ccbf930f834df90b92  affxparser_1.74.0.tar.gz
379f1364d7159e2456c9310556e990e6  oligoClasses_1.64.0.tar.gz
3ba99dedbe0ebca2fcaa81a400777469  oligo_1.66.0.tar.gz
8e13f85ece49c38da73eaf7b2247f5f0  pd.hta.2.0_3.12.2.tar.gz
9debbc190cc5f9f0ad97fbc452fe7025  hta20transcriptcluster.db_8.8.0.tar.gz
CHECKSUMS
)

export R_LIBS_USER="${LIB_DIR}"
# This server rejects preprocessCore's default pthread creation during RMA.
# Keep a project-local, explicitly non-threaded build ahead of all other
# libraries. Checking lib.loc prevents the system threaded build from being
# mistaken for the required server-compatible package.
if ! Rscript -e "quit(status=if(length(find.package('preprocessCore', lib.loc='${LIB_DIR}', quiet=TRUE))) 0 else 1)"; then
  R CMD INSTALL \
    --configure-args="--disable-threading" \
    --library="${LIB_DIR}" \
    "${PREPROCESSCORE_SOURCE}"
fi

install_if_missing() {
  local package="$1"
  local source="$2"
  if ! Rscript -e "quit(status=if(requireNamespace('${package}', quietly=TRUE)) 0 else 1)"; then
    R CMD INSTALL --library="${LIB_DIR}" "${SOURCE_DIR}/${source}"
  fi
}

install_if_missing ff ff_4.5.3.tar.gz
install_if_missing affxparser affxparser_1.74.0.tar.gz
install_if_missing oligoClasses oligoClasses_1.64.0.tar.gz
install_if_missing oligo oligo_1.66.0.tar.gz
install_if_missing pd.hta.2.0 pd.hta.2.0_3.12.2.tar.gz
install_if_missing hta20transcriptcluster.db hta20transcriptcluster.db_8.8.0.tar.gz

Rscript - <<'RS'
required <- c(
  "preprocessCore", "ff", "affxparser", "oligoClasses", "oligo", "pd.hta.2.0",
  "hta20transcriptcluster.db", "AnnotationDbi", "Biobase", "matrixStats"
)
ok <- vapply(required, requireNamespace, logical(1), quietly = TRUE)
print(data.frame(package = required, installed = ok), row.names = FALSE)
if (!all(ok)) stop("One or more Stage 9A R packages are unavailable")
RS
