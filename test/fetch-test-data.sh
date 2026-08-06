#!/usr/bin/env bash
set -euo pipefail

readonly REVISION="1d08e47c586fca21163c4e7362d409e62b1c1943"
readonly REPOSITORY="https://github.com/nih-megcore/TEST_ctf_data.git"
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly DATA_DIR="${TEST_DATA_DIR:-${SCRIPT_DIR}/../.test-data}"

verify_data() {
    [[ -d "${DATA_DIR}/.git" ]] || return 1
    [[ "$(git -C "${DATA_DIR}" rev-parse HEAD)" == "${REVISION}" ]] || return 1
    [[ "$(git -C "${DATA_DIR}" rev-parse 'HEAD:MRI/ABABABAB_refaced+orig.HEAD')" == "ae97b40cf210985ca935d5532415e171a37649a0" ]] || return 1
    [[ "$(git -C "${DATA_DIR}" rev-parse 'HEAD:MRI/ABABABAB_refaced+orig.BRIK.gz')" == "eb18f9a7015e9095d67865b7ca74852125865280" ]] || return 1
    [[ "$(git -C "${DATA_DIR}" rev-parse 'HEAD:20010101/ABABABAB_airpuff_20010101_001.ds')" == "d497adbd07caf94c3eb0824b6667d2e32da1c687" ]] || return 1
    [[ -f "${DATA_DIR}/MRI/ABABABAB_refaced+orig.HEAD" ]] || return 1
    [[ -f "${DATA_DIR}/MRI/ABABABAB_refaced+orig.BRIK.gz" ]] || return 1
    [[ -f "${DATA_DIR}/20010101/ABABABAB_airpuff_20010101_001.ds/ABABABAB_airpuff_20010101_001.meg4" ]] || return 1
}

if verify_data; then
    printf 'Test data already verified at %s\n' "${DATA_DIR}"
    exit 0
fi

if [[ -e "${DATA_DIR}" && ! -d "${DATA_DIR}/.git" ]]; then
    printf 'Refusing to replace non-git path: %s\n' "${DATA_DIR}" >&2
    exit 1
fi

if [[ ! -d "${DATA_DIR}/.git" ]]; then
    git clone --filter=blob:none --no-checkout "${REPOSITORY}" "${DATA_DIR}"
fi

git -C "${DATA_DIR}" sparse-checkout init --no-cone
git -C "${DATA_DIR}" sparse-checkout set --no-cone \
    '/MRI/ABABABAB_refaced+orig.HEAD' \
    '/MRI/ABABABAB_refaced+orig.BRIK.gz' \
    '/20010101/ABABABAB_airpuff_20010101_001.ds/'
git -C "${DATA_DIR}" fetch --depth 1 origin "${REVISION}"
git -C "${DATA_DIR}" checkout --detach FETCH_HEAD

if ! verify_data; then
    printf 'Downloaded test data failed verification\n' >&2
    exit 1
fi
printf 'Downloaded and verified test data at %s\n' "${DATA_DIR}"
