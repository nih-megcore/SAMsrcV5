#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
readonly DATA_DIR="${TEST_DATA_DIR:-${ROOT_DIR}/.test-data}"
readonly RESULTS_DIR="${TEST_RESULTS_DIR:-${ROOT_DIR}/.test-results}"
readonly CCACHE_HOST_DIR="${CCACHE_HOST_DIR:-${ROOT_DIR}/.ccache/afni}"
readonly AFNI_IMAGE="${AFNI_IMAGE:-ghcr.io/jstout211/afni:latest}"
readonly TEST_IMAGE="${SAM_TEST_IMAGE:-sam2multi-afni-test:local}"
readonly PYTEST_MARKER="${SAM_PYTEST_MARKER:-integration}"
readonly TEST_REPORT="${SAM_TEST_REPORT:-integration-junit.xml}"

if [[ -n "${CONTAINER_ENGINE:-}" ]]; then
    ENGINE="${CONTAINER_ENGINE}"
elif command -v podman >/dev/null 2>&1; then
    ENGINE="podman"
elif command -v docker >/dev/null 2>&1; then
    ENGINE="docker"
else
    printf 'Neither podman nor docker was found in PATH\n' >&2
    exit 1
fi

TEST_DATA_DIR="${DATA_DIR}" "${SCRIPT_DIR}/fetch-test-data.sh"
mkdir -p "${RESULTS_DIR}" "${CCACHE_HOST_DIR}"

if [[ "${BUILD_TEST_IMAGE:-1}" == "1" ]]; then
    "${ENGINE}" build \
        --build-arg "AFNI_IMAGE=${AFNI_IMAGE}" \
        --file "${SCRIPT_DIR}/container/test.Dockerfile" \
        --tag "${TEST_IMAGE}" \
        "${ROOT_DIR}"
fi

RUN_OPTIONS=(--rm)
if [[ "$(basename "${ENGINE}")" == "podman" ]]; then
    RUN_OPTIONS+=(--security-opt label=disable)
fi

"${ENGINE}" run "${RUN_OPTIONS[@]}" \
    --env CCACHE_DIR=/ccache \
    --env TEST_DATA_DIR=/test-data \
    --env TEST_RESULTS_DIR=/results \
    --env "SAM_PYTEST_MARKER=${PYTEST_MARKER}" \
    --env "SAM_TEST_REPORT=${TEST_REPORT}" \
    --volume "${DATA_DIR}:/test-data:ro" \
    --volume "${RESULTS_DIR}:/results" \
    --volume "${CCACHE_HOST_DIR}:/ccache" \
    "${TEST_IMAGE}" \
    bash -lc 'set -euo pipefail
cd /work
make -j"$(nproc)" CC="ccache gcc"
make -C test integration-build CC="ccache gcc"
python3 -m pytest -m "$SAM_PYTEST_MARKER" -vv \
    --junitxml="/results/$SAM_TEST_REPORT" \
    test/integration'
