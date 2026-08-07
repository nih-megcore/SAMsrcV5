#!/usr/bin/env bash
set -euo pipefail

: "${SAM_DEPS_PREFIX:?Set SAM_DEPS_PREFIX to the dependency install prefix}"
: "${SAM_DEPS_CACHE:?Set SAM_DEPS_CACHE to the dependency archive cache}"

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly JOBS="${SAM_BUILD_JOBS:-2}"
readonly WORK_DIR="$(mktemp -d)"
trap 'rm -rf "${WORK_DIR}"' EXIT

case "$(uname -s)" in
    MINGW*|MSYS*) readonly PIC_FLAG="" ;;
    *) readonly PIC_FLAG="-fPIC" ;;
esac

"${SCRIPT_DIR}/fetch_native_deps.sh"
mkdir -p "${SAM_DEPS_PREFIX}"

tar -xzf "${SAM_DEPS_CACHE}/gsl-2.8.tar.gz" -C "${WORK_DIR}"
(
    cd "${WORK_DIR}/gsl-2.8"
    ./configure \
        --prefix="${SAM_DEPS_PREFIX}" \
        --disable-shared \
        --enable-static \
        CFLAGS="${CFLAGS:--O2} ${PIC_FLAG}"
    make -j"${JOBS}"
    make install
    mkdir -p "${SAM_DEPS_PREFIX}/share/samsrcv5-licenses"
    cp COPYING "${SAM_DEPS_PREFIX}/share/samsrcv5-licenses/GSL-GPL-3.0.txt"
)

tar -xzf "${SAM_DEPS_CACHE}/fftw-3.3.11.tar.gz" -C "${WORK_DIR}"
(
    cd "${WORK_DIR}/fftw-3.3.11"
    ./configure \
        --prefix="${SAM_DEPS_PREFIX}" \
        --disable-fortran \
        --disable-shared \
        --enable-static \
        CFLAGS="${CFLAGS:--O2} ${PIC_FLAG}"
    make -j"${JOBS}"
    make install
    mkdir -p "${SAM_DEPS_PREFIX}/share/samsrcv5-licenses"
    cp COPYING "${SAM_DEPS_PREFIX}/share/samsrcv5-licenses/FFTW-GPL-2.0-or-later.txt"
)
