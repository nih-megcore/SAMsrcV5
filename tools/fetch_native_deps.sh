#!/usr/bin/env bash
set -euo pipefail

: "${SAM_DEPS_CACHE:?Set SAM_DEPS_CACHE to a writable download directory}"
mkdir -p "${SAM_DEPS_CACHE}"

fetch() {
    local url="$1"
    local filename="$2"
    local expected="$3"
    local destination="${SAM_DEPS_CACHE}/${filename}"
    if [[ ! -f "${destination}" ]]; then
        curl --fail --location --retry 3 --output "${destination}" "${url}"
    fi
    local actual
    if command -v sha512sum >/dev/null 2>&1; then
        # Avoid GNU-only --check/--status flags: macOS runners provide a
        # BSD-style sha512sum with a smaller command-line interface.
        actual="$(sha512sum "${destination}" | awk '{print $1}')"
    else
        actual="$(shasum -a 512 "${destination}" | awk '{print $1}')"
    fi
    [[ "${actual}" == "${expected}" ]]
}

fetch \
    "https://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz" \
    "gsl-2.8.tar.gz" \
    "4427f6ce59dc14eabd6d31ef1fcac1849b4d7357faf48873aef642464ddf21cc9b500d516f08b410f02a2daa9a6ff30220f3995584b0a6ae2f73c522d1abb66b"
fetch \
    "https://www.fftw.org/fftw-3.3.11.tar.gz" \
    "fftw-3.3.11.tar.gz" \
    "ca1bf80490dc6955a0ab49b1af05d6658c2ecc0968b3bde5b4af22271e47d30cd38f6f8347e8e6124091b6a17717447942bee95f94ca29574ba71c4d167af351"
