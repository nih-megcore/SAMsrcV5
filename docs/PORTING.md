# Portable packaging design

SAMsrcV5 now uses `pyproject.toml`, scikit-build-core, and CMake to produce a
platform wheel. Compilation happens in GitHub Actions; installing that wheel
with pip only unpacks Python files, data, launchers, and prebuilt executables.

| Component | Original project | Portable package | Action |
| --- | --- | --- | --- |
| Python tools | Scripts installed from `Mains/` by Make | `samsrcv5._legacy` plus console entry points | Preserve command names and resolve data relative to the installed package |
| Native programs | Make targets in `Mains/` | Executables in `samsrcv5/_bin` | Build explicitly listed sources with CMake and invoke them through Python launchers |
| Native libraries | Host GSL and FFTW installations | GSL 2.8 and FFTW 3.3.11 statically linked in CI | Download pinned sources, verify SHA-512, and retain corresponding source artifacts |
| Package data | Make copy/install rules | `samsrcv5/data` | Install AFNI templates, GTK UI, and atlas data explicitly |
| CLI surface | Programs and `.py` scripts on `PATH` | `[project.scripts]` | Preserve the supported Make-default command set and common `.py` aliases |
| Tests | Make unit tests and AFNI/CTF pipeline | Existing tests plus installed-wheel tests | Test repaired wheels outside the checkout on every target platform |

The wheel matrix covers Linux x86_64 and aarch64, macOS x86_64 and arm64, and
Windows AMD64. The native payload has no CPython ABI dependency, so one
`py3-none-<platform>` wheel per platform supports all declared Python versions
(3.10 and newer).

GSL and FFTW are GPL libraries. Their code is statically included in the wheel,
so the binary distribution is licensed as GPL-3.0-or-later; the wheel also
contains the project license, both dependency license texts, and
`THIRD_PARTY_NOTICES.md`. Each CI run uploads the project sdist together with
the exact GSL and FFTW source archives used for the build.

AFNI, FreeSurfer, qhull, GTK3, and PyGObject are workflow/runtime tools rather
than link-time dependencies. They remain external and are detected when the
corresponding commands are used.
