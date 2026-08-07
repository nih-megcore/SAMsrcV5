# Third-party notices

SAMsrcV5 wheels contain statically linked copies of the following libraries:

- GNU Scientific Library (GSL) 2.8, licensed under GPL-3.0-or-later.
  Source: <https://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz>
- FFTW 3.3.11, licensed under GPL-2.0-or-later.
  Source: <https://www.fftw.org/fftw-3.3.11.tar.gz>

Linux and Windows wheels may contain GCC OpenMP/runtime components, licensed
under GPL-3.0-or-later with the GCC Runtime Library Exception 3.1. The NIFTI-1
header in `include/nifti1.h` is public domain. Several CTF-format compatibility
headers retain their original CTF Systems or Biomagnetic Technologies
copyright notices.

The build workflow publishes the exact dependency source archives and build
scripts beside binary wheel artifacts as corresponding source.
