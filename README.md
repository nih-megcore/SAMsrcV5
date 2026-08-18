# SAMsrcV5
Synthetic Aperture Magnetometry Suite Version 5

This is a public release branch of the SAMsrcV5 suite of tools, used for
source localization of CTF format MEG data.

These tools are designed to work with AFNI.

The main webpage for the SAM suite is here:

https://megcore.nih.gov/index.php?title=Source_Localization_-_SAM

## Installation

SAMsrcV5 5.0.0 is packaged as a platform wheel containing the compiled SAM
programs and their FFTW/GSL runtime code. Installing a wheel does not require a
C compiler, Make, FFTW, or GSL on the user's computer:

```sh
python -m venv .venv
. .venv/bin/activate                # Windows: .venv\Scripts\activate
python -m pip install samsrcv5
```

Until the wheel artifacts are published to PyPI, download the wheel matching
your platform from the `Build wheel artifacts` GitHub Actions run and pass its
filename to `python -m pip install`.

Wheels are built for Linux x86_64/aarch64, macOS x86_64/arm64, and Windows
AMD64. The computational commands (`sam_cov`, `sam_wts`, `sam_3d`, `sam_4d`,
`sam_power`, `sam_simulate`, and `OPMsim`) are native programs. AFNI workflows
such as `3dNormalize` and `orthohull` still require AFNI commands on `PATH`;
`orthohull` also uses qhull. FreeSurfer tools and `ROIbuilder` likewise require
FreeSurfer and PyGObject/GTK respectively.

Source installations are intended for developers and still require a C17
compiler plus FFTW and GSL development files. Set `SAM_DEPS_ROOT` to a common
dependency prefix when they are outside standard search paths:

```sh
SAM_DEPS_ROOT=/path/to/deps python -m pip install .
```

The legacy `make` build remains supported on Unix-like development systems.

## SAM input and output directories

SAM products no longer have to live inside the MEG dataset directory. The
pipeline commands accept independent roots for existing inputs and newly
generated outputs:

```sh
sam_cov -r FILENAME.ds -m analysis.param \
    -o_SAMdir /work/covariances
sam_wts -r FILENAME.ds -m weights.param \
    -i_SAMdir /work/covariances -o_SAMdir /work/weights
sam_3d -r FILENAME.ds -m image.param \
    -i_SAMdir /work/weights -o_SAMdir /work/images
```

Each option names the SAM root itself. If either option is omitted, that side
independently defaults to `FILENAME.ds/SAM`. Output directories and missing
parents are created automatically. `InputSAMDirectory` and
`OutputSAMDirectory` provide the same settings in parameter files, while an
explicit `ImageDirectory` continues to control final image placement.

## Tests

The test suite has a fast host-side unit/CLI layer and an AFNI/CTF integration
layer. The integration tests download only the required MRI and MEG fixture
paths from `nih-megcore/TEST_ctf_data`, pinned to commit
`1d08e47c586fca21163c4e7362d409e62b1c1943`, and verify the Git object IDs
before use.

Run the legacy unit tests on a system with GSL, FFTW, Python, and pytest installed:

```sh
make test-unit
```

Run the complete pipeline locally with Podman (selected automatically when it
is installed):

```sh
make test-integration
```

Run the opt-in full-brain beamformer at 5 mm resolution locally with:

```sh
make test-slow
```

The slow test derives its volume from the MRI brain hull and writes viewable
mean and variance NIfTI images to `.test-results/slow/`. It is deliberately
excluded from `make test` and GitHub Actions.

Set `CONTAINER_ENGINE=docker` to use Docker explicitly. Integration results,
logs, JUnit XML, representative NIfTI files, and the generated hull are written
to `.test-results/`. GitHub Actions runs both layers for every push and pull
request; its AFNI image layers, fixture checkout, and compiler output are
cached between runs. See [test/README.md](test/README.md) for the full test plan
and expected numerical outputs.
