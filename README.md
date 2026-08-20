# SAMsrcV5
Synthetic Aperture Magnetometry Suite Version 5

This is a public release branch of the SAMsrcV5 suite of tools, used for
source localization of CTF format MEG data.

These tools are designed to work with AFNI.

The main webpage for the SAM suite is here:

https://megcore.nih.gov/index.php?title=Source_Localization_-_SAM

## Installation

SAMsrcV5 5.1.0 is packaged as a platform wheel containing the compiled SAM
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
`sam_epi`, `sam_ers`, `sam_power`, `sam_simulate`, and `OPMsim`) are native programs. AFNI workflows
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

## BIDS fiducials to AFNI

Convert a BIDS T1w NIfTI image and the `AnatomicalLandmarkCoordinates` from its
matching JSON sidecar into an AFNI dataset with Nasion, Left Ear, and Right Ear
tags:

```sh
convert_json_fids_to_head sub-01_T1w.nii.gz
convert_json_fids_to_head sub-01_T1w.nii
convert_json_fids_to_head sub-01_T1w.nii.gz --output-dir afni --overwrite
```

The landmark voxel coordinates are transformed with the affine stored in the
input NIfTI before conversion to AFNI's LPS coordinates. The command requires
AFNI's `3dcopy` on `PATH`. The same operation is available from Python:

```python
from samsrcv5.fiducials import convert_json_fids_to_head

brik, head = convert_json_fids_to_head("sub-01_T1w.nii.gz")
```

The packaged `orthohull` and `orthohull.py` commands perform this conversion
automatically for `.nii` and `.nii.gz` inputs. The matching JSON sidecar must
contain valid NAS, LPA, and RPA coordinates. The generated AFNI HEAD/BRIK pair
is kept in a private temporary directory and removed when `orthohull` exits,
including after an error.

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

Parameter files are optional for `sam_cov`, `sam_wts`, `sam_3d`, `sam_4d`,
`sam_ers`, and `sam_power` when all required analysis parameters are supplied
on the command line. Parameters with multiple values retain their parameter-file
spelling as long options:

```sh
sam_3d -r FILENAME.ds \
    --InputSAMDirectory /work/weights \
    --OutputSAMDirectory /work/images \
    --Marker stim -0.1 0.3 TRUE \
    --CovBand 13 35 --ImageBand 13 35 \
    --CovType GLOBAL --ImageMetric Power
```

Without `-m`, `cmdline` is used in place of the parameter-file basename for
input and output names. Use `--OutName`, `--CovName`, or `--WtsName` when the
run should use another name. Command-line values continue to override values
from a supplied parameter file.

## SAM parameter editor

Run `sam_param_gui` to create or edit parameter files for `sam_cov`, `sam_wts`,
`sam_3d`, and `sam_ers`. The Tkinter interface groups the current supported
parameters by function, highlights those used by the selected program, checks
their syntax, and previews the generated plain-text file.

The editor creates `~/samparams` when necessary and saves all `.param` files
there. Existing files may be opened from elsewhere; comments, blank lines, and
parameters unknown to the editor are retained when the managed copy is saved.
Tk must be provided by the host Python installation (often through an operating
system package named `python3-tk`).

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
