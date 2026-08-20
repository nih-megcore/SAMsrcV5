# Test plan

The suite separates fast build/unit coverage from a real-data integration
pipeline. AFNI-dependent work always runs in an image derived from
`ghcr.io/jstout211/afni:latest`; the local runner prefers Podman and GitHub
Actions uses Docker/BuildKit.

## Test data

`fetch-test-data.sh` performs a partial, sparse clone of
`nih-megcore/TEST_ctf_data` at commit
`1d08e47c586fca21163c4e7362d409e62b1c1943`. Before reuse it verifies these Git
objects:

| Path | Git object |
| --- | --- |
| `MRI/ABABABAB_refaced+orig.HEAD` | `ae97b40cf210985ca935d5532415e171a37649a0` |
| `MRI/ABABABAB_refaced+orig.BRIK.gz` | `eb18f9a7015e9095d67865b7ca74852125865280` |
| `20010101/ABABABAB_airpuff_20010101_001.ds` | `d497adbd07caf94c3eb0824b6667d2e32da1c687` |

## Unit and CLI layer

`make test-unit` builds the package and checks:

- Cartesian/spherical coordinate and FIELD/GSL conversions;
- rotation and orthogonal-plane invariants;
- demeaning, detrending, power, Hanning, and Kendall calculations;
- pseudoinverse and SAM beamformer solutions with known matrices;
- BIDS landmark validation, affine-to-AFNI coordinate conversion, output
  protection, and fiducial HEAD attributes;
- the expected installed command set, help output, argument rejection, and
  `1dstats` output for a known sample.

## Installed-wheel layer

`test/wheel` checks distribution metadata, package data, the complete console
script set, native help/error behavior, `1dstats`, and the Python replacement
for `3dNormalize`. Cibuildwheel installs each repaired wheel into a clean test
environment and runs this layer automatically. For a local wheel, run it from
outside the checkout so the installed package cannot be shadowed:

```sh
cd /tmp
/path/to/venv/bin/python -m pytest \
    /path/to/SAM2MULTI/test/wheel -q --import-mode=importlib
```

## AFNI and CTF integration layer

`make test-integration` runs the following pipeline in the AFNI test image:

1. Verify MRI dimensions `208×256×256×1`, signed voxel spacing
   `[-1,-1,1]`, origin, LPI orientation, ORIG space, short datatype, range
   `0..3697`, and fiducial labels.
2. Run `3dNormalize -z`; require nonzero-voxel mean approximately zero and
   standard deviation approximately one.
3. Run `orthohull -i0`; require the aligned/mask datasets, AFNI surface files,
   a 9,002-vertex finite hull with unit normals, and multisphere outputs.
4. Parse the CTF dataset through SAMLIB; require one epoch, 272 primary
   channels, 18,000 samples at 300 Hz, and 120 markers (103 `stim`, 17
   `missingstim`).
5. Run `sam_cov`; validate binary headers, channel count, 5–70 Hz band,
   symmetry, finite values, positive diagonals, covariance types, and exact
   segment/sample counts for Global, Orient, Sum, and stim matrices.
6. Run `sam_wts` with the generated hull and an eighth-order Nolte model over
   a 64-voxel ROI. Require a `4×4×4×272` NIfTI weight volume and 64 finite
   condition-number values.
7. Repeat covariance, weights, `sam_epi`, `sam_ers`, and `sam_3d` with separate
   SAM input and output roots. Require a finite 3D epilepsy image, the `cmdline`
   naming fallback for command-line-only workflows, a 3D+time ERS NIfTI image,
   companion noise files, reconstructed run parameters, and no dataset-local
   `SAM` directory.
8. Run `sam_3d`; require `4×4×4×1` IRP mean and variance NIfTI images with
   `[10,10,-10]` mm spacing. Aggregate nonzero statistics must remain within
   5% of the following baselines:

| Image | Mean | Standard deviation |
| --- | ---: | ---: |
| Mean | 4.30361 | 1.16041 |
| Variance | 2.76794 | 1.73154 |

The exact expectations are stored in `fixtures/expected.json`.

## Optional full-brain 5 mm test

`make test-slow` runs a local, opt-in extension of the real-data pipeline. It
generates the MRI hull, computes the hull's complete X/Y/Z bounds, and runs
`sam_wts` on a 5 mm (`ImageStep 0.5` cm) grid using all 272 MEG channels and an
eighth-order Nolte model. It then runs `sam_3d` for the `stim` marker and
requires finite, nonzero mean and variance images with the same full-volume
dimensions and 5 mm voxel spacing.

The generated NIfTI images, hull, numerical summary, logs, and JUnit report are
written under `.test-results/slow/` (with logs and the report rooted in
`.test-results/`). This test is not part of `make test` or GitHub Actions
because it is substantially more computationally expensive than the bounded
integration ROI.

## Automation and outputs

`.github/workflows/test.yml` runs both layers on pushes, pull requests, and
manual dispatch. It caches the AFNI image layers, pinned fixture checkout, and
host/container compiler results.

Integration logs, JUnit XML, `pipeline-summary.json`, the hull, and the two
representative NIfTI outputs are written to `.test-results/` locally and
uploaded as the `afni-integration-results` Actions artifact. Unit pytest output
is uploaded separately as `unit-test-report`.
