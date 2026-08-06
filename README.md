# SAMsrcV5
Synthetic Aperture Magnetometry Suite Version 5

This is a public release branch of the SAMsrcV5 suite of tools, used for
source localization of CTF format MEG data.

These tools are designed to work with AFNI.

The main webpage for the SAM suite is here:

https://megcore.nih.gov/index.php?title=Source_Localization_-_SAM

## Tests

The test suite has a fast host-side unit/CLI layer and an AFNI/CTF integration
layer. The integration tests download only the required MRI and MEG fixture
paths from `nih-megcore/TEST_ctf_data`, pinned to commit
`1d08e47c586fca21163c4e7362d409e62b1c1943`, and verify the Git object IDs
before use.

Run the unit tests on a system with GSL, FFTW, Python, and pytest installed:

```sh
make test-unit
```

Run the complete pipeline locally with Podman (selected automatically when it
is installed):

```sh
make test-integration
```

Set `CONTAINER_ENGINE=docker` to use Docker explicitly. Integration results,
logs, JUnit XML, representative NIfTI files, and the generated hull are written
to `.test-results/`. GitHub Actions runs both layers for every push and pull
request; its AFNI image layers, fixture checkout, and compiler output are
cached between runs. See [test/README.md](test/README.md) for the full test plan
and expected numerical outputs.
