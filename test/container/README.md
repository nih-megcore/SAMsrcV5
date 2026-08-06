# AFNI containers

```
neurodocker generate docker \
    --pkg-manager yum \
    --base-image fedora:43 \
    --afni method=binaries version=latest \
> afni-binaries.Dockerfile

```

```
podman build --tag afni:latest --file afni-binaries.Dockerfile .
```

This is available on ghcr.io/jstout211/afni:latest

`test.Dockerfile` extends that image with SAM2MULTI's build and Python test
dependencies. Normally it should be invoked through `../run-container-tests.sh`:

```sh
CONTAINER_ENGINE=podman ../run-container-tests.sh
```

The runner defaults to Podman when both Podman and Docker are installed. The
GitHub Actions workflow uses Docker/BuildKit so that the large AFNI base and
dependency layers can use GitHub's shared build cache.
