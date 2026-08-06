# Use neurdocker to containerize

```
neurodocker generate docker \
    --pkg-manager yum \
    --base-image fedora:43 \
    --afni method=binaries version=latest \
> afni-binaries.Dockerfile

```

```
docker build --tag afni:latest --file afni-binaries.Dockerfile .
```

This is available on ghcr.io/jstout211/afni:latest
