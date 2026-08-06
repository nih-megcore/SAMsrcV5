ARG AFNI_IMAGE=ghcr.io/jstout211/afni:latest
FROM ${AFNI_IMAGE}

RUN dnf install -y \
        ccache \
        fftw-devel \
        findutils \
        git \
        gsl-devel \
        gzip \
        python3-numpy \
        python3-pytest \
    && dnf clean all \
    && rm -rf /var/cache/dnf

ENV PATH="/work/bin:/opt/afni-latest:${PATH}" \
    OMP_NUM_THREADS=2

WORKDIR /work

COPY . /work

RUN cd config && ./configure
