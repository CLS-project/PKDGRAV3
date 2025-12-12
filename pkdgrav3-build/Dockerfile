# Create the build/test image for pkdgrav3.
# Make sure the required files are present: euclid_flagship_500.hdf5
# docker image build --memory=4g --memory-swap=4g -t dpotter/pkdgrav3-build:latest pkdgrav3-build
FROM nvidia/cuda:12.1.0-devel-ubuntu22.04 AS pkdgrav3-build-mpich
LABEL maintainer="douglas.potter@uzh.ch"

# Base system + build deps
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
        tzdata \
        autoconf automake pkg-config \
        libfftw3-dev libfftw3-mpi-dev \
        gcc g++ gfortran make \
        git wget \
        libgsl0-dev libboost-all-dev \
        libhdf5-dev libmemkind-dev libhwloc-dev \
        python3-dev python3-pip python3-numpy \
        python3-tomli python3-xmlrunner \
        python3-ddt python3-nose python3-h5py \
    && rm -rf /var/lib/apt/lists/*

# Python build tools
RUN python3 -m pip install --no-cache-dir --upgrade \
        Cython cmake dill clang-format==21.1.7

# Prepare build directory
WORKDIR /build

COPY euclid_flagship_500.hdf5 \
     b0-final.std \
     b0-final-np-asym-k1.acc.npy \
     b0-final-p0.10-asym-k1.acc.npy \
     b0-final-ref.grp.npy \
     ./

COPY sph-data ./sph-data
COPY ANEOSmaterial ./ANEOSmaterial
COPY EOSlib ./EOSlib

# MPI config
ENV OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMPI_MCA_btl_vader_single_copy_mechanism=none

ENTRYPOINT ["/bin/bash"]