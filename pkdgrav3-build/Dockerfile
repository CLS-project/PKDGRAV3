# Create the build/test image for pkdgrav3.
# Make sure the required files are present: euclid_flagship_500.hdf5
# docker image build --memory=4g --memory-swap=4g -t dpotter/pkdgrav3-build:latest pkdgrav3-build
FROM nvidia/cuda:12.1.0-devel-ubuntu22.04 as pkdgrav3-build-mpich
LABEL maintainer="douglas.potter@uzh.ch"
RUN apt update && DEBIAN_FRONTEND="noninteractive" apt -y install tzdata && apt install -y autoconf automake libfftw3-dev libfftw3-mpi-dev pkg-config cmake gcc g++ make gfortran wget git libgsl0-dev libboost-all-dev libhdf5-dev libmemkind-dev libhwloc-dev && apt-get clean all
RUN apt-get update && apt-get install -y python3-dev cython3 python3-pip python3-numpy python3-tomli python3-xmlrunner python3-ddt python3-nose
RUN mkdir /build
WORKDIR /build
COPY euclid_flagship_500.hdf5 b0-final.std b0-final-np-asym-k1.acc.npy b0-final-p0.10-asym-k1.acc.npy b0-final-ref.grp.npy ./
COPY sph-data ./sph-data
COPY ANEOSmaterial ./ANEOSmaterial
COPY EOSlib ./EOSlib
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
ENV OMPI_MCA_btl_vader_single_copy_mechanism=none
ENTRYPOINT /bin/bash
