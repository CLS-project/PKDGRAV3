# Create the build/test image for pkdgrav3.
# Make sure the required files are present: euclid_flagship_500.hdf5
# docker image build --memory=1g --memory-swap=1g -t dpotter/pkdgrav3-build:latest pkdgrav3-build
FROM nvidia/cuda:10.1-devel-ubuntu18.04 as pkdgrav3-build-mpich
LABEL maintainer="douglas.potter@uzh.ch"
ARG MPICH_VERSION=3.1.4
ARG FFTW_VERSION=3.3.8
RUN apt-get update && apt-get install -y autoconf automake pkg-config cmake gcc g++ make gfortran wget git libgsl0-dev libhdf5-dev libmemkind-dev libhwloc-dev && apt-get clean all
RUN mkdir /build
RUN cd /build && wget http://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz\
  && tar xvzf mpich-${MPICH_VERSION}.tar.gz && cd /build/mpich-${MPICH_VERSION} \
  && ./configure && make && make install && make clean && rm -rf /build/mpich-${MPICH_VERSION}.tar.gz /build/mpich-${MPICH_VERSION}
RUN cd /build && wget http://www.fftw.org/fftw-${FFTW_VERSION}.tar.gz\
  && tar xvzf fftw-${FFTW_VERSION}.tar.gz && cd /build/fftw-${FFTW_VERSION} \
  && ./configure --enable-float --enable-mpi --enable-threads\
  && make && make install && make clean && rm -rf /build/fftw-${FFTW_VERSION}.tar.gz /build/fftw-${FFTW_VERSION}
RUN apt-get install -y python3-dev python3-pip python3-numpy && apt-get clean all && pip3 install xmlrunner ddt
WORKDIR /build
COPY euclid_flagship_500.hdf5 b0-final.std b0-final-np-asym-k1.acc.npy b0-final-p0.10-asym-k1.acc.npy ./
ENTRYPOINT /bin/bash
