FROM ubuntu as builder
MAINTAINER Douglas Potter "douglas.potter@uzh.ch"

RUN apt-get update
RUN apt-get install -y openmpi-bin libopenmpi-dev libgsl0-dev autoconf libfftw3-dev libfftw3-mpi-dev python3-pandas git cmake
RUN git clone --depth=1 https://bitbucket.org/dpotter/pkdgrav3.git
RUN cd tmp && cmake /pkdgrav3 && make install

FROM ubuntu
ENTRYPOINT ["/usr/local/bin/pkdgrav3"]
RUN apt-get update
RUN apt-get install -y openmpi-bin libgsl2 libfftw3-mpi3
COPY --from=builder /usr/local/bin/pkdgrav3 /usr/local/bin/pkdgrav3
