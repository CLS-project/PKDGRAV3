# Creates a production image for pkdgrav3. Requires the build image (Dockerfile.build).
# docker image build --build-arg TARGET_ARCHITECTURE=skylake --build-arg BRANCH=develop -t dpotter/pkdgrav3:3.0.4-skylake - <Dockerfile
# docker run -it --rm -v ${PWD}:/pkdgrav3 -u $(id -u ${USER}):$(id -g ${USER}) dpotter/pkdgrav3:3.0.4-skylake cosmology.par
FROM dpotter/pkdgrav3-build as builder
LABEL maintainer="douglas.potter@uzh.ch"
ARG BRANCH=develop
ARG TARGET_ARCHITECTURE=auto
RUN cd /build && git clone --depth=1 --branch ${BRANCH} https://bitbucket.org/dpotter/pkdgrav3.git\
  && cd /tmp && cmake -DTARGET_ARCHITECTURE=${TARGET_ARCHITECTURE} /build/pkdgrav3 && make install

FROM nvidia/cuda:10.1-base-ubuntu18.04
LABEL maintainer="douglas.potter@uzh.ch"
WORKDIR /pkdgrav3
ENTRYPOINT ["/usr/local/bin/pkdgrav3"]
COPY --from=builder /usr/local/bin/pkdgrav3 /usr/local/bin/
COPY --from=builder /usr/local/lib/libmpi.so.12.1.1 /usr/local/lib/
RUN ln -s libmpi.so.12.1.1 /usr/local/lib/libmpi.so.12 && ln -s libmpi.so.12.1.1 /usr/local/lib/libmpi.so && ldconfig
RUN apt-get update && apt-get install -y --no-install-recommends libgsl23 libhdf5-100 libmemkind0 libhwloc5 libpython3.6 && apt-get clean all
