# pkdgrav3

Douglas Potter <douglas.potter@uzh.ch>  
Joachim Stadel <stadel@physik.uzh.ch>  

The User Guide for pkdgrav3 can be found at: [https://pkdgrav3.readthedocs.io/](https://pkdgrav3.readthedocs.io/)

## Quick Start

The pkdgrav3 code uses the "cmake" build system. It is recommended to
use an "out of source" build. The easiest way to accomplish this is to
create a subdirectory in the pkdgrav3 source directory:

```
cd /path/to/pkdgrav3
cmake -S . -B build
cmake --build build
```

This will build a single executable ``pkdgrav3`` which will be found in the ``build`` directory as well as other utility programs.

## Prerequisites

Detailed instructions for building pkdgrav3, as well as advice for specific systems can be found in the
[online installation guide](https://pkdgrav3.readthedocs.io/en/latest/install.html).

### CMake - cmake build system

Most modern systems already have cmake installed.
Pkdgrav3 requires version 3.14 or newer of cmake.
You can check with:

```
pkdgrav3:~> cmake --version
cmake version 3.27.8
```

If you need a more recent version is can be found at: [https://cmake.org/](https://cmake.org/)

### Python3 - Python Library and Interpreter

Most modern distributions have Python3 installed.
It should be automatically detected, but if not refer to: [https://cmake.org/cmake/help/latest/module/FindPython3.html](https://cmake.org/cmake/help/latest/module/FindPython3.html)

Versions of Python that are "end-of-life" are not necessarily supported, and in particular versions before 3.8 will not work. Support status of Python releases can be checked here: [https://devguide.python.org/versions/](https://devguide.python.org/versions/)

### GSL - The GNU Scientific Library
This library is usually available on HPC systems, but if not it must be downloaded and compiled, and can be found at this URL: [https://www.gnu.org/software/gsl/](https://www.gnu.org/software/gsl/)

pkdgrav3 will locate the GSL installation by invoking gsl-config, so make sure that it is in your PATH.
Alternatively, you can tell CMake where to find it by defining ``GSL_ROOT_ROOT``:

```
cmake -DGSL_ROOT_DIR=/opt/gsl/2.5
```

### Boost - C++ Template Lirary
This library is usually available on HPC systems, but if not it must be downloaded and compiled, and can be found at this URL: [https://www.boost.org](https://www.boost.org)

pkdgrav3 will locate the Boost installation by using the CMake module: [https://cmake.org/cmake/help/latest/module/FindBoost.html](https://cmake.org/cmake/help/latest/module/FindBoost.html)

Generally setting ``BOOST_ROOT`` is required for user installations:

```
cmake -DBOOST_ROOT=/path/to/boost ...
```

### FFTW3 - Fast Fourier Transform Library

If FFTW3 is is not available on your system it can be obtained from: [http://www.fftw.org/](http://www.fftw.org/)

If CMake does not automatically find FFTW then you can define ``FFTW_ROOT``:

```
cmake -DFFTW_ROOT=/path/to/fftw ...
```

If you build it yourself you must enable support for float and MPI transforms.
You will need to have a working MPI and provide ``--enable-mpi`` to the ``configure`` script to compile and build MPI support.
Additionally you will need to provide ``--enable-single`` to the ``configure`` script to support single precision transforms. It is common to build both libraries by running configure/make twice, with and without ``--enable-single``.

### CUDA (optional)

If your system has a CUDA capable GPU then pkdgrav3 can use it.
The necessary toolkits can be downloaded from nVidia:
[https://developer.nvidia.com/cuda-downloads](https://developer.nvidia.com/cuda-downloads)

## Configuration

### Potentials in Lightcone particle output


The potential for each particle can be output for the lightcone.
Enable by setting ``POTENTIAL_IN_LIGHTCONE`` when running cmake:

```
cmake -DPOTENTIAL_IN_LIGHTCONE=True ...
```

## Build

Once CMake has been run to produce a Makefile and associated files,
the "make" command is used to build the program, as in:

```
cmake --build build
```

The build can be done in parallel so if you are on, for example,
a 16 core machine, the build process can be sped up with::

```
cmake --build -j 16 build
```

## Running

This version is run using the MPI system on the cluster in question.
Normally this involves a special command (often "mpirun" or "mpiexec"),
for example:

```
mpiexec pkdgrav3 simfile.par
```

Consult your cluster documentation on how to run MPI programs.

## IA annotations

For the CUDA compiler, the GNU compiler must be older than 6.0
This is the command used for preparing and compiling my version, on virus:

```
mkdir build
cd build
CC=/scratch/isaacaa/opt/gcc53/bin/gcc CXX=/scratch/isaacaa/opt/gcc53/bin/g++ cmake -DFFTW_ROOT=/scratch/isaacaa/opt/fftw3  ..
make
```

The newly added compile-time flags are described in README.hydro
