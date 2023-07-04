============
Installation
============

The pkdgrav3 code requires the following packages:

CMake_ - cmake build system
    Most modern systems already have cmake installed. Pkdgrav3 requires
    version 3.1 or newer of cmake. You can check with "cmake --version"::

        pkdgrav3:~> cmake --version
        cmake version 3.26.0

    .. _CMake: https://cmake.org/

Message Passing Interface (MPI)
  Most distributions or Supercomputer centers provide an MPI implementation.
  If you need to download a version, OpenMPI_ has been well tested.

  .. _OpenMPI: https://www.open-mpi.org

FFTW_
    If FFTW_ is available then two advanced features are enabled in pkdgrav3.
      #. Initial Condition Generation, and,
      #. Power spectrum measurement

    If it is not available on your system you can download it from the FFTW_ website and compile it yourself.

    .. _FFTW: http://www.fftw.org/

    If CMake does not automatically find FFTW then you can define FFTW_ROOT::

        cmake -DFFTW_ROOT=/path/to/fftw

GNU Scientific Library (GSL_)
    This library is usually available on HPC systems, but if not it must be
    downloaded and compiled.

    .. _GSL: https://www.gnu.org/software/gsl/

    pkdgrav3 will locate the GSL installation by invoking gsl-config, so make
    sure that it is in your PATH. Alternatively, you can tell CMake where to
    find it by defining GSL_ROOT_ROOT::

        cmake -DGSL_ROOT_DIR=/opt/gsl/2.5

Python3 - Python Library and Interpreter
    Most modern distributions have Python3 installed. It should be automatically
    detected, but if not refer to the FindPython3_ module of ``CMake``

        .. _FindPython3: https://cmake.org/cmake/help/latest/module/FindPython3.html

    You will need to install the following Python modules:
        * numpy (often installed globally)
        * Cython
        * tomli (if using Python versions before 3.11)

Boost_ C++ Library
    This library is usually available on HPC systems, but if not it must be
    downloaded and compiled.

        .. _Boost: https://www.boost.org

    pkdgrav3 will locate the Boost installation by using the CMake module FindBoost_.

        .. _FindBoost: https://cmake.org/cmake/help/latest/module/FindBoost.html

    Generally setting ``BOOST_ROOT`` is required for user installations.

CUDA_ (optional)
    If your system has a CUDA_ capable GPU then pkdgrav3 can use it.
    The necessary toolkits can be downloaded from nVidia.

        .. _CUDA: https://developer.nvidia.com/cuda-downloads

+++++++++++++++++++++++++++
Python Virtual Environments
+++++++++++++++++++++++++++

There are multiple different tools available to create and manage
Python virtual environments. The following works with the build-in system.

    python -m venv /path/to/pkdgrav3/.venv

When you want to compile or run the code you need to activate the environment.

    source /path/to/pkdgrav3/.venv/bin/activate

The location can be anywhere you want. To install packages use ``pip``.
A requirements.txt is provided that will install the packages necessary
to compile, run and test the code.

    python -m pip install -r requirements.txt

++++++++++++++++++++++++
Compilation Instructions
++++++++++++++++++++++++

The code uses ``CMake`` to detect all dependencies and to compile the code. If the packages are installed correctly (see below),
the one only need run ``cmake`` to configure the build environment, and again to build it::

    cd /path/to/pkdgrav3
    cmake -S . -B build
    cmake --build build

++++++++++++++++++++++
Directory Organization
++++++++++++++++++++++

Both ``CMake`` and ``Python`` allow for a flexible layout of your files. The choice of
how to organise your files is a personal one, but you may consider doing so on a
per-project basis. Consider the case where your source files are in a subdirectory
of your home directors called ``sources`` while your project specific files are
located in ``project``.

+------------------------------+---------------------------------+
| What                         | Location                        |
+==============================+=================================+
| pkdgrav3 source code         | $HOME/sources/pkdgrav3          |
+------------------------------+---------------------------------+
| project1                     | $HOME/project/project1          |
+------------------------------+---------------------------------+

You need to choose a ``build`` directory for pkdgrav3, and virtual environment
directory for Python. Since both of these are in some sense associated with the
project you could put then both in the project directory, for example:

+------------------------------+---------------------------------+
| What                         | Location                        |
+==============================+=================================+
| pkdgrav3 build directory     | $HOME/project/project1/build    |
+------------------------------+---------------------------------+
| virtual environment directory| $HOME/project/project1/.venv    |
+------------------------------+---------------------------------+

This can be easily achieved with::

    cd $HOME/project/project1

    python -m venv .venv
    source .venv/bin/activate
    python -m pip install -r $HOME/sources/pkdgrav3/requirements.txt

    cmake -S $HOME/sources/pkdgrav3 -B build
    cmake --build build
    cmake --install build --prefix $HOME/project/project1

You can then file the pkdgrav3 executable here::

    $HOME/project/project1/bin/pkdgrav3

The ``build`` directory is no longer needed and can be removed.
Obviously if you expect to make changes to ``pkdgrav3`` and recompile
you would leave it intact. It is also not necessary to run the install
phase as you can run ``pkdgrav3`` directly from the build directory,
or copy it somewhere more convenient.

++++++++++++++++++++++++++++++++++++++++++
Swiss National Supercomputer Center (CSCS)
++++++++++++++++++++++++++++++++++++++++++

The necessary libraries can be selected using the modules subsystem with the exception
of some of the Python packages which can be installed using ``pip`` and ``requirements.txt``.

-----------------------
System Specific Modules
-----------------------

Eiger
-----

The GNU programming environment is recommended, thus issue the following commands::

    module swap PrgEnv-cray PrgEnv-gnu
    module load cpeGNU GSL Boost cray-hdf5 cray-fftw CMake cray-python hwloc

Piz Daint
---------

You can compile for the GPU partition and run on the multi-core partition, thus issue the following commands::

    module load daint-gpu cudatoolkit cray-hdf5 GSL Boost cray-fftw CMake cray-python

---------------
Python Packages
---------------

Some additional Python package are required to compile the code. **After** you load the
system modules above, use the following to install the necessary packages::

    cd /path/to/pkdgrav3
    python -m venv --system-site-packages .venv
    source .venv/bin/activate
    python -m pip install -r requirements.txt

This should allow you to build the code and run it in the future without activating the
virtual environment. If you need to recompile, then you must still activate the
virtual environment.

++++++
Ubuntu
++++++

Recent versions of Ubuntu are well supported. The code has been tested on Ubuntu 22.04 (Jammy Jellyfish) as well
as Ubuntu 20.04 (Focal Fossa). It is not recommended to use Ubuntu 18.04 as there is a bug_ with the packaged version
of OpenMPI that was never fixed. If you use Ubuntu prior to 20.04 you need to compile your own MPI (and FFTW).

.. _bug: https://bugs.launchpad.net/ubuntu/+source/openmpi/+bug/1731938

For other versions of Ubuntu (specifically 20.04 and 22.04), the following packages should be sufficient::

    sudo apt update
    sudo apt install -y autoconf automake pkg-config cmake gcc g++ make gfortran git
    sudo apt install -y libfftw3-dev libfftw3-mpi-dev libgsl0-dev libboost-all-dev libhdf5-dev libmemkind-dev libhwloc-dev
    sudo apt install -y python3-dev cython3 python3-pip python3-numpy python3-ddt python3-nose python3-tomli

If you intend to run the test suite, you need the ``xmlrunner`` package. On Ubuntu 22.04 you can install this with::

    sudo apt install -y python3-xmlrunner

For Ubuntu 20.03 you need to create a virtual environment and install it with ``pip``::

    pip3 install xmlrunner

+++++
MacOS
+++++

You will need Xcode_ which can be install from the App Store. You need to launch it at least once.

.. _Xcode: https://apps.apple.com/us/app/xcode/

The following instructions assume that you have installed and are using
Homebrew_ as your package manager.

.. _Homebrew: https://brew.sh

Then install the following packages::

    brew install cmake boost fftw git gsl open-mpi hdf5-mpi python pyenv pyenv-virtualenv

