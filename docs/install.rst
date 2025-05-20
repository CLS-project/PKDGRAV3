============
Installation
============

The pkdgrav3 code requires the following packages:

CMake_ - cmake build system
    Most modern systems already have cmake installed. Pkdgrav3 requires
    version 3.24 or newer of cmake. You can check with "cmake --version"::

        pkdgrav3:~> cmake --version
        cmake version 3.26.0

    .. _CMake: https://cmake.org/

    The CMake package will be installed in the Python virtul environment (see below).

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

+++++++++++++++++++++
spack package manager
+++++++++++++++++++++

You can use ``spack`` to install the necessary dependencies and compiler pkdgrav3.
To start with ``spack`` you simply clone the repository to a directory of your choice,
and source the supplied environment file.

    git clone -c feature.manyFiles=true --depth=2 https://github.com/spack/spack.git
    source spack/share/spack/setup-env.sh

You then need to add the pkdgrav3 repository.

    spack repo add spack_repo
    ==> Added repo with namespace 'pkdgrav3'.

You can check how the dependencies will be handled with:

    spack spec pkdgrav3

For local development you should build in a spack environment:

    cd /path/to/pkdgrav3
    spack env activate  .
    spack install pkdgrav3



++++++++++++++++++++++++++++++++++++++++++
Swiss National Supercomputer Center (CSCS)
++++++++++++++++++++++++++++++++++++++++++


The Swiss National Supercomputer Center has moved towards their own tool (``uenv``) for
package management. On "Daint" the necessary tool is already installed,
while on "Eiger" it may need to be installed.

----------------
Installing uenv2
----------------

Eiger
-----


To install ``uenv2`` on Eiger first follow the instructions found at `Installing uv`_.
Next follow the instructions for `Installing uenv`_.

  .. _Installing uv: https://docs.astral.sh/uv/getting-started/installation/
  .. _Installing uenv: https://confluence.cscs.ch/spaces/KB/pages/701726829/uenv+user+environments#uenvuserenvironments-Installing


If setup correctly you should see version 8.1.0 (or perhaps later)::

    [eiger]$ uenv --version
    8.1.0

If you use both Eiger and Daint then you should take care to only change uenv if you are on Eiger.
To do this the change to ``~/.bashrc`` should be the following::

    export PATH=$HOME/.local/$(uname -m)/bin:$PATH
    if test ${HOST##eiger} != $HOST ; then
        unset -f uenv
    fi

Piz Daint
---------

On Daint the tool is already installed::

    [daint]$ uenv --version
    8.1.0

-----------
Using uenv2
-----------

To find available user environments use the ``uenv find`` command::


    [daint][dpotter@daint-ln003 ~]$ uenv image find service::pkdgrav3
    uenv                     arch   system  id                size(MB)  date
    pkdgrav3/3.4:1824353160  gh200  daint   ee775978c79587fb   4,044    2025-05-19
    pkdgrav3/3.4:1825038404  gh200  daint   5e0b0ac41d73e4e5   3,866    2025-05-19

Your will need to pull the version you want, for example ``pkdgrav3/3.4``.
If there are multiple images with the same version then you should choose the most recent::

    [daint]$ uenv image pull service::pkdgrav3/3.4:1825038404
    pulling 5e0b0ac41d73e4e5 100.00% ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 3867/3867 (0.00 MB/s)
    updating pkdgrav3/3.4:1825038404@daint%gh200

The same holds for Eiger::

    [eiger]$ uenv image find service::pkdgrav3
    uenv                     arch  system  id                size(MB)  date
    pkdgrav3/3.4:1825268321  zen2  eiger   09cbc44a4026a1f6     643    2025-05-19
    pkdgrav3/3.4:1826427762  zen2  eiger   c3268965cf74074c     643    2025-05-20
    [eiger]$ uenv image pull service::pkdgrav3/3.4:1826427762
    pulling c3268965cf74074c 100.00% ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 644/644 (0.00 MB/s)
    updating pkdgrav3/3.4:1826427762@eiger%zen2

----------------
Running pkdgrav3
----------------

To run ``pkdgrav3`` from the image (the specific version) you can use the ``pkdgrav3`` view::


    uenv start --view=pkdgrav3 pkdgrav3
    pkdgrav3 cosmology.par

If running in the queue then use the ``slurm`` flags::

    #SBATCH --uenv=pkdgrav3 --view=pkdgrav3

-------------------
Developing pkdgrav3
-------------------

If you want to compile your own version of ``pkdgrav3`` then you can still use a user environment
and it will provide the necessary tools. In this case use the ``modules`` view::

    [daint]$ uenv start --view=modules pkdgrav3
    [daint]$ module load boost cmake cray-mpich fftw gcc gsl hdf5 pkdgrav3-python

Note that the ``pkdgrav3-python`` module loads ``python`` and sets up a virtual environment
will the python modules necessary to compile and run ``pkdgrav3``.

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

If you intend to use ``spack`` you can instruct it to use the Homebrew versions::

    spack external find --not-buildable cmake git gmake gsl hdf5 openmpi apple-clang gcc

This will create the file ``~/.spack/packages.yaml``. In addtion you will want to add
the following packages::

      python:
        externals:
        - spec: python@3.12
          prefix: /opt/homebrew/opt/python@3.12
        buildable: false
      fftw:
        externals:
        - spec: fftw@3.3.10
          prefix: /opt/homebrew/Cellar/fftw/3.3.10_2
        buildable: false
      boost:
        externals:
        - spec: boost@1.88.0
          prefix: /opt/homebrew/Cellar/boost/1.88.0
        buildable: false

You can verify their version and location with::

    brew --prefix python@3.12 fftw boost

+++++++++++++++++++++
OpenMPI Based Systems
+++++++++++++++++++++

The code uses a Hybrid MPI/threads mechanism for parallization. Some versions of OpenMPI
will bind each MPI rank to a single core in an effort to be helpful. This will result
in very poor performance. To disable this "feature", define ``OMPI_MCA_hwloc_base_binding_policy``
to be ``none`` when running::

    export OMPI_MCA_hwloc_base_binding_policy=none
    srun ./pkdgrav3 cosmology.par
