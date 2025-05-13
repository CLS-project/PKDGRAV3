from spack.package import *
from llnl.util import tty
import os
import sys
import subprocess

# Versions of pkdgrav3 are found in a separate file
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from pkdgrav3_common.pkdgrav3_versions import pkdgrav3_versions

class Pkdgrav3(CMakePackage):
    """PKDGRAV3: Cosmological N-body code for large-scale simulations."""

    homepage = "https://dpotter@bitbucket.org/dpotter/pkdgrav3"
    git      = "https://dpotter@bitbucket.org/dpotter/pkdgrav3.git"

    maintainers = ['dpotter']

    # Register all of the versions
    for v, info in pkdgrav3_versions.items():
        kwargs = {}
        if 'tag' in info:
            kwargs['tag'] = info['tag']
        if 'branch' in info:
            kwargs['branch'] = info['branch']
        if info.get('preferred', False):
            kwargs['preferred'] = True
        version(v, **kwargs)

    # We also want the matching python
    for v in pkdgrav3_versions:
        depends_on('pkdgrav3-python@{0}'.format(v), when='@{0}'.format(v), type=('build', 'run'))

    variant('build_type', default='RelWithDebInfo',
            description='CMake build type',
            values=('Debug', 'Release', 'RelWithDebInfo', 'MinSizeRel'))

    variant('cuda', default=False, description='Enable CUDA support')
    variant('cuda_arch', default='80', description='CUDA architecture',
            values=('60', '70', '75', '80', '86', '90'), multi=False)

    depends_on('cmake@3.15:', type='build')
    depends_on('mpi')
    depends_on('fftw')
    depends_on('gsl')
    depends_on('hdf5+hl')
    depends_on('git')
    depends_on('boost+thread')
    depends_on('cuda', when='+cuda')

    def cmake_args(self):
        boost_spec = self.spec['boost']
        if not boost_spec.variants['thread'].value:
            raise InstallError('This package requires boost+thread')
        spec = self.spec
        python_prefix = self.spec['pkdgrav3-python'].prefix
        venv_python = os.path.join(python_prefix, '.venv', 'bin', 'python')
        args = []

        args.append(self.define('Python3_ROOT_DIR', python_prefix))
        args.append(self.define('Python3_EXECUTABLE', venv_python))
        args.append(self.define('Python_EXECUTABLE', venv_python))
        args.append(self.define('PYTHON_EXECUTABLE', venv_python))
        args.append(self.define('PKDGRAV3_USE_MPI', True))
        args.append(self.define('MPI_C_COMPILER', spec['mpi'].mpicc))
        args.append(self.define('MPI_CXX_COMPILER', spec['mpi'].mpicxx))
        args.append(self.define('FFTW_ROOT', spec['fftw'].prefix))
        args.append(self.define('GSL_ROOT_DIR', spec['gsl'].prefix))
        args.append(self.define('HDF5_ROOT', spec['hdf5'].prefix))
        args.append(self.define('Boost_ROOT', spec['boost'].prefix))

        if '+cuda' in spec:
            args.append(self.define('PKDGRAV3_USE_CUDA', True))
            cuda_arch = spec.variants['cuda_arch'].value
            args.append(self.define('CMAKE_CUDA_ARCHITECTURES', cuda_arch))
        else:
            args.append(self.define('PKDGRAV3_USE_CUDA', False))

        return args
