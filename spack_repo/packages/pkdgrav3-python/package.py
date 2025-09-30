from spack.package import *
from llnl.util import tty
import os
import sys
import subprocess

# Versions of pkdgrav3 are found in a separate file
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from pkdgrav3_common.pkdgrav3_versions import pkdgrav3_versions

class Pkdgrav3Python(Package):
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

    depends_on('python@3.12', type=('build', 'run'))

    def install(self, spec, prefix):

        venv_dir = os.path.join(prefix, '.venv')
        python_version = self.spec['python'].version.up_to(2)  # '3.12'
        python_exe = os.path.join(self.spec['python'].prefix.bin, f'python{python_version}')

        tty.msg(f"[install] Creating runtime venv in {venv_dir}")
        subprocess.check_call([python_exe, '-m', 'venv', '--system-site-packages', venv_dir])

        # Install runtime Python requirements
        source_root = self.stage.source_path
        req_file = os.path.join(source_root, 'requirements.txt')

        if not os.path.isfile(req_file):
            tty.error(f"[install] requirements.txt not found at {req_file}")
            raise FileNotFoundError(f"requirements.txt not found at {req_file}")

        tty.msg(f"[install] Installing Python requirements into runtime venv from {req_file}")
        pip_exe = os.path.join(venv_dir, 'bin', 'pip')
        subprocess.check_call([pip_exe, 'install', '-r', req_file])

    def setup_run_environment(self, env):
        # Set PATH and PYTHONPATH to activate venv behavior on spack load
        venv_dir = os.path.join(self.prefix, '.venv')
        env.prepend_path('PATH', os.path.join(venv_dir, 'bin'))

        python_version = '{0}.{1}'.format(self.spec['python'].version[0], self.spec['python'].version[1])
        site_packages = os.path.join(venv_dir, 'lib', f'python{python_version}', 'site-packages')

        if os.path.isdir(site_packages):
            env.prepend_path('PYTHONPATH', site_packages)
            tty.msg(f"[runenv] Added PYTHONPATH: {site_packages}")
        else:
            tty.warn(f"[runenv] site-packages not found at {site_packages}, PYTHONPATH not set.")
