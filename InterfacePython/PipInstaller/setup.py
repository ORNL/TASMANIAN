
import sys, site

# do standard skbuild setup
from packaging.version import LegacyVersion
from skbuild.exceptions import SKBuildError
from skbuild.cmaker import get_cmake_version
from skbuild import setup  # This line replaces 'from setuptools import setup'

# Add CMake as a build requirement if cmake is not installed or too old
setup_requires = []
try:
    if LegacyVersion(get_cmake_version()) < LegacyVersion("3.10"):
        setup_requires.append('cmake>=3.10')
except SKBuildError:
    setup_requires.append('cmake>=3.10')
setup_requires.append('numpy>=1.10')

with open('README.md', 'r') as fh:
     readme_file = fh.readlines()

long_description = ""
for line in readme_file[3:]:
    if line.rstrip() == "Quick Install":
        break
    else:
        long_description += line

long_description += "### Quick Install\n Tasmanian supports `--user` and venv install only, see the on-line documentation for details.\n"

# find out whether this is a virtual environment, real_prefix is an older test, base_refix is the newer one
if hasattr(sys, 'real_prefix') or (hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix):
    final_install_path = sys.prefix # sys.prefix points to the virtual environment root
    isvirtual = True
else:
    isvirtual = False
    try:
        final_install_path = site.getuserbase()
    except:
        import os
        # some implementations do not provide compatible 'site' package, assume default Linux behavior
        final_install_path = os.getenv('HOME') + "/.local/"

# check if using OSX Framework environment
isosxframework = False
if sys.platform == 'darwin':
    try:
        if 'python/site-packages' in site.getusersitepackages():
            # appears to be Mac Framework using Library/Python/X.Y/lib/python/site-packages
            isosxframework = True
    except:
        # cannot determine if using Mac Framework
        pass

# setup cmake arguments
cmake_args=[
        '-DCMAKE_BUILD_TYPE=Release',
        '-DBUILD_SHARED_LIBS=ON',
        '-DTasmanian_ENABLE_RECOMMENDED:BOOL=ON',
        '-DTasmanian_ENABLE_PYTHON:BOOL=ON',
        '-DPYTHON_EXECUTABLE:PATH={0:1s}'.format(sys.executable),
        '-DTasmanian_python_pip_final:PATH={0:1s}/'.format(final_install_path),
        ]
if isvirtual:
    cmake_args.append('-DTasmanian_windows_virtual:BOOL=ON')
if isosxframework:
    cmake_args.append('-DTasmanian_osx_framework:BOOL=ON')


# call the actual package setup command
setup(
    name='Tasmanian',
    version='7.7rc1',
    author='Miroslav Stoyanov',
    author_email='stoyanovmk@ornl.gov',
    description='UQ library for sparse grids and Bayesian inference',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://tasmanian.ornl.gov',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: C++',
        'Operating System :: OS Independent',
    ],
    install_requires=['numpy>=1.10'],
    ### cmake portion of the setup, specific to skbuild ###
    setup_requires=setup_requires,
    cmake_args=cmake_args
)
