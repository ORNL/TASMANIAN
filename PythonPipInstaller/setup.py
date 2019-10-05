# First catch Tasmanian specific options so CUDA and manual selection
# of the BLAS libraries would be possible
import sys
enable_cuda = False
cuda_path = ""
blas_libs = ""
for opt in sys.argv:
    if opt.startswith("-cuda"):
        # remove the option to avoid confusion with the standard options
        sys.argv.remove(opt)
        enable_cuda = True
        cuda_path = opt.split("=")[1]
    elif opt.startswith("-blas"):
        sys.argv.remove(opt)
        blas_libs = opt.split("=")[1]


# do standard skbuild setup
from packaging.version import LegacyVersion
from skbuild.exceptions import SKBuildError
from skbuild.cmaker import get_cmake_version

import os

# Add CMake as a build requirement if cmake is not installed or too old
setup_requires = []
try:
    if LegacyVersion(get_cmake_version()) < LegacyVersion("3.10"):
        setup_requires.append('cmake')
except SKBuildError:
    setup_requires.append('cmake')


from skbuild import setup  # This line replaces 'from setuptools import setup'


with open('README.md', 'r') as fh:
     readme_file = fh.readlines()

long_description = ""
for line in readme_file:
    if line.rstrip() == "Quick Install":
        break
    else:
        long_description += line


final_install_path = os.getenv('VIRTUAL_ENV') if os.getenv('VIRTUAL_ENV') is not None else os.getenv('HOME') + "/.local/"

cmake_args=[
        '-DCMAKE_BUILD_TYPE=Release',
        '-DBUILD_SHARED_LIBS=ON',
        '-DTasmanian_ENABLE_RECOMMENDED:BOOL=ON',
        '-DPYTHON_EXECUTABLE:PATH={0:1s}'.format(sys.executable),
        '-DTasmanian_python_pip_final:PATH={0:1s}/'.format(final_install_path),
        '-DTasmanian_cinfo:STRING="{0:1s}"'.format('None'),
        '-DTasmanian_ENABLE_CUDA:BOOL={0:1s}'.format('ON' if enable_cuda else 'OFF'),
        ]
if cuda_path != "":
    cmake_args.append('-DCMAKE_CUDA_COMPILER:PATH={0:1s}'.format(cuda_path))
if blas_libs != "":
    cmake_args.append('-DBLAS_LIBRARIES={0:1s}'.format(blas_libs))

setup(
    name='Tasmanian',
    version='7.0',
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
    ### cmake portion of the setup, specific to skbuild ###
    setup_requires=setup_requires,
    cmake_args=cmake_args
)
