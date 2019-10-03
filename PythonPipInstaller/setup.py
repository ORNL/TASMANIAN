
from packaging.version import LegacyVersion
from skbuild.exceptions import SKBuildError
from skbuild.cmaker import get_cmake_version

import skbuild
import setuptools
import os

# Add CMake as a build requirement if cmake is not installed or too old
setup_requires = []
try:
    if LegacyVersion(get_cmake_version()) < LegacyVersion("3.10"):
        setup_requires.append('cmake')
except SKBuildError:
    setup_requires.append('cmake')

import sys, site
from skbuild import setup  # This line replaces 'from setuptools import setup'

with open('README.md', 'r') as fh:
     readme_file = fh.readlines()

long_description = ""
for line in readme_file:
    if line.rstrip() == "Quick Install":
        break
    else:
        long_description += line

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
    cmake_args=[
        '-DCMAKE_BUILD_TYPE=Release',
        '-DBUILD_SHARED_LIBS=ON',
        '-DTasmanian_ENABLE_RECOMMENDED:BOOL=ON',
        '-DPYTHON_EXECUTABLE:PATH={0:1s}'.format(sys.executable),
        '-DTasmanian_python_pip_final:PATH={0:1s}/.local/'.format(os.getenv("HOME"))
    ],
)
