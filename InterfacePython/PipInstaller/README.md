# Python PIP

The latest stable release of Tasmanian is included in Pip:
[https://pypi.org/project/Tasmanian/](https://pypi.org/project/Tasmanian/)

This folder contains the setup files needed to build the Tasmanian pip archive,
i.e., the source distribution `sdist` package.

### Quick Install

Tasmanian supports installs under Linux, Max OSX, and Windows using Python 3 and
either user or virtual environment installation.

* if using an pip before 1.10 then must first install the build dependencies:
```
python3 -m pip install scikit-build cmake packaging numpy
```
* then perform either a user install
```
python3 -m pip install Tasmanian --user
```
* or virtual environment install
```
python3 -m pip install Tasmanian
```
If the installer is run more than once, the `--no-cache-dir` can be helpful to force repeating
the CMake configuration process and thus clean any errors or obsolete settings in cache.

### Requirements

The following Python modules are required:
```
    setuptools     version compatible with scikit-build
    wheel          version compatible with scikit-build
    packaging      version containing packaging.version.LegacyVersion
    scikit-build   version >= 0.10.0
    cmake          version >= 3.10 (see the top of CMakeLists.txt)
    numpy          version >= 1.10
```
The `scikit-build` package is responsible for turning the CMake build system
into one compatible with `pip`. The rest of the dependencies are there
to facilitate the build process and only `numpy` is required at runtime.
However, `pip` cannot automatically install build dependencies and hence
these have to be installed beforehand.

The process will require a compatible C++ compiler, make sure to have installed
* Linux: `build-essential` with `g++` or the equivalent for the corresponding distribution
* MacOSX: `xcode` or `gcc`
* Windows: MS Visual Studio

### Quick Build Offline

Commands to build the package:
```
cd <repo root, e.g., TASMANIAN>
bash ./InterfacePython/PipInstaller/make_tarball.sh
python3 setup.py sdist
```
At this point, the `dist/Tasmanian-<version>.tar.gz` package file will be created.
xcode
The package can be installed off-line, e.g.,
```
python3 -m pip install dist/Tasmanian-<version>.tar.gz --user
```

### Dependencies

* According to [PEP518](https://www.python.org/dev/peps/pep-0518/), the `pyproject.toml` file must specify the dependencies that are needed to run the `setup.py` script, e.g., scikit-build and packaging.
* The `setup_requires` list in the `setup.py` script specifies the dependencies for the build process, e.g., make.
* The `install_requires` list in the `setup.py` script specifies the run-time dependencies, e.g., numpy.

### Options

Pip accepts both `--global-option=` and `--install-option=` parameters and scikit-build will automatically translate some of those into corresponding CMake commands. However, some python distributions (e.g., Ubuntu) deliberately ignore the options, not sure about the reason behind that, perhaps to limit breaks due to incorrect user-provided install paths. As an alternative, when using Pip, Tasmanian will accept environmental variables to setup BLAS, CUDA, MAGAM, MPI, and MATLAB options.

### Behind the Curtains

The `pip` installer will execute the following steps:
* un-tar the files into a temporary folder
* build Tasmanian (out-of-source) using CMake
    * sci-kit will attempt to download and/or compile missing dependencies, e.g., compilers, CUDA and CMake
* install Tasmanian to another temporary folder
* copy all installed files to <home>/.local or another relevant location
* delete the temporary folders

Currently Tasmanian has three shared libraries and they find each other through `rpath`,
i.e., without the need to be added to the global library path.
This necessitates the usage of `Tasmanian_final_install_path` which is usually
the same as `CMAKE_INSTALL_PREFIX` except when working with `pip`,
in which case the `rpath` is set to the final install path, e.g., `<home>/.local`.

The `setup.py` file sets the CMake options and predicts the install path.
General information is also set here.

The `MANIFEST.in` describes the source files that need to be added to
the package tarball. Some files, such as `.gitignore`, `install` and `Makefile`
are simply not needed by `pip` and are omitted deliberately.
