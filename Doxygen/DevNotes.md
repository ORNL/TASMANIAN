# Developer Notes

Notes and checklists for Tasmanian developers to help procedures such as release checklists and version bumps.

#### Version Bump

When increasing the version, e.g., 7.5 to 7.6:
* update the project directive project() in the main CMakeLists.txt file
* set the `Tasmanian_version_comment` variable to develop, release candidate, or release (empty)
* update the Makefile
* update `Config/AltBuildSystems/TasmanianConfig.hpp`
* update the `setup.py` in `InterfacePython\PipInstaller`

#### Release Checklist

When making a new release:
* test CMake installer on all OS, Linux, Mac and Windows
* test the pip-installer on both Linux and OS (if possible also test on Windows)
* test CUDA, ROCm, OneAPI and MPI on Linux
* check the installation documentation and the currently tested and supported compilers and libraries
* update the change log file

#### After Release Checklist

* update the pip-installer
* update the spack package
    * deprecate one older version for each new version enabled
