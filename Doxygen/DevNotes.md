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
TODO