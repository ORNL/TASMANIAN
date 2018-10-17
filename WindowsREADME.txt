Starting with Tasmanian 6.0, CMake is the recommended way of building
Tasmanian under MS Windows. The old WindowsMake script is no longer
supported. To use CMake:

1. Launch the CMake GUI
2. Select source and build folders
3. Hit configure
4. Adjust the options as desired
5. Hit configure and generate

After the CMake project is configured and generated, Tasmanian can be
build from the command prompt (cmd) using the following commands:
```
cd <cmake-binary-directory>
cmake --build . --config Release
ctest -C Release
cmake --build . --config Release --target install
```
Both Release and Debug build types are supported.
