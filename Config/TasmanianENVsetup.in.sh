#!/usr/bin/env bash
# source this file to execute tasgrid and tasdream after install
# also sets the python path
export PATH=$PATH:@CMAKE_INSTALL_PREFIX@/bin/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:@CMAKE_INSTALL_PREFIX@/lib/
if [[ "@Tasmanian_ENABLE_PYTHON@" == "ON" ]]; then
    export PYTHONPATH=$PYTHONPATH:@Tasmanian_PYTHONPATH@
fi
