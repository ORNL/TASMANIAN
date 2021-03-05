#!/usr/bin/env bash
# source this file to execute tasgrid after install
# also sets the python path
export PATH=$PATH:"@Tasmanian_final_install_path@"/bin/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"@Tasmanian_final_install_path@"/lib/
if [[ "@Tasmanian_ENABLE_HIP@" == "ON" ]]; then
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"@Tasmanian_hipomp_rpath@"
fi
if [[ "@Tasmanian_ENABLE_PYTHON@" == "ON" ]]; then
    export PYTHONPATH=$PYTHONPATH:"@Tasmanian_PYTHONPATH@"
fi

# export old and new style cmake search paths
export Tasmanian_DIR="@Tasmanian_final_install_path@"
export Tasmanian_ROOT="@Tasmanian_final_install_path@"
