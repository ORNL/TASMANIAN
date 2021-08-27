# libEnsemble Integration

Tasmanian has integration with the [libEnsemble](https://github.com/Libensemble/libensemble) project, which is a Python library to coordinate the concurrent evaluation of dynamic ensembles of calculations. The workflow management system of libEnsemble can be wrapped around the Tasmanian sparse grid methods for constructing surrogate models and libEnsemble can serve as an intermediary between Tasmanian and the model of interest.

The functionality provided by libEnsemble is similar to that of `TasGrid::constructSurrogate()` and `TasGrid::mpiConstructSurrogate()`. While the native Tasmanian methods are functional, libEnsemble provides much more sophisticated workflow algorithms and much wider support for hardware architectures. One of the main advantages of libEnsemble is that the code is portable on shared memory multi-threaded environments, as well as numerous distributed cluster systems.


### Requirements

Using Tasmanian and libEnsemble together requires:

* Tasmanian 7.0 or newer
* Tasmanian Python bindings are enabled, e.g., with `Tasmanian_ENABLE_PYTHON=ON` or the pip-installer
* libEnsemble from the latest development branch
* see the [libEnsemble](https://github.com/Libensemble/libensemble) documentation for any additional requirements


### Installation

The easiest way to install the latest development version of libEnsemble is to clone the repo and use a local pip-installer:
```
    git clone https://github.com/Libensemble/libensemble.git
    cd libensemble
    git checkout develop
    python3 setup.py sdist
    python3 -m pip uninstall libensemble
    python3 -m pip install dist/*.tar.gz
```
See the [Installation](Doxygen/Installation.md) instructions for Tasmanian.

Notes:
* the `pip` approach required both `git` and the Python `pip` modules to be installed
* see also the libEnsemble instructions for more details
* both Tasmanian and libEnsemble works well within a Python `venv` environment


### Basic Usage

An example of the Tasmanian-libEnsemble integration can be found in the [libEnsemble testing suite](https://github.com/Libensemble/libensemble/blob/develop/libensemble/tests/regression_tests/test_persistent_tasmanian.py). The Tasmanian sparse grid construction algorithm is implemented as the `sparse_grid_batched` libEnsemble native `gen_f` function:
```
from libensemble.gen_funcs.persistent_tasmanian import sparse_grid_batched as gen_f_batched
```
The construction method accepts the following user parameters:
```
gen_specs['user']['tasmanian_init'] -> a callable function with no inputs that returns an initialized object of Tasmanian.SparseGrid
gen_specs['user']['tasmanian_checkpoint_file'] -> a file for intermediary check pointing as well as the final grid
```
No additional parameters are required for static grid construction, e.g., similar to `TasGrid::loadNeededPoints()`. Iterative refinement can be enabled with the additional options for either anisotropic strategy
```
    gen_specs['user']['refinement'] = 'setAnisotropicRefinement'
    gen_specs['user']['sType'] = 'iptotal'
    gen_specs['user']['iMinGrowth'] = 10
    gen_specs['user']['iOutput'] = 0
```
or surplus strategy
```
    gen_specs['user']['refinement'] = 'setSurplusRefinement'
    gen_specs['user']['fTolerance'] = 1.E-2
    gen_specs['user']['sCriteria'] = 'classic'
    gen_specs['user']['iOutput'] = 0
```
See `TasGrid::TasmanianSparseGrid::setAnisotropicRefinement()` and `TasGrid::TasmanianSparseGrid::setSurplusRefinement()` for details.



