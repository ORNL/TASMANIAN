# Tasmanian

The Toolkit for Adaptive Stochastic Modeling and Non-Intrusive ApproximatioN is a collection of robust libraries for high dimensional integration and interpolation as well as parameter calibration. The code consists of several modules that can be used individually or conjointly.

Visit us at: [http://tasmanian.ornl.gov/](http://tasmanian.ornl.gov/)

Sparse Grids
--------------

Sparse Grids is a family of algorithms for constructing multidimensional quadrature and interpolation rules
using multiple tensor products of one dimensional rules with varying degree of precision.
The Tasmanian Sparse Grids Module implements a variety of grids that fall into five major categories:
* **Global grids** use polynomials supported over the entire domain of integration or interpolation.
    Such grids are suitable for approximating functions with smooth response.
* **Sequence grids** work much like Global grids, but use optimized internal data-structures for rules that
    are based on sequences of points formed from solving a greedy optimization problem
* **Local polynomial grids** use hierarchical piece-wise polynomials with compact support.
    Such grids are suitable for approximating functions with localized sharp behavior.
* **Wavelet grids** use special functions that form a Riesz basis which allows for more accurate
    local error estimations. Compared to Local polynomial grids, the wavelet basis can provide
    similar accuracy with significantly fewer points, but the advantage of the Riesz basis could also
    be negated from the higher Lebesgue constant near the domain boundary.
* **Fourier grids** use trigonometric functions with varying frequency and rely on Discrete (complex)
    Fourier Transforms. Such grids are suitable for approximating periodic functions,
    since periodicity if preserved in the approximation.

DREAM
--------------

The DiffeRential Evolution Adaptive Metropolis is a method to draw samples from an arbitrary probability
distribution defined by an arbitrary non-negative function (not necessarily normalized to integrate to 1).
The DREAM approach is similar to the classical Markov Chain Monte Carlo, but it evolves a large number
of chains simultaneously which leads to better parallelization and (potentially) faster convergence.
In addition, multiple chains allow for better exploration of the probability domain, which is often
advantageous when working with multi-modal distributions.

One of the main applications of DREAM is in the field of Bayesian inference, where samples are drawn
from a posterior distribution comprised from a data-informed likelihood and an arbitrary model.
The DREAM module of Tasmanian can use Tasmanian Sparse Grids approximation to either the model
or the likelihood.


### Please cite us
If you use Tasmanian for your research, please cite the Manual and our work on global and locally adaptive grids.

[http://tasmanian.ornl.gov/documents/Tasmanian.bib](http://tasmanian.ornl.gov/documents/Tasmanian.bib)

Quick Install
--------------

* the basic way: using GNU Make, `g++` and optionally `gfortran` and `/usr/bin/python`
```
  make
  make test     (will fail if /usr/bin/env python is missing numpy or ctypes modules)
  make matlab   (optional: sets matlab work folder to ./tsgMatlabWorkFolder/)
  make python3  (optional: sets #!/usr/bin/env python3)
  make fortran  (optional: compile Fortran libraries)
  make examples
  make clean
```
* the easy way: using cmake and the `install` script
```
  ./install <install-path> <optional: matlab-work-folder> <additional options>
  ./install --help  (lists all options)
```
* the CMake way: see the detailed instruction for a full list of commands
```
  mkdir Build
  cd Build
  cmake <options> <path-to-Tasmanian-source>
  make
  make test
  make install
  make test_install
```
* under MS Windows: use the cmake gui to set the folders and options then use the command prompt (`cmd.exe`) to enter the build folder
```
  cd <cmake build folder>
  cmake --build . --config Release
  ctest -C Release
  cmake --build . --config Release --target install
  # both Debug and Release are the supported config modes above
```
* Tasmanian is also included in Spack: [https://spack.io/](https://spack.io/)

Basic Usage of Tasmanian
--------------

* see the Examples in the install prefix:
```
  <install-prefix>/share/Tasmanian/examples/
```
* the Examples source code is in:
```
  <source-root>/<module>/Examples/
```
