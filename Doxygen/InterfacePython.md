# Python Interface

The Python interface uses `c_types` and links to the C interface of Tasmanian.
The C interface uses a series of functions that take a `void*` pointer that references
an instance of the C++ class TasGrid::TasmanianSparseGrid.
The Python module takes a hold of the C void pointer and encapsulates it into a Python class
with interface near identical to that of C++.
The vectors are replaced by 1D and 2D `numpy.ndarray` structures and the enumerated types
are replaced by strings.


### Requirements

Both Python 2 and 3 are supported and the installed module is fully compatible with both
regardless of the Python version used during install and testing.
Tasmanian uses the following modules:
* **required** `c_types`
* **required** `numpy`
* **optional** `matplotlib.pyplot`


### Module

Starting with Tasmanian 7.0 all Python capabilities are included in a single module:
```
    import Tasmanian
```
Older version of Tasmanian have Python bindings only for the Sparse Grid class
and the module was called `TasmanianSG`. The API is backwards compatible
there is no need to change the old code.

Note that the Tasmanian install path must be included in the Python search path:
```
    <install-prefix>/lib/pythonX.Y/site-packages/ # following the Python convention
    <install-prefix>/share/Tasmanian/python/      # Tasmanian version independent path
```
Either path is sufficient if it is added to environment variable `PYTHONPATH`
(e.g., via the `TasmanianENVsetup.sh` script) or with `sys.path.append()` command.

Example for creating an instance of a sparse grid:
```
    import Tasmanian
    grid = Tasmanian.SparseGrid()
    grid.makeGlobalGrid(...)
```
See the Python examples in `<install-prefix>/share/Tasmanian/examples/`.


### Additional Documentation

Every class and method in Tasmanian comes with additional comments designed to be
accessed with the Python `help()` function, e.g.,
```
    help(Tasmanian.SparseGrid.makeGlobalGrid)
```
Note that the names of the class methods match the C++ interface.


### Error Handling

The Tasmanian Python module comes with a comprehensive error checking mechanisms
where the validity of the inputs to every method are checked and `TasmanianInputError`
exception is raised whenever invalid input is encountered.
The `TasmanianInputError` class has two variables:
* **sVariable** indicates the variable with invalid input or the name of the method
  if the method was called out of order, e.g., `evaluate()` cannot be called before `loadNeededPoints()`
* **sMessage** is a human readable string giving specific description to of the error.


### Acceleration and Thread Safety

All acceleration modes are available through the Python interface and the same
rules for thread safety apply as in the C++ interface.
One point of difference is the `evaluate()` function, which on the Python side
will link to `SparseGrid.evaluateFast()` and therefore not safe to call from different threads.
Safe evaluations can be performed with `SparseGrid.evaluateThreadSafe()`


### Basic Plotting

If `matplotlib.pyplot` module is available, `Tasmanian.SparseGrid` will enable two plotting methods:
* **plotPoints2D()** plots the nodes associated with the grid
* **plotResponse2D()** plots a color image corresponding to the response surface
Both methods work only with two dimensional grids.

Note that `matplotlib.pyplot` doesn't have to be available during build time,
it can be installed later and still used with a current installation of Tasmanian.
