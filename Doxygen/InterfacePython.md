# Python Interface

The Python interface uses `c_types` and links to the C interface of Tasmanian
as well as several additional C-style functions, e.g., when dealing with Addon templates.
The C++ classes are converted to `void*` pointers and then wrapped inside Python
classes with near identical interface.
The C++ vectors are replaced by 1D and 2D `numpy.ndarray` structures,
the enumerated types are replaced by strings,
and the C++ `std::function` and `lambda` expressions by Python callable objects and lambdas.
See the Fortran 2003 comments about the row major matrix format.


### Requirements

Both Python 2 and 3 are supported and the installed module is fully compatible with both
regardless of the Python version used during install and testing.
Tasmanian uses the following modules:
* **required** `sys`
* **required** `c_types`
* **required** `numpy`
* **optional** `matplotlib.pyplot`

Tasmanian is also available through PIP: [https://pypi.org/project/Tasmanian/](https://pypi.org/project/Tasmanian/);
however, the PIP installation comes with additional requirements for facilitating the CMake build process.
See the PIP section on the installation instructions page.

### Module

Starting with Tasmanian 7.0 all Python capabilities are included in a single module:
```
    import Tasmanian
```
Older version of Tasmanian have Python bindings only for the Sparse Grid class
and the module was called `TasmanianSG`. The API is backwards compatible
and there is no need to change existing codes.

Note that the Tasmanian install path must be included in the Python search path:
```
    <install-prefix>/lib/pythonX.Y/site-packages/ # following the Python convention
    <install-prefix>/share/Tasmanian/python/      # Tasmanian version independent path
```
Either path is sufficient if it is added to the environment variable `PYTHONPATH`
(e.g., via the `TasmanianENVsetup.sh` script) or with the `sys.path.append()` command.

Example for creating an instance of a sparse grid:
```
    import Tasmanian
    grid = Tasmanian.SparseGrid()
    grid.makeGlobalGrid(...)
```
Example for using DREAM sampling:
```
    import Tasmanian
    state = Tasmanian.DREAM.State(...)
    state.setState(Tasmanian.DREAM.genUniformSamples(...))
    Tasmanian.DREAM.Sample(iNumBurnupIterations, iNumCollectIterations,
                           lambda x : numpy.exp( - x[:,0]**2 ),
                           Tasmanian.DREAM.Domain("unbounded"),
                           state,
                           Tasmanian.DREAM.IndependentUpdate("uniform", 0.5),
                           Tasmanian.DREAM.DifferentialUpdate(90))
```
See the Python examples in `<install-prefix>/share/Tasmanian/examples/`.


### Additional Documentation

Every class and method in Tasmanian comes with additional comments designed to be
accessed with the Python `help()` function, e.g.,
```
    help(Tasmanian.SparseGrid.makeGlobalGrid)
```

* Sparse grid module and C++ equivalent:
```
    Tasmanian.SparseGrid (class)                  -> TasGrid::TasmanianSparseGrid (class)
    Tasmanian.TasmanianSimpleSparseMatrix (class) -> triplet vectors (pntr, indx, vals)
```
* DREAM module and C++ equivalent:
```
         import Tasmanian.DREAM as DREAM
(class)  DREAM.State                        -> TasDREAM::TasmanianDREAM              (class)
(class)  DREAM.LikelihoodGaussIsotropic     -> TasDREAM::LikelihoodGaussIsotropic    (class)
(class)  DREAM.LikelihoodGaussAnisotropic   -> TasDREAM::LikelihoodGaussAnisotropic  (class)
(class)  DREAM.Posterior                    -> TasDREAM::posterior                  (method)
(method) DREAM.PosteriorEmbeddedLikelihood  -> TasDREAM::posterior         (method-overload)
(method) DREAM.genUniformSamples            -> TasDREAM::genUniformSamples          (method)
(method) DREAM.genGaussianSamples           -> TasDREAM::genGaussianSamples         (method)
(method) DREAM.Sample                       -> TasDREAM::SampleDREAM                (method)
```
* Addon module methods and C++ equivalent:
```
    Tasmanian.loadNeededValues              -> TasGrid::loadNeededValues <overwrite_loaded=false>
    Tasmanian.reloadLoadedValues            -> TasGrid::loadNeededValues <overwrite_loaded=true>
    Tasmanian.constructAnisotropicSurrogate -> TasGrid::constructSurrogate (anistropic overload)
    Tasmanian.constructSurplusSurrogate     -> TasGrid::constructSurrogate (surplus overload)
```


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
One point of difference is the `SparseGrid.evaluate()` function, which on the Python side
will link to `SparseGrid.evaluateFast()` and therefore not safe to call from different threads.
Safe evaluations can be performed with `SparseGrid.evaluateThreadSafe()`.

The model callable inputs to load/reload and construct methods will be executed in
parallel on most Python implementations.


### Basic Plotting

If `matplotlib.pyplot` module is available, `Tasmanian.SparseGrid` will enable two plotting methods:
* **plotPoints2D()** plots the nodes associated with the grid
* **plotResponse2D()** plots a color image corresponding to the response surface
Both methods work only with two dimensional grids.

Note that `matplotlib.pyplot` doesn't have to be available during build time,
it can be installed later and still used with a current installation of Tasmanian.
