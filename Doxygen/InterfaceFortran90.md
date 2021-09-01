# Interface Fortran 90/95

*note that the Fortran 90 interface is deprecated, the object oriented Fortran 2003 interface is the currently recommended one*

Tasmanian 7.0 comes with a mostly-stable Fortran 90/95 module that wraps around the C++ sparse grids library and gives access to all of the core functionality. The interface uses a user defined type that is effectively a pointer to a **TasmanainSparseGrid** objects. A new grid has to be initialized with **tsgAllocateGrid()**  which creates an empty grid, then the grid can be used in all follow on calls to the Fortran module. Once work is done with the grid, **tsgDeallocateGrid()** must be called to free all memory used by this grid. For example:
```
program foo
  use TasmanianSG
  type(TasmanianSparseGrid) :: grid
  double precision, pointer :: points(:,:), weights(:)

  call tsgAllocateGrid(grid)
  call tsgMakeGlobalGrid(grid, 2, 0, 1, tsg_level, tsg_leja)
  points => tsgGetPoints(grid)
  weights => tsgGetQuadratureWeights(grid)
  ... use the quadrature ...
  call tsgDeallocateGrid(grid)
end program foo
```
The enumerate types used by the C++ interface are replaced by hard-coded constant integers, i.e., the module defines a parameter integer **tsg_clenshaw_curtis**, which can be used to call **tsgMakeGlobalGrid()**. All vectors are replaced by Fortran vectors and 2D matrixes, similar to the calls to the command line tool **tasgrid**. However, Fortran is using column-major storage format, hence the row-column format is reversed to what you see in Python or **tasgrid** calls, i.e., when the C++ interface refers to the first point which is stored in the first set of* **getNumDimensions()** entries, that means the first column of the Fortran matrix.

All defined integer parameters start wit **tsg_** and all functions/subroutines start with **tsg**. A full list is provided at the top of TasmanianSG.f90 file. Examples are provided in the examples_sparse_grids.f90 file.

The Fourier grids use complex numbers and hence special functions are provides:
```
  tsgEvaluateComplexHierarchicalFunctions(...)
  tsgGetComplexHierarchicalCoefficients(...)
  tsgGetComplexHierarchicalCoefficientsStatic(...)
```

The interface comes with negligible overhead, except for the case of sparse hierarchical functions, where extra overhead is needed since there is no *clean* way to pass pointer to dynamic memory created in C++ to Fortran 90/95.
