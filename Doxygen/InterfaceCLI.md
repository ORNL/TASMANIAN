# Command Line Interface

The **tasgrid** executable is a command line interface to the Tasmanian Sparse Grid module. It provides the ability to create and manipulate sparse grids, save and load them into files and optionally interface with another program via text files. For the most part, **tasgrid** reads a grid from a file, calls one or more of the functions described in the previous section and then saves the resulting grid.

The commands for **tasgrid** correspond to calls to the C++ API, where scalar inputs are given as command line arguments and vector/array parameters are given as matrix files, see the end of this section for the matrix file format.

```
  ./tasgrid <command> <option1> <value1> <option2> <value2> ....
```

The first input to the executable is the command that specifies the action that needs to be taken. The command is followed by options and values.

Every command is associated with a set of options, extra options are ignored. See the help subsection on how to find the right options for each command.

### Example commands

```
  ./tasgrid -mq -dim 4 -depth 2 -type qptotal -1d gauss-legendre -p
```

Make quadrature rule in 4 dimensions that can integrate exactly all quadratic polynomials, print the result to the screen. Note that the first column is the weight.

```
  ./tasgrid -mg -dim 3 -out 2 -depth 4 -type iptotal -1d clenshaw-curtis -gf example_grid_file
  ./tasgrid -l -gf example_grid_file -vf file_with_values
  ./tasgrid -e -gf example_grid_file -xf file_with_points -of result_file
```
The first command creates a global grid with 3 dimensions and clenshaw-curtis points that interpolate exactly all polynomials of order 4. The grid is stored in *example_grid_file*. In the second command, model values are read from the *file_with_values* and loaded into the grid. In the final command, the interpolant is evaluated at the points specified in *file_with_points* and the result is stored in the last file.


### Command: -h, help, -help, --help

```
  ./tasgrid --help
  ./tasgrid -makequadrature help
```

Prints information about the usage of **tasgrid**.  In addition, writing **help** after any command will print information specific to that command; effectively, **help** is a universal option.


### Commands and corresponding C++ functions

```
  ./tasgrid -makeglobal         ->  makeGlobalGrid()
  ./tasgrid -makesequence       ->  makeSequenceGrid()
  ./tasgrid -makelocalpoly      ->  makeLocalPolynomialGrid()
  ./tasgrid -makewavelet        ->  makeWaveletGrid()
  ./tasgrid -makefourier        ->  makeFourierGrid()
  ./tasgrid -makequadrature     ->  (one of the grids above, see comments)
  ./tasgrid -makeupdate         ->  updateGlobalGrid()/updateSequenceGrid()
  ./tasgrid -setconformal       ->  setConformalTransformASIN()
  ./tasgrid -getquadrature      ->  getQuadratureWeights()/getPoints()
  ./tasgrid -getinterweights    ->  getInterpolationWeights()
  ./tasgrid -getpoints          ->  getPoints()
  ./tasgrid -getneededpoints    ->  getNeededPoints()
  ./tasgrid -loadvalues         ->  loadNeededPoints()
  ./tasgrid -evaluate           ->  evaluateBatch()
  ./tasgrid -evalhierarchyd     ->  evaluateHierarchicalFunctions()
  ./tasgrid -evalhierarchys     ->  evaluateSparseHierarchicalFunctions()
  ./tasgrid -integrate          ->  integrate()
  ./tasgrid -getanisotropy      ->  estimateAnisotropicCoefficients()
  ./tasgrid -refineaniso        ->  setAnisotropicRefinement()
  ./tasgrid -refinesurp         ->  setSurplusRefinement()
  ./tasgrid -refine             ->  setAnisotropicRefinement()/setSurplusRefinement()
  ./tasgrid -cancelrefine       ->  clearRefinement()
  ./tasgrid -mergerefine        ->  mergeRefinement()
  ./tasgrid -getcoefficients    ->  getHierarchicalCoefficients()
  ./tasgrid -setcoefficients    ->  setHierarchicalCoefficients()
  ./tasgrid -getpoly            ->  getGlobalPolynomialSpace()
  ./tasgrid -summary            ->  printStats()
  ./tasgrid <command> help      -> show more info for this command
```

Additional notes:
* The domain types for all grids are set during the *make* command, domains cannot be changed with the **tasgrid** executable since domain changes always change the nodes and effectively generates a new grid.
* Make quadrature creates a grid with zero outputs and type that is based on the one dimensional rule, e.g., the **-makeupdate** grid will automatically detect sequence or global grids.
* The **-getquadrature** command will generate larger matrix, where the first column is the weights and the rest correspond to the points.
* The **-getinterweights** command can work with multiple points at a time, the call will use OpenMP (if available).
* The **-refine** command will call anisotropic refinement on Global, Sequence, and Fourier grids, and surplus refinement for Local Polynomial and Wavelet grids.
* The coefficients and hierarchical functions for Fourier grids work with complex numbers, meaning that each pair of consecutive numbers correspond to one complex number (real and complex parts). This the matrices have twice as many columns. Note that this also applies to the coefficients as inputs and outputs (which is different from the C++ API).
* The **-evaluate** command accepts **-gpuid** options, which allows to select a CUDA device to use for acceleration. If the option is omitted, GPU acceleration will not be used.


### Command: -listtypes

```
  ./tasgrid -listtypes
```

List the available one dimensional quadrature and interpolation rules as well as the different types of grids, refinement and conformal mapping types. Use this command to see the correct spelling of all string options.

### Command: -version or -info

```
  ./tasgrid -version
  ./tasgrid -v
  ./tasgrid -info
```
Prints the version of the library, the available acceleration options, and (if CUDA is enabled) the visible CUDA devices.


### Command: -test

```
./tasgrid -test
./tasgrid -test random
./tasgrid -test verbose
```

Since Tasmanian 6.0 the sparse grids testing is moved to a different executable, i.e., **gridtest**. The test method of **tasgrid** is still included but it covers only a sub-set of the tests. Also, the tests take longer, especially when CUDA is enables, since large reference solutions have to be computed on the slower CPU. The **gridtest** executable takes the **random** and **verbose** switches, but does not need the **-test** command.

The tests rely on random number generation to estimate the accuracy of computed interpolants. If the test fails, this may be indication of a problem with the hard-coded random seed. Using the **random** option will reset the seed on every run and will provide more statistically significant results.

The **verbose** will print more detailed output. This affects only the successful tests, failed tests always print verbose information.


### Matrix File Format

The matrix files have two formats, binary and ASCII. The simple text file describes a two dimensional array of real (double-precision) numbers. The file contains two integers on the first line indicating the number of rows and columns. Those are followed by the actual entries of the matrix one row at a time.

The file containing

```
3 4
1.0  2.0  3.0  4.0
5.0  6.0  7.0  8.0
9.0 10.0 11.0 12.0
```
represents the matrix
\f$
	\left( \begin{array}{rrrr}
	1 & 2 & 3 & 4 \\
	5 & 6 & 7 & 8 \\
	9 & 10 & 11 & 12 \\
	\end{array} \right)
\f$


A matrix file may contain only one row or column, e.g.,
```
1 2
13.0 14.0
```

In binary format, the file starts with three characters *TSG* indicating that this is a binary Tasmanian file. The characters are followed by two integers and the double-precision numbers that correspond to the matrix being read left-to-right top-to-bottom.

All files used by **tasgrid** have the above format with three exceptions. The **-gridfile** option contains saved sparse grids and it is not intended for editing outside of the **tasgrid** calls (it is OK to modify using other Tasmanian API calls, e.g., C++ or Python). The **-anisotropyfile** option requires a matrix with one column and it should contain double-precision numbers that have integer values. The **-customrulefile** has special format is described the Custom Rule File Format section.

The default mode is to use binary files for all calls to **tasgrid**, but ASCII files are easier to debug and potentially easier to import to external codes. This **tasgrid** has an option **-ascii** that can be added to any command and will force the resulting output to be written in ASCII format.
