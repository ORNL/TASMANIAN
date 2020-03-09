/*!
 * \file TasmanianSparseGrid.i
 *
 * Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 * Distributed under a BSD 3-Clause License: see LICENSE for details.
 */

%{
#include "TasmanianSparseGrid.hpp"
%}

// Allow native fortran arrays to be passed to pointer/arrays
%include <typemaps.i>
%apply SWIGTYPE ARRAY[] {
    int*,
    float*,
    double*,
    int[],
    float[],
    double[]
};

namespace TasGrid {
// The array-like typemap can't be applied to return values
%ignore TasmanianSparseGrid::getHierarchicalCoefficients;
%ignore TasmanianSparseGrid::getLoadedValues;

%ignore TasmanianSparseGrid::setSurplusRefinement(double, int);
%ignore TasmanianSparseGrid::setSurplusRefinement(double, int, int const*);

%ignore TasmanianSparseGrid::getPointsIndexes;
%ignore TasmanianSparseGrid::getNeededIndexes;

%ignore readGrid(std::string const&);

%ignore TasmanianSparseGrid::operator=;
%ignore TasmanianSparseGrid::operator EvaluateCallable;
%ignore TasmanianSparseGrid::TasmanianSparseGrid(TasmanianSparseGrid &&);
%ignore TasmanianSparseGrid::copyGrid(const TasmanianSparseGrid *, int, int);
%ignore TasmanianSparseGrid::copyGrid(const TasmanianSparseGrid *, int);
%ignore TasmanianSparseGrid::copyGrid(const TasmanianSparseGrid *);
}

// Rename the factory methods to account for the lack of namespac
%rename(TasmanianGlobalGrid) TasGrid::makeGlobalGrid;
%rename(TasmanianSequenceGrid) TasGrid::makeSequenceGrid;
%rename(TasmanianLocalPolynomialGrid) TasGrid::makeLocalPolynomialGrid;
%rename(TasmanianFourierGrid) TasGrid::makeFourierGrid;
%rename(TasmanianWaveletGrid) TasGrid::makeWaveletGrid;
%rename(TasmanianReadGrid) TasGrid::readGrid;
%rename(TasmanianCopyGrid) TasGrid::copyGrid;

// Additional functional-style of helper methods
// Signature, makes them public for the module
%insert("fdecl") %{
    public :: tsgGetLoadedPoints
    public :: tsgGetNeededPoints
    public :: tsgGetPoints
    public :: tsgGetQuadratureWeights
%}
// implementation
%insert("fsubprograms") %{

function tsgGetLoadedPoints(grid) result(fresult)
    use, intrinsic :: ISO_C_BINDING
    type(TasmanianSparseGrid), intent(in) :: grid
    real(C_DOUBLE), dimension(:,:), pointer :: fresult

    allocate(fresult(grid%getNumDimensions(), grid%getNumLoaded()))
    call grid%getLoadedPoints(fresult(:,1))
end function

function tsgGetNeededPoints(grid) result(fresult)
    use, intrinsic :: ISO_C_BINDING
    type(TasmanianSparseGrid), intent(in) :: grid
    real(C_DOUBLE), dimension(:,:), pointer :: fresult

    allocate(fresult(grid%getNumDimensions(), grid%getNumNeeded()))
    call grid%getNeededPoints(fresult(:,1))
end function

function tsgGetPoints(grid) result(fresult)
    use, intrinsic :: ISO_C_BINDING
    type(TasmanianSparseGrid), intent(in) :: grid
    real(C_DOUBLE), dimension(:,:), pointer :: fresult

    allocate(fresult(grid%getNumDimensions(), grid%getNumPoints()))
    call grid%getPoints(fresult(:,1))
end function

function tsgGetQuadratureWeights(grid) result(fresult)
    use, intrinsic :: ISO_C_BINDING
    type(TasmanianSparseGrid), intent(in) :: grid
    real(C_DOUBLE), dimension(:), pointer :: fresult

    allocate(fresult(grid%getNumPoints()))
    call grid%getQuadratureWeights(fresult)
end function

%}





// Wrap clases/methods found in this file
%include "TasmanianSparseGrid.hpp"
