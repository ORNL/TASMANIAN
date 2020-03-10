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

// Using custom factory methods to take advantage of the Fortran optional parameters
// Factory methods have many-many parameters, better optional ones helps a lot
// makeEmpty() is equivalent to TasmanianSparseGrid(), can ignore
%ignore TasGrid::makeGlobalGrid;
%ignore TasGrid::makeSequenceGrid;
%ignore TasGrid::makeFourierGrid;
%ignore TasGrid::makeLocalPolynomialGrid;
%ignore TasGrid::makeWaveletGrid;
%ignore TasGrid::makeEmpty;
%rename(TasmanianReadGrid) TasGrid::readGrid;
%rename(TasmanianCopyGrid) TasGrid::copyGrid;

// Additional functional-style of helper methods
// Signature, makes them public for the module
%insert("fdecl") %{
    public :: TasmanianGlobalGrid
    public :: TasmanianSequenceGrid
    public :: TasmanianFourierGrid
    public :: TasmanianLocalPolynomialGrid
    public :: TasmanianWaveletGrid
    public :: tsgGetLoadedPoints
    public :: tsgGetNeededPoints
    public :: tsgGetPoints
    public :: tsgGetQuadratureWeights
%}
// implementation
%insert("fsubprograms") %{

function TasmanianGlobalGrid(dims, outs, depth, gtype, rule, &
                             aweights, alpha, beta, custom_filename, level_limits) result(grid)
    type(TasmanianSparseGrid) :: grid
    integer(C_INT), intent(in) :: dims, outs, depth, gtype, rule
    integer(C_INT), optional, target  :: aweights(:)
    double precision, optional :: alpha, beta
    character(len=*), target, optional :: custom_filename
    integer(C_INT), optional :: level_limits(dims)
    real(C_DOUBLE) :: al, be
    integer(C_INT), dimension(:), pointer  :: aw(:)

    if ( present(aweights) ) then
        aw => aweights
    else
        allocate(aw(2 * dims))
        aw(1:dims) = 1;
        aw(dims+1:2*dims) = 0;
    endif
    if ( present(alpha) ) then
        al = alpha
    else
        al = 0.0D-0
    endif
    if ( present(beta) ) then
        be = beta
    else
        be = 0.0D-0
    endif

    grid%swigdata = swigc_new_TasmanianSparseGrid__SWIG_0()
    if ( present(custom_filename) ) then
        if ( present(level_limits) ) then
            call grid%makeGlobalGrid(dims, outs, depth, gtype, rule, aw, al, be, custom_filename, level_limits)
        else
            call grid%makeGlobalGrid(dims, outs, depth, gtype, rule, aw, al, be, custom_filename)
        endif
    else
        if ( present(level_limits) ) then
            call grid%makeGlobalGrid(dims, outs, depth, gtype, rule, aw, al, be, "", level_limits)
        else
            call grid%makeGlobalGrid(dims, outs, depth, gtype, rule, aw, al, be)
        endif
    endif
    if ( .not. present(aweights) ) then
        deallocate(aw)
    endif
end function

function TasmanianSequenceGrid(dims, outs, depth, gtype, rule, aweights, level_limits) result(grid)
    type(TasmanianSparseGrid) :: grid
    integer(C_INT) :: dims, outs, depth, gtype, rule
    integer(C_INT), optional, target  :: aweights(:)
    integer(C_INT), optional :: level_limits(dims)
    integer(C_INT) :: aw(2 * dims)

    grid%swigdata = swigc_new_TasmanianSparseGrid__SWIG_0()
    if ( present(aweights) ) then
        if ( present(level_limits) ) then
            call grid%makeSequenceGrid(dims, outs, depth, gtype, rule, aweights, level_limits)
        else
            call grid%makeSequenceGrid(dims, outs, depth, gtype, rule, aweights)
        endif
    else
        if ( present(level_limits) ) then
            aw(1:dims) = 1;
            aw(dims+1:2*dims) = 0;
            call grid%makeSequenceGrid(dims, outs, depth, gtype, rule, aw, level_limits)
        else
            call grid%makeSequenceGrid(dims, outs, depth, gtype, rule)
        endif
    endif
end function

function TasmanianLocalPolynomialGrid(dims, outs, depth, order, rule, level_limits) result(grid)
    type(TasmanianSparseGrid) :: grid
    integer(C_INT) :: dims, outs, depth, ord, ru
    integer(C_INT), optional :: order, rule
    integer(C_INT), optional :: level_limits(dims)

    if ( present(order) ) then
        ord = order
    else
        ord = 1
    endif
    if ( present(rule) ) then
        ru = rule
    else
        ru = tsg_rule_localp
    endif

    grid%swigdata = swigc_new_TasmanianSparseGrid__SWIG_0()
    if ( present(level_limits) ) then
        call grid%makeLocalPolynomialGrid(dims, outs, depth, ord, ru, level_limits)
    else
        call grid%makeLocalPolynomialGrid(dims, outs, depth, ord, ru)
    endif
end function

function TasmanianWaveletGrid(dims, outs, depth, order, level_limits) result(grid)
    type(TasmanianSparseGrid) :: grid
    integer(C_INT) :: dims, outs, depth, ord
    integer(C_INT), optional :: order
    integer(C_INT), optional :: level_limits(dims)

    if ( present(order) ) then
        ord = order
    else
        ord = 1
    endif

    grid%swigdata = swigc_new_TasmanianSparseGrid__SWIG_0()
    if ( present(level_limits) ) then
        call grid%makeWaveletGrid(dims, outs, depth, ord, level_limits)
    else
        call grid%makeWaveletGrid(dims, outs, depth, ord)
    endif
end function

function TasmanianFourierGrid(dims, outs, depth, gtype, aweights, level_limits) result(grid)
    type(TasmanianSparseGrid) :: grid
    integer(C_INT) :: dims, outs, depth, gtype
    integer(C_INT), dimension(:), optional :: aweights(:)
    integer(C_INT), optional :: level_limits(dims)
    integer(C_INT) :: aw(2 * dims)

    grid%swigdata = swigc_new_TasmanianSparseGrid__SWIG_0()
    if ( present(aweights) ) then
        if ( present(level_limits) ) then
            call grid%makeFourierGrid(dims, outs, depth, gtype, aweights, level_limits)
        else
            call grid%makeFourierGrid(dims, outs, depth, gtype, aweights)
        endif
    else
        if ( present(level_limits) ) then
            aw(1:dims) = 1;
            aw(dims+1:2*dims) = 0;
            call grid%makeFourierGrid(dims, outs, depth, gtype, aw, level_limits)
        else
            call grid%makeFourierGrid(dims, outs, depth, gtype)
        endif
    endif
end function

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
