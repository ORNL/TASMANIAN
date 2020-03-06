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

// Wrap clases/methods found in this file
%include "TasmanianSparseGrid.hpp"
