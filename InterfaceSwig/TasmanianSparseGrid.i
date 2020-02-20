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

// The array-like typemap can't be applied to return values
%ignore TasGrid::TasmanianSparseGrid::getHierarchicalCoefficients;
%ignore TasGrid::TasmanianSparseGrid::getLoadedValues;

// Wrap clases/methods found in this file
%include "TasmanianSparseGrid.hpp"
