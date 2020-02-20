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
    const int*,
    const double[],
    double[]
};

// Wrap clases/methods found in this file
%include "TasmanianSparseGrid.hpp"
