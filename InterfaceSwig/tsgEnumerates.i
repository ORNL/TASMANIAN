/*!
 * \file tsgEnumerates.i
 *
 * Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 * Distributed under a BSD 3-Clause License: see LICENSE for details.
 */

%{
#include "tsgEnumerates.hpp"
%}

// Wrap types found in this file as native fortran enums
%fortranconst;
%ignore accel_gpu_hip;
%ignore accel_gpu_rocblas;
%include "tsgEnumerates.hpp"
%nofortranconst;
