/*!
 * \file tasmanian.i
 *
 * Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 * Distributed under a BSD 3-Clause License: see LICENSE for details.
 */

%module "tasmanian_mpi_swig"

/* -------------------------------------------------------------------------
 * Header definition macros
 * ------------------------------------------------------------------------- */

%define %tasmanian_add_header
%insert("fbegin") %{
! TASMANIAN project, https://github.com/ORNL/TASMANIAN
! Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
! Distributed under an MIT open source license: see LICENSE for details.
%}
%insert("begin") %{
/*
 * TASMANIAN project, https://github.com/ORNL/TASMANIAN
 * Copyright (c) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 * Distributed under an MIT open source license: see LICENSE for details.
 */
%}
%enddef

%tasmanian_add_header

/* -------------------------------------------------------------------------
 * Exception handling
 * ------------------------------------------------------------------------- */

// Rename the error variables' internal C symbols
#define SWIG_FORTRAN_ERROR_INT tasmanian_ierr
#define SWIG_FORTRAN_ERROR_STR tasmanian_get_serr

%include <extern_exception.i>

/* -------------------------------------------------------------------------
 * Data types and instantiation
 * ------------------------------------------------------------------------- */

// Note: stdint.i inserts #include <stdint.h>
%include <stdint.i>

// Support std::string
%include <std_string.i>

%import tasmanian_swig.i
