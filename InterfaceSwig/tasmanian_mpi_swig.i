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

%include "mpi.i"

%import tasmanian_swig.i

%ignore std::basic_streambuf<char, std::char_traits<char>>;
%ignore TasGrid::VectorToStreamBuffer;
%ignore getMPIRank(MPI_Comm comm);

#define Tasmanian_ENABLE_MPI
%{
#include <tsgMPIScatterGrid.hpp>
%}

%include "tsgMPIScatterGrid.hpp"

%rename(tsgMPIGridSend) TasGrid::MPIGridSend<true>;
%template(tsgMPIGridSendBin) TasGrid::MPIGridSend<true>;

%rename(tsgMPIGridRecv) TasGrid::MPIGridRecv<true>;
%template(tsgMPIGridRecvBin) TasGrid::MPIGridRecv<true>;

%rename(tsgMPIGridBcast) TasGrid::MPIGridBcast<true>;
%template(tsgMPIGridBcastBin) TasGrid::MPIGridBcast<true>;

%rename(tsgMPIGridScatterOutputs) TasGrid::MPIGridScatterOutputs<true>;
%template(tsgMPIGridScatterOutputsBin) TasGrid::MPIGridScatterOutputs<true>;
