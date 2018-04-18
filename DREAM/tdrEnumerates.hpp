#ifndef __TASMANIAN_DREAM_ENUMERATES_HPP
#define __TASMANIAN_DREAM_ENUMERATES_HPP

#include <iostream>
#include <cstdlib>
#include <math.h>

#include "TasmanianConfig.hpp"

#ifdef Tasmanian_ENABLE_MPI
#include <mpi.h>
#endif

namespace TasDREAM{

using std::cout; // for debugging purposes
using std::endl; // for debugging purposes

using std::cerr; // for error messages

enum TypeDistribution{
    dist_uniform,
    dist_gaussian,
    dist_truncated_gaussian,
    dist_weibull,
    dist_exponential,
    dist_beta,
    dist_gamma,
    dist_custom
};

enum TypeLikelihood{
    likely_gauss_scale, // scale, diagonal and dense refer to the type of covariance
    likely_gauss_diagonal,
    likely_gauss_dense
};

}

#endif
