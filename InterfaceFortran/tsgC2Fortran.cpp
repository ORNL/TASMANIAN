/*
 * Copyright (c) 2017, Miroslav Stoyanov
 *
 * This file is part of
 * Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
 *    and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND IMPLIED.
 * THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT,
 * COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.
 * THE USER ASSUMES RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR ARISING OUT OF,
 * IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
 */

#ifndef __TSG_C2FORT_CPP
#define __TSG_C2FORT_CPP

#include "TasmanianAddons.hpp"

using std::cerr;
using std::endl;

using namespace TasGrid;

extern "C" void tsgc2fstr_(int *length, const char *str);
extern "C" void tsgc2fvec_(int *length, double *vect);
extern "C" void tsgc2fmat_(int *rows, int *cols, double *mat);

extern "C"{

////////////////////////////////////////////////////////////////////////
//   MAIN INTERFACE
////////////////////////////////////////////////////////////////////////

void tsggvm_(int *ver){ *ver = TasmanianSparseGrid::getVersionMajor(); }
void tsggvn_(int *ver){ *ver = TasmanianSparseGrid::getVersionMinor(); }
void tsggli_(){
    // this is clumsy, may have to think of something else
    const char *lic = TasmanianSparseGrid::getLicense();
    int l = 0;
    while(lic[l] != '\0') l++;
    tsgc2fstr_(&l, lic);
}

// read/write
void tsgwri_(TasmanianSparseGrid **grid, const char *filename, int *binary){ (*grid)->write(filename, (*binary != 0)); }
void tsgrea_(TasmanianSparseGrid **grid, const char *filename, int *status){
    try{
        (*grid)->read(filename);
        *status = 1;
    }catch(std::runtime_error &e){
        cerr << e.what() << endl;
        *status = 0;
    }
}

void tsgalloc_(TasmanianSparseGrid **grid){
    (*grid) = new TasmanianSparseGrid();
}
void tsgfree_(TasmanianSparseGrid **grid){
    delete (*grid);
}

// create
void tsgmg_(TasmanianSparseGrid **grid, int *dimensions, int *outputs, int *depth, int *type, int *rule, int *opt_flags,
            const int *aniso_weights, double *alpha, double *beta, const char *customRuleFilename, const int *llimits){
    const int   *aw = (opt_flags[0] != 0 ) ? aniso_weights      : nullptr;
    double       al = (opt_flags[1] != 0 ) ? *alpha             : 0.0;
    double       be = (opt_flags[2] != 0 ) ? *beta              : 0.0;
    const char *cfn = (opt_flags[3] != 0 ) ? customRuleFilename : nullptr;
    const int   *ll = (opt_flags[4] != 0 ) ? llimits            : nullptr;

    (*grid)->makeGlobalGrid(*dimensions, *outputs, *depth, IO::getDepthTypeInt(*type), IO::getRuleInt(*rule), aw, al, be, cfn, ll);
}
void tsgms_(TasmanianSparseGrid **grid, int *dimensions, int *outputs, int *depth, int *type, int *rule, int *opt_flags,
            const int *aniso_weights, const int *llimits){
    const int *aw  = (opt_flags[0] != 0 ) ? aniso_weights      : nullptr;
    const int *ll  = (opt_flags[1] != 0 ) ? llimits            : nullptr;

    (*grid)->makeSequenceGrid(*dimensions, *outputs, *depth, IO::getDepthTypeInt(*type), IO::getRuleInt(*rule), aw, ll);
}
void tsgml_(TasmanianSparseGrid **grid, int *dimensions, int *outputs, int *depth, int *opt_flags,
            int *order, int *rule, const int *llimits){
    int ru         = (opt_flags[0] != 0) ? *rule   : 1;
    int ord        = (opt_flags[1] != 0) ? *order  : 1;
    const int *ll  = (opt_flags[2] != 0) ? llimits : nullptr;

    (*grid)->makeLocalPolynomialGrid(*dimensions, *outputs, *depth, ord, IO::getRuleInt(ru), ll);
}
void tsgmw_(TasmanianSparseGrid **grid, int *dimensions, int *outputs, int *depth, int *opt_flags, int *order, const int *llimits){
    int ord       = (opt_flags[0] != 0) ? *order  : 1;
    const int *ll = (opt_flags[1] != 0) ? llimits : nullptr;

    (*grid)->makeWaveletGrid(*dimensions, *outputs, *depth, ord, ll);
}
void tsgmf_(TasmanianSparseGrid **grid, int *dimensions, int *outputs, int *depth, int *type, int *opt_flags, const int *aniso_weights, const int *llimits){
    const int *aw = (opt_flags[0] != 0 ) ? aniso_weights : nullptr;
    const int *ll = (opt_flags[1] != 0 ) ? llimits       : nullptr;

    (*grid)->makeFourierGrid(*dimensions, *outputs, *depth, IO::getDepthTypeInt(*type), aw, ll);
}

// copy/updare
void tsgcp_(TasmanianSparseGrid **grid, TasmanianSparseGrid **source){ (*grid)->copyGrid(*source); }

void tsgug_(TasmanianSparseGrid **grid, int *depth, int *type, int *opt_flags, const int *anisotropic_weights){
    const int *aw = (opt_flags[0] != 0) ? anisotropic_weights : nullptr;
    (*grid)->updateGlobalGrid(*depth, IO::getDepthTypeInt(*type), aw);
}
void tsgus_(TasmanianSparseGrid **grid, int *depth, int *type, int *opt_flags, const int *anisotropic_weights){
    const int *aw = (opt_flags[0] != 0) ? anisotropic_weights : nullptr;
    (*grid)->updateSequenceGrid(*depth, IO::getDepthTypeInt(*type), aw);
}

// getAlpha/Beta/Order/Dims/Outs/Rule //
void tsggal_(TasmanianSparseGrid **grid, double *alpha){ *alpha = (*grid)->getAlpha(); }
void tsggbe_(TasmanianSparseGrid **grid, double *beta){ *beta = (*grid)->getBeta(); }
void tsggor_(TasmanianSparseGrid **grid, int *order){ *order = (*grid)->getOrder(); }
void tsggnd_(TasmanianSparseGrid **grid, int *dims){ *dims = (*grid)->getNumDimensions(); }
void tsggno_(TasmanianSparseGrid **grid, int *outs){ *outs = (*grid)->getNumOutputs(); }
void tsggru_(TasmanianSparseGrid **grid, int *rule){ *rule = IO::getRuleInt((*grid)->getRule()); }

// getNumNeeded/Loaded/Points
void tsggnn_(TasmanianSparseGrid **grid, int *num_points){ *num_points = (*grid)->getNumNeeded(); }
void tsggnl_(TasmanianSparseGrid **grid, int *num_points){ *num_points = (*grid)->getNumLoaded(); }
void tsggnp_(TasmanianSparseGrid **grid, int *num_points){ *num_points = (*grid)->getNumPoints(); }

// getLoaded/Needed/Points
void tsgglp_(TasmanianSparseGrid **grid, double *points){ (*grid)->getLoadedPoints(points); }
void tsggdp_(TasmanianSparseGrid **grid, double *points){ (*grid)->getNeededPoints(points); }
void tsggpp_(TasmanianSparseGrid **grid, double *points){ (*grid)->getPoints(points); }

// getQuadrature/InterpolationWeights
void tsggqw_(TasmanianSparseGrid **grid, double *weights){ (*grid)->getQuadratureWeights(weights); }
void tsggiw_(TasmanianSparseGrid **grid, const double *x, double *weights){ (*grid)->getInterpolationWeights(x,weights); }

// set/is/clear/getDomainTransform
void tsgsdt_(TasmanianSparseGrid **grid, const double *transform_a, const double *transform_b){ (*grid)->setDomainTransform(transform_a, transform_b); }
void tsgidt_(TasmanianSparseGrid **grid, int *result){ *result = ((*grid)->isSetDomainTransfrom()) ? 1 : 0; }
void tsgcdt_(TasmanianSparseGrid **grid){ (*grid)->clearDomainTransform(); }
void tsggdt_(TasmanianSparseGrid **grid, double *transform_a, double *transform_b){ (*grid)->getDomainTransform(transform_a, transform_b); }

// loadNeededPoints
void tsglnp_(TasmanianSparseGrid **grid, const double *vals){ (*grid)->loadNeededPoints(vals); }

// evaluate/Fast/Batch/integrate
void tsgeva_(TasmanianSparseGrid **grid, const double *x, double *y){ (*grid)->evaluate(x, y); }
void tsgevf_(TasmanianSparseGrid **grid, const double *x, double *y){ (*grid)->evaluateFast(x, y); }
void tsgevb_(TasmanianSparseGrid **grid, const double *x, int *num_x, double *y){ (*grid)->evaluateBatch(x, (*num_x), y); }
void tsgint_(TasmanianSparseGrid **grid, double *q){ (*grid)->integrate(q); }

// hierarchical functions/coefficients
void tsgehf_(TasmanianSparseGrid **grid, const double *x, int *num_x, double *y){
    (*grid)->evaluateHierarchicalFunctions(x, (*num_x), y);}
void tsgehs_(TasmanianSparseGrid **grid, const double *x, int *num_x, int *pntr, int *indx, double *vals){
    (*grid)->evaluateSparseHierarchicalFunctionsStatic(x, *num_x, pntr, indx, vals);}
void tsgehz_(TasmanianSparseGrid **grid, const double *x, int *num_x, int *num_nz){
    *num_nz = (*grid)->evaluateSparseHierarchicalFunctionsGetNZ(x, *num_x);}
void tsgghc_(TasmanianSparseGrid **grid, double *c){
    const double *cc = (*grid)->getHierarchicalCoefficients();
    (*grid)->isFourier() ? std::copy(cc, cc + 2 * (*grid)->getNumPoints() * (*grid)->getNumOutputs(), c)
                                     : std::copy(cc, cc + (*grid)->getNumPoints() * (*grid)->getNumOutputs(), c);
}
void tsgshc_(TasmanianSparseGrid **grid, double *c){(*grid)->setHierarchicalCoefficients(c);}
void tsghsu_(TasmanianSparseGrid **grid, double *c){
    std::vector<double> sup = (*grid)->getHierarchicalSupport();
    std::copy(sup.begin(), sup.end(), c);
}

// setAnisotropic/Surplus/Refinement
void tsgsar_(TasmanianSparseGrid **grid, int *type, int *min_growth, int *output, int *opt_flags, const int *llimits){
    const int *ll = (opt_flags[0] != 0) ? llimits : nullptr;
    (*grid)->setAnisotropicRefinement(IO::getDepthTypeInt(*type), *min_growth, *output, ll);
}
void tsgeac_(TasmanianSparseGrid **grid, int *type, int *output, int *result){
    auto coeff = (*grid)->estimateAnisotropicCoefficients(IO::getDepthTypeInt(*type), *output);
    int num_coeff = (*grid)->getNumDimensions();
    if ((*type == 2) || (*type == 4) || (*type == 6)) num_coeff *= 2;
    for(int i=0; i<num_coeff; i++) result[i] = coeff[i];
}
void tsgssr_(TasmanianSparseGrid **grid, double *tol, int *output, int *opt_flags, const int *llimits){
    const int *ll;
    ll = (opt_flags[0] != 0) ? llimits : nullptr;
    (*grid)->setSurplusRefinement(*tol, *output, ll); }
void tsgshr_(TasmanianSparseGrid **grid, double *tol, int *type, int *opt_flags, int *output, const int *llimits){
    int theout    = (opt_flags[0] != 0) ? *output : -1;
    const int *ll = (opt_flags[1] != 0) ? llimits : nullptr;

    (*grid)->setSurplusRefinement(*tol, IO::getTypeRefinementInt(*type), theout, ll);
}
void tsgcre_(TasmanianSparseGrid **grid){ (*grid)->clearRefinement(); }
void tsgmre_(TasmanianSparseGrid **grid){ (*grid)->mergeRefinement(); }

// set/is/clear/getConformalTransform
void tsgsca_(TasmanianSparseGrid **grid, int *trunc){ (*grid)->setConformalTransformASIN(Utils::copyArray(trunc, (*grid)->getNumDimensions())); }
void tsgica_(TasmanianSparseGrid **grid, int *result){ *result = ((*grid)->isSetConformalTransformASIN()) ? 1 : 0; }
void tsgcct_(TasmanianSparseGrid **grid){ (*grid)->clearConformalTransform(); }
void tsggca_(TasmanianSparseGrid **grid, int *trunc){
    auto truncation_vector = (*grid)->getConformalTransformASIN();
    std::copy(truncation_vector.begin(), truncation_vector.end(), trunc);
}

// isGlobal/Sequence/LocalPolynomial/Wavelet/Fourier
void tsgisg_(TasmanianSparseGrid **grid, int *status){*status = ((*grid)->isGlobal()) ? 1 : 0;}
void tsgiss_(TasmanianSparseGrid **grid, int *status){*status = ((*grid)->isSequence()) ? 1 : 0;}
void tsgisl_(TasmanianSparseGrid **grid, int *status){*status = ((*grid)->isLocalPolynomial()) ? 1 : 0;}
void tsgisw_(TasmanianSparseGrid **grid, int *status){*status = ((*grid)->isWavelet()) ? 1 : 0;}
void tsgisf_(TasmanianSparseGrid **grid, int *status){*status = ((*grid)->isFourier()) ? 1 : 0;}

// print stats
void tsgpri_(TasmanianSparseGrid **grid){ (*grid)->printStats(); }

// get/enableAcceleration
void tsgacc_(TasmanianSparseGrid **grid, int *acc){
    (*grid)->enableAcceleration(AccelerationMeta::getIOIntAcceleration(*acc));
}
void tsggac_(TasmanianSparseGrid **grid, int *acc){
    TypeAcceleration accel = (*grid)->getAccelerationType();
    *acc = AccelerationMeta::getIOAccelerationInt(accel);
}
void tsgsgi_(TasmanianSparseGrid **grid, int *gpuID){ (*grid)->setGPUID(*gpuID); }
void tsgggi_(TasmanianSparseGrid **grid, int *gpuID){ *gpuID = (*grid)->getGPUID(); }
void tsggng_(int *gpus){ *gpus = TasmanianSparseGrid::getNumGPUs(); }
void tsgggm_(int *gpuID, int *mem){ *mem = TasmanianSparseGrid::getGPUMemory(*gpuID); }

// MPI send/recv/bcast
#ifdef Tasmanian_ENABLE_MPI
void tsgmpigsend_(TasmanianSparseGrid **grid, int *dest, int *tag, MPI_Fint *comm, int *ierr){
    *ierr = MPIGridSend(**grid, *dest, *tag, MPI_Comm_f2c(*comm));
}
void tsgmpigrecv_(TasmanianSparseGrid **grid, int *src, int *tag, MPI_Fint *comm, MPI_Fint *stat, int *ierr){
    MPI_Status status;
    *ierr = MPI_Status_f2c(stat, &status);
    if (*ierr != MPI_SUCCESS) return;
    *ierr = MPIGridRecv(**grid, *src, *tag, MPI_Comm_f2c(*comm), &status);
}
void tsgmpigcast_(TasmanianSparseGrid **grid, int *root, MPI_Fint *comm, int *ierr){
    *ierr = MPIGridBcast(**grid, *root, MPI_Comm_f2c(*comm));
}
#else
void tsgmpigsend_(TasmanianSparseGrid**, int*, int*, int*, int*){}
void tsgmpigrecv_(TasmanianSparseGrid**, int*, int*, int*, int*, int*){}
void tsgmpigcast_(TasmanianSparseGrid**, int*, int*, int*){}
#endif

}
#endif
