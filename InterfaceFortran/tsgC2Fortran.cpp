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

#include "TasmanianSparseGrid.hpp"

using namespace TasGrid;

//extern "C" void tsgc2fint_(int *i);
//extern "C" void tsgc2fdouble_(double *i);
extern "C" void tsgc2fstr_(int *length, const char *str);
extern "C" void tsgc2fvec_(int *length, double *vect);
extern "C" void tsgc2fmat_(int *rows, int *cols, double *mat);


TasmanianSparseGrid **_tsg_grid_list;
int _tsg_num_grids;

extern "C"{

void tsgbeg_(){
    _tsg_num_grids = 4; // assume we are working with 4 grids
    _tsg_grid_list = new TasmanianSparseGrid*[_tsg_num_grids];
    for(int i=0; i<_tsg_num_grids; i++) _tsg_grid_list[i] = 0;
}
void tsgnew_(int *returnID){
    int id = 0;
    while((id < _tsg_num_grids) && (_tsg_grid_list[id] != 0)) id++;
    if (id == _tsg_num_grids){
        TasmanianSparseGrid **_tsg_grid_list_old = _tsg_grid_list;
        _tsg_num_grids *= 2;
        _tsg_grid_list = new TasmanianSparseGrid*[_tsg_num_grids];
        for(int i=0; i<_tsg_num_grids/2; i++) _tsg_grid_list[i] = _tsg_grid_list_old[i];
        for(int i=_tsg_num_grids/2; i<_tsg_num_grids; i++) _tsg_grid_list_old[i] = 0;
        delete[] _tsg_grid_list_old;
        _tsg_grid_list_old = 0;
    }
    _tsg_grid_list[id] = new TasmanianSparseGrid();
    *returnID = id;
}
void tsgfre_(int *id){
    if (*id < _tsg_num_grids){
        if (_tsg_grid_list[*id] != 0){
            delete _tsg_grid_list[*id];
            _tsg_grid_list[*id] = 0;
        }
    }
}
void tsgend_(){
    for(int i=0; i<_tsg_num_grids; i++){
        if (_tsg_grid_list[i] != 0) delete _tsg_grid_list[i];
    }
    delete[] _tsg_grid_list;
    _tsg_grid_list = 0;
    _tsg_num_grids = 0;
}
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

void tsgcer_(int *id){ _tsg_grid_list[*id]->setErrorLog(&cerr); }

// read/write
void tsgwri_(int *id, const char *filename, int *binary){ _tsg_grid_list[*id]->write(filename, (*binary != 0)); }
void tsgrea_(int *id, const char *filename){ _tsg_grid_list[*id]->read(filename); }

// create
void tsgmg_(int *id, int *dimensions, int *outputs, int *depth, int *type, int *rule, const int *anisotropic_weights, double *alpha, double *beta, const int *llimits){
    _tsg_grid_list[*id]->makeGlobalGrid(*dimensions, *outputs, *depth,
        OneDimensionalMeta::getIOTypeInt(*type), OneDimensionalMeta::getIORuleInt(*rule),
        anisotropic_weights, *alpha, *beta, 0, llimits);
}
void tsgms_(int *id, int *dimensions, int *outputs, int *depth, int *type, int *rule, const int *anisotropic_weights){
    _tsg_grid_list[*id]->makeSequenceGrid(*dimensions, *outputs, *depth,
        OneDimensionalMeta::getIOTypeInt(*type), OneDimensionalMeta::getIORuleInt(*rule), anisotropic_weights);
}
void tsgml_(int *id, int *dimensions, int *outputs, int *depth, int *order, int *rule){
    _tsg_grid_list[*id]->makeLocalPolynomialGrid(*dimensions, *outputs, *depth, *order, OneDimensionalMeta::getIORuleInt(*rule));
}
void tsgmw_(int *id, int *dimensions, int *outputs, int *depth, int *order){
    _tsg_grid_list[*id]->makeWaveletGrid(*dimensions, *outputs, *depth, *order);
}

// copy/updare
void tsgcp_(int *id, int *source){ _tsg_grid_list[*id]->copyGrid(_tsg_grid_list[*source]); }

void tsgug_(int *id, int *depth, int *type, const int *anisotropic_weights){
    _tsg_grid_list[*id]->updateGlobalGrid(*depth, OneDimensionalMeta::getIOTypeInt(*type), anisotropic_weights);
}
void tsgus_(int *id, int *depth, int *type, const int *anisotropic_weights){
    _tsg_grid_list[*id]->updateSequenceGrid(*depth, OneDimensionalMeta::getIOTypeInt(*type), anisotropic_weights);
}

// getAlpha/Beta/Order/Dims/Outs/Rule //
void tsggal_(int *id, double *alpha){ *alpha = _tsg_grid_list[*id]->getAlpha(); }
void tsggbe_(int *id, double *beta){ *beta = _tsg_grid_list[*id]->getBeta(); }
void tsggor_(int *id, int *order){ *order = _tsg_grid_list[*id]->getOrder(); }
void tsggnd_(int *id, int *dims){ *dims = _tsg_grid_list[*id]->getNumDimensions(); }
void tsggno_(int *id, int *outs){ *outs = _tsg_grid_list[*id]->getNumOutputs(); }
void tsggru_(int *id, int *rule){ *rule = OneDimensionalMeta::getIORuleInt(_tsg_grid_list[*id]->getRule()); }

// getNumNeeded/Loaded/Points
void tsggnn_(int *id, int *num_points){ *num_points = _tsg_grid_list[*id]->getNumNeeded(); }
void tsggnl_(int *id, int *num_points){ *num_points = _tsg_grid_list[*id]->getNumLoaded(); }
void tsggnp_(int *id, int *num_points){ *num_points = _tsg_grid_list[*id]->getNumPoints(); }

// getLoaded/Needed/Points
void tsgglp_(int *id, double *points){ _tsg_grid_list[*id]->getLoadedPoints(points); }
void tsggdp_(int *id, double *points){ _tsg_grid_list[*id]->getNeededPoints(points); }
void tsggpp_(int *id, double *points){ _tsg_grid_list[*id]->getPoints(points); }

// getQuadrature/InterpolationWeights
void tsggqw_(int *id, double *weights){ _tsg_grid_list[*id]->getQuadratureWeights(weights); }
void tsggiw_(int *id, double *weights){ _tsg_grid_list[*id]->getQuadratureWeights(weights); }

// set/is/clear/getDomainTransform
void tsgsdt_(int *id, const double *transform_a, const double *transform_b){ _tsg_grid_list[*id]->setDomainTransform(transform_a, transform_b); }
void tsgidt_(int *id, int *result){ *result = (_tsg_grid_list[*id]->isSetDomainTransfrom()) ? 1 : 0; }
void tsgcdt_(int *id){ _tsg_grid_list[*id]->clearDomainTransform(); }
void tsggdt_(int *id, double *transform_a, double *transform_b){ _tsg_grid_list[*id]->getDomainTransform(transform_a, transform_b); }

// loadNeededPoints
void tsglnp_(int *id, const double *vals){ _tsg_grid_list[*id]->loadNeededPoints(vals); }

// evaluate/Fast/Batch/integrate
void tsgeva_(int *id, const double *x, double *y){ _tsg_grid_list[*id]->evaluate(x, y); }
void tsgevf_(int *id, const double *x, double *y){ _tsg_grid_list[*id]->evaluateFast(x, y); }
void tsgevb_(int *id, const double *x, int *num_x, double *y){ _tsg_grid_list[*id]->evaluateBatch(x, (*num_x), y); }
void tsgint_(int *id, double *q){ _tsg_grid_list[*id]->integrate(q); }

// setAnisotropic/Surplus/Refinement
void tsgsar_(int *id, int *type, int *min_growth, int *output){
    _tsg_grid_list[*id]->setAnisotropicRefinement(OneDimensionalMeta::getIOTypeInt(*type), *min_growth, *output);
}
void tsgeac_(int *id, int *type, int *output, int *result){
    int *coeff = _tsg_grid_list[*id]->estimateAnisotropicCoefficients(OneDimensionalMeta::getIOTypeInt(*type), *output);
    int num_coeff = _tsg_grid_list[*id]->getNumDimensions();
    if ((*type == 2) || (*type == 4) || (*type == 6)) num_coeff *= 2;
    for(int i=0; i<num_coeff; i++) result[i] = coeff[i];
    delete[] coeff;
}
void tsgssr_(int *id, double *tol, int *output){ _tsg_grid_list[*id]->setSurplusRefinement(*tol, *output, 0); }
void tsgshr_(int *id, double *tol, int *type, int *output){
    _tsg_grid_list[*id]->setSurplusRefinement(*tol, OneDimensionalMeta::getIOTypeRefinementInt(*type), *output);
}
void tsgcre_(int *id){ _tsg_grid_list[*id]->clearRefinement(); }

// set/is/clear/getConformalTransform
void tsgsca_(int *id, int *trunc){ _tsg_grid_list[*id]->setConformalTransformASIN(trunc); }
void tsgica_(int *id, int *result){ *result = (_tsg_grid_list[*id]->isSetConformalTransformASIN()) ? 1 : 0; }
void tsgcct_(int *id){ _tsg_grid_list[*id]->clearConformalTransform(); }
void tsggca_(int *id, int *trunc){ _tsg_grid_list[*id]->getConformalTransformASIN(trunc); }

// print stats
void tsgpri_(int *id){ _tsg_grid_list[*id]->printStats(); }

// get/enableAcceleration
void tsgacc_(int *id, int *acc){
    TypeAcceleration accel = accel_none;
    if (*acc == 1){ accel = accel_cpu_blas;   }else
    if (*acc == 2){ accel = accel_gpu_cublas; }else
    if (*acc == 3){ accel = accel_gpu_cuda;   }else
    if (*acc == 4){ accel = accel_gpu_magma;  }
    _tsg_grid_list[*id]->enableAcceleration(accel);
}
void tsggac_(int *id, int *acc){
    TypeAcceleration accel = _tsg_grid_list[*id]->getAccelerationType();
    if (accel == accel_none){       *acc = 0; }else
    if (accel == accel_cpu_blas){   *acc = 1; }else
    if (accel == accel_gpu_cublas){ *acc = 2; }else
    if (accel == accel_gpu_cuda){   *acc = 3; }else
    if (accel == accel_gpu_magma){  *acc = 4; }
}
void tsgsgi_(int *id, int *gpuID){ _tsg_grid_list[*id]->setGPUID(*gpuID); }
void tsgggi_(int *id, int *gpuID){ *gpuID = _tsg_grid_list[*id]->getGPUID(); }
void tsggng_(int *gpus){ *gpus = TasmanianSparseGrid::getNumGPUs(); }
void tsgggm_(int *gpuID, int *mem){ *mem = TasmanianSparseGrid::getGPUMemory(*gpuID); }

}
#endif
