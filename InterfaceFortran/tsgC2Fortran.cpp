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

extern "C" void tsgc2fint_(int *i);
extern "C" void tsgc2fdouble_(double *i);
extern "C" void tsgc2fvec_(int *length, double *vect);
extern "C" void tsgc2fmat_(int *rows, int *cols, double *mat);


TasmanianSparseGrid **_tsg_grid_list;
int _num_grids;

extern "C"{

void tsgbeg_(){
    _num_grids = 5; // assume we are working with 5 grids
    _tsg_grid_list = new TasmanianSparseGrid*[_num_grids];
    for(int i=0; i<_num_grids; i++) _tsg_grid_list[i] = 0;
}
void tsgnew_(){
    int id = 0;
    while((id < _num_grids) && (_tsg_grid_list[id] != 0)) id++;
    if (id == _num_grids){
        TasmanianSparseGrid **_tsg_grid_list_old = _tsg_grid_list;
        _num_grids += 5;
        _tsg_grid_list = new TasmanianSparseGrid*[_num_grids];
        for(int i=0; i<_num_grids-5; i++) _tsg_grid_list[i] = _tsg_grid_list_old[i];
        for(int i=_num_grids-5; i<_num_grids; i++) _tsg_grid_list_old[i] = 0;
        delete[] _tsg_grid_list_old;
        _tsg_grid_list_old = 0;
    }
    //cout << "Creating grid: " << id << endl;
    _tsg_grid_list[id] = new TasmanianSparseGrid();
    tsgc2fint_(&id);
}
void tsgfre_(int *id){
    if (*id < _num_grids){
        if (_tsg_grid_list[*id] != 0){
            delete _tsg_grid_list[*id];
            _tsg_grid_list[*id] = 0;
        }
    }
}
void tsgend_(){
    for(int i=0; i<_num_grids; i++){
        if (_tsg_grid_list[i] != 0) delete _tsg_grid_list[i];
    }
    delete[] _tsg_grid_list;
    _tsg_grid_list = 0;
    _num_grids = 0;
}
////////////////////////////////////////////////////////////////////////
// MAIN INTERFACE
////////////////////////////////////////////////////////////////////////

void tsggvm_(){
    int ver = TasmanianSparseGrid::getVersionMajor(); 
    tsgc2fint_(&ver);
}
void tsggvn_(){ 
    int ver = TasmanianSparseGrid::getVersionMinor(); 
    tsgc2fint_(&ver);
}

void tsgcer_(int id){ _tsg_grid_list[id]->setErrorLog(&cerr); }

// read/write

// create
void tsgmg_(int *id, int *dimensions, int *outputs, int *depth, int *type, int *rule, const int *anisotropic_weights, double *alpha, double *beta){
    _tsg_grid_list[*id]->makeGlobalGrid(*dimensions, *outputs, *depth,
        OneDimensionalMeta::getIOTypeInt(*type), OneDimensionalMeta::getIORuleInt(*rule), 
        anisotropic_weights, *alpha, *beta);
}
void tsgms_(int *id, int *dimensions, int *outputs, int *depth, int *type, int *rule, const int *anisotropic_weights){
    _tsg_grid_list[*id]->makeSequenceGrid(*dimensions, *outputs, *depth,
        OneDimensionalMeta::getIOTypeInt(*type), OneDimensionalMeta::getIORuleInt(*rule), anisotropic_weights);
}
void tsgml_(int *id, int *dimensions, int *outputs, int *depth, int *order, int *rule){
    _tsg_grid_list[*id]->makeLocalPolynomialGrid(*dimensions, *outputs, *depth, OneDimensionalMeta::getIORuleInt(*rule));
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
void tsggal_(int *id){
    double alpha = _tsg_grid_list[*id]->getAlpha();
    tsgc2fdouble_(&alpha);
}
void tsggbe_(int *id){ 
    double beta = _tsg_grid_list[*id]->getBeta(); 
    tsgc2fdouble_(&beta);
}
void tsggor_(int *id){ 
    int order = _tsg_grid_list[*id]->getOrder(); 
    tsgc2fint_(&order);
}
void tsggnd_(int *id){ 
    int dims = _tsg_grid_list[*id]->getNumDimensions(); 
    tsgc2fint_(&dims);
}
void tsggno_(int *id){
    int outs = _tsg_grid_list[*id]->getNumOutputs(); 
    tsgc2fint_(&outs);
}
void tsggru_(int *id){
    int rule = OneDimensionalMeta::getIORuleInt(_tsg_grid_list[*id]->getRule()); 
    tsgc2fint_(&rule);
}

// getNumNeeded/Loaded/Points
void tsggnn_(int *id){ 
    int num_points = _tsg_grid_list[*id]->getNumNeeded(); 
    tsgc2fint_(&num_points);
}
void tsggnl_(int *id){ 
    int num_points = _tsg_grid_list[*id]->getNumLoaded(); 
    tsgc2fint_(&num_points);
}
void tsggnp_(int *id, int *num_points){
    *num_points = _tsg_grid_list[*id]->getNumPoints();
    //tsgc2fint_(&num_points);
}

// getLoaded/Needed/Points
void tsgglp_(int *id){
    int rows = _tsg_grid_list[*id]->getNumDimensions();
    int cols = _tsg_grid_list[*id]->getNumPoints();
    double *points = (double*) malloc(rows * cols * sizeof(double));
    _tsg_grid_list[*id]->getLoadedPoints(points);
    tsgc2fmat_(&rows, &cols, points);
    points = 0;
}
void tsggls_(int *id, double *points){ _tsg_grid_list[*id]->getLoadedPoints(points); }
void tsggdp_(int *id){
    int rows = _tsg_grid_list[*id]->getNumDimensions();
    int cols = _tsg_grid_list[*id]->getNumPoints();
    double *points = (double*) malloc(rows * cols * sizeof(double));
    _tsg_grid_list[*id]->getNeededPoints(points);
    tsgc2fmat_(&rows, &cols, points);
    points = 0;
}
void tsggds_(int *id, double *points){ _tsg_grid_list[*id]->getNeededPoints(points); }
void tsggpp_(int *id){
    int rows = _tsg_grid_list[*id]->getNumDimensions();
    int cols = _tsg_grid_list[*id]->getNumPoints();
    double *points = (double*) malloc(rows * cols * sizeof(double));
    _tsg_grid_list[*id]->getPoints(points);
    tsgc2fmat_(&rows, &cols, points);
    points = 0;
}
void tsggps_(int *id, double *points){ _tsg_grid_list[*id]->getPoints(points); }

// getQuadrature/InterpolationWeights
void tsggqw_(int *id){
    int num_points = _tsg_grid_list[*id]->getNumPoints();
    double *weights = (double*) malloc(num_points * sizeof(double));
    _tsg_grid_list[*id]->getQuadratureWeights(weights);
    tsgc2fvec_(&num_points, weights);
    weights = 0;
}
void tsggqs_(int *id, double *weights){ _tsg_grid_list[*id]->getQuadratureWeights(weights); }
void tsggiw_(int *id){
    int num_points = _tsg_grid_list[*id]->getNumPoints();
    double *weights = (double*) malloc(num_points * sizeof(double));
    _tsg_grid_list[*id]->getQuadratureWeights(weights);
    tsgc2fvec_(&num_points, weights);
    weights = 0;
}
void tsggis_(int *id, double *weights){ _tsg_grid_list[*id]->getQuadratureWeights(weights); }

}

#endif

