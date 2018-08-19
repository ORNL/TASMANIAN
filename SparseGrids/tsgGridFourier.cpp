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

#ifndef __TASMANIAN_SPARSE_GRID_FOURIER_CPP
#define __TASMANIAN_SPARSE_GRID_FOURIER_CPP

#include "tsgGridFourier.hpp"
#include "tsgHiddenExternals.hpp"

namespace TasGrid{

GridFourier::GridFourier() : num_dimensions(0), num_outputs(0), wrapper(0), tensors(0), active_tensors(0), active_w(0),
    max_levels(0), points(0), needed(0), exponents(0), fourier_coefs(0), exponent_refs(0), tensor_refs(0), values(0), accel(0)
{}

GridFourier::GridFourier(const GridFourier &fourier) : num_dimensions(0), num_outputs(0), wrapper(0), tensors(0), active_tensors(0),
    active_w(0), max_levels(0), points(0), needed(0), exponents(0), fourier_coefs(0), exponent_refs(0), tensor_refs(0), values(0), accel(0){
    copyGrid(&fourier);
}

GridFourier::~GridFourier(){ reset(); }

void GridFourier::write(std::ofstream &ofs) const{
    using std::endl;

    ofs << std::scientific; ofs.precision(17);
    ofs << num_dimensions << " " << num_outputs << endl;
    if (num_dimensions > 0){
        tensors->write(ofs);
        active_tensors->write(ofs);
        ofs << active_w[0];
        for(int i=1; i<active_tensors->getNumIndexes(); i++){
            ofs << " " << active_w[i];
        }
        ofs << endl;
        if (points == 0){
            ofs << "0" << endl;
        }else{
            ofs << "1 ";
            points->write(ofs);
        }
        if (needed == 0){
            ofs << "0" << endl;
        }else{
            ofs << "1 ";
            needed->write(ofs);
        }
        ofs << max_levels[0];
        for(int j=1; j<num_dimensions; j++){
            ofs << " " << max_levels[j];
        }
        ofs << endl;
        if (num_outputs > 0){
            values->write(ofs);
            if (fourier_coefs != 0){
                ofs << "1";
                for(int i=0; i < 2*num_outputs*getNumPoints(); i++){
                    ofs << " " << fourier_coefs[i];
                }
            }else{
                ofs << "0";
            }
        }

        /* not needed right now; will need later for refinement
        if (updated_tensors != 0){
            ofs << "1" << endl;
            updated_tensors->write(ofs);
            updated_active_tensors->write(ofs);
            ofs << updated_active_w[0];
            for(int i=1; i<updated_active_tensors->getNumIndexes(); i++){
                ofs << " " << updated_active_w[i];
            }
        }else{
            ofs << "0";
        }
        */

        ofs << endl;
    }
}

void GridFourier::read(std::ifstream &ifs){
    reset();
    ifs >> num_dimensions >> num_outputs;
    if (num_dimensions > 0){
        int flag;

        tensors = new IndexSet(num_dimensions);  tensors->read(ifs);
        active_tensors = new IndexSet(num_dimensions);  active_tensors->read(ifs);
        active_w = new int[active_tensors->getNumIndexes()];  for(int i=0; i<active_tensors->getNumIndexes(); i++){ ifs >> active_w[i]; }
        ifs >> flag; if (flag == 1){ points = new IndexSet(num_dimensions); points->read(ifs); }
        ifs >> flag; if (flag == 1){ needed = new IndexSet(num_dimensions); needed->read(ifs); }
        max_levels.resize(num_dimensions); for(auto &m : max_levels){ ifs >> m; }

        IndexSet *work = (points != 0) ? points : needed;

        if (num_outputs > 0){
            values = new StorageSet(0, 0); values->read(ifs);
            ifs >> flag;
            if (flag == 1){
                fourier_coefs = new double[2 * num_outputs * work->getNumIndexes()];
                for(int i=0; i<2*num_outputs*work->getNumIndexes(); i++){
                    ifs >> fourier_coefs[i];
                }
            }else{
                fourier_coefs = 0;
            }
        }

        IndexManipulator IM(num_dimensions);
        int oned_max_level = max_levels[0];
        int nz_weights = active_tensors->getNumIndexes();
        for(int j=1; j<num_dimensions; j++){ if (oned_max_level < max_levels[j]) oned_max_level = max_levels[j]; }

        OneDimensionalMeta meta(0);
        wrapper = new OneDimensionalWrapper(&meta, oned_max_level, rule_fourier, 0.0, 0.0);

        UnsortedIndexSet *exponents_unsorted = new UnsortedIndexSet(num_dimensions, work->getNumIndexes());
        int *exponent = new int[num_dimensions];

        for (int i=0; i<work->getNumIndexes(); i++){
            for(int j=0; j<work->getNumDimensions(); j++){
                exponent[j] = (work->getIndex(i)[j] % 2 == 0 ? -work->getIndex(i)[j]/2 : (work->getIndex(i)[j]+1)/2);
            }
            exponents_unsorted->addIndex(exponent);
        }

        exponents = new IndexSet(exponents_unsorted);
        delete[] exponent;
        delete exponents_unsorted;

        exponent_refs = new int*[nz_weights];
        tensor_refs = new int*[nz_weights];
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<nz_weights; i++){
            exponent_refs[i] = referenceExponents(active_tensors->getIndex(i), exponents);
            tensor_refs[i] = IM.referenceNestedPoints(active_tensors->getIndex(i), wrapper, work);
        }

        int dummy;
        IM.getMaxLevels(work, max_power, dummy);
    }
}

void GridFourier::writeBinary(std::ofstream &ofs) const{
    int num_dim_out[2];
    num_dim_out[0] = num_dimensions;
    num_dim_out[1] = num_outputs;
    ofs.write((char*) num_dim_out, 2*sizeof(int));
    if (num_dimensions > 0){
        tensors->writeBinary(ofs);
        active_tensors->writeBinary(ofs);
        ofs.write((char*) active_w, active_tensors->getNumIndexes() * sizeof(int));
        char flag;
        if (points == 0){
            flag = 'n'; ofs.write(&flag, sizeof(char));
        }else{
            flag = 'y'; ofs.write(&flag, sizeof(char));
            points->writeBinary(ofs);
        }
        if (needed == 0){
            flag = 'n'; ofs.write(&flag, sizeof(char));
        }else{
            flag = 'y'; ofs.write(&flag, sizeof(char));
            needed->writeBinary(ofs);
        }
        ofs.write((char*) max_levels.data(), num_dimensions * sizeof(int));

        if (num_outputs > 0){
            values->writeBinary(ofs);
            if (fourier_coefs != 0){
                flag = 'y'; ofs.write(&flag, sizeof(char));
                ofs.write((char*) fourier_coefs, 2 * getNumPoints() * num_outputs * sizeof(double));
            }else{
                flag = 'n'; ofs.write(&flag, sizeof(char));
            }
        }

        /* don't need this right now; will need later when refinement is added
        if (updated_tensors != 0){
            flag = 'y'; ofs.write(&flag, sizeof(char));
            updated_tensors->writeBinary(ofs);
            updated_active_tensors->writeBinary(ofs);
            ofs.write((char*) updated_active_w, updated_active_tensors->getNumIndexes() * sizeof(int));
        }else{
            flag = 'n'; ofs.write(&flag, sizeof(char));
        }
        */
    }
}

void GridFourier::readBinary(std::ifstream &ifs){
    reset();
    int num_dim_out[2];
    ifs.read((char*) num_dim_out, 2*sizeof(int));
    num_dimensions = num_dim_out[0];
    num_outputs = num_dim_out[1];

    if (num_dimensions > 0){
        tensors = new IndexSet(num_dimensions);  tensors->readBinary(ifs);
        active_tensors = new IndexSet(num_dimensions);  active_tensors->readBinary(ifs);
        active_w = new int[active_tensors->getNumIndexes()];
        ifs.read((char*) active_w, active_tensors->getNumIndexes() * sizeof(int));

        char flag;
        ifs.read((char*) &flag, sizeof(char)); if (flag == 'y'){ points = new IndexSet(num_dimensions); points->readBinary(ifs); }
        ifs.read((char*) &flag, sizeof(char)); if (flag == 'y'){ needed = new IndexSet(num_dimensions); needed->readBinary(ifs); }

        IndexSet *work = (points != 0) ? points : needed;

        max_levels.resize(num_dimensions);
        ifs.read((char*) max_levels.data(), num_dimensions * sizeof(int));

        if (num_outputs > 0){
            values = new StorageSet(0, 0); values->readBinary(ifs);
            ifs.read((char*) &flag, sizeof(char));
            if (flag == 'y'){
                fourier_coefs = new double[2* num_outputs * work->getNumIndexes()];
                ifs.read((char*) fourier_coefs, 2 * num_outputs * work->getNumIndexes() * sizeof(double));
            }else{
                fourier_coefs = 0;
            }
        }

        IndexManipulator IM(num_dimensions);
        int nz_weights = active_tensors->getNumIndexes();
        int oned_max_level;
        oned_max_level = max_levels[0];
        for(int j=1; j<num_dimensions; j++) if (oned_max_level < max_levels[j]) oned_max_level = max_levels[j];

        OneDimensionalMeta meta(0);
        wrapper = new OneDimensionalWrapper(&meta, oned_max_level, rule_fourier, 0.0, 0.0);

        UnsortedIndexSet* exponents_unsorted = new UnsortedIndexSet(num_dimensions, work->getNumIndexes());
        int *exponent = new int[num_dimensions];

        for (int i=0; i<work->getNumIndexes(); i++){
            for(int j=0; j<work->getNumDimensions(); j++){
                exponent[j] = (work->getIndex(i)[j] % 2 == 0 ? -work->getIndex(i)[j]/2 : (work->getIndex(i)[j]+1)/2);
            }
            exponents_unsorted->addIndex(exponent);
        }
        exponents = new IndexSet(exponents_unsorted);
        delete[] exponent;
        delete exponents_unsorted;

        exponent_refs = new int*[nz_weights];
        tensor_refs = new int*[nz_weights];
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<nz_weights; i++){
            exponent_refs[i] = referenceExponents(active_tensors->getIndex(i), exponents);
            tensor_refs[i] = IM.referenceNestedPoints(active_tensors->getIndex(i), wrapper, work);
        }

        int dummy;
        IM.getMaxLevels(work, max_power, dummy);
    }
}

void GridFourier::reset(){
    clearAccelerationData();
    if (exponent_refs != 0){ for(int i=0; i<active_tensors->getNumIndexes(); i++){ delete[] exponent_refs[i]; exponent_refs[i] = 0; } delete[] exponent_refs; exponent_refs = 0; }
    if (tensor_refs != 0){ for(int i=0; i<active_tensors->getNumIndexes(); i++){ delete[] tensor_refs[i]; tensor_refs[i] = 0; } delete[] tensor_refs; tensor_refs = 0; }
    if (wrapper != 0){ delete wrapper; wrapper = 0; }
    if (tensors != 0){ delete tensors; tensors = 0; }
    if (active_tensors != 0){ delete active_tensors; active_tensors = 0; }
    if (active_w != 0){ delete[] active_w; active_w = 0; }
    if (points != 0){ delete points; points = 0; }
    if (needed != 0){ delete needed; needed = 0; }
    if (exponents != 0){ delete exponents; exponents = 0; }
    if (values != 0){ delete values; values = 0; }
    if (fourier_coefs != 0){ delete[] fourier_coefs; fourier_coefs = 0; }
    num_dimensions = 0;
    num_outputs = 0;
}

void GridFourier::makeGrid(int cnum_dimensions, int cnum_outputs, int depth, TypeDepth type, const std::vector<int> &anisotropic_weights, const std::vector<int> &level_limits){
    IndexManipulator IM(cnum_dimensions);
    IndexSet *tset = IM.selectTensors(depth, type, anisotropic_weights, rule_fourier);
    if (!level_limits.empty()){
        IndexSet *limited = IM.removeIndexesByLimit(tset, level_limits);
        if (limited != 0){
            delete tset;
            tset = limited;
        }
    }

    setTensors(tset, cnum_outputs);
}

void GridFourier::copyGrid(const GridFourier *fourier){
    IndexSet *tset = new IndexSet(fourier->tensors);
    setTensors(tset, fourier->num_outputs);
    if ((num_outputs > 0) && (fourier->points != 0)){ // if there are values inside the source object
        loadNeededPoints(fourier->values->getValues(0));
    }
}

void GridFourier::setTensors(IndexSet* &tset, int cnum_outputs){
    reset();
    num_dimensions = tset->getNumDimensions();
    num_outputs = cnum_outputs;

    tensors = tset;
    tset = 0;

    IndexManipulator IM(num_dimensions);

    OneDimensionalMeta meta(0);
    int max_level; IM.getMaxLevels(tensors, max_levels, max_level);
    wrapper = new OneDimensionalWrapper(&meta, max_level, rule_fourier);

    std::vector<int> tensors_w;
    IM.makeTensorWeights(tensors, tensors_w);
    active_tensors = IM.nonzeroSubset(tensors, tensors_w);

    int nz_weights = active_tensors->getNumIndexes();

    active_w = new int[nz_weights];
    tensor_refs = new int*[nz_weights];
    int count = 0;
    for(int i=0; i<tensors->getNumIndexes(); i++){ if (tensors_w[i] != 0) active_w[count++] = tensors_w[i]; }

    needed = IM.generateNestedPoints(tensors, wrapper); // nested grids exploit nesting

    UnsortedIndexSet* exponents_unsorted = new UnsortedIndexSet(num_dimensions, needed->getNumIndexes());
    int *exponent = new int[num_dimensions];

    for (int i=0; i<needed->getNumIndexes(); i++){
        for(int j=0; j<needed->getNumDimensions(); j++){
            exponent[j] = (needed->getIndex(i)[j] % 2 == 0 ? -needed->getIndex(i)[j]/2 : (needed->getIndex(i)[j]+1)/2);
        }
        exponents_unsorted->addIndex(exponent);
    }

    exponents = new IndexSet(exponents_unsorted);
    delete[] exponent;
    delete exponents_unsorted;

    exponent_refs = new int*[nz_weights];
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<nz_weights; i++){
        exponent_refs[i] = referenceExponents(active_tensors->getIndex(i), exponents);
        tensor_refs[i] = IM.referenceNestedPoints(active_tensors->getIndex(i), wrapper, needed);
    }

    if (num_outputs == 0){
        points = needed;
        needed = 0;
    }else{
        values = new StorageSet(num_outputs, needed->getNumIndexes());
    }

    int dummy;
    IM.getMaxLevels(((points != 0) ? points : needed), max_power, dummy);
}

int* GridFourier::referenceExponents(const int levels[], const IndexSet *ilist){

    /*
     * This function ensures the correct match-up between Fourier coefficients and basis functions.
     * The basis exponents are stored in an IndexSet whose ordering is different from the needed/points
     * IndexSet (since exponents may be negative).
     *
     * The 1D Fourier coefficients are \hat{f}^l_j, where j = -(3^l-1)/2, ..., 0, ..., (3^l-1)/2 and
     *
     * \hat{f}^l_j = \sum_{n=0}^{3^l-1} f(x_n) exp(-2 * pi * sqrt(-1) * j * x_n)
     *             = \sum_{n=0}^{3^l-1} f(x_n) exp(-2 * pi * sqrt(-1) * j * n/3^l)
     *
     * Now, \hat{f}^l_j is the coefficient of the basis function exp(2 * pi * sqrt(-1) * j * x).
     * However, the Fourier transform code returns the coefficients with the indexing j' = 0, 1,
     * ..., 3^l-1.  For j' = 0, 1, ..., (3^l-1)/2, the corresponding basis function has exponent j'.
     * For j' = (3^l-1)/2 + 1, ..., 3^l-1, we march clockwise around the unit circle to see
     * the corresponding exponent is -3^l + j'. Algebraically, for j' > (3^l-1)/2,
     *
     * \hat{f}^l_{j'} = \sum_{n=0}^{3^l-1} f(x_n) [exp(-2 * pi * sqrt(-1) * j' /3^l)]^n
     *                = \sum_{n=0}^{3^l-1} f(x_n) [exp(-2 * pi * sqrt(-1) * (j' - 3^l) /3^l) * exp(-2 * pi * sqrt(-1))]^n
     *                = \sum_{n=0}^{3^l-1} f(x_n) [exp(-2 * pi * sqrt(-1) * (j' - 3^l) /3^l)]^n
     *                = \hat{f}^l_{j' - 3^l}
     */
    int *num_points = new int[num_dimensions];
    int num_total = 1;
    for(int j=0; j<num_dimensions; j++){  num_points[j] = wrapper->getNumPoints(levels[j]); num_total *= num_points[j];  }

    int* refs = new int[num_total];
    int *p = new int[num_dimensions];

    for(int i=0; i<num_total; i++){
        int t = i;
        for(int j=num_dimensions-1; j>=0; j--){
            int tmp = t % num_points[j];
            p[j] = (tmp <= (num_points[j]-1)/2 ? tmp : -num_points[j] + tmp);
            t /= num_points[j];
        }
        refs[i] = ilist->getSlot(p);
    }

    delete[] p;
    delete[] num_points;

    return refs;
}

int GridFourier::getNumDimensions() const{ return num_dimensions; }
int GridFourier::getNumOutputs() const{ return num_outputs; }
TypeOneDRule GridFourier::getRule() const{ return rule_fourier; }

int GridFourier::getNumLoaded() const{ return (((points == 0) || (num_outputs == 0)) ? 0 : points->getNumIndexes()); }
int GridFourier::getNumNeeded() const{ return ((needed == 0) ? 0 : needed->getNumIndexes()); }
int GridFourier::getNumPoints() const{ return ((points == 0) ? getNumNeeded() : points->getNumIndexes()); }

void GridFourier::loadNeededPoints(const double *vals, TypeAcceleration){
    if (accel != 0) accel->resetGPULoadedData();

    if (points == 0){ //setting points for the first time
        values->setValues(vals);
        points = needed;
        needed = 0;
    }else{ //resetting the points
        values->setValues(vals);
    }
    //if we add anisotropic or surplus refinement, I'll need to add a third case here

    int dummy;
    IndexManipulator IM(num_dimensions);
    IM.getMaxLevels(points, max_power, dummy);

    calculateFourierCoefficients();
}

void GridFourier::getLoadedPoints(double *x) const{
    int num_points = points->getNumIndexes();
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = points->getIndex(i);
        for(int j=0; j<num_dimensions; j++){
            x[i*num_dimensions + j] = wrapper->getNode(p[j]);
        }
    }
}
void GridFourier::getNeededPoints(double *x) const{
    int num_points = needed->getNumIndexes();
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = needed->getIndex(i);
        for(int j=0; j<num_dimensions; j++){
            x[i*num_dimensions + j] = wrapper->getNode(p[j]);
        }
    }
}
void GridFourier::getPoints(double *x) const{
    if (points == 0){ getNeededPoints(x); }else{ getLoadedPoints(x); };
}

void GridFourier::generateIndexingMap(std::vector<std::vector<int>> &index_map) const{
    // The internal point-indexing of Tasmanian goes 0, 1/3, 2/3, 1/9, 2/9, 4/9 ....
    // Fourier transform (and coefficients) need spacial order 0, 1/9, 2/9, 3/9=1/3, ...
    // Create a map, where at level 0: 0 -> 0, level 1: 0 1 2 -> 0 1 2, level 2: 0 1 2 3 4 5 6 7 8 -> 0 3 4 1 5 6 2 7 8
    // The map takes a point from previous map and adds two more points ...
    // Thus, a spacial point i on level l is Tasmanian point index_map[l][i]
    IndexManipulator IM(num_dimensions);
    int maxl = IM.getMaxLevel(active_tensors) + 1;
    index_map.resize(maxl);
    index_map[0].resize(1, 0);
    int c = 1;
    for(int l=1; l<maxl; l++){
        index_map[l].resize(3*c); // next level is 3 times the size of the previous
        auto im = index_map[l].begin();
        for(auto i: index_map[l-1]){
            *im++ = i; // point from old level
            *im++ = c++; // two new points
            *im++ = c++;
        }
    }
}

void GridFourier::calculateFourierCoefficients(){
    // There are three indexing schemes needed here
    // First, is the way nodes are indexed and stored in the "work" IndexSet (same as all other nested grids)
    //       the order is contiguous in level, meaning that node indexed by i is always the same node regardless of level
    //       and all nodes on level l are indexed from 0 till 3^l
    // Second, is the spacial order needed by the Fourier transform
    //       nodes have to appear left to right in order for the transform to work so reindexing has to be done
    //       see generateIndexingMap() for index detail
    //       reindexing is done when all data to for the tensor is put into one data structure
    // Third, the exponents of the basis functions have to be indexed and some functions will have negative exponents
    //       following the same idea as the first indexing, the exponents are best used contiguously
    //       reorder the exponents so that 0, 1, 2, 3, 4 ... map to exponents 0, -1, 1, -2, 2, -3, 3 ...
    //       the formula is exponent = (point + 1) / 2, if point is odd, make exponent negative
    // The (point + 1) / 2 maps First indexing to Third indexing, generateIndexingMap() takes care of First -> Second
    // The Second -> Third indexing is done per-tensor, where the non-negative exponents are associated with coefficients
    //     going from left to right (in the order of the Fourier coefficients), while the negative coefficients go
    //     in reverse right to left order "int rj = (p[j] % 2 == 0) ? (p[j]+1) / 2 : num_oned_points[j] - (p[j]+1) / 2;"
    int num_points = getNumPoints();

    IndexSet *work = (points == 0 ? needed : points);
    std::vector<std::vector<int>> index_map;
    generateIndexingMap(index_map);

    if (fourier_coefs != 0){ delete[] fourier_coefs; fourier_coefs = 0; }
    fourier_coefs = new double[2 * num_outputs * num_points];
    std::fill(fourier_coefs, fourier_coefs + 2*num_outputs*num_points, 0.0);

    for(int n=0; n<active_tensors->getNumIndexes(); n++){
        const int* levels = active_tensors->getIndex(n);
        int num_tensor_points = 1;
        std::vector<int> num_oned_points(num_dimensions);
        std::vector<int> cnum_oned_points(num_dimensions, 1); // cumulative number of points
        for(int j=0; j<num_dimensions; j++){
            num_oned_points[j] = wrapper->getNumPoints(levels[j]);
            num_tensor_points *= num_oned_points[j];
        }
        for(int j=num_dimensions-2; j>=0; j--) cnum_oned_points[j] = num_oned_points[j+1] * cnum_oned_points[j+1];

        std::vector<std::vector<std::complex<double>>> tensor_data(num_tensor_points);
        std::vector<int> p(num_dimensions);
        for(int i=0; i<num_tensor_points; i++){
            // We interpret this "i" as running through the spatial indexing; convert to internal
            int t=i;
            for(int j=num_dimensions-1; j>=0; j--){
                // here p[] is the Tasmanian index of index from real space i (t % num_oned_points[j])
                p[j] = index_map[levels[j]][t % num_oned_points[j]];
                t /= num_oned_points[j];
            }
            //refs[i] = ; // refs[i] is the index of Tasmanian (index set) corresponding to real index "i"
            const double *v = values->getValues(work->getSlot(p));
            tensor_data[i].resize(num_outputs);
            std::copy(v, v + num_outputs, tensor_data[i].data());
        }

        TasmanianFourierTransform::fast_fourier_transform(tensor_data, num_oned_points);

        double tensorw = ((double) active_w[n]) / ((double) num_tensor_points);
        for(int i=0; i<num_tensor_points; i++){
            int t = i;
            int r = 0; // holds the real index corresponding to the power
            for(int j=num_dimensions-1; j>=0; j--){
                // here rj is the real multi-index corresponding to "i"
                p[j] = t % num_oned_points[j];
                int rj = (p[j] % 2 == 0) ? (p[j]+1) / 2 : num_oned_points[j] - (p[j]+1) / 2; // +/- index
                r += rj * cnum_oned_points[j];
                t /= num_oned_points[j];
            }
            t = work->getSlot(p); // holds the Tasmanian index corresponding to real index p

            // Combine with tensor weights
            double *fc_real = &(fourier_coefs[num_outputs*t]);
            double *fc_imag = &(fourier_coefs[num_outputs*(t + num_points)]);

             for(auto d : tensor_data[r]){
                 *fc_real++ += tensorw * d.real();
                 *fc_imag++ += tensorw * d.imag();
             }
        }
    }
}

void GridFourier::getInterpolationWeights(const double x[], double weights[]) const {
    /*
    I[f](x) = c^T * \Phi(x) = (U*P*f)^T * \Phi(x)           (U represents normalized forward Fourier transform; P represents reordering of f_i before going into FT)
                            = f^T * (P^T * U^T * \Phi(x))   (P^T = P^(-1) since P is a permutation matrix)

    Note that U is the DFT operator (complex) and the transposes are ONLY REAL transposes, so U^T = U.
    */

    IndexSet *work = (points == 0 ? needed : points);
    std::vector<std::vector<int>> index_map;
    generateIndexingMap(index_map);

    std::fill(weights, weights+getNumPoints(), 0.0);
    double *basisFuncs = new double[2* getNumPoints()];
    computeExponentials<true>(x, basisFuncs);

    if (true){ // use the new fft algorithm

    for(int n=0; n<active_tensors->getNumIndexes(); n++){
        const int *levels = active_tensors->getIndex(n);
        int num_tensor_points = 1;
        std::vector<int> num_oned_points(num_dimensions);

        for(int j=0; j<num_dimensions; j++){
            num_oned_points[j] = wrapper->getNumPoints(levels[j]);
            num_tensor_points *= num_oned_points[j];
        }

        std::vector<std::vector<std::complex<double>>> tensor_data(num_tensor_points);
        auto v = tensor_data.begin();
        for(int i=0; i<num_tensor_points; i++) // v++->resize() means apply resize to the vector that v references and then increment v
            v++->resize(1, std::complex<double> (basisFuncs[2*exponent_refs[n][i]], basisFuncs[2*exponent_refs[n][i] + 1]));

        TasmanianFourierTransform::fast_fourier_transform(tensor_data, num_oned_points);

        std::vector<int> p(num_dimensions);
        v = tensor_data.begin();
        double tensorw = ((double) active_w[n]) / ((double) num_tensor_points);
        for(int i=0; i<num_tensor_points; i++){
            // We interpret this "i" as running through the spatial indexing; convert to internal
            int t=i;
            for(int j=num_dimensions-1; j>=0; j--){
                p[j] = index_map[levels[j]][t % num_oned_points[j]];
                t /= num_oned_points[j];
            }
            weights[work->getSlot(p)] += tensorw * v++->data()->real();
        }
    }

    }else{ // use the old dft algorithm

    for(int n=0; n<active_tensors->getNumIndexes(); n++){
        const int *levels = active_tensors->getIndex(n);
        int num_tensor_points = 1;
        std::vector<int> num_oned_points(num_dimensions);

        for(int j=0; j<num_dimensions; j++){
            num_oned_points[j] = wrapper->getNumPoints(levels[j]);
            num_tensor_points *= num_oned_points[j];
        }

        std::complex<double> *in = new std::complex<double>[num_tensor_points];
        std::complex<double> *out = new std::complex<double>[num_tensor_points];

        for(int i=0; i<num_tensor_points; i++){
            in[i] = std::complex<double>(basisFuncs[2*exponent_refs[n][i]], basisFuncs[2*exponent_refs[n][i] + 1]);
        }

        TasmanianFourierTransform::discrete_fourier_transform(num_dimensions, num_oned_points.data(), in, out);

        std::vector<int> p(num_dimensions);
        for(int i=0; i<num_tensor_points; i++){
            // We interpret this "i" as running through the spatial indexing; convert to internal
            int t=i;
            for(int j=num_dimensions-1; j>=0; j--){
                p[j] = index_map[levels[j]][t % num_oned_points[j]];
                t /= num_oned_points[j];
            }
            weights[work->getSlot(p)] += ((double) active_w[n]) * out[i].real()/((double) num_tensor_points);
        }
        delete[] in;
        delete[] out;
    }

    }
    delete[] basisFuncs;
}

void GridFourier::getQuadratureWeights(double weights[]) const{

    /*
     * When integrating the Fourier series on a tensored grid, all the
     * nonzero modes vanish, and we're left with the normalized Fourier
     * coeff for e^0 (sum of the data divided by number of points)
     */

    int num_points = getNumPoints();
    std::fill(weights, weights+num_points, 0.0);
    for(int n=0; n<active_tensors->getNumIndexes(); n++){
        const int *levels = active_tensors->getIndex(n);
        int num_tensor_points = 1;
        for(int j=0; j<num_dimensions; j++){
            num_tensor_points *= wrapper->getNumPoints(levels[j]);
        }
        for(int i=0; i<num_tensor_points; i++){
            weights[tensor_refs[n][i]] += ((double) active_w[n])/((double) num_tensor_points);
        }
    }
}

void GridFourier::evaluate(const double x[], double y[]) const{
    int num_points = points->getNumIndexes();
    TasBLAS::setzero(num_outputs, y);
    std::vector<double> wreal(num_points);
    std::vector<double> wimag(num_points);
    computeBasis<double, false>(points, x, wreal.data(), wimag.data());
    for(int i=0; i<num_points; i++){
        const double *fcreal = &(fourier_coefs[i*num_outputs]);
        const double *fcimag = &(fourier_coefs[(i+num_points) * num_outputs]);
        double wr = wreal[i];
        double wi = wimag[i];
        for(int k=0; k<num_outputs; k++) y[k] += wr * fcreal[k] - wi * fcimag[k];
    }
}
void GridFourier::evaluateBatch(const double x[], int num_x, double y[]) const{
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        evaluate(&(x[((size_t) i) * ((size_t) num_dimensions)]), &(y[((size_t) i) * ((size_t) num_outputs)]));
    }
}

void GridFourier::evaluateFastCPUblas(const double x[], double y[]) const{
    evaluate(x,y);
}
void GridFourier::evaluateFastGPUcublas(const double x[], double y[]) const{
    evaluate(x,y);
}
void GridFourier::evaluateFastGPUcuda(const double x[], double y[]) const{
    evaluate(x,y);
}
void GridFourier::evaluateFastGPUmagma(int, const double x[], double y[]) const{
    evaluate(x,y);
}

void GridFourier::evaluateBatchCPUblas(const double x[], int num_x, double y[]) const{
    evaluateBatch(x, num_x, y);
}
void GridFourier::evaluateBatchGPUcublas(const double x[], int num_x, double y[]) const {
    evaluateBatch(x, num_x, y);
}
void GridFourier::evaluateBatchGPUcuda(const double x[], int num_x, double y[]) const {
    evaluateBatch(x, num_x, y);
}
void GridFourier::evaluateBatchGPUmagma(int, const double x[], int num_x, double y[]) const {
    evaluateBatch(x, num_x, y);
}

void GridFourier::integrate(double q[], double *conformal_correction) const{
    std::fill(q, q+num_outputs, 0.0);
    if (conformal_correction == 0){
        // everything vanishes except the Fourier coeff of e^0
        int *zeros_num_dim = new int[num_dimensions];
        std::fill(zeros_num_dim, zeros_num_dim + num_dimensions, 0);
        int idx = exponents->getSlot(zeros_num_dim);
        for(int k=0; k<num_outputs; k++){
            q[k] = fourier_coefs[num_outputs * idx + k];
        }
        delete[] zeros_num_dim;
    }else{
        // Do the expensive computation if we have a conformal map
        double *w = new double[getNumPoints()];
        getQuadratureWeights(w);
        for(int i=0; i<points->getNumIndexes(); i++){
            w[i] *= conformal_correction[i];
            const double *v = values->getValues(i);
            for(int k=0; k<num_outputs; k++){
                q[k] += w[i] * v[k];
            }
        }
        delete[] w;
    }
}

void GridFourier::evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const{
    // y must be of size 2*num_x*num_points
    int num_points = getNumPoints();
    Data2D<double> xx; xx.cload(num_dimensions, num_x, x);
    Data2D<double> yy; yy.load(2 * num_points, num_x, y);
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        computeBasis<double, true>(((points == 0) ? needed : points), xx.getCStrip(i), yy.getStrip(i), 0);
    }
}
void GridFourier::evaluateHierarchicalFunctionsInternal(const double x[], int num_x, double M_real[], double M_imag[]) const{
    // y must be of size num_x * num_nodes * 2
    int num_points = getNumPoints();
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        double *w = new double[2 * num_points];
        computeExponentials<false>(&(x[((size_t) i) * ((size_t) num_dimensions)]), w);
        for(int m=0; m<num_points; m++){
            M_real[i*num_points+m] = w[m];
            M_imag[i*num_points+m] = w[m + num_points];
        }
    }
}

void GridFourier::setHierarchicalCoefficients(const double c[], TypeAcceleration){
    // takes c to be length 2*num_outputs*num_points
    // first num_points*num_outputs are the real part; second num_points*num_outputs are the imaginary part

    if (accel != 0) accel->resetGPULoadedData();
    if (points == 0){
        points = needed;
        needed = 0;
    }
    if (fourier_coefs != 0) delete[] fourier_coefs;
    fourier_coefs = new double[2 * getNumPoints() * num_outputs];
    for(int i=0; i<2 * getNumPoints() * num_outputs; i++) fourier_coefs[i] = c[i];
}

void GridFourier::clearAccelerationData(){
    if (accel != 0){
        delete accel;
        accel = 0;
    }
}
void GridFourier::clearRefinement(){ return; }     // to be expanded later
void GridFourier::mergeRefinement(){ return; }     // to be expanded later

const int* GridFourier::getPointIndexes() const{
    return ((points == 0) ? needed->getIndex(0) : points->getIndex(0));
}
const IndexSet* GridFourier::getExponents() const{
    return exponents;
}
const double* GridFourier::getFourierCoefs() const{
    return ((double*) fourier_coefs);
}

} // end TasGrid

#endif
