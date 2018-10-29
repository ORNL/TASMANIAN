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

GridFourier::GridFourier() : num_dimensions(0), num_outputs(0), wrapper(0), tensors(0), active_tensors(0),
    max_levels(0), points(0), needed(0)
{}

GridFourier::GridFourier(const GridFourier &fourier) : num_dimensions(0), num_outputs(0), wrapper(0), tensors(0), active_tensors(0),
    max_levels(0), points(0), needed(0){
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
            values.write(ofs);
            if (fourier_coefs.getNumStrips() > 0){
                ofs << "1";
                for(auto c : *fourier_coefs.getVector()) ofs << " " << c;
            }else{
                ofs << "0";
            }
        }
        ofs << endl;
        ofs << "0"; // indicate no refinement, will implement in the future

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
        active_w.resize(active_tensors->getNumIndexes());  for(auto &w : active_w) ifs >> w;
        ifs >> flag; if (flag == 1){ points = new IndexSet(num_dimensions); points->read(ifs); }
        ifs >> flag; if (flag == 1){ needed = new IndexSet(num_dimensions); needed->read(ifs); }
        max_levels.resize(num_dimensions); for(auto &m : max_levels){ ifs >> m; }

        IndexSet *work = (points != 0) ? points : needed;

        if (num_outputs > 0){
            values.read(ifs);
            ifs >> flag;
            if (flag == 1){
                fourier_coefs.resize(num_outputs, 2 * work->getNumIndexes());
                for(auto &c : *fourier_coefs.getVector()) ifs >> c;
            }
        }
        ifs >> flag;

        IndexManipulator IM(num_dimensions);
        int oned_max_level = max_levels[0];
        for(int j=1; j<num_dimensions; j++){ if (oned_max_level < max_levels[j]) oned_max_level = max_levels[j]; }

        OneDimensionalMeta meta(0);
        wrapper = new OneDimensionalWrapper(&meta, oned_max_level, rule_fourier, 0.0, 0.0);

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
        ofs.write((char*) active_w.data(), active_tensors->getNumIndexes() * sizeof(int));
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
            values.writeBinary(ofs);
            if (fourier_coefs.getNumStrips() > 0){
                flag = 'y'; ofs.write(&flag, sizeof(char));
                ofs.write((char*) fourier_coefs.getCStrip(0), 2 * ((size_t) getNumPoints()) * ((size_t) num_outputs) * sizeof(double));
            }else{
                flag = 'n'; ofs.write(&flag, sizeof(char));
            }
        }
        flag = 'n'; ofs.write(&flag, sizeof(char)); // no refinement, will implement in the future

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
        active_w.resize(active_tensors->getNumIndexes());
        ifs.read((char*) active_w.data(), active_tensors->getNumIndexes() * sizeof(int));

        char flag;
        ifs.read((char*) &flag, sizeof(char)); if (flag == 'y'){ points = new IndexSet(num_dimensions); points->readBinary(ifs); }
        ifs.read((char*) &flag, sizeof(char)); if (flag == 'y'){ needed = new IndexSet(num_dimensions); needed->readBinary(ifs); }

        IndexSet *work = (points != 0) ? points : needed;

        max_levels.resize(num_dimensions);
        ifs.read((char*) max_levels.data(), num_dimensions * sizeof(int));

        if (num_outputs > 0){
            values.readBinary(ifs);
            ifs.read((char*) &flag, sizeof(char));
            if (flag == 'y'){
                fourier_coefs.resize(num_outputs, 2 * work->getNumIndexes());
                ifs.read((char*) fourier_coefs.getStrip(0), 2 * ((size_t) num_outputs) * ((size_t) work->getNumIndexes()) * sizeof(double));
            }
        }
        ifs.read((char*) &flag, sizeof(char));

        IndexManipulator IM(num_dimensions);
        int oned_max_level;
        oned_max_level = max_levels[0];
        for(int j=1; j<num_dimensions; j++) if (oned_max_level < max_levels[j]) oned_max_level = max_levels[j];

        OneDimensionalMeta meta(0);
        wrapper = new OneDimensionalWrapper(&meta, oned_max_level, rule_fourier, 0.0, 0.0);

        int dummy;
        IM.getMaxLevels(work, max_power, dummy);
    }
}

void GridFourier::reset(){
    clearAccelerationData();
    if (wrapper != 0){ delete wrapper; wrapper = 0; }
    if (tensors != 0){ delete tensors; tensors = 0; }
    if (active_tensors != 0){ delete active_tensors; active_tensors = 0; }
    if (points != 0){ delete points; points = 0; }
    if (needed != 0){ delete needed; needed = 0; }
    values = StorageSet();
    active_w.clear();
    fourier_coefs.clear();
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
        loadNeededPoints(fourier->values.getValues(0));
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
    wrapper = new OneDimensionalWrapper(&meta, max_level, rule_fourier, 0.0, 0.0);

    std::vector<int> tensors_w;
    IM.makeTensorWeights(tensors, tensors_w);
    active_tensors = IM.nonzeroSubset(tensors, tensors_w);

    int nz_weights = active_tensors->getNumIndexes();

    active_w.reserve(nz_weights);
    for(int i=0; i<tensors->getNumIndexes(); i++){ if (tensors_w[i] != 0) active_w.push_back(tensors_w[i]); }

    needed = IM.generateNestedPoints(tensors, wrapper); // nested grids exploit nesting

    if (num_outputs == 0){
        points = needed;
        needed = 0;
    }else{
        values.resize(num_outputs, needed->getNumIndexes());
    }

    int dummy;
    IM.getMaxLevels(((points != 0) ? points : needed), max_power, dummy);
}

int GridFourier::getNumDimensions() const{ return num_dimensions; }
int GridFourier::getNumOutputs() const{ return num_outputs; }
TypeOneDRule GridFourier::getRule() const{ return rule_fourier; }

int GridFourier::getNumLoaded() const{ return (((points == 0) || (num_outputs == 0)) ? 0 : points->getNumIndexes()); }
int GridFourier::getNumNeeded() const{ return ((needed == 0) ? 0 : needed->getNumIndexes()); }
int GridFourier::getNumPoints() const{ return ((points == 0) ? getNumNeeded() : points->getNumIndexes()); }

void GridFourier::loadNeededPoints(const double *vals, TypeAcceleration){
    clearAccelerationData();
    if (points == 0){ //setting points for the first time
        values.setValues(vals);
        points = needed;
        needed = 0;
    }else{ //resetting the points
        values.setValues(vals);
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

    fourier_coefs.resize(num_outputs, 2 * num_points, 0.0);

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
            const double *v = values.getValues(work->getSlot(p));
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
            double *fc_real = fourier_coefs.getStrip(t);
            double *fc_imag = fourier_coefs.getStrip(t + num_points);

             for(auto d : tensor_data[r]){
                 *fc_real++ += tensorw * d.real();
                 *fc_imag++ += tensorw * d.imag();
             }
        }
    }
}

void GridFourier::getInterpolationWeights(const double x[], double weights[]) const {
    // if Fourier coefficient are c, Data from the target function is f, and values of the basis functions are b
    // then we have c = A * f (where A is both the Fourier transform and the reindexing)
    // and the value of the interpolant is result = <c, b> = <A f, b> = <f, A^* b>, where A^* is conjugate transpose
    // However, we consider only the real values (the complex ones add-up to zero), thus we really need A^T (regular transpose)
    // The Fourier transform is symmetric with respect to regular transpose, which leaves only the indexing
    // Take the basis functions, reindex and reorder to a data strucutre, take FFT, reindex and reorder into the weights

    IndexSet *work = (points == 0 ? needed : points);
    std::vector<std::vector<int>> index_map;
    generateIndexingMap(index_map);

    std::fill(weights, weights+getNumPoints(), 0.0);

    std::vector<std::complex<double>> basisFuncs(work->getNumIndexes());
    computeBasis<double, true>(work, x, (double*) basisFuncs.data(), 0);

    for(int n=0; n<active_tensors->getNumIndexes(); n++){
        const int *levels = active_tensors->getIndex(n);
        int num_tensor_points = 1;
        std::vector<int> num_oned_points(num_dimensions);
        std::vector<int> cnum_oned_points(num_dimensions, 1);
        for(int j=0; j<num_dimensions; j++){
            num_oned_points[j] = wrapper->getNumPoints(levels[j]);
            num_tensor_points *= num_oned_points[j];
        }
        for(int j=num_dimensions-2; j>=0; j--) cnum_oned_points[j] = num_oned_points[j+1] * cnum_oned_points[j+1];

        std::vector<std::vector<std::complex<double>>> tensor_data(num_tensor_points);
        std::vector<int> p(num_dimensions);
        for(int i=0; i<num_tensor_points; i++){ // v++->resize() means apply resize to the vector that v references and then increment v
            int r = 0;
            int t = i;
            // here rj is the real multi-index corresponding to "i"
            for(int j=num_dimensions-1; j>=0; j--){
                p[j] = t % num_oned_points[j]; // index of the exponent is p, convert to spacial index rj (cumulative r)
                t /= num_oned_points[j];
                int rj = (p[j] % 2 == 0) ? (p[j]+1) / 2 : num_oned_points[j] - (p[j]+1) / 2; // +/- index
                r += rj * cnum_oned_points[j];
            }
            tensor_data[r].resize(1, std::complex<double> (basisFuncs[work->getSlot(p)]));
        }

        TasmanianFourierTransform::fast_fourier_transform(tensor_data, num_oned_points);

        auto v = tensor_data.begin(); // v++->data() means: get the 0-th entry from current v and move to the next v
        double tensorw = ((double) active_w[n]) / ((double) num_tensor_points);
        for(int i=0; i<num_tensor_points; i++){
            // We interpret this "i" as running through the spatial indexing; convert to internal
            int t=i;
            for(int j=num_dimensions-1; j>=0; j--){ // here p is the index of the spacial point in Tasmanian indexing
                p[j] = index_map[levels[j]][t % num_oned_points[j]];
                t /= num_oned_points[j];
            }
            weights[work->getSlot(p)] += tensorw * v++->data()->real();
        }
    }
}

void GridFourier::getQuadratureWeights(double weights[]) const{
    // When integrating the Fourier series on a tensored grid, all the
    // nonzero modes vanish, and we're left with the normalized Fourier
    // coeff for e^0 (sum of the data divided by number of points)

    IndexSet *work = (points != 0) ? points : needed;
    int num_points = work->getNumIndexes();
    IndexManipulator IM(num_dimensions);
    std::fill(weights, weights+num_points, 0.0);
    for(int n=0; n<active_tensors->getNumIndexes(); n++){
        const int *levels = active_tensors->getIndex(n);
        int num_tensor_points = 1;
        for(int j=0; j<num_dimensions; j++){
            num_tensor_points *= wrapper->getNumPoints(levels[j]);
        }
        std::vector<int> refs;
        IM.referencePoints<true>(levels, wrapper, work, refs);
        double tensorw = ((double) active_w[n]) / ((double) num_tensor_points);
        for(int i=0; i<num_tensor_points; i++){
            weights[refs[i]] += tensorw;
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
        const double *fcreal = fourier_coefs.getCStrip(i);
        const double *fcimag = fourier_coefs.getCStrip(i + num_points);
        double wr = wreal[i];
        double wi = wimag[i];
        for(int k=0; k<num_outputs; k++) y[k] += wr * fcreal[k] - wi * fcimag[k];
    }
}
#ifdef Tasmanian_ENABLE_BLAS
void GridFourier::evaluateFastCPUblas(const double x[], double y[]) const{
    int num_points = points->getNumIndexes();
    TasBLAS::setzero(num_outputs, y);
    std::vector<double> wreal(num_points);
    std::vector<double> wimag(num_points);
    computeBasis<double, false>(points, x, wreal.data(), wimag.data());
    TasBLAS::dgemv(num_outputs, num_points, fourier_coefs.getCStrip(0), wreal.data(), y, 1.0, 0.0);
    TasBLAS::dgemv(num_outputs, num_points, fourier_coefs.getCStrip(num_points), wimag.data(), y, -1.0, 1.0);
}
#else
void GridFourier::evaluateFastCPUblas(const double[], double[]) const{}
#endif

#ifdef Tasmanian_ENABLE_CUDA
void GridFourier::evaluateFastGPUcublas(const double x[], double y[]) const{
    prepareCudaData();
    int num_points = points->getNumIndexes();
    std::vector<double> wreal(num_points);
    std::vector<double> wimag(num_points);
    cudaDoubles gpuY;
    computeBasis<double, false>(points, x, wreal.data(), wimag.data());
    cuda_engine.cublasDGEMM(num_outputs, 1, num_points,  1.0, cuda_real, wreal, 0.0, gpuY);
    cuda_engine.cublasDGEMM(num_outputs, 1, num_points, -1.0, cuda_imag, wimag, 1.0, gpuY);
    gpuY.unload(y);
}
#else
void GridFourier::evaluateFastGPUcublas(const double[], double[]) const{}
#endif

void GridFourier::evaluateFastGPUcuda(const double x[], double y[]) const{
    evaluateFastGPUcublas(x,y);
}

#ifdef Tasmanian_ENABLE_MAGMA
void GridFourier::evaluateFastGPUmagma(int gpuID, const double x[], double y[]) const{
    prepareCudaData();
    int num_points = points->getNumIndexes();
    std::vector<double> wreal(num_points);
    std::vector<double> wimag(num_points);
    cudaDoubles gpuY;
    computeBasis<double, false>(points, x, wreal.data(), wimag.data());
    cuda_engine.magmaCudaDGEMM(gpuID, num_outputs, 1, num_points,  1.0, cuda_real, wreal, 0.0, gpuY);
    cuda_engine.magmaCudaDGEMM(gpuID, num_outputs, 1, num_points, -1.0, cuda_imag, wimag, 1.0, gpuY);
    gpuY.unload(y);
}
#else
void GridFourier::evaluateFastGPUmagma(int, const double[], double[]) const{}
#endif

void GridFourier::evaluateBatch(const double x[], int num_x, double y[]) const{
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        evaluate(&(x[((size_t) i) * ((size_t) num_dimensions)]), &(y[((size_t) i) * ((size_t) num_outputs)]));
    }
}

#ifdef Tasmanian_ENABLE_BLAS
void GridFourier::evaluateBatchCPUblas(const double x[], int num_x, double y[]) const{
    int num_points = points->getNumIndexes();
    Data2D<double> wreal;
    Data2D<double> wimag;
    evaluateHierarchicalFunctionsInternal(x, num_x, wreal, wimag);
    TasBLAS::dgemm(num_outputs, num_x, num_points, 1.0, fourier_coefs.getCStrip(0), wreal.getStrip(0), 0.0, y);
    TasBLAS::dgemm(num_outputs, num_x, num_points, -1.0, fourier_coefs.getCStrip(num_points), wimag.getStrip(0), 1.0, y);
}
#else
void GridFourier::evaluateBatchCPUblas(const double[], int, double[]) const{}
#endif

#ifdef Tasmanian_ENABLE_CUDA
void GridFourier::evaluateBatchGPUcublas(const double x[], int num_x, double y[]) const{
    int num_points = points->getNumIndexes();
    prepareCudaData();
    Data2D<double> wreal;
    Data2D<double> wimag;
    evaluateHierarchicalFunctionsInternal(x, num_x, wreal, wimag);
    cudaDoubles gpuY;
    cuda_engine.cublasDGEMM(num_outputs, num_x, num_points,  1.0, cuda_real, *(wreal.getVector()), 0.0, gpuY);
    cuda_engine.cublasDGEMM(num_outputs, num_x, num_points, -1.0, cuda_imag, *(wimag.getVector()), 1.0, gpuY);
    gpuY.unload(y);
}
void GridFourier::evaluateBatchGPUcuda(const double x[], int num_x, double y[]) const{
    prepareCudaData();
    int num_points = points->getNumIndexes();
    cudaDoubles wreal, wimag, gpu_x(num_dimensions, num_x, x);
    evaluateHierarchicalFunctionsInternalGPU(gpu_x.data(), num_x, wreal, wimag);
    cudaDoubles gpuY;
    cuda_engine.cublasDGEMM(num_outputs, num_x, num_points,  1.0, cuda_real, wreal, 0.0, gpuY);
    cuda_engine.cublasDGEMM(num_outputs, num_x, num_points, -1.0, cuda_imag, wimag, 1.0, gpuY);
    gpuY.unload(y);
}
#else
void GridFourier::evaluateBatchGPUcublas(const double[], int, double[]) const{}
void GridFourier::evaluateBatchGPUcuda(const double[], int, double[]) const{}
#endif



#ifdef Tasmanian_ENABLE_MAGMA
void GridFourier::evaluateBatchGPUmagma(int gpuID, const double x[], int num_x, double y[]) const{
    prepareCudaData();
    int num_points = points->getNumIndexes();
    cudaDoubles wreal, wimag, gpu_x(num_dimensions, num_x, x);
    evaluateHierarchicalFunctionsInternalGPU(gpu_x.data(), num_x, wreal, wimag);
    cudaDoubles gpuY;
    cuda_engine.magmaCudaDGEMM(gpuID, num_outputs, num_x, num_points,  1.0, cuda_real, wreal, 0.0, gpuY);
    cuda_engine.magmaCudaDGEMM(gpuID, num_outputs, num_x, num_points, -1.0, cuda_imag, wimag, 1.0, gpuY);
    gpuY.unload(y);
}
#else
void GridFourier::evaluateBatchGPUmagma(int, const double[], int, double[]) const{}
#endif

void GridFourier::integrate(double q[], double *conformal_correction) const{
    std::fill(q, q+num_outputs, 0.0);
    if (conformal_correction == 0){
        // everything vanishes except the Fourier coeff of e^0
        std::copy(fourier_coefs.getCStrip(0), fourier_coefs.getCStrip(0) + num_outputs, q);
    }else{
        // Do the expensive computation if we have a conformal map
        std::vector<double> w(getNumPoints());
        getQuadratureWeights(w.data());
        for(int i=0; i<points->getNumIndexes(); i++){
            w[i] *= conformal_correction[i];
            const double *v = values.getValues(i);
            for(int k=0; k<num_outputs; k++){
                q[k] += w[i] * v[k];
            }
        }
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
void GridFourier::evaluateHierarchicalFunctionsInternal(const double x[], int num_x, Data2D<double> &wreal, Data2D<double> &wimag) const{
    // when performing internal evaluations, split the matrix into real and complex components
    // thus only two real gemm() operations can be used (as opposed to one complex gemm)
    IndexSet *work = (points != 0) ? points : needed;
    int num_points = getNumPoints();
    Data2D<double> xx; xx.cload(num_dimensions, num_x, x);
    wreal.resize(num_points, num_x);
    wimag.resize(num_points, num_x);
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        computeBasis<double, false>(work, xx.getCStrip(i), wreal.getStrip(i), wimag.getStrip(i));
    }
}

void GridFourier::setHierarchicalCoefficients(const double c[], TypeAcceleration){
    // takes c to be length 2*num_outputs*num_points
    // first num_points*num_outputs are the real part; second num_points*num_outputs are the imaginary part

    clearAccelerationData();
    if (points == 0){
        points = needed;
        needed = 0;
    }
    fourier_coefs.resize(num_outputs, 2 * getNumPoints());
    std::copy(c, c + 2 * ((size_t) num_outputs) * ((size_t) getNumPoints()), fourier_coefs.getStrip(0));
}

#ifdef Tasmanian_ENABLE_CUDA
void GridFourier::evaluateHierarchicalFunctionsGPU(const double gpu_x[], int num_x, double gpu_y[]) const{
    loadCudaPoints();
    TasCUDA::devalfor(num_dimensions, num_x, max_levels, gpu_x, cuda_num_nodes, cuda_points, gpu_y, 0);
}
void GridFourier::evaluateHierarchicalFunctionsInternalGPU(const double gpu_x[], int num_x, cudaDoubles &wreal, cudaDoubles &wimag) const{
    size_t num_weights = ((size_t) points->getNumIndexes()) * ((size_t) num_x);
    if (wreal.size() != num_weights) wreal.resize(num_weights);
    if (wimag.size() != num_weights) wimag.resize(num_weights);
    loadCudaPoints();
    TasCUDA::devalfor(num_dimensions, num_x, max_levels, gpu_x, cuda_num_nodes, cuda_points, wreal.data(), wimag.data());
}
#endif

void GridFourier::clearAccelerationData(){
    #ifdef Tasmanian_ENABLE_CUDA
    cuda_real.clear();
    cuda_imag.clear();
    cuda_num_nodes.clear();
    cuda_points.clear();
    cuda_engine.reset();
    #endif
}
void GridFourier::clearRefinement(){ return; }     // to be expanded later
void GridFourier::mergeRefinement(){ return; }     // to be expanded later

const int* GridFourier::getPointIndexes() const{
    return ((points == 0) ? needed->getIndex(0) : points->getIndex(0));
}
const double* GridFourier::getFourierCoefs() const{
    return fourier_coefs.getCStrip(0);
}

} // end TasGrid

#endif
