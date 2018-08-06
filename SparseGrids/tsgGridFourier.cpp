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
                for(int i=0; i < num_outputs*getNumPoints(); i++){
                    ofs << " " << fourier_coefs[i].real() << " " << fourier_coefs[i].imag();
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

void GridFourier::read(std::ifstream &ifs, std::ostream *logstream){
    reset();
    ifs >> num_dimensions >> num_outputs;
    if (num_dimensions > 0){
        int flag;

        tensors = new IndexSet(num_dimensions);  tensors->read(ifs);
        active_tensors = new IndexSet(num_dimensions);  active_tensors->read(ifs);
        active_w = new int[active_tensors->getNumIndexes()];  for(int i=0; i<active_tensors->getNumIndexes(); i++){ ifs >> active_w[i]; }
        ifs >> flag; if (flag == 1){ points = new IndexSet(num_dimensions); points->read(ifs); }
        ifs >> flag; if (flag == 1){ needed = new IndexSet(num_dimensions); needed->read(ifs); }
        max_levels = new int[num_dimensions];  for(int j=0; j<num_dimensions; j++){ ifs >> max_levels[j]; }

        IndexSet *work = (points != 0) ? points : needed;

        if (num_outputs > 0){
            values = new StorageSet(0, 0); values->read(ifs);
            ifs >> flag;
            if (flag == 1){
                fourier_coefs = new std::complex<double>[num_outputs * work->getNumIndexes()];
                double fourier_real; double fourier_imag;
                for(int i=0; i<num_outputs*work->getNumIndexes(); i++){
                    ifs >> fourier_real >> fourier_imag;
                    fourier_coefs[i] = std::complex<double>(fourier_real, fourier_imag);
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
        wrapper = new OneDimensionalWrapper(&meta, oned_max_level, rule_fourier, 0.0, 0.0, logstream);

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
        work = 0;
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
        ofs.write((char*) max_levels, num_dimensions * sizeof(int));

        if (num_outputs > 0){
            values->writeBinary(ofs);
            if (fourier_coefs != 0){
                flag = 'y'; ofs.write(&flag, sizeof(char));
                const double *fcoefs = getFourierCoefs();     // real double array of length 2*num_outputs*num_nodes
                ofs.write((char*) fcoefs, 2 * getNumPoints() * num_outputs * sizeof(double));
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

void GridFourier::readBinary(std::ifstream &ifs, std::ostream *logstream){
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

        max_levels = new int[num_dimensions];
        ifs.read((char*) max_levels, num_dimensions * sizeof(int));

        if (num_outputs > 0){
            values = new StorageSet(0, 0); values->readBinary(ifs);
            ifs.read((char*) &flag, sizeof(char));
            if (flag == 'y'){
                double *fcoefs = new double[2 * num_outputs * work->getNumIndexes()];
                fourier_coefs = new std::complex<double>[num_outputs * work->getNumIndexes()];
                ifs.read((char*) fcoefs, 2 * num_outputs * work->getNumIndexes() * sizeof(double));
                for(int i=0; i<num_outputs * work->getNumIndexes(); i++){
                    fourier_coefs[i] = std::complex<double>(fcoefs[2*i], fcoefs[2*i+1]);
                }
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
        wrapper = new OneDimensionalWrapper(&meta, oned_max_level, rule_fourier, 0.0, 0.0, logstream);

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
        work = 0;
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
    if (max_levels != 0){ delete[] max_levels; max_levels = 0; }
    if (points != 0){ delete points; points = 0; }
    if (needed != 0){ delete needed; needed = 0; }
    if (exponents != 0){ delete exponents; exponents = 0; }
    if (values != 0){ delete values; values = 0; }
    if (fourier_coefs != 0){ delete[] fourier_coefs; fourier_coefs = 0; }
    num_dimensions = 0;
    num_outputs = 0;
}

void GridFourier::makeGrid(int cnum_dimensions, int cnum_outputs, int depth, TypeDepth type, const int* anisotropic_weights, const int* level_limits){
    IndexManipulator IM(cnum_dimensions);
    IndexSet *tset = IM.selectTensors(depth, type, anisotropic_weights, rule_fourier);
    if (level_limits != 0){
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
    max_levels = new int[num_dimensions];
    int max_level; IM.getMaxLevels(tensors, max_levels, max_level);
    wrapper = new OneDimensionalWrapper(&meta, max_level, rule_fourier);

    int* tensors_w = IM.makeTensorWeights(tensors);
    active_tensors = IM.nonzeroSubset(tensors, tensors_w);

    int nz_weights = active_tensors->getNumIndexes();

    active_w = new int[nz_weights];
    tensor_refs = new int*[nz_weights];
    int count = 0;
    for(int i=0; i<tensors->getNumIndexes(); i++){ if (tensors_w[i] != 0) active_w[count++] = tensors_w[i]; }

    delete[] tensors_w;

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

}

int* GridFourier::referenceExponents(const int levels[], const IndexSet *list){

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
        refs[i] = list->getSlot(p);
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

int GridFourier::convertIndexes(const int i, const int levels[]) const {

    /*
     * This routine ensures that function values are loaded into the Fourier transform in order of
     * increasing SPATIAL x-values.
     *
     * Take the example of the level-2 1D grid. Tasmanian's internal indexing reports the points as
     *
     * [0, 1, 2, ..., 8]
     *
     * which corresponds spatially to
     *
     * [0, 1/3, 2/3, 1/9, 2/9, 4/9, 5/9, 7/9, 8/9].
     *
     * That is, the new points added on level-2 appear after the level-1 points, even though spatially
     * they are interwoven. We now order the grid spatially, in terms of increasing x-values:
     *
     * [0, 1/9, ..., 8/9]
     *
     * The input parameter i is the i-th entry in the spatially ordered list. The output is the corresponding
     * internal index used by Tasmanian. For example, with p[0] = 2, convertIndexes(1, p) would return 4 (the
     * index of 1/9 in the internally ordered list of points).
     */

    IndexSet *work = (points == 0 ? needed : points);

    int* cnum_oned_points = new int[num_dimensions];
    for(int j=0; j<num_dimensions; j++){
        cnum_oned_points[j] = wrapper->getNumPoints(levels[j]);
    }

    // This i is spatial indexing
    int t=i;
    int* p=new int[num_dimensions];
    for(int j=num_dimensions-1; j>=0; j--){
        p[j] = t % cnum_oned_points[j];
        t /= cnum_oned_points[j];
    }
    // p[] stores the spatial index as a tensor address

    // Now we move from spatial to internal indexing
    int *p_internal = new int[num_dimensions];
    for(int j=0; j<num_dimensions; j++){
        int tmp = p[j];
        if (tmp == 0){
            p_internal[j] = 0;
        }else{
            int go_back = 0;
            while(tmp % 3 == 0 && tmp != 0){
                go_back += 1;
                tmp /= 3;
            }

            int num_prev_level = wrapper->getNumPoints(levels[j]-go_back-1);
            p_internal[j] = num_prev_level + (tmp/3)*2 + tmp % 3 - 1;
        }
    }

    int result = work->getSlot(p_internal);
    delete[] cnum_oned_points;
    delete[] p_internal;
    delete[] p;
    return result;
}

void GridFourier::calculateFourierCoefficients(){
    int num_nodes = getNumPoints();

    if (fourier_coefs != 0){ delete[] fourier_coefs; fourier_coefs = 0; }
    fourier_coefs = new std::complex<double>[num_outputs * num_nodes];
    std::fill(fourier_coefs, fourier_coefs + num_outputs*num_nodes, std::complex<double>(0,0));
    for(int k=0; k<num_outputs; k++){
        for(int n=0; n<active_tensors->getNumIndexes(); n++){
            const int* levels = active_tensors->getIndex(n);
            int num_tensor_points = 1;
            int* num_oned_points = new int[num_dimensions];
            for(int j=0; j<num_dimensions; j++){
                num_oned_points[j] = wrapper->getNumPoints(levels[j]);
                num_tensor_points *= num_oned_points[j];
            }

            std::complex<double> *in = new std::complex<double>[num_tensor_points];
            std::complex<double> *out = new std::complex<double>[num_tensor_points];
            for(int i=0; i<num_tensor_points; i++){
                // We interpret this "i" as running through the spatial indexing; convert to internal
                int key = convertIndexes(i, levels);
                const double *v = values->getValues(key);
                in[i] = v[k];
            }

            // Execute FFT
            TasmanianFourierTransform::discrete_fourier_transform(num_dimensions, num_oned_points, in, out);

            for(int i=0; i<num_tensor_points; i++){
                // Combine with tensor weights
                fourier_coefs[num_outputs*(exponent_refs[n][i]) + k] += ((double) active_w[n]) * out[i] / ((double) num_tensor_points);
            }
            delete[] in;
            delete[] out;
            delete[] num_oned_points;
        }
    }
}

std::complex<double>* GridFourier::getBasisFunctions(const double x[]) const {
    std::complex<double> *weights = new std::complex<double>[exponents->getNumIndexes()];
    getBasisFunctions(x,weights);
    return weights;
}
void GridFourier::getBasisFunctions(const double x[], std::complex<double> weights[]) const {
    std::complex<double> unit_imag(0.0, 1.0);    // this is sqrt(-1)
    for(int i=0; i<exponents->getNumIndexes(); i++){
        weights[i] = exp(2 * M_PI * unit_imag * ((double) exponents->getIndex(i)[0]) * x[0]);
        for(int j=1; j<exponents->getNumDimensions(); j++){
            weights[i] *= exp(2 * M_PI * unit_imag * ((double) exponents->getIndex(i)[j]) * x[j]);
        }
    }
}
void GridFourier::getBasisFunctions(const double x[], double weights[]) const {
    // weights has length 2*num_nodes here
    std::complex<double> *tmp = getBasisFunctions(x);
    for(int i=0; i<exponents->getNumIndexes(); i++){
        weights[2*i] = tmp[i].real();
        weights[2*i+1] = tmp[i].imag();
    }
    delete[] tmp;
}

void GridFourier::getInterpolationWeights(const double x[], double weights[]) const {
    /*
    I[f](x) = c^T * \Phi(x) = (U*P*f)^T * \Phi(x)           (U represents normalized forward Fourier transform; P represents reordering of f_i before going into FT)
                            = f^T * (P^T * U^T * \Phi(x))   (P^T = P^(-1) since P is a permutation matrix)

    Note that U is the DFT operator (complex) and the transposes are ONLY REAL transposes, so U^T = U.
    */

    std::fill(weights, weights+getNumPoints(), 0.0);
    std::complex<double> *basisFuncs = getBasisFunctions(x);

    for(int n=0; n<active_tensors->getNumIndexes(); n++){
        const int *levels = active_tensors->getIndex(n);
        int num_tensor_points = 1;
        int *num_oned_points = new int[num_dimensions];

        for(int j=0; j<num_dimensions; j++){
            num_oned_points[j] = wrapper->getNumPoints(levels[j]);
            num_tensor_points *= num_oned_points[j];
        }

        std::complex<double> *in = new std::complex<double>[num_tensor_points];
        std::complex<double> *out = new std::complex<double>[num_tensor_points];

        for(int i=0; i<num_tensor_points; i++){
            in[i] = basisFuncs[exponent_refs[n][i]];
        }

        TasmanianFourierTransform::discrete_fourier_transform(num_dimensions, num_oned_points, in, out);

        for(int i=0; i<num_tensor_points; i++){
            int key = convertIndexes(i, levels);
            weights[key] += ((double) active_w[n]) * out[i].real()/((double) num_tensor_points);
        }
        delete[] in;
        delete[] out;
        delete[] num_oned_points;
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
    std::complex<double> *w = getBasisFunctions(x);
    TasBLAS::setzero(num_outputs, y);
    for(int k=0; k<num_outputs; k++){
        for(int i=0; i<points->getNumIndexes(); i++){
            y[k] += (w[i] * fourier_coefs[i*num_outputs+k]).real();
        }
    }
    delete[] w;
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
void GridFourier::evaluateFastGPUcublas(const double x[], double y[], std::ostream*) const{
    evaluate(x,y);
}
void GridFourier::evaluateFastGPUcuda(const double x[], double y[], std::ostream*) const{
    evaluate(x,y);
}
void GridFourier::evaluateFastGPUmagma(int, const double x[], double y[], std::ostream*) const{
    evaluate(x,y);
}

void GridFourier::evaluateBatchCPUblas(const double x[], int num_x, double y[]) const{
    evaluateBatch(x, num_x, y);
}
void GridFourier::evaluateBatchGPUcublas(const double x[], int num_x, double y[], std::ostream*) const {
    evaluateBatch(x, num_x, y);
}
void GridFourier::evaluateBatchGPUcuda(const double x[], int num_x, double y[], std::ostream*) const {
    evaluateBatch(x, num_x, y);
}
void GridFourier::evaluateBatchGPUmagma(int, const double x[], int num_x, double y[], std::ostream*) const {
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
            q[k] = fourier_coefs[num_outputs * idx + k].real();
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
    // y must be of size num_x * num_nodes * 2
    int num_points = getNumPoints();
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        getBasisFunctions(&(x[((size_t) i) * ((size_t) num_dimensions)]), &(y[((size_t) i) * ((size_t) 2*num_points)]));
    }
}
void GridFourier::setHierarchicalCoefficients(const double c[], TypeAcceleration, std::ostream*){
    // takes c to be length 2*num_outputs*num_nodes
    // first two entries are real and imag parts of the Fourier coef for the first basis function and first output dimension
    // second two entries are real/imag parts of the Fourier coef for the second basis function and first output dim
    // and so on
    if (accel != 0) accel->resetGPULoadedData();
    if (points == 0){
        points = needed;
        needed = 0;
    }
    if (fourier_coefs != 0) delete[] fourier_coefs;
    fourier_coefs = new std::complex<double>[num_outputs * getNumPoints()];
    for(int i=0; i<num_outputs*getNumPoints(); i++){
        fourier_coefs[i] = std::complex<double>(c[2*i], c[2*i+1]);
    }
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
