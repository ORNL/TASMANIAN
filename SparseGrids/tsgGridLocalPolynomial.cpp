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

#ifndef __TASMANIAN_SPARSE_GRID_LPOLY_CPP
#define __TASMANIAN_SPARSE_GRID_LPOLY_CPP

#include "tsgGridLocalPolynomial.hpp"

#include "tsgHiddenExternals.hpp"
#include "tsgCudaMacros.hpp" // remove later

namespace TasGrid{

GridLocalPolynomial::GridLocalPolynomial() : num_dimensions(0), num_outputs(0), order(1), top_level(0),
                         surpluses(0), points(0), needed(0), values(0), parents(0), num_roots(0), roots(0), pntr(0), indx(0), rule(0),
                         accel(0), backend_flavor(flavor_auto), force_sparse(false)  {}
GridLocalPolynomial::GridLocalPolynomial(const GridLocalPolynomial &pwpoly) : num_dimensions(0), num_outputs(0), order(1), top_level(0),
                         surpluses(0), points(0), needed(0), values(0), parents(0), num_roots(0), roots(0), pntr(0), indx(0), rule(0),
                         accel(0), backend_flavor(flavor_auto), force_sparse(false)
{
    copyGrid(&pwpoly);
}

GridLocalPolynomial::~GridLocalPolynomial(){
    reset();
}
void GridLocalPolynomial::reset(bool clear_rule){
    clearAccelerationData();
    num_dimensions = num_outputs = top_level = 0;
    if (surpluses != 0){ delete[] surpluses; surpluses = 0; }
    if (points != 0){ delete points; points = 0; }
    if (needed != 0){ delete needed; needed = 0; }
    if (values != 0){ delete values; values = 0; }
    if (parents != 0){ delete[] parents; parents = 0; }
    if (roots != 0){ delete[] roots;  roots = 0; }
    if (pntr != 0){ delete[] pntr;  pntr = 0; }
    if (indx != 0){ delete[] indx;  indx = 0; }
    if (clear_rule){ rule = 0; order = 1; }
    backend_flavor = flavor_auto;
    force_sparse = false;
}

void GridLocalPolynomial::write(std::ofstream &ofs) const{
    ofs << std::scientific; ofs.precision(17);
    ofs << num_dimensions << " " << num_outputs << " " << order << " " << top_level << endl;
    if (num_dimensions > 0){
        ofs << OneDimensionalMeta::getIORuleString(rule->getType()) << endl;
        if (points == 0){
            ofs << "0" << endl;
        }else{
            ofs << "1 ";
            points->write(ofs);
        }
        if (surpluses == 0){
            ofs << "0" << endl;
        }else{
            ofs << "1 ";
            for(size_t i=0; i<((size_t) points->getNumIndexes()) * ((size_t) num_outputs); i++){ ofs << " " << surpluses[i]; } ofs << endl;
        }
        if (needed == 0){
            ofs << "0" << endl;
        }else{
            ofs << "1 ";
            needed->write(ofs);
        }
        if (parents == 0){
            ofs << "0" << endl;
        }else{
            ofs << "1";
            for(int i=0; i<rule->getMaxNumParents()*(points->getNumIndexes()); i++){ ofs << " " << parents[i]; } ofs << endl;
        }
        int num_points = (points == 0) ? needed->getNumIndexes() : points->getNumIndexes();
        ofs << num_roots; for(int i=0; i<num_roots; i++){  ofs << " " << roots[i];  } ofs << endl;
        ofs << pntr[0]; for(int i=1; i<=num_points; i++){  ofs << " " << pntr[i];  }  ofs << endl;
        ofs << indx[0]; for(int i=1; i<pntr[num_points]; i++){  ofs << " " << indx[i];  }  ofs << endl;
        if (num_outputs > 0) values->write(ofs);
    }
}
void GridLocalPolynomial::writeBinary(std::ofstream &ofs) const{
    int dims[4];
    dims[0] = num_dimensions;
    dims[1] = num_outputs;
    dims[2] = order;
    dims[3] = top_level;
    ofs.write((char*) dims, 4 * sizeof(int));
    if (num_dimensions > 0){
        dims[0] = OneDimensionalMeta::getIORuleInt(rule->getType());
        ofs.write((char*) dims, sizeof(int));

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
        if (surpluses == 0){
            flag = 'n'; ofs.write(&flag, sizeof(char));
        }else{
            flag = 'y'; ofs.write(&flag, sizeof(char));
            ofs.write((char*) surpluses, ((size_t) num_outputs) * ((size_t) points->getNumIndexes()) * sizeof(double));
        }
        if (parents == 0){
            flag = 'n'; ofs.write(&flag, sizeof(char));
        }else{
            flag = 'y'; ofs.write(&flag, sizeof(char));
            ofs.write((char*) parents, rule->getMaxNumParents()*(points->getNumIndexes()) * sizeof(int));
        }

        int num_points = (points == 0) ? needed->getNumIndexes() : points->getNumIndexes();
        ofs.write((char*) &num_roots, sizeof(int));
        ofs.write((char*) roots, num_roots * sizeof(int));
        ofs.write((char*) pntr, (num_points+1) * sizeof(int));
        ofs.write((char*) indx, pntr[num_points] * sizeof(int));

        if (num_outputs > 0) values->writeBinary(ofs);
    }
}
void GridLocalPolynomial::read(std::ifstream &ifs){
    reset();
    ifs >> num_dimensions >> num_outputs >> order >> top_level;
    if (num_dimensions > 0){
        int flag;
        std::string T;
        ifs >> T;
        TypeOneDRule crule = OneDimensionalMeta::getIORuleString(T.c_str());
        if (crule == rule_localp){
            rule = &rpoly;
        }else if (crule == rule_semilocalp){
            rule = &rsemipoly;
        }else if (crule == rule_localp0){
            rule = &rpoly0;
        }
        if (order == 0){
            rule = &rpolyc;
        }
        rule->setMaxOrder(order);

        ifs >> flag; if (flag == 1){ points = new IndexSet(num_dimensions); points->read(ifs); }
        ifs >> flag;
        if (flag == 1){
            surpluses = new double[((size_t) points->getNumIndexes()) * ((size_t) num_outputs)];
            for(size_t i=0; i<((size_t) points->getNumIndexes()) * ((size_t) num_outputs); i++){ ifs >> surpluses[i]; }
        }
        ifs >> flag;
        if (flag == 1){
            needed = new IndexSet(num_dimensions);
            needed->read(ifs);
        }
        ifs >> flag;
        if (flag == 1){
            int num_parents = rule->getMaxNumParents() * points->getNumIndexes();
            parents = new int[num_parents];
            for(int i=0; i<num_parents; i++){ ifs >> parents[i]; }
        }

        int num_points = (points == 0) ? needed->getNumIndexes() : points->getNumIndexes();

        ifs >> num_roots;
        roots = new int[num_roots];
        for(int i=0; i<num_roots; i++) ifs >> roots[i];
        pntr  = new int[num_points + 1];
        for(int i=0; i<=num_points; i++) ifs >> pntr[i];
        if (pntr[num_points] > 0){
            indx  = new int[pntr[num_points]];
            for(int i=0; i<pntr[num_points]; i++) ifs >> indx[i];
        }else{
            indx  = new int[1];  ifs >> indx[0]; // there is a special case when the grid has only one point without any children
        }

        if (num_outputs > 0){ values = new StorageSet(0, 0); values->read(ifs); }
    }
}
void GridLocalPolynomial::readBinary(std::ifstream &ifs){
    reset();
    int dims[4];
    ifs.read((char*) dims, 4 * sizeof(int));
    num_dimensions = dims[0];
    num_outputs = dims[1];
    order = dims[2];
    top_level = dims[3];
    if (num_dimensions > 0){
        ifs.read((char*) dims, sizeof(int));
        TypeOneDRule crule = OneDimensionalMeta::getIORuleInt(dims[0]);
        if (crule == rule_localp){
            rule = &rpoly;
        }else if (crule == rule_semilocalp){
            rule = &rsemipoly;
        }else if (crule == rule_localp0){
            rule = &rpoly0;
        }
        if (order == 0){
            rule = &rpolyc;
        }
        rule->setMaxOrder(order);

        char flag;
        ifs.read((char*) &flag, sizeof(char)); if (flag == 'y'){ points = new IndexSet(num_dimensions); points->readBinary(ifs); }
        ifs.read((char*) &flag, sizeof(char)); if (flag == 'y'){ needed = new IndexSet(num_dimensions); needed->readBinary(ifs); }

        ifs.read((char*) &flag, sizeof(char));
        if (flag == 'y'){
            surpluses = new double[((size_t) num_outputs) * ((size_t) points->getNumIndexes())];
            ifs.read((char*) surpluses, ((size_t) num_outputs) * ((size_t) points->getNumIndexes()) * sizeof(double));
        }

        ifs.read((char*) &flag, sizeof(char));
        if (flag == 'y'){
            parents = new int[rule->getMaxNumParents() * points->getNumIndexes()];
            ifs.read((char*) parents, rule->getMaxNumParents() * points->getNumIndexes() * sizeof(int));
        }

        int num_points = (points == 0) ? needed->getNumIndexes() : points->getNumIndexes();
        ifs.read((char*) &num_roots, sizeof(int));
        roots = new int[num_roots];
        ifs.read((char*) roots, num_roots * sizeof(int));

        pntr  = new int[num_points + 1];
        ifs.read((char*) pntr, (num_points+1) * sizeof(int));

        if (pntr[num_points] > 0){
            indx  = new int[pntr[num_points]];
            ifs.read((char*) indx, pntr[num_points] * sizeof(int));
        }else{
            indx  = new int[1];
            ifs.read((char*) indx, sizeof(int));
        }

        if (num_outputs > 0){ values = new StorageSet(0, 0); values->readBinary(ifs); }
    }
}

void GridLocalPolynomial::makeGrid(int cnum_dimensions, int cnum_outputs, int depth, int corder, TypeOneDRule crule, const int *level_limits){
    reset();
    num_dimensions = cnum_dimensions;
    num_outputs = cnum_outputs;
    order = corder;

    TypeOneDRule effective_rule = (crule == rule_localp0) ? rule_localp0 : (((crule == rule_semilocalp) && ((order == -1) || (order > 1))) ? rule_semilocalp : rule_localp);
    if (effective_rule == rule_localp){
        rule = &rpoly;
    }else if (effective_rule == rule_semilocalp){
        rule = &rsemipoly;
    }else if (effective_rule == rule_localp0){
        rule = &rpoly0;
    }
    if (order == 0){
        rule = &rpolyc;
    }
    rule->setMaxOrder(order);

    IndexManipulator IM(num_dimensions);
    UnsortedIndexSet* deltas = IM.getToalDegreeDeltas(depth);

    // Limits come here
    if (level_limits != 0){
        UnsortedIndexSet *limited = IM.removeIndexesByLimit(deltas, level_limits);
        if (limited != 0){
            delete deltas;
            deltas = limited;
        }
    }

    needed = IM.generatePointsFromDeltas(deltas, rule);
    delete deltas;

    buildTree();

    if (num_outputs == 0){
        points = needed;
        needed = 0;
        parents = IM.computeDAGupLocal(points, rule);
    }else{
        values = new StorageSet(num_outputs, needed->getNumIndexes());
    }
}

void GridLocalPolynomial::copyGrid(const GridLocalPolynomial *pwpoly){
    reset();
    num_dimensions = pwpoly->num_dimensions;
    num_outputs = pwpoly->num_outputs;
    order = pwpoly->order;

    if (pwpoly->rule->getType() == rule_localp){
        rule = &rpoly;
    }else if (pwpoly->rule->getType() == rule_semilocalp){
        rule = &rsemipoly;
    }else if (pwpoly->rule->getType() == rule_localp0){
        rule = &rpoly0;
    }
    if (pwpoly->rule->getMaxOrder() == 0){
        rule = &rpolyc;
    }
    rule->setMaxOrder(pwpoly->rule->getMaxOrder());

    if (pwpoly->points != 0) points = new IndexSet(pwpoly->points);
    if (pwpoly->needed != 0) needed = new IndexSet(pwpoly->needed);

    buildTree();

    if (pwpoly->values != 0) values = new StorageSet(pwpoly->values);

    if ((points != 0) && (num_outputs > 0)){ // points are loaded
        surpluses = new double[((size_t) points->getNumIndexes()) * ((size_t) num_outputs)];
        std::copy(pwpoly->surpluses, pwpoly->surpluses + ((size_t) points->getNumIndexes()) * ((size_t) num_outputs), surpluses);
    }
}

int GridLocalPolynomial::getNumDimensions() const{ return num_dimensions; }
int GridLocalPolynomial::getNumOutputs() const{ return num_outputs; }
TypeOneDRule GridLocalPolynomial::getRule() const{ return rule->getType(); }
int GridLocalPolynomial::getOrder() const{ return order; }

int GridLocalPolynomial::getNumLoaded() const{ return (((points == 0) || (num_outputs == 0)) ? 0 : points->getNumIndexes()); }
int GridLocalPolynomial::getNumNeeded() const{ return ((needed == 0) ? 0 : needed->getNumIndexes()); }
int GridLocalPolynomial::getNumPoints() const{ return ((points == 0) ? getNumNeeded() : points->getNumIndexes()); }

double* GridLocalPolynomial::getLoadedPoints() const{
    if (points == 0) return 0;
    int num_points = points->getNumIndexes();
    if (num_points == 0) return 0;
    double *x = new double[num_dimensions * num_points];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = points->getIndex(i);
        for(int j=0; j<num_dimensions; j++){
            x[i*num_dimensions + j] = rule->getNode(p[j]);
        }
    }
    return x;
}
void GridLocalPolynomial::getLoadedPoints(double *x) const{
    int num_points = points->getNumIndexes();
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = points->getIndex(i);
        for(int j=0; j<num_dimensions; j++){
            x[i*num_dimensions + j] = rule->getNode(p[j]);
        }
    }
}
double* GridLocalPolynomial::getNeededPoints() const{
    if (needed == 0) return 0;
    int num_points = needed->getNumIndexes();
    if (num_points == 0) return 0;
    double *x = new double[num_dimensions * num_points];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = needed->getIndex(i);
        for(int j=0; j<num_dimensions; j++){
            x[i*num_dimensions + j] = rule->getNode(p[j]);
        }
    }
    return x;
}
void GridLocalPolynomial::getNeededPoints(double *x) const{
    int num_points = needed->getNumIndexes();
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = needed->getIndex(i);
        for(int j=0; j<num_dimensions; j++){
            x[i*num_dimensions + j] = rule->getNode(p[j]);
        }
    }
}
double* GridLocalPolynomial::getPoints() const{
    return ((points == 0) ? getNeededPoints() : getLoadedPoints());
}
void GridLocalPolynomial::getPoints(double *x) const{
    if (points == 0){ getNeededPoints(x); }else{ getLoadedPoints(x); }
}

void GridLocalPolynomial::evaluate(const double x[], double y[]) const{
    int *monkey_count = new int[top_level+1];
    int *monkey_tail = new int[top_level+1];

    bool isSupported;
    size_t offset;

    std::fill(y, y + num_outputs, 0.0);

    for(int r=0; r<num_roots; r++){
        double basis_value = evalBasisSupported(points->getIndex(roots[r]), x, isSupported);

        if (isSupported){
            offset = roots[r] * num_outputs;
            for(int k=0; k<num_outputs; k++) y[k] += basis_value * surpluses[offset + k];

            int current = 0;
            monkey_tail[0] = roots[r];
            monkey_count[0] = pntr[roots[r]];

            while(monkey_count[0] < pntr[monkey_tail[0]+1]){
                if (monkey_count[current] < pntr[monkey_tail[current]+1]){
                    offset = indx[monkey_count[current]];
                    basis_value = evalBasisSupported(points->getIndex(offset), x, isSupported);
                    if (isSupported){
                        offset *= num_outputs;
                        for(int k=0; k<num_outputs; k++) y[k] += basis_value * surpluses[offset + k];
                        offset /= num_outputs;

                        monkey_tail[++current] = offset;
                        monkey_count[current] = pntr[offset];
                    }else{
                        monkey_count[current]++;
                    }
                }else{
                    monkey_count[--current]++;
                }
            }
        }
    }

    delete[] monkey_count;
    delete[] monkey_tail;
}

void GridLocalPolynomial::evaluateFastCPUblas(const double x[], double y[]) const{ evaluate(x, y); } // standard BLAS cannot accelerate dense matrix times a sparse vector
#if defined(Tasmanian_ENABLE_CUBLAS) || defined(Tasmanian_ENABLE_CUDA)
void GridLocalPolynomial::evaluateFastGPUcublas(const double x[], double y[], std::ostream *os) const{
    if (num_outputs < 64){
        evaluate(x, y);
        return;
    }

    int num_points = points->getNumIndexes();
    makeCheckAccelerationData(accel_gpu_cublas, os);
    checkAccelerationGPUValues();
    AccelerationDataGPUFull *gpu_acc = (AccelerationDataGPUFull*) accel;

    int num_nz, *sindx = 0;
    double *svals = 0;
    buildSparseVector<false>(x, num_nz, 0, 0);
    sindx = new int[num_nz];
    svals = new double[num_nz];
    buildSparseVector<true>(x, num_nz, sindx, svals);
    int *gpu_sindx = TasCUDA::cudaSend<int>(num_nz, sindx, os);
    double *gpu_svals = TasCUDA::cudaSend<double>(num_nz, svals, os);
    double *gpu_y = TasCUDA::cudaNew<double>(num_outputs, os);

    #ifdef Tasmanian_ENABLE_CUBLAS
    gpu_acc->cusparseMatveci(num_outputs, num_points, num_nz, gpu_sindx, gpu_svals, gpu_y);
    #else
    TasCUDA::cudaSparseVecDenseMat(num_outputs, num_points, num_nz, gpu_acc->getGPUValues(), gpu_sindx, gpu_svals, gpu_y);
    #endif

    TasCUDA::cudaRecv<double>(num_outputs, gpu_y, y, os);

    TasCUDA::cudaDel<int>(gpu_sindx, os);
    TasCUDA::cudaDel<double>(gpu_svals, os);
    TasCUDA::cudaDel<double>(gpu_y, os);
    delete[] sindx;
    delete[] svals;
}
#else
void GridLocalPolynomial::evaluateFastGPUcublas(const double x[], double y[], std::ostream*) const{ evaluate(x, y); }
#endif
// evaluation of a single x cannot be accelerated with a gpu (not parallelizable), do that on the CPU and use the GPU only for the case of many outputs
void GridLocalPolynomial::evaluateFastGPUcuda(const double x[], double y[], std::ostream* os) const{ evaluateFastGPUcublas(x, y, os); }
void GridLocalPolynomial::evaluateFastGPUmagma(const double x[], double y[], std::ostream*) const{ evaluate(x, y); }

void GridLocalPolynomial::evaluateBatch(const double x[], int num_x, double y[]) const{
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        evaluate(&(x[((size_t) i) * ((size_t) num_dimensions)]), &(y[((size_t) i) * ((size_t) num_outputs)]));
    }
}
void GridLocalPolynomial::evaluateBatchCPUblas(const double x[], int num_x, double y[]) const{
    if (force_sparse || (((backend_flavor == flavor_auto) || (backend_flavor == flavor_sparse_sparse)) && (num_outputs <= TSG_LOCALP_BLAS_NUM_OUTPUTS))){
        evaluateBatch(x, num_x, y);
        return;
    }

    int *sindx, *spntr;
    double *svals;
    buildSpareBasisMatrix(x, num_x, 32, spntr, sindx, svals); // build sparse matrix corresponding to x

    int num_points = (points == 0) ? needed->getNumIndexes() : points->getNumIndexes();
    double nnz = (double) spntr[num_x];
    double total_size = ((double) num_x) * ((double) num_points);

    if ((backend_flavor == flavor_sparse_dense) || (backend_flavor == flavor_dense_dense) || ((backend_flavor == flavor_auto) && (nnz / total_size > 0.1))){
        // potentially wastes a lot of memory
        double *A = new double[((size_t) num_x) * ((size_t) num_points)];
        std::fill(A, A + ((size_t) num_x) * ((size_t) num_points), 0.0);
        for(int i=0; i<num_x; i++){
            for(int j=spntr[i]; j<spntr[i+1]; j++){
                A[((size_t) i) * ((size_t) num_points) + ((size_t) sindx[j])] = svals[j];
            }
        }
        TasBLAS::dgemm(num_outputs, num_x, num_points, 1.0, surpluses, A, 0.0, y);
        delete[] A;
    }else{
        #pragma omp parallel for
        for(int i=0; i<num_x; i++){
            double *this_y = &(y[i*num_outputs]);
            std::fill(this_y, this_y + num_outputs, 0.0);
            for(int j=spntr[i]; j<spntr[i+1]; j++){
                double v = svals[j];
                double *this_surp = &(surpluses[((size_t) sindx[j]) * ((size_t) num_outputs)]);
                for(int k=0; k<num_outputs; k++) this_y[k] += v * this_surp[k];
            }
        }
    }

    delete[] sindx;
    delete[] spntr;
    delete[] svals;
}
#ifdef Tasmanian_ENABLE_CUBLAS
void GridLocalPolynomial::evaluateBatchGPUcublas(const double x[], int num_x, double y[], std::ostream *os) const{
    int num_points = points->getNumIndexes();
    makeCheckAccelerationData(accel_gpu_cublas, os);
    checkAccelerationGPUValues();
    AccelerationDataGPUFull *gpu_acc = (AccelerationDataGPUFull*) accel;

    int *sindx, *spntr;
    double *svals;
    buildSpareBasisMatrix(x, num_x, 32, spntr, sindx, svals); // build sparse matrix corresponding to x

    gpu_acc->cusparseMatmul(true, num_points, num_outputs, num_x, spntr, sindx, svals, 0, y); // true for pointers residing on the cpu

    delete[] svals;
    delete[] sindx;
    delete[] spntr;
}
#else
void GridLocalPolynomial::evaluateBatchGPUcublas(const double x[], int num_x, double y[], std::ostream *) const{ evaluateBatchCPUblas(x, num_x, y); }
#endif // Tasmanian_ENABLE_CUDA

void GridLocalPolynomial::evaluateBatchGPUcuda(const double x[], int num_x, double y[], std::ostream *os) const{
    #ifdef Tasmanian_ENABLE_CUDA
    if ((order == -1) || (order > 2)){ // GPU evaluations are availabe only for order 0, 1, and 2. Cubic will come later, but higher order will not be supported
        evaluateBatchGPUcublas(x, num_x, y, os);
        return;
    }
    int num_points = points->getNumIndexes();
    makeCheckAccelerationData(accel_gpu_cuda, os);
    checkAccelerationGPUValues();
    checkAccelerationGPUNodes();
    AccelerationDataGPUFull *gpu_acc = (AccelerationDataGPUFull*) accel;

    TypeLocalPolynomialBackendFlavor flv = backend_flavor;
    if (backend_flavor == flavor_auto){
        flv = (num_points > 1024) ? flavor_sparse_sparse : flavor_dense_dense;
    }
    if (force_sparse) flv = flavor_sparse_sparse;
    //flv = flavor_sparse_dense;
    if (flv == flavor_dense_dense){
        double *gpu_x = TasCUDA::cudaSend<double>(num_x * num_dimensions, x, os);
        double *gpu_weights = TasCUDA::cudaNew<double>(((size_t) num_x) * ((size_t) points->getNumIndexes()), os);
        double *gpu_result = TasCUDA::cudaNew<double>(((size_t) num_x) * ((size_t) values->getNumOutputs()), os);

        buildDenseBasisMatrixGPU(gpu_x, num_x, gpu_weights, os);
        #ifdef Tasmanian_ENABLE_CUBLAS
        gpu_acc->cublasDGEMM(values->getNumOutputs(), num_x, num_points, gpu_weights, gpu_result);
        #else
        TasCUDA::cudaDgemm(values->getNumOutputs(), num_x, num_points, gpu_acc->getGPUValues(), gpu_weights, gpu_result);
        #endif // Tasmanian_ENABLE_CUBLAS

        TasCUDA::cudaRecv<double>(num_x * values->getNumOutputs(), gpu_result, y, os);

        TasCUDA::cudaDel<double>(gpu_result, os);
        TasCUDA::cudaDel<double>(gpu_weights, os);
        TasCUDA::cudaDel<double>(gpu_x, os);
    }else if (flv == flavor_sparse_sparse){
        double *gpu_x = TasCUDA::cudaSend<double>(num_x * num_dimensions, x, os);
        double *gpu_y = TasCUDA::cudaNew<double>(((size_t) num_x) * ((size_t) num_outputs), os);

        int *gpu_spntr, *gpu_sindx, num_nz = 0;
        double *gpu_svals;
        buildSparseBasisMatrixGPU(gpu_x, num_x, gpu_spntr, gpu_sindx, gpu_svals, num_nz, os);

        #ifdef Tasmanian_ENABLE_CUBLAS
        if (num_outputs == 1){
            gpu_acc->cusparseMatvec(num_points, num_x, gpu_spntr, gpu_sindx, gpu_svals, num_nz, gpu_y);
        }else{
            gpu_acc->cusparseMatmul(false, num_points, num_outputs, num_x, gpu_spntr, gpu_sindx, gpu_svals, num_nz, gpu_y); // false for pointers already living on the gpu
        }
        #else
        TasCUDA::cudaSparseMatmul(num_x, num_outputs, num_nz, gpu_spntr, gpu_sindx, gpu_svals, gpu_acc->getGPUValues(), gpu_y);
        #endif // Tasmanian_ENABLE_CUBLAS
        TasCUDA::cudaRecv<double>(((size_t) num_x) * ((size_t) num_outputs), gpu_y, y, os);

        TasCUDA::cudaDel<int>(gpu_spntr, os);
        TasCUDA::cudaDel<int>(gpu_sindx, os);
        TasCUDA::cudaDel<double>(gpu_svals, os);
        TasCUDA::cudaDel<double>(gpu_y, os);
        TasCUDA::cudaDel<double>(gpu_x, os);
    }else if (flv == flavor_sparse_dense){
        double *gpu_x = TasCUDA::cudaSend<double>(num_x * num_dimensions, x, os);
        double *gpu_weights = TasCUDA::cudaNew<double>(((size_t) num_x) * ((size_t) points->getNumIndexes()), os);
        double *gpu_result = TasCUDA::cudaNew<double>(((size_t) num_x) * ((size_t) values->getNumOutputs()), os);

        checkAccelerationGPUHierarchy();
        TasCUDA::devalpwpoly_sparse_dense(order, rule->getType(), num_dimensions, num_x, num_points, gpu_x, gpu_acc->getGPUNodes(), gpu_acc->getGPUSupport(),
                                gpu_acc->getGPUpntr(), gpu_acc->getGPUindx(), num_roots, gpu_acc->getGPUroots(), gpu_weights);

        #ifdef Tasmanian_ENABLE_CUBLAS
        gpu_acc->cublasDGEMM(values->getNumOutputs(), num_x, num_points, gpu_weights, gpu_result);
        #else
        TasCUDA::cudaDgemm(values->getNumOutputs(), num_x, num_points, gpu_acc->getGPUValues(), gpu_weights, gpu_result);
        #endif // Tasmanian_ENABLE_CUBLAS

        TasCUDA::cudaRecv<double>(((size_t) num_x) * ((size_t) values->getNumOutputs()), gpu_result, y, os);

        TasCUDA::cudaDel<double>(gpu_result, os);
        TasCUDA::cudaDel<double>(gpu_weights, os);
        TasCUDA::cudaDel<double>(gpu_x, os);
    }
    #else
    evaluateBatchGPUcublas(x, num_x, y, os);
    #endif // Tasmanian_ENABLE_CUDA
}
void GridLocalPolynomial::evaluateBatchGPUmagma(const double x[], int num_x, double y[], std::ostream *os) const{
    evaluateBatchGPUcublas(x, num_x, y, os);
}

void GridLocalPolynomial::loadNeededPoints(const double *vals, TypeAcceleration acc){
    if (points == 0){
        values->setValues(vals);
        points = needed;
        needed = 0;
    }else if (needed == 0){
        values->setValues(vals);
    }else{
        values->addValues(points, needed, vals);
        points->addIndexSet(needed);
        delete needed; needed = 0;
        buildTree();
    }
    if (accel != 0) accel->resetGPULoadedData();
    if (acc == accel_gpu_cublas){
        recomputeSurplusesGPUcublas();
    }else if (acc == accel_gpu_cuda){
        recomputeSurplusesGPUcuda();
    }else{
        recomputeSurpluses();
    }
}
void GridLocalPolynomial::mergeRefinement(){
    if (needed == 0) return; // nothing to do
    int num_all_points = getNumLoaded() + getNumNeeded();
    size_t num_vals = ((size_t) num_all_points) * ((size_t) num_outputs);
    double *vals = new double[num_vals];
    std::fill(vals, vals + num_vals, 0.0);
    values->setValuesPointer(vals, num_all_points);
    if (points == 0){
        points = needed;
        needed = 0;
    }else{
        points->addIndexSet(needed);
        delete needed; needed = 0;
        buildTree();
    }
    if (surpluses != 0) delete[] surpluses;
    surpluses = new double[num_vals];
    std::fill(surpluses, surpluses + num_vals, 0.0);
}


double* GridLocalPolynomial::getInterpolationWeights(const double x[]) const{
    IndexSet *work = (points == 0) ? needed : points;
    double *weights = new double[work->getNumIndexes()];
    getInterpolationWeights(x, weights);
    return weights;
}

void GridLocalPolynomial::getInterpolationWeights(const double x[], double *weights) const{
    IndexSet *work = (points == 0) ? needed : points;

    int *active_points = new int[work->getNumIndexes()];
    std::fill(weights, weights + work->getNumIndexes(), 0.0);
    int num_active_points = 0;

    int *monkey_count = new int[top_level+1];
    int *monkey_tail = new int[top_level+1];

    double basis_value;
    bool isSupported;
    int offset;

    for(int r=0; r<num_roots; r++){

        basis_value = evalBasisSupported(work->getIndex(roots[r]), x, isSupported);

        if (isSupported){
            active_points[num_active_points++] = roots[r];
            weights[roots[r]] = basis_value;

            int current = 0;
            monkey_tail[0] = roots[r];
            monkey_count[0] = pntr[roots[r]];

            while(monkey_count[0] < pntr[monkey_tail[0]+1]){
                if (monkey_count[current] < pntr[monkey_tail[current]+1]){
                    offset = indx[monkey_count[current]];

                    basis_value = evalBasisSupported(work->getIndex(offset), x, isSupported);

                    if (isSupported){
                        active_points[num_active_points++] = offset;
                        weights[offset] = basis_value;

                        monkey_tail[++current] = offset;
                        monkey_count[current] = pntr[offset];
                    }else{
                        monkey_count[current]++;
                    }
                }else{
                    monkey_count[--current]++;
                }
            }
        }
    }

    // apply the transpose of the surplus transformation
    const int *dagUp;
    if (num_outputs == 0){
        dagUp = parents;
    }else{
        IndexManipulator IM(num_dimensions);
        dagUp = IM.computeDAGupLocal(work, rule);
    }

    int *level = new int[num_active_points];
    int active_top_level = 0;
    for(int i=0; i<num_active_points; i++){
        const int *p = work->getIndex(active_points[i]);
        level[i] = rule->getLevel(p[0]);
        for(int j=1; j<num_dimensions; j++){
            level[i] += rule->getLevel(p[j]);
        }
        if (active_top_level < level[i]) active_top_level = level[i];
    }

    bool *used = new bool[work->getNumIndexes()];
    double *node = new double[num_dimensions];
    int max_parents = rule->getMaxNumParents() * num_dimensions;

    for(int l=active_top_level; l>0; l--){
        for(int i=0; i<num_active_points; i++){
            if (level[i] == l){
                const int* p = work->getIndex(active_points[i]);
                for(int j=0; j<num_dimensions; j++) node[j] = rule->getNode(p[j]);

                std::fill(used, used + work->getNumIndexes(), false);

                monkey_count[0] = 0;
                monkey_tail[0] = active_points[i];
                int current = 0;

                while(monkey_count[0] < max_parents){
                    if (monkey_count[current] < max_parents){
                        int branch = dagUp[ monkey_tail[current] * max_parents + monkey_count[current] ];
                        if ((branch == -1) || used[branch]){
                            monkey_count[current]++;
                        }else{
                            const int *func = work->getIndex(branch);
                            basis_value = rule->evalRaw(func[0], node[0]);
                            for(int j=1; j<num_dimensions; j++) basis_value *= rule->evalRaw(func[j], node[j]);
                            weights[branch] -= weights[active_points[i]] * basis_value;
                            used[branch] = true;

                            monkey_count[++current] = 0;
                            monkey_tail[current] = branch;
                        }
                    }else{
                        monkey_count[--current]++;
                    }
                }
            }
        }
    }

    delete[] used;
    delete[] node;

    delete[] level;

    delete[] monkey_count;
    delete[] monkey_tail;

    delete[] active_points;

    if (num_outputs > 0){
        delete[] dagUp;
    }
}

void GridLocalPolynomial::evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const{
    IndexSet *work = (points == 0) ? needed : points;
    int num_points = work->getNumIndexes();
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        const double *this_x = &(x[i*num_dimensions]);
        double *this_y = &(y[((size_t) i) * ((size_t) num_points)]);
        bool dummy;
        for(int j=0; j<num_points; j++){
            this_y[j] = evalBasisSupported(work->getIndex(j), this_x, dummy);
        }
    }
}


void GridLocalPolynomial::recomputeSurplusesGPUcublas(){
#ifdef Tasmanian_ENABLE_CUBLAS
    int num_points = points->getNumIndexes();
    if (surpluses != 0) delete[] surpluses;
    surpluses = new double[num_points * num_outputs];

    if ((accel != 0) && (!accel->isCompatible(accel_gpu_cublas))){
        delete accel;
        accel = 0;
    }
    if (accel == 0){ accel = (BaseAccelerationData*) (new AccelerationDataGPUFull()); }
    AccelerationDataGPUFull *gpu = (AccelerationDataGPUFull*) accel;
    gpu->setLogStream(&cerr); // debug mode

    double *x = getPoints();

    int *sindx, *spntr;
    double *svals;
    buildSpareBasisMatrix(x, num_points, 32, spntr, sindx, svals); // build sparse matrix corresponding to x
    delete[] x;

    int *rindx, *rpntr;
    double *rvals;
    int c = 0;
    for(int i=0; i<spntr[num_points]; i++) if (svals[i] != 0.0) c++;
    rindx = new int[c];
    rvals = new double[c];
    rpntr = new int[num_points + 1];
    rpntr[0] = 0;
    c = 0;
    for(int i=0; i<num_points; i++){
        for(int j=spntr[i]; j<spntr[i+1]; j++){
            if ((svals[j] != 0.0) && (sindx[j] != i)){
                rindx[c] = sindx[j];
                rvals[c] = svals[j];
                c++;
            }
        }
        rindx[c] = i;
        rvals[c] = 1.0;
        c++;
        rpntr[i+1] = c;
    }

    gpu->cusparseDCRSMM(num_points, num_outputs, rpntr, rindx, rvals, values->getValues(0), surpluses);

    delete[] sindx;
    delete[] svals;
    delete[] spntr;
    delete[] rindx;
    delete[] rvals;
    delete[] rpntr;
#else
    recomputeSurpluses();
#endif // Tasmanian_ENABLE_CUBLAS
}
void GridLocalPolynomial::recomputeSurplusesGPUcuda(){
#ifdef Tasmanian_ENABLE_CUDA
    int num_points = points->getNumIndexes();
    if (surpluses != 0) delete[] surpluses;
    surpluses = new double[((size_t) num_points) * ((size_t) num_outputs)];

    int *level = new int[num_points];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = points->getIndex(i);
        level[i] = rule->getLevel(p[0]);
        for(int j=1; j<num_dimensions; j++){
            level[i] += rule->getLevel(p[j]);
        }
    }

    double *x = getPoints();

    int *sindx, *spntr;
    double *svals;
    buildSpareBasisMatrix(x, num_points, 32, spntr, sindx, svals); // build sparse matrix corresponding to x
    delete[] x;

    int *rindx, *rpntr;
    double *rvals;
    int c = 0;
    for(int i=0; i<spntr[num_points]; i++) if (svals[i] != 0.0) c++;
    rindx = new int[c];
    rvals = new double[c];
    rpntr = new int[num_points + 1];
    rpntr[0] = 0;
    c = 0;
    for(int i=0; i<num_points; i++){
        for(int j=spntr[i]; j<spntr[i+1]; j++){
            if ((svals[j] != 0.0) && (sindx[j] != i)){
                rindx[c] = sindx[j];
                rvals[c] = svals[j];
                c++;
            }
        }
        rpntr[i+1] = c;
    }

    TasCUDA::d3gecss(num_outputs, num_points, level, top_level, rpntr, rindx, rvals, values->getValues(0), surpluses, &cerr);

    delete[] sindx;
    delete[] svals;
    delete[] spntr;
    delete[] rindx;
    delete[] rvals;
    delete[] rpntr;
    delete[] level;
#else
    recomputeSurpluses();
#endif // Tasmanian_ENABLE_CUDA
}

void GridLocalPolynomial::recomputeSurpluses(){
    int num_points = points->getNumIndexes();
    if (surpluses != 0) delete[] surpluses;
    surpluses = new double[((size_t) num_points) * ((size_t) num_outputs)];

    const double* v = values->getValues(0);
    std::copy(v, v + ((size_t) num_points) * ((size_t) num_outputs), surpluses);

    IndexManipulator IM(num_dimensions);
    int *dagUp = IM.computeDAGupLocal(points, rule);

    //int max_parents = (rule->isSemiLocal()) ? 2*num_dimensions : num_dimensions;
    int max_parents = rule->getMaxNumParents() * num_dimensions;

    int *level = new int[num_points];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = points->getIndex(i);
        level[i] = rule->getLevel(p[0]);
        for(int j=1; j<num_dimensions; j++){
            level[i] += rule->getLevel(p[j]);
        }
    }

    for(int l=1; l<=top_level; l++){
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<num_points; i++){
            if (level[i] == l){
                const int* p = points->getIndex(i);
                double *x = new double[num_dimensions];
                double *surpi = &(surpluses[((size_t) i) * ((size_t) num_outputs)]);
                for(int j=0; j<num_dimensions; j++) x[j] = rule->getNode(p[j]);

                int *monkey_count = new int[top_level + 1];
                int *monkey_tail = new int[top_level + 1];
                bool *used = new bool[num_points];
                std::fill(used, used + num_points, false);

                int current = 0;

                monkey_count[0] = 0;
                monkey_tail[0] = i;

                while(monkey_count[0] < max_parents){
                    if (monkey_count[current] < max_parents){
                        int branch = dagUp[monkey_tail[current] * max_parents + monkey_count[current]];
                        if ((branch == -1) || (used[branch])){
                            monkey_count[current]++;
                        }else{
                            const double *branch_surp = &(surpluses[((size_t) branch) * ((size_t) num_outputs)]);
                            double basis_value = evalBasisRaw(points->getIndex(branch), x);
                            for(int k=0; k<num_outputs; k++){
                                surpi[k] -= basis_value * branch_surp[k];
                            }
                            used[branch] = true;

                            monkey_count[++current] = 0;
                            monkey_tail[current] = branch;
                        }
                    }else{
                        monkey_count[--current]++;
                    }
                }

                delete[] x;
                delete[] used;
                delete[] monkey_tail;
                delete[] monkey_count;
            }
        }
    }

    delete[] dagUp;
    delete[] level;
}

double GridLocalPolynomial::evalBasisRaw(const int point[], const double x[]) const{
    double f = rule->evalRaw(point[0], x[0]);
    for(int j=1; j<num_dimensions; j++) f *= rule->evalRaw(point[j], x[j]);
    return f;
}
double GridLocalPolynomial::evalBasisSupported(const int point[], const double x[], bool &isSupported) const{
    double f = rule->evalSupport(point[0], x[0], isSupported);
    if (!isSupported) return 0.0;
    for(int j=1; j<num_dimensions; j++){
        f *= rule->evalSupport(point[j], x[j], isSupported);
        if (!isSupported) return 0.0;
    }
    return f;
}

void GridLocalPolynomial::buildSpareBasisMatrix(const double x[], int num_x, int num_chunk, int* &spntr, int* &sindx, double* &svals) const{
    int num_blocks, num_last, stripe_size;

    int *stripes, *last_stripe_size, **tpntr, ***tindx;
    double ***tvals;

    buildSparseMatrixBlockForm(x, num_x, num_chunk, num_blocks, num_last, stripe_size, stripes, last_stripe_size, tpntr, tindx, tvals);

    // assemble the matrix from the local junks
    spntr = new int[num_x+1];
    spntr[0] = 0;
    int c = 0, block_offset = 0;
    for(int b=0; b<num_blocks; b++){
        int this_size = (b == num_blocks-1) ? num_last : num_chunk;
        for(int i=0; i<this_size; i++){
            spntr[c+1] = block_offset + tpntr[b][i+1];
            c++;
        }
        block_offset = spntr[c];
        delete[] tpntr[b];
    }
    delete[] tpntr;

    int num_nz = spntr[num_x];
    sindx = new int[num_nz];
    svals = new double[num_nz];
    c = 0;
    for(int b=0; b<num_blocks; b++){
        for(int s=0; s<stripes[b]-1; s++){
            std::copy(tindx[b][s], tindx[b][s] + stripe_size, &(sindx[c]));
            std::copy(tvals[b][s], tvals[b][s] + stripe_size, &(svals[c]));
            c += stripe_size;
            delete[] tindx[b][s];
            delete[] tvals[b][s];
        }
        std::copy(tindx[b][stripes[b]-1], tindx[b][stripes[b]-1] + last_stripe_size[b], &(sindx[c]));
        std::copy(tvals[b][stripes[b]-1], tvals[b][stripes[b]-1] + last_stripe_size[b], &(svals[c]));
        c += last_stripe_size[b];
        delete[] tindx[b][stripes[b]-1];
        delete[] tvals[b][stripes[b]-1];
        delete[] tindx[b];
        delete[] tvals[b];
    }
    delete[] tindx;
    delete[] tvals;
    delete[] stripes;
    delete[] last_stripe_size;
}
void GridLocalPolynomial::buildSpareBasisMatrixStatic(const double x[], int num_x, int num_chunk, int *spntr, int *sindx, double *svals) const{
    int num_blocks, num_last, stripe_size;

    int *stripes, *last_stripe_size, **tpntr, ***tindx;
    double ***tvals;

    buildSparseMatrixBlockForm(x, num_x, num_chunk, num_blocks, num_last, stripe_size, stripes, last_stripe_size, tpntr, tindx, tvals);

    // assemble the matrix from the local junks
    spntr[0] = 0;
    int c = 0, block_offset = 0;
    for(int b=0; b<num_blocks; b++){
        int this_size = (b == num_blocks-1) ? num_last : num_chunk;
        for(int i=0; i<this_size; i++){
            spntr[c+1] = block_offset + tpntr[b][i+1];
            c++;
        }
        block_offset = spntr[c];
        delete[] tpntr[b];
    }
    delete[] tpntr;

    c = 0;
    for(int b=0; b<num_blocks; b++){
        for(int s=0; s<stripes[b]-1; s++){
            std::copy(tindx[b][s], tindx[b][s] + stripe_size, &(sindx[c]));
            std::copy(tvals[b][s], tvals[b][s] + stripe_size, &(svals[c]));
            c += stripe_size;
            delete[] tindx[b][s];
            delete[] tvals[b][s];
        }
        std::copy(tindx[b][stripes[b]-1], tindx[b][stripes[b]-1] + last_stripe_size[b], &(sindx[c]));
        std::copy(tvals[b][stripes[b]-1], tvals[b][stripes[b]-1] + last_stripe_size[b], &(svals[c]));
        c += last_stripe_size[b];
        delete[] tindx[b][stripes[b]-1];
        delete[] tvals[b][stripes[b]-1];
        delete[] tindx[b];
        delete[] tvals[b];
    }
    delete[] tindx;
    delete[] tvals;
    delete[] stripes;
    delete[] last_stripe_size;
}
int GridLocalPolynomial::getSpareBasisMatrixNZ(const double x[], int num_x, int num_chunk) const{
    IndexSet *work = (points == 0) ? needed : points;
    int num_blocks = (num_x / num_chunk) + ((num_x % num_chunk > 0) ? 1 : 0); // number of blocks
    int num_last = ((num_x % num_chunk > 0) ? num_x % num_chunk : num_chunk);

    int total_nz = 0;

    for(int b=0; b<num_blocks; b++){
        int c = 0; // count the current entry

        int *monkey_count = new int[top_level+1]; // monkey business (traversing threes)
        int *monkey_tail = new int[top_level+1];
        bool isSupported;
        int offset;

        int this_size = (b == num_blocks-1) ? num_last : num_chunk;
        for(int i=0; i<this_size; i++){ // for each point
            const double *this_x = &(x[(b*num_chunk + i) * num_dimensions]);

            for(int r=0; r<num_roots; r++){
                evalBasisSupported(work->getIndex(roots[r]), this_x, isSupported);

                if (isSupported){
                    // add point to this sparse column
                    c++;

                    int current = 0;
                    monkey_tail[0] = roots[r];
                    monkey_count[0] = pntr[roots[r]];

                    while(monkey_count[0] < pntr[monkey_tail[0]+1]){
                        if (monkey_count[current] < pntr[monkey_tail[current]+1]){
                            offset = indx[monkey_count[current]];
                            evalBasisSupported(work->getIndex(offset), this_x, isSupported);
                            if (isSupported){
                                c++;
                                monkey_tail[++current] = offset;
                                monkey_count[current] = pntr[offset];
                            }else{
                                monkey_count[current]++;
                            }
                        }else{
                            monkey_count[--current]++;
                        }
                    }
                }
            }
        }

        total_nz = c;

        delete[] monkey_count;
        delete[] monkey_tail;
    }

    return total_nz;
}
void GridLocalPolynomial::buildSparseMatrixBlockForm(const double x[], int num_x, int num_chunk, int &num_blocks, int &num_last, int &stripe_size,
                                                    int* &stripes, int* &last_stripe_size, int** &tpntr, int*** &tindx, double*** &tvals) const{
    IndexSet *work = (points == 0) ? needed : points;
    num_blocks = (num_x / num_chunk) + ((num_x % num_chunk > 0) ? 1 : 0); // number of blocks
    num_last = ((num_x % num_chunk > 0) ? num_x % num_chunk : num_chunk);
    stripe_size = work->getNumIndexes();

    stripes = new int[num_blocks]; std::fill(stripes, stripes + num_blocks, 1);
    last_stripe_size = new int[num_blocks];
    tpntr = new int*[num_blocks];
    for(int i=0; i<num_blocks-1; i++) tpntr[i] = new int[num_chunk+1];
    tpntr[num_blocks-1] = new int[num_last+1];
    tindx = new int**[num_blocks];
    tvals = new double**[num_blocks];
    for(int b=0; b<num_blocks; b++){
        tindx[b] = new int*[1]; // start with only one stripe
        tindx[b][0] = new int[stripe_size];
        tvals[b] = new double*[1]; // start with only one stripe
        tvals[b][0] = new double[stripe_size];
    }

    #pragma omp parallel for
    for(int b=0; b<num_blocks; b++){
        int c = 0; // count the current entry
        int s = 0; // index the current stripe
        tpntr[b][0] = 0;

        int *monkey_count = new int[top_level+1]; // monkey business (traversing threes)
        int *monkey_tail = new int[top_level+1];
        bool isSupported;
        int offset;

        int this_size = (b == num_blocks-1) ? num_last : num_chunk;
        for(int i=0; i<this_size; i++){ // for each point
            tpntr[b][i+1] = tpntr[b][i];
            const double *this_x = &(x[(b*num_chunk + i) * num_dimensions]);

            for(int r=0; r<num_roots; r++){
                double basis_value = evalBasisSupported(work->getIndex(roots[r]), this_x, isSupported);

                if (isSupported){
                    // add point to this sparse column
                    if (c == stripe_size){ // add a stripe
                        stripes[b]++;
                        int **save_tindx = tindx[b];
                        double **save_tvals = tvals[b];
                        tindx[b] = new int*[stripes[b]];
                        tvals[b] = new double*[stripes[b]];
                        for(int j=0; j<stripes[b]-1; j++){
                            tindx[b][j] = save_tindx[j];
                            tvals[b][j] = save_tvals[j];
                        }
                        tindx[b][stripes[b]-1] = new int[stripe_size];
                        tvals[b][stripes[b]-1] = new double[stripe_size];
                        c = 0;
                        s++;

                        delete[] save_tindx;
                        delete[] save_tvals;
                    }
                    // add an entry
                    tindx[b][s][c] = roots[r];
                    tvals[b][s][c] = basis_value;
                    tpntr[b][i+1]++;
                    c++;

                    int current = 0;
                    monkey_tail[0] = roots[r];
                    monkey_count[0] = pntr[roots[r]];

                    while(monkey_count[0] < pntr[monkey_tail[0]+1]){
                        if (monkey_count[current] < pntr[monkey_tail[current]+1]){
                            offset = indx[monkey_count[current]];
                            basis_value = evalBasisSupported(work->getIndex(offset), this_x, isSupported);
                            if (isSupported){
                                if (c == stripe_size){ // add a stripe
                                    stripes[b]++;
                                    int **save_tindx = tindx[b];
                                    double **save_tvals = tvals[b];
                                    tindx[b] = new int*[stripes[b]];
                                    tvals[b] = new double*[stripes[b]];
                                    for(int j=0; j<stripes[b]-1; j++){
                                        tindx[b][j] = save_tindx[j];
                                        tvals[b][j] = save_tvals[j];
                                    }
                                    tindx[b][stripes[b]-1] = new int[stripe_size];
                                    tvals[b][stripes[b]-1] = new double[stripe_size];
                                    c = 0;
                                    s++;

                                    delete[] save_tindx;
                                    delete[] save_tvals;
                                }
                                // add an entry
                                tindx[b][s][c] = offset;
                                tvals[b][s][c] = basis_value;
                                tpntr[b][i+1]++;
                                c++;

                                monkey_tail[++current] = offset;
                                monkey_count[current] = pntr[offset];
                            }else{
                                monkey_count[current]++;
                            }
                        }else{
                            monkey_count[--current]++;
                        }
                    }
                }
            }

        }
        last_stripe_size[b] = c;

        delete[] monkey_count;
        delete[] monkey_tail;
    }
}

#ifdef Tasmanian_ENABLE_CUDA
void GridLocalPolynomial::buildDenseBasisMatrixGPU(const double gpu_x[], int cpu_num_x, double gpu_y[], std::ostream *os) const{
    int num_points = getNumPoints();
    makeCheckAccelerationData(accel_gpu_cuda, os);
    checkAccelerationGPUNodes();
    AccelerationDataGPUFull *gpu_acc = (AccelerationDataGPUFull*) accel;
    TasCUDA::devalpwpoly(order, rule->getType(), num_dimensions, cpu_num_x, num_points, gpu_x, gpu_acc->getGPUNodes(), gpu_acc->getGPUSupport(), gpu_y);
}
void GridLocalPolynomial::buildSparseBasisMatrixGPU(const double gpu_x[], int cpu_num_x, int* &gpu_spntr, int* &gpu_sindx, double* &gpu_svals, int &num_nz, std::ostream *os) const{
    int num_points = getNumPoints();
    makeCheckAccelerationData(accel_gpu_cuda, os);
    checkAccelerationGPUNodes();
    checkAccelerationGPUHierarchy();
    AccelerationDataGPUFull *gpu_acc = (AccelerationDataGPUFull*) accel;
    TasCUDA::devalpwpoly_sparse(order, rule->getType(), num_dimensions, cpu_num_x, num_points, gpu_x, gpu_acc->getGPUNodes(), gpu_acc->getGPUSupport(),
                                gpu_acc->getGPUpntr(), gpu_acc->getGPUindx(), num_roots, gpu_acc->getGPUroots(),
                                gpu_spntr, gpu_sindx, gpu_svals, num_nz, os);
}

#else
void GridLocalPolynomial::buildDenseBasisMatrixGPU(const double*, int, double*, std::ostream*) const{}
void GridLocalPolynomial::buildSparseBasisMatrixGPU(const double*, int, int*&, int*&, double*&, int&, std::ostream*) const{}
#endif // Tasmanian_ENABLE_CUDA

void GridLocalPolynomial::buildTree(){
    if (roots != 0) delete[] roots;
    if (pntr != 0) delete[] pntr;
    if (indx != 0) delete[] indx;

    IndexSet *work = (points == 0) ? needed : points;
    int num_points = work->getNumIndexes();

    int *level = new int[num_points];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = work->getIndex(i);
        level[i] = rule->getLevel(p[0]);
        for(int j=1; j<num_dimensions; j++){
            level[i] += rule->getLevel(p[j]);
        }
    }

    top_level = 0;
    for(int i=0; i<num_points; i++) if (top_level < level[i]) top_level = level[i];

    int max_1d_kids = rule->getMaxNumKids();
    int max_kids = max_1d_kids*num_dimensions;
    int* monkey_count = new int[top_level + 1];
    int* monkey_tail  = new int[top_level + 1];

    int  *tree = new int[max_kids * num_points];  std::fill(tree, tree + max_kids * num_points, -1);
    bool *free = new bool[num_points];        std::fill(free, free + num_points, true);
    int  *kid  = new int[num_dimensions];

    int next_root = 0;

    int crts = 12;
    int *rts = new int[crts];
    int num_rts = 0;

    while(next_root != -1){
        if (num_rts == crts){
            int *tmp = rts;
            crts *= 2;
            rts = new int[crts];
            std::copy(tmp, tmp + num_rts, rts);
            delete[] tmp;
        }
        rts[num_rts++] = next_root;
        free[next_root] = false;

        monkey_tail[0] = next_root;
        monkey_count[0] = 0;
        int current = 0;

        while(monkey_count[0] < max_kids){
            if (monkey_count[current] < max_kids){
                const int *p = work->getIndex(monkey_tail[current]);

                int dir = monkey_count[current] / max_1d_kids;
                int ikid = rule->getKid(p[dir], monkey_count[current] % max_1d_kids);

                if (ikid == -1){
                    monkey_count[current]++;
                }else{
                    std::copy(p, p + num_dimensions, kid);
                    kid[dir] = ikid;
                    int t = work->getSlot(kid);
                    if ((t == -1) || (!free[t])){
                        monkey_count[current]++;
                    }else{
                        tree[monkey_tail[current] * max_kids + monkey_count[current]] =  t;
                        monkey_count[++current] = 0;
                        monkey_tail[current] = t;
                        free[t] = false;
                    }
                }
            }else{
                monkey_count[--current]++;
            }
        }

        next_root = -1;
        int next_level = top_level+1;
        for(int i=0; i<num_points; i++){
            if (free[i] && (level[i] < next_level)){
                next_root = i;
                next_level = level[i];
            }
        }
    }

    num_roots = num_rts;
    roots = new int[num_roots];
    std::copy(rts, rts + num_roots, roots);

    pntr = new int[num_points + 1];
    pntr[0] = 0;
    for(int i=0; i<num_points; i++){
        pntr[i+1] = pntr[i];
        for(int j=0; j<max_kids; j++) if (tree[i * max_kids + j] > -1) pntr[i+1]++;
    }

    indx = (pntr[num_points] > 0) ? new int[pntr[num_points]] : new int[1]; indx[0] = 0;
    int count = 0;
    for(int i=0; i<num_points; i++){
        for(int j=0; j<max_kids; j++){
            int t = tree[i * max_kids + j];
            if (t > -1) indx[count++] = t;
        }
    }

    delete[] kid;
    delete[] free;
    delete[] tree;
    delete[] monkey_tail;
    delete[] monkey_count;
    delete[] rts;
    delete[] level;
}

void GridLocalPolynomial::getBasisIntegrals(double *integrals) const{
    IndexSet *work = (points == 0) ? needed : points;

    int n = 0;
    double *w = 0, *x = 0;

    if ((rule->getMaxOrder() == -1) || (rule->getMaxOrder() > 3) ){
        n = top_level / 2 + 1;
        OneDimensionalNodes::getGaussLegendre(n, w, x);
    }

    for(int i=0; i<work->getNumIndexes(); i++){
        const int* p = work->getIndex(i);
        integrals[i] = rule->getArea(p[0], n, w, x);
        for(int j=1; j<num_dimensions; j++){
            integrals[i] *= rule->getArea(p[j], n, w, x);
        }
    }

    if ((rule->getMaxOrder() == -1) || (rule->getMaxOrder() > 3)){
        delete[] x;
        delete[] w;
    }
}
double* GridLocalPolynomial::getQuadratureWeights() const{
    IndexSet *work = (points == 0) ? needed : points;
    double *weights = new double[work->getNumIndexes()];
    getQuadratureWeights(weights);
    return weights;
}
void GridLocalPolynomial::getQuadratureWeights(double *weights) const{
    IndexSet *work = (points == 0) ? needed : points;
    getBasisIntegrals(weights);

    int *monkey_count = new int[top_level+1];
    int *monkey_tail = new int[top_level+1];

    double basis_value;

    const int *dagUp;
    if (num_outputs == 0){
        dagUp = parents;
    }else{
        IndexManipulator IM(num_dimensions);
        dagUp = IM.computeDAGupLocal(work, rule);
    }

    int num_points = work->getNumIndexes();
    int *level = new int[num_points];
    for(int i=0; i<num_points; i++){
        const int *p = work->getIndex(i);
        level[i] = rule->getLevel(p[0]);
        for(int j=1; j<num_dimensions; j++){
            level[i] += rule->getLevel(p[j]);
        }
    }

    bool *used = new bool[work->getNumIndexes()];
    double *node = new double[num_dimensions];
    int max_parents = rule->getMaxNumParents() * num_dimensions;

    for(int l=top_level; l>0; l--){
        for(int i=0; i<num_points; i++){
            if (level[i] == l){
                const int* p = work->getIndex(i);
                for(int j=0; j<num_dimensions; j++) node[j] = rule->getNode(p[j]);

                std::fill(used, used + work->getNumIndexes(), false);

                monkey_count[0] = 0;
                monkey_tail[0] = i;
                int current = 0;

                while(monkey_count[0] < max_parents){
                    if (monkey_count[current] < max_parents){
                        int branch = dagUp[ monkey_tail[current] * max_parents + monkey_count[current] ];
                        if ((branch == -1) || used[branch]){
                            monkey_count[current]++;
                        }else{
                            const int *func = work->getIndex(branch);
                            basis_value = rule->evalRaw(func[0], node[0]);
                            for(int j=1; j<num_dimensions; j++) basis_value *= rule->evalRaw(func[j], node[j]);
                            weights[branch] -= weights[i] * basis_value;
                            used[branch] = true;

                            monkey_count[++current] = 0;
                            monkey_tail[current] = branch;
                        }
                    }else{
                        monkey_count[--current]++;
                    }
                }
            }
        }
    }

    delete[] used;
    delete[] node;

    delete[] level;

    delete[] monkey_count;
    delete[] monkey_tail;

    if (num_outputs > 0)  delete[] dagUp;
}

void GridLocalPolynomial::integrate(double q[], double *conformal_correction) const{
    int num_points = points->getNumIndexes();
    std::fill(q, q + num_outputs, 0.0);

    if (conformal_correction == 0){
        double *integrals = new double[num_points];
        getBasisIntegrals(integrals);
        for(int i=0; i<num_points; i++){
            for(int k=0; k<num_outputs; k++){
                q[k] += integrals[i] * surpluses[i*num_outputs + k];
            }
        }
        delete[] integrals;
    }else{
        double *w = getQuadratureWeights();
        for(int i=0; i<num_points; i++){
            w[i] *= conformal_correction[i];
            const double *vals = values->getValues(i);
            for(int k=0; k<num_outputs; k++){
                q[k] += w[i] * vals[k];
            }
        }
        delete[] w;
    }
}

double* GridLocalPolynomial::getNormalization() const{
    double* norm = new double[num_outputs];  std::fill(norm, norm + num_outputs, 0.0);
    for(int i=0; i<points->getNumIndexes(); i++){
        const double *v = values->getValues(i);
        for(int j=0; j<num_outputs; j++){
            if (norm[j] < fabs(v[j])) norm[j] = fabs(v[j]);
        }
    }
    return norm;
}

int* GridLocalPolynomial::buildUpdateMap(double tolerance, TypeRefinement criteria, int output, const double *scale_correction) const{
    int num_points = points->getNumIndexes();
    int *pmap = new int[num_points * num_dimensions];  std::fill(pmap, pmap + num_points * num_dimensions, 0);

    double *norm = getNormalization();

    const double *scale = scale_correction;
    double *temp_scale = 0;
    if (scale == 0){
        size_t size_scale = ((size_t) num_points) * ((size_t) num_outputs);
        temp_scale = new double[size_scale]; std::fill(temp_scale, temp_scale + size_scale, 1.0);
        scale = temp_scale;
    }

    if (tolerance == 0.0){
        std::fill(pmap, pmap + num_points * num_dimensions, 1);
    }else if ((criteria == refine_classic) || (criteria == refine_parents_first)){
        #pragma omp parallel for
        for(int i=0; i<num_points; i++){
            bool small = true;
            if (output == -1){
                for(size_t k=0; k<((size_t) num_outputs); k++){
                    if (small && ((scale[((size_t)i) * ((size_t) num_outputs) + k] * fabs(surpluses[((size_t)i) * ((size_t) num_outputs) + k]) / norm[k]) > tolerance)) small = false;
                }
            }else{
                small = !((scale[i] * fabs(surpluses[((size_t)i) * ((size_t) num_outputs) + output]) / norm[output]) > tolerance);
            }
            if (!small){
                std::fill(&(pmap[i*num_dimensions]), &(pmap[i*num_dimensions]) + num_dimensions, 1);
            }
        }
    }else{
        // construct a series of 1D interpolants and use a refinement criteria that is a combination of the two hierarchical coefficients
        IndexManipulator IM(num_dimensions);
        int *dagUp = IM.computeDAGupLocal(points, rule);

        int max_1D_parents = rule->getMaxNumParents();
        int max_parents = max_1D_parents * num_dimensions;

        SplitDirections split(points);

        #pragma omp parallel for
        for(int s=0; s<split.getNumJobs(); s++){ // split.getNumJobs() gives the number of 1D interpolants to construct
            int d = split.getJobDirection(s);
            int nump = split.getJobNumPoints(s);
            const int *pnts = split.getJobPoints(s);

            int *global_to_pnts = new int[num_points];
            int *levels = new int[nump];
            int max_level = 0;

            int active_outputs = (output == -1) ? num_outputs : 1;

            double *vals = new double[((size_t) nump) * ((size_t) active_outputs)];

            for(int i=0; i<nump; i++){
                const double* v = values->getValues(pnts[i]);
                const int *p = points->getIndex(pnts[i]);
                if (output == -1){
                    std::copy(v, v + num_outputs, &(vals[((size_t) i) * ((size_t) num_outputs)]));
                }else{
                    vals[i] = v[output];
                }
                global_to_pnts[pnts[i]] = i;
                levels[i] = rule->getLevel(p[d]);
                if (max_level < levels[i]) max_level = levels[i];
            }

            int *monkey_count = new int[max_level + 1];
            int *monkey_tail = new int[max_level + 1];
            bool *used = new bool[nump];

            for(int l=1; l<=max_level; l++){
                for(int i=0; i<nump; i++){
                    if (levels[i] == l){
                        const int *p = points->getIndex(pnts[i]);
                        double x = rule->getNode(p[d]);
                        double *valsi = &(vals[((size_t) i) * ((size_t) active_outputs)]);

                        int current = 0;
                        monkey_count[0] = d * max_1D_parents;
                        monkey_tail[0] = pnts[i]; // uses the global indexes
                        std::fill(used, used + nump, false);

                        while(monkey_count[0] < (d+1) * max_1D_parents){
                            if (monkey_count[current] < (d+1) * max_1D_parents){
                                int branch = dagUp[monkey_tail[current] * max_parents + monkey_count[current]];
                                if ((branch == -1) || (used[global_to_pnts[branch]])){
                                    monkey_count[current]++;
                                }else{
                                    const int *branch_point = points->getIndex(branch);
                                    double basis_value = rule->evalRaw(branch_point[d], x);
                                    const double *branch_vals = &(vals[((size_t) global_to_pnts[branch]) * ((size_t) active_outputs)]);
                                    for(int k=0; k<active_outputs; k++){
                                        valsi[k] -= basis_value * branch_vals[k];
                                    }

                                    used[global_to_pnts[branch]] = true;
                                    monkey_count[++current] = d * max_1D_parents;
                                    monkey_tail[current] = branch;
                                }
                            }else{
                                monkey_count[--current]++;
                            }
                        }
                    }
                }
            }

            delete[] used;
            delete[] monkey_tail;
            delete[] monkey_count;

            // at this point, vals contains the one directional surpluses
            for(size_t i=0; i<((size_t) nump); i++){
                bool small = true;
                if (output == -1){
                    for(size_t k=0; k<((size_t) num_outputs); k++){
                        if (small && ((scale[i*num_outputs + k] * fabs(surpluses[pnts[i]*num_outputs + k]) / norm[k]) > tolerance) && ((scale[i*num_outputs + k] * fabs(vals[i*num_outputs + k]) / norm[k]) > tolerance)) small = false;
                    }
                }else{
                    if (((scale[i] * fabs(surpluses[pnts[i]*num_outputs + output]) / norm[output]) > tolerance) && ((scale[i] * fabs(vals[i]) / norm[output]) > tolerance)) small = false;
                }
                pmap[pnts[i]*num_dimensions + d] = (small) ? 0 : 1;;
            }

            delete[] vals;
            delete[] global_to_pnts;
            delete[] levels;
        }

        delete[] dagUp;
    }

    delete[] norm;
    if (temp_scale != 0) delete[] temp_scale;

    return pmap;
}

bool GridLocalPolynomial::addParent(const int point[], int direction, GranulatedIndexSet *destination, IndexSet *exclude) const{
    int *dad = new int[num_dimensions];  std::copy(point, point + num_dimensions, dad);
    bool added = false;
    dad[direction] = rule->getParent(point[direction]);
    if ((dad[direction] != -1) && (exclude->getSlot(dad) == -1)){
        destination->addIndex(dad);
        added = true;
    }
    dad[direction] = rule->getStepParent(point[direction]);
    if ((dad[direction] != -1) && (exclude->getSlot(dad) == -1)){
        destination->addIndex(dad);
        added = true;
    }
    delete[] dad;
    return added;
}
void GridLocalPolynomial::addChild(const int point[], int direction, GranulatedIndexSet *destination, IndexSet *exclude)const{
    int *kid = new int[num_dimensions];  std::copy(point, point + num_dimensions, kid);
    int max_1d_kids = rule->getMaxNumKids();
    for(int i=0; i<max_1d_kids; i++){
        kid[direction] = rule->getKid(point[direction], i);
        if ((kid[direction] != -1) && (exclude->getSlot(kid) == -1)){
            destination->addIndex(kid);
        }
    }
    delete[] kid;
}
void GridLocalPolynomial::addChildLimited(const int point[], int direction, GranulatedIndexSet *destination, IndexSet *exclude, const int *level_limits) const{
    int *kid = new int[num_dimensions];  std::copy(point, point + num_dimensions, kid);
    int max_1d_kids = rule->getMaxNumKids();
    for(int i=0; i<max_1d_kids; i++){
        kid[direction] = rule->getKid(point[direction], i);
        if ((kid[direction] != -1) && (rule->getLevel(kid[direction]) <= level_limits[direction]) && (exclude->getSlot(kid) == -1)){
            destination->addIndex(kid);
        }
    }
    delete[] kid;
}

void GridLocalPolynomial::clearRefinement(){
    if (needed != 0){ delete needed; needed = 0; }
}
const double* GridLocalPolynomial::getSurpluses() const{
    return surpluses;
}
const int* GridLocalPolynomial::getPointIndexes() const{
    return ((points == 0) ? needed->getIndex(0) : points->getIndex(0));
}
const int* GridLocalPolynomial::getNeededIndexes() const{
    //cout << "HERE " << needed->getNumIndexes() << endl;
    return ((needed != 0) ? needed->getIndex(0) : 0);
}
void GridLocalPolynomial::setSurplusRefinement(double tolerance, TypeRefinement criteria, int output, const int *level_limits, const double *scale_correction){
    clearRefinement();

    int *pmap = buildUpdateMap(tolerance, criteria, output, scale_correction);

    bool useParents = (criteria == refine_fds) || (criteria == refine_parents_first);

    GranulatedIndexSet *refined = new GranulatedIndexSet(num_dimensions);

    int num_points = points->getNumIndexes();

    if (level_limits == 0){
        for(int i=0; i<num_points; i++){
            for(int j=0; j<num_dimensions; j++){
                if (pmap[i*num_dimensions+j] == 1){ // if this dimension needs to be refined
                    if (!(useParents && addParent(points->getIndex(i), j, refined, points))){
                        addChild(points->getIndex(i), j, refined, points);
                    }
                }
            }
        }
    }else{
        for(int i=0; i<num_points; i++){
            for(int j=0; j<num_dimensions; j++){
                if (pmap[i*num_dimensions+j] == 1){ // if this dimension needs to be refined
                    if (!(useParents && addParent(points->getIndex(i), j, refined, points))){
                        addChildLimited(points->getIndex(i), j, refined, points, level_limits);
                    }
                }
            }
        }
    }

    //////////// TEST CODE /////////////////////
    // complete the set of points to ensure no missing parents
    // often times it results in real excess of points,
    // but I had to test
//    if (refined->getNumIndexes() > 0){
//        IndexManipulator IM(num_dimensions);
//        IndexSet* total = new IndexSet(refined);
//        IndexSet* completion = IM.getLowerCompletion(total);
//        cout << "total = " << total->getNumIndexes() << "  " << completion->getNumIndexes() << endl;
//        if ((completion != 0) && (completion->getNumIndexes() > 0)){
//            total->addIndexSet(completion);
//            delete completion;
//        }
//        needed = total->diffSets(points);
//        cout << " needed = " << needed->getNumIndexes() << endl;
//    }
    //////// END TEST CODE /////////////////
    if (refined->getNumIndexes() > 0){
        needed = new IndexSet(refined);
    }
    delete refined;

    delete[] pmap;
}
int GridLocalPolynomial::removePointsByHierarchicalCoefficient(double tolerance, int output, const double *scale_correction){
    clearRefinement();
    int num_points = points->getNumIndexes();
    bool *pmap = new bool[num_points]; // point map, set to true if the point is to be kept, false otherwise

    double *norm = getNormalization();

    const double *scale = scale_correction;
    double *temp_scale = 0;
    if (scale == 0){
        temp_scale = new double[num_points * num_outputs]; std::fill(temp_scale, temp_scale + num_points * num_outputs, 1.0);
        scale = temp_scale;
    }

    for(int i=0; i<num_points; i++){
        bool small = true;
        if (output == -1){
            for(int k=0; k<num_outputs; k++){
                if (small && ((scale[i*num_outputs + k] * fabs(surpluses[i*num_outputs + k]) / norm[k]) > tolerance)) small = false;
            }
        }else{
            small = !((scale[i] * fabs(surpluses[i*num_outputs + output]) / norm[output]) > tolerance);
        }
        pmap[i] = !small;
    }

    int num_kept = 0; for(int i=0; i<num_points; i++) num_kept += (pmap[i]) ? 1 : 0;
    if (num_kept == 0){
        delete[] norm;
        delete[] pmap;
        return 0;
    }

    if (num_kept == num_points){
        delete[] norm;
        delete[] pmap;
        return num_points;
    }

    // save a copy of the points and the values
    int *point_kept = new int[num_kept * num_dimensions];
    double *values_kept = new double[num_kept * num_outputs];

    num_kept = 0;
    for(int i=0; i<num_points; i++){
        if (pmap[i]){
            std::copy(points->getIndex(i), points->getIndex(i) + num_dimensions, &(point_kept[num_kept*num_dimensions]));
            std::copy(values->getValues(i), values->getValues(i) + num_outputs, &(values_kept[num_kept*num_outputs]));
            num_kept++;
        }
    }

    int dims = num_dimensions, outs = num_outputs;

    reset(false);
    num_dimensions = dims; num_outputs = outs;

    points = new IndexSet(num_dimensions, num_kept, point_kept);
    values = new StorageSet(num_outputs, num_kept);
    values->setValues(values_kept);
    delete[] values_kept;

    buildTree();
    recomputeSurpluses();

    delete[] norm;
    delete[] pmap;
    if (temp_scale == 0) delete[] temp_scale;

    return points->getNumIndexes();
}

void GridLocalPolynomial::setHierarchicalCoefficients(const double c[], TypeAcceleration acc, std::ostream *os){
    if (accel != 0) accel->resetGPULoadedData();
    double *vals = 0;
    bool aliased = false;
    if (points != 0){
        clearRefinement();
        vals = values->aliasValues();
        aliased = true;
    }else{
        points = needed;
        needed = 0;
        vals = new double[((size_t) points->getNumIndexes()) * ((size_t) num_outputs)];
    }
    int num_ponits = points->getNumIndexes();
    if (surpluses != 0) delete[] surpluses;
    surpluses = new double[((size_t) num_ponits) * ((size_t) num_outputs)];
    std::copy(c, c + ((size_t) num_ponits) * ((size_t) num_outputs), surpluses);
    double *x = getPoints();
    if (acc == accel_cpu_blas){
        evaluateBatchCPUblas(x, points->getNumIndexes(), vals);
    }else if (acc == accel_gpu_cublas){
        evaluateBatchGPUcublas(x, points->getNumIndexes(), vals, os);
    }else if (acc == accel_gpu_cuda){
        evaluateBatchGPUcuda(x, points->getNumIndexes(), vals, os);
    }else{
        evaluateBatch(x, points->getNumIndexes(), vals);
    }
    delete[] x;
    if (! aliased){
        values->setValuesPointer(vals, num_ponits);
    }
}

#if defined(Tasmanian_ENABLE_CUBLAS) || defined(Tasmanian_ENABLE_CUDA)
void GridLocalPolynomial::makeCheckAccelerationData(TypeAcceleration acc, std::ostream *os) const{
    if (AccelerationMeta::isAccTypeFullMemoryGPU(acc)){
        if ((accel != 0) && (!accel->isCompatible(acc))){
            delete accel;
            accel = 0;
        }
        if (accel == 0){ accel = (BaseAccelerationData*) (new AccelerationDataGPUFull()); }
        ((AccelerationDataGPUFull*) accel)->setLogStream(os);
    }
}
void GridLocalPolynomial::checkAccelerationGPUValues() const{
    AccelerationDataGPUFull *gpu = (AccelerationDataGPUFull*) accel;
    double *gpu_values = gpu->getGPUValues();
    if (gpu_values == 0) gpu->loadGPUValues(points->getNumIndexes() * values->getNumOutputs(), surpluses);
}
void GridLocalPolynomial::checkAccelerationGPUNodes() const{
    AccelerationDataGPUFull *gpu = (AccelerationDataGPUFull*) accel;
    if (gpu->getGPUNodes() == 0){
        IndexSet *work = (points == 0) ? needed : points;
        int num_entries = num_dimensions * work->getNumIndexes();
        double *cpu_nodes = getPoints();
        double *cpu_support = new double[num_entries];
        if (rule->getType() == rule_localp){
            switch(order){
            case 0: encodeSupportForGPU<0, rule_localp>(work, cpu_support); break;
            case 2: encodeSupportForGPU<2, rule_localp>(work, cpu_support); break;
            default:
                encodeSupportForGPU<1, rule_localp>(work, cpu_support);
            }
        }else if (rule->getType() == rule_semilocalp){
            encodeSupportForGPU<2, rule_semilocalp>(work, cpu_support);
        }else{
            switch(order){
            case 2: encodeSupportForGPU<2, rule_localp0>(work, cpu_support); break;
            default:
                encodeSupportForGPU<1, rule_localp0>(work, cpu_support);
            }
        }
        gpu->loadGPUNodesSupport(num_entries, cpu_nodes, cpu_support);
        delete[] cpu_nodes;
        delete[] cpu_support;
    }
}
void GridLocalPolynomial::checkAccelerationGPUHierarchy() const{
    AccelerationDataGPUFull *gpu = (AccelerationDataGPUFull*) accel;
    int num_points = getNumPoints();
    gpu->loadGPUHierarchy(num_points, pntr, indx, num_roots, roots);
}
#else
void GridLocalPolynomial::makeCheckAccelerationData(TypeAcceleration, std::ostream *) const{}
void GridLocalPolynomial::checkAccelerationGPUValues() const{}
void GridLocalPolynomial::checkAccelerationGPUNodes() const{}
void GridLocalPolynomial::checkAccelerationGPUHierarchy() const{}
#endif // Tasmanian_ENABLE_CUBLAS

void GridLocalPolynomial::clearAccelerationData(){
    if (accel != 0){
        delete accel;
        accel = 0;
    }
}

}

#endif
