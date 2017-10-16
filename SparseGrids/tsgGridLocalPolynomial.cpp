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

namespace TasGrid{

GridLocalPolynomial::GridLocalPolynomial() : num_dimensions(0), num_outputs(0), order(1), top_level(0),
                         surpluses(0), points(0), needed(0), values(0), parents(0), num_roots(0), roots(0), pntr(0), indx(0), rule(0),
                         accel(0)  {}
GridLocalPolynomial::GridLocalPolynomial(const GridLocalPolynomial &pwpoly) : num_dimensions(0), num_outputs(0), order(1), top_level(0),
                         surpluses(0), points(0), needed(0), values(0), parents(0), num_roots(0), roots(0), pntr(0), indx(0), rule(0),
                         accel(0)
{
    copyGrid(&pwpoly);
}

GridLocalPolynomial::~GridLocalPolynomial(){
    reset();
}
void GridLocalPolynomial::reset(bool clear_rule){
    clearAccelerationData();
    num_dimensions = num_outputs = top_level = 0;
    if (surpluses != 0){  delete[] surpluses; surpluses = 0;  }
    if (points != 0){  delete points; points = 0;  }
    if (needed != 0){  delete needed; needed = 0;  }
    if (values != 0){  delete values; values = 0;  }
    if (parents != 0){  delete[] parents; parents = 0;  }
    if (roots != 0){  delete[] roots;  roots = 0;  }
    if (pntr != 0){  delete[] pntr;  pntr = 0;  }
    if (indx != 0){  delete[] indx;  indx = 0;  }
    if (clear_rule){ rule = 0; order = 1; }
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
            for(int i=0; i<points->getNumIndexes() * num_outputs; i++){  ofs << " " << surpluses[i];  } ofs << endl;
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
            //if (rule->isSemiLocal()){
            //    for(int i=0; i<2*(points->getNumIndexes()); i++){ ofs << " " << parents[i]; } ofs << endl;
            //}else{
            //    for(int i=0; i<points->getNumIndexes(); i++){ ofs << " " << parents[i]; } ofs << endl;
            //}
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
            ofs.write((char*) surpluses, num_outputs * points->getNumIndexes() * sizeof(double));
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

        ifs >> flag;  if (flag == 1){  points    = new IndexSet(num_dimensions);   points->read(ifs);  }
        ifs >> flag;  if (flag == 1){  surpluses = new double[points->getNumIndexes() * num_outputs]; for(int i=0; i<points->getNumIndexes() * num_outputs; i++){  ifs >> surpluses[i];  }  }
        ifs >> flag;  if (flag == 1){  needed    = new IndexSet(num_dimensions);   needed->read(ifs);  }
        //ifs >> flag;  if (flag == 1){  int num_parents = ((rule->isSemiLocal()) ? 2 : 1) * points->getNumIndexes();  parents = new int[num_parents];  for(int i=0; i<num_parents; i++){  ifs >> parents[i];  }  }
        ifs >> flag;  if (flag == 1){  int num_parents = rule->getMaxNumParents() * points->getNumIndexes();  parents = new int[num_parents];  for(int i=0; i<num_parents; i++){  ifs >> parents[i];  }  }

        int num_points = (points == 0) ? needed->getNumIndexes() : points->getNumIndexes();

        ifs >> num_roots;
        roots = new int[num_roots];     for(int i=0; i<num_roots; i++) ifs >> roots[i];
        pntr  = new int[num_points + 1];    for(int i=0; i<=num_points; i++) ifs >> pntr[i];
        if (pntr[num_points] > 0){
            indx  = new int[pntr[num_points]];  for(int i=0; i<pntr[num_points]; i++) ifs >> indx[i];
        }else{
            indx  = new int[1];  ifs >> indx[0]; // there is a special case when the grid has only one point without any children
        }

        if (num_outputs > 0){  values = new StorageSet(0, 0); values->read(ifs);  }
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

        ifs.read((char*) &flag, sizeof(char)); if (flag == 'y'){ surpluses = new double[num_outputs * points->getNumIndexes()]; ifs.read((char*) surpluses, num_outputs * points->getNumIndexes() * sizeof(double)); }
        ifs.read((char*) &flag, sizeof(char)); if (flag == 'y'){ parents = new int[rule->getMaxNumParents() * points->getNumIndexes()]; ifs.read((char*) parents, rule->getMaxNumParents() * points->getNumIndexes() * sizeof(int)); }

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

void GridLocalPolynomial::makeGrid(int cnum_dimensions, int cnum_outputs, int depth, int corder, TypeOneDRule crule){
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
        recomputeSurpluses();
    }
}

int GridLocalPolynomial::getNumDimensions() const{  return num_dimensions;  }
int GridLocalPolynomial::getNumOutputs() const{  return num_outputs;  }
TypeOneDRule GridLocalPolynomial::getRule() const{  return rule->getType();  }
int GridLocalPolynomial::getOrder() const{  return order;  }

int GridLocalPolynomial::getNumLoaded() const{  return ((points == 0) ? 0 : points->getNumIndexes());  }
int GridLocalPolynomial::getNumNeeded() const{  return ((needed == 0) ? 0 : needed->getNumIndexes());  }
int GridLocalPolynomial::getNumPoints() const{  return ((points == 0) ? getNumNeeded() : getNumLoaded());  }

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
    int offset;

    std::fill(y, y + num_outputs, 0.0);

    //int count_support = 0;
    //int count_tested = 0;

    for(int r=0; r<num_roots; r++){
        double basis_value = evalBasisSupported(points->getIndex(roots[r]), x, isSupported);

        if (isSupported){
            offset = roots[r] * num_outputs;
            for(int k=0; k<num_outputs; k++) y[k] += basis_value * surpluses[offset + k];
            //cout << roots[r] << "   " << basis_value << endl;

            int current = 0;
            monkey_tail[0] = roots[r];
            monkey_count[0] = pntr[roots[r]];

            while(monkey_count[0] < pntr[monkey_tail[0]+1]){
                if (monkey_count[current] < pntr[monkey_tail[current]+1]){
                    //count_tested++;
                    offset = indx[monkey_count[current]];
                    basis_value = evalBasisSupported(points->getIndex(offset), x, isSupported);
                    if (isSupported){
                        //count_support++;
                        offset *= num_outputs;
                        for(int k=0; k<num_outputs; k++) y[k] += basis_value * surpluses[offset + k];
                        offset /= num_outputs;
                        //cout << offset << "   " << basis_value << endl;

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

    //cout << count_tested << "   " << count_support << endl;

    delete[] monkey_count;
    delete[] monkey_tail;
}

void GridLocalPolynomial::evaluateFastCPUblas(const double x[], double y[]) const{ evaluate(x, y); }
void GridLocalPolynomial::evaluateFastGPUcublas(const double x[], double y[], std::ostream *os) const{ evaluate(x, y); }
void GridLocalPolynomial::evaluateFastGPUcuda(const double x[], double y[], std::ostream *os) const{ evaluate(x, y); }
void GridLocalPolynomial::evaluateFastGPUmagma(const double x[], double y[], std::ostream *os) const{ evaluate(x, y); }

void GridLocalPolynomial::evaluateBatch(const double x[], int num_x, double y[]) const{
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        evaluate(&(x[i*num_dimensions]), &(y[i*num_outputs]));
    }
}
void GridLocalPolynomial::evaluateBatchCPUblas(const double x[], int num_x, double y[]) const{
    evaluateBatch(x, num_x, y);
    return;

    int *sindx, *spntr;
    double *svals;
    buildSpareBasisMatrix(x, num_x, 32, spntr, sindx, svals); // build sparse matrix corresponding to x

    // how do you optimize sparse BLAS? This is slower
    for(int i=0; i<num_x; i++){
        double *this_y = &(y[i*num_outputs]);
        std::fill(this_y, this_y + num_outputs, 0.0);
        for(int j=spntr[i]; j<spntr[i+1]; j++){
            TasBLAS::daxpy(num_outputs, svals[j], &(surpluses[sindx[j] * num_outputs]), this_y);
        }
    }
}
void GridLocalPolynomial::evaluateBatchGPUcublas(const double x[], int num_x, double y[], std::ostream *os) const{
    //evaluateBatchBLAS(x, num_x, y);
    #ifdef TASMANIAN_CUBLAS
    int num_points = points->getNumIndexes();
    makeCheckAccelerationData(accel_gpu_cublas, os);
    AccelerationDataGPUFull *gpu_acc = (AccelerationDataGPUFull*) accel;

    int *sindx, *spntr;
    double *svals;
    buildSpareBasisMatrix(x, num_x, 32, spntr, sindx, svals); // build sparse matrix corresponding to x

    //TasCUDA::d3gecs(num_outputs, num_x, num_points, gpu_acc->getGPUValues(), spntr, sindx, svals, y, &cerr);
    gpu_acc->cusparseDCRMM2(num_points, num_outputs, num_x, spntr, sindx, svals, y);

    delete[] svals;
    delete[] sindx;
    delete[] spntr;
    #else
    evaluateBatchCPUblas(x, num_x, y);
    #endif // TASMANIAN_CUDA
}
void GridLocalPolynomial::evaluateBatchGPUcuda(const double x[], int num_x, double y[], std::ostream *os) const{
    #ifdef TASMANIAN_CUDA
    int num_points = points->getNumIndexes();
    makeCheckAccelerationData(accel_gpu_cuda, os);
    AccelerationDataGPUFull *gpu_acc = (AccelerationDataGPUFull*) accel;

    int *sindx, *spntr;
    double *svals;
    buildSpareBasisMatrix(x, num_x, 32, spntr, sindx, svals); // build sparse matrix corresponding to x

    TasCUDA::d3gecs(num_outputs, num_x, num_points, gpu_acc->getGPUValues(), spntr, sindx, svals, y, &cerr);
    //gpu_acc->cusparseDCRMM2(num_points, num_outputs, num_x, spntr, sindx, svals, y);

    delete[] svals;
    delete[] sindx;
    delete[] spntr;
    #else
    evaluateBatchGPUcublas(x, num_x, y, os);
    #endif // TASMANIAN_CUDA
}
void GridLocalPolynomial::evaluateBatchGPUmagma(const double x[], int num_x, double y[], std::ostream *os) const{
    evaluateBatchGPUcublas(x, num_x, y, os);
}

void GridLocalPolynomial::loadNeededPoints(const double *vals){
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
    recomputeSurpluses();
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
    //int max_parents = (rule->isSemiLocal()) ? 2*num_dimensions : num_dimensions;
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

void GridLocalPolynomial::recomputeSurpluses(){
    int num_ponits = points->getNumIndexes();
    if (surpluses != 0) delete[] surpluses;
    surpluses = new double[num_ponits * num_outputs];

    if (num_outputs > 2){
        #pragma omp parallel for schedule(static)
        for(int i=0; i<num_ponits; i++){
            const double* v = values->getValues(i);
            std::copy(v, v + num_outputs, &(surpluses[i*num_outputs]));
        }
    }else{
        const double* v = values->getValues(0);
        std::copy(v, v + num_ponits * num_outputs, surpluses);
    }

    IndexManipulator IM(num_dimensions);
    int *dagUp = IM.computeDAGupLocal(points, rule);

    //int max_parents = (rule->isSemiLocal()) ? 2*num_dimensions : num_dimensions;
    int max_parents = rule->getMaxNumParents() * num_dimensions;

    int *level = new int[num_ponits];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_ponits; i++){
        const int *p = points->getIndex(i);
        level[i] = rule->getLevel(p[0]);
        for(int j=1; j<num_dimensions; j++){
            level[i] += rule->getLevel(p[j]);
        }
    }

    for(int l=1; l<=top_level; l++){
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<num_ponits; i++){
            if (level[i] == l){
                const int* p = points->getIndex(i);
                double *x = new double[num_dimensions];
                for(int j=0; j<num_dimensions; j++) x[j] = rule->getNode(p[j]);

                int *monkey_count = new int[top_level + 1];
                int *monkey_tail = new int[top_level + 1];
                bool *used = new bool[num_ponits];
                std::fill(used, used + num_ponits, false);

                int current = 0;

                monkey_count[0] = 0;
                monkey_tail[0] = i;

                while(monkey_count[0] < max_parents){
                    if (monkey_count[current] < max_parents){
                        int branch = dagUp[monkey_tail[current] * max_parents + monkey_count[current]];
                        if ((branch == -1) || (used[branch])){
                            monkey_count[current]++;
                        }else{
                            double basis_value = evalBasisRaw(points->getIndex(branch), x);
                            for(int k=0; k<num_outputs; k++){
                                surpluses[i*num_outputs + k] -= basis_value * surpluses[branch * num_outputs + k];
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
    int num_blocks = (num_x / num_chunk) + ((num_x % num_chunk > 0) ? 1 : 0); // number of blocks
    int num_last = ((num_x % num_chunk > 0) ? num_x % num_chunk : num_chunk);
    int stripe_size = points->getNumIndexes();

    int *stripes = new int[num_blocks]; std::fill(stripes, stripes + num_blocks, 1);
    int *last_stripe_size = new int[num_blocks];
    int ** tpntr = new int*[num_blocks];
    for(int i=0; i<num_blocks-1; i++) tpntr[i] = new int[num_chunk+1];
    tpntr[num_blocks-1] = new int[num_last+1];
    int ***tindx = new int**[num_blocks];
    double ***tvals = new double**[num_blocks];
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
                double basis_value = evalBasisSupported(points->getIndex(roots[r]), this_x, isSupported);

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
                            //count_tested++;
                            offset = indx[monkey_count[current]];
                            basis_value = evalBasisSupported(points->getIndex(offset), this_x, isSupported);
                            if (isSupported){
                                //count_support++;
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

    // assemble the matrix from the local junks
    spntr = new int[num_x+1];
    spntr[0] = 0;
    int c = 0, block_offset = 0;
    for(int b=0; b<num_blocks; b++){
        int this_size = (b == num_blocks-1) ? num_last : num_chunk;
        //cout << "tpn 0 = " << tpntr[b][0] << endl;
        for(int i=0; i<this_size; i++){
            spntr[c+1] = block_offset + tpntr[b][i+1];
            //cout << "tpn 0 = " << tpntr[b][i+1] << endl;
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

    //int max_kids = 2*num_dimensions;
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

                //int dir = monkey_count[current] / 2;
                //int ikid = (monkey_count[current] % 2 == 0) ? rule->getKidLeft(p[dir]) : rule->getKidRight(p[dir]);
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
        OneDimensionalNodes GL;
        n = top_level / 2 + 1;
        GL.getGaussLegendre(n, w, x);
    }

    //double *integrals = new double[work->getNumIndexes()];
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
    //int max_parents = (rule->isSemiLocal()) ? 2*num_dimensions : num_dimensions;
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

int* GridLocalPolynomial::buildUpdateMap(double tolerance, TypeRefinement criteria, int output) const{
    int num_points = points->getNumIndexes();
    int *map = new int[num_points * num_dimensions];  std::fill(map, map + num_points * num_dimensions, 0);

    double *norm = getNormalization();

    // scaling for L^1 and L^2 norms
    //double *supports = new double[num_points];
    //for(int i=0; i<num_points; i++){
    //    const int *p = points->getIndex(i);
    //    supports[i] = rule->getSupport(p[0]);
    //    for(int j=1; j<num_dimensions; j++) supports[i] *= rule->getSupport(p[j]);
    //}

    if (tolerance == 0.0){
        std::fill(map, map + num_points * num_dimensions, 1);
    }else if ((criteria == refine_classic) || (criteria == refine_parents_first)){
        #pragma omp parallel for
        for(int i=0; i<num_points; i++){
            bool small = true;
            if (output == -1){
                for(int k=0; k<num_outputs; k++){
                    if (small && ((fabs(surpluses[i*num_outputs + k]) / norm[k]) > tolerance)) small = false;
                }
            }else{
                small = !((fabs(surpluses[i*num_outputs + output]) / norm[output]) > tolerance);
            }
            //small = false; // MIRO: FIX THIS!
            if (!small){
                std::fill(&(map[i*num_dimensions]), &(map[i*num_dimensions]) + num_dimensions, 1);
            }
        }
    }else{
        IndexManipulator IM(num_dimensions);
        int *dagUp = IM.computeDAGupLocal(points, rule);

        int max_1D_parents = rule->getMaxNumParents();
        int max_parents = max_1D_parents * num_dimensions;

        SplitDirections split(points);

        #pragma omp parallel for
        for(int s=0; s<split.getNumJobs(); s++){
            int d = split.getJobDirection(s);
            int nump = split.getJobNumPoints(s);
            const int *pnts = split.getJobPoints(s);

            int *global_to_pnts = new int[num_points];
            int *levels = new int[nump];
            int max_level = 0;

            int active_outputs = (output == -1) ? num_outputs : 1;

            double *vals = new double[nump * active_outputs];

            for(int i=0; i<nump; i++){
                const double* v = values->getValues(pnts[i]);
                const int *p = points->getIndex(pnts[i]);
                if (output == -1){
                    std::copy(v, v + num_outputs, &(vals[i*num_outputs]));
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
                                    for(int k=0; k<active_outputs; k++){
                                        vals[i * active_outputs + k] -= basis_value * vals[global_to_pnts[branch] * active_outputs + k];
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
            for(int i=0; i<nump; i++){
                bool small = true;
                if (output == -1){
                    for(int k=0; k<num_outputs; k++){
                        if (small && ((fabs(surpluses[pnts[i]*num_outputs + k]) / norm[k]) > tolerance) && ((fabs(vals[i*num_outputs + k]) / norm[k]) > tolerance)) small = false;
                    }
                }else{
                    if (((fabs(surpluses[pnts[i]*num_outputs + output]) / norm[output]) > tolerance) && ((fabs(vals[i]) / norm[output]) > tolerance)) small = false;
                }
                map[pnts[i]*num_dimensions + d] = (small) ? 0 : 1;;
            }

            delete[] vals;
            delete[] global_to_pnts;
            delete[] levels;
        }

        delete[] dagUp;
    }

    delete[] norm;

    return map;
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
//    kid[direction] = rule->getKidLeft(point[direction]);
//    if ((kid[direction] != -1) && (exclude->getSlot(kid) == -1)){
//        destination->addIndex(kid);
//    }
//    kid[direction] = rule->getKidRight(point[direction]);
//    if ((kid[direction] != -1) && (exclude->getSlot(kid) == -1)){
//        destination->addIndex(kid);
//    }
    delete[] kid;
}

void GridLocalPolynomial::clearRefinement(){
    if (needed != 0){  delete needed;  needed = 0;  }
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
void GridLocalPolynomial::setSurplusRefinement(double tolerance, TypeRefinement criteria, int output){
    clearRefinement();

    int *map = buildUpdateMap(tolerance, criteria, output);

    bool useParents = (criteria == refine_fds) || (criteria == refine_parents_first);

    GranulatedIndexSet *refined = new GranulatedIndexSet(num_dimensions);

    int num_points = points->getNumIndexes();

    for(int i=0; i<num_points; i++){
        for(int j=0; j<num_dimensions; j++){
            if (map[i*num_dimensions+j] == 1){ // if this dimension needs to be refined
                if (!(useParents && addParent(points->getIndex(i), j, refined, points))){
                    addChild(points->getIndex(i), j, refined, points);
                }
            }
        }
    }

    if (refined->getNumIndexes() > 0){
        needed = new IndexSet(refined);
    }
    delete refined;

    delete[] map;
}
int GridLocalPolynomial::removePointsBySurplus(double tolerance, int output){
    int num_points = points->getNumIndexes();
    bool *map = new bool[num_points]; //std::fill(map, map + num_points, false);

    double *norm = getNormalization();

    //#pragma omp parallel for
    for(int i=0; i<num_points; i++){
        bool small = true;
        if (output == -1){
            for(int k=0; k<num_outputs; k++){
                if (small && ((fabs(surpluses[i*num_outputs + k]) / norm[k]) > tolerance)) small = false;
            }
        }else{
            small = !((fabs(surpluses[i*num_outputs + output]) / norm[output]) > tolerance);
        }
        map[i] = !small;
        //if (!small) cout << "keep " << i << endl;
    }

    int num_kept = 0; for(int i=0; i<num_points; i++) num_kept += (map[i]) ? 1 : 0;
    //cout << "num kept = " << num_kept << endl;
    if (num_kept == 0){
        delete[] norm;
        delete[] map;
        return 0;
    }

    // save a copy of the points and the values
    int *point_kept = new int[num_kept * num_dimensions];
    double *values_kept = new double[num_kept * num_outputs];

    num_kept = 0;
    for(int i=0; i<num_points; i++){
        if (map[i]){
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
    delete[] map;

    return points->getNumIndexes();
}

double* GridLocalPolynomial::evalHierarchicalFunctions(const double x[]) const{
    IndexSet *work = (points == 0) ? needed : points;
    int num_points = work->getNumIndexes();
    double *vals = new double[num_points];
    bool dummy;
    for(int i=0; i<num_points; i++){
        //vals[i] = evalBasisRaw(work->getIndex(i), x);
        vals[i] = evalBasisSupported(work->getIndex(i), x, dummy);
    }
    return vals;
}
void GridLocalPolynomial::setHierarchicalCoefficients(const double c[]){
    int num_ponits = points->getNumIndexes();
    if (surpluses != 0) delete[] surpluses;
    surpluses = new double[num_ponits * num_outputs];
    std::copy(c, c + num_ponits * num_outputs, surpluses);
}

void GridLocalPolynomial::makeCheckAccelerationData(TypeAcceleration acc, std::ostream *os) const{
    #if defined(TASMANIAN_CUBLAS) || defined(TASMANIAN_CUDA)
    if (AccelerationMeta::isAccTypeFullMemoryGPU(acc)){
        if ((accel != 0) && (!accel->isCompatible(acc))){
            delete accel;
            accel = 0;
        }
        if (accel == 0){ accel = (BaseAccelerationData*) (new AccelerationDataGPUFull()); }
        AccelerationDataGPUFull *gpu = (AccelerationDataGPUFull*) accel;
        gpu->setLogStream(os);
        double *gpu_values = gpu->getGPUValues();
        if (gpu_values == 0) gpu->loadGPUValues(points->getNumIndexes() * values->getNumOutputs(), surpluses);
    }
    #endif // TASMANIAN_CUBLAS
}

void GridLocalPolynomial::clearAccelerationData(){
    if (accel != 0){
        delete accel;
        accel = 0;
    }
}

}

#endif
