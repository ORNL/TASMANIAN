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

#ifndef __TASMANIAN_SPARSE_GRID_LPOLY_HPP
#define __TASMANIAN_SPARSE_GRID_LPOLY_HPP

#include <vector>
#include <memory>

#include "tsgIOHelpers.hpp"
#include "tsgEnumerates.hpp"
#include "tsgIndexSets.hpp"
#include "tsgIndexManipulator.hpp"
#include "tsgGridCore.hpp"
#include "tsgRuleLocalPolynomial.hpp"

#include "tsgAcceleratedDataStructures.hpp"

//! \internal
//! \file tsgGridLocalPolynomial.hpp
//! \brief Algorithms for manipulating sets of multi-indexes.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianLocalPolynomialGrids
//!
//! Contains the implementation of the Local Polynomial Grids

//! \internal
//! \defgroup TasmanianLocalPolynomialGrids Local Polynomial Grid
//!
//! \par Local Polynomial Grid
//! A grid constructed from a hierarchy of one dimensional points and
//! basis with compact support.

namespace TasGrid{

class GridLocalPolynomial : public BaseCanonicalGrid{
public:
    GridLocalPolynomial();
    ~GridLocalPolynomial();

    bool isLocalPolynomial() const{ return true; }

    template<bool useAscii> void write(std::ostream &os) const;
    template<bool useAscii> void read(std::istream &is);

    void makeGrid(int cnum_dimensions, int cnum_outputs, int depth, int corder, TypeOneDRule crule, const std::vector<int> &level_limits);
    void copyGrid(const GridLocalPolynomial *pwpoly);

    int getNumDimensions() const;
    int getNumOutputs() const;
    TypeOneDRule getRule() const;
    int getOrder() const;

    int getNumLoaded() const;
    int getNumNeeded() const;
    int getNumPoints() const;

    void getLoadedPoints(double *x) const;
    void getNeededPoints(double *x) const;
    void getPoints(double *x) const; // returns the loaded points unless no points are loaded, then returns the needed points

    void getQuadratureWeights(double weights[]) const;
    void getInterpolationWeights(const double x[], double weights[]) const;

    void loadNeededPoints(const double *vals, TypeAcceleration acc = accel_none);

    void evaluate(const double x[], double y[]) const;
    void integrate(double q[], double *conformal_correction) const;

    void evaluateBatch(const double x[], int num_x, double y[]) const;

    #ifdef Tasmanian_ENABLE_BLAS
    void evaluateFastCPUblas(const double x[], double y[]) const;
    void evaluateBatchCPUblas(const double x[], int num_x, double y[]) const;
    #endif

    #ifdef Tasmanian_ENABLE_CUDA
    void evaluateFastGPUcublas(const double x[], double y[]) const;
    void evaluateFastGPUcuda(const double x[], double y[]) const;
    void evaluateBatchGPUcublas(const double x[], int num_x, double y[]) const;
    void evaluateBatchGPUcuda(const double x[], int num_x, double y[]) const;
    #endif

    #ifdef Tasmanian_ENABLE_MAGMA
    void evaluateFastGPUmagma(int gpuID, const double x[], double y[]) const;
    void evaluateBatchGPUmagma(int gpuID, const double x[], int num_x, double y[]) const;
    #endif

    void setSurplusRefinement(double tolerance, TypeRefinement criteria, int output, const std::vector<int> &level_limits, const double *scale_correction);
    void clearRefinement();
    void mergeRefinement();
    int removePointsByHierarchicalCoefficient(double tolerance, int output, const double *scale_correction); // returns the number of points kept

    void evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const;
    void setHierarchicalCoefficients(const double c[], TypeAcceleration acc);

    void clearAccelerationData();
    void setFavorSparse(bool favor);

    const double* getSurpluses() const;
    const int* getPointIndexes() const;
    const int* getNeededIndexes() const;

    void buildSpareBasisMatrix(const double x[], int num_x, int num_chunk, int* &spntr, int* &sindx, double* &svals) const;
    void buildSpareBasisMatrix(const double x[], int num_x, int num_chunk, std::vector<int> &spntr, std::vector<int> &sindx, std::vector<double> &svals) const;
    void buildSpareBasisMatrixStatic(const double x[], int num_x, int num_chunk, int *spntr, int *sindx, double *svals) const;
    int getSpareBasisMatrixNZ(const double x[], int num_x) const;

    void buildDenseBasisMatrixGPU(const double gpu_x[], int cpu_num_x, double gpu_y[]) const;
    void buildSparseBasisMatrixGPU(const double gpu_x[], int cpu_num_x, int* &gpu_spntr, int* &gpu_sindx, double* &gpu_svals, int &num_nz) const;

    #ifdef Tasmanian_ENABLE_CUDA
    void buildDenseBasisMatrixGPU(const double gpu_x[], int cpu_num_x, cudaDoubles &gpu_y) const;
    void buildSparseBasisMatrixGPU(const double gpu_x[], int cpu_num_x, cudaInts &gpu_spntr, cudaInts &gpu_sindx, cudaDoubles &gpu_svals) const;
    #endif

protected:
    void reset(bool clear_rule = true);

    //! \internal
    //! \brief makes the unique pointer associated with this rule, assuming that **order** is already set
    //! \ingroup TasmanianLocalPolynomialGrids
    void makeRule(TypeOneDRule trule);

    void buildTree();

    void recomputeSurpluses();

    void buildSparseMatrixBlockForm(const double x[], int num_x, int num_chunk, std::vector<int> &numnz,
                                    std::vector<std::vector<int>> &tindx, std::vector<std::vector<double>> &tvals) const;

    template<int mode>
    void buildSparseVector(const MultiIndexSet &work, const double x[], int &num_nz, std::vector<int> &sindx, std::vector<double> &svals) const{
        // This operates under several modes, no need to make separate enumerates for internal API, just use integer codes
        // mode 0, count the non-zeros, sindx and svals will never be accessed can pass dummy empty vectors
        // mode 1, num_nz is input (already pre-computed), sindx and svals are resized and filled
        // mode 2, num_nz is output (counted and filled), sindx and svals are NOT resized, std::vector::push_back() is used
        std::vector<int> monkey_count(top_level+1);
        std::vector<int> monkey_tail(top_level+1);

        if (mode == 1){
            sindx.resize(num_nz);
            svals.resize(num_nz);
        }

        bool isSupported;
        int p;

        num_nz = 0;

        for(const auto &r : roots){
            double basis_value = evalBasisSupported(work.getIndex(r), x, isSupported);

            if (isSupported){
                if (mode == 1){
                    sindx[num_nz] = r;
                    svals[num_nz] = basis_value;
                }else if (mode == 2){
                    sindx.push_back(r);
                    svals.push_back(basis_value);
                }
                num_nz++;

                int current = 0;
                monkey_tail[0] = r;
                monkey_count[0] = pntr[r];

                while(monkey_count[0] < pntr[monkey_tail[0]+1]){
                    if (monkey_count[current] < pntr[monkey_tail[current]+1]){
                        p = indx[monkey_count[current]];
                        basis_value = evalBasisSupported(work.getIndex(p), x, isSupported);
                        if (isSupported){
                            if (mode == 1){
                                sindx[num_nz] = p;
                                svals[num_nz] = basis_value;
                            }else if (mode == 2){
                                sindx.push_back(p);
                                svals.push_back(basis_value);
                            }
                            num_nz++;

                            monkey_tail[++current] = p;
                            monkey_count[current] = pntr[p];
                        }else{
                            monkey_count[current]++;
                        }
                    }else{
                        monkey_count[--current]++;
                    }
                }
            }
        }

        // according to https://docs.nvidia.com/cuda/cusparse/index.html#sparse-format
        // "... it is assumed that the indices are provided in increasing order and that each index appears only once."
        // This may not be a requirement for cusparseDgemvi(), but it may be that I have not tested it sufficiently
        // Also, see AccelerationDataGPUFull::cusparseMatveci() for inaccuracies in Nvidia documentation
        if (mode == 1){
            bool isNotSorted = false;
            for(int i=0; i<num_nz-1; i++) if (sindx[i] > sindx[i+1]) isNotSorted = true;
            if (isNotSorted){ // sort the vector
                std::vector<int> idx1(num_nz);
                std::vector<int> idx2(num_nz);
                std::vector<double> vls1(num_nz);
                std::vector<double> vls2(num_nz);

                int loop_end = (num_nz % 2 == 1) ? num_nz - 1 : num_nz;
                for(int i=0; i<loop_end; i+=2){
                    if (sindx[i] < sindx[i+1]){
                        idx1[i] = sindx[i];  idx1[i+1] = sindx[i+1];
                        vls1[i] = svals[i];  vls1[i+1] = svals[i+1];
                    }else{
                        idx1[i] = sindx[i+1];  idx1[i+1] = sindx[i];
                        vls1[i] = svals[i+1];  vls1[i+1] = svals[i];
                    }
                }
                if (num_nz % 2 == 1){
                    idx1[num_nz - 1] = sindx[num_nz - 1];
                    vls1[num_nz - 1] = svals[num_nz - 1];
                }

                int stride = 2;
                while(stride < num_nz){
                    int c = 0;
                    while(c < num_nz){
                        int acurrent = c;
                        int bcurrent = c + stride;
                        int aend = (acurrent + stride < num_nz) ? acurrent + stride : num_nz;
                        int bend = (bcurrent + stride < num_nz) ? bcurrent + stride : num_nz;
                        while((acurrent < aend) || (bcurrent  < bend)){
                            bool picka;
                            if (bcurrent >= bend){
                                picka = true;
                            }else if (acurrent >= aend){
                                picka = false;
                            }else{
                                picka = (idx1[acurrent] < idx1[bcurrent]);
                            }
                            if (picka){
                                idx2[c] = idx1[acurrent];
                                vls2[c] = vls1[acurrent];
                                acurrent++;
                            }else{
                                idx2[c] = idx1[bcurrent];
                                vls2[c] = vls1[bcurrent];
                                bcurrent++;
                            }
                            c++;
                        }
                    }
                    stride *= 2;
                    std::swap(idx1, idx2);
                    std::swap(vls1, vls2);
                }
                sindx = idx1;
                svals = vls1;
            }
        }
    }

    double evalBasisRaw(const int point[], const double x[]) const;
    double evalBasisSupported(const int point[], const double x[], bool &isSupported) const;

    void getBasisIntegrals(double *integrals) const;

    void getNormalization(std::vector<double> &norms) const;

    void buildUpdateMap(double tolerance, TypeRefinement criteria, int output, const double *scale_correction, Data2D<int> &map2) const;

    bool addParent(const int point[], int direction, const MultiIndexSet &exclude, Data2D<int> &destination) const;
    void addChild(const int point[], int direction, const MultiIndexSet &exclude, Data2D<int> &destination) const;
    void addChildLimited(const int point[], int direction, const MultiIndexSet &exclude, const std::vector<int> &level_limits, Data2D<int> &destination) const;

    #ifdef Tasmanian_ENABLE_CUDA
    // synchronize with tasgpu_devalpwpoly_feval
    template<int order, TypeOneDRule crule>
    void encodeSupportForGPU(const MultiIndexSet &work, double *cpu_support) const{
        for(int i=0; i<work.getNumIndexes(); i++){
            const int* p = work.getIndex(i);
            for(int j=0; j<num_dimensions; j++){
                cpu_support[i*num_dimensions + j] = rule->getSupport(p[j]);
                if (order != 0){
                    if (order == 2) cpu_support[i*num_dimensions + j] *= cpu_support[i*num_dimensions + j];
                    if ((crule == rule_localp) || (crule == rule_semilocalp)) if (p[j] == 0) cpu_support[i*num_dimensions + j] = -1.0; // constant function
                    if ((crule == rule_localp) && (order == 2)){
                        if (p[j] == 1) cpu_support[i*num_dimensions + j] = -2.0;
                        else if (p[j] == 2) cpu_support[i*num_dimensions + j] = -3.0;
                    }
                    if ((crule == rule_semilocalp) && (order == 2)){
                        if (p[j] == 1) cpu_support[i*num_dimensions + j] = -4.0;
                        else if (p[j] == 2) cpu_support[i*num_dimensions + j] = -5.0;
                    }
                    if ((crule == rule_localpb) && (order == 2)){
                        if (p[j] < 2) cpu_support[i*num_dimensions + j] = -2.0; // linear functions on level 0
                    }
                }
            }
        }
    }
    // synchronize with tasgpu_devalpwpoly_feval
    template<int order, TypeOneDRule crule>
    void encodeSupportForGPU(const MultiIndexSet &work, Data2D<double> &cpu_support) const{
        cpu_support.resize(num_dimensions, work.getNumIndexes());
        for(int i=0; i<work.getNumIndexes(); i++){
            const int* p = work.getIndex(i);
            double *s = cpu_support.getStrip(i);
            for(int j=0; j<num_dimensions; j++){
                s[j] = rule->getSupport(p[j]);
                if (order != 0){
                    if (order == 2) s[j] *= s[j];
                    if ((crule == rule_localp) || (crule == rule_semilocalp)) if (p[j] == 0) s[j] = -1.0; // constant function
                    if ((crule == rule_localp) && (order == 2)){
                        if (p[j] == 1) s[j] = -2.0;
                        else if (p[j] == 2) s[j] = -3.0;
                    }
                    if ((crule == rule_semilocalp) && (order == 2)){
                        if (p[j] == 1) s[j] = -4.0;
                        else if (p[j] == 2) s[j] = -5.0;
                    }
                    if ((crule == rule_localpb) && (order == 2)){
                        if (p[j] < 2) s[j] = -2.0; // linear functions on level 0
                    }
                }
            }
        }
    }
    void loadCudaData() const{
        cuda_surpluses.load(surpluses.getVector());
        std::vector<double> cpu_nodes(((size_t) getNumPoints()) * ((size_t) num_dimensions));
        getPoints(cpu_nodes.data());
        cuda_nodes.load(cpu_nodes);
        Data2D<double> cpu_support;
        const MultiIndexSet &work = (points.empty()) ? needed : points;
        if (rule->getType() == rule_localp){
            switch(order){
            case 0: encodeSupportForGPU<0, rule_localp>(work, cpu_support); break;
            case 2: encodeSupportForGPU<2, rule_localp>(work, cpu_support); break;
            default:
                encodeSupportForGPU<1, rule_localp>(work, cpu_support);
            }
        }else if (rule->getType() == rule_semilocalp){
            encodeSupportForGPU<2, rule_semilocalp>(work, cpu_support);
        }else if (rule->getType() == rule_localpb){
            switch(order){
            case 2: encodeSupportForGPU<2, rule_localpb>(work, cpu_support); break;
            default:
                encodeSupportForGPU<1, rule_localpb>(work, cpu_support);
            }
        }else{
            switch(order){
            case 2: encodeSupportForGPU<2, rule_localp0>(work, cpu_support); break;
            default:
                encodeSupportForGPU<1, rule_localp0>(work, cpu_support);
            }
        }
        cuda_support.load(cpu_support.getVector());
        cuda_pntr.load(pntr);
        cuda_indx.load(indx);
        cuda_roots.load(roots);
    }
    void clearCudaLoadedData(){
        cuda_surpluses.clear();
        cuda_nodes.clear();
        cuda_support.clear();
        cuda_pntr.clear();
        cuda_indx.clear();
        cuda_roots.clear();
    }
    #endif

private:
    int num_dimensions, num_outputs, order, top_level;

    Data2D<double> surpluses;

    MultiIndexSet points;
    MultiIndexSet needed;

    StorageSet values;
    Data2D<int> parents;

    // tree for evaluation
    std::vector<int> roots;
    std::vector<int> pntr;
    std::vector<int> indx;

    std::unique_ptr<BaseRuleLocalPolynomial> rule;

    int sparse_affinity;

    #ifdef Tasmanian_ENABLE_CUDA
    mutable LinearAlgebraEngineGPU cuda_engine;
    mutable cudaDoubles cuda_surpluses, cuda_nodes, cuda_support;
    mutable cudaInts cuda_pntr, cuda_indx, cuda_roots;
    #endif
};

}

#endif
