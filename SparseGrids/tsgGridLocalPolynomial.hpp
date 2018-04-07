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

#include "tsgEnumerates.hpp"
#include "tsgIndexSets.hpp"
#include "tsgIndexManipulator.hpp"
#include "tsgGridCore.hpp"
#include "tsgRuleLocalPolynomial.hpp"

#include "tsgAcceleratedDataStructures.hpp"

namespace TasGrid{

class GridLocalPolynomial : public BaseCanonicalGrid{
public:
    GridLocalPolynomial();
    GridLocalPolynomial(const GridLocalPolynomial &pwpoly);
    ~GridLocalPolynomial();

    void write(std::ofstream &ofs) const;
    void read(std::ifstream &ifs);

    void writeBinary(std::ofstream &ofs) const;
    void readBinary(std::ifstream &ifs);

    void makeGrid(int cnum_dimensions, int cnum_outputs, int depth, int corder, TypeOneDRule crule, const int *level_limits = 0);
    void copyGrid(const GridLocalPolynomial *pwpoly);

    int getNumDimensions() const;
    int getNumOutputs() const;
    TypeOneDRule getRule() const;
    int getOrder() const;

    int getNumLoaded() const;
    int getNumNeeded() const;
    int getNumPoints() const;

    double* getLoadedPoints() const;
    void getLoadedPoints(double *x) const;
    double* getNeededPoints() const;
    void getNeededPoints(double *x) const;
    double* getPoints() const;
    void getPoints(double *x) const; // returns the loaded points unless no points are loaded, then returns the needed points

    double* getQuadratureWeights() const;
    void getQuadratureWeights(double weights[]) const;
    double* getInterpolationWeights(const double x[]) const;
    void getInterpolationWeights(const double x[], double weights[]) const;

    void loadNeededPoints(const double *vals, TypeAcceleration acc = accel_none);

    void evaluate(const double x[], double y[]) const;
    void integrate(double q[], double *conformal_correction) const;

    void evaluateFastCPUblas(const double x[], double y[]) const;
    void evaluateFastGPUcublas(const double x[], double y[], std::ostream *os) const;
    void evaluateFastGPUcuda(const double x[], double y[], std::ostream *os) const;
    void evaluateFastGPUmagma(const double x[], double y[], std::ostream *os) const;

    void evaluateBatch(const double x[], int num_x, double y[]) const;
    void evaluateBatchCPUblas(const double x[], int num_x, double y[]) const;
    void evaluateBatchGPUcublas(const double x[], int num_x, double y[], std::ostream *os) const;
    void evaluateBatchGPUcuda(const double x[], int num_x, double y[], std::ostream *os) const;
    void evaluateBatchGPUmagma(const double x[], int num_x, double y[], std::ostream *os) const;

    void setSurplusRefinement(double tolerance, TypeRefinement criteria, int output, const int *level_limits, const double *scale_correction);
    void clearRefinement();
    void mergeRefinement();
    int removePointsByHierarchicalCoefficient(double tolerance, int output, const double *scale_correction); // returns the number of points kept

    void evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const;
    void setHierarchicalCoefficients(const double c[], TypeAcceleration acc, std::ostream *os);

    void clearAccelerationData();

    const double* getSurpluses() const;
    const int* getPointIndexes() const;
    const int* getNeededIndexes() const;

    void buildSpareBasisMatrix(const double x[], int num_x, int num_chunk, int* &spntr, int* &sindx, double* &svals) const;
    void buildSpareBasisMatrixStatic(const double x[], int num_x, int num_chunk, int *spntr, int *sindx, double *svals) const;
    int getSpareBasisMatrixNZ(const double x[], int num_x, int num_chunk) const;

    // EXPERIMENTAL: GPU evaluateHierarchicalFunctionsGPU()
    void buildDenseBasisMatrixGPU(const double gpu_x[], int cpu_num_x, double gpu_y[], std::ostream *os) const;
    void buildSparseBasisMatrixGPU(const double gpu_x[], int cpu_num_x, int* &gpu_spntr, int* &gpu_sindx, double* &gpu_svals, int &num_nz, std::ostream *os) const;

    // EXPERIMENTAL: mostly for tuning and testing purposes, force certain backend behavior
    inline void setBackendFlavor(TypeLocalPolynomialBackendFlavor new_flavor){ backend_flavor = new_flavor; }
    inline void setForceSparse(){ force_sparse = true; }

protected:
    void reset(bool clear_rule = true);

    void buildTree();

    void recomputeSurpluses();
    void recomputeSurplusesGPUcublas();
    void recomputeSurplusesGPUcuda();

    void buildSparseMatrixBlockForm(const double x[], int num_x, int num_chunk, int &num_blocks, int &num_last, int &stripe_size,
                                    int* &stripes, int* &last_stripe_size, int** &tpntr, int*** &tindx, double*** &tvals) const;

    template<bool fill>
    void buildSparseVector(const double x[], int &num_nz, int *sindx, double *svals) const{
        int *monkey_count = new int[top_level+1];
        int *monkey_tail = new int[top_level+1];

        bool isSupported;
        size_t offset;

        num_nz = 0;

        for(int r=0; r<num_roots; r++){
            double basis_value = evalBasisSupported(points->getIndex(roots[r]), x, isSupported);

            if (isSupported){
                offset = roots[r] * num_outputs;
                if (fill){
                    sindx[num_nz] = roots[r];
                    svals[num_nz] = basis_value;
                }
                num_nz++;

                int current = 0;
                monkey_tail[0] = roots[r];
                monkey_count[0] = pntr[roots[r]];

                while(monkey_count[0] < pntr[monkey_tail[0]+1]){
                    if (monkey_count[current] < pntr[monkey_tail[current]+1]){
                        offset = indx[monkey_count[current]];
                        basis_value = evalBasisSupported(points->getIndex(offset), x, isSupported);
                        if (isSupported){
                            if (fill){
                                sindx[num_nz] = offset;
                                svals[num_nz] = basis_value;
                            }
                            num_nz++;

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

    double evalBasisRaw(const int point[], const double x[]) const;
    double evalBasisSupported(const int point[], const double x[], bool &isSupported) const;

    void getBasisIntegrals(double *integrals) const;

    double* getNormalization() const;

    int* buildUpdateMap(double tolerance, TypeRefinement criteria, int output, const double *scale_correction) const;

    bool addParent(const int point[], int direction, GranulatedIndexSet *destination, IndexSet *exclude) const;
    void addChild(const int point[], int direction, GranulatedIndexSet *destination, IndexSet *exclude) const;
    void addChildLimited(const int point[], int direction, GranulatedIndexSet *destination, IndexSet *exclude, const int *level_limits) const;

    void makeCheckAccelerationData(TypeAcceleration acc, std::ostream *os) const;
    void checkAccelerationGPUValues() const;
    void checkAccelerationGPUNodes() const;
    void checkAccelerationGPUHierarchy() const;

    // synchronize with tasgpu_devalpwpoly_feval
    template<int order, TypeOneDRule crule>
    void encodeSupportForGPU(const IndexSet *work, double *cpu_support) const{
        for(int i=0; i<work->getNumIndexes(); i++){
            const int* p = work->getIndex(i);
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
                }
            }
        }
    }

private:
    int num_dimensions, num_outputs, order, top_level;

    double *surpluses;

    IndexSet *points;
    IndexSet *needed;

    StorageSet *values;
    int *parents;

    // tree for evaluation
    int num_roots, *roots;
    int *pntr, *indx;

    BaseRuleLocalPolynomial *rule;

    templRuleLocalPolynomial<rule_localp, false> rpoly;
    templRuleLocalPolynomial<rule_semilocalp, false> rsemipoly;
    templRuleLocalPolynomial<rule_localp0, false> rpoly0;
    templRuleLocalPolynomial<rule_localp, true> rpolyc;

    mutable BaseAccelerationData *accel;
    TypeLocalPolynomialBackendFlavor backend_flavor;
    bool force_sparse;
};

}

#endif
