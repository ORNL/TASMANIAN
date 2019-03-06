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

GridLocalPolynomial::GridLocalPolynomial() : num_dimensions(0), num_outputs(0), order(1), top_level(0), sparse_affinity(0)  {}
GridLocalPolynomial::~GridLocalPolynomial(){}

void GridLocalPolynomial::reset(bool clear_rule){
    clearAccelerationData();
    num_dimensions = num_outputs = top_level = 0;
    points = MultiIndexSet();
    needed = MultiIndexSet();
    values = StorageSet();
    if (clear_rule){ rule = std::unique_ptr<BaseRuleLocalPolynomial>(); order = 1; }
    parents.load(0, 0, 0);
    sparse_affinity = 0;
    surpluses.clear();
}
template<class T> std::unique_ptr<T> make_unique_ptr(){ return std::unique_ptr<T>(new T()); } // in C++14 this is called std::make_unique()
void GridLocalPolynomial::makeRule(TypeOneDRule crule){
    if (crule == rule_localp){
        rule = make_unique_ptr<templRuleLocalPolynomial<rule_localp, false>>();
    }else if (crule == rule_semilocalp){
        rule = make_unique_ptr<templRuleLocalPolynomial<rule_semilocalp, false>>();
    }else if (crule == rule_localp0){
        rule = make_unique_ptr<templRuleLocalPolynomial<rule_localp0, false>>();
    }else if (crule == rule_localpb){
        rule = make_unique_ptr<templRuleLocalPolynomial<rule_localpb, false>>();
    }
    if (order == 0){
        rule = make_unique_ptr<templRuleLocalPolynomial<rule_localp, true>>();
    }
    rule->setMaxOrder(order);
}

template<bool useAscii> void GridLocalPolynomial::write(std::ostream &os) const{
    if (useAscii){ os << std::scientific; os.precision(17); }
    IO::writeNumbers<useAscii, IO::pad_line>(os, num_dimensions, num_outputs, order, top_level);
    IO::writeRule<useAscii>(rule->getType(), os);
    IO::writeFlag<useAscii, IO::pad_auto>(!points.empty(), os);
    if (!points.empty()) points.write<useAscii>(os);
    if (useAscii){ // backwards compatible: surpluses and needed, or needed and surpluses
        IO::writeFlag<useAscii, IO::pad_auto>((surpluses.getNumStrips() != 0), os);
        if (surpluses.getNumStrips() != 0) IO::writeVector<useAscii, IO::pad_line>(surpluses.getVector(), os);
        IO::writeFlag<useAscii, IO::pad_auto>(!needed.empty(), os);
        if (!needed.empty()) needed.write<useAscii>(os);
    }else{
        IO::writeFlag<useAscii, IO::pad_auto>(!needed.empty(), os);
        if (!needed.empty()) needed.write<useAscii>(os);
        IO::writeFlag<useAscii, IO::pad_auto>((surpluses.getNumStrips() != 0), os);
        if (surpluses.getNumStrips() != 0) IO::writeVector<useAscii, IO::pad_line>(surpluses.getVector(), os);
    }
    IO::writeFlag<useAscii, IO::pad_auto>((parents.getNumStrips() != 0), os);
    if (parents.getNumStrips() != 0) IO::writeVector<useAscii, IO::pad_line>(parents.getVector(), os);

    IO::writeNumbers<useAscii, IO::pad_rspace>(os, (int) roots.size());
    IO::writeVector<useAscii, IO::pad_line>(roots, os);
    IO::writeVector<useAscii, IO::pad_line>(pntr, os);
    IO::writeVector<useAscii, IO::pad_line>(indx, os);

    if (num_outputs > 0) values.write<useAscii>(os);
}

template<bool useAscii> void GridLocalPolynomial::read(std::istream &is){
    reset();
    num_dimensions = IO::readNumber<useAscii, int>(is);
    num_outputs = IO::readNumber<useAscii, int>(is);
    order = IO::readNumber<useAscii, int>(is);
    top_level = IO::readNumber<useAscii, int>(is);
    TypeOneDRule crule = IO::readRule<useAscii>(is);
    makeRule(crule);

    if (IO::readFlag<useAscii>(is)) points.read<useAscii>(is);
    if (useAscii){ // backwards compatible: surpluses and needed, or needed and surpluses
        if (IO::readFlag<useAscii>(is)){
            surpluses.resize(num_outputs, points.getNumIndexes());
            IO::readVector<useAscii>(is, surpluses.getVector());
        }
        if (IO::readFlag<useAscii>(is)) needed.read<useAscii>(is);
    }else{
        if (IO::readFlag<useAscii>(is)) needed.read<useAscii>(is);
        if (IO::readFlag<useAscii>(is)){
            surpluses.resize(num_outputs, points.getNumIndexes());
            IO::readVector<useAscii>(is, surpluses.getVector());
        }
    }
    if (IO::readFlag<useAscii>(is)){
        parents.resize(rule->getMaxNumParents() * num_dimensions, points.getNumIndexes());
        IO::readVector<useAscii>(is, parents.getVector());
    }

    size_t num_points = (size_t) ((points.empty()) ? needed.getNumIndexes() : points.getNumIndexes());
    roots.resize((size_t) IO::readNumber<useAscii, int>(is));
    IO::readVector<useAscii>(is, roots);
    pntr.resize(num_points + 1);
    IO::readVector<useAscii>(is, pntr);
    if (pntr[num_points] > 0){
        indx.resize((size_t) pntr[num_points]);
        IO::readVector<useAscii>(is, indx);
    }else{
        indx.resize(1);
        indx[0] = IO::readNumber<useAscii, int>(is); // there is a special case when the grid has only one point without any children
    }

    if (num_outputs > 0) values.read<useAscii>(is);
}

template void GridLocalPolynomial::write<true>(std::ostream &) const;
template void GridLocalPolynomial::write<false>(std::ostream &) const;
template void GridLocalPolynomial::read<true>(std::istream &);
template void GridLocalPolynomial::read<false>(std::istream &);

void GridLocalPolynomial::makeGrid(int cnum_dimensions, int cnum_outputs, int depth, int corder, TypeOneDRule crule, const std::vector<int> &level_limits){
    reset();
    num_dimensions = cnum_dimensions;
    num_outputs = cnum_outputs;
    order = corder;

    TypeOneDRule effective_rule = ((crule == rule_semilocalp) && ((order == -1) || (order > 1))) ? rule_semilocalp : rule_localp; // semi-localp of order 1 is identical to localp
    if (crule == rule_localp0) effective_rule = rule_localp0;
    if (crule == rule_localpb) effective_rule = rule_localpb;

    makeRule(effective_rule);

    MultiIndexSet tensors;
    tensors.setNumDimensions(num_dimensions);
    MultiIndexManipulations::selectTensors(depth, type_level, [&](int i) -> long long{ return i; }, std::vector<int>(), tensors);

    if (!level_limits.empty()) MultiIndexManipulations::removeIndexesByLimit(level_limits, tensors);

    MultiIndexManipulations::generateNestedPoints(tensors, [&](int l) -> int{ return rule->getNumPoints(l); }, needed);

    buildTree();

    if (num_outputs == 0){
        points = std::move(needed);
        needed = MultiIndexSet();
        MultiIndexManipulations::computeDAGup(points, rule.get(), parents);
    }else{
        values.resize(num_outputs, needed.getNumIndexes());
    }
}

void GridLocalPolynomial::copyGrid(const GridLocalPolynomial *pwpoly){
    reset();
    num_dimensions = pwpoly->num_dimensions;
    num_outputs = pwpoly->num_outputs;
    order = pwpoly->order;

    makeRule(pwpoly->rule->getType());

    points = pwpoly->points;
    needed = pwpoly->needed;
    parents = pwpoly->parents;

    buildTree();

    values = pwpoly->values;

    if ((!points.empty()) && (num_outputs > 0)){ // points are loaded
        surpluses.resize(num_outputs, points.getNumIndexes());
        surpluses.getVector() = pwpoly->surpluses.getVector(); // copy assignment
    }
}

int GridLocalPolynomial::getNumDimensions() const{ return num_dimensions; }
int GridLocalPolynomial::getNumOutputs() const{ return num_outputs; }
TypeOneDRule GridLocalPolynomial::getRule() const{ return rule->getType(); }
int GridLocalPolynomial::getOrder() const{ return order; }

int GridLocalPolynomial::getNumLoaded() const{ return (num_outputs == 0) ? 0 : points.getNumIndexes(); }
int GridLocalPolynomial::getNumNeeded() const{ return needed.getNumIndexes(); }
int GridLocalPolynomial::getNumPoints() const{ return ((points.empty()) ? needed.getNumIndexes() : points.getNumIndexes()); }

void GridLocalPolynomial::getLoadedPoints(double *x) const{
    int num_points = points.getNumIndexes();
    Data2D<double> split;
    split.load(num_dimensions, num_points, x);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = points.getIndex(i);
        double *xx = split.getStrip(i);
        for(int j=0; j<num_dimensions; j++){
            xx[j] = rule->getNode(p[j]);
        }
    }
}
void GridLocalPolynomial::getNeededPoints(double *x) const{
    int num_points = needed.getNumIndexes();
    Data2D<double> split;
    split.load(num_dimensions, num_points, x);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = needed.getIndex(i);
        double *xx = split.getStrip(i);
        for(int j=0; j<num_dimensions; j++){
            xx[j] = rule->getNode(p[j]);
        }
    }
}
void GridLocalPolynomial::getPoints(double *x) const{
    if (points.empty()){ getNeededPoints(x); }else{ getLoadedPoints(x); }
}

void GridLocalPolynomial::evaluate(const double x[], double y[]) const{
    std::fill_n(y, num_outputs, 0.0);
    std::vector<int> sindx; // dummy variables, never references in mode 0 below
    std::vector<double> svals;
    walkTree<0>(points, x, sindx, svals, y);
}
void GridLocalPolynomial::evaluateBatch(const double x[], int num_x, double y[]) const{
    if (num_x == 1){ evaluate(x, y); return; }
    Data2D<double> xx; xx.cload(num_dimensions, num_x, x);
    Data2D<double> yy; yy.load(num_outputs, num_x, y);
    #pragma omp parallel for
    for(int i=0; i<num_x; i++)
        evaluate(xx.getCStrip(i), yy.getStrip(i));
}

#ifdef Tasmanian_ENABLE_BLAS
void GridLocalPolynomial::evaluateBlas(const double x[], int num_x, double y[]) const{
    if ((sparse_affinity == 1) || ((sparse_affinity == 0) && (num_outputs <= TSG_LOCALP_BLAS_NUM_OUTPUTS))){
        evaluateBatch(x, num_x, y);
        return;
    }

    std::vector<int> sindx, spntr;
    std::vector<double> svals;
    buildSpareBasisMatrix(x, num_x, 32, spntr, sindx, svals); // build sparse matrix corresponding to x

    int num_points = points.getNumIndexes();
    double nnz = (double) spntr[num_x];
    double total_size = ((double) num_x) * ((double) num_points);

    if ((sparse_affinity == -1) || ((sparse_affinity == 0) && (nnz / total_size > 0.1))){
        // potentially wastes a lot of memory
        Data2D<double> A;
        A.resize(num_points, num_x, 0.0);
        for(int i=0; i<num_x; i++){
            double *row = A.getStrip(i);
            for(int j=spntr[i]; j<spntr[i+1]; j++) row[sindx[j]] = svals[j];
        }
        TasBLAS::denseMultiply(num_outputs, num_x, num_points, 1.0, surpluses.getCStrip(0), A.getCStrip(0), 0.0, y);
    }else{
        Data2D<double> yy; yy.load(num_outputs, num_x, y);
        #pragma omp parallel for
        for(int i=0; i<num_x; i++){
            double *this_y = yy.getStrip(i);
            std::fill(this_y, this_y + num_outputs, 0.0);
            for(int j=spntr[i]; j<spntr[i+1]; j++){
                double v = svals[j];
                const double *s = surpluses.getCStrip(sindx[j]);
                for(int k=0; k<num_outputs; k++) this_y[k] += v * s[k];
            }
        }
    }
}
#endif

#ifdef Tasmanian_ENABLE_CUDA
void GridLocalPolynomial::evaluateCudaMixed(CudaEngine *engine, const double x[], int num_x, double y[]) const{
    loadCudaSurpluses();

    std::vector<int> sindx, spntr;
    std::vector<double> svals;

    if (num_x > 1){
        buildSpareBasisMatrix(x, num_x, 32, spntr, sindx, svals);
    }else{
        walkTree<2>(points, x, sindx, svals, nullptr);
    }
    engine->sparseMultiply(num_outputs, num_x, points.getNumIndexes(), 1.0, cuda_cache->surpluses, spntr, sindx, svals, y);
}
void GridLocalPolynomial::evaluateCuda(CudaEngine *engine, const double x[], int num_x, double y[]) const{
    if ((order == -1) || (order > 2) || (num_x == 1)){
        // GPU evaluations are available only for order 0, 1, and 2. Cubic will come later, but higher order will not be supported.
        // cannot use GPU to accelerate the evaluation of a single vector
        evaluateCudaMixed(engine, x, num_x, y);
        return;
    }
    loadCudaSurpluses();
    int num_points = points.getNumIndexes();

    bool useDense = (sparse_affinity == -1) || ((sparse_affinity == 0) && (num_dimensions > 6)); // dimension is the real criteria here

    CudaVector<double> gpu_x;
    gpu_x.load(((size_t) num_dimensions) * ((size_t) num_x), x);

    CudaVector<double> gpu_result(num_x, num_outputs);
    if (useDense){
        CudaVector<double> gpu_basis(num_x, num_points);
        buildDenseBasisMatrixGPU(gpu_x.data(), num_x, gpu_basis);

        engine->denseMultiply(num_outputs, num_x, num_points, 1.0, cuda_cache->surpluses, gpu_basis, 0.0, gpu_result);
    }else{
        CudaVector<int> gpu_spntr, gpu_sindx;
        CudaVector<double> gpu_svals;
        buildSparseBasisMatrixGPU(gpu_x.data(), num_x, gpu_spntr, gpu_sindx, gpu_svals);

        engine->sparseMultiply(num_outputs, num_x, num_points, 1.0, cuda_cache->surpluses, gpu_spntr, gpu_sindx, gpu_svals, 0.0, gpu_result);
    }
    gpu_result.unload(y);
}
#endif

void GridLocalPolynomial::loadNeededPoints(const double *vals, TypeAcceleration){
    #ifdef Tasmanian_ENABLE_CUDA
    clearCudaSurpluses();
    #endif
    if (needed.empty()){
        values.setValues(vals);
    }else{
        #ifdef Tasmanian_ENABLE_CUDA
        clearCudaBasisHierarchy();
        #endif
        if (points.empty()){ // initial grid, just relabel needed as points (loaded)
            values.setValues(vals);
            points = std::move(needed);
            needed = MultiIndexSet();
        }else{ // merge needed and points
            values.addValues(points, needed, vals);
            points.addSortedInsexes(needed.getVector());
            needed = MultiIndexSet();
            buildTree();
        }
    }
    recomputeSurpluses();
}
void GridLocalPolynomial::mergeRefinement(){
    if (needed.empty()) return; // nothing to do
    #ifdef Tasmanian_ENABLE_CUDA
    clearCudaSurpluses();
    #endif
    int num_all_points = getNumLoaded() + getNumNeeded();
    size_t num_vals = ((size_t) num_all_points) * ((size_t) num_outputs);
    std::vector<double> vals(num_vals, 0.0);
    values.setValues(vals);
    if (points.empty()){
        points = std::move(needed);
        needed = MultiIndexSet();
    }else{
        points.addMultiIndexSet(needed);
        needed = MultiIndexSet();
        buildTree();
    }
    surpluses.resize(num_outputs, num_all_points);
    surpluses.fill(0.0);
}

void GridLocalPolynomial::beginConstruction(){}
void GridLocalPolynomial::writeConstructionDataBinary(std::ofstream &ofs) const{}
void GridLocalPolynomial::writeConstructionData(std::ofstream &ofs) const{}
void GridLocalPolynomial::readConstructionDataBinary(std::ifstream &ifs){}
void GridLocalPolynomial::readConstructionData(std::ifstream &ifs){}
void GridLocalPolynomial::getCandidateConstructionPoints(double tolerance, TypeRefinement criteria, int output,
                                                         std::vector<int> const &level_limits, double const *scale_correction, std::vector<double> &x){}
void GridLocalPolynomial::loadConstructedPoint(const double x[], const std::vector<double> &y){}
void GridLocalPolynomial::finishConstruction(){}

void GridLocalPolynomial::getInterpolationWeights(const double x[], double *weights) const{
    const MultiIndexSet &work = (points.empty()) ? needed : points;

    std::vector<int> active_points;
    std::vector<double> hbasis_values;
    std::fill_n(weights, work.getNumIndexes(), 0.0);

    // construct a sparse vector and apply transpose surplus transformation
    walkTree<1>(work, x, active_points, hbasis_values, nullptr);
    auto ibasis = hbasis_values.begin();
    for(auto i : active_points) weights[i] = *ibasis++;

    // apply the transpose of the surplus transformation
    Data2D<int> lparents;
    if (parents.getNumStrips() != work.getNumIndexes()) // if the current dag loaded in parents does not reflect the indexes in work
        MultiIndexManipulations::computeDAGup(work, rule.get(), lparents);

    const Data2D<int> &dagUp = (parents.getNumStrips() != work.getNumIndexes()) ? lparents : parents;

    std::vector<int> level(active_points.size());
    int active_top_level = 0;
    for(size_t i=0; i<active_points.size(); i++){
        const int *p = work.getIndex(active_points[i]);
        int current_level = rule->getLevel(p[0]);
        for(int j=1; j<num_dimensions; j++){
            current_level += rule->getLevel(p[j]);
        }
        if (active_top_level < current_level) active_top_level = current_level;
        level[i] = current_level;
    }

    std::vector<int> monkey_count(top_level+1);
    std::vector<int> monkey_tail(top_level+1);
    std::vector<bool> used(work.getNumIndexes());
    std::vector<double> node(num_dimensions);
    int max_parents = rule->getMaxNumParents() * num_dimensions;

    for(int l=active_top_level; l>0; l--){
        for(size_t i=0; i<active_points.size(); i++){
            if (level[i] == l){
                const int* p = work.getIndex(active_points[i]);
                for(int j=0; j<num_dimensions; j++) node[j] = rule->getNode(p[j]);

                std::fill(used.begin(), used.end(), false);

                monkey_count[0] = 0;
                monkey_tail[0] = active_points[i];
                int current = 0;

                while(monkey_count[0] < max_parents){
                    if (monkey_count[current] < max_parents){
                        int branch = dagUp.getCStrip(monkey_tail[current])[monkey_count[current]];
                        if ((branch == -1) || used[branch]){
                            monkey_count[current]++;
                        }else{
                            const int *func = work.getIndex(branch);
                            double basis_value = rule->evalRaw(func[0], node[0]);
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
}

void GridLocalPolynomial::evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const{
    const MultiIndexSet &work = (points.empty()) ? needed : points;
    int num_points = work.getNumIndexes();
    Data2D<double> yy; yy.load(num_points, num_x, y);
    Data2D<double> xx; xx.cload(num_dimensions, num_x, x);
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        const double *this_x = xx.getCStrip(i);
        double *this_y = yy.getStrip(i);
        bool dummy;
        for(int j=0; j<num_points; j++){
            this_y[j] = evalBasisSupported(work.getIndex(j), this_x, dummy);
        }
    }
}

void GridLocalPolynomial::recomputeSurpluses(){
    int num_points = points.getNumIndexes();

    surpluses.resize(num_outputs, num_points);
    surpluses.getVector() = values.aliasValues(); // copy assignment

    Data2D<int> dagUp;
    MultiIndexManipulations::computeDAGup(points, rule.get(), dagUp);

    std::vector<int> level = MultiIndexManipulations::computeLevels(points, rule.get());

    updateSurpluses(points, top_level, level, dagUp);
}

void GridLocalPolynomial::updateSurpluses(MultiIndexSet const &work, int max_level, std::vector<int> const &level, Data2D<int> const &dagUp){
    int num_points = work.getNumIndexes();
    int max_parents = num_dimensions * rule->getMaxNumParents();
    for(int l=1; l<=max_level; l++){
        #pragma omp parallel for schedule(dynamic)
        for(int i=0; i<num_points; i++){
            if (level[i] == l){
                int const *p = work.getIndex(i);
                std::vector<double> x(num_dimensions);
                std::transform(p, p + num_dimensions, x.begin(), [&](int k)->double{ return rule->getNode(k); });
                double *surpi = surpluses.getStrip(i);

                std::vector<int> monkey_count(max_level + 1);
                std::vector<int> monkey_tail(max_level + 1);
                std::vector<bool> used(num_points, false);

                int current = 0;

                monkey_count[0] = 0;
                monkey_tail[0] = i;

                while(monkey_count[0] < max_parents){
                    if (monkey_count[current] < max_parents){
                        int branch = dagUp.getCStrip(monkey_tail[current])[monkey_count[current]];
                        if ((branch == -1) || (used[branch])){
                            monkey_count[current]++;
                        }else{
                            const double *branch_surp = surpluses.getCStrip(branch);
                            double basis_value = evalBasisRaw(work.getIndex(branch), x.data());
                            for(int k=0; k<num_outputs; k++)
                                surpi[k] -= basis_value * branch_surp[k];
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
    std::vector<std::vector<int>> tindx;
    std::vector<std::vector<double>> tvals;
    std::vector<int> numnz;
    buildSparseMatrixBlockForm(x, num_x, num_chunk, numnz, tindx, tvals);

    spntr = new int[num_x + 1];

    int nz = 0;
    for(int i=0; i<num_x; i++){
        spntr[i] = nz;
        nz += numnz[i];
    }
    spntr[num_x] = nz;

    sindx = new int[nz];
    svals = new double[nz];

    int c = 0;
    for(auto &idx : tindx) for(auto i: idx) sindx[c++] = i;
    c = 0;
    for(auto &vls : tvals) for(auto v: vls) svals[c++] = v;
}
void GridLocalPolynomial::buildSpareBasisMatrix(const double x[], int num_x, int num_chunk, std::vector<int> &spntr, std::vector<int> &sindx, std::vector<double> &svals) const{
    std::vector<std::vector<int>> tindx;
    std::vector<std::vector<double>> tvals;
    std::vector<int> numnz;
    buildSparseMatrixBlockForm(x, num_x, num_chunk, numnz, tindx, tvals);

    spntr.resize(num_x + 1);

    int nz = 0;
    auto inz = numnz.begin();
    for(auto &s: spntr){
        s = nz;
        nz += *inz++;
    }

    sindx.resize(nz);
    svals.resize(nz);

    auto ii = sindx.begin();
    for(auto &idx : tindx) for(auto i: idx) *ii++ = i;
    auto iv = svals.begin();
    for(auto &vls : tvals) for(auto v: vls) *iv++ = v;
}
void GridLocalPolynomial::buildSpareBasisMatrixStatic(const double x[], int num_x, int num_chunk, int *spntr, int *sindx, double *svals) const{
    std::vector<std::vector<int>> tindx;
    std::vector<std::vector<double>> tvals;
    std::vector<int> numnz;
    buildSparseMatrixBlockForm(x, num_x, num_chunk, numnz, tindx, tvals);

    int nz = 0;
    for(int i=0; i<num_x; i++){
        spntr[i] = nz;
        nz += numnz[i];
    }
    spntr[num_x] = nz;
    int c = 0;
    for(auto &idx : tindx) for(auto i: idx) sindx[c++] = i;
    c = 0;
    for(auto &vls : tvals) for(auto v: vls) svals[c++] = v;
}
int GridLocalPolynomial::getSpareBasisMatrixNZ(const double x[], int num_x) const{
    const MultiIndexSet &work = (points.empty()) ? needed : points;

    Data2D<double> xx; xx.cload(num_dimensions, num_x, x);
    std::vector<int> num_nz(num_x);

    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        std::vector<int> sindx;
        std::vector<double> svals;
        walkTree<1>(work, xx.getCStrip(i), sindx, svals, nullptr);
        num_nz[i] = (int) sindx.size();
    }

    return std::accumulate(num_nz.begin(), num_nz.end(), 0);
}

void GridLocalPolynomial::buildSparseMatrixBlockForm(const double x[], int num_x, int num_chunk, std::vector<int> &numnz, std::vector<std::vector<int>> &tindx, std::vector<std::vector<double>> &tvals) const{
    // numnz will be resized to (num_x + 1) with the last entry set to a dummy zero
    numnz.resize(num_x + 1);
    int num_blocks = num_x / num_chunk + ((num_x % num_chunk != 0) ? 1 : 0);

    tindx.resize(num_blocks);
    tvals.resize(num_blocks);

    const MultiIndexSet &work = (points.empty()) ? needed : points;
    Data2D<double> xx; xx.cload(num_dimensions, num_x, x);

    #pragma omp parallel for
    for(int b=0; b<num_blocks; b++){
        tindx[b].resize(0);
        tvals[b].resize(0);
        int chunk_size = (b < num_blocks - 1) ? num_chunk : (num_x - (num_blocks - 1) * num_chunk);
        for(int i = b * num_chunk; i < b * num_chunk + chunk_size; i++){
            numnz[i] = (int) tindx[b].size();
            walkTree<1>(work, xx.getCStrip(i), tindx[b], tvals[b], nullptr);
            numnz[i] = (int) tindx[b].size() - numnz[i];
        }
    }
    numnz[num_x] = 0;
}

#ifdef Tasmanian_ENABLE_CUDA
void GridLocalPolynomial::buildDenseBasisMatrixGPU(const double gpu_x[], int cpu_num_x, CudaVector<double> &gpu_y) const{
    loadCudaBasis();
    int num_points = getNumPoints();
    gpu_y.resize(((size_t) cpu_num_x) * ((size_t) num_points));
    TasCUDA::devalpwpoly(order, rule->getType(), num_dimensions, cpu_num_x, num_points, gpu_x, cuda_cache->nodes.data(), cuda_cache->support.data(), gpu_y.data());
}
void GridLocalPolynomial::buildSparseBasisMatrixGPU(const double gpu_x[], int cpu_num_x, CudaVector<int> &gpu_spntr, CudaVector<int> &gpu_sindx, CudaVector<double> &gpu_svals) const{
    loadCudaBasis();
    loadCudaHierarchy();
    TasCUDA::devalpwpoly_sparse(order, rule->getType(), num_dimensions, cpu_num_x, getNumPoints(), gpu_x, cuda_cache->nodes, cuda_cache->support,
                                cuda_cache->hpntr, cuda_cache->hindx, cuda_cache->hroots, gpu_spntr, gpu_sindx, gpu_svals);
}
#endif // Tasmanian_ENABLE_CUDA

void GridLocalPolynomial::buildTree(){
    const MultiIndexSet &work = (points.empty()) ? needed : points;
    int num_points = work.getNumIndexes();

    std::vector<int> level = MultiIndexManipulations::computeLevels(work, rule.get());

    top_level = *std::max_element(level.begin(), level.end());

    int max_1d_kids = rule->getMaxNumKids();
    int max_kids = max_1d_kids*num_dimensions;
    std::vector<int> monkey_count(top_level + 1);
    std::vector<int> monkey_tail(top_level + 1);

    std::vector<int>  tree(max_kids * num_points, -1);
    std::vector<bool> free(num_points, true);
    std::vector<int>  kid(num_dimensions);

    int next_root = 0;
    roots.resize(0);

    while(next_root != -1){
        roots.push_back(next_root);
        free[next_root] = false;

        monkey_tail[0] = next_root;
        monkey_count[0] = 0;
        int current = 0;

        while(monkey_count[0] < max_kids){
            if (monkey_count[current] < max_kids){
                const int *p = work.getIndex(monkey_tail[current]);

                int dir = monkey_count[current] / max_1d_kids;
                int ikid = rule->getKid(p[dir], monkey_count[current] % max_1d_kids);

                if (ikid == -1){
                    monkey_count[current]++;
                }else{
                    std::copy(p, p + num_dimensions, kid.data());
                    kid[dir] = ikid;
                    int t = work.getSlot(kid.data());
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
        int next_level = top_level + 1;
        for(int i=0; i<num_points; i++){
            if (free[i] && (level[i] < next_level)){
                next_root = i;
                next_level = level[i];
            }
        }
    }

    pntr.resize(num_points + 1);
    pntr[0] = 0;
    for(int i=0; i<num_points; i++){
        pntr[i+1] = pntr[i];
        for(int j=0; j<max_kids; j++) if (tree[i * max_kids + j] > -1) pntr[i+1]++;
    }

    indx.resize((pntr[num_points] > 0) ? pntr[num_points] : 1);
    indx[0] = 0;
    int count = 0;
    for(int i=0; i<num_points; i++){
        for(int j=0; j<max_kids; j++){
            int t = tree[i * max_kids + j];
            if (t > -1) indx[count++] = t;
        }
    }
}

void GridLocalPolynomial::getBasisIntegrals(double *integrals) const{
    const MultiIndexSet &work = (points.empty()) ? needed : points;

    int n = 0;
    std::vector<double> w, x;

    if ((rule->getMaxOrder() == -1) || (rule->getMaxOrder() > 3) ){
        n = top_level / 2 + 1;
        OneDimensionalNodes::getGaussLegendre(n, w, x);
    }

    for(int i=0; i<work.getNumIndexes(); i++){
        const int* p = work.getIndex(i);
        integrals[i] = rule->getArea(p[0], n, w.data(), x.data());
        for(int j=1; j<num_dimensions; j++){
            integrals[i] *= rule->getArea(p[j], n, w.data(), x.data());
        }
    }
}
void GridLocalPolynomial::getQuadratureWeights(double *weights) const{
    const MultiIndexSet &work = (points.empty()) ? needed : points;
    getBasisIntegrals(weights);

    std::vector<int> monkey_count(top_level+1);
    std::vector<int> monkey_tail(top_level+1);

    double basis_value;

    Data2D<int> lparents;
    if (parents.getNumStrips() != work.getNumIndexes())
        MultiIndexManipulations::computeDAGup(work, rule.get(), lparents);

    const Data2D<int> &dagUp = (parents.getNumStrips() != work.getNumIndexes()) ? lparents : parents;

    int num_points = work.getNumIndexes();
    std::vector<int> level = MultiIndexManipulations::computeLevels(work, rule.get());

    std::vector<double> node(num_dimensions);
    int max_parents = rule->getMaxNumParents() * num_dimensions;

    for(int l=top_level; l>0; l--){
        for(int i=0; i<num_points; i++){
            if (level[i] == l){
                const int* p = work.getIndex(i);
                for(int j=0; j<num_dimensions; j++) node[j] = rule->getNode(p[j]);

                std::vector<bool> used(work.getNumIndexes(), false);

                monkey_count[0] = 0;
                monkey_tail[0] = i;
                int current = 0;

                while(monkey_count[0] < max_parents){
                    if (monkey_count[current] < max_parents){
                        int branch = dagUp.getCStrip(monkey_tail[current])[monkey_count[current]];
                        if ((branch == -1) || used[branch]){
                            monkey_count[current]++;
                        }else{
                            const int *func = work.getIndex(branch);
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
}

void GridLocalPolynomial::integrate(double q[], double *conformal_correction) const{
    int num_points = points.getNumIndexes();
    std::fill(q, q + num_outputs, 0.0);

    if (conformal_correction == 0){
        std::vector<double> integrals(num_points);
        getBasisIntegrals(integrals.data());
        for(int i=0; i<num_points; i++){
            const double *s = surpluses.getCStrip(i);
            double wi = integrals[i];
            for(int k=0; k<num_outputs; k++) q[k] += wi * s[k];
        }
    }else{
        std::vector<double> w(num_points);
        getQuadratureWeights(w.data());
        for(int i=0; i<num_points; i++){
            double wi = w[i] * conformal_correction[i];
            const double *vals = values.getValues(i);
            for(int k=0; k<num_outputs; k++) q[k] += wi * vals[k];
        }
    }
}

void GridLocalPolynomial::getNormalization(std::vector<double> &norms) const{
    norms.resize(num_outputs);
    std::fill(norms.begin(), norms.end(), 0.0);
    for(int i=0; i<points.getNumIndexes(); i++){
        const double *v = values.getValues(i);
        for(int j=0; j<num_outputs; j++){
            if (norms[j] < fabs(v[j])) norms[j] = fabs(v[j]);
        }
    }
}

void GridLocalPolynomial::buildUpdateMap(double tolerance, TypeRefinement criteria, int output, const double *scale_correction, Data2D<int> &map2) const{
    int num_points = points.getNumIndexes();
    map2.resize(num_dimensions, num_points);
    if (tolerance == 0.0){
        map2.fill(1); // if tolerance is 0, refine everything
        return;
    }else{
        map2.fill(0);
    }

    std::vector<double> norm;
    getNormalization(norm);

    int active_outputs = (output == -1) ? num_outputs : 1;
    Data2D<double> scale;
    if (scale_correction == 0){
        scale.resize(active_outputs, num_points, 1.0);
    }else{
        scale.cload(active_outputs, num_points, scale_correction);
    }

    if ((criteria == refine_classic) || (criteria == refine_parents_first)){
        #pragma omp parallel for
        for(int i=0; i<num_points; i++){
            bool small = true;
            const double *s = surpluses.getCStrip(i);
            const double *c = scale.getCStrip(i);
            if (output == -1){
                for(int k=0; k<num_outputs; k++) small = small && ((c[k] * fabs(s[k]) / norm[k]) <= tolerance);
            }else{
                small = ((c[0] * fabs(s[output]) / norm[output]) <= tolerance);
            }
            if (!small){
                int *m = map2.getStrip(i);
                std::fill(m, m + num_dimensions, 1);
            }
        }
    }else{
        // construct a series of 1D interpolants and use a refinement criteria that is a combination of the two hierarchical coefficients
        Data2D<int> dagUp;
        MultiIndexManipulations::computeDAGup(points, rule.get(), dagUp);

        int max_1D_parents = rule->getMaxNumParents();

        SplitDirections split(points);

        #pragma omp parallel for
        for(int j=0; j<split.getNumJobs(); j++){ // split.getNumJobs() gives the number of 1D interpolants to construct
            int d = split.getJobDirection(j);
            int nump = split.getJobNumPoints(j);
            const int *pnts = split.getJobPoints(j);

            std::vector<int> global_to_pnts(num_points);
            std::vector<int> levels(nump);

            int max_level = 0;

            Data2D<double> vals;
            vals.resize(active_outputs, nump);

            for(int i=0; i<nump; i++){
                const double* v = values.getValues(pnts[i]);
                const int *p = points.getIndex(pnts[i]);
                if (output == -1){
                    std::copy(v, v + num_outputs, vals.getStrip(i));
                }else{
                    vals.getStrip(i)[0] = v[output];
                }
                global_to_pnts[pnts[i]] = i;
                levels[i] = rule->getLevel(p[d]);
                if (max_level < levels[i]) max_level = levels[i];
            }

            std::vector<int> monkey_count(max_level + 1);
            std::vector<int> monkey_tail(max_level + 1);

            for(int l=1; l<=max_level; l++){
                for(int i=0; i<nump; i++){
                    if (levels[i] == l){
                        const int *p = points.getIndex(pnts[i]);
                        double x = rule->getNode(p[d]);
                        double *valsi = vals.getStrip(i);

                        int current = 0;
                        monkey_count[0] = d * max_1D_parents;
                        monkey_tail[0] = pnts[i]; // uses the global indexes
                        std::vector<bool> used(nump, false);

                        while(monkey_count[0] < (d+1) * max_1D_parents){
                            if (monkey_count[current] < (d+1) * max_1D_parents){
                                int branch = dagUp.getStrip(monkey_tail[current])[monkey_count[current]];
                                if ((branch == -1) || (used[global_to_pnts[branch]])){
                                    monkey_count[current]++;
                                }else{
                                    const int *branch_point = points.getIndex(branch);
                                    double basis_value = rule->evalRaw(branch_point[d], x);
                                    const double *branch_vals = vals.getCStrip(global_to_pnts[branch]);
                                    for(int k=0; k<active_outputs; k++) valsi[k] -= basis_value * branch_vals[k];

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

            // at this point, vals contains the one directional surpluses
            for(int i=0; i<nump; i++){
                const double *s = surpluses.getCStrip(pnts[i]);
                const double *c = scale.getCStrip(pnts[i]);
                const double *v = vals.getCStrip(i);
                bool small = true;
                if (output == -1){
                    for(int k=0; k<num_outputs; k++){
                        small = small && (((c[k] * fabs(s[k]) / norm[k]) <= tolerance) || ((c[k] * fabs(v[k]) / norm[k]) <= tolerance));
                    }
                }else{
                    small = ((c[0] * fabs(s[output]) / norm[output]) <= tolerance) || ((c[0] * fabs(v[0]) / norm[output]) <= tolerance);
                }
                map2.getStrip(pnts[i])[d] = (small) ? 0 : 1;;
            }
        }
    }
}
MultiIndexSet GridLocalPolynomial::getRefinementCanidates(double tolerance, TypeRefinement criteria, int output, const std::vector<int> &level_limits, const double *scale_correction) const{
    Data2D<int> pmap;
    buildUpdateMap(tolerance, criteria, output, scale_correction, pmap);

    bool useParents = (criteria == refine_fds) || (criteria == refine_parents_first);

    Data2D<int> refined;
    refined.resize(num_dimensions, 0);

    int num_points = points.getNumIndexes();

    if (level_limits.empty()){
        for(int i=0; i<num_points; i++){
            const int *map = pmap.getCStrip(i);
            for(int j=0; j<num_dimensions; j++){
                if (map[j] == 1){ // if this dimension needs to be refined
                    if (!(useParents && addParent(points.getIndex(i), j, points, refined))){
                        addChild(points.getIndex(i), j, points, refined);
                    }
                }
            }
        }
    }else{
        for(int i=0; i<num_points; i++){
            const int *map = pmap.getCStrip(i);
            for(int j=0; j<num_dimensions; j++){
                if (map[j] == 1){ // if this dimension needs to be refined
                    if (!(useParents && addParent(points.getIndex(i), j, points, refined))){
                        addChildLimited(points.getIndex(i), j, points, level_limits, refined);
                    }
                }
            }
        }
    }

    MultiIndexSet candidates;
    if (refined.getNumStrips() > 0){
        candidates.setNumDimensions(num_dimensions);
        candidates.addData2D(refined);
    }

    return candidates;
}

bool GridLocalPolynomial::addParent(const int point[], int direction, const MultiIndexSet &exclude, Data2D<int> &destination) const{
    std::vector<int> dad(num_dimensions);
    std::copy_n(point, num_dimensions, dad.data());
    bool added = false;
    dad[direction] = rule->getParent(point[direction]);
    if ((dad[direction] != -1) && exclude.missing(dad)){
        destination.appendStrip(dad);
        added = true;
    }
    dad[direction] = rule->getStepParent(point[direction]);
    if ((dad[direction] != -1) && exclude.missing(dad)){
        destination.appendStrip(dad);
        added = true;
    }
    return added;
}
void GridLocalPolynomial::addChild(const int point[], int direction, const MultiIndexSet &exclude, Data2D<int> &destination) const{
    std::vector<int> kid(num_dimensions);
    std::copy_n(point, num_dimensions, kid.data());
    int max_1d_kids = rule->getMaxNumKids();
    for(int i=0; i<max_1d_kids; i++){
        kid[direction] = rule->getKid(point[direction], i);
        if ((kid[direction] != -1) && exclude.missing(kid)){
            destination.appendStrip(kid);
        }
    }
}
void GridLocalPolynomial::addChildLimited(const int point[], int direction, const MultiIndexSet &exclude, const std::vector<int> &level_limits, Data2D<int> &destination) const{
    std::vector<int> kid(num_dimensions);
    std::copy_n(point, num_dimensions, kid.data());
    int max_1d_kids = rule->getMaxNumKids();
    for(int i=0; i<max_1d_kids; i++){
        kid[direction] = rule->getKid(point[direction], i);
        if ((kid[direction] != -1)
            && ((level_limits[direction] == -1) || (rule->getLevel(kid[direction]) <= level_limits[direction]))
            && exclude.missing(kid)){
            destination.appendStrip(kid);
        }
    }
}

void GridLocalPolynomial::clearRefinement(){ needed = MultiIndexSet(); }
const double* GridLocalPolynomial::getSurpluses() const{
    return surpluses.getVector().data();
}
const int* GridLocalPolynomial::getPointIndexes() const{
    return ((points.empty()) ? needed.getIndex(0) : points.getIndex(0));
}
const int* GridLocalPolynomial::getNeededIndexes() const{
    return (needed.empty()) ? 0 : needed.getIndex(0);
}
void GridLocalPolynomial::setSurplusRefinement(double tolerance, TypeRefinement criteria, int output, const std::vector<int> &level_limits, const double *scale_correction){
    clearRefinement();

    needed = getRefinementCanidates(tolerance, criteria, output, level_limits, scale_correction);
}
int GridLocalPolynomial::removePointsByHierarchicalCoefficient(double tolerance, int output, const double *scale_correction){
    clearRefinement();
    int num_points = points.getNumIndexes();
    std::vector<bool> pmap(num_points); // point map, set to true if the point is to be kept, false otherwise

    std::vector<double> norm;
    getNormalization(norm);

    Data2D<double> scale;
    int active_outputs = (output == -1) ? num_outputs : 1;
    if (scale_correction == 0){
        scale.resize(active_outputs, num_points, 1.0);
    }else{
        scale.cload(active_outputs, num_points, scale_correction);
    }

    for(int i=0; i<num_points; i++){
        bool small = true;
        const double *s = surpluses.getCStrip(i);
        const double *c = scale.getCStrip(i);
        if (output == -1){
            for(int k=0; k<num_outputs; k++) small = small && ((c[k] * fabs(s[k]) / norm[k]) <= tolerance);
        }else{
            small = ((c[0] * fabs(s[output]) / norm[output]) <= tolerance);
        }
        pmap[i] = !small;
    }

    int num_kept = 0; for(int i=0; i<num_points; i++) num_kept += (pmap[i]) ? 1 : 0;

    if (num_kept == num_points) return num_points; // trivial case, remove nothing

    // save a copy of the points and the values
    std::vector<int> point_kept(((size_t) num_kept) * ((size_t) num_dimensions));
    Data2D<int> pp;
    pp.load(num_dimensions, num_kept, point_kept.data()); // to handle double indexing

    StorageSet values_kept;
    values_kept.resize(num_outputs, num_kept);
    values_kept.aliasValues().resize(((size_t) num_kept) * ((size_t) num_outputs));

    num_kept = 0;
    for(int i=0; i<num_points; i++){
        if (pmap[i]){
            std::copy(points.getIndex(i), points.getIndex(i) + num_dimensions, pp.getStrip(num_kept));
            std::copy_n(values.getValues(i), num_outputs, values_kept.getValues(num_kept));
            num_kept++;
        }
    }

    int dims = num_dimensions, outs = num_outputs;

    reset(false);
    num_dimensions = dims; num_outputs = outs;
    if (num_kept == 0) return 0; // trivial case, remove all

    points = MultiIndexSet();
    points.setNumDimensions(num_dimensions);
    points.addSortedInsexes(point_kept);

    values = std::move(values_kept);

    buildTree();
    recomputeSurpluses();

    return points.getNumIndexes();
}

void GridLocalPolynomial::setHierarchicalCoefficients(const double c[], TypeAcceleration){
    #ifdef Tasmanian_ENABLE_CUDA
    clearCudaSurpluses();
    #endif
    if (points.empty()){
        points = std::move(needed);
        needed = MultiIndexSet();
    }else{
        clearRefinement();
    }
    surpluses.resize(num_outputs, getNumPoints());
    std::copy_n(c, surpluses.getTotalEntries(), surpluses.getVector().data());

    std::vector<double> &vals = values.aliasValues();
    vals.resize(surpluses.getTotalEntries());

    std::vector<double> x(((size_t) getNumPoints()) * ((size_t) num_dimensions));
    getPoints(x.data());
    evaluateBatch(x.data(), points.getNumIndexes(), vals.data());
}

void GridLocalPolynomial::clearAccelerationData(){
    #ifdef Tasmanian_ENABLE_CUDA
    cuda_cache.reset();
    #endif
}

void GridLocalPolynomial::setFavorSparse(bool favor){
    // sparse_affinity == -1: use dense algorithms
    // sparse_affinity ==  1: use sparse algorithms
    // sparse_affinity ==  0: let Tasmanian decide
    // favor true/false bumps you upper or lower on the scale
    if (favor && (sparse_affinity < 1)) sparse_affinity++;
    if (!favor && (sparse_affinity > -1)) sparse_affinity--;
}

}

#endif
