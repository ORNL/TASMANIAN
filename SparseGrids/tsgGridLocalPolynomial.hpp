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

#include "tsgGridCore.hpp"

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

#ifndef __TASMANIAN_DOXYGEN_SKIP
class GridLocalPolynomial : public BaseCanonicalGrid{
public:
    GridLocalPolynomial(AccelerationContext const *acc) : BaseCanonicalGrid(acc), order(1), top_level(0){}
    friend struct GridReaderVersion5<GridLocalPolynomial>;
    GridLocalPolynomial(AccelerationContext const *acc, const GridLocalPolynomial *pwpoly, int ibegin, int iend);
    GridLocalPolynomial(AccelerationContext const *acc, int cnum_dimensions, int cnum_outputs, int depth, int corder, TypeOneDRule crule, const std::vector<int> &level_limits);
    ~GridLocalPolynomial() = default;

    bool isLocalPolynomial() const override{ return true; }

    void write(std::ostream &os, bool iomode) const override{ if (iomode == mode_ascii) write<mode_ascii>(os); else write<mode_binary>(os); }

    template<bool iomode> void write(std::ostream &os) const;

    TypeOneDRule getRule() const override{ return rule->getType(); }
    int getOrder() const{ return order; }

    void getLoadedPoints(double *x) const override;
    void getNeededPoints(double *x) const override;
    void getPoints(double *x) const override; // returns the loaded points unless no points are loaded, then returns the needed points

    void getQuadratureWeights(double weights[]) const override;
    void getInterpolationWeights(const double x[], double weights[]) const override;

    void loadNeededValues(const double *vals) override;

    void evaluate(const double x[], double y[]) const override;
    void integrate(double q[], double *conformal_correction) const override;

    void evaluateBatchOpenMP(const double x[], int num_x, double y[]) const;
    void evaluateBatch(const double x[], int num_x, double y[]) const override;

    void loadNeededValuesGPU(const double *vals);
    void evaluateGpuMixed(const double x[], int num_x, double y[]) const;
    void evaluateBatchGPU(const double gpu_x[], int cpu_num_x, double gpu_y[]) const override;
    void evaluateBatchGPU(const float gpu_x[], int cpu_num_x, float gpu_y[]) const override;
    template<typename T> void evaluateBatchGPUtempl(const T gpu_x[], int cpu_num_x, T gpu_y[]) const;
    void evaluateHierarchicalFunctionsGPU(const double gpu_x[], int cpu_num_x, double *gpu_y) const override;
    void buildSparseBasisMatrixGPU(const double gpu_x[], int cpu_num_x, GpuVector<int> &gpu_spntr, GpuVector<int> &gpu_sindx, GpuVector<double> &gpu_svals) const;
    void evaluateHierarchicalFunctionsGPU(const float gpu_x[], int cpu_num_x, float *gpu_y) const override;
    void buildSparseBasisMatrixGPU(const float gpu_x[], int cpu_num_x, GpuVector<int> &gpu_spntr, GpuVector<int> &gpu_sindx, GpuVector<float> &gpu_svals) const;

    void setSurplusRefinement(double tolerance, TypeRefinement criteria, int output, const std::vector<int> &level_limits, const double *scale_correction);
    void clearRefinement() override;
    void mergeRefinement() override;
    /*!
     * \brief Returns a vector with the normalized and rescaled coefficients.
     *
     * \param output indicates which output to use or can be -1 to indicate the use of all outputs (returns the max for each output)
     * \param scale_correction user provided corrections that rescale the coefficients by spacial or other importance
     */
    std::vector<double> getScaledCoefficients(int output, const double *scale_correction);
    int removePointsByHierarchicalCoefficient(double tolerance, int output, const double *scale_correction); // returns the number of points kept
    void removePointsByHierarchicalCoefficient(int new_num_points, int output, const double *scale_correction);
    /*!
     * \brief Remove all points marked as \b false in the \b pmap.
     *
     * \param pmap is a vector with size equal to the loaded points and flags as \b false each point that should be removed.
     */
    int removeMappedPoints(std::vector<bool> const &pmap);

    void beginConstruction() override;
    void writeConstructionData(std::ostream &os, bool) const override;
    void readConstructionData(std::istream &is, bool) override;
    std::vector<double> getCandidateConstructionPoints(double tolerance, TypeRefinement criteria, int output, std::vector<int> const &level_limits, double const *scale_correction);
    void loadConstructedPoint(const double x[], const std::vector<double> &y) override;
    void loadConstructedPoint(const double x[], int numx, const double y[]) override;
    void finishConstruction() override;

    void evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const override;
    std::vector<double> getSupport() const override final;
    void setHierarchicalCoefficients(const double c[]) override;
    void integrateHierarchicalFunctions(double integrals[]) const override;

    void updateAccelerationData(AccelerationContext::ChangeType change) const override;

    const double* getSurpluses() const;
    const int* getNeededIndexes() const;

    void buildSpareBasisMatrix(const double x[], int num_x, int num_chunk, std::vector<int> &spntr, std::vector<int> &sindx, std::vector<double> &svals) const;
    void buildSpareBasisMatrixStatic(const double x[], int num_x, int num_chunk, int *spntr, int *sindx, double *svals) const;
    int getSpareBasisMatrixNZ(const double x[], int num_x) const;

protected:
    //! \brief Create a new grid with given parameters and moving the data out of the vectors and sets.
    GridLocalPolynomial(AccelerationContext const *acc, int cnum_dimensions, int cnum_outputs, int corder, TypeOneDRule crule, std::vector<int> &&pnts, std::vector<double> &&vals, std::vector<double> &&surps);

    //! \brief Used as part of the loadNeededPoints() algorithm, updates the values and cuda cache, but does not touch the surpluses.
    void updateValues(double const *vals);

    void buildTree();

    //! \brief Returns a list of indexes of the nodes in \b points that are descendants of the \b point.
    std::vector<int> getSubGraph(std::vector<int> const &point) const;

    //! \brief Add the \b point to the grid using the \b values.
    void expandGrid(std::vector<int> const &point, std::vector<double> const &value);

    //! \brief Return the multi-index of canonical point \b x.
    std::vector<int> getMultiIndex(const double x[]);

    //! \brief Looks for a batch of constructed points and processes all that will result in a connected graph.
    void loadConstructedPoints();

    void recomputeSurpluses();

    /*!
     * \brief Update the surpluses for a portion of the graph.
     *
     * Updates the surpluses (i.e., hierarchical coefficients) for a portion of the DAG graph.
     * - \b work is the point set to consider
     * - \b max_level is the largest element in \b level
     * - \b level is the hierarchical level of the points in \b work (i.e., not the sum of multi-index entries);
     *   points where \b level is zero will not be computed
     * - \b dagUp must have been computed using \b MultiIndexManipulations::computeDAGup(\b work, \b rule, \b dagUp)
     *
     * Note: adjusting the \b level vector allows to update the surpluses for only a portion of the graph.
     */
    void updateSurpluses(MultiIndexSet const &work, int max_level, std::vector<int> const &level, Data2D<int> const &dagUp);

    void buildSparseMatrixBlockForm(const double x[], int num_x, int num_chunk, std::vector<int> &numnz,
                                    std::vector<std::vector<int>> &tindx, std::vector<std::vector<double>> &tvals) const;

    /*!
     * \brief Walk through all the nodes of the tree and touches only the nodes supported at \b x.
     *
     * The template is instantiated in different modes:
     * - \b mode \b 0, ignore \b sindx and \b svals, find the non-zero basis functions multiply them by the surpluses and add to \b y
     * - \b mode \b 1, ignore \b y, form a sparse vector by std::vector::push_back() to the \b sindx and \b svals
     * - \b mode \b 2, same as \b mode \b 1 but it also sorts the entries within the vector (requirement of Nvidia cusparseDgemvi)
     *
     * In all cases, \b work is the \b points or \b needed set that has been used to construct the tree.
     */
    template<int mode>
    void walkTree(const MultiIndexSet &work, const double x[], std::vector<int> &sindx, std::vector<double> &svals, double *y) const{
        std::vector<int> monkey_count(top_level+1); // traverse the tree, counts the branches of the current node
        std::vector<int> monkey_tail(top_level+1); // traverse the tree, keeps track of the previous node (history)

        for(const auto &r : roots){
            bool isSupported;
            double basis_value = evalBasisSupported(work.getIndex(r), x, isSupported);

            if (isSupported){
                if (mode == 0){
                    double const *s = surpluses.getStrip(r);
                    for(int k=0; k<num_outputs; k++) y[k] += basis_value * s[k];
                }else{
                    sindx.push_back(r);
                    svals.push_back(basis_value);
                }

                int current = 0;
                monkey_tail[0] = r;
                monkey_count[0] = pntr[r];

                while(monkey_count[0] < pntr[monkey_tail[0]+1]){
                    if (monkey_count[current] < pntr[monkey_tail[current]+1]){
                        int p = indx[monkey_count[current]];
                        basis_value = evalBasisSupported(work.getIndex(p), x, isSupported);
                        if (isSupported){
                            if (mode == 0){
                                double const *s = surpluses.getStrip(p);
                                for(int k=0; k<num_outputs; k++) y[k] += basis_value * s[k];
                            }else{
                                sindx.push_back(p);
                                svals.push_back(basis_value);
                            }

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
        if (mode == 2){
            std::vector<int> map(sindx);
            std::iota(map.begin(), map.end(), 0);
            std::sort(map.begin(), map.end(), [&](int a, int b)->bool{ return (sindx[a] < sindx[b]); });

            std::vector<int> idx = sindx;
            std::vector<double> vls = svals;
            std::transform(map.begin(), map.end(), sindx.begin(), [&](int i)->int{ return idx[i]; });
            std::transform(map.begin(), map.end(), svals.begin(), [&](int i)->double{ return vls[i]; });
        }
    }

    double evalBasisRaw(const int point[], const double x[]) const;
    double evalBasisSupported(const int point[], const double x[], bool &isSupported) const;

    std::vector<double> getNormalization() const;

    Data2D<int> buildUpdateMap(double tolerance, TypeRefinement criteria, int output, const double *scale_correction) const;
    MultiIndexSet getRefinementCanidates(double tolerance, TypeRefinement criteria, int output, const std::vector<int> &level_limits, const double *scale_correction) const;

    bool addParent(const int point[], int direction, const MultiIndexSet &exclude, Data2D<int> &destination) const;
    void addChild(const int point[], int direction, const MultiIndexSet &exclude, Data2D<int> &destination) const;
    void addChildLimited(const int point[], int direction, const MultiIndexSet &exclude, const std::vector<int> &level_limits, Data2D<int> &destination) const;

    void clearGpuSurpluses();
    void clearGpuBasisHierarchy();

private:
    int order, top_level;

    Data2D<double> surpluses;

    Data2D<int> parents;

    // tree for evaluation
    std::vector<int> roots;
    std::vector<int> pntr;
    std::vector<int> indx;

    std::unique_ptr<BaseRuleLocalPolynomial> rule;

    std::unique_ptr<SimpleConstructData> dynamic_values;

    // synchronize with tasgpu_devalpwpoly_feval
    template<int ord, TypeOneDRule crule, typename T>
    Data2D<T> encodeSupportForGPU(const MultiIndexSet &work) const{
        Data2D<T> cpu_support(num_dimensions, work.getNumIndexes());
        for(int i=0; i<work.getNumIndexes(); i++){
            const int* p = work.getIndex(i);
            T *s = cpu_support.getStrip(i);
            for(int j=0; j<num_dimensions; j++){
                s[j] = static_cast<T>(rule->getSupport(p[j]));
                if (ord != 0){
                    if (ord == 2) s[j] *= s[j];
                    if ((crule == rule_localp) || (crule == rule_semilocalp)) if (p[j] == 0) s[j] = static_cast<T>(-1.0); // constant function
                    if ((crule == rule_localp) && (ord == 2)){
                        if (p[j] == 1) s[j] = static_cast<T>(-2.0);
                        else if (p[j] == 2) s[j] = static_cast<T>(-3.0);
                    }
                    if ((crule == rule_semilocalp) && (ord == 2)){
                        if (p[j] == 1) s[j] = static_cast<T>(-4.0);
                        else if (p[j] == 2) s[j] = static_cast<T>(-5.0);
                    }
                    if ((crule == rule_localpb) && (ord == 2)){
                        if (p[j] < 2) s[j] = static_cast<T>(-2.0); // linear functions on level 0
                    }
                }
            }
        }
        return cpu_support;
    }
    std::unique_ptr<CudaLocalPolynomialData<double>>& getGpuCacheOverload(double) const{ return gpu_cache; }
    std::unique_ptr<CudaLocalPolynomialData<float>>& getGpuCacheOverload(float) const{ return gpu_cachef; }
    template<typename T> std::unique_ptr<CudaLocalPolynomialData<T>>& getGpuCache() const{
        return getGpuCacheOverload(static_cast<T>(0.0));
    }
    template<typename T> void loadGpuBasis() const;
    template<typename T> void loadGpuHierarchy() const;
    template<typename T> void loadGpuSurpluses() const;
    mutable std::unique_ptr<CudaLocalPolynomialData<double>> gpu_cache;
    mutable std::unique_ptr<CudaLocalPolynomialData<float>> gpu_cachef;
};

// Old version reader
template<> struct GridReaderVersion5<GridLocalPolynomial>{
    template<typename iomode> static std::unique_ptr<GridLocalPolynomial> read(AccelerationContext const *acc, std::istream &is){
        std::unique_ptr<GridLocalPolynomial> grid = Utils::make_unique<GridLocalPolynomial>(acc);

            grid->num_dimensions = IO::readNumber<iomode, int>(is);
            grid->num_outputs = IO::readNumber<iomode, int>(is);
            grid->order = IO::readNumber<iomode, int>(is);
            grid->top_level = IO::readNumber<iomode, int>(is);
            TypeOneDRule rule = IO::readRule<iomode>(is);
            grid->rule = makeRuleLocalPolynomial(rule, grid->order);

            if (IO::readFlag<iomode>(is)) grid->points = MultiIndexSet(is, iomode());
            if (std::is_same<iomode, IO::mode_ascii_type>::value){ // backwards compatible: surpluses and needed, or needed and surpluses
                if (IO::readFlag<iomode>(is))
                    grid->surpluses = IO::readData2D<iomode, double>(is, grid->num_outputs, grid->points.getNumIndexes());
                if (IO::readFlag<iomode>(is)) grid->needed = MultiIndexSet(is, iomode());
            }else{
                if (IO::readFlag<iomode>(is)) grid->needed = MultiIndexSet(is, iomode());
                if (IO::readFlag<iomode>(is))
                    grid->surpluses = IO::readData2D<iomode, double>(is, grid->num_outputs, grid->points.getNumIndexes());
            }
            if (IO::readFlag<iomode>(is))
                grid->parents = IO::readData2D<iomode, int>(is, grid->rule->getMaxNumParents() * grid->num_dimensions, grid->points.getNumIndexes());

            size_t num_points = (size_t) ((grid->points.empty()) ? grid->needed.getNumIndexes() : grid->points.getNumIndexes());
            grid->roots = std::vector<int>((size_t) IO::readNumber<iomode, int>(is));
            if (grid->roots.size() > 0){
                IO::readVector<iomode>(is, grid->roots);
                grid->pntr = IO::readVector<iomode, int>(is, num_points + 1);
                if (grid->pntr[num_points] > 0){
                    grid->indx = IO::readVector<iomode, int>(is, grid->pntr[num_points]);
                }else{
                    grid->indx = IO::readVector<iomode, int>(is, 1);
                }
            }

            if (grid->num_outputs > 0) grid->values = StorageSet(is, iomode());

        return grid;
    }
};
#endif // __TASMANIAN_DOXYGEN_SKIP
}

#endif
