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
#ifndef __TASGRID_TEST_HELPERS_HPP
#define __TASGRID_TEST_HELPERS_HPP

#include <iostream>
#include "TasmanianSparseGrid.hpp"

using std::cout;
using std::endl;

using namespace TasGrid;

template<typename T> struct wrap_array{
    template<typename U> wrap_array(T *in_raw, U in_num) : raw(in_raw), num(size_t(in_num)){}
    T *raw;
    size_t num;
    size_t size() const{ return num; }
    T const& operator[](size_t i) const{ return raw[i]; }
};

template<typename VectorLike1, typename VectorLike2>
double err1(size_t num, VectorLike1 const &x, VectorLike2 const &y){
    if ((x.size() < num) || (y.size() < num)) throw std::runtime_error("vector size is insufficient");
    double err = 0.0;
    for(size_t i=0; i<num; i++) err = std::max(err, std::abs(x[i] - y[i]));
    return err;
}
template<typename VectorLike1, typename VectorLike2>
double err1(VectorLike1 const &x, VectorLike2 const &y){
    if (x.size() != y.size()) throw std::runtime_error("vector size mismatch");
    return err1(x.size(), x, y);
}

template<typename VectorLike>
std::vector<std::pair<int, typename VectorLike::value_type>>
extractRow(size_t i, std::vector<int> const &pntr, std::vector<int> const &indx, VectorLike const &vals){
    std::vector<std::pair<int, typename VectorLike::value_type>> result;
    result.reserve(pntr[i+1] - pntr[i]);
    for(int j=pntr[i]; j<pntr[i+1]; j++)
        result.push_back(std::pair<int, typename VectorLike::value_type>(indx[j], vals[j]));
    return result;
}

template<typename T>
void sortRow(std::vector<std::pair<int, T>> &row){
    std::sort(row.begin(), row.end(), [&](std::pair<int, T> &a, std::pair<int, T> &b)->
                                      bool{ return a.first < b.first; });
}

template<typename VectorLike1, typename VectorLike2>
double err1_sparse(std::vector<int> const &xpntr, std::vector<int> const &xindx, VectorLike1 const &xvals,
                   std::vector<int> const &ypntr, std::vector<int> const &yindx, VectorLike2 const &yvals){
    if (xpntr.size() != ypntr.size()) throw std::runtime_error("mismatch in number of sparse rows");
    double err = 0.0;
    for(size_t i=0; i<xpntr.size()-1; i++){
        auto xrow = extractRow(i, xpntr, xindx, xvals); sortRow(xrow);
        auto yrow = extractRow(i, ypntr, yindx, yvals); sortRow(yrow);
        auto ix = xrow.begin();
        auto iy = yrow.begin();
        while( ix < xrow.end() or iy < yrow.end() ){
            double xv = 0.0, yv = 0.0;
            if (ix == xrow.end()){ // read all x, using only y
                xv = 0.0;
                yv = iy++ ->second;
            }else if (iy == yrow.end()){ // read all y, using only x
                xv = ix++ ->second;
                yv = 0.0;
            }else if (ix->first == iy->first){ // index pair match
                xv = ix++ ->second;
                yv = iy++ ->second;
            }else if (ix->first < iy->first){ // y skips an index
                xv = ix++ ->second;
                yv = 0.0;
            }else{ // since sorted, x must have skipped an index
                xv = 0.0;
                yv = iy++ ->second;
            }
            err = std::max(err, std::abs(xv - yv));
        }
    }
    return err;
}

inline bool testPass(double err, double tol, std::string message, TasmanianSparseGrid const &grid = TasmanianSparseGrid()){
    if (err > tol){
        cout << std::scientific; cout.precision(16);
        cout << "ERROR: " << message << "\n  observed: " << err << "  expected: " << tol << "\n";
        if (!grid.empty()) grid.printStats();
        return false;
    }else{
        return true;
    }
}

inline std::vector<double> const& getVector(std::vector<double> const &x, std::vector<double>&, size_t){ return x; }
inline std::vector<float> const& getVector(std::vector<double> const &x, std::vector<float> &t, size_t num){
    t = std::vector<float>(num);
    std::transform(x.begin(), x.begin() + num, t.begin(), [](double const&v)->float{ return static_cast<float>(v); });
    return t;
}

struct GridMethodBatch{};
struct GridMethodFast{};
template<typename T, typename EvalMethod>
bool testAccEval(std::vector<double> const &x, std::vector<double> const &y, int numx, double tolerance, TasmanianSparseGrid &grid, std::string message){
    std::vector<T> test_y, temp_x;
    std::vector<T> const &test_x = getVector(x, temp_x, Utils::size_mult(grid.getNumDimensions(), numx));

    if (std::is_same<EvalMethod, GridMethodBatch>::value){
        grid.evaluateBatch(test_x, test_y);
    }else if (std::is_same<EvalMethod, GridMethodFast>::value){
        test_y = std::vector<T>(Utils::size_mult(numx, grid.getNumOutputs()));
        Utils::Wrapper2D<T const> wx(grid.getNumDimensions(), test_x.data());
        Utils::Wrapper2D<T> wy(grid.getNumOutputs(), test_y.data());
        for(int i=0; i<numx; i++)
            grid.evaluateFast(wx.getStrip(i), wy.getStrip(i));
    }
    return testPass(err1(Utils::size_mult(numx, grid.getNumOutputs()), test_y, y), tolerance, message, grid);
}

inline bool canUseCudaKernels(TasmanianSparseGrid const &grid){
    if (!(grid.getAccelerationType() == accel_gpu_cuda || grid.getAccelerationType() == accel_gpu_magma)) return false;
    if (grid.isLocalPolynomial() && (grid.getOrder() > 2 || grid.getOrder() == -1)) return false;
    if (grid.isWavelet() && grid.getOrder() == 3) return false;
    return true;
}

#ifdef Tasmanian_ENABLE_GPU
struct GridMethodHierBasisGPU{};
struct GridMethodEvalBatchGPU{};
template<typename T, typename GridMethod>
bool testDenseGPU(std::vector<double> const &x, std::vector<double> const &y, int numx, double tolerance, TasmanianSparseGrid const &grid, std::string message){
    AccelerationContext const *acceleration = grid.getAccelerationContext();
    GpuVector<T> gpux;
    gpux.load(acceleration, x);
    GpuVector<T> gpuy(acceleration, ((grid.isFourier() && std::is_same<GridMethod, GridMethodHierBasisGPU>::value) ? 2 : 1) * numx,
                       (std::is_same<GridMethod, GridMethodHierBasisGPU>::value) ? grid.getNumPoints() : grid.getNumOutputs());
    if (std::is_same<GridMethod, GridMethodEvalBatchGPU>::value){
        grid.evaluateBatchGPU(gpux.data(), numx, gpuy.data());
    }else if (std::is_same<GridMethod, GridMethodHierBasisGPU>::value){
        grid.evaluateHierarchicalFunctionsGPU(gpux.data(), numx, gpuy.data());
    }
    auto cpuy = gpuy.unload(acceleration);
    return testPass(err1(cpuy, y), tolerance, message, grid);
}
template<typename T>
bool testHBasisGPUSparse(std::vector<double> const &x,
                         std::vector<int> const &pntr, std::vector<int> const &indx, std::vector<double> const &vals,
                         double tolerance, TasmanianSparseGrid const &grid, std::string message){
    if (grid.empty()){ cout << "ERROR: cannot test an empty grid\n"; return false; }
    AccelerationContext const *acceleration = grid.getAccelerationContext();
    GpuVector<T> gpux;
    gpux.load(acceleration, x);

    int nump = (int) x.size() / grid.getNumDimensions();
    int *gpu_indx = 0, *gpu_pntr = 0, num_nz = 0;
    T *gpu_vals = 0;

    grid.evaluateSparseHierarchicalFunctionsGPU(gpux.data(), nump, gpu_pntr, gpu_indx, gpu_vals, num_nz);

    std::vector<int> cpntr; AccelerationMeta::recvGpuArray(acceleration, nump + 1, gpu_pntr, cpntr);
    std::vector<int> cindx; AccelerationMeta::recvGpuArray(acceleration, num_nz, gpu_indx, cindx);
    std::vector<T> cvals;   AccelerationMeta::recvGpuArray(acceleration, num_nz, gpu_vals, cvals);

    AccelerationMeta::delGpuArray(acceleration, gpu_pntr);
    AccelerationMeta::delGpuArray(acceleration, gpu_indx);
    AccelerationMeta::delGpuArray(acceleration, gpu_vals);

    return testPass(err1_sparse(pntr, indx, vals, cpntr, cindx, cvals), tolerance, message, grid);
}
#endif


#endif
