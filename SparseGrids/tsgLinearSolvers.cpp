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

#ifndef __TASMANIAN_LINEAR_SOLVERS_CPP
#define __TASMANIAN_LINEAR_SOLVERS_CPP

#include "tsgLinearSolvers.hpp"
#include "tsgTPLWrappers.hpp"

namespace TasGrid{

void TasmanianDenseSolver::solveLeastSquares(int n, int m, const double A[], double b[], double *x){
    // form Ar = A' * A
    // form x = A' * b
    std::vector<double> Ar(m * m);
    #pragma omp parallel for
    for(int i=0; i<m; i++){
        for(int j=0; j<m; j++){
            double sum = 0.0;
            for(int k=0; k<n; k++){
                sum += A[i*n + k] * A[j*n + k];
            }
            Ar[i*m + j] = sum;
        }
        x[i] = 0.0;
        for(int k=0; k<n; k++){
            x[i] += A[i*n+k] * b[k];
        }
    }

    // factorize Ar
    for(int i=0; i<m; i++){
        Ar[i*m + i] = std::sqrt(Ar[i*m + i]);

        for(int j=i+1; j<m; j++){
            Ar[i*m + j] /= Ar[i*m + i];
        }

        for(int k=i+1; k<m; k++){
            for(int j=i+1; j <= k; j++){
                Ar[j*m + k] -= Ar[i*m + j] * Ar[i*m + k];
            }
        }
    }

    // solve L^(-1) x
    for(int i=0; i<m; i++){
        x[i] /= Ar[i*m + i];
        for(int j=i+1; j<m; j++){
            x[j] -= x[i] * Ar[i*m + j];
        }
    }

    // solve L^(-T) x
    for(int i=m-1; i>=0; i--){
        for(int j=i+1; j<m; j++){
            x[i] -= x[j] * Ar[i*m + j];
        }
        x[i] /= Ar[i*m + i];
    }
}

void TasmanianDenseSolver::solveLeastSquares(AccelerationContext const *acceleration, int n, int m, double A[], double b[], double *x){
    if (acceleration->blasCompatible()){
        TasBLAS::solveLS('N', n, m, A, b);
        std::copy_n(b, m, x);
        return;
    }else{
        solveLeastSquares(n, m, A, b, x);
    }
}

template<typename scalar_type>
void TasmanianDenseSolver::solvesLeastSquares(AccelerationContext const *acceleration, int n, int m, scalar_type A[], int nrhs, scalar_type B[]){
    switch(acceleration->mode){
        case accel_gpu_magma:
            TasGpu::solveLSmultiOOC(acceleration, n, m, A, nrhs, B);
            break;
        case accel_gpu_cuda:
        case accel_gpu_cublas:
            acceleration->setDevice();
            TasGpu::solveLSmulti(acceleration, n, m, A, nrhs, B);
            break;
        case accel_cpu_blas:
            TasBLAS::solveLSmulti(n, m, A, nrhs, B);
            break;
        default:
            throw std::runtime_error("Dense least-squares solve attempted without BLAS or CUDA acceleration enabled.");
    };
}
template<typename scalar_type>
void TasmanianDenseSolver::solvesLeastSquaresGPU(AccelerationContext const *acceleration, int n, int m, scalar_type A[], int nrhs, scalar_type B[]){
    if (not acceleration->on_gpu())
        throw std::runtime_error("solvesLeastSquaresGPU() requires a GPU mode to be enabled.");
    TasGpu::solveLSmultiGPU(acceleration, n, m, A, nrhs, B);
}

template void TasmanianDenseSolver::solvesLeastSquares<double>(AccelerationContext const*, int, int, double[], int, double[]);
template void TasmanianDenseSolver::solvesLeastSquares<std::complex<double>>(AccelerationContext const*, int, int, std::complex<double>[], int, std::complex<double>[]);
template void TasmanianDenseSolver::solvesLeastSquaresGPU<double>(AccelerationContext const*, int, int, double[], int, double[]);
template void TasmanianDenseSolver::solvesLeastSquaresGPU<std::complex<double>>(AccelerationContext const*, int, int, std::complex<double>[], int, std::complex<double>[]);

namespace Utils{ // implements block transpose method
template<typename scalar_type>
inline void copy_transpose(long long M, long long N, scalar_type const A[], long long lda, scalar_type B[], long long ldb){
    for(int i=0; i<M; i++)
        for(int j=0; j<N; j++)
            B[i * ldb + j] = A[j * lda + i];
}
template<typename scalar_type>
void transpose(long long M, long long N, scalar_type const A[], scalar_type B[]){
    constexpr long long bsize = 64;
    long long bM = (M / bsize) + ((M % bsize == 0) ? 0 : 1); // number of blocks in M and N
    long long bN = (N / bsize) + ((N % bsize == 0) ? 0 : 1);
    #pragma omp parallel for
    for(long long t =0; t < bM * bN; t++){
        long long i = t / bN;
        long long j = t % bN;
        copy_transpose(std::min(bsize, M - i * bsize), std::min(bsize, N - j * bsize),
                        A + i * bsize + j * bsize * M, M,
                        B + i * bsize * N + j * bsize, N);
    }
}

template void transpose(long long, long long, float const[], float[]);
template void transpose(long long, long long, double const[], double[]);
template void transpose(long long, long long, std::complex<double> const[], std::complex<double>[]);
}

std::vector<double> TasmanianTridiagonalSolver::getSymmetricEigenvalues(int n, std::vector<double> const &diag, std::vector<double> const &offdiag){
    #ifndef Tasmanian_ENABLE_BLAS
    throw std::runtime_error("getSymmetricEigenvalues() requires BLAS acceleration to be enabled!");
    #endif
    std::vector<double> result = diag;
    std::vector<double> dummy = offdiag;
    TasBLAS::sterf(n, result.data(), dummy.data());
    return result;
}

void TasmanianTridiagonalSolver::decompose(std::vector<double> &diag, std::vector<double> &off_diag, const double mu0,
                                           std::vector<double> &nodes, std::vector<double> &weights) {
    switch(TasmanianTridiagonalSolver::decompose_version) {
        case 1 :
            weights = std::vector<double>(diag.size(), 0.0);
            weights[0] = sqrt(mu0);
            nodes = diag;
            off_diag.push_back(0.0);
            decompose1(static_cast<int>(diag.size()), nodes, off_diag, weights);
            break;
        case 2 :
            decompose2(diag, off_diag, mu0, nodes, weights);
            break;
        default :
            throw std::invalid_argument("ERROR: decompose_version must be a valid number!");
    }
}

void TasmanianTridiagonalSolver::decompose1(int n, std::vector<double> &d, std::vector<double> &e, std::vector<double> &z){
    const double tol = Maths::num_tol;
    if (n == 1){ z[0] = z[0]*z[0]; return; }

    for(int l=0; l<n-1; l++){
        int m = l;
        while((m < n-1) && (std::abs(e[m]) > tol)) m++;

        while (m != l){
            double p = d[l];
            double g = (d[l+1] - p) / (2.0 * e[l]);
            double r = std::sqrt(g*g + 1.0);

            g = d[m] - p + e[l] / (g + Maths::sign(g) * r); // sign function here may be unstable

            double s = 1.0;
            double c = 1.0;
            p = 0.0;

            for(int i=m-1; i>=l; i--){
                double f = s * e[i];
                double b = c * e[i];

                if (std::abs(f) >= std::abs(g)){
                    c = g / f;
                    r = std::sqrt(c*c + 1.0);
                    e[i+1] = f*r;
                    s = 1.0 / r;
                    c *= s;
                }else{
                    s = f / g;
                    r =  std::sqrt(s*s + 1.0);
                    e[i+1] = g * r;
                    c = 1.0 / r;
                    s *= c;
                }

                g = d[i+1] - p;
                r = (d[i] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i+1] = g + p;
                g = c * r - b;
                f = z[i+1];
                z[i+1] = s * z[i] + c * f;
                z[i] = c * z[i] - s * f;
            }

            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;

            m = l;
            while((m < n-1) && (std::abs(e[m]) > tol)) m++;
        }
    }

    for(int i=1; i<n; i++){
        for(int j=0; j<n-1; j++){
            if (d[j] > d[j+1]){
                double p = d[j];
                d[j] = d[j+1];
                d[j+1] = p;
                p = z[j];
                z[j] = z[j+1];
                z[j+1] = p;
            }
        }
    }
    for(int i=0; i<n; i++){
           z[i] = z[i]*z[i];
    }
}

void TasmanianTridiagonalSolver::decompose2(std::vector<double> &diag, std::vector<double> &off_diag, const double mu0,
                                            std::vector<double> &nodes, std::vector<double> &weights) {

    // Ensure compatibility with the ALGOL implementation.
    // NOTE: diag and off_diag are 1-indexed, while nodes and weights are 0-indexed after this step.
    size_t n = diag.size();
    assert(off_diag.size() == n-1);
    nodes.resize(n);
    weights.resize(n);
    diag.insert(diag.begin(), 0.0);
    off_diag.insert(off_diag.begin(), 0.0);

    // SETUP block from ALGOL code.
    // ALGOL COMMENT: Find the maximum row sum norm and initialize weights.
    off_diag[0] = 0.0;
    double norm = 0.0;
    for (size_t i=1; i<=n-1; i++) {
        norm = std::max(norm, std::fabs(off_diag[i-1]) + std::fabs(diag[i]) + std::fabs(off_diag[i-1]));
        weights[i-1] = 0.0;
    }
    norm = std::max(norm, std::fabs(diag[n]) + std::fabs(off_diag[n-1]));
    weights[n-1] = 0.0;
    weights[0] = 1.0; // Fix the bug in the ALGOL code.
    double eps = norm * std::numeric_limits<double>::epsilon();
    size_t m = n;
    double lambda{norm}, lambda1{norm}, lambda2{norm}, rho{norm};

    // INSPECT block from ALGOL code.
    // ALGOL COMMENT: Look for convergence of lower diagonal element.
    while (m > 0) {
        if (std::fabs(off_diag[m-1]) <= eps) {
            nodes[m-1] = diag[m];
            weights[m-1] = mu0 * weights[m-1] * weights[m-1];
            rho = lambda1 < lambda1 ? lambda1 : lambda2;
            m--;
            continue;
        }
        // ALGOL COMMENT: Small off diagonal element means matrix can be split.
        size_t k = m-1;
        while (std::fabs(off_diag[k-1]) > eps) k--;
        // ALGOL COMMENT: Find eigenvalues of lower 2-by-2 and select accelerating shift.
        double b2 = off_diag[m-1] * off_diag[m-1];
        double det = std::sqrt((diag[m-1] - diag[m]) * (diag[m-1] - diag[m]) + 4.0 * b2);
        double aa = diag[m-1] + diag[m];
        lambda2 = 0.5 * (aa >= 0 ? aa + det : aa - det);
        lambda1 = (diag[m-1] * diag[m] - b2) / lambda2;
        double eigmax = std::max(lambda1, lambda2);
        if (std::fabs(eigmax - rho) <= 0.125 * std::fabs(eigmax)) {
            lambda = eigmax;
        }
        rho = eigmax;
        // ALGOL COMMENT: Transform block from k to m.
        double cj = off_diag[k];
        off_diag[k-1] = diag[k] - lambda;
        for (size_t j=k; j<=m-1; j++) {
            double r = std::sqrt(cj * cj + off_diag[j-1] * off_diag[j-1]);
            double st = cj / r;            double st2 = st * st;
            double ct = off_diag[j-1] / r; double ct2 = ct * ct;
            double sc = st * ct;           double aj = diag[j];
            double bj = off_diag[j];       double wj = weights[j-1];
            // Order below is important!
            diag[j] = aj * ct2 + 2.0 * bj * sc + diag[j+1] * st2;
            off_diag[j] = (aj - diag[j+1]) * sc + bj * (st2 - ct2);
            diag[j+1] = aj * st2 - 2.0 * bj * sc + diag[j+1] * ct2;
            cj = off_diag[j+1] * st;
            off_diag[j+1] = -off_diag[j+1] * ct;
            off_diag[j-1] = r;
            // Account for the offset of the indices.
            weights[j-1] = wj * ct + weights[j] * st;
            weights[j] = wj * st - weights[j] * ct;
        }
        off_diag[k-1] = 0.0;
    }

    // SORT block from ALGOL code.
    // ALGOL COMMENT: Arrange abscissas in ascending order.
    // NOTE: the original code used an exchange sort, which is O(n^2).
    std::vector<size_t> I(n);
    for (size_t i=0; i<n; i++) I[i] = i;
    std::sort(I.begin(), I.end(), [&nodes](size_t i, size_t j){return nodes[i] < nodes[j];});
    for (size_t i=0; i<n; i++) {
        while (I[i] != i) {
            std::swap(nodes[i], nodes[I[i]]);
            std::swap(weights[i], weights[I[i]]);
            std::swap(I[i], I[I[i]]);
        }
    }
}

void TasmanianFourierTransform::fast_fourier_transform(std::vector<std::vector<std::complex<double>>> &data, std::vector<int> &num_points){
    int num_dimensions = (int) num_points.size();
    int num_total = 1;
    for(auto n: num_points) num_total *= n;
    std::vector<int> cumulative_points(num_dimensions);
    for(int k=0; k<num_dimensions; k++){
        // split the data into vectors of 1-D transforms
        std::vector<std::vector<int>> maps1d(num_total / num_points[k]); // total number of 1-D transforms
        for(auto &v: maps1d) v.reserve(num_points[k]); // each 1-D transform will have size num_points[k]

        std::fill(cumulative_points.begin(), cumulative_points.end(), 1);
        for(int j=num_dimensions-2; j>=0; j--) cumulative_points[j] = cumulative_points[j+1] * ((j+1 != k) ? num_points[j+1] : 1);

        // convert i to tuple, i.e., i -> p0, p1, p2 .. pd, where d = num_dimensions
        // all tuples that differ only in the k-th entry will belong to the same 1-D transform
        // the index of the 1D transform is i1d = \sum_{j \neq k} pj * cumulative_points[j]
        for(int i=0; i<num_total; i++){
            int t = i; // used to construct the tensor index
            int i1d = 0; // index of the corresponding 1D transform
            for(int j=num_dimensions-1; j>=0; j--){
                int pj = t % num_points[j]; // tensor index j
                if (j != k) i1d += cumulative_points[j] * pj;
                t /= num_points[j];
            }
            maps1d[i1d].push_back(i);
        }

        #pragma omp parallel for // perform the 1D transforms
        for(int i=0; i<(int) maps1d.size(); i++){
            fast_fourier_transform1D(data, maps1d[i]);
        }
    }
}

void TasmanianFourierTransform::fast_fourier_transform1D(std::vector<std::vector<std::complex<double>>> &data, std::vector<int> &indexes){
    //
    // Given vector x_n with size N, the Fourier transform F_k is defined as: F_k = \sum_{n=0}^{N-1} \exp(- 2 \pi k n / N) x_n
    // Assuming that N = 3^l for some l, we can sub-divide the transform into strips of 3
    // let j = 0, 1, 2; let k = 0, .., N/3; let x_{0,m} = x_{3m}; let x_{1,m} = x_{3m+1}; and let x_{2,m} = x_{3m+2}
    // F_{k + j N / 3} =                                         \sum_{m=0}^{N/3 - 1} x_{0,m} \exp(-2 \pi k m / (N / 3))
    //                   + \exp(-2 \pi k / N) \exp(-2 \pi j / 3) \sum_{m=0}^{N/3 - 1} x_{1,m} \exp(-2 \pi k m / (N / 3))
    //                   + \exp(-4 \pi k / N) \exp(-4 \pi j / 3) \sum_{m=0}^{N/3 - 1} x_{2,m} \exp(-2 \pi k m / (N / 3))
    // The three sums are the Fourier coefficients of x_{0, m}, x_{1, m}, and x_{2, m}
    // The terms \exp(-2 \pi k / N) \exp(-2 \pi j / 3), and \exp(-4 \pi k / N) \exp(-4 \pi j / 3) are the twiddle factors
    // The procedure is recursive splitting the transform into small sets, all the way to size 3
    //
    int num_outputs = (int) data[0].size(); // get the problem dimensions, num outputs and num entries for the 1D transform
    int num_entries = (int) indexes.size(); // the size of the 1D problem, i.e., N
    if (num_entries == 1) return; // nothing to do for size 1
    // a copy of the data is needed to swap back and forth, thus we make two copies and swap between them
    std::vector<std::vector<std::complex<double>>> V(num_entries);
    auto v = V.begin();
    for(auto i: indexes) *v++ = data[i]; // copy from the data only the indexes needed for the 1D transform
    std::vector<std::vector<std::complex<double>>> W(num_entries);
    for(auto &w : W) w.resize(num_outputs); // allocate storage for the second data set

    // the radix-3 FFT algorithm uses two common twiddle factors from known angles +/- 2 pi/3
    std::complex<double> twidlep(-0.5, -std::sqrt(3.0) / 2.0); // angle of -2 pi/3
    std::complex<double> twidlem(-0.5,  std::sqrt(3.0) / 2.0); // angle of  2 pi/3 = -4 pi/3

    int stride = num_entries / 3; // the jump between entries, e.g., in one level of split stride is 3, split again and stride is 9 ... up to N / 3
    int length = 3;               // the number of entries in the sub-sequences, i.e., how large k can be (see above), smallest sub-sequence uses length 3

    for(int i=0; i<stride; i++){ // do the 3 transform, multiply by 3 by 3 matrix
        auto x1 = V[i].begin();
        auto x2 = V[i + stride].begin();
        auto x3 = V[i + 2 * stride].begin(); // x1, x2, x3 are the three entries of a sub-sequence

        auto y1 = W[i].begin();
        auto y2 = W[i + stride].begin();
        auto y3 = W[i + 2 * stride].begin(); // y1, y2, y3 are the resulting Fourier coefficients

        for(int k=0; k<num_outputs; k++){
            *y1++ = *x1 + *x2 + *x3;
            *y2++ = *x1 + twidlep * *x2 + twidlem * *x3;
            *y3++ = *x1 + twidlem * *x2 + twidlep * *x3;
            x1++; x2++; x3++;
        }
    }

    std::swap(V, W); // swap, now V contains the computed transform of the sub-sequences with size 3, W will be used for scratch space

    // merge smaller sequences, do the recursion
    while(stride / 3 > 0){ // when the stride that we just computed is equal to 1, then stop the recursion
        int biglength = 3 * length; // big sequence, i.e., F_k has this total length
        int bigstride = stride / 3;

        double theta = -2.0 * Maths::pi / ((double) biglength);
        std::complex<double> expstep(std::cos(theta), std::sin(theta)); // initialize the twiddle factors common for this level of sub-sequences
        std::complex<double> expstep2 = expstep * expstep;

        // merge sets of 3 sub-sequences (coefficients of x_{i,m}) into 3 pieces of one sequence F_{k + j N / 3}
        for(int i=0; i<bigstride; i++){ // total number of triples of sequences is bigstride
            std::complex<double> t01(1.0, 0.0);
            std::complex<double> t02(1.0, 0.0);

            std::complex<double> t11 = twidlep;
            std::complex<double> t12 = twidlem;

            std::complex<double> t21 = twidlem;
            std::complex<double> t22 = twidlep; // the twiddle factors form a 3 by 3 matrix [1, 1, 1; 1, t11, t12; 1, t21, t22;]

            for(int k=0; k<length; k++){ // number of entries in the sub-sequences
                auto x1 = V[i + k * stride].begin();
                auto x2 = V[i + k * stride + bigstride].begin();
                auto x3 = V[i + k * stride + 2 * bigstride].begin(); // x1, x2, x3 are the next entries of the sub-sequence (i.e., the sums)

                auto y1 = W[i + k * bigstride].begin();
                auto y2 = W[i + (k + length) * bigstride].begin();
                auto y3 = W[i + (k + 2 * length) * bigstride].begin(); // y1, y2, y3 are the F_{k + j N / 3}

                for(int o=0; o<num_outputs; o++){ // traverse through all the outputs
                    *y1++ = *x1 + t01 * *x2 + t02 * *x3;
                    *y2++ = *x1 + t11 * *x2 + t12 * *x3;
                    *y3++ = *x1 + t21 * *x2 + t22 * *x3;
                    x1++; x2++; x3++;
                }

                // update the twiddle factors for the next index k
                t01 *= expstep;
                t11 *= expstep;
                t21 *= expstep;
                t02 *= expstep2;
                t12 *= expstep2;
                t22 *= expstep2;
            }
        }

        std::swap(V, W); // swap the data, V holds the current set of indexes and W is the next set

        stride = bigstride;
        length = biglength;
    }

    // copy back the solution into the data structure
    v = V.begin();
    for(auto i : indexes) data[i] = *v++;
}

namespace TasSparse{

WaveletBasisMatrix::WaveletBasisMatrix(AccelerationContext const *acceleration,
                                       const std::vector<int> &lpntr, const std::vector<std::vector<int>> &lindx, const std::vector<std::vector<double>> &lvals) : tol(Maths::num_tol), num_rows(static_cast<int>(lpntr.size())){

    if (num_rows == 0) return; // make an empty matrix

    // hip doesn't have rocsolver yet, so BLAS is required for dense operations
    if (acceleration->mode != accel_none and useDense(acceleration, num_rows)){ // dense mode
        dense = std::vector<double>(Utils::size_mult(num_rows, num_rows));
        Utils::Wrapper2D<double> dense_rows(num_rows, dense.data());

        auto idx = lindx.begin();
        auto vls = lvals.begin();
        auto ii = idx->begin();
        auto vv = vls->begin();
        for(int i=0; i<num_rows; i++){
            if (ii == idx->end()){
                idx++;
                vls++;
                if (idx != lindx.end()){
                    ii = idx->begin();
                    vv = vls->begin();
                }
            }

            double *r = dense_rows.getStrip(i);
            for(int j=0; j<lpntr[i]; j++)
                r[*ii++] = *vv++;
        }

        if (acceleration->mode != accel_cpu_blas){ // using GPU
            acceleration->setDevice();
            gpu_dense = GpuVector<double>(acceleration, dense);
            dense = std::vector<double>();
        }
    }else{ // sparse mode
        pntr = std::vector<int>(num_rows+1, 0);
        for(int i=0; i<num_rows; i++)
            pntr[i+1] = pntr[i] + lpntr[i];

        int num_nz = pntr.back();
        indx.reserve(num_nz);
        vals.reserve(num_nz);

        for(const auto &idx1 : lindx)
            indx.insert(indx.end(), idx1.begin(), idx1.end());
        for(const auto &vls1 : lvals)
            vals.insert(vals.end(), vls1.begin(), vls1.end());
    }
    factorize(acceleration);
}

void WaveletBasisMatrix::factorize(AccelerationContext const *acceleration){
    if (not gpu_dense.empty()){
        gpu_ipiv = GpuVector<int_gpu_lapack>(acceleration, num_rows);
        TasGpu::factorizePLU(acceleration, num_rows, gpu_dense.data(), gpu_ipiv.data());
    }else if (not dense.empty()){
        ipiv = std::vector<int>(num_rows);
        TasBLAS::getrf(num_rows, num_rows, dense.data(), num_rows, ipiv.data());
    }else{
        computeILU();
    }
}

void WaveletBasisMatrix::computeILU(){
    indxD.resize(num_rows);
    ilu.resize(pntr[num_rows]);
    for(int i=0; i<num_rows; i++){
        int j = pntr[i];
        while(indx[j] < i){ j++; };
        indxD[i] = j;
    }

    ilu = vals;

    for(int i=0; i<num_rows-1; i++){
        double u = ilu[indxD[i]];
        #pragma omp parallel for
        for(int j=i+1; j<num_rows; j++){ // update the rest of the matrix, each row can be done in parallel
            int jc = pntr[j];
            while(indx[jc] < i){ jc++; }
            if (indx[jc] == i){
                ilu[jc] /= u;
                double l = ilu[jc];
                int ik = indxD[i]+1;
                int jk = jc+1;
                while((ik<pntr[i+1]) && (jk<pntr[j+1])){
                    if (indx[ik] == indx[jk]){
                        ilu[jk] -= l * ilu[ik];
                        ik++; jk++;
                    }else if (indx[ik] < indx[jk]){
                        ik++;
                    }else{
                        jk++;
                    }
                }
            }
        }
    }
}

void WaveletBasisMatrix::invertTransposed(AccelerationContext const *acceleration, double b[]) const{
    if (not gpu_dense.empty()){
        GpuVector<double> gpu_b(acceleration, num_rows, 1, b);
        TasGpu::solvePLU(acceleration, 'N', num_rows, gpu_dense.data(), gpu_ipiv.data(), gpu_b.data());
        gpu_b.unload(acceleration, b);
    }else if (not dense.empty()){
        TasBLAS::getrs('N', num_rows, 1, dense.data(), num_rows, ipiv.data(), b, num_rows);
    }else{
        // using sparse algorithm
        if (acceleration->blasCompatible())
            solve<use_transpose, use_blas>(std::vector<double>(b, b + num_rows).data(), b);
        else
            solve<use_transpose, no_blas>(std::vector<double>(b, b + num_rows).data(), b);
    }
}

void WaveletBasisMatrix::invert(AccelerationContext const *acceleration, int num_colums, double B[]){
    if (not gpu_dense.empty()){
        GpuVector<double> gpu_b(acceleration, num_rows, num_colums, B);
        if (num_colums == 1){
            TasGpu::solvePLU(acceleration, 'T', num_rows, gpu_dense.data(), gpu_ipiv.data(), gpu_b.data());
        }else{
            TasGpu::solvePLU(acceleration, 'T', num_rows, gpu_dense.data(), gpu_ipiv.data(), num_colums, gpu_b.data());
        }
        gpu_b.unload(acceleration, B);
    }else if (not dense.empty()){
        if (num_colums == 1){
            TasBLAS::getrs('T', num_rows, 1, dense.data(), num_rows, ipiv.data(), B, num_rows);
        }else{
            TasBLAS::trsm('R', 'U', 'N', 'N', num_colums, num_rows, 1.0, dense.data(), num_rows, B, num_colums);
            TasBLAS::trsm('R', 'L', 'N', 'U', num_colums, num_rows, 1.0, dense.data(), num_rows, B, num_colums);
            // permute
            Utils::Wrapper2D<double> rows(num_colums, B);
            for(int i=0; i<num_rows; i++){
                if (ipiv[i] - 1 != i){
                    TasBLAS::vswap(num_colums, rows.getStrip(i), 1, rows.getStrip(ipiv[i] - 1), 1);
                }
            }
        }
    }else{
        if (num_colums == 1){
            std::vector<double> b(B, B + num_rows);
            if (acceleration->blasCompatible())
                solve<no_transpose, use_blas>(b.data(), B);
            else
                solve<no_transpose, no_blas>(b.data(), B);
            return;
        }
        Utils::Wrapper2D<double> wrapb(num_colums, B);
        std::vector<double> b(num_rows);
        std::vector<double> x(num_rows);
        for(int k=0; k<num_colums; k++){
            for(int i=0; i<num_rows; i++){
                b[i] = wrapb.getStrip(i)[k];
                x[i] = b[i];
            }

            if (acceleration->blasCompatible())
                solve<no_transpose, use_blas>(b.data(), x.data());
            else
                solve<no_transpose, no_blas>(b.data(), x.data());
            for(int i=0; i<num_rows; i++)
                wrapb.getStrip(i)[k] = x[i];
        }
    }
}

template<bool transpose> void WaveletBasisMatrix::applyILU(double x[]) const{
    if (transpose){
        for(int i=0; i<num_rows; i++){
            x[i] /= ilu[indxD[i]];
            for(int j=indxD[i]+1; j<pntr[i+1]; j++)
                x[indx[j]] -= ilu[j] * x[i];
        }
        for(int i=num_rows-2; i>=0; i--)
            for(int j=pntr[i]; j<indxD[i]; j++)
                x[indx[j]] -= ilu[j] * x[i];
    }else{
        for(int i=1; i<num_rows; i++){
            for(int j=pntr[i]; j<indxD[i]; j++){
                x[i] -= ilu[j] * x[indx[j]];
            }
        }
        for(int i=num_rows-1; i>=0; i--){
            for(int j=indxD[i]+1; j<pntr[i+1]; j++){
                x[i] -= ilu[j] * x[indx[j]];
            }
            x[i] /= ilu[indxD[i]];
        }
    }
}
template<bool transpose> void WaveletBasisMatrix::apply(double const x[], double r[]) const{
    if (transpose){
        std::fill_n(r, num_rows, 0.0);
        for(int i=0; i<num_rows; i++){
            for(int j=pntr[i]; j<pntr[i+1]; j++)
                r[indx[j]] += vals[j] * x[i];
        }
    }else{
        for(int i=0; i<num_rows; i++){
            double sum = 0.0;
            for(int j=pntr[i]; j<pntr[i+1]; j++)
                sum += vals[j] * x[indx[j]];
            r[i] = sum;
        }
    }
}

void WaveletBasisMatrix::residual(double const x[], double const b[], double r[]) const{
    for(int i=0; i<num_rows; i++){
        double sum = 0.0;
        for(int j=pntr[i]; j<pntr[i+1]; j++)
            sum += vals[j] * x[indx[j]];
        r[i] = b[i] - sum;
    }
}

// compute the norm, scale x by 1 / norm, return the norm
inline double rescale(int num_entries, double x[]){
    double nrm = 0.0;
    for(int i=0; i<num_entries; i++) nrm += x[i] * x[i];
    nrm = std::sqrt(nrm);
    double s = 1.0 / nrm;
    if (nrm > 0.0) for(int i=0; i<num_entries; i++) x[i] *= s;
    return nrm;
}
inline double rescale_blas(int num_entries, double x[]){
    if (AccelerationMeta::isAvailable(accel_cpu_blas)){
        double nrm = TasBLAS::norm2(num_entries, x, 1);
        if (nrm > 0.0) TasBLAS::scal(num_entries, 1.0 / nrm, x, 1);
        return nrm;
    }else{
        return rescale(num_entries, x);
    }
}

// project the krylov basis
inline void projectKrylov(int inner_itr, int max_inner, int num_rows, std::vector<double> &W, std::vector<double> &H){
    #pragma omp parallel for
    for(int i=0; i<inner_itr; i++){
        H[i*max_inner + inner_itr-1] = 0.0; for(int j=0; j<num_rows; j++){  H[i*max_inner + inner_itr-1] += W[inner_itr*num_rows+j] * W[i*num_rows+j];  };
    }

    #pragma omp parallel for
    for(int j=0; j<num_rows; j++){
        for(int i=0; i<inner_itr; i++){
            W[inner_itr*num_rows+j] -= H[i*max_inner + inner_itr-1] * W[i*num_rows+j];
        }
    }
}
inline void projectKrylov_blas(int inner_itr, int max_inner, int num_rows, std::vector<double> &W, std::vector<double> &H){
    if (AccelerationMeta::isAvailable(accel_cpu_blas)){
        TasBLAS::gemv('T', num_rows, inner_itr,  1.0, W.data(), num_rows, &W[inner_itr * num_rows], 1, 0.0, &H[inner_itr-1], max_inner);
        TasBLAS::gemv('N', num_rows, inner_itr, -1.0, W.data(), num_rows, &H[inner_itr-1], max_inner, 1.0, &W[inner_itr * num_rows], 1);
    }else{
        projectKrylov(inner_itr, max_inner, num_rows, W, H);
    }
}
// reconstruct the krylov solution from the basis
inline void reconstructKrylov(int inner_itr, int max_inner, int num_rows, std::vector<double> const &W, std::vector<double> const &H,
                              std::vector<double> &Z, double x[]){
    Z[inner_itr] /= H[inner_itr * max_inner + inner_itr];
    for(int i=inner_itr-1; i>-1; i--){
        double beta = 0.0;
        for(int j=i+1; j<=inner_itr; j++){
            beta += H[i*max_inner + j] * Z[j];
        };
        Z[i] = (Z[i] - beta) / H[i * max_inner + i];
    }

    for(int i=0; i<=inner_itr; i++){
        for(int j=0; j<num_rows; j++){
            x[j] += Z[i] * W[i*num_rows+j];
        }
    }
}
inline void reconstructKrylov_blas(int inner_itr, int max_inner, int num_rows, std::vector<double> const &W, std::vector<double> const &H,
                                   std::vector<double> &Z, double x[]){
    if (AccelerationMeta::isAvailable(accel_cpu_blas)){
        TasBLAS::trsv('L', 'T', 'N', inner_itr, H.data(), max_inner, Z.data(), 1);
        TasBLAS::gemv('N', num_rows, inner_itr, 1.0, W.data(), num_rows, Z.data(), 1, 1.0, x, 1);
    }else{
        reconstructKrylov(inner_itr, max_inner, num_rows, W, H, Z, x);
    }
}

template<bool transpose, bool blas>
void WaveletBasisMatrix::solve(const double b[], double x[]) const{
    int max_inner = 30;
    int max_outer = 80;
    std::vector<double> W((max_inner+1) * num_rows); // Krylov basis

    std::vector<double> H(max_inner * (max_inner+1)); // holds the transformation for the normalized basis
    std::vector<double> S(max_inner); // std::sin and std::cos of the Givens rotations
    std::vector<double> C(max_inner+1);
    std::vector<double> Z(max_inner); // holds the coefficients of the solution

    double alpha, beta; // temp variables

    double outer_res = tol + 1.0; // outer and inner residual
    int outer_itr = 0; // counts the inner and outer iterations

    std::vector<double> temp;

    std::fill_n(x, num_rows, 0.0); // zero initial guess, I wonder if we can improve this

    while (outer_res > tol and outer_itr < max_outer){
        if (transpose){
            temp = std::vector<double>(x, x + num_rows);
            applyILU<transpose>(temp.data());
            apply<transpose>(temp.data(), W.data());
            for(int i=0; i<num_rows; i++){
                W[i] = b[i] - W[i];
            }
        }else{
            residual(x, b, W.data());
            applyILU<transpose>(W.data());
        }

        double inner_res = (blas) ? rescale_blas(num_rows, W.data()) : rescale(num_rows, W.data());
        Z[0] = inner_res;

        int inner_itr = 0; // counts the size of the basis

        while ((inner_res > tol) && (inner_itr < max_inner-1)){
            inner_itr++;

            if (transpose){
                std::copy_n(&(W[num_rows*(inner_itr-1)]), num_rows, temp.data());
                applyILU<transpose>(temp.data());
                apply<transpose>(temp.data(), &W[inner_itr*num_rows]);
            }else{
                apply<transpose>(&W[num_rows*(inner_itr-1)], &W[inner_itr*num_rows]);
                applyILU<transpose>(&W[inner_itr*num_rows]);
            }

            if (blas)
                projectKrylov_blas(inner_itr, max_inner, num_rows, W, H);
            else
                projectKrylov(inner_itr, max_inner, num_rows, W, H);

            beta = (blas) ? rescale_blas(num_rows, &W[inner_itr*num_rows]) : rescale(num_rows, &W[inner_itr*num_rows]);

            for (int i=0; i<inner_itr-1; i++){ // form the next row of the transformation
                alpha = H[i*max_inner + inner_itr-1];
                H[   i*max_inner + inner_itr-1] = C[i] * alpha + S[i] * H[(i+1)*max_inner + inner_itr-1];
                H[(i+1)*max_inner + inner_itr-1] = S[i] * alpha - C[i] * H[(i+1)*max_inner + inner_itr-1];
            };

            alpha = std::sqrt(beta * beta  +  H[(inner_itr-1)*max_inner + inner_itr-1] * H[(inner_itr-1)*max_inner + inner_itr-1]);

            // set the next set of Givens rotations
            S[inner_itr-1] = beta / alpha;
            C[inner_itr-1] = H[(inner_itr-1)*max_inner + inner_itr-1] / alpha;

            H[(inner_itr-1) * max_inner + inner_itr-1] = alpha;

            // Z is used to reconstruct the solution in the end
            Z[inner_itr] = S[inner_itr-1]*Z[inner_itr-1];
            Z[inner_itr-1] = C[inner_itr-1]*Z[inner_itr-1]; // apply it on z

            inner_res = std::abs(Z[inner_itr]);
        }

        inner_itr--;

        if (inner_itr > -1){ // if the first guess was not within TOL of the true solution
            if (blas)
                reconstructKrylov_blas(inner_itr, max_inner, num_rows, W, H, Z, x);
            else
                reconstructKrylov(inner_itr, max_inner, num_rows, W, H, Z, x);

            if (transpose)
                applyILU<transpose>(x);

        }

        outer_res = inner_res;
        outer_itr++;
    }
}

} /* namespace TasSparse */

}
#endif
