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

#ifndef __TASMANIAN_LINEAR_SOLVERS_HPP
#define __TASMANIAN_LINEAR_SOLVERS_HPP

#include <complex>

#include "tsgAcceleratedDataStructures.hpp"

//! \internal
//! \file tsgLinearSolvers.hpp
//! \brief Linear solvers.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianLinearSolvers
//!
//! Many sparse grids methods rely on various linear solvers, sparse iterative, eigenvalue,
//! fast-fourier-transform, and least-squares. The solvers are tuned for the specific
//! problems and can perform as good (or even better) than general third-party methods.
//! In addition, this reduces the external dependencies.

/*!
 * \internal
 * \ingroup TasmanianSG
 * \addtogroup TasmanianLinearSolvers Linear solvers
 *
 * \par Linear Solvers
 * Many sparse grids methods rely on various linear solvers, sparse iterative, eigenvalue,
 * fast-fourier-transform, and least-squares. The solvers are tuned for the specific
 * problems and can perform as good (or even better) than general third-party methods.
 * In addition, this reduces the external dependencies.
 */

namespace TasGrid{

//! \internal
//! \brief Methods for dense linear algebra.
//! \ingroup TasmanianLinearSolvers
namespace TasmanianDenseSolver{
    //! \internal
    //! \brief Least squares solver, used to infer anisotropic coefficients and thus rarely exceeds 10 - 20.
    //! \ingroup TasmanianLinearSolvers

    //! Matrix \b A is \b n by \b m (with \b n much larger), vector \b b has length \b n,
    //! the result is given in \b x (length \b m). The regularized term \b reg is added to the diagonal of
    //! the product `transpose(A) * A` to stabilize the Cholesky decomposition.
    void solveLeastSquares(int n, int m, const double A[], const double b[], double reg, double *x);
    //! \brief The same solver, but can take advantage of TPL.
    void solveLeastSquares(AccelerationContext const *acceleration, int n, int m, double A[], const double b[], double reg, double *x);
}

//! \internal
//! \brief Methods for tridiagonal eigenvalue problems.
//! \ingroup TasmanianLinearSolvers
namespace TasmanianTridiagonalSolver{
    //! \internal
    //! \brief Methods for tridiagonal eigenvalue problems, used to compute nodes and weights for Gaussian rules.
    //! \ingroup TasmanianLinearSolvers
    void decompose(int n, std::vector<double> &d, std::vector<double> &e, std::vector<double> &z);
}

//! \internal
//! \brief Methods for Fast-Fourier-Transform.
//! \ingroup TasmanianLinearSolvers
namespace TasmanianFourierTransform{
    //! \internal
    //! \brief Transfrom the data for a multi-dimensional tensor.
    //! \ingroup TasmanianLinearSolvers

    //! The \b num_points vector defines the number of points used by the tensor in different directions.
    //! The \b data holds the values for each point in the tensor, each point has a vector of values
    //! with size equal to the number of model outputs.
    //!
    //! The data is split into one-dimensional \a lines and each is tackled with \b fast_fourier_transform1D()
    void fast_fourier_transform(std::vector<std::vector<std::complex<double>>> &data, std::vector<int> &num_points);

    //! \internal
    //! \brief Perform one dimensional fast-fourier-transform (using radix-3).
    //! \ingroup TasmanianLinearSolvers

    //! Consider the line which is a subset of \b data defined by \b indexes and perform the radix-3 fast-fourier-transform.
    //! Called from \b fast_fourier_transform().
    void fast_fourier_transform1D(std::vector<std::vector<std::complex<double>>> &data, std::vector<int> &indexes);
}

//! \internal
//! \brief Methods for sparse linear algebra.
//! \ingroup TasmanianLinearSolvers
namespace TasSparse{

//! \internal
//! \brief Used to manipulate the wavelet values and solve for the wavelet coefficients.
//! \ingroup TasmanianLinearSolvers
class WaveletBasisMatrix{
public:
    //! \brief Default constructor, create an empty matrix.
    WaveletBasisMatrix() : tol(Maths::num_tol), num_rows(0){}
    //! \brief Initialize the matrix with the given set of indexes.
    WaveletBasisMatrix(AccelerationContext const *acceleration, const std::vector<int> &lpntr, const std::vector<std::vector<int>> &lindx, const std::vector<std::vector<double>> &lvals);
    //! \brief Default destructor.
    ~WaveletBasisMatrix() = default;

    //! \brief Return true if using the sparse mode and false is using dense mode or empty.
    bool isSparse() const{ return (num_rows > 0) and dense.empty(); }
    //! \brief Return true if using the dense mode and false is using sparse mode or empty.
    bool isDense() const{ return (num_rows > 0) and not dense.empty(); }

    //! \brief Return the number of rows in the matrix.
    int getNumRows() const{ return num_rows; }

    //! \brief Overwrites vector \b b with \b x that solves \f$ A^T x = b \f$.
    void invertTransposed(AccelerationContext const *acceleration, double b[]) const;

    //! \brief Overwrites the row-major matrix \b B with \b X that solves \f$ A X = B \f$.
    void invert(AccelerationContext const *acceleration, int num_colums, double B[]);

    //! \brief Solve `op(A) x = b` where `op` is either identity (find the coefficients) or transpose (find the interpolation weights).
    template<bool transpose, bool blas> void solve(const double b[], double x[]) const;

protected:
    //! \brief Compute the incomplete lower-upper decomposition of the matrix (zero extra fill).
    void computeILU();
    //! \brief Apply the preconditioner on a vector.
    template<bool transpose> void applyILU(double x[]) const;
    //! \brief Sets to be the product of the sparse matrix times x.
    template<bool transpose> void apply(double const x[], double r[]) const;
    //! \brief Computes the residual.
    void residual(double const x[], double const b[], double r[]) const;
    //! \brief Expressive call to transpose method.
    static constexpr bool use_transpose = true;
    //! \brief Expressive call to no transpose method.
    static constexpr bool no_transpose = false;
    //! \brief Expressive call to blas variant method.
    static constexpr bool use_blas = true;
    //! \brief Expressive call to no-blas variant method.
    static constexpr bool no_blas = false;

private:
    double tol;
    int num_rows;
    std::vector<int> pntr, indx, indxD;
    std::vector<double> vals, ilu;
    std::vector<double> dense; // if using the dense format
    std::vector<int> ipiv; // pivots for the dense factorize
};

}

}

#endif
