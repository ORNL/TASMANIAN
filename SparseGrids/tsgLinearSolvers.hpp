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

#include <list>
#include <vector>
#include <fstream>
#include <math.h>

#include "tsgEnumerates.hpp"

namespace TasGrid{

namespace TasmanianDenseSolver{
    void solveLeastSquares(int n, int m, const double A[], const double b[], double reg, double *x);
}

namespace TasmanianTridiagonalSolver{
    void decompose(int n, double d[], double e[], double z[]);
}

}

// Sparse linear algebra is used by the wavelet grids only
// tsgEnumerates.hpp is needed by TasSparse for the hardcoded constants, i.e., TSG_DROP_TOL
namespace TasSparse{

using std::endl;

// the enumerates defined here are used only by the sparse code
enum TsgSparseSorting{
	unsorted,
	csc_sort,
	csr_sort
};

enum TsgSparseType{
	coo,
	csc,
	csr,
	none
};

enum TsgCgStatus{
	converged,
	max_iter_reached
};

// (I,J,V) triplet used to represent a sparse matrix entry where A[i,j] += v.
typedef struct COO_Triplet {
	unsigned int i;
	unsigned int j;
	double v;
	inline bool operator==(COO_Triplet &other) const{
		return (i == other.i) && (j == other.j);
	}
} COO_Triplet;

// Sparse matrix represented by above triplets.
class TsgSparseCOO {
public:
	TsgSparseCOO();
	TsgSparseCOO(int in_m, int in_n);
	~TsgSparseCOO();

	void addPoint(unsigned int i, unsigned int j, double v);
	void clear();
	int getNumRows();
	int getNumCols();
	void combine(TsgSparseCOO &m);
	friend class TsgSparseCSC;
	friend class TsgSparseCSR;
	friend class SparseMatrix;

protected:
	std::list<COO_Triplet> values;
	static bool csrComp(COO_Triplet a, COO_Triplet b);
	static bool cscComp(COO_Triplet a, COO_Triplet b);
	void sort(TsgSparseSorting type);
	void coalesce();
	int m, n, nnz;
private:
	TsgSparseSorting sorted;

}; // End TsgSparseCOO

// Make a new matrix class that would store the ILU preconditioner and solve either the regular or adjoined problem all in one class
class SparseMatrix{
public:
    SparseMatrix(TsgSparseCOO &M);
    ~SparseMatrix();

    int getNumRows() const;

	void solve(const double b[], double x[], bool transposed = false) const;

protected:
    void clear();

    void computeILU();

private:
    const double tol;
    int num_rows;
    int *pntr, *indx, *indxD;
    double *vals, *ilu;
};

}

#endif
