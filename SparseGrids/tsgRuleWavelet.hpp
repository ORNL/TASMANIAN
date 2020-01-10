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

#ifndef __TASMANIAN_SPARSE_GRID_WAVELET_RULE_HPP
#define __TASMANIAN_SPARSE_GRID_WAVELET_RULE_HPP

#include "tsgGridCore.hpp"

namespace TasGrid{
// These macros are used in accessing coarse and fine level coefficients for the cascade algorithm.

#ifndef __TASMANIAN_DOXYGEN_SKIP
class RuleWavelet{
public:
    RuleWavelet(int corder, int iter_depth) : order(0), iteration_depth(iter_depth), num_data_points((1 << iteration_depth) + 1){
        updateOrder(corder);
    }
    ~RuleWavelet() = default;

    // Interface to Extend
    int getNumPoints(int level) const; // get the number of points associated with level (also the index of the first point on the level+1)

    const char* getDescription() const;

    double getNode(int point) const; // returns the x-value of a point
    int getOrder() const;
    void updateOrder(int new_order); // Sets the order of the underlying wavelet rule. Involves recalculating approximations if order==3.

    double getWeight(int point) const; // get the quadrature weight associated with the point
    double eval(int point, double x) const; // returns the value of point at location x (there is assumed 1-1 correspondence between points and functions)

    int getLevel(int point) const; // returns the hierarchical level of a point
    void getChildren(int point, int &first, int &second) const; // Given a point, return the children (if any)
    int getParent(int point) const; // Returns the parent of the given node

    void getShiftScale(int point, double &scale, double &shift) const; // encodes for GPU purposes

    double getSupport(int point) const; // return the support of the point, for reporting purposes (not used in eval)
protected:
    inline double eval_linear(int pt, double x) const;
    inline double eval_cubic(int pt, double x) const;
    inline double linear_boundary_wavelet(double x) const;
    inline double linear_central_wavelet(double x) const;
    int order;
    int iteration_depth;
    int num_data_points;
    static void cubic_cascade(double *y, int starting_level, int iteration_depth);

    inline int find_index(double x) const;
    inline double interpolate(const double *y, double x) const;

    std::vector<std::vector<double>> data;
    std::vector<double> cachexs;
};
#endif // __TASMANIAN_DOXYGEN_SKIP

} // namespace TasGrid
#endif // __TASMANIAN_SPARSE_GRID_WAVELET_RULE_HPP
