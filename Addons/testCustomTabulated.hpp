/*
 * Copyright (c) 2021, Miroslav Stoyanov & William Kong
 *
 * This file is part of
 * Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following
 * conditions are met:
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
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND
 * IMPLIED. THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE
 * OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL
 * ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE. THE USER ASSUMES
 * RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING
 * FROM OR ARISING OUT OF, IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
 */


#include "TasmanianAddons.hpp"

// Test getSubrules().
inline bool testGetSubrules() {
    bool passed = true;
    // Level == 0
    TasGrid::CustomTabulated empty_ct;
    TasGrid::CustomTabulated empty_sub_ct = TasGrid::getSubrules(empty_ct, 0, 2, "Empty rules");
    if (empty_sub_ct.getNumLevels() != 0) {
        passed = false;
    }
    // Level >= 1
    std::vector<int> n_vec = {1, 2, 10};
    std::vector<int> start_index_vec = {0, 1, 4};
    std::vector<int> stride_vec = {1, 2, 3};
    for (auto n : n_vec) {
        // Create the full CustomTabulated instance with n levels.
        std::vector<int> num_nodes(n), precision(n);
        std::vector<std::vector<double>> nodes(n), weights(n);
        for (int i=0; i<n; i++) {
            num_nodes[i] = i + 1;
            precision[i] = i;
            OneDimensionalNodes::getGaussLegendre(num_nodes[i], weights[i], nodes[i]);
        }
        TasGrid::CustomTabulated ct = TasGrid::CustomTabulated(std::move(num_nodes), std::move(precision), std::move(nodes),
                                                               std::move(weights), "Gauss-Legendre rules");
        // Create the subset and test different start_index and stride combinations.
        for (auto start_index : start_index_vec) {
            for (auto stride : stride_vec) {
                std::string descr = "subset_" + std::to_string(start_index) + "_" + std::to_string(stride); 
                TasGrid::CustomTabulated sub_ct = TasGrid::getSubrules(ct, start_index, stride, descr);
                if (sub_ct.getDescription() != descr) {
                    passed = false;
                }
                int sub_level = 0;
                for (int level=start_index; level<ct.getNumLevels(); level+=stride) {
                    if (ct.getNumPoints(level) != sub_ct.getNumPoints(sub_level)) {
                        passed = false;
                    }
                    if (ct.getQExact(level) != sub_ct.getQExact(sub_level)) {
                        passed = false;
                    }
                    std::vector<double> w, x, sub_w, sub_x;
                    ct.getWeightsNodes(level, w, x);
                    sub_ct.getWeightsNodes(sub_level, sub_w, sub_x);
                    for (int i=0; i < ct.getNumPoints(level); i++) {
                        if (w[i] != sub_w[i] || x[i] != sub_x[i]) {
                            passed = false;
                        }
                    }
                    sub_level++;
                }
                if (sub_ct.getNumLevels() != sub_level) {
                    passed = false;
                }
            }
        }
    }
    return passed;
}

// Main wrapper
inline bool testCustomTabulated() {
    bool passed = true;
    passed = passed && testGetSubrules() ;
    return passed;
}
