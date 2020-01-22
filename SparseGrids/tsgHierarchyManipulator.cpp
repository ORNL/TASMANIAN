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

#ifndef __TSG_HIERARCHY_MANIPULATOR_CPP
#define __TSG_HIERARCHY_MANIPULATOR_CPP

#include "tsgHierarchyManipulator.hpp"

namespace TasGrid{

namespace HierarchyManipulations{

Data2D<int> computeDAGup(MultiIndexSet const &mset, const BaseRuleLocalPolynomial *rule){
    size_t num_dimensions = mset.getNumDimensions();
    int num_points = mset.getNumIndexes();
    if (rule->getMaxNumParents() > 1){ // allow for multiple parents and level 0 may have more than one node
        int max_parents = rule->getMaxNumParents() * (int) num_dimensions;
        Data2D<int> parents(max_parents, num_points, -1);
        int level0_offset = rule->getNumPoints(0);
        #pragma omp parallel for schedule(static)
        for(int i=0; i<num_points; i++){
            const int *p = mset.getIndex(i);
            std::vector<int> dad(num_dimensions);
            std::copy_n(p, num_dimensions, dad.data());
            int *pp = parents.getStrip(i);
            for(size_t j=0; j<num_dimensions; j++){
                if (dad[j] >= level0_offset){
                    int current = p[j];
                    dad[j] = rule->getParent(current);
                    pp[2*j] = mset.getSlot(dad);
                    while ((dad[j] >= level0_offset) && (pp[2*j] == -1)){
                        current = dad[j];
                        dad[j] = rule->getParent(current);
                        pp[2*j] = mset.getSlot(dad);
                    }
                    dad[j] = rule->getStepParent(current);
                    if (dad[j] != -1){
                        pp[2*j + 1] = mset.getSlot(dad);
                    }
                    dad[j] = p[j];
                }
            }
        }
        return parents;
    }else{ // this assumes that level zero has only one node
        Data2D<int> parents((int) num_dimensions, num_points);
        #pragma omp parallel for schedule(static)
        for(int i=0; i<num_points; i++){
            const int *p = mset.getIndex(i);
            std::vector<int> dad(num_dimensions);
            std::copy_n(p, num_dimensions, dad.data());
            int *pp = parents.getStrip(i);
            for(size_t j=0; j<num_dimensions; j++){
                if (dad[j] == 0){
                    pp[j] = -1;
                }else{
                    dad[j] = rule->getParent(dad[j]);
                    pp[j] = mset.getSlot(dad.data());
                    while((dad[j] != 0) && (pp[j] == -1)){
                        dad[j] = rule->getParent(dad[j]);
                        pp[j] = mset.getSlot(dad);
                    }
                    dad[j] = p[j];
                }
            }
        }
        return parents;
    }
}

Data2D<int> computeDAGDown(MultiIndexSet const &mset, const BaseRuleLocalPolynomial *rule){
    size_t num_dimensions = mset.getNumDimensions();
    int max_1d_kids = rule->getMaxNumKids();
    int num_points = mset.getNumIndexes();
    Data2D<int> kids(Utils::size_mult(max_1d_kids, num_dimensions), num_points);

    #pragma omp parallel for
    for(int i=0; i<num_points; i++){
        std::vector<int> kid(num_dimensions);
        std::copy_n(mset.getIndex(i), num_dimensions, kid.data());
        auto family = kids.getIStrip(i);

        for(size_t j=0; j<num_dimensions; j++){
            int current = kid[j];
            for(int k=0; k<max_1d_kids; k++){
                kid[j] = rule->getKid(current, k);
                *family++ = (kid[j] == -1) ? -1 : mset.getSlot(kid);
            }
            kid[j] = current;
        }
    }

    return kids;
}

std::vector<int> computeLevels(MultiIndexSet const &mset, BaseRuleLocalPolynomial const *rule){
    size_t num_dimensions = mset.getNumDimensions();
    int num_points = mset.getNumIndexes();
    std::vector<int> level((size_t) num_points);
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = mset.getIndex(i);
        int current_level = rule->getLevel(p[0]);
        for(size_t j=1; j<num_dimensions; j++){
            current_level += rule->getLevel(p[j]);
        }
        level[i] = current_level;
    }
    return level;
}

void completeToLower(MultiIndexSet const &mset, MultiIndexSet &refined, BaseRuleLocalPolynomial const *rule){
    size_t num_dimensions = mset.getNumDimensions();
    size_t num_added = 1; // set to 1 to start the loop
    while(num_added > 0){
        Data2D<int> addons(num_dimensions, 0);
        int num_points = refined.getNumIndexes();

        for(int i=0; i<num_points; i++){
            std::vector<int> parent(refined.getIndex(i), refined.getIndex(i) + num_dimensions);
            for(auto &p : parent){
                int r = p;
                p = rule->getParent(r);
                if ((p != -1) && refined.missing(parent) && mset.missing(parent))
                    addons.appendStrip(parent);
                p = rule->getStepParent(r);
                if ((p != -1) && refined.missing(parent) && mset.missing(parent))
                    addons.appendStrip(parent);
                p = r;
            }
        }

        num_added = addons.getNumStrips();
        if (num_added > 0) refined += addons;
    }
}

SplitDirections::SplitDirections(const MultiIndexSet &points){
    // split the points into "jobs", where each job represents a batch of
    // points that lay on a line in some direction
    // int    job_directions[i] gives the direction of the i-th job
    // vector job_pnts[i] gives a list of the points (using indexes within the set)
    int num_points = points.getNumIndexes();
    size_t num_dimensions = points.getNumDimensions();

    auto doesBelongSameLine = [&](const int a[], const int b[], size_t direction)->
                                    bool{
                                        for(size_t i=0; i<num_dimensions; i++)
                                            if ((i != direction) && (a[i] != b[i])) return false;
                                        return true;
                                    };

    for(size_t d=0; d<num_dimensions; d++){
        // working with direction d
        // sort all points but ignore index d
        std::vector<int> map(num_points);
        std::iota(map.begin(), map.end(), 0);
        std::sort(map.begin(), map.end(), [&](int a, int b)->bool{
            const int * idxa = points.getIndex(a);
            const int * idxb = points.getIndex(b);
            // lexigographical order ignoring dimension d
            for(size_t j=0; j<num_dimensions; j++)
                if (j != d){
                    if (idxa[j] < idxb[j]) return true;
                    if (idxa[j] > idxb[j]) return false;
                }
            return false;
        });

        auto imap = map.begin();
        while(imap != map.end()){
            // new job, get reference index
            const int *p = points.getIndex(*imap);
            job_directions.push_back((int) d);
            job_pnts.emplace_back(std::vector<int>(1, *imap++));
            // while the points are in the same direction as the reference, add to the same job
            while((imap != map.end()) && doesBelongSameLine(p, points.getIndex(*imap), d))
                job_pnts.back().push_back(*imap++);
        }
    }
}

} // HierarchyManipulations

} // TasGrid

#endif
