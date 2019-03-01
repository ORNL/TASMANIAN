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

#ifndef __TSG_BASE_CLASS_CPP
#define __TSG_BASE_CLASS_CPP

#include "tsgGridCore.hpp"

namespace TasGrid{

BaseCanonicalGrid::BaseCanonicalGrid(){}
BaseCanonicalGrid::~BaseCanonicalGrid(){}

SplitDirections::SplitDirections(const MultiIndexSet &points){
    // split the points into "jobs", where each job represents a batch of
    // points that lay on a line in some direction
    // int    job_directions[i] gives the direction of the i-th job
    // vector job_pnts[i] gives a list of the points (using indexes within the set)
    int num_points = points.getNumIndexes();
    num_dimensions = points.getNumDimensions();

    for(int d=0; d<num_dimensions; d++){
        // working with direction d
        // sort all points but ignore index d
        std::vector<int> map(num_points);
        std::iota(map.begin(), map.end(), 0);
        std::sort(map.begin(), map.end(), [&](int a, int b)->bool{
            const int * idxa = points.getIndex(a);
            const int * idxb = points.getIndex(b);
            // lexigographical order ignoring dimension d
            for(int j=0; j<num_dimensions; j++)
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
            job_directions.push_back(d);
            std::vector<int> pnts = {*imap++};
            // while the points are in the same direction as the reference, add to the same job
            while((imap != map.end()) && doesBelongSameLine(p, points.getIndex(*imap), d))
                pnts.push_back(*imap++);
            job_pnts.push_back(std::move(pnts));
        }
    }
}
SplitDirections::~SplitDirections(){}

int SplitDirections::getNumJobs() const{ return (int) job_pnts.size(); }
int SplitDirections::getJobDirection(int job) const{ return job_directions[job]; }
int SplitDirections::getJobNumPoints(int job) const{ return (int) job_pnts[job].size(); }
const int* SplitDirections::getJobPoints(int job) const{ return job_pnts[job].data(); }

bool SplitDirections::doesBelongSameLine(const int a[], const int b[], int direction) const{
    for(int i=0; i<num_dimensions; i++){
        if ((i != direction) && (a[i] != b[i])) return false;
    }
    return true;
}

}

#endif
