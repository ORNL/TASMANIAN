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

SplitDirections::SplitDirections(const IndexSet *points) : num_dimensions(0), num_allocated_jobs(0), num_jobs(0), job_directions(0), num_job_pnts(0), job_pnts(0) {
    int num_points = points->getNumIndexes();
    num_dimensions = points->getNumDimensions();
    int *map = new int[num_points * num_dimensions];  std::fill(map, map + num_points * num_dimensions, 0);

    num_allocated_jobs = 256;
    job_directions = new int[num_allocated_jobs];
    num_job_pnts   = new int[num_allocated_jobs];
    job_pnts       = new int*[num_allocated_jobs];

    // This section is taking too long! Find a way to parallelize it!
    for(int s=0; s<num_points; s++){
        const int *p = points->getIndex(s);
        for(int d=0; d<num_dimensions; d++){
            if (map[s*num_dimensions + d] == 0){
                if (num_jobs == num_allocated_jobs){
                    num_allocated_jobs *= 2;
                    int *tmp = job_directions;
                    job_directions = new int[num_allocated_jobs];  std::copy(tmp, tmp + num_jobs, job_directions);
                    delete[] tmp;
                    tmp = num_job_pnts;
                    num_job_pnts = new int[num_allocated_jobs];  std::copy(tmp, tmp + num_jobs, num_job_pnts);
                    delete[] tmp;
                    int **ttmp = job_pnts;
                    job_pnts = new int*[num_allocated_jobs];  std::copy(ttmp, ttmp + num_jobs, job_pnts);
                    delete[] ttmp;
                }

                job_directions[num_jobs] = d;

                int nump_allocated = 128;
                job_pnts[num_jobs] = new int[nump_allocated];
                num_job_pnts[num_jobs] = 0;
                for(int i=0; i<num_points; i++){
                    if ((map[i*num_dimensions + d] == 0) && (doesBelongSameLine(p, points->getIndex(i), d))){
                        // adding a new point
                        if (num_job_pnts[num_jobs] == nump_allocated){
                            nump_allocated *= 2;
                            int *tmp = job_pnts[num_jobs];
                            job_pnts[num_jobs] = new int[nump_allocated];  std::copy(tmp, tmp + num_job_pnts[num_jobs], job_pnts[num_jobs]);
                            delete[] tmp;
                        }

                        job_pnts[num_jobs][num_job_pnts[num_jobs]] = i;
                        num_job_pnts[num_jobs]++;
                        map[i*num_dimensions + d] = -1;
                    }
                }

                num_jobs++;
            }
        }
    }

    delete[] map;
}
SplitDirections::~SplitDirections(){
    for(int s=0; s<num_jobs; s++){
        delete[] job_pnts[s];
    }
    delete[] job_pnts;
    delete[] num_job_pnts;
    delete[] job_directions;
}

int SplitDirections::getNumJobs() const{  return num_jobs;  }
int SplitDirections::getJobDirection(int job) const{  return job_directions[job];  }
int SplitDirections::getJobNumPoints(int job) const{  return num_job_pnts[job];  }
const int* SplitDirections::getJobPoints(int job) const{  return job_pnts[job];  }

bool SplitDirections::doesBelongSameLine(const int a[], const int b[], int direction) const{
    for(int i=0; i<num_dimensions; i++){
        if ((i != direction) && (a[i] != b[i])) return false;
    }
    return true;
}

}

#endif
