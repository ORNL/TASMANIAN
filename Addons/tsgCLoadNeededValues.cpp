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
#ifndef __TASMANIAN_ADDONS_CLOADNEEDEDVALS_CPP
#define __TASMANIAN_ADDONS_CLOADNEEDEDVALS_CPP

#include "TasmanianAddons.hpp"

extern "C"{

// x.size(), x.data(), y.size(), y.data(), thread_id
using tsg_lnp_model = void (*)(int, double const*, int, double*, int, int*);

void tsgLoadNeededValues(int overwrite, tsg_lnp_model pymodel, void *grid_pntr, int num_threads, int *err){
    *err = 1; // produce an error code, if not ended early
    try{
        TasGrid::TasmanianSparseGrid &grid = *reinterpret_cast<TasGrid::TasmanianSparseGrid*>(grid_pntr);

        int const num_dimensions = grid.getNumDimensions();
        int const num_outputs    = grid.getNumOutputs();

        auto cpp_model = [&](double const x[], double y[], size_t thread_id)->
            void{
                int error_code = 0;
                pymodel(num_dimensions, x, num_outputs, y, (int) thread_id, &error_code);
                if (error_code != 0) throw std::runtime_error("The Python callback returned an error in tsgLoadNeededValues()");
            };

        constexpr bool needed = true;
        constexpr bool loaded = false;

        if (num_threads > 1){
            if (overwrite != 0){
                TasGrid::loadNeededValues<TasGrid::mode_parallel, needed>(cpp_model, grid, num_threads);
            }else{
                TasGrid::loadNeededValues<TasGrid::mode_parallel, loaded>(cpp_model, grid, num_threads);
            }
        }else{
            if (overwrite != 0){
                TasGrid::loadNeededValues<TasGrid::mode_sequential, needed>(cpp_model, grid, 1);
            }else{
                TasGrid::loadNeededValues<TasGrid::mode_sequential, loaded>(cpp_model, grid, 1);
            }
        }
        *err = 0; // got here, no exceptions were encountered
    }catch(std::runtime_error &){} // keep *err set to 1
}

} // extern "C"
#endif
