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
#ifndef __TASMANIAN_ADDONS_CCONSTRUCTSURROGATE_CPP
#define __TASMANIAN_ADDONS_CCONSTRUCTSURROGATE_CPP

#include "TasmanianAddons.hpp"

extern "C"{

// num_samples, num_dimensions, x.data(), num_outputs, y.data(), thread_id, error_code
using tsg_scs_model = void (*)(int, int, double const*, int, double*, int, int*);

void tsgConstructSurrogateNoIGSurplus
           (tsg_scs_model pymodel,
            int max_num_points, int num_parallel_jobs, int max_samples_per_job,
            void *grid_pntr,
            double tolerance, const char* s_criteria, int output,
            int *llimits, const char *checkpoint_filename, int *err){

    *err = 1;
    TasGrid::TasmanianSparseGrid &grid = *reinterpret_cast<TasGrid::TasmanianSparseGrid*>(grid_pntr);

    int const num_dimensions = grid.getNumDimensions();
    int const num_outputs    = grid.getNumOutputs();

    TasGrid::TypeRefinement criteria = TasGrid::IO::getTypeRefinementString(s_criteria);

    std::vector<int> level_limits = TasGrid::Utils::copyArray(llimits, num_dimensions);
    std::string cfname = (checkpoint_filename != nullptr) ? std::string(checkpoint_filename) : std::string();

    auto cpp_model = [&](std::vector<double> const &x, std::vector<double> &y, size_t thread_id)->
        void{
            int sample_size = (int) x.size() / num_dimensions;
            int error_code = 0;
            pymodel(sample_size, num_dimensions, x.data(), num_outputs, y.data(), (int) thread_id, &error_code);
            if (error_code != 0) throw std::runtime_error("The Python callback returned an error in tsgConstructSurrogateNoIGSurplus()");
        };

    try{
        if (num_parallel_jobs > 1){
            TasGrid::constructSurrogate<TasGrid::mode_parallel, TasGrid::no_initial_guess>
                (cpp_model, max_num_points, num_parallel_jobs, max_samples_per_job,
                grid, tolerance, criteria, output, level_limits, cfname);
        }else{
            TasGrid::constructSurrogate<TasGrid::mode_sequential, TasGrid::no_initial_guess>
                (cpp_model, max_num_points, num_parallel_jobs, max_samples_per_job,
                grid, tolerance, criteria, output, level_limits, cfname);
        }
        *err = 0; // success
    }catch(std::runtime_error &){} // *err will remain 1
}

void tsgConstructSurrogateNoIGAniso
           (tsg_scs_model pymodel,
            int max_num_points, int num_parallel_jobs, int max_samples_per_job, void *grid_pntr,
            const char* s_type, int output, int *llimits, const char *checkpoint_filename, int *err){

    *err = 1;
    TasGrid::TasmanianSparseGrid &grid = *reinterpret_cast<TasGrid::TasmanianSparseGrid*>(grid_pntr);

    int const num_dimensions = grid.getNumDimensions();
    int const num_outputs    = grid.getNumOutputs();

    TasGrid::TypeDepth dtype = TasGrid::IO::getDepthTypeString(s_type);

    std::vector<int> level_limits = TasGrid::Utils::copyArray(llimits, num_dimensions);
    std::string cfname = (checkpoint_filename != nullptr) ? std::string(checkpoint_filename) : std::string();

    auto cpp_model = [&](std::vector<double> const &x, std::vector<double> &y, size_t thread_id)->
        void{
            int sample_size = (int) x.size() / num_dimensions;
            int error_code = 0;
            pymodel(sample_size, num_dimensions, x.data(), num_outputs, y.data(), (int) thread_id, &error_code);
            if (error_code != 0) throw std::runtime_error("The Python callback returned an error in tsgConstructSurrogateNoIGAniso()");
        };

    try{
        if (num_parallel_jobs > 1){
            TasGrid::constructSurrogate<TasGrid::mode_parallel, TasGrid::no_initial_guess>
                (cpp_model, max_num_points, num_parallel_jobs, max_samples_per_job,
                grid, dtype, output, level_limits, cfname);
        }else{
            TasGrid::constructSurrogate<TasGrid::mode_sequential, TasGrid::no_initial_guess>
                (cpp_model, max_num_points, num_parallel_jobs, max_samples_per_job,
                grid, dtype, output, level_limits, cfname);
        }
        *err = 0; // success
    }catch(std::runtime_error &){} // *err will remain 1
}

void tsgConstructSurrogateNoIGAnisoFixed
           (tsg_scs_model pymodel,
            int max_num_points, int num_parallel_jobs, int max_samples_per_job, void *grid_pntr,
            const char* s_type, int *aweights, int *llimits, const char *checkpoint_filename, int *err){

    *err = 1;
    TasGrid::TasmanianSparseGrid &grid = *reinterpret_cast<TasGrid::TasmanianSparseGrid*>(grid_pntr);

    int const num_dimensions = grid.getNumDimensions();
    int const num_outputs    = grid.getNumOutputs();

    TasGrid::TypeDepth dtype = TasGrid::IO::getDepthTypeString(s_type);

    std::vector<int> anisotropic_weights = TasGrid::Utils::copyArray(aweights, num_dimensions *
                                        ((TasGrid::OneDimensionalMeta::isTypeCurved(dtype)) ? 2 : 1));
    std::vector<int> level_limits = TasGrid::Utils::copyArray(llimits, num_dimensions);
    std::string cfname = (checkpoint_filename != nullptr) ? std::string(checkpoint_filename) : std::string();

    auto cpp_model = [&](std::vector<double> const &x, std::vector<double> &y, size_t thread_id)->
        void{
            int sample_size = (int) x.size() / num_dimensions;
            int error_code = 0;
            pymodel(sample_size, num_dimensions, x.data(), num_outputs, y.data(), (int) thread_id, &error_code);
            if (error_code != 0) throw std::runtime_error("The Python callback returned an error in tsgConstructSurrogateNoIGAnisoFixed()");
        };

    try{
        if (num_parallel_jobs > 1){
            TasGrid::constructSurrogate<TasGrid::mode_parallel, TasGrid::no_initial_guess>
                (cpp_model, max_num_points, num_parallel_jobs, max_samples_per_job,
                grid, dtype, anisotropic_weights, level_limits, cfname);
        }else{
            TasGrid::constructSurrogate<TasGrid::mode_sequential, TasGrid::no_initial_guess>
                (cpp_model, max_num_points, num_parallel_jobs, max_samples_per_job,
                grid, dtype, anisotropic_weights, level_limits, cfname);
        }
        *err = 0; // success
    }catch(std::runtime_error &){} // *err will remain 1
}

// num_samples, num_dimensions, x.data(), has_init_guess, num_outputs, y.data(), thread_id, error_code
using tsg_ics_model = void (*)(int, int, double const*, int, int, double*, int, int*);

void tsgConstructSurrogateWiIGSurplus
           (tsg_ics_model pymodel,
            int max_num_points, int num_parallel_jobs, int max_samples_per_job,
            void *grid_pntr,
            double tolerance, const char* s_criteria, int output,
            int *llimits, const char *checkpoint_filename, int *err){

    *err = 1;
    TasGrid::TasmanianSparseGrid &grid = *reinterpret_cast<TasGrid::TasmanianSparseGrid*>(grid_pntr);

    int const num_dimensions = grid.getNumDimensions();
    int const num_outputs    = grid.getNumOutputs();

    TasGrid::TypeRefinement criteria = TasGrid::IO::getTypeRefinementString(s_criteria);

    std::vector<int> level_limits = TasGrid::Utils::copyArray(llimits, num_dimensions);
    std::string cfname = (checkpoint_filename != nullptr) ? std::string(checkpoint_filename) : std::string();

    auto cpp_model = [&](std::vector<double> const &x, std::vector<double> &y, size_t thread_id)->
        void{
            int sample_size = (int) x.size() / num_dimensions;

            int ihas_guess = 1; // assume there is a guess
            if (y.empty()){
                ihas_guess = 0; // no guess
                y.resize(TasGrid::Utils::size_mult(sample_size, num_outputs));
            }

            int error_code = 0;
            pymodel(sample_size, num_dimensions, x.data(), ihas_guess, num_outputs, y.data(), (int) thread_id, &error_code);
            if (error_code != 0) throw std::runtime_error("The Python callback returned an error in tsgConstructSurrogateWiIGSurplus()");
        };

    try{
        if (num_parallel_jobs > 1){
            TasGrid::constructSurrogate<TasGrid::mode_parallel, TasGrid::with_initial_guess>
                (cpp_model, max_num_points, num_parallel_jobs, max_samples_per_job,
                grid, tolerance, criteria, output, level_limits, cfname);
        }else{
            TasGrid::constructSurrogate<TasGrid::mode_sequential, TasGrid::with_initial_guess>
                (cpp_model, max_num_points, num_parallel_jobs, max_samples_per_job,
                grid, tolerance, criteria, output, level_limits, cfname);
        }
        *err = 0; // success
    }catch(std::runtime_error &){} // *err will remain 1
}

void tsgConstructSurrogateWiIGAniso
           (tsg_ics_model pymodel,
            int max_num_points, int num_parallel_jobs, int max_samples_per_job, void *grid_pntr,
            const char* s_type, int output, int *llimits, const char *checkpoint_filename, int *err){

    *err = 1;
    TasGrid::TasmanianSparseGrid &grid = *reinterpret_cast<TasGrid::TasmanianSparseGrid*>(grid_pntr);

    int const num_dimensions = grid.getNumDimensions();
    int const num_outputs    = grid.getNumOutputs();

    TasGrid::TypeDepth dtype = TasGrid::IO::getDepthTypeString(s_type);

    std::vector<int> level_limits = TasGrid::Utils::copyArray(llimits, num_dimensions);
    std::string cfname = (checkpoint_filename != nullptr) ? std::string(checkpoint_filename) : std::string();

    auto cpp_model = [&](std::vector<double> const &x, std::vector<double> &y, size_t thread_id)->
        void{
            int sample_size = (int) x.size() / num_dimensions;

            int ihas_guess = 1; // assume there is a guess
            if (y.empty()){
                ihas_guess = 0; // no guess
                y.resize(TasGrid::Utils::size_mult(sample_size, num_outputs));
            }

            int error_code = 0;
            pymodel(sample_size, num_dimensions, x.data(), ihas_guess, num_outputs, y.data(), (int) thread_id, &error_code);
            if (error_code != 0) throw std::runtime_error("The Python callback returned an error in tsgConstructSurrogateWiIGAniso()");
        };

    try{
        if (num_parallel_jobs > 1){
            TasGrid::constructSurrogate<TasGrid::mode_parallel, TasGrid::with_initial_guess>
                (cpp_model, max_num_points, num_parallel_jobs, max_samples_per_job,
                grid, dtype, output, level_limits, cfname);
        }else{
            TasGrid::constructSurrogate<TasGrid::mode_sequential, TasGrid::with_initial_guess>
                (cpp_model, max_num_points, num_parallel_jobs, max_samples_per_job,
                grid, dtype, output, level_limits, cfname);
        }
        *err = 0; // success
    }catch(std::runtime_error &){} // *err will remain 1
}

void tsgConstructSurrogateWiIGAnisoFixed
           (tsg_ics_model pymodel,
            int max_num_points, int num_parallel_jobs, int max_samples_per_job, void *grid_pntr,
            const char* s_type, int *aweights, int *llimits, const char *checkpoint_filename, int *err){

    *err = 1;
    TasGrid::TasmanianSparseGrid &grid = *reinterpret_cast<TasGrid::TasmanianSparseGrid*>(grid_pntr);

    int const num_dimensions = grid.getNumDimensions();
    int const num_outputs    = grid.getNumOutputs();

    TasGrid::TypeDepth dtype = TasGrid::IO::getDepthTypeString(s_type);

    std::vector<int> anisotropic_weights = TasGrid::Utils::copyArray(aweights, num_dimensions *
                                        ((TasGrid::OneDimensionalMeta::isTypeCurved(dtype)) ? 2 : 1));
    std::vector<int> level_limits = TasGrid::Utils::copyArray(llimits, num_dimensions);
    std::string cfname = (checkpoint_filename != nullptr) ? std::string(checkpoint_filename) : std::string();

    auto cpp_model = [&](std::vector<double> const &x, std::vector<double> &y, size_t thread_id)->
        void{
            int sample_size = (int) x.size() / num_dimensions;

            int ihas_guess = 1; // assume there is a guess
            if (y.empty()){
                ihas_guess = 0; // no guess
                y.resize(TasGrid::Utils::size_mult(sample_size, num_outputs));
            }

            int error_code = 0;
            pymodel(sample_size, num_dimensions, x.data(), ihas_guess, num_outputs, y.data(), (int) thread_id, &error_code);
            if (error_code != 0) throw std::runtime_error("The Python callback returned an error in tsgConstructSurrogateWiIGAnisoFixed()");
        };

    try{
        if (num_parallel_jobs > 1){
            TasGrid::constructSurrogate<TasGrid::mode_parallel, TasGrid::with_initial_guess>
                (cpp_model, max_num_points, num_parallel_jobs, max_samples_per_job,
                grid, dtype, anisotropic_weights, level_limits, cfname);
        }else{
            TasGrid::constructSurrogate<TasGrid::mode_sequential, TasGrid::with_initial_guess>
                (cpp_model, max_num_points, num_parallel_jobs, max_samples_per_job,
                grid, dtype, anisotropic_weights, level_limits, cfname);
        }
        *err = 0; // success
    }catch(std::runtime_error &){} // *err will remain 1
}

} // extern "C"
#endif
