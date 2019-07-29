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

#ifndef __TASMANIAN_ADDONS_MPICONSTRUCTGRID_HPP
#define __TASMANIAN_ADDONS_MPICONSTRUCTGRID_HPP

/*!
 * \internal
 * \file tsgMPIConstructGrid.hpp
 * \brief Sparse Grids construction through sampling and MPI.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianAddonsCommon
 *
 * Templates that construct a sparse grid with sampling distributed across an MPI communicator.
 * \endinternal
 */

#include "tsgMPIScatterGrid.hpp"
#include "tsgConstructSurrogate.hpp"

/*!
 * \ingroup TasmanianAddons
 * \addtogroup TasmanianAddonsMPIConstruct MPI Automated Surrogate Construction Procedure
 *
 * Distributed sampling procedure, construction similar to TasGrid::constructSurrogate()
 * but the parallelism is associated with different ranks in an MPI communicator.
 * These templates require Tasmanian_ENABLE_MPI=ON.
 */

#ifdef Tasmanian_ENABLE_MPI

namespace TasGrid{

/*!
 * \internal
 * \ingroup TasmanianAddonsConstruct
 * \brief MPI construction algorithm using generic candidates procedure.
 *
 * Implements the MPI algorithm with generic refinement procedure given to TasGrid::constructCommon().
 * \endinternal
 */
template<bool use_initial_guess>
void mpiConstructCommon(std::function<void(std::vector<double> const &x, std::vector<double> &y)> model,
                        int num_dimensions, int num_outputs,
                        size_t max_num_points, size_t max_samples_per_job,
                        size_t max_num_ranks, int tagx, int tagy, int root, MPI_Comm comm,
                        TasmanianSparseGrid &grid,
                        std::function<std::vector<double>(TasmanianSparseGrid &)> candidates,
                        std::string const &checkpoint_filename){
    int num_ranks;
    MPI_Comm_size(comm, &num_ranks);
    // force max_num_ranks between 1 and num_ranks
    max_num_ranks = std::max(std::min(max_num_ranks, size_t(num_ranks)), size_t(1));
    int me; // my rank within the comm
    MPI_Comm_rank(comm, &me);

    std::vector<int> thread_to_rank(num_ranks); // map the thread id to the rank, root is used last
    std::iota(thread_to_rank.begin(), thread_to_rank.begin() + root, 0);
    std::iota(thread_to_rank.begin() + root, thread_to_rank.end(), root + 1);
    thread_to_rank.back() = root;

    std::mutex seq_threads;

    size_t max_message = sizeof(double) * Utils::size_mult(num_dimensions, max_samples_per_job) + sizeof(size_t);
    if (use_initial_guess == with_initial_guess)
        max_message += sizeof(double) * Utils::size_mult(num_outputs, max_samples_per_job) + sizeof(size_t);

    if (me == root){ // main sampler
        MPI_Barrier(comm);
        constructCommon<TasGrid::mode_parallel, use_initial_guess>(
            [&](std::vector<double> const &x, std::vector<double> &y, size_t thread_id)->void{
                int rankid = thread_to_rank[thread_id];
                if (rankid == root){
                    model(x, y); // thread belongs to this rank, do work locally
                }else{
                    std::stringstream ss; // pack x and y into a message
                    IO::writeNumbers<mode_binary, IO::pad_none>(ss, x.size());
                    IO::writeVector<mode_binary, IO::pad_none>(x, ss);
                    if (use_initial_guess == with_initial_guess){ // pack y only if used
                        IO::writeNumbers<mode_binary, IO::pad_none>(ss, y.size());
                        IO::writeVector<mode_binary, IO::pad_none>(y, ss);
                    }
                    MPI_Send(ss.str().c_str(), (int) (ss.str().size()), MPI_BYTE, rankid, tagx, comm); // send x
                    MPI_Status status;
                    y.resize(Utils::size_mult(num_outputs, x.size() / num_dimensions));
                    MPI_Recv(y.data(), (int) y.size(), MPI_DOUBLE, rankid, tagy, comm, &status); // receive y
                }
            }, max_num_points, max_num_ranks, max_samples_per_job, grid, candidates, checkpoint_filename);

        std::stringstream stop_message;
        IO::writeNumbers<mode_binary, IO::pad_none>(stop_message, size_t(0));
        if (use_initial_guess == with_initial_guess) IO::writeNumbers<mode_binary, IO::pad_none>(stop_message, size_t(0)); // pass 0 for y
        for(size_t id=0; id<max_num_ranks; id++){
            if (thread_to_rank[id] != root){
                MPI_Send(stop_message.str().c_str(), (int) (stop_message.str().size()), MPI_BYTE, thread_to_rank[id], tagx, comm); // send x
            }
        }
    }else{
        size_t thread_id = 0; // find the thread id that should communicate with my rank
        while(thread_to_rank[thread_id] != me) thread_id++;
        MPI_Barrier(comm);
        if (thread_id < max_num_ranks){ // if such thread would exist do work, otherwise do nothing
            std::vector<char> message_buffer(max_message);
            std::vector<double> x, y;
            do{
                MPI_Status status;
                MPI_Recv(message_buffer.data(), (int) message_buffer.size(), MPI_BYTE, root, tagx, comm, &status); // receive x

                // unpack the message
                VectorToStreamBuffer data_buffer(message_buffer); // do not modify buff after this point
                std::istream is(&data_buffer);
                x.resize(IO::readNumber<mode_binary, size_t>(is));
                IO::readVector<mode_binary>(is, x);
                if (use_initial_guess == with_initial_guess){ // unpack y, if necessary
                    y.resize(IO::readNumber<mode_binary, size_t>(is));
                    IO::readVector<mode_binary>(is, y);
                }else{
                    y.resize(Utils::size_mult(num_outputs, x.size() / num_dimensions));
                }

                if (!x.empty()){
                    model(x, y); // call the model
                    MPI_Send(y.data(), (int) (y.size()), MPI_DOUBLE, root, tagy, comm); // send y
                }
            }while(!x.empty());
        }
    }
}

template<bool use_initial_guess = no_initial_guess>
void mpiConstructSurrogate(std::function<void(std::vector<double> const &x, std::vector<double> &y)> model,
                           int num_dimensions, int num_outputs,
                           size_t max_num_points, size_t max_samples_per_job,
                           size_t max_num_ranks, int tagx, int tagy, int root, MPI_Comm comm,
                           TasmanianSparseGrid &grid,
                           double tolerance, TypeRefinement criteria, int output = -1,
                           std::vector<int> const &level_limits = std::vector<int>(),
                           std::vector<double> const &scale_correction = std::vector<double>(),
                           std::string const &checkpoint_filename = std::string()){

    mpiConstructCommon<use_initial_guess>(model, num_dimensions, num_outputs, max_num_points, max_samples_per_job,
                                          max_num_ranks, tagx, tagy, root, comm, grid,
                                          [&](TasmanianSparseGrid &g)->std::vector<double>{
                                              return g.getCandidateConstructionPoints(tolerance, criteria, output, level_limits, scale_correction);
                                          },
                                          checkpoint_filename);
}

template<bool use_initial_guess = no_initial_guess>
void mpiConstructSurrogate(std::function<void(std::vector<double> const &x, std::vector<double> &y)> model,
                           int num_dimensions, int num_outputs,
                           size_t max_num_points, size_t max_samples_per_job,
                           size_t max_num_ranks, int tagx, int tagy, int root, MPI_Comm comm,
                           TasmanianSparseGrid &grid,
                           TypeDepth type, std::vector<int> const &anisotropic_weights = std::vector<int>(),
                           std::vector<int> const &level_limits = std::vector<int>(),
                           std::string const &checkpoint_filename = std::string()){

    mpiConstructCommon<use_initial_guess>(model, num_dimensions, num_outputs, max_num_points, max_samples_per_job,
                                          max_num_ranks, tagx, tagy, root, comm, grid,
                                          [&](TasmanianSparseGrid &g)->std::vector<double>{
                                               return g.getCandidateConstructionPoints(type, anisotropic_weights, level_limits);
                                          },
                                          checkpoint_filename);
}

template<bool use_initial_guess = no_initial_guess>
void mpiConstructSurrogate(std::function<void(std::vector<double> const &x, std::vector<double> &y)> model,
                           int num_dimensions, int num_outputs,
                           size_t max_num_points, size_t max_samples_per_job,
                           size_t max_num_ranks, int tagx, int tagy, int root, MPI_Comm comm,
                           TasmanianSparseGrid &grid,
                           TypeDepth type, int output, std::vector<int> const &level_limits = std::vector<int>(),
                           std::string const &checkpoint_filename = std::string()){

    mpiConstructCommon<use_initial_guess>(model, num_dimensions, num_outputs, max_num_points, max_samples_per_job,
                                          max_num_ranks, tagx, tagy, root, comm, grid,
                                          [&](TasmanianSparseGrid &g)->std::vector<double>{
                                               return g.getCandidateConstructionPoints(type, output, level_limits);
                                           },
                                          checkpoint_filename);
}

}

#endif // Tasmanian_ENABLE_MPI

#endif
