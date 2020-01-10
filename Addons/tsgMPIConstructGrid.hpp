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

#include "tsgMPIScatterDream.hpp"
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
 * \ingroup TasmanianAddonsMPIConstruct
 * \brief Signature of a model function to be used in the MPI construction procedures.
 *
 * The logic is identical to the TasGrid::ModelSignature used in the TasGrid::constructSurrogate()
 * templates, with the exception of the removal of the thread-id.
 * In an MPI context, the thread-id is replaced by the MPI rank which can be queried using MPI_Comm_rank().
 */
using ModelSignatureMPI = std::function<void(std::vector<double> const &x, std::vector<double> &y)>;

/*!
 * \internal
 * \ingroup TasmanianAddonsMPIConstruct
 * \brief MPI construction algorithm using generic candidates procedure.
 *
 * Implements the MPI algorithm with generic refinement procedure given to TasGrid::constructCommon().
 * \endinternal
 */
template<bool use_initial_guess>
void mpiConstructCommon(ModelSignatureMPI model,
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
    int const me = getMPIRank(comm); // my rank within the comm

    std::vector<int> thread_to_rank(num_ranks); // map the thread id to the rank, root is used last
    std::iota(thread_to_rank.begin(), thread_to_rank.begin() + root, 0);
    std::iota(thread_to_rank.begin() + root, thread_to_rank.end(), root + 1);
    thread_to_rank.back() = root;

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
                    MPI_Request request;
                    y.resize(Utils::size_mult(num_outputs, x.size() / num_dimensions));
                    MPI_Irecv(y.data(), (int) y.size(), MPI_DOUBLE, rankid, tagy, comm, &request); // receive y
                    MPI_Wait(&request, &status); // wait for the result from the model
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

                x = IO::readVector<IO::mode_binary_type, double>(is, IO::readNumber<IO::mode_binary_type, size_t>(is));

                if (use_initial_guess == with_initial_guess){ // unpack y, if necessary
                    y = IO::readVector<IO::mode_binary_type, double>(is, IO::readNumber<IO::mode_binary_type, size_t>(is));
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

/*!
 * \ingroup TasmanianAddonsMPIConstruct
 * \brief Construct a sparse grid surrogate to the model defined by the lambda, MPI version.
 *
 * The logic is the same as in TasGrid::constructSurrogate(), construct adaptive sparse grid
 * surrogate to the mode defined by the lambda using parallel computations.
 * The model will be executed with different inputs \b x on different ranks within the MPI communicator.
 * All ranks of across the communicator must make an identical call with identical parameters
 * before the \b grid. The grid and the refinement and checkpoint parameters will be used
 * only by the \b root rank and will not be addressed by the rest.
 *
 * \par MPI and Multi-Threading
 * This template will call TasGrid::constructSurrogate() and it will spawn a separate thread
 * for each active ranks within the communicator. Thus, MPI must support multi-threading
 * and the multi-threaded initialization call must be used:
 * \code
 *   // instead of calling MPI_Init() call the following
 *   int threads_available;
 *   MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &threads_available);
 *   if (threads_available != MPI_THREAD_MULTIPLE)
 *      std::cerr << "MPI threading not available, TasGrid::mpiConstructSurrogate() will not work.\n";
 * \endcode
 * Note that the MPIGridSend() and other similar methods do not require multi-threading.
 *
 * The list of inputs is very similar to TasGrid::constructSurrogate(),
 * here we outline the main differences.
 *
 * \tparam use_initial_guess defines whether to compute an initial guess for the model outputs
 *          and send that over the communicator. The MPI procedure always runs in parallel mode
 *          hence there is no additional template parameter.
 *
 * \param model the inputs \b x and outputs \b y of the model are the same as in
 *      TasGrid::constructSurrogate(); however, there is no thread_id.
 *      Only one call of the model will be executed on each active rank and thus
 *      the MPI native MPI_Comm_rank() can be used in place of the thread_id.
 *
 * \param num_dimensions must match \b grid.getNumDimensions().
 *      This is a required parameter since only the \b root rank is assumed to have
 *      access to the grid, the other ranks need to know the correct dimension without the grid.
 * \param num_outputs must match \b grid.getNumOutputs(), see \b num_dimensions.
 *
 * \param max_num_points defines the computational budget, same as TasGrid::constructSurrogate().
 *
 * \param max_samples_per_job defines the number of samples that will be given to the model
 *      in each call, see TasGrid::constructSurrogate().
 *
 * \param max_num_ranks defines the number of ranks to use, if it exceeds the size of the
 *      communicator then all ranks will be used in the process. If less than the size,
 *      the \b root process will not make calls to the model but it will be working
 *      on the internal grid methods. This is useful when dealing with many points where
 *      the cost of loading points is not negligible.
 *      If the \b max_num_ranks is less than the communicator size minus one, then some of
 *      the ranks will be completely idle (usually the ranks with higher number).
 *
 * \param tagx is the MPI_Send() tag to use when communicating the inputs \b x to the worker rank.
 * \param tagy is the MPI_Send() tag to use when communicating the outputs \b y to the root rank.
 *
 * \param root is the MPI rank within the communicator that has the initial copy of the \b grid
 *      as well as the refinement parameters that are identical to TasGrid::constructSurrogate().
 *      The final result from the construction will be held by the \b root rank only,
 *      the grids for the rest of other ranks will not be accessed.
 *
 * \param comm is the MPI communicator where sampling will take place, e.g., MPI_COMM_WORLD.
 *
 * \param grid on rank \b root is the grid to be constructed, ignored when the rank is different.
 *
 * \param tolerance refinement parameter, see TasGrid::constructSurrogate().
 * \param criteria refinement parameter, see TasGrid::constructSurrogate().
 * \param output refinement parameter, see TasGrid::constructSurrogate().
 * \param level_limits refinement parameter, see TasGrid::constructSurrogate().
 * \param scale_correction refinement parameter, see TasGrid::constructSurrogate().
 *
 * \param checkpoint_filename is either empty (skip checkpoints) or the filename to use
 *      to save the intermediate results. Note that only the \b root rank will issue
 *      read/write commands to the files.
 *
 * \par Example 1: Model that uses no MPI
 * Construct a surrogate to a model described by a parametrized system of equations
 * using a local polynomial gird and arbitrary number of MPI ranks. Rank 0 will be
 * responsible for the grid construction and computing the initial guess for all
 * solutions, the rest of the ranks will be calling a solver which makes no MPI calls.
 * The computational budget is 10,000 points and the solver works with one equation
 * at a time.
 * \code
 *  // solves the equation using the initial guess and makes no MPI calls
 *  void solver(std::vector<double> const &inputs,
 *              std::vector<double> &solution,
 *              std::vector<double> &initial guess);
 *  ...
 *  int me, num_ranks, tagx = 11, tagy = 12;
 *  MPI_Comm_rank(MPI_COMM_WORLD, &me); // rank of this process
 *  MPI_Comm_size(comm, &num_ranks); // all ranks
 *
 *  int num_inputs = ...;  // must match the solver
 *  int num_outputs = ...; // must match the solver
 *
 *  TasGrid::TasmanianSparseGrid grid;
 *  if (me == 0)
 *      grid = TasGrid::makeLocalPolynomialGrid(num_inputs, num_outputs, ...);
 *
 *  TasGrid::mpiConstructSurrogate<TasGrid::with_initial_guess>(
 *      [&](std::vector<double> const &x, std::vector<double> &y)
 *      ->void{
 *          auto initial_guess = y; // on input y is the initial guess
 *          std::vector<double> solution;
 *          solver(x, solution, initial_guess);
 *          y = solution; // on output y is the solution
 *      },
 *      num_dimensions, num_outputs,
 *      10000, 1, num_ranks - 1, tagx, tagy, 0, MPI_COMM_WORLD,
 *      grid, tolerance, TasGrid::refine_fds, ... );
 *
 *  // at this point, grid on rank 0 will hold the surrogate to solver()
 *  // grid on all other processes will be empty()
 * \endcode
 *
 * \par Example 2: Model that works on sub-communicator
 * Suppose that the model requires multiple MPI ranks for efficient computing
 * and it can operate on an arbitrary MPI communicator.
 * In this example, we take the MPI world communicator and we sub-divide it
 * into many smaller comms, then each sub-communicator will be used to handle
 * an independent sample.
 * In this example, all ranks participate in the process, assuming that
 * the cost of the grid manipulations is negligible compared to the model simulations.
 * There is no initial guess in the complicated model.
 * \code
 *  // model that computes the outputs using all ranks on an arbitrary communicator
 *  void model(std::vector<double> const &inputs,
 *             std::vector<double> &outputs,
 *             MPI_Comm comm);
 *  ...
 *  // number of ranks to assign to the communicators for each call to model()
 *  int ranks_per_model = ...;
 *  int world_me, num_ranks, tagx = 11, tagy = 12;
 *  MPI_Comm_rank(MPI_COMM_WORLD, &world_me); // rank of this process in world comm
 *  MPI_Comm comm_worker; // local comm for the work done by the model
 *  // split the world communicator into local ones
 *  MPI_Comm_split(MPI_COMM_WORLD, world_me / ranks_per_model, world_me % ranks_per_model, &comm_worker);
 *  int worker_me;
 *  MPI_Comm_rank(comm_worker, &worker_me); // this procees rank within the worker comm
 *  // create a sub-comm for the processes that would commonicate with Tasmanian
 *  MPI_Comm comm_tasmanian; // communicator for Tasmanian operations
 *  MPI_Comm_split(MPI_COMM_WORLD, (worker_me == 0) ? 0 : MPI_UNDEFINED, world_me, &comm_tasmanian);
 *  int num_tasmanian_ranks;
 *  if (worker_me == 0) MPI_Comm_size(comm_tasmanian, &num_tasmanian_ranks);
 *  int tasmanian_root = 0; // world_me = zero is the same as tasmanian_root = zero
 *
 *  int num_inputs = ...;  // must match the solver
 *  int num_outputs = ...; // must match the solver
 *
 *  if (worker_me == 0){
 *      mpiConstructSurrogate(
 *          [&](...)->void{
 *              MPI_Bcast(... sent x to all ranks in comm_worker ...);
 *              model(..., comm_worker); // participate in this model evaluation
 *              MPI_Gather(... collect all results in this process ...);
 *          },
 *      num_inputs, num_outputs, 10000, 1, num_tasmanian_ranks, tagx, tagy,
 *      tasmanian_root, comm_tasmanian,
 *      grid, ...);
 *  }else{
 *      MPI_Bcast(... get x from rank 0 on comm_worker ...);
 *      model(x, outputs, comm_worker);
 *      MPI_Gather(... gather the results in rank 0 on comm_worker ...);
 *  }
 * \endcode
 *  \b Note: the code above demonstrates the idea of the split world communicator,
 *  an actual Implementation would differ in details and specifics.
 */
template<bool use_initial_guess = no_initial_guess>
void mpiConstructSurrogate(ModelSignatureMPI model,
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

/*!
 * \ingroup TasmanianAddonsMPIConstruct
 * \brief Construct a sparse grid surrogate to the model defined by the lambda, MPI version.
 *
 * Uses the user provided \b anisotropic_weights to order the samples by importance,
 * see\n TasmanianSparseGrid::getCandidateConstructionPoints().
 * and the overloads to TasGrid::constructSurrogate().
 */
template<bool use_initial_guess = no_initial_guess>
void mpiConstructSurrogate(ModelSignatureMPI model,
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

/*!
 * \ingroup TasmanianAddonsMPIConstruct
 * \brief Construct a sparse grid surrogate to the model defined by the lambda, MPI version.
 *
 * Uses anisotropic weights to order the samples by importance,
 * starts with a fully isotropic grid until enough points are loaded to allow to estimate the weights,
 * see TasmanianSparseGrid::getCandidateConstructionPoints().
 * and the overloads to TasGrid::constructSurrogate().
 */
template<bool use_initial_guess = no_initial_guess>
void mpiConstructSurrogate(ModelSignatureMPI model,
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

/*!
 * \ingroup TasmanianAddonsMPIConstruct
 * \brief Executes the worker (i.e., non-root) portion of the MPI sampling.
 *
 * Using any of the TasGrid::mpiConstructSurrogate() overloads, the grid, the budget and the refinement
 * parameters are accessed only by the root rank, the rest is not used by the workers.
 * This template can be instantiated for the non-root rank with only the relevant inputs.
 *
 * - The call to this method \b must be paired with a call to TasGrid::mpiConstructSurrogate() on the root rank.
 * - The inputs (including the template parameter) \b must match the ones in the call to TasGrid::mpiConstructSurrogate().
 *
 */
template<bool use_initial_guess = no_initial_guess>
void mpiConstructWorker(ModelSignatureMPI model,
                        int num_dimensions, int num_outputs, size_t max_samples_per_job,
                        size_t max_num_ranks, int tagx, int tagy, int root, MPI_Comm comm){

    TasmanianSparseGrid grid; // workers use an empty grid
    mpiConstructCommon<use_initial_guess>(model, num_dimensions, num_outputs, 0, max_samples_per_job,
                                          max_num_ranks, tagx, tagy, root, comm, grid,
                                          [&](TasmanianSparseGrid&)->std::vector<double>{
                                               return std::vector<double>();
                                          }, std::string());
}

}

#endif // Tasmanian_ENABLE_MPI

#endif
