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

#ifndef __TASMANIAN_ADDONS_MPIGRIDSCATTER_HPP
#define __TASMANIAN_ADDONS_MPIGRIDSCATTER_HPP

/*!
 * \internal
 * \file tsgMPIScatterGrid.hpp
 * \brief Sparse Grids send/receive through MPI.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianAddonsCommon
 *
 * Templates that communicate sparse grids through MPI commands.
 * \endinternal
 */

#include "tsgAddonsCommon.hpp"

/*!
 * \ingroup TasmanianAddons
 * \addtogroup TasmanianAddonsMPIGridSend MPI Send/Receive/Bcast/Scatter for Sparse Grids
 *
 * Methods to send/receive TasGrid::TasmanianSparseGrid objects.
 * The syntax mimics the raw MPI_Send and MPI_Recv calls,
 * and the templates require Tasmanian_ENABLE_MPI=ON.
 */

#ifdef Tasmanian_ENABLE_MPI

namespace TasGrid{

/*!
 * \internal
 * \ingroup TasmanianAddonsMPIGridSend
 * \brief Coverts a vector to basic stream-buffer.
 *
 * The stream-buffer is the data portion of the stream, this class converts
 * a char-vector to such buffer.
 * The stream buffer will simply assume the begin, end and data-type pointers
 * from the vector, therefore, using iterator invalidating operations on
 * the vector that invalidate will also invalidate the stream-buffer.
 * \endinternal
 */
class VectorToStreamBuffer : public std::basic_streambuf<char, std::char_traits<char>>{
public:
    //! \brief Make a stream-buffer from the \b data vector.
    VectorToStreamBuffer(std::vector<char> &data){
        setg(&*data.begin(), data.data(), &*data.end());
    }
};

/*!
 * \internal
 * \ingroup TasmanianAddonsMPIGridSend
 * \brief Utility to return the rank within the given comm.
 *
 * \endinternal
 */
inline int getMPIRank(MPI_Comm comm){ int rank; MPI_Comm_rank(comm, &rank); return rank; }

/*!
 * \ingroup TasmanianAddonsMPIGridSend
 * \brief Send a grid to another process in the MPI comm.
 *
 * Send the grid to the destination rank across an MPI comm.
 * The grid is frist writtein to a stream in either binary or ASCII
 * format and then send with a single MPI_Send() call.
 * Binary results in smaller messages and less computational overhead;
 * thus, ASCII is provided mostly for debugging purposes.
 *
 * \tparam binary defines whether to use binary (\b true) or ASCII (\b false) mode.
 *                Recommended use with constexpr constant TasGrid::mode_binary
 *                and TasGrid::mode_ascii.
 *
 * \param grid        is the grid to send.
 * \param destination is the rank of the recipient MPI process.
 * \param tag         is the tag to use for the MPI message.
 * \param comm        is the MPI comm where the source and destination reside.
 *
 * \return the error code of the MPI_Send() command.
 *
 * \b Note: this call must be mirrored by TasGrid::MPIGridRecv() on the
 *          destination process.
 *
 * Example usage, process 0 creates a grid and sends it to process 1:
 * \code
 *   int me, tag = 0;
 *   MPI_Comm_rank(MPI_COMM_WORLD, &me);
 *   MPI_Status status;
 *   TasGrid::TasmanianSparseGrid grid;
 *   if (me == 0) grid.makeGlobalGrid(3, 4, 5, TasGrid::type_level, TasGrid::rule_clenshawcurtis);
 *   if (me == 0) MPIGridSend(grid, 1, tag, MPI_COMM_WORLD);
 *   else if (me == 1) MPIGridRecv(grid, 0, tag, MPI_COMM_WORLD, &status);
 *   // at this line, process 1 has a grid equivalent to that of process 0
 *   // processes with rank 2 and above do nothing, i.e., they have an empty grid
 * \endcode
 */
template<bool binary = TasGrid::mode_binary>
int MPIGridSend(TasmanianSparseGrid const &grid, int destination, int tag, MPI_Comm comm){
    std::stringstream ss;
    grid.write(ss, binary);
    while(ss.str().size() % 16 != 0) ss << " ";
    return MPI_Send(ss.str().c_str(), (int) (ss.str().size() / 16), MPI_LONG_DOUBLE, destination, tag, comm);
}

/*!
 * \ingroup TasmanianAddonsMPIGridSend
 * \brief Receive a grid from another process in the MPI comm.
 *
 * Receive a grid that has been send with TasGrid::MPIGridSend().
 * This call intercepts both messages and compiles them into a sparse grid object.
 *
 * \tparam binary defines whether to use binary (\b true) or ASCII (\b false) mode.
 *                Recommended use with constexpr constant TasGrid::mode_binary
 *                and TasGrid::mode_ascii.
 *
 * \param grid is the output grid, it will be overwritten with grid send by
 *             the source rank similar to TasGrid::TasmanianSparseGrid::read().
 * \param source   is the rank of the process in the MPI comm that issued the send command.
 * \param tag      is the tag used in the MPI send command.
 * \param comm     is the MPI comm where the source and destination reside.
 * \param status is the status of the MPI_Recv() command.
 *
 * \return the error code of the MPI_Recv() command.
 */
template<bool binary = TasGrid::mode_binary>
int MPIGridRecv(TasmanianSparseGrid &grid, int source, int tag, MPI_Comm comm, MPI_Status *status = MPI_STATUS_IGNORE){
    MPI_Status internal_status;
    if (status == MPI_STATUS_IGNORE) status = &internal_status;

    int short_data_size;
    MPI_Probe(source, tag, comm, status);
    MPI_Get_count(status, MPI_LONG_DOUBLE, &short_data_size);

    size_t data_size = Utils::size_mult(short_data_size, 16);

    std::vector<char> buff(data_size);
    auto result = MPI_Recv(buff.data(), (int) (data_size / 16), MPI_LONG_DOUBLE, source, tag, comm, status);

    VectorToStreamBuffer data_buffer(buff); // do not modify buff after this point
    std::istream is(&data_buffer);
    grid.read(is, binary);
    return result;
}

/*!
 * \ingroup TasmanianAddonsMPIGridSend
 * \brief Broadcast a grid to all processes in an MPI comm.
 *
 * Make all \b grid variables for all process in the \b comm match the grid
 * on the \b root process.
 * This call uses two MPI_Bcast() calls, the grid size (in memory units)
 * and the actual grid data.
 * The transfer can be done in either binary or ASCII format, but binary results
 * in smaller messages and less computational overhead;
 * thus, ASCII is provided mostly for debugging purposes.
 *
 * \tparam binary defines whether to use binary (\b true) or ASCII (\b false) mode.
 *                Recommended use with constexpr constant TasGrid::mode_binary
 *                and TasGrid::mode_ascii.
 *
 * \param grid is the grid to broadcast across the MPI comm, the grid on the \b root
 *             process will not be modified (i.e., treat as const),
 *             in all other cases, the grid will be overwritten similar to
 *             the TasGrid::TasmanianSparseGrid::read().
 * \param root is the process that holds the data that needs to be send across.
 * \param comm is the MPI comm of all process that need to share the grid.
 *
 * \return the error code of the fist failed MPI_Bcast() command
 *         (corresponding to either the size of the data message),
 *         if MPI_SUCCESS is returned then both messages were successful.
 *
 * Example usage, process 0 reads a grid from a file and sends it to all processes:
 * \code
 *   int me;
 *   MPI_Comm_rank(MPI_COMM_WORLD, &me);
 *   auto grid = (me == 0) ? TasGrid::readGrid("foo") : TasGrid::TasmanianSparseGrid();
 *   MPIGridBcast(grid, 0, MPI_COMM_WORLD);
 *   // at this line, every process has the same grid as if they all read it from "foo"
 * \endcode
 */
template<bool binary = TasGrid::mode_binary>
int MPIGridBcast(TasmanianSparseGrid &grid, int root, MPI_Comm comm){
    int me = getMPIRank(comm); // my rank within the comm
    if (me == root){ // sends the grid
        std::stringstream ss;
        grid.write(ss, binary);

        while(ss.str().size() % 16 != 0) ss << " "; // pad with empty chars to align to 16 bytes, i.e., long double

        unsigned long long data_size = (unsigned long long) ss.str().size();
        auto result = MPI_Bcast(&data_size, 1, MPI_UNSIGNED_LONG_LONG, me, comm);
        if (result != MPI_SUCCESS) return result;

        return MPI_Bcast(const_cast<char*>(ss.str().c_str()), (int) (data_size / 16), MPI_LONG_DOUBLE, me, comm); // Bcast root does not modify the buffer, this is const-correct
    }else{ // receives the grid
        unsigned long long data_size;

        auto result = MPI_Bcast(&data_size, 1, MPI_UNSIGNED_LONG_LONG, root, comm);
        if (result != MPI_SUCCESS) return result;

        std::vector<char> buff((size_t) data_size);
        result = MPI_Bcast(buff.data(), (int) (buff.size() / 16), MPI_LONG_DOUBLE, root, comm);

        VectorToStreamBuffer data_buffer(buff); // do not modify buff after this point
        std::istream is(&data_buffer);
        grid.read(is, binary);
        return result;
    }
}

/*!
 * \ingroup TasmanianAddonsMPIGridSend
 * \brief Split the grid across the comm where each rank receives an equal portion of the total outputs.
 *
 * Split the grid across the ranks in the comm so that each rank receives a grid with the same type, rule,
 * points, etc., but with a subset of the total outputs. The distribution is as close to even as possible,
 * if there are less outputs than ranks, some ranks will receive an empty grid.
 *
 * \b Note: this does not use MPI_Scatter(), instead it makes multiple calls to MPIGridSend() and MPIGridRecv().
 *
 * \tparam binary defines whether to use binary (\b true) or ASCII (\b false) mode, see TasGrid::MPIGridSend().
 *
 * \param source grid located on the \b root rank is the grid to be distributed across,
 *               for all other ranks the source will not be referenced.
 * \param destination is the grid where the local portion of the scatter will be stored,
 *                    the existing grid will be overwritten.
 *                    If the \b source outputs are less than the number of \b comm ranks,
 *                    then some of the destination grids will be empty.
 * \param root is the rank that will hold the source sparse grid.
 * \param tag  same as in TasGrid::MPIGridSend().
 * \param comm is the MPI comm of all process that need to share a portion of the grid.
 *
 * Example usage, rank 0 creates a large grid and scatters is across comm:
 * \code
 *   int me;
 *   MPI_Comm_rank(MPI_COMM_WORLD, &me);
 *   auto source = (me == 0) ? TasGrid::readGrid("grid_with_many_outputs") : TasGrid::TasmanianSparseGrid();
 *   TasGrid::TasmanianSparseGrid grid;
 *   MPIGridScatterOutputs(source, grid, 0, 1, 2, comm);
 *   // at this line, every process has a portion of the source at grid
 *   // if the comm has 3 ranks,
 *   // then 7 outputs will be split into 3 2 2
 *   // and 2 outputs will become 1 1 empty
 * \endcode
 */
template<bool binary = TasGrid::mode_binary>
int MPIGridScatterOutputs(TasmanianSparseGrid const &source, TasmanianSparseGrid &destination, int root, int tag, MPI_Comm comm){
    int me = getMPIRank(comm); // my rank within the comm

    if (me == root){ // splitting and sending the grid
        int num_ranks; MPI_Comm_size(comm, &num_ranks);
        int num_effective_ranks = std::min(num_ranks, source.getNumOutputs());

        int stride = source.getNumOutputs() / num_effective_ranks;
        int extras = source.getNumOutputs() % num_effective_ranks;

        // return the starting offset for the given rank
        auto offset = [&](int rank)->int{ return rank * stride + std::min(rank, extras); };

        for(int rank=0; rank<num_effective_ranks; rank++){
            if (rank == root){ // this is me, take my own copy of the grid
                destination = copyGrid(source, offset(rank), offset(rank+1));
            }else{ // send the grid out
                auto result = MPIGridSend<binary>(copyGrid(source, offset(rank), offset(rank+1)) , rank, tag, comm);
                if (result != MPI_SUCCESS) return result;
            }
        }
        for(int rank=num_effective_ranks; rank<num_ranks; rank++){ // if there are any grid remaining, set those to empty
            if (rank == root){ // this is me, take my own copy of the grid
                destination = makeEmpty();
            }else{ // send the grid out
                auto result = MPIGridSend<binary>(makeEmpty() , rank, tag, comm);
                if (result != MPI_SUCCESS) return result;
            }
        }
        return MPI_SUCCESS; // if we got here, all was successful
    }else{ // receiving a grid
        return MPIGridRecv<binary>(destination, root, tag, comm);
    }
}

}

#endif // Tasmanian_ENABLE_MPI

#endif
