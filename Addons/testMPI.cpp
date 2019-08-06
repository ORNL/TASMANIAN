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

#include "testMPI.hpp"
#include "testMPIDream.hpp"

using std::cout;
using std::setw;

int main(int argc, char ** argv){

    //MPI_Init(&argc, &argv); // MPI_THREAD_MULTIPLE requires MPI_Init_thread()
    int threads_available;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &threads_available);

    int me = TasGrid::getMPIRank(MPI_COMM_WORLD);
    if (me == 0) cout << "\n";

    // --------------- Send/Recv <ascii> ----------------- //
    if (!testSendReceive<TasGrid::mode_ascii>()) throw std::runtime_error("failed Send/Recv ascii");
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0)
        cout << "    MPI Send/Recv     <ascii>    Pass\n";

    // --------------- Send/Recv <binary> ----------------- //
    if (!testSendReceive<TasGrid::mode_binary>()) throw std::runtime_error("failed Send/Recv binary");
    MPI_Barrier(MPI_COMM_WORLD);
    if (!testLikelySendRecv()) throw std::runtime_error("failed Send/Recv DREAM");
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0)
        cout << "    MPI Send/Recv    <binary>    Pass\n";

    // ----------------- Bcast <ascii> ------------------ //
    if (!testBcast<TasGrid::mode_ascii>()) throw std::runtime_error("failed Bcast ascii");
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0)
        cout << "        MPI Bcast     <ascii>    Pass\n";

    // ----------------- Bcast <binary> ----------------- //
    if (!testBcast<TasGrid::mode_binary>()) throw std::runtime_error("failed Bcast binary");
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0)
        cout << "        MPI Bcast    <binary>    Pass\n";

    // ----------------- Scatter <ascii> --------------- //
    if (!testScatterOutputs<TasGrid::mode_ascii>()) throw std::runtime_error("failed Scatter ascii");
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0)
        cout << "      MPI Scatter     <ascii>    Pass\n";

    // ----------------- Scatter <binary> --------------- //
    if (!testScatterOutputs<TasGrid::mode_binary>()) throw std::runtime_error("failed Scatter binary");
    MPI_Barrier(MPI_COMM_WORLD);
    if (!testLikelyScatter()) throw std::runtime_error("failed Scatter DREAM");
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0)
        cout << "      MPI Scatter    <binary>    Pass\n";

    if (threads_available == MPI_THREAD_MULTIPLE){
        // ----------------- Construct <no guess> ----------- //
        testMPIconstruct<no_initial_guess>();
        MPI_Barrier(MPI_COMM_WORLD);
        if (me == 0)
            cout << "    MPI Construct   <no-init>    Pass\n";

        // ----------------- Construct <with guess> --------- //
        testMPIconstruct<with_initial_guess>();
        MPI_Barrier(MPI_COMM_WORLD);
        if (me == 0)
            cout << "    MPI Construct <with-init>    Pass\n";

        // ----------------- Construct <with guess> --------- //
        testMPIconstructStrict();
        MPI_Barrier(MPI_COMM_WORLD);
        if (me == 0)
            cout << "    MPI Construct    <strict>    Pass\n";

    }else{
        if (me == 0)
            cout << "\n Skipping MPI construction since this version of MPI\n does not seem to support MPI_THREAD_MULTIPLE\n\n";
    }

    testMPIDream();
    MPI_Barrier(MPI_COMM_WORLD);
    if (me == 0)
        cout << "        MPI Dream  <sampling>    Pass\n";

    // --------------- Finalize ------------------------- //
    MPI_Finalize();

    if (me == 0) cout << endl;
    return 0;
}
