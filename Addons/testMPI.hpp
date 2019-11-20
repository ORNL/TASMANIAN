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

#include "TasmanianAddons.hpp"
#include "tasgridCLICommon.hpp"

/*!
 * \brief Check if the points in the two grids match bit-wise.
 */
inline bool checkPoints(TasGrid::TasmanianSparseGrid const &gridA, TasGrid::TasmanianSparseGrid const &gridB){
    if (gridA.getNumPoints()     != gridB.getNumPoints())     return false;
    if (gridA.getNumDimensions() != gridB.getNumDimensions()) return false;
    auto pA = gridA.getPoints();
    auto pB = gridB.getPoints();
    double err = 0.0;
    for(auto x = pA.begin(), y = pB.begin(); x != pA.end(); x++, y++) err += std::abs(*x - *y);
    return (err == 0.0); // bit-wise match is reasonable to expect here, but use with caution for some grids
}

/*!
 * \brief Simple test of MPI Send/Recv of sparse grids, binary and ascii formats.
 */
template<bool use_binary>
bool testSendReceive(){
    MPI_Comm comm = MPI_COMM_WORLD;
    auto true_grid = TasGrid::makeGlobalGrid(5, 3, 4, TasGrid::type_level, TasGrid::rule_clenshawcurtis);

    int tag = 1;
    int const me = TasGrid::getMPIRank(comm);

    if (me == 0){
        return (TasGrid::MPIGridSend<use_binary>(true_grid, 1, tag, comm) == MPI_SUCCESS);
    }else if (me == 1){
        MPI_Status status;
        TasGrid::TasmanianSparseGrid grid;
        auto result = TasGrid::MPIGridRecv<use_binary>(grid, 0, tag, comm, &status);
        if (result != MPI_SUCCESS) return false;
        return checkPoints(true_grid, grid);
    }else{
        return true;
    }
}

/*!
 * \brief Simple test of MPI Send/Recv of sparse grids, binary and ascii formats.
 */
template<bool use_binary>
bool testBcast(){
    MPI_Comm comm = MPI_COMM_WORLD;
    auto true_grid = TasGrid::makeGlobalGrid(5, 1, 4, TasGrid::type_level, TasGrid::rule_clenshawcurtis);

    int const me = TasGrid::getMPIRank(comm);

    if (me == 1){ // using proc 1 to Bcast the grid
        return (TasGrid::MPIGridBcast<use_binary>(true_grid, 1, comm) == MPI_SUCCESS);
    }else{
        TasGrid::TasmanianSparseGrid grid;
        auto result = TasGrid::MPIGridBcast<use_binary>(grid, 1, comm);
        if (result != MPI_SUCCESS) return false;
        return checkPoints(true_grid, grid);
    }
}

/*!
 * \brief Simple test of MPI Scatter Outputs of sparse grids, binary and ascii formats.
 */
template<bool use_binary>
bool testScatterOutputs(){
    // grid has 7 outputs split between 3 ranks gives (3 2 2)
    MPI_Comm comm = MPI_COMM_WORLD;
    int const me = TasGrid::getMPIRank(comm);

    auto reference_grid = TasGrid::makeGlobalGrid(3, (me == 0) ? 3 : 2, 4, TasGrid::type_level, TasGrid::rule_clenshawcurtis);

    loadNeededPoints<false, false>([&](double const x[], double y[], size_t)->void{
                                        double expval = std::exp(x[0] + x[1] + x[2]); // 3 inputs
                                        if (me == 0){
                                            y[0] = expval;
                                            y[1] = 2.0 * expval;
                                            y[2] = 3.0 * expval;
                                        }else if (me == 1){
                                            y[0] = 4.0 * expval;
                                            y[1] = 5.0 * expval;
                                        }else{
                                            y[0] = 6.0 * expval;
                                            y[1] = 7.0 * expval;
                                        }
                                    }, reference_grid, 0);

    TasmanianSparseGrid grid; // received grid

    // use rank 1 for the root
    if (me == 1){
        auto full_grid = TasGrid::makeGlobalGrid(3, 7, 4, TasGrid::type_level, TasGrid::rule_clenshawcurtis);
        loadNeededPoints<false, false>([&](double const x[], double y[], size_t)->void{
                                            double expval = std::exp(x[0] + x[1] + x[2]); // 3 inputs
                                            for(size_t i=0; i<7; i++)
                                                y[i] = double(i+1) * expval;
                                        }, full_grid, 0);
        MPIGridScatterOutputs<use_binary>(full_grid, grid, 1, 2, comm);
    }else{
        MPIGridScatterOutputs<use_binary>(TasmanianSparseGrid(), grid, 1, 2, comm);
    }

    std::minstd_rand park_miller(99);
    std::uniform_real_distribution<double> unif(-1.0, 1.0);
    std::vector<double> test_points(3 * 1000);
    for(auto &t : test_points) t = unif(park_miller);

    auto match = [&](TasmanianSparseGrid const &a, TasmanianSparseGrid const &b)->bool{
        std::vector<double> resa, resb; // reference and actual result
        a.evaluateBatch(test_points, resa);
        b.evaluateBatch(test_points, resb);
        double err = 0.0;
        for(auto ia = resa.begin(), ib = resb.begin(); ia != resa.end(); ia++, ib++)
            err = std::max(err, std::abs(*ia - *ib));
        return (err < 1.E-13);
    };

    if (!match(grid, reference_grid)) throw std::runtime_error("ERROR: first iteration of MPIGridScatterOutputs() failed.");

    MPIGridScatterOutputs<use_binary>(copyGrid(grid), grid, 1, 2, comm);
    if (me == 2){
        if (!grid.empty()) throw std::runtime_error("ERROR: second iteration of MPIGridScatterOutputs() failed.");
    }else{
        reference_grid = TasGrid::makeGlobalGrid(3, 1, 4, TasGrid::type_level, TasGrid::rule_clenshawcurtis);
        loadNeededPoints<false, false>([&](double const x[], double y[], size_t)->void{
                                            double expval = std::exp(x[0] + x[1] + x[2]); // 3 inputs
                                            y[0] = ((me == 0) ? 4.0 : 5.0) * expval;
                                        }, reference_grid, 0);
        if (!match(grid, reference_grid)) throw std::runtime_error("ERROR: second iteration of MPIGridScatterOutputs() failed.");
    }

    MPIGridScatterOutputs<use_binary>(copyGrid(grid), grid, 1, 2, comm);
    if (me == 0){
        reference_grid = TasGrid::makeGlobalGrid(3, 1, 4, TasGrid::type_level, TasGrid::rule_clenshawcurtis);
        loadNeededPoints<false, false>([&](double const x[], double y[], size_t)->void{
                                            y[0] = 5.0 * std::exp(x[0] + x[1] + x[2]);
                                        }, reference_grid, 0);
        if (!match(grid, reference_grid)) throw std::runtime_error("ERROR: third iteration of MPIGridScatterOutputs() failed.");
    }else{
        if (!grid.empty()) throw std::runtime_error("ERROR: third iteration of MPIGridScatterOutputs() failed.");
    }

    return true;
}

template<bool use_initial_guess>
void testMPIconstruct(){
    MPI_Comm comm = MPI_COMM_WORLD;
    int const me = TasGrid::getMPIRank(comm);

    std::minstd_rand park_miller(99);
    std::uniform_real_distribution<double> unif(-1.0, 1.0);
    std::vector<double> test_points(3 * 1000);
    for(auto &t : test_points) t = unif(park_miller);

    auto match = [&](TasmanianSparseGrid const &a, TasmanianSparseGrid const &b)->bool{
        std::vector<double> resa, resb; // reference and actual result
        a.evaluateBatch(test_points, resa);
        b.evaluateBatch(test_points, resb);
        double err = 0.0;
        for(auto ia = resa.begin(), ib = resb.begin(); ia != resa.end(); ia++, ib++)
            err = std::max(err, std::abs(*ia - *ib));
        constexpr double tolerance = 1.E-2;
        if (err >= tolerance) std::cout << "error = " << err << "  expected " << tolerance << std::endl;
        return (err < tolerance);
    };

    auto model = [&](std::vector<double> const &x, std::vector<double> &y)->void{
        size_t num_samples = x.size() / 3;
        if (use_initial_guess == with_initial_guess)
            y.resize(num_samples * 2); // y can be empty
        for(size_t i=0; i<num_samples; i++){ // for each sample
            y[2*i + 0] = std::exp(x[3*i + 0] + x[3*i + 1] + x[3*i + 2]);
            y[2*i + 1] = std::sin(x[3*i + 0]) * std::cos(x[3*i + 1]) + std::sin(x[3*i + 2]) * std::cos(x[3*i + 1]);
        }
    };
    auto modelt = [&](std::vector<double> const &x, std::vector<double> &y, size_t)->void{
        model(x, y);
    };

    auto grid = TasGrid::makeLocalPolynomialGrid(3, 2, 3);

    if (me == 0){
        mpiConstructSurrogate<use_initial_guess>(model, 3, 2, 1000, 2, 3, 11, 22, 0, comm,
                                                 grid, 1.E-5, refine_classic, -1);
    }else{
        mpiConstructWorker<use_initial_guess>(model, 3, 2, 2, 3, 11, 22, 0, comm);
    }

    if (me == 0){
        auto reference_grid = TasGrid::makeLocalPolynomialGrid(3, 2, 3);
        constructSurrogate<mode_sequential>(modelt, 1000, 0, 1, reference_grid, 1.E-5, refine_classic, -1);
        if (!match(grid, reference_grid)) throw std::runtime_error("testMPIconstruct() grids mismatch.");
    }
}

// this test must produce grids that match to within numeric precision
// no matter the order of samples or any other considerations
void testMPIconstructStrict(){
    MPI_Comm comm = MPI_COMM_WORLD;
    int const me = TasGrid::getMPIRank(comm);

    auto match = [](TasmanianSparseGrid const &a, TasmanianSparseGrid const &b)->bool{
        if (a.getNumLoaded() != b.getNumLoaded()) return false;
        auto pa = a.getLoadedPoints();
        auto pb = b.getLoadedPoints();
        double err = 0.0;
        for(auto ia = pa.begin(), ib = pb.begin(); ia != pa.end(); ia++, ib++)
            err = std::max(err, std::abs(*ia - *ib));
        if (err >= Maths::num_tol)
            cout << "points mismatch: " << err << endl;
        return (err < Maths::num_tol);
    };

    auto modelt = [&](std::vector<double> const &x, std::vector<double> &y, size_t)->void{
        size_t num_samples = x.size() / 2;
        for(size_t i=0; i<num_samples; i++) // for each sample
            y[i] = std::exp(x[2*i] + x[2*i + 1]);
    };
    auto model = [&](std::vector<double> const &x, std::vector<double> &y)->void{
        modelt(x, y, 0);
        if (me == 0) throw std::runtime_error("ERROR: rank 0 should not participate in this.");
        //std::this_thread::sleep_for(std::chrono::milliseconds(1100));
    };

    std::vector<int> aweights = {1, 2};
    auto reference_grid = TasGrid::makeSequenceGrid(2, 1, 6, TasGrid::type_level, TasGrid::rule_leja, aweights);
    auto grid = copyGrid(reference_grid);
    loadNeededPoints<mode_sequential>(modelt, reference_grid, 0);

    mpiConstructSurrogate<no_initial_guess>
        (model, 2, 1, reference_grid.getNumLoaded(), 1, 2, 11, 22, 0, comm, grid, TasGrid::type_iptotal, aweights);

    if (me == 0){
        if (!match(grid, reference_grid)) throw std::runtime_error("ERROR: CV construction failed.");
    }

}
