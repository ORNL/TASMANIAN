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

inline bool testLikelySendRecv(){
    int me = TasGrid::getMPIRank(MPI_COMM_WORLD);
    TasDREAM::LikelihoodGaussIsotropic ref_isolike(10.0, {1.0, 2.0, 3.0});

    if (me == 0){
        if (TasDREAM::MPILikelihoodSend(ref_isolike, 1, 11, MPI_COMM_WORLD) != MPI_SUCCESS) return false;
    }else if (me == 1){
        TasDREAM::LikelihoodGaussIsotropic isolike;
        if (TasDREAM::MPILikelihoodRecv(isolike, 0, 11, MPI_COMM_WORLD) != MPI_SUCCESS) return false;
        std::vector<double> model = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}; // full rank matrix to cover all entries
        std::vector<double> result(3), true_result(3);
        isolike.getLikelihood(TasDREAM::logform, model, result);
        ref_isolike.getLikelihood(TasDREAM::logform, model, true_result);
        for(size_t i=0; i<3; i++) if (std::abs(result[i] - true_result[i]) > TasGrid::Maths::num_tol) return false;
    }

    TasDREAM::LikelihoodGaussAnisotropic ref_alike({4.0, 5.0, 6.0}, {1.0, 2.0, 3.0});

    if (me == 1){
        if (TasDREAM::MPILikelihoodSend(ref_alike, 2, 12, MPI_COMM_WORLD) != MPI_SUCCESS) return false;
    }else if (me == 2){
        TasDREAM::LikelihoodGaussAnisotropic alike;
        if (TasDREAM::MPILikelihoodRecv(alike, 1, 12, MPI_COMM_WORLD) != MPI_SUCCESS) return false;
        std::vector<double> model = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}; // full rank matrix to cover all entries
        std::vector<double> result(3), true_result(3);
        alike.getLikelihood(TasDREAM::logform, model, result);
        ref_alike.getLikelihood(TasDREAM::logform, model, true_result);
        for(size_t i=0; i<3; i++) if (std::abs(result[i] - true_result[i]) > TasGrid::Maths::num_tol) return false;
    }

    return true;
}

inline bool testLikelyScatter(){
    int me = TasGrid::getMPIRank(MPI_COMM_WORLD);
    int tag = 11;
    TasDREAM::LikelihoodGaussIsotropic source(10.0, {1.0, 2.0, 3.0});

    TasDREAM::LikelihoodGaussIsotropic reference;
    TasDREAM::LikelihoodGaussIsotropic destination;
    if (me == 0){
        reference.setData(10.0, {1.0});
        MPILikelihoodScatter(source, destination, 0, tag, MPI_COMM_WORLD);
    }else if (me == 1){
        reference.setData(10.0, {2.0});
        MPILikelihoodScatter(TasDREAM::LikelihoodGaussIsotropic(), destination, 0, tag, MPI_COMM_WORLD);
    }else{
        reference.setData(10.0, {3.0});
        MPILikelihoodScatter(TasDREAM::LikelihoodGaussIsotropic(), destination, 0, tag, MPI_COMM_WORLD);
    }

    std::vector<double> model = {1.0}; // full rank matrix to cover all entries
    std::vector<double> result(1), true_result(1);
    destination.getLikelihood(TasDREAM::logform, model, result);
    reference.getLikelihood(TasDREAM::logform, model, true_result);
    if (std::abs(result[0] - true_result[0]) > TasGrid::Maths::num_tol) return false;

    TasDREAM::LikelihoodGaussAnisotropic asource({14.0, 15.0}, {1.0, 2.0});

    TasDREAM::LikelihoodGaussAnisotropic areference;
    TasDREAM::LikelihoodGaussAnisotropic adestination;
    if (me == 0){
        areference.setData({14.0}, {1.0});
        MPILikelihoodScatter(asource, adestination, 0, tag, MPI_COMM_WORLD);
    }else if (me == 1){
        areference.setData({15.0}, {2.0});
        MPILikelihoodScatter(TasDREAM::LikelihoodGaussAnisotropic(), adestination, 0, tag, MPI_COMM_WORLD);
    }else{
        MPILikelihoodScatter(TasDREAM::LikelihoodGaussAnisotropic(), adestination, 0, tag, MPI_COMM_WORLD);
    }

    if (me != 2){
        result = {0.0};
        true_result = {11.0};
        adestination.getLikelihood(TasDREAM::logform, model, result);
        areference.getLikelihood(TasDREAM::logform, model, true_result);
        if (std::abs(result[0] - true_result[0]) > TasGrid::Maths::num_tol) return false;
    }else{
        if (adestination.getNumOutputs() != 0){
            std::cout << "last rank did not receive empty likelihood." << std::endl;
            return false;
        }
    }

    return true;
}

void testMPIDream(){
    int num_chains = 10;

    int me = TasGrid::getMPIRank(MPI_COMM_WORLD);
    auto full_grid = TasGrid::makeSequenceGrid(2, 7, 2, TasGrid::type_level, TasGrid::rule_rleja);
    TasGrid::loadNeededPoints<mode_sequential>([&](std::vector<double> const &x, std::vector<double> &y, size_t)->void{
            y = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
            for(size_t i=0; i<y.size(); i++)
                y[i] += 0.1 * x[i%2]; // add perturbation to y
        }, full_grid, 0);
    auto grid = TasGrid::makeEmpty();

    TasDREAM::LikelihoodGaussAnisotropic full_likelihood({0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7}, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0});
    TasDREAM::LikelihoodGaussAnisotropic likely;

    TasDREAM::TasmanianDREAM full_state(num_chains, 2);
    TasDREAM::TasmanianDREAM state = (me == 0) ? full_state : TasDREAM::TasmanianDREAM();

    TasGrid::MPIGridScatterOutputs(full_grid, grid, 0, 11, MPI_COMM_WORLD);
    TasDREAM::MPILikelihoodScatter(full_likelihood, likely, 0, 13, MPI_COMM_WORLD);

    std::minstd_rand park_miller_init(42), park_miller1(77), park_miller2(77);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    std::vector<double> initial_state;
    TasDREAM::genGaussianSamples({0.0, 0.0}, {0.2, 0.2}, num_chains, initial_state, [&]()->double{ return unif(park_miller_init); });
    if (me == 0) state.setState(initial_state);
    full_state.setState(initial_state);

    TasDREAM::SampleDREAM(10, 10,
        TasDREAM::DistributedPosterior<TasDREAM::regform>(grid, likely, TasDREAM::uniform_prior, 2, num_chains, 0, MPI_COMM_WORLD),
        grid.getDomainInside(),
        state,
        TasDREAM::dist_uniform, 0.05,
        TasDREAM::const_percent<50>,
        [&]()->double{ return unif(park_miller1); }
    );

    TasDREAM::SampleDREAM(10, 10,
        TasDREAM::posterior<TasDREAM::regform>(full_grid, full_likelihood, TasDREAM::uniform_prior),
        grid.getDomainInside(),
        full_state,
        TasDREAM::dist_uniform, 0.05,
        TasDREAM::const_percent<50>,
        [&]()->double{ return unif(park_miller2); }
    );

    if (me == 0){
        std::vector<double> mean, variance;
        state.getHistoryMeanVariance(mean, variance);
        std::vector<double> ref_mean, ref_variance;
        full_state.getHistoryMeanVariance(ref_mean, ref_variance);
        if (((std::abs(mean[0] - ref_mean[0]) + std::abs(mean[1] - ref_mean[1])) > 1.E-9) ||
            ((std::abs(variance[0] - ref_variance[0]) + std::abs(variance[1] - ref_variance[1])) > 1.E-9))
            throw std::runtime_error("ERROR: mismatch in sampling between reference and computed DREAM.");
    }
}
