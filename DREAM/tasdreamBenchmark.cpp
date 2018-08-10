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

#ifndef __TASMANIAN_TASDREAM_BENCHMARKS_CPP
#define __TASMANIAN_TASDREAM_BENCHMARKS_CPP

#include "tasdreamBenchmark.hpp"

TasmanianSparseGrid* prepareGrid(int num_outputs, int depth, int mpi_me, int mpi_all){
    TasmanianSparseGrid *grid = new TasmanianSparseGrid();

    double delta_x = 1.0 / ((double) (num_outputs * mpi_all)), delta_x2 = delta_x / 2.0;

    grid->makeGlobalGrid(2, num_outputs, depth, type_iptotal, rule_clenshawcurtis);
    //grid->makeLocalPolynomialGrid(2, num_outputs, depth, 1, rule_localp);

    double x_root = ((double) mpi_me) / ((double) mpi_all);

    double a[2] = {1.0, 1.0};
    double b[2] = {5.0, 5.0};
    grid->setDomainTransform(a,b);

    double *grid_points = grid->getPoints();
    int num_sg_points = grid->getNumPoints();

    double *simulations = new double[num_sg_points * num_outputs];

    for(int i=0; i<num_sg_points; i++){
        double x = x_root + delta_x2; // * ((double) (2*j+1));
        for(int j=0; j<num_outputs; j++){
            simulations[i*num_outputs + j] = sin(grid_points[2*i] * M_PI * x) + sin(grid_points[2*i+1] * M_PI * x);
            x += delta_x;
        }
    }

    grid->loadNeededPoints(simulations);
    delete[] simulations;
    delete[] grid_points;

    return grid;
}
double* getData(int num_outputs, int mpi_me, int mpi_all){
    double delta_x = 1.0 / ((double) (num_outputs * mpi_all)), delta_x2 = delta_x / 2.0;
    double x = delta_x2 + ((double) mpi_me) / ((double) mpi_all);
    double *data = new double[num_outputs];
    for(int j=0; j<num_outputs; j++){
        data[j] = sin(2.0 * M_PI * x) + sin(3.0 * M_PI * x);
        x += delta_x;
    }
    return data;
}
void writeMatrix(const char *filename, int rows, int cols, const double mat[]){
    std::ofstream ofs;
    ofs.open(filename);
    ofs << rows << " " << cols << endl;
    ofs << std::scientific << std::setprecision(17);
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
            ofs << setw(25) << mat[i*cols + j] << " ";
        }
        ofs << endl;
    }
    ofs.close();
}

void sharedBenchmarkBasicAlpha(int num_outputs, int depth, int num_chains, int num_burnup, int num_mcmc, int gpuID, const char* outfilename){
    // num_outputs = 100; depth = 8; num_chains = 2000; num_burnup = 10; num_mcmc = 5;

//    miro@mantis32:~/RAM/TasBuild$ ./tasdream -bench basic-alpha 524288 36 1024 128 64
//    1409
//    Elapsed time: 53 (GPU 0)
//    miro@mantis32:~/RAM/TasBuild$ ./tasdream -bench basic-alpha 524288 36 1024 128 64
//    1409
//    Elapsed time: 419 (GPU 1)
//    miro@mantis32:~/RAM/TasBuild$ ./tasdream -bench basic-alpha 524288 36 1024 128 64
//    1409
//    Elapsed time: 2765
//    miro@mantis32:~/RAM/TasBuild$ ./tasdream -bench basic-alpha 65536 36 1024 128 64
//    1409
//    Elapsed time: 321
//    miro@mantis32:~/RAM/TasBuild$ ./tasdream -bench basic-alpha 65536 36 1024 128 64
//    1409
//    Elapsed time: 162 (GPU 0)
//    miro@mantis32:~/RAM/TasBuild$ ./tasdream -bench basic-alpha 65536 36 1024 128 64
//    1409
//    Elapsed time: 55 (GPU 1)


    TasmanianSparseGrid *grid = prepareGrid(num_outputs, depth);
    //cout << grid->getNumPoints() << endl;
    cout << "basic-alpha, out=" << num_outputs
                        << ", nodes=" << grid->getNumPoints()
                        << ", chains=" << num_chains
                        << ", burnup=" << num_burnup
                        << ", mcmc=" << num_mcmc
                        << ", gpu=";
    if (gpuID == -1){
        cout << "none";
    }else{
        cout << gpuID;
        //if (!grid->isCudaEnabled()) cout << " (cuda disabled)";
    }
    if (outfilename != 0){
        cout << ", outfile=" << endl;
    }else{
        cout << endl;
    }

    //if (grid->isCudaEnabled() && (gpuID > -1)){
    if (gpuID > -1){
        grid->enableAcceleration(accel_gpu_cublas);
        if (gpuID+1 > grid->getNumGPUs()){
            cout << "GPU " << gpuID << " is not available";
            gpuID = grid->getNumGPUs()-1;
            cout << ", using GPU " << gpuID << endl;
        }
        grid->setGPUID(gpuID);
    }

    PosteriorFromModel *post = new PosteriorFromModel(grid);

    double *data = getData(num_outputs);

    double scale = ((double) num_outputs) / 10.0;
    GaussianLikelihood *likely = new GaussianLikelihood(num_outputs, likely_gauss_scale, &scale, 1, data);
    post->setLikelihood(likely);

    TasmanianDREAM *dream = new TasmanianDREAM();
    dream->setProbabilityWeightFunction(post);

    dream->setNumChains(num_chains);
    GaussianPDF gauss(0.0, 0.01);
    dream->setCorrectionAll(&gauss);

    int start_time = (int) time(0);
    double *mcmc = dream->collectSamples(num_burnup, num_mcmc);
    int end_time = (int) time(0);

    cout << "Elapsed time: " << end_time - start_time << endl;

    if (outfilename != 0){
        writeMatrix(outfilename, num_mcmc*num_chains, grid->getNumDimensions(), mcmc);
    }

    delete[] mcmc;
    delete dream;
    delete likely;
    delete[] data;
    delete post;
    delete grid;
}


#ifdef MPI_VERSION
void mpiBenchmarkBasicAlpha(int num_outputs, int depth, int num_chains, int num_burnup, int num_mcmc, const char* outfilename){

    MPI_Init(NULL, NULL);

    int mpi_me, mpi_all;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_me);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_all);

    TasmanianSparseGrid *grid = prepareGrid(num_outputs, depth, mpi_me, mpi_all);
    if (mpi_me == 0){
        cout << "basic-alpha-mpi, out="    << num_outputs
                            << ", nodes="  << grid->getNumPoints()
                            << ", chains=" << num_chains
                            << ", burnup=" << num_burnup
                            << ", mcmc="   << num_mcmc
                            << ", mpi="    << mpi_all << endl;
    }

    PosteriorFromModel *post = new PosteriorFromModel(grid);

    double *data = getData(num_outputs, mpi_me, mpi_all);

    //ScaledGaussianLikelihoood *like = new ScaledGaussianLikelihoood(num_outputs, 10.0 / ((double) (num_outputs * mpi_all)));
    //post->setLikelihood(like);

    double scale = ((double) (num_outputs * mpi_all)) / 10.0;
    GaussianLikelihood *likely = new GaussianLikelihood(num_outputs, likely_gauss_scale, &scale, 1, data);
    post->setLikelihood(likely);

    DistributedPosteriorTSGModel *dist = new DistributedPosteriorTSGModel(MPI_COMM_WORLD, post);
    dist->setNumChanis(num_chains);
    bool useLogForm = true;

    if (mpi_me == 0){

        TasmanianDREAM dream;
        dream.setProbabilityWeightFunction(dist);

        dream.setNumChains(num_chains);
        GaussianPDF gauss(0.0, 0.01);
        dream.setCorrectionAll(&gauss);

        int start_time = (int) time(0);
        double *mcmc = dream.collectSamples(num_burnup, num_mcmc, useLogForm);
        int end_time = (int) time(0);

        cout << "Elapsed time: " << end_time - start_time << endl;

        if (outfilename != 0){
            writeMatrix(outfilename, num_mcmc*num_chains, grid->getNumDimensions(), mcmc);
        }

        delete[] mcmc;

    }else{
        dist->workerLoop(useLogForm); // main runs the chains, the rest wait in a loop only to evaluate the likelihood
    }

    dist->endWorkerLoop();


    delete dist;
    delete likely;
    delete[] data;
    delete post;
    delete grid;

    MPI_Finalize();

}
#endif


#endif // __TASMANIAN_TASDREAM_BENCHMARKS_CPP
