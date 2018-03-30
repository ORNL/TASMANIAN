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

#ifndef __TASMANIAN_DREAM_HPP
#define __TASMANIAN_DREAM_HPP

#include "TasmanianSparseGrid.hpp"

#include "tdrEnumerates.hpp"
#include "tdrCorePDF.hpp"

namespace TasDREAM{

class CustomModelWrapper{ // use this class for inheritance purposes only
public:
    CustomModelWrapper();
    virtual ~CustomModelWrapper();
    virtual int getNumDimensions() const = 0;
    virtual int getNumOutputs() const = 0;
    virtual void evaluate(const double x[], int num_points, double y[]) const = 0;
};


class ProbabilityWeightFunction{ // use this class for inheritance purposes only
public:
    ProbabilityWeightFunction();
    virtual ~ProbabilityWeightFunction();

    virtual int getNumDimensions() const = 0;

    virtual void evaluate(int num_points, const double x[], double y[], bool useLogForm) = 0;
    // in most cases evaluate should be const, but for caching purposes you may want it to be not a const function
    // WARNING: TasmanianDREAM class holds an alias to an instance of this class,
    //      modifying the class after calling setProbabilityWeightFunction could result in undefined behavior
    //      TasmanianDREAM class does not call delete on the pointer

    virtual void getDomainBounds(bool* lower_bound, bool* upper_bound) = 0;
    virtual void getDomainBounds(double* lower_bound, double* upper_bound) = 0;
    virtual void getInitialSample(double x[]) = 0;
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Link between Sparse Grids and DREAM modules
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class PosteriorFromModel : public ProbabilityWeightFunction{
public:
    PosteriorFromModel(const TasGrid::TasmanianSparseGrid *model, std::ostream *os = 0);
    PosteriorFromModel(const CustomModelWrapper *model, std::ostream *os = 0);
    ~PosteriorFromModel();
    void overwritePDF(int dimension, BasePDF* pdf);

    void setErrorLog(std::ostream *os);

    int getNumDimensions() const;

    void evaluate(int num_points, const double x[], double y[], bool useLogForm);

    void getInitialSample(double x[]);
    void setLikelihood(BaseLikelihood *likelihood);
    void setData(int num_data_samples, const double *posterior_data);

    void getDomainBounds(bool* lower_bound, bool* upper_bound);
    void getDomainBounds(double* lower_bound, double* upper_bound);

private:
    const TasGrid::TasmanianSparseGrid *grid;
    const CustomModelWrapper *cmodel;

    int num_dimensions, num_outputs;
    BasePDF **priors;
    bool *priors_created_here;

    int num_data;
    const double *data;

    BaseLikelihood *likely;

    int num_cache;
    double *model_cache;

    std::ostream *logstream;
};


#ifdef MPI_VERSION
// Assumes grids split by outputs, each chunk sits on a different node
// Assume likelihood multiply across nodes (so we communicate only scalars)
class DistributedPosteriorTSGModel : public ProbabilityWeightFunction{
public:
    DistributedPosteriorTSGModel(MPI_Comm in_comm, PosteriorFromModel *local_posterior, std::ostream *os = 0);
    ~DistributedPosteriorTSGModel();

    void setErrorLog(std::ostream *os);

    int getNumDimensions() const;
    void evaluate(int num_points, const double x[], double y[], bool useLogForm);

    void getInitialSample(double x[]);

    void setNumChanis(int num_dream_chains); // needed for MPI communication purposes

    void getDomainBounds(bool* lower_bound, bool* upper_bound);
    void getDomainBounds(double* lower_bound, double* upper_bound);

    void workerLoop(bool useLogForm);
    void endWorkerLoop();

private:
    PosteriorFromModel *posterior;

    MPI_Comm comm;
    int comm_me, comm_size;

    int num_dimensions, num_chains;

    std::ostream *logstream;
};
#endif // MPI_VERSION

// use a sparse grid as interpolated likelihood
class LikelihoodTSG : public ProbabilityWeightFunction{
public:
    LikelihoodTSG(const TasGrid::TasmanianSparseGrid *likely, bool savedLogForm, std::ostream *os = 0);
    ~LikelihoodTSG();
    void setPDF(int dimension, BasePDF* &pdf);

    void setErrorLog(std::ostream *os);

    int getNumDimensions() const;

    void evaluate(int num_points, const double x[], double y[], bool useLogForm);

    void getInitialSample(double x[]);

    void getDomainBounds(bool* lower_bound, bool* upper_bound);
    void getDomainBounds(double* lower_bound, double* upper_bound);

private:
    const TasGrid::TasmanianSparseGrid *grid;
    bool savedLogarithmForm;

    int num_dimensions;
    BasePDF **priors;

    std::ostream *logstream;
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  DREAM Samples
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class TasmanianDREAM{
public:
    TasmanianDREAM(std::ostream *os = 0);
    ~TasmanianDREAM();
    void overwriteBaseUnifrom(const BaseUniform *new_uniform);

    static const char* getVersion();
    static int getVersionMajor();
    static int getVersionMinor();
    static const char* getLicense();

    void setErrorLog(std::ostream *os);

    void setNumChains(int num_dream_chains);
    int getNumChains() const;

    void setJumpScale(double jump_scale);
    double getJumpScale();

    void setCorrectionAll(BasePDF *correct);
    void setCorrection(int dim, BasePDF *correct);
    const BasePDF* getCorrection(int dim);

    int getNumDimensions() const;

    // read/write chain state
    void setChainState(const double* state);
    //void clearPDFValues(); // delete cached pdf values, use when the state has been saved from a different pdf

    void setProbabilityWeightFunction(ProbabilityWeightFunction *probability_weight);
    // in most cases evaluate should be const, but for caching purposes you may want it to be not a const function
    // WARNING: TasmanianDREAM class holds an alias to an instance of the ProbabilityWeightFunction class,
    //      modifying the class after calling setProbabilityWeightFunction could result in undefined behavior
    //      TasmanianDREAM class does not call delete on the pointer

    double* collectSamples(int num_burnup, int num_samples, bool useLogForm = false);
    double* getPDFHistory() const;

protected:
    void clear();

    void advanceMCMCDREAM(bool useLogForm);

private:
    int num_dimensions, num_chains;
    ProbabilityWeightFunction *pdf;

    double jump;
    BasePDF **corrections;

    double *chain_state, *pdf_values, *old_state, *new_pdf_values;

    bool *isBoudnedBelow, *isBoudnedAbove;
    double *boundBelow, *boundAbove;

    const BaseUniform *core;
    CppUniformSampler unifrom_cpp;

    int num_pdf_history;
    double *pdf_history;

    std::ostream *logstream;
};

}

#endif
