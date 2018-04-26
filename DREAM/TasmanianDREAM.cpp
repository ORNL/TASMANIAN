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

#ifndef __TASMANIAN_DREAM_CPP
#define __TASMANIAN_DREAM_CPP

#include "TasmanianDREAM.hpp"

namespace TasDREAM{

CustomModelWrapper::CustomModelWrapper(){}
CustomModelWrapper::~CustomModelWrapper(){}
ProbabilityWeightFunction::ProbabilityWeightFunction(){}
ProbabilityWeightFunction::~ProbabilityWeightFunction(){}

PosteriorFromModel::PosteriorFromModel(const TasGrid::TasmanianSparseGrid *model, std::ostream *os) :
    grid(model), cmodel(0), num_dimensions(0), num_outputs(0), priors(0), priors_created_here(0),
    num_data(0), data(0), likely(0), num_cache(0), model_cache(0), logstream(os)
{
    #ifndef TASMANIAN_XSDK
    if (logstream == 0) logstream = &cerr;
    #endif // TASMANIAN_XSDK
    num_dimensions = grid->getNumDimensions();
    if (num_dimensions == 0){
        if (logstream != 0) (*logstream) << "ERROR: PosteriorFromModel cannot load a grid with no information" << endl;
    }
    num_outputs = grid->getNumOutputs();
    if (num_outputs < 1){
        if (logstream != 0) (*logstream) << "ERROR: cannot work with a grid with no outputs" << endl;
    }
    priors = new BasePDF*[num_dimensions];
    SparseGridDomainToPDF::assumeDefaultPDF(model, priors);
    priors_created_here = new bool[num_dimensions];
    for(int j=0; j<num_dimensions; j++) priors_created_here[j] = true;
}
PosteriorFromModel::PosteriorFromModel(const CustomModelWrapper *model, std::ostream *os) :
    grid(0), cmodel(model), num_dimensions(0), num_outputs(0), priors(0),  priors_created_here(0),
    num_data(0), data(0), likely(0), num_cache(0), model_cache(0), logstream(os)
{
    #ifndef TASMANIAN_XSDK
    if (logstream == 0) logstream = &cerr;
    #endif // TASMANIAN_XSDK
    num_dimensions = cmodel->getNumDimensions();
    if (num_dimensions == 0){
        if (logstream != 0) (*logstream) << "ERROR: PosteriorFromModel cannot load a model with no information" << endl;
    }
    num_outputs = cmodel->getNumOutputs();
    if (num_outputs < 1){
        if (logstream != 0) (*logstream) << "ERROR: cannot work with a model with no outputs" << endl;
    }
    priors = new BasePDF*[num_dimensions];
    for(int j=0; j<num_dimensions; j++) priors[j] = 0;
    priors_created_here = new bool[num_dimensions];
    for(int j=0; j<num_dimensions; j++) priors_created_here[j] = false;
}
PosteriorFromModel::~PosteriorFromModel(){
    if (priors != 0){
        for(int j=0; j<num_dimensions; j++) if ((priors[j] != 0) && priors_created_here[j]) delete priors[j];
        delete[] priors;
    }
    if (priors_created_here != 0) delete[] priors_created_here;
    if (model_cache != 0) delete[] model_cache;
}
void PosteriorFromModel::overwritePDF(int dimension, BasePDF* pdf){
    if ((dimension < 0) || (dimension >= num_dimensions)) if (logstream != 0) (*logstream) << "ERROR: attempt to overwritePDF for dimension outside of range" << endl;
    if ((priors[dimension] != 0) && priors_created_here[dimension]) delete priors[dimension];
    priors[dimension] = pdf;
    priors_created_here[dimension] = false;
    //pdf = 0;
}
void PosteriorFromModel::setErrorLog(std::ostream *os){ logstream = os; }
int PosteriorFromModel::getNumDimensions() const{ return num_dimensions; }

void PosteriorFromModel::evaluate(int num_points, const double x[], double y[], bool useLogForm){
    if (num_points != num_cache){
        if (model_cache != 0) delete[] model_cache;
        num_cache = num_points;
        model_cache = new double[num_cache * num_outputs];
    }

    if (grid != 0){
        grid->evaluateBatch(x, num_points, model_cache); // fastest
//        for(int i=0; i<num_points; i++){
//            //grid->evaluate(&(x[i*num_dimensions]), &(model_cache[i*num_outputs])); // slow, thread safe
//            grid->fastEvaluate(&(x[i*num_dimensions]), &(model_cache[i*num_outputs])); // faster, not thread safe
//        }
    }else{
        cmodel->evaluate(x, num_points, model_cache); // fastest
    }

    likely->getLikelihood(num_points, model_cache, num_data, data, y, useLogForm);

    if (useLogForm){
        for(int i=0; i<num_points; i++){
            for(int j=0; j<num_dimensions; j++) y[i] += priors[j]->getDensityLog(x[i*num_dimensions +j]);
        }
    }else{
        for(int i=0; i<num_points; i++){
            for(int j=0; j<num_dimensions; j++) y[i] *= priors[j]->getDensity(x[i*num_dimensions +j]);
        }
    }
}

void PosteriorFromModel::setLikelihood(BaseLikelihood *likelihood){ likely = likelihood; }
void PosteriorFromModel::setData(int num_data_samples, const double *posterior_data){
    num_data = num_data_samples;
    data = posterior_data;
}

void PosteriorFromModel::getInitialSample(double x[]){ for(int j=0; j<num_dimensions; j++) x[j] = priors[j]->getSample(); }

void PosteriorFromModel::getDomainBounds(bool* lower_bound, bool* upper_bound){
    for(int j=0; j<num_dimensions; j++){
        lower_bound[j] = priors[j]->isBoundedBelow();
        upper_bound[j] = priors[j]->isBoundedAbove();
    }
}
void PosteriorFromModel::getDomainBounds(double* lower_bound, double* upper_bound){
    for(int j=0; j<num_dimensions; j++){
        lower_bound[j] = priors[j]->getBoundBelow();
        upper_bound[j] = priors[j]->getBoundAbove();
    }
}


#ifdef MPI_VERSION
DistributedPosteriorTSGModel::DistributedPosteriorTSGModel(MPI_Comm in_comm, PosteriorFromModel *local_posterior, std::ostream *os) :
    posterior(local_posterior), comm(in_comm), num_chains(0), logstream(os){

    #ifndef TASMANIAN_XSDK
    if (logstream == 0) logstream = &cerr;
    #endif // TASMANIAN_XSDK

    MPI_Comm_rank(comm, &comm_me);
    MPI_Comm_size(comm, &comm_size);

    num_dimensions = local_posterior->getNumDimensions();
}
DistributedPosteriorTSGModel::~DistributedPosteriorTSGModel(){}

void DistributedPosteriorTSGModel::setErrorLog(std::ostream *os){ logstream = os; }

int DistributedPosteriorTSGModel::getNumDimensions() const{ return posterior->getNumDimensions(); }

void DistributedPosteriorTSGModel::evaluate(int num_points, const double x[], double y[], bool useLogForm){
    // MPI witchcraft
    double *local_y = new double[num_chains];
    double *x_extended = new double[num_chains*num_dimensions + 2]; std::copy(x, x+num_dimensions*num_points, x_extended);
    x_extended[num_dimensions*num_chains] = 1.0; // send message: keep working
    x_extended[num_dimensions*num_chains+1] = (double) num_points; // send message: with number of points
    MPI_Bcast((void*) x_extended, num_dimensions*num_chains+2, MPI_DOUBLE, 0, comm);

    posterior->evaluate(num_points, x, local_y, useLogForm);

    for(int i=num_points; i<num_chains; i++) local_y[i] = 0.0;

    if (useLogForm){
        MPI_Reduce(local_y, y, num_points, MPI_DOUBLE, MPI_SUM, 0, comm);
    }else{
        MPI_Reduce(local_y, y, num_points, MPI_DOUBLE, MPI_PROD, 0, comm);
    }
}

void DistributedPosteriorTSGModel::getInitialSample(double x[]){ posterior->getInitialSample(x); }

void DistributedPosteriorTSGModel::setNumChanis(int num_dream_chains){ num_chains = num_dream_chains; }

void DistributedPosteriorTSGModel::getDomainBounds(bool* lower_bound, bool* upper_bound){
    posterior->getDomainBounds(lower_bound, upper_bound);
}
void DistributedPosteriorTSGModel::getDomainBounds(double* lower_bound, double* upper_bound){
    posterior->getDomainBounds(lower_bound, upper_bound);
}

void DistributedPosteriorTSGModel::workerLoop(bool useLogForm){
    double *x = new double[num_dimensions*num_chains+2];
    double *local_y = new double[num_chains];

    bool keep_working = true;

    while(keep_working){
        MPI_Bcast((void*) x, num_dimensions*num_chains+2, MPI_DOUBLE, 0, comm);

        if (x[num_dimensions*num_chains] == 1.0){ // received message: keep working

            int num_points = (int) x[num_dimensions*num_chains+1];

            posterior->evaluate(num_points, x, local_y, useLogForm);

            for(int i=num_points; i<num_chains; i++) local_y[i] = 0.0;

            if (useLogForm){
                MPI_Reduce(local_y, 0, num_points, MPI_DOUBLE, MPI_SUM, 0, comm);
            }else{
                MPI_Reduce(local_y, 0, num_points, MPI_DOUBLE, MPI_PROD, 0, comm);
            }
        }else{
            keep_working = false;
        }
    }
    delete[] local_y;
    delete[] x;
}
void DistributedPosteriorTSGModel::endWorkerLoop(){
    if (comm_me == 0){
        double *x = new double[num_dimensions*num_chains+2];
        std::fill(x, x + num_dimensions*num_chains+2, 0.0); // send message: stop working
        MPI_Bcast((void*) x, num_dimensions*num_chains+2, MPI_DOUBLE, 0, comm);
        delete[] x;
    }
}
#endif // MPI_VERSION


LikelihoodTSG::LikelihoodTSG(const TasGrid::TasmanianSparseGrid *likely, bool savedLogForm, std::ostream *os) :
    grid(likely), savedLogarithmForm(savedLogForm), num_dimensions(0), priors(0), logstream(os)
{
    #ifndef TASMANIAN_XSDK
    if (logstream == 0) logstream = &cerr;
    #endif // TASMANIAN_XSDK
    num_dimensions = grid->getNumDimensions();
    if (num_dimensions < 1){
        if (logstream != 0) (*logstream) << "ERROR: cannot work with an empty Grid" << endl;
    }
    priors = new BasePDF*[num_dimensions];
    SparseGridDomainToPDF::assumeDefaultPDF(likely, priors);
}
LikelihoodTSG::~LikelihoodTSG(){
    if (priors != 0){
        for(int i=0; i<num_dimensions; i++) delete priors[i];
        delete[] priors;
    }
}
void LikelihoodTSG::setPDF(int dimension, BasePDF* &pdf){
    if ((dimension < 0) || (dimension >= num_dimensions)) return;
    if (priors[dimension] != 0) delete priors[dimension];
    priors[dimension] = pdf;
    pdf = 0;
}
void LikelihoodTSG::setErrorLog(std::ostream *os){ logstream = os; }
int LikelihoodTSG::getNumDimensions() const{ return num_dimensions; }

void LikelihoodTSG::evaluate(int num_points, const double x[], double y[], bool useLogForm){

    grid->evaluateBatch(x, num_points, y); // fastest
//    for(int i=0; i<num_points; i++){
//        //grid->evaluate(&(x[i*num_dimensions]), &(y[i])); // slow, thread safe
//        grid->fastEvaluate(&(x[i*num_dimensions]), &(y[i])); // faster, not thread safe
//    }
    if (savedLogarithmForm ^ useLogForm){ // using form other than the one saved
        if (savedLogarithmForm){
            for(int i=0; i<num_points; i++) y[i] = exp(y[i]);
        }else{
            for(int i=0; i<num_points; i++) y[i] = log(y[i]);
        }
    }

    if (useLogForm){
        for(int i=0; i<num_points; i++){
            for(int j=0; j<num_dimensions; j++) y[i] += priors[j]->getDensityLog(x[i*num_dimensions +j]);
        }
    }else{
        for(int i=0; i<num_points; i++){
            for(int j=0; j<num_dimensions; j++) y[i] *= priors[j]->getDensity(x[i*num_dimensions +j]);
        }
    }
}

void LikelihoodTSG::getInitialSample(double x[]){ for(int j=0; j<num_dimensions; j++) x[j] = priors[j]->getSample(); }

void LikelihoodTSG::getDomainBounds(bool* lower_bound, bool* upper_bound){
    for(int j=0; j<num_dimensions; j++){
        lower_bound[j] = priors[j]->isBoundedBelow();
        upper_bound[j] = priors[j]->isBoundedAbove();
    }
}
void LikelihoodTSG::getDomainBounds(double* lower_bound, double* upper_bound){
    for(int j=0; j<num_dimensions; j++){
        lower_bound[j] = priors[j]->getBoundBelow();
        upper_bound[j] = priors[j]->getBoundAbove();
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Actual DREAM Sampler
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TasmanianDREAM::TasmanianDREAM(std::ostream *os): num_dimensions(-1), num_chains(-1),
    pdf(0), jump(1.0), corrections(0), chain_state(0), pdf_values(0), old_state(0), new_pdf_values(0),
    isBoudnedBelow(0), isBoudnedAbove(0), boundBelow(0), boundAbove(0), num_pdf_history(0),
    pdf_history(0), logstream(os)
{
    #ifndef TASMANIAN_XSDK
    if (logstream == 0) logstream = &cerr;
    #endif // TASMANIAN_XSDK
    core = &unifrom_cpp;
}
TasmanianDREAM::~TasmanianDREAM(){
    num_dimensions = -1;
    clear();
}
void TasmanianDREAM::overwriteBaseUnifrom(const BaseUniform *new_uniform){
    core = new_uniform;
}

const char* TasmanianDREAM::getVersion(){ return TASMANIAN_VERSION_STRING; }
const char* TasmanianDREAM::getLicense(){ return TASMANIAN_LICENSE; }
int TasmanianDREAM::getVersionMajor(){ return TASMANIAN_VERSION_MAJOR; }
int TasmanianDREAM::getVersionMinor(){ return TASMANIAN_VERSION_MINOR; }

void TasmanianDREAM::clear(){
    if (corrections != 0){
        for(int j=0; j<num_dimensions; j++) corrections[j] = 0;
        delete[] corrections;
        corrections = 0;
    }
    num_dimensions = -1;
    num_chains = -1;
    pdf = 0;
    jump = 0.0;
    if (chain_state != 0){ delete[] chain_state; chain_state = 0; }
    if (pdf_values != 0){ delete[] pdf_values; pdf_values = 0; }
    if (isBoudnedBelow != 0){ delete[] isBoudnedBelow; isBoudnedBelow = 0; }
    if (isBoudnedAbove != 0){ delete[] isBoudnedAbove; isBoudnedAbove = 0; }
    if (boundBelow != 0){ delete[] boundBelow; boundBelow = 0; }
    if (boundAbove != 0){ delete[] boundAbove; boundAbove = 0; }
    if (pdf_history != 0){ delete[] pdf_history; pdf_history = 0; }
}

void TasmanianDREAM::setErrorLog(std::ostream *os){ logstream = os; }

void TasmanianDREAM::setNumChains(int num_dream_chains){
    num_chains = num_dream_chains;
    if (chain_state != 0){ delete[] chain_state; chain_state = 0; }
    if (pdf_values != 0){ delete[] pdf_values; pdf_values = 0; }
    if (old_state != 0){ delete[] old_state; old_state = 0; }
    if (new_pdf_values != 0){ delete[] new_pdf_values; new_pdf_values = 0; }
}
int TasmanianDREAM::getNumChains() const{ return num_chains; }
int TasmanianDREAM::getNumDimensions() const{ return num_dimensions; }

void TasmanianDREAM::setJumpScale(double jump_scale){ jump = jump_scale; }
double TasmanianDREAM::getJumpScale(){ return jump; }

void TasmanianDREAM::setCorrectionAll(BasePDF *correct){
    if (num_dimensions == 0) if (logstream != 0) (*logstream) << "ERROR: cannot set correction before the pdf" << endl;
    for(int j=0; j<num_dimensions; j++) corrections[j] = correct;
}
void TasmanianDREAM::setCorrection(int dim, BasePDF *correct){
    if (num_dimensions == 0) if (logstream != 0) (*logstream) << "ERROR: cannot set correction before the pdf" << endl;
    if ((dim < 0) || (dim >= num_dimensions)) if (logstream != 0) (*logstream) << "ERROR: incorrect dimension" << endl;
    corrections[dim] = correct;
}
const BasePDF* TasmanianDREAM::getCorrection(int dim){
    if ((dim < 0) || (dim >= num_dimensions)){
        if (logstream != 0) (*logstream) << "ERROR: incorrect dimension" << endl;
        return 0;
    }
    return corrections[dim];
}

void TasmanianDREAM::setProbabilityWeightFunction(ProbabilityWeightFunction *probability_weight){
    clear();
    pdf = probability_weight;
    num_dimensions = probability_weight->getNumDimensions();
    isBoudnedBelow = new bool[num_dimensions];
    isBoudnedAbove = new bool[num_dimensions];
    boundBelow = new double[num_dimensions];
    boundAbove = new double[num_dimensions];
    probability_weight->getDomainBounds(isBoudnedBelow, isBoudnedAbove);
    probability_weight->getDomainBounds(boundBelow, boundAbove);
    corrections = new BasePDF*[num_dimensions];
    for(int j=0; j<num_dimensions; j++) corrections[j] = 0;
}

double* TasmanianDREAM::collectSamples(int num_burnup, int num_samples, bool useLogForm){
    if (num_chains < 1){ if (logstream != 0) (*logstream) << "No chains specified, cannot collect samples" << endl; return 0; } // no chains specified
    if (chain_state == 0){
        chain_state = new double[num_chains * num_dimensions];
        for(int j=0; j<num_chains; j++){
            pdf->getInitialSample(&(chain_state[j*num_dimensions]));
        }
    }
    if (pdf_values == 0){
        pdf_values = new double[num_chains];
        pdf->evaluate(num_chains, chain_state, pdf_values, useLogForm);
    }

    old_state = new double[num_chains * num_dimensions];
    new_pdf_values = new double[num_chains];

    for(int i=0; i<num_burnup; i++){
        advanceMCMCDREAM(useLogForm);
    }

    double *samples = new double[num_samples * num_dimensions * num_chains];
    if (pdf_history != 0){ delete[] pdf_history; }
    pdf_history = new double[num_samples * num_chains];
    num_pdf_history = num_samples * num_chains;

    for(int i=0; i<num_samples; i++){
        advanceMCMCDREAM(useLogForm);
        std::copy(chain_state, chain_state + num_dimensions*num_chains, &(samples[i*num_dimensions*num_chains]));
        std::copy(pdf_values, pdf_values + num_chains, &(pdf_history[i*num_chains]));
    }

    delete[] new_pdf_values;
    new_pdf_values = 0;
    delete[] old_state;
    old_state = 0;

    return samples;
}

void TasmanianDREAM::advanceMCMCDREAM(bool useLogForm){
    std::copy(chain_state, chain_state + num_chains*num_dimensions, old_state);
    bool *valid = new bool[num_chains];
    int num_need_evaluation = num_chains;

    bool allValid = true; //, savedGaussian = false;
    double unilength = 1.0 / ((double) num_chains);

    for(int i=0; i<num_chains; i++){
        int index1 = (int) (core->getSample01() / unilength);
        int index2 = (int) (core->getSample01() / unilength);
        if (index1 >= num_chains) index1 = num_chains; // this is needed in case core->getSample01() returns 1.0
        if (index2 >= num_chains) index2 = num_chains;


        valid[i] = true;

        for(int j=0; (j<num_dimensions) && valid[i]; j++){
            chain_state[i*num_dimensions + j] += old_state[index1*num_dimensions + j] - old_state[index2*num_dimensions + j];
            if (corrections[j] != 0) chain_state[i*num_dimensions + j] += corrections[j]->getSample();
            if ((isBoudnedBelow[j] && (chain_state[i*num_dimensions + j] < boundBelow[j]))
                || (isBoudnedAbove[j] && (chain_state[i*num_dimensions + j] > boundAbove[j]))){
                valid[i] = false;
                allValid = false;
                num_need_evaluation--;
            }
        }
    }

    // extract only the valid entries to evaluate
    double *need_evaluation;
    if (allValid){
        num_need_evaluation = num_chains;
        need_evaluation = chain_state;
    }else{
        need_evaluation = new double[num_dimensions * num_need_evaluation];
        num_need_evaluation = 0;
        for(int i=0; i<num_chains; i++){
            if (valid[i])
                std::copy(&(chain_state[i*num_dimensions]), &(chain_state[i*num_dimensions]) + num_dimensions, &(need_evaluation[num_dimensions * num_need_evaluation++]));
        }
    }

    pdf->evaluate(num_need_evaluation, need_evaluation, new_pdf_values, useLogForm);

    // clean memory and reorder the pdf values putting 0 in the invalid spots (for log case 0 should be -infty, hence using valid for accept/reject too
    if (!allValid){
        delete[] need_evaluation;

        int c = num_chains - num_need_evaluation;
        for(int i=num_chains-1; (i>=0) && (c>0) ; i--){
            if (valid[i]){
                new_pdf_values[i] = new_pdf_values[i-c];
            }else{
                new_pdf_values[i] = 0.0;
                c--;
            }
        }
    }

    // accept/reject
    for(int i=0; i<num_chains; i++){
        int decided = (valid[i]) ? -1 : 0;
        if (decided == -1){
            if (new_pdf_values[i] > pdf_values[i]){
                decided = 1;
            }else{
                if (useLogForm){
                    decided = ((new_pdf_values[i] - pdf_values[i]) > log(core->getSample01())) ? 1 : 0;
                }else{
                    decided = ((new_pdf_values[i] / pdf_values[i]) > core->getSample01()) ? 1 : 0;
                }
            }
        }

        if (decided == 1){
            pdf_values[i] = new_pdf_values[i];
        }else{
            std::copy(&(old_state[i*num_dimensions]), &(old_state[i*num_dimensions]) + num_dimensions, &(chain_state[i*num_dimensions]));
        }
    }

    delete[] valid;
}

double* TasmanianDREAM::getPDFHistory() const{
    double *hist = new double[num_pdf_history];
    std::copy(pdf_history, pdf_history + num_pdf_history, hist);
    return hist;
}

}

#endif
