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

void CustomModelWrapper::evaluate(const double*, int, double*) const{} // kept for backwards compatibility
void CustomModelWrapper::evaluate(const std::vector<double> &, std::vector<double> &y) const{ y.clear(); }

void ProbabilityWeightFunction::evaluate(int, const double*, double*, bool){} // kept for backwards compatibility
void ProbabilityWeightFunction::evaluate(const std::vector<double> &, std::vector<double> &y, bool){ y.clear(); }

void ProbabilityWeightFunction::getDomainBounds(bool*, bool*){} // kept for backwards compatibility
void ProbabilityWeightFunction::getDomainBounds(double*, double*){} // kept for backwards compatibility
void ProbabilityWeightFunction::getDomainBounds(std::vector<bool> &lower, std::vector<bool> &upper){ lower.clear(); upper.clear(); }
void ProbabilityWeightFunction::getDomainBounds(std::vector<double> &lower, std::vector<double> &upper){ lower.clear(); upper.clear(); }

PosteriorFromModel::PosteriorFromModel(const TasGrid::TasmanianSparseGrid *model) :
    grid(model), cmodel(0), num_dimensions(0), num_outputs(0),
    num_data(0), data(0), likely(0)
{
    num_dimensions = grid->getNumDimensions();
    if (num_dimensions == 0) throw std::runtime_error("ERROR: PosteriorFromModel() cannot load a grid with no dimension");
    num_outputs = grid->getNumOutputs();
    if (num_outputs < 1) throw std::runtime_error("ERROR:PosteriorFromModel() cannot work with a grid with no outputs");

    SparseGridDomainToPDF::assumeDefaultPDF(model, internal_priors);
    active_priors = internal_priors; // copy assignment
}
PosteriorFromModel::PosteriorFromModel(const CustomModelWrapper *model) :
    grid(0), cmodel(model), num_dimensions(0), num_outputs(0),
    num_data(0), data(0), likely(0)
{
    num_dimensions = cmodel->getNumDimensions();
    if (num_dimensions < 1) throw std::runtime_error("ERROR: PosteriorFromModel() cannot load a model with dimension < 1");
    num_outputs = cmodel->getNumOutputs();
    if (num_outputs < 1) throw std::runtime_error("ERROR:PosteriorFromModel() cannot work with a model with no outputs");

    internal_priors.resize(num_dimensions, 0);
    active_priors.resize(num_dimensions, 0);
}
PosteriorFromModel::~PosteriorFromModel(){
    for(auto &p : internal_priors) if (p != 0) delete p;
}
void PosteriorFromModel::overwritePDF(int dimension, BasePDF* pdf){
    if ((dimension < 0) || (dimension >= num_dimensions)) throw std::invalid_argument("ERROR: overwritePDF() called with incorrect dimension");
    active_priors[dimension] = pdf;
}
int PosteriorFromModel::getNumDimensions() const{ return num_dimensions; }

void PosteriorFromModel::evaluate(const std::vector<double> &x, std::vector<double> &y, bool useLogForm){
    size_t num_points = x.size() / num_dimensions;

    std::vector<double> model_output;
    if (grid != 0){
        grid->evaluateBatch(x, model_output); // fastest
    }else{
        cmodel->evaluate(x, model_output); // default vector API
        if ((num_points > 0) && (model_output.size() == 0)){
            model_output.resize(num_points * num_outputs);
            cmodel->evaluate(x.data(), (int) num_points, model_output.data());
        }
    }

    likely->getLikelihood((int) num_points, model_output.data(), y, num_data, data, useLogForm);

    if (useLogForm){
        for(size_t i=0; i<num_points; i++){
            for(int j=0; j<num_dimensions; j++) y[i] += active_priors[j]->getDensityLog(x[i*num_dimensions +j]);
        }
    }else{
        for(size_t i=0; i<num_points; i++){
            for(int j=0; j<num_dimensions; j++) y[i] *= active_priors[j]->getDensity(x[i*num_dimensions +j]);
        }
    }
}

void PosteriorFromModel::setLikelihood(BaseLikelihood *likelihood){ likely = likelihood; }
void PosteriorFromModel::setData(int num_data_samples, const double *posterior_data){
    num_data = num_data_samples;
    data = posterior_data;
}

void PosteriorFromModel::getInitialSample(double x[]){ for(int j=0; j<num_dimensions; j++) x[j] = active_priors[j]->getSample(); }

void PosteriorFromModel::getDomainBounds(std::vector<bool> &lower, std::vector<bool> &upper){
    lower.resize(num_dimensions);
    upper.resize(num_dimensions);
    for(int j=0; j<num_dimensions; j++){
        lower[j] = active_priors[j]->isBoundedBelow();
        upper[j] = active_priors[j]->isBoundedAbove();
    }
}
void PosteriorFromModel::getDomainBounds(std::vector<double> &lower, std::vector<double> &upper){
    lower.resize(num_dimensions);
    upper.resize(num_dimensions);
    for(int j=0; j<num_dimensions; j++){
        lower[j] = active_priors[j]->getBoundBelow();
        upper[j] = active_priors[j]->getBoundAbove();
    }
}


#ifdef MPI_VERSION
DistributedPosteriorTSGModel::DistributedPosteriorTSGModel(MPI_Comm in_comm, PosteriorFromModel *local_posterior) :
    posterior(local_posterior), comm(in_comm), num_chains(0){
    MPI_Comm_rank(comm, &comm_me);
    MPI_Comm_size(comm, &comm_size);

    num_dimensions = local_posterior->getNumDimensions();
}
DistributedPosteriorTSGModel::~DistributedPosteriorTSGModel(){}

int DistributedPosteriorTSGModel::getNumDimensions() const{ return posterior->getNumDimensions(); }

void DistributedPosteriorTSGModel::evaluate(const std::vector<double> &x, std::vector<double> &y, bool useLogForm){
    int num_points = x.size() / num_dimensions;
    // MPI witchcraft
    std::vector<double> local_y(num_points);
    std::vector<double> x_extended(num_chains*num_dimensions + 2);
    std::copy(x.begin(), x.end(), x_extended.data());
    x_extended[num_dimensions*num_chains] = 1.0; // send message: keep working
    x_extended[num_dimensions*num_chains+1] = (double) num_points; // send message: with number of points
    MPI_Bcast((void*) x_extended.data(), num_dimensions*num_chains+2, MPI_DOUBLE, 0, comm);

    posterior->evaluate(x, local_y, useLogForm);

    y.resize(num_points);
    MPI_Reduce(local_y.data(), y.data(), num_points, MPI_DOUBLE, (useLogForm) ? MPI_SUM : MPI_PROD, 0, comm);
}

void DistributedPosteriorTSGModel::getInitialSample(double x[]){ posterior->getInitialSample(x); }

void DistributedPosteriorTSGModel::setNumChanis(int num_dream_chains){ num_chains = num_dream_chains; }

void DistributedPosteriorTSGModel::getDomainBounds(std::vector<bool> &lower, std::vector<bool> &upper){
    posterior->getDomainBounds(lower, upper);
}
void DistributedPosteriorTSGModel::getDomainBounds(std::vector<double> &lower, std::vector<double> &upper){
    posterior->getDomainBounds(lower, upper);
}

void DistributedPosteriorTSGModel::workerLoop(bool useLogForm){
    std::vector<double> x(num_dimensions*num_chains+2);
    std::vector<double> local_y(num_chains);

    bool keep_working = true;

    while(keep_working){
        MPI_Bcast((void*) x.data(), num_dimensions*num_chains+2, MPI_DOUBLE, 0, comm);

        if (x[num_dimensions*num_chains] == 1.0){ // received message: keep working

            int num_points = (int) x[num_dimensions*num_chains+1];

            posterior->evaluate(x, local_y, useLogForm);

            MPI_Reduce(local_y.data(), 0, num_points, MPI_DOUBLE,((useLogForm) ? MPI_SUM : MPI_PROD), 0, comm);
        }else{
            keep_working = false;
        }
    }
}
void DistributedPosteriorTSGModel::endWorkerLoop(){
    if (comm_me == 0){
        std::vector<double> x(num_dimensions*num_chains+2, 0.0); // send message: stop working
        MPI_Bcast((void*) x.data(), num_dimensions*num_chains+2, MPI_DOUBLE, 0, comm);
    }
}
#endif // MPI_VERSION


LikelihoodTSG::LikelihoodTSG(const TasGrid::TasmanianSparseGrid *likely, bool savedLogForm) :
    grid(likely), savedLogarithmForm(savedLogForm), num_dimensions(0)
{
    num_dimensions = grid->getNumDimensions();
    if (num_dimensions < 1) throw std::runtime_error("ERROR: in LikelihoodTSG() the likelihood specified has dimension less than 1");
    SparseGridDomainToPDF::assumeDefaultPDF(likely, internal_priors);
    active_priors = internal_priors;
}
LikelihoodTSG::~LikelihoodTSG(){
    for(auto &p : internal_priors) delete p;
}
void LikelihoodTSG::setPDF(int dimension, BasePDF* pdf){
    if ((dimension < 0) || (dimension >= num_dimensions)) return;
    active_priors[dimension] = pdf;
}
int LikelihoodTSG::getNumDimensions() const{ return num_dimensions; }

void LikelihoodTSG::evaluate(const std::vector<double> &x, std::vector<double> &y, bool useLogForm){
    size_t num_points = x.size() / num_dimensions;

    grid->evaluateBatch(x, y); // fastest

    if (savedLogarithmForm ^ useLogForm){ // using form other than the one saved
        if (savedLogarithmForm){
            for(size_t i=0; i<num_points; i++) y[i] = exp(y[i]);
        }else{
            for(size_t i=0; i<num_points; i++) y[i] = log(y[i]);
        }
    }

    if (useLogForm){
        for(size_t i=0; i<num_points; i++){
            for(int j=0; j<num_dimensions; j++) y[i] += active_priors[j]->getDensityLog(x[i*num_dimensions +j]);
        }
    }else{
        for(size_t i=0; i<num_points; i++){
            for(int j=0; j<num_dimensions; j++) y[i] *= active_priors[j]->getDensity(x[i*num_dimensions +j]);
        }
    }
}

void LikelihoodTSG::getInitialSample(double x[]){ for(int j=0; j<num_dimensions; j++) x[j] = active_priors[j]->getSample(); }

void LikelihoodTSG::getDomainBounds(std::vector<bool> &lower, std::vector<bool> &upper){
    lower.resize(num_dimensions);
    upper.resize(num_dimensions);
    for(int j=0; j<num_dimensions; j++){
        lower[j] = active_priors[j]->isBoundedBelow();
        upper[j] = active_priors[j]->isBoundedAbove();
    }
}
void LikelihoodTSG::getDomainBounds(std::vector<double> &lower, std::vector<double> &upper){
    lower.resize(num_dimensions);
    upper.resize(num_dimensions);
    for(int j=0; j<num_dimensions; j++){
        lower[j] = active_priors[j]->getBoundBelow();
        upper[j] = active_priors[j]->getBoundAbove();
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Actual DREAM Sampler
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TasmanianDREAM::TasmanianDREAM(): num_dimensions(-1), num_chains(-1),
    pdf(0), jump(1.0), state_initialized(false), values_initialized(false), values_logform(false)
{
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
    corrections.resize(0); // ensures I don't keep old unwanted pointers
    state_initialized = false;
    values_initialized = false;
    num_dimensions = -1;
    num_chains = -1;
    pdf = 0;
    jump = 0.0;
}

void TasmanianDREAM::setNumChains(int num_dream_chains){
    num_chains = num_dream_chains;
    state_initialized = false;
    values_initialized = false;
}
int TasmanianDREAM::getNumChains() const{ return num_chains; }
int TasmanianDREAM::getNumDimensions() const{ return num_dimensions; }

void TasmanianDREAM::setJumpScale(double jump_scale){ jump = jump_scale; }
double TasmanianDREAM::getJumpScale(){ return jump; }

void TasmanianDREAM::setCorrectionAll(BasePDF *correct){
    if (num_dimensions == 0) throw std::runtime_error("ERROR: cannot setCorrectionAll() before setProbabilityWeightFunction()");
    for(auto &p : corrections) p = correct;
}
void TasmanianDREAM::setCorrection(int dim, BasePDF *correct){
    if (num_dimensions == 0) throw std::runtime_error("ERROR: cannot setCorrection() before the setProbabilityWeightFunction()");
    if ((dim < 0) || (dim >= num_dimensions)) throw std::invalid_argument("ERROR: setCorrection() called with incorrect dimension");
    corrections[dim] = correct;
}
const BasePDF* TasmanianDREAM::getCorrection(int dim){
    if ((dim < 0) || (dim >= num_dimensions)){
        throw std::invalid_argument("ERROR: getCorrection() called with incorrect dimension");
    }
    return corrections[dim];
}

void TasmanianDREAM::setProbabilityWeightFunction(ProbabilityWeightFunction *probability_weight){
    clear();
    pdf = probability_weight;
    num_dimensions = probability_weight->getNumDimensions();

    corrections.resize(num_dimensions, 0); // clear sets this vector to zero size

    isBoudnedBelow.resize(num_dimensions);
    isBoudnedAbove.resize(num_dimensions);
    boundBelow.resize(num_dimensions);
    boundAbove.resize(num_dimensions);

    probability_weight->getDomainBounds(isBoudnedBelow, isBoudnedAbove);
    probability_weight->getDomainBounds(boundBelow, boundAbove);

    if (boundBelow.size() == 0){ // vector API not present, use the array version
        isBoudnedBelow.resize(num_dimensions);
        isBoudnedAbove.resize(num_dimensions);
        boundBelow.resize(num_dimensions);
        boundAbove.resize(num_dimensions);

        // std::vector<bool> is a special class without .data() member
        bool *bbelow = new bool[num_dimensions];
        bool *babove = new bool[num_dimensions];
        probability_weight->getDomainBounds(bbelow, babove);
        for(int i=0; i<num_dimensions; i++){
            isBoudnedBelow[i] = bbelow[i];
            isBoudnedAbove[i] = babove[i];
        }
        delete[] bbelow;
        delete[] babove;

        probability_weight->getDomainBounds(boundBelow.data(), boundAbove.data());
    }
}

void TasmanianDREAM::setChainState(const double* state){
    chain_state.resize(num_dimensions * num_chains);
    std::copy(state, state + num_dimensions * num_chains, chain_state.data());
    state_initialized = true;
    values_initialized = false;
}
void TasmanianDREAM::setChainState(const std::vector<double> &state){
    chain_state = state; // copy assignment
    if (chain_state.size() != (size_t) (num_dimensions * num_chains)) num_chains = (int) (chain_state.size() / num_dimensions);
    state_initialized = true;
    values_initialized = false;
}

void TasmanianDREAM::collectSamples(int num_burnup, int num_samples, double *samples, bool useLogForm){
    if (num_chains < 1) throw std::runtime_error("ERROR: must call setNumChains() before collectSamples()");
    if (!state_initialized){
        chain_state.resize(num_chains * num_dimensions);
        for(int j=0; j<num_chains; j++){
            pdf->getInitialSample(&(chain_state[j*num_dimensions]));
        }
        state_initialized = true;
    }
    if (!values_initialized || (values_logform != useLogForm)){ // if log form is changed mid evaluations
        pdf_values.clear();
        pdf->evaluate(chain_state, pdf_values, useLogForm);
        if (pdf_values.size() == 0){
            pdf_values.resize(num_dimensions * num_chains);
            pdf->evaluate(num_chains, chain_state.data(), pdf_values.data(), useLogForm);
        }
        values_initialized = true;
        values_logform = useLogForm;
    }

    for(int i=0; i<num_burnup; i++){
        advanceMCMCDREAM(useLogForm);
    }

    pdf_history.resize(num_samples * num_chains);

    auto iter_hist = pdf_history.begin();
    for(int i=0; i<num_samples; i++){
        advanceMCMCDREAM(useLogForm);
        std::copy(chain_state.begin(), chain_state.end(), &(samples[i*num_dimensions*num_chains]));
        std::copy(pdf_values.begin(), pdf_values.end(), iter_hist);
        advance(iter_hist, num_chains);
    }
}
double* TasmanianDREAM::collectSamples(int num_burnup, int num_samples, bool useLogForm){
    if (num_chains < 1) throw std::runtime_error("ERROR: must call setNumChains() before collectSamples()");
    double *samples = new double[num_samples * num_dimensions * num_chains];
    collectSamples(num_burnup, num_samples, samples, useLogForm);
    return samples;
}
void TasmanianDREAM::collectSamples(int num_burnup, int num_samples, std::vector<double> &samples, bool useLogForm){
    if (num_chains < 1) throw std::runtime_error("ERROR: must call setNumChains() before collectSamples()");
    samples.resize(num_samples * num_dimensions * num_chains);
    collectSamples(num_burnup, num_samples, samples.data(), useLogForm);
}

void TasmanianDREAM::advanceMCMCDREAM(bool useLogForm){
    std::vector<double> old_state = chain_state; // copy assignment
    std::vector<bool> valid(num_chains);
    int num_need_evaluation = num_chains;

    bool allValid = true; //, savedGaussian = false;
    double unilength = ((double) num_chains);

    for(int i=0; i<num_chains; i++){
        int index1 = (int) (core->getSample01() * unilength);
        int index2 = (int) (core->getSample01() * unilength);
        if (index1 >= num_chains) index1 = num_chains - 1; // this is needed in case core->getSample01() returns 1.0
        if (index2 >= num_chains) index2 = num_chains - 1;

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
    std::vector<double>* need_evaluation;
    std::vector<double> selected_for_evaluation;
    if (allValid){
        num_need_evaluation = num_chains;
        need_evaluation = &chain_state;
    }else{
        selected_for_evaluation.resize(num_dimensions * num_need_evaluation);
        num_need_evaluation = 0;
        for(int i=0; i<num_chains; i++){
            if (valid[i])
                std::copy(&(chain_state[i*num_dimensions]), &(chain_state[i*num_dimensions]) + num_dimensions, &(selected_for_evaluation[num_dimensions * num_need_evaluation++]));
        }
        need_evaluation = &selected_for_evaluation;
    }

    std::vector<double> new_pdf_values(num_chains);
    std::vector<double> computed_values;

    pdf->evaluate(*need_evaluation, computed_values, useLogForm);
    if (computed_values.size() == 0){
        pdf->evaluate(num_need_evaluation, need_evaluation->data(), new_pdf_values.data(), useLogForm);
    }else{
        std::copy(computed_values.begin(), computed_values.end(), new_pdf_values.data());
    }

    // clean memory and reorder the pdf values putting 0 in the invalid spots (for log case 0 should be -infty, hence using valid for accept/reject too
    if (!allValid){
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
}

double* TasmanianDREAM::getPDFHistory() const{
    double *hist = new double[pdf_history.size()];
    std::copy(pdf_history.begin(), pdf_history.end(), hist);
    return hist;
}

void TasmanianDREAM::getPDFHistory(std::vector<double> &history) const{
    history = pdf_history; // copy assignment
}

}

#endif
