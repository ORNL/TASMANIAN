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

#ifndef __TASGRID_WRAPPER_CPP
#define __TASGRID_WRAPPER_CPP

#include "tsgExoticQuadrature.hpp"
#include "tasgridWrapper.hpp"

TasgridWrapper::TasgridWrapper() : command(command_none), num_dimensions(0), num_outputs(-1), depth(-1), order(1),
    depth_type(type_none), rule(rule_none),
    conformal(conformal_none), alpha(0.0), beta(0.0), set_alpha(false), set_beta(false), tolerance(0.0), set_tolerance(false),
    ref_output(-1), min_growth(-1), tref(refine_fds), set_tref(false),
    printCout(false), useASCII(false), set_gpuid(-1), shift(0.0), set_shift(false)
{}
TasgridWrapper::~TasgridWrapper(){}

TypeCommand TasgridWrapper::hasCommand(std::string const &s){
    std::map<std::string, TypeCommand> commands = {
            {"-makeglobal",     command_makeglobal},     {"-mg", command_makeglobal},
            {"-makesequence",   command_makesequence},   {"-ms", command_makesequence},
            {"-makelocalpoly",  command_makelocalp},     {"-mp", command_makelocalp},
            {"-makewavelet",    command_makewavelet},    {"-mw", command_makewavelet},
            {"-makefourier",    command_makefourier},    {"-mf", command_makefourier},
            {"-makequadrature", command_makequadrature}, {"-mq", command_makequadrature},
            {"-makeexoquad",    command_makeexoquad},    {"-meq", command_makeexoquad},
            {"-makeupdate",      command_update},          {"-mu",   command_update},
            {"-setconformal",    command_setconformal},    {"-sc",   command_setconformal},
            {"-getquadrature",   command_getquadrature},   {"-gq",   command_getquadrature},
            {"-getinterweights", command_getinterweights}, {"-gi",   command_getinterweights},
            {"-getpoints",  command_getpoints},  {"-gp", command_getpoints},
            {"-getneeded",  command_getneeded},  {"-gn", command_getneeded},
            {"-loadvalues", command_loadvalues}, {"-l",  command_loadvalues},
            {"-evaluate",   command_evaluate},   {"-e",  command_evaluate},
            {"-integrate",  command_integrate},  {"-i",  command_integrate},
            {"-evalhierarchyd", command_evalhierarchical_dense},  {"-ehd", command_evalhierarchical_dense},
            {"-evalhierarchys", command_evalhierarchical_sparse}, {"-ehs", command_evalhierarchical_sparse},
            {"-gethsupport", command_gethsupport}, {"-ghsup", command_gethsupport},
            {"-getanisotropy", command_getanisocoeff}, {"-ga", command_getanisocoeff},
            {"-refinesurp",    command_refine_surp},   {"-rs", command_refine_surp},
            {"-refineaniso",   command_refine_aniso},  {"-ra", command_refine_aniso},
            {"-refine",        command_refine},        {"-r",  command_refine},
            {"-cancelrefine",  command_refine_clear},  {"-cr",   command_refine_clear},
            {"-mergerefine",   command_refine_merge},  {"-mr",   command_refine_merge},
            {"-summary", command_summary}, {"-s",   command_summary},
            {"-getcoefficients", command_getcoefficients}, {"-gc", command_getcoefficients},
            {"-setcoefficients", command_setcoefficients}, {"-sc", command_setcoefficients},
            {"-getpoly",          command_getpoly},
            {"-getpointsindexes", command_getpointsindex},
            {"-getneededindexes", command_getneededindex}
        };

    try{
        return commands.at(s);
    }catch(std::out_of_range &){
        return command_none;
    }
}

TypeConformalMap TasgridWrapper::getConfromalType(const char* name){
    if (std::string(name) == "asin"){
        return conformal_asin;
    }else{
        return conformal_none;
    }
}

bool TasgridWrapper::isCreateCommand(TypeCommand com){
    return ( (com == command_makeglobal) || (com == command_makesequence) || (com == command_makelocalp) ||
             (com == command_makewavelet) || (com == command_makefourier) || (com == command_makequadrature) ||
             (com == command_makeexoquad));
}

bool TasgridWrapper::checkSane() const{
    bool pass = true;
    if (command == command_none){
        cerr << "ERROR: checkSane(), no command specified" << endl;  return false;
    }else if (command == command_makeglobal){
        if (num_dimensions < 1){ cerr << "ERROR: must specify number of dimensions" << endl; pass = false; }
        if (num_outputs < 0){ cerr << "ERROR: must specify number of outputs (could use zero)" << endl; pass = false; }
        if (depth < 0){ cerr << "ERROR: must specify depth (e.g., level or polynomial degree)" << endl; pass = false; }
        if (depth_type == type_none){ cerr << "ERROR: must specify depth_type (e.g., select levels or polynomial basis)" << endl; pass = false; }
        if (rule == rule_none){ cerr << "ERROR: must specify rule to use (e.g., clenshaw-curtis)" << endl; pass = false; }
        if (!(OneDimensionalMeta::isGlobal(rule))){ cerr << "ERROR: cannot use global grids with rule: " << IO::getRuleString(rule) << endl; pass = false; }
        if ((rule == rule_gaussgegenbauer) || (rule == rule_gausslaguerre) || (rule == rule_gausshermite) || (rule == rule_gaussgegenbauerodd) || (rule == rule_gausshermiteodd) ){
            if (!set_alpha){ cerr << "ERROR: one dimensional rule " << IO::getRuleString(rule) << " requires alpha parameter" << endl; pass = false; }
        }else if (rule == rule_gaussjacobi){
            if (!set_alpha){ cerr << "ERROR: one dimensional rule " << IO::getRuleString(rule) << " requires alpha parameter" << endl; pass = false; }
            if (!set_beta ){ cerr << "ERROR: one dimensional rule " << IO::getRuleString(rule) << " requires beta parameter" << endl; pass = false; }
        }else{
            if (set_alpha){ cerr << "WARNING: alpha parameter set, but one dimensional rule " << IO::getRuleString(rule) << " doesn't depend on alpha" << endl; }
            if (set_beta ){ cerr << "WARNING: beta parameter set, but one dimensional rule " << IO::getRuleString(rule) << " doesn't depend on beta" << endl; }
        }
        if (customfilename.empty()){
            if (rule == rule_customtabulated){
                cerr << "ERROR: custom-tabulated rule specified, but no -customflile given" << endl; pass = false;
            }
        }else{
            if (rule != rule_customtabulated){
                cerr << "WARNING: -customflile given, but is only valid for the custom-tabulated rule and not rule " << IO::getRuleString(rule) << endl;
            }
        }
        if (gridfilename.empty() && outfilename.empty() && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -gridfile, -outfile or -print" << endl; pass = false;
        }
        if ((conformal == conformal_none) != conformalfilename.empty()){
            cerr << "WARNING: conformal transform requires both -conformaltype and -conformalfile, ignoring conformal mapping" << endl;
        }
        return pass;
    }else if (command == command_makesequence){
        if (num_dimensions < 1){ cerr << "ERROR: must specify number of dimensions" << endl; pass = false; }
        if (num_outputs < 0){ cerr << "ERROR: must specify number of outputs (could use zero)" << endl; pass = false; }
        if (depth < 0){ cerr << "ERROR: must specify depth (e.g., level or polynomial degree)" << endl; pass = false; }
        if (depth_type == type_none){ cerr << "ERROR: must specify depth_type (e.g., select levels or polynomial basis)" << endl; pass = false; }
        if (rule == rule_none){ cerr << "ERROR: must specify rule to use (e.g., rleja)" << endl; pass = false; }
        if (!(OneDimensionalMeta::isSequence(rule))){ cerr << "ERROR: rule is set to " << IO::getRuleString(rule) << " which is not a sequence rule (e.g., leja, rleja, min/max-lebesgue)" << endl; pass = false; }
        if (gridfilename.empty() && outfilename.empty() && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -gridfile, -outfile or -print" << endl; pass = false;
        }
        if ((conformal == conformal_none) != conformalfilename.empty()){
            cerr << "WARNING: conformal transform requires both -conformaltype and -conformalfile, ignoring conformal mapping" << endl;
        }
        return pass;
    }else if (command == command_makelocalp){
        if (num_dimensions < 1){ cerr << "ERROR: must specify number of dimensions" << endl; pass = false; }
        if (num_outputs < 0){ cerr << "ERROR: must specify number of outputs (could use zero)" << endl; pass = false; }
        if (depth < 0){ cerr << "ERROR: must specify depth (e.g., level or polynomial degree)" << endl; pass = false; }
        if (order < -1){ cerr << "ERROR: the maximum order cannot be less than -1"; pass = false; }
        if (rule == rule_none){ cerr << "ERROR: must specify rule to use (e.g., localp)" << endl; pass = false; }
        if (!(OneDimensionalMeta::isLocalPolynomial(rule))){ cerr << "ERROR: cannot use a local polynomial grid with rule: " << IO::getRuleString(rule) << endl; pass = false; }
        if (gridfilename.empty() && outfilename.empty() && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -gridfile, -outfile or -print" << endl; pass = false;
        }
        if ((conformal == conformal_none) != conformalfilename.empty()){
            cerr << "WARNING: conformal transform requires both -conformaltype and -conformalfile, ignoring conformal mapping" << endl;
        }
        return pass;
    }else if (command == command_makewavelet){
        if (num_dimensions < 1){ cerr << "ERROR: must specify number of dimensions" << endl; pass = false; }
        if (num_outputs < 0){ cerr << "ERROR: must specify number of outputs (could use zero)" << endl; pass = false; }
        if (depth < 0){ cerr << "ERROR: must specify depth (e.g., level or polynomial degree)" << endl; pass = false; }
        if ((order != 1) && (order != 3)){ cerr << "ERROR: the order must be either 1 or 3"; pass = false; }
        if (gridfilename.empty() && outfilename.empty() && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -gridfile, -outfile or -print" << endl; pass = false;
        }
        if ((conformal == conformal_none) != conformalfilename.empty()){
            cerr << "WARNING: conformal transform requires both -conformaltype and -conformalfile, ignoring conformal mapping" << endl;
        }
        return pass;
    }else if (command == command_makefourier){
        if (num_dimensions < 1){ cerr << "ERROR: must specify number of dimensions" << endl; pass = false; }
        if (num_outputs < 0){ cerr << "ERROR: must specify number of outputs (could use zero)" << endl; pass = false; }
        if (depth < 0){ cerr << "ERROR: must specify depth" << endl; pass = false; }
        if (depth_type == type_none){ cerr << "ERROR: must specify depth_type" << endl; pass = false; }
        if (gridfilename.empty() && outfilename.empty() && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -gridfile, -outfile or -print" << endl; pass = false;
        }
        if ((conformal == conformal_none) != conformalfilename.empty()){
            cerr << "WARNING: conformal transform requires both -conformaltype and -conformalfile, ignoring conformal mapping" << endl;
        }
        return pass;
    }else if (command == command_makequadrature){
        if (num_dimensions < 1){ cerr << "ERROR: must specify number of dimensions" << endl; pass = false; }
        if (num_outputs != -1){ cerr << "WARNING: ignoring the -outputs specified for the -makequadrature command" << endl; }
        if (depth < 0){ cerr << "ERROR: must specify depth (e.g., level or polynomial degree)" << endl; pass = false; }
        if (rule == rule_none){ cerr << "ERROR: must specify rule to use (e.g., clenshaw-curtis)" << endl; pass = false; }
        if (OneDimensionalMeta::isGlobal(rule)){ // global quadrature
            if (depth_type == type_none){ cerr << "ERROR: must specify depth_type (e.g., select levels or polynomial basis)" << endl; pass = false; }

            if ((rule == rule_gaussgegenbauer) || (rule == rule_gausslaguerre) || (rule == rule_gausshermite) || (rule == rule_gaussgegenbauerodd) || (rule == rule_gausshermiteodd)){
                if (!set_alpha){ cerr << "ERROR: one dimensional rule " << IO::getRuleString(rule) << " requires alpha parameter" << endl; pass = false; }
            }
            if (rule == rule_gaussjacobi){
                if (!set_alpha){ cerr << "ERROR: one dimensional rule " << IO::getRuleString(rule) << " requires alpha parameter" << endl; pass = false; }
                if (!set_beta ){ cerr << "ERROR: one dimensional rule " << IO::getRuleString(rule) << " requires beta parameter" << endl; pass = false; }
            }
        }else if (OneDimensionalMeta::isLocalPolynomial(rule)){
            if (order < -1){ cerr << "ERROR: the maximum order cannot be less than -1"; pass = false; }
        }else if (OneDimensionalMeta::isWavelet(rule)){
            if ((order != 1) && (order != 3)){ cerr << "ERROR: the order must be either 1 or 3"; pass = false; }
        }else if (OneDimensionalMeta::isFourier(rule)){
            if (depth_type == type_none){ cerr << "ERROR: must specify depth_type (e.g., select levels or polynomial basis)" << endl; pass = false; }
        }else{
            if (rule == rule_none){
                cerr << "ERROR: must specify rule to use (e.g., clenshaw-curtis or localp)" << endl; pass = false;
            }else{
                cerr << "ERROR: cannot make a quadrature with rule " << IO::getRuleString(rule) << endl; pass = false;
            }
        }
        if (outfilename.empty() && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -outfile or -print" << endl; pass = false;
        }
        if (!gridfilename.empty()){  cerr << "WARNING: quadrature does not output a -gridfile, if you need a gridfile use -makeglobal/-makelocalpoly commands followed by -getquadrature" << endl; }
        if ((conformal == conformal_none) != conformalfilename.empty()){
            cerr << "WARNING: conformal transform requires both -conformaltype and -conformalfile, ignoring conformal mapping" << endl;
        }
        return pass;
    }else if (command == command_makeexoquad){
        if (depth < 0){ cerr << "ERROR: must specify depth (e.g., level or polynomial degree)" << endl; pass = false; }
        if (!set_shift){ cerr << "ERROR: must specify shift parameter" << endl; pass = false; }
        if (weightfilename.empty()){ cerr << "ERROR: must specify filename of weight function surrogate/interpolant" << endl; pass = false; }
        if (description.empty()){ cerr << "ERROR: must specify description" << endl; pass = false; }
        if (outfilename.empty() && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -outfile or -print" << endl; pass = false;
        }
        return pass;

    }else if (command == command_update){
        if (gridfilename.empty()){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        if (depth < 0){ cerr << "ERROR: must specify depth (e.g., level or polynomial degree)" << endl; pass = false; }
        if (depth_type == type_none){ cerr << "ERROR: must specify depth_type (e.g., select levels or polynomial basis)" << endl; pass = false; }
    }else if (command == command_setconformal){
        if (gridfilename.empty()){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        if (conformal == conformal_none){  cerr << "ERROR: must specify valid -conformaltype" << endl; pass = false;  }
        if (conformalfilename.empty()){ cerr << "ERROR: must specify valid -conformalfile" << endl; pass = false; }
    }else if ((command == command_getquadrature)   || (command == command_getpoints) || (command == command_getneeded) ||
              (command == command_getcoefficients) || (command == command_gethsupport)){
        if (gridfilename.empty()){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        if (outfilename.empty() && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -outfile or -print" << endl; pass = false;
        }
        return pass;
    }else if ((command == command_loadvalues) || (command == command_setcoefficients)){
        if (gridfilename.empty()){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        if (valsfilename.empty()){ cerr << "ERROR: must specify valid -valsfile" << endl; pass = false; }
        return pass;
    }else if ((command == command_getinterweights) || (command == command_evaluate)
              || (command == command_evalhierarchical_dense) || (command == command_evalhierarchical_sparse)){
        if (gridfilename.empty()){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        if (xfilename.empty()){ cerr << "ERROR: must specify valid -pointsfile" << endl; pass = false; }
        if (outfilename.empty() && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -outfile or -print" << endl; pass = false;
        }
    }else if (command == command_integrate){
        if (gridfilename.empty()){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        if (outfilename.empty() && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -outfile or -print" << endl; pass = false;
        }
    }else if (command == command_getanisocoeff){
        if (depth_type == type_none){
            cerr << "ERROR: must specify type of coefficients with valid -type!" << endl;
            return false;
        }
        if (gridfilename.empty()){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        if (outfilename.empty() && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -outfile or -print" << endl; pass = false;
        }
    }else if (command == command_refine){
        if (gridfilename.empty()){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        // additional checks in refineGrid() since checks depend on the type of grid
    }else if (command == command_refine_aniso){
        if (gridfilename.empty()){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        // additional checks in refineGrid() since checks depend on the type of grid
    }else if (command == command_refine_surp){
        if (gridfilename.empty()){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        // additional checks in refineGrid() since checks depend on the type of grid
    }else if ((command == command_refine_clear) || (command == command_refine_merge)){
        if (gridfilename.empty()){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
    }else if (command == command_getrefcoeff){
        if (gridfilename.empty()){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
    }else if (command == command_getpoly){
        if (gridfilename.empty()){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        if (depth_type == type_none){ cerr << "ERROR: must specify depth_type (e.g., select levels or polynomial basis)" << endl; pass = false; }
        if ((depth_type != type_iptotal) && (depth_type != type_ipcurved) && (depth_type != type_iptensor) && (depth_type != type_iphyperbolic) &&
             (depth_type != type_qptotal) && (depth_type != type_qpcurved) && (depth_type != type_qptensor) && (depth_type != type_qphyperbolic)){
            cerr << "ERROR: the type here must start with either i or q indicating whether we seek the polynomils for integration or interpolation." << endl; pass = false;
        }
        if (outfilename.empty() && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -outfile or -print" << endl; pass = false;
        }
    }else if (command == command_summary){
        if (gridfilename.empty()){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
    }else if (command == command_getcoefficients){
        if (gridfilename.empty()){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
    }

    return pass;
}

void TasgridWrapper::createGlobalGird(){
    auto weights = readAnisotropicFile((OneDimensionalMeta::getControurType(depth_type) == type_curved) ? 2*num_dimensions : num_dimensions);
    auto llimits = readLevelLimits(num_dimensions);
    grid.makeGlobalGrid(num_dimensions, num_outputs, depth, depth_type, rule, weights, alpha, beta, customfilename.c_str(), llimits);
    if (!transformfilename.empty()){
        auto transforms = readTransform();
        grid.setDomainTransform(transforms.first, transforms.second);
    }
    if ((conformal != conformal_none) && (!conformalfilename.empty())){
        if (!setConformalTransformation()){
            cerr << "ERROR: could not set conformal transform" << endl;
        }
    }
}
void TasgridWrapper::createSequenceOrFourierGird(){
    auto weights = readAnisotropicFile((OneDimensionalMeta::getControurType(depth_type) == type_curved) ? 2*num_dimensions : num_dimensions);
    auto llimits = readLevelLimits(num_dimensions);
    if (command == command_makefourier || rule == rule_fourier){    // rule condition is for quadrature
        grid.makeFourierGrid(num_dimensions, num_outputs, depth, depth_type, weights, llimits);
    }else{
        grid.makeSequenceGrid(num_dimensions, num_outputs, depth, depth_type, rule, weights, llimits);
    }
    if (!transformfilename.empty()){
        auto transforms = readTransform();
        grid.setDomainTransform(transforms.first, transforms.second);
    }
    if ((conformal != conformal_none) && (!conformalfilename.empty())){
        if (!setConformalTransformation()){
            cerr << "ERROR: could not set conformal transform" << endl;
        }
    }
}
void TasgridWrapper::createLocalPolynomialGird(){
    auto llimits = readLevelLimits(num_dimensions);
    grid.makeLocalPolynomialGrid(num_dimensions, num_outputs, depth, order, rule, llimits);
    if (!transformfilename.empty()){
        auto transforms = readTransform();
        grid.setDomainTransform(transforms.first, transforms.second);
    }
    if ((conformal != conformal_none) && (!conformalfilename.empty())){
        if (!setConformalTransformation()){
            cerr << "ERROR: could not set conformal transform" << endl;
        }
    }
}
void TasgridWrapper::createWaveletGird(){
    auto llimits = readLevelLimits(num_dimensions);
    grid.makeWaveletGrid(num_dimensions, num_outputs, depth, order, llimits);
    if (!transformfilename.empty()){
        auto transforms = readTransform();
        grid.setDomainTransform(transforms.first, transforms.second);
    }
    if ((conformal != conformal_none) && (!conformalfilename.empty())){
        if (!setConformalTransformation()){
            cerr << "ERROR: could not set conformal transform" << endl;
        }
    }
}
void TasgridWrapper::createQuadrature(){
    if (num_outputs != 0) num_outputs = 0;
    if (OneDimensionalMeta::isGlobal(rule)){
        createGlobalGird();
    }else if (OneDimensionalMeta::isLocalPolynomial(rule)){
        createLocalPolynomialGird();
    }else if (OneDimensionalMeta::isWavelet(rule)){
        createWaveletGird();
    }else if (OneDimensionalMeta::isFourier(rule)){
        createSequenceOrFourierGird();
    }else{
        cerr << "ERROR: createQuadrature" << endl;
    }
}
void TasgridWrapper::createExoticQuadrature(){
    TasGrid::TasmanianSparseGrid weight_surrogate;
    weight_surrogate.read(weightfilename.c_str());
    ct = TasGrid::getExoticQuadrature(depth, shift, weight_surrogate, description.c_str(), is_symmetric_weight_function);
}
bool TasgridWrapper::updateGrid(){
    if (!(grid.isGlobal() || grid.isSequence() || grid.isFourier())){
        cerr << "ERROR: -makeupdate can be called only for Global and Sequence grids" << endl;
        return false;
    }
    auto weights = readAnisotropicFile((OneDimensionalMeta::getControurType(depth_type) == type_curved) ? 2*num_dimensions : num_dimensions);
    if (grid.isGlobal()){
        grid.updateGlobalGrid(depth, depth_type, weights);
    }else if (grid.isSequence()){
        grid.updateSequenceGrid(depth, depth_type, weights);
    }else{
        grid.updateFourierGrid(depth, depth_type, weights);
    }
    return true;
}
void TasgridWrapper::writeGrid() const{ grid.write(gridfilename.c_str(), (useASCII) ? mode_ascii : mode_binary); }
bool TasgridWrapper::readGrid(){
    try{
        grid.read(gridfilename.c_str());
        return true;
    }catch(std::runtime_error &e){
        cerr << e.what() << endl;
        return false;
    }
}
void TasgridWrapper::outputPoints(bool useNeeded) const{
    int num_p, num_d = grid.getNumDimensions();
    std::vector<double> points;
    if (outfilename.empty() && (printCout == false)) return;
    if (useNeeded){
        num_p = grid.getNumNeeded();
        points = grid.getNeededPoints();
    }else{
        num_p = grid.getNumPoints();
        points = grid.getPoints();
    }
    writeMatrix(outfilename, num_p, num_d, points.data());
    printMatrix(num_p, num_d, points.data());
}
void TasgridWrapper::outputQuadrature() const{
    if (outfilename.empty() && (printCout == false)) return;
    int num_p = grid.getNumPoints();
    int num_d = grid.getNumDimensions();
    int offset = num_d + 1;
    auto points  = grid.getPoints();
    auto weights = grid.getQuadratureWeights();
    Data2D<double> combined(num_d + 1, num_p);
    auto ip = points.begin();
    auto iw = weights.begin();
    for(int i=0; i<num_p; i++){
        double *c = combined.getStrip(i);
        c[0] = *iw++;
        std::copy_n(ip, size_t(num_d), &(c[1]));
        std::advance(ip, num_d);
    }
    writeMatrix(outfilename, num_p, offset, combined.getStrip(0));
    printMatrix(num_p, offset, combined.getStrip(0));
}
void TasgridWrapper::outputExoticQuadrature() const{
    if (!outfilename.empty()){
        std::ofstream ofs(outfilename, std::ios::out | std::ios::trunc);
        ct.write<mode_ascii>(ofs);
    }
    if (printCout){
        ct.write<mode_ascii>(cout);
    }
    return;
}
void TasgridWrapper::outputHierarchicalCoefficients() const{
    const double *coeff = grid.getHierarchicalCoefficients();
    int num_pnts = grid.getNumPoints();
    int num_outs = grid.getNumOutputs();
    if (grid.isFourier()){
        // use interwoven format for complex coefficients
        Data2D<double> coeff_fourier(2 * num_outs, num_pnts);
        Utils::Wrapper2D<const double> real(num_outs, coeff);
        Utils::Wrapper2D<const double> imag(num_outs, real.getStrip(num_pnts));
        for(int p=0; p<num_pnts; p++){
            double *c = coeff_fourier.getStrip(p);
            double const *r = real.getStrip(p);
            double const *i = imag.getStrip(p);

            for(size_t j=0; j<size_t(num_outs); j++){
                c[2*j    ] = r[j];
                c[2*j + 1] = i[j];
            }
        }
        writeMatrix(outfilename, num_pnts, 2 * num_outs, coeff_fourier.getStrip(0));
        printMatrix(num_pnts, num_outs, coeff_fourier.getStrip(0), true);
    }else{
        writeMatrix(outfilename, num_pnts, num_outs, coeff);
        printMatrix(num_pnts, num_outs, coeff);
    }
}
void TasgridWrapper::outputHierachicalSupport() const{
    std::vector<double> supp = grid.getHierarchicalSupport();
    writeMatrix(outfilename, grid.getNumPoints(), grid.getNumDimensions(), supp.data());
    printMatrix(grid.getNumPoints(), grid.getNumDimensions(), supp.data());
}
bool TasgridWrapper::setConformalTransformation(){
    if (conformal == conformal_asin){
        auto mat = readMatrix(conformalfilename);
        if (mat.getNumStrips() != 1){
            cerr << "ERROR: the conformal file for asin should contain only one row" << endl;
            return false;
        }
        if (mat.getStride() != (size_t) grid.getNumDimensions()){
            cerr << "ERROR: the conformal file for asin should contain " << grid.getNumDimensions() << " columns, instead it has " << mat.getStride() << endl;
            return false;
        }
        std::vector<int> coeff(mat.getTotalEntries());
        std::transform(mat.begin(), mat.end(), coeff.begin(), [](double x)->int{ return static_cast<int>(x); });

        grid.setConformalTransformASIN(coeff);
        return true;
    }else{
        cerr << "ERROR: conformal type not implemented" << endl;
        return false;
    }
}
bool TasgridWrapper::loadValues(){
    auto vals = readMatrix(valsfilename);
    if (grid.getNumNeeded() == 0){
        if (vals.getNumStrips() != grid.getNumLoaded()){
            cerr << "ERROR: grid has " << grid.getNumLoaded() << " new values, but " << valsfilename << " specifies " << vals.getNumStrips() << endl;
            return false;
        }
    }else{
        if (vals.getNumStrips() != grid.getNumNeeded()){
            cerr << "ERROR: grid is awaiting " << grid.getNumNeeded() << " new values, but " << valsfilename << " specifies " << vals.getNumStrips() << endl;
            return false;
        }
    }
    if (vals.getStride() != (size_t) grid.getNumOutputs()){
        cerr << "ERROR: grid is set for " << grid.getNumOutputs() << " outputs, but " << valsfilename << " specifies " << vals.getStride() << endl;
        return false;
    }
    grid.loadNeededPoints(vals.data());
    return true;
}
bool TasgridWrapper::getInterWeights(){
    auto x = readMatrix(xfilename);
    if (x.getStride() != (size_t) grid.getNumDimensions()){
        cerr << "ERROR: grid is set for " << grid.getNumDimensions() << " dimensions, but " << xfilename << " specifies " << x.getStride() << endl;
        return false;
    }
    if (x.empty()){
        cerr << "ERROR: no points specified in " << xfilename << endl;
        return false;
    }
    int numx = x.getNumStrips();
    int num_p = grid.getNumPoints();
    Data2D<double> result(grid.getNumPoints(), numx);

    #pragma omp parallel for
    for(int i=0; i<numx; i++) // in windows OpenMP loop counters must be signed ??
        grid.getInterpolationWeights(x.getStrip(i), result.getStrip(i));

    writeMatrix(outfilename, numx, (int) num_p, result.getStrip(0));
    printMatrix(numx, (int) num_p, result.getStrip(0));

    return true;
}
bool TasgridWrapper::getEvaluate(){
    if (grid.getNumLoaded() == 0){
        cerr << "ERROR: no values loaded in the grid, cannot evaluate!" << endl;
        return false;
    }
    if (grid.getNumOutputs() == 0){
        cerr << "ERROR: no outputs set for the grid, nothing to evaluate!" << endl;
        return false;
    }
    auto x = readMatrix(xfilename);
    if (x.getStride() != (size_t) grid.getNumDimensions()){
        cerr << "ERROR: grid is set for " << grid.getNumDimensions() << " dimensions, but " << xfilename << " specifies " << x.getStride() << endl;
        return false;
    }
    if (x.empty()){
        cerr << "ERROR: no points specified in " << xfilename << endl;
        return false;
    }
    int num_out = grid.getNumOutputs();

    if (set_gpuid > -1){
        if (set_gpuid < grid.getNumGPUs()){
            grid.enableAcceleration(accel_gpu_cuda);
            grid.setGPUID(set_gpuid);
        }else{
            cerr << "WARNING: invalud GPU specified " << set_gpuid << endl;
        }
    }

    int num_points = x.getNumStrips();
    std::vector<double> result;
    grid.evaluateBatch(x.release(), result);

    writeMatrix(outfilename, num_points, (int) num_out, result.data());
    printMatrix(num_points, (int) num_out, result.data());

    return true;
}
bool TasgridWrapper::getIntegrate(){
    if (grid.getNumLoaded() == 0){
        cerr << "ERROR: no values loaded in the grid, cannot evaluate!" << endl;
        return false;
    }
    if (grid.getNumOutputs() == 0){
        cerr << "ERROR: no outputs set for the grid, nothing to evaluate!" << endl;
        return false;
    }
    int num_out = grid.getNumOutputs();
    std::vector<double> q;
    grid.integrate(q);

    writeMatrix(outfilename, 1, num_out, q.data());
    printMatrix(1, num_out, q.data());

    return true;
}
bool TasgridWrapper::getAnisoCoeff(){
    if (grid.getNumOutputs() == 0){
        cerr << "ERROR: cannot estimate coefficients with no outputs!" << endl;
        return false;
    }
    if (grid.getNumLoaded() == 0){
        cerr << "ERROR: cannot estimate coefficients for a grid with no loaded values!" << endl;
        return false;
    }
    std::vector<int> ab;
    if (grid.isSequence()){
        ab = grid.estimateAnisotropicCoefficients(depth_type, ref_output);
    }else{
        if (ref_output >= grid.getNumOutputs()){
            cerr << "ERROR: -ref_output " << ref_output << " is specified, however, the grid has only " << grid.getNumOutputs() << " outputs!" << endl;
            cerr << " HINT: the outputs are indexed starting at zero!" << endl;
            return false;
        }
        if ((ref_output == -1) && (grid.getNumOutputs() > 1)){
            cerr << "ERROR: must specify a refinement output with -ref_output option!" << endl;
            return false;
        }else if (ref_output == -1) ref_output = 0;
        ab = grid.estimateAnisotropicCoefficients(depth_type, ref_output);
    }
    size_t num_coeff = (size_t) grid.getNumDimensions();
    if (OneDimensionalMeta::getControurType(depth_type) == type_curved)
        num_coeff *= 2;

    std::vector<double> coeff(num_coeff);
    std::transform(ab.begin(), ab.end(), coeff.begin(), [](int i)->double{ return double(i); });

    writeMatrix(outfilename, (int) num_coeff, 1, coeff.data());
    printMatrix((int) num_coeff, 1, coeff.data());

    return true;
}

bool TasgridWrapper::refineGrid(){
    // put the sanity check code here, since many of the parameters of the refinement depend on the type of grid being used
    if (grid.getNumOutputs() == 0){
        cerr << "ERROR: cannot refine a grid with no outputs!" << endl;
        return false;
    }
    if (grid.getNumLoaded() == 0){
        cerr << "ERROR: cannot refine a grid with no loaded values!" << endl;
        return false;
    }
    auto llimits = readLevelLimits(grid.getNumDimensions());
    TypeCommand effective_command = command;
    if (command == command_refine){
        if (grid.isGlobal() || grid.isSequence()){
            effective_command = command_refine_aniso;
        }else{
            effective_command = command_refine_surp;
        }
    }
    if (effective_command == command_refine_aniso){
        if ((!grid.isGlobal()) && (!grid.isSequence()) && (!grid.isFourier())){
            cerr << "ERROR: anisotropic refinement can be used only for global/sequence/Fourier grids!" << endl;
            return false;
        }
        if (set_tolerance){
            cerr << "WARNING: anisotropic refinement ignores the -tolerance option!" << endl;
        }
        if (depth_type == type_none){
            cerr << "ERROR: anisotropic refinement requires -type!" << endl;
            return false;
        }
        if (set_tref){
            cerr << "WARNING: anisotropic refinement ignores the -reftype option!" << endl;
        }
        if (min_growth < 1) min_growth = 1;
        if (grid.isSequence()){
            grid.setAnisotropicRefinement(depth_type, min_growth, ref_output, llimits);
        }else{
            if (ref_output >= grid.getNumOutputs()){
                cerr << "ERROR: -ref_output " << ref_output << " is specified, however, the grid has only " << grid.getNumOutputs() << " outputs!" << endl;
                cerr << " HINT: the outputs are indexed starting at zero!" << endl;
                return false;
            }
            if ((ref_output == -1) && (grid.getNumOutputs() > 1)){
                cerr << "ERROR: must specify a refinement output with -ref_output option!" << endl;
                return false;
            }else if (ref_output == -1) ref_output = 0;
            grid.setAnisotropicRefinement(depth_type, min_growth, ref_output, llimits);
        }
    }else{ // using surplus refinement
        if (!set_tolerance){
            cerr << "ERROR: must specify -tolerance for surplus refinement!" << endl;
            return false;
        }
        if ((grid.isLocalPolynomial()) || (grid.isWavelet())){
            if (!set_tref){
                cerr << "ERROR: must specify -reftype option!" << endl;
                return false;
            }
            Data2D<double> scale;
            if (!valsfilename.empty()){
                scale = readMatrix(valsfilename);
                if (scale.getNumStrips() != grid.getNumPoints()){
                    cerr << "ERROR: the number of weights must match the number of points." << endl;
                    return false;
                }
                if ((ref_output == -1) && (scale.getStride() != (size_t) grid.getNumOutputs())){
                    cerr << "ERROR: the number of weights must match the number of outputs." << endl;
                    return false;
                }
                if ((ref_output > -1) && (scale.getStride() != 1)){
                    cerr << "ERROR: there must be one weight per output." << endl;
                    return false;
                }
            }
            grid.setSurplusRefinement(tolerance, tref, ref_output, llimits, scale.release());
        }else if (grid.isSequence()){
            grid.setSurplusRefinement(tolerance, ref_output, llimits);
        }else{
            if (ref_output >= grid.getNumOutputs()){
                cerr << "ERROR: -ref_output " << ref_output << " is specified, however, the grid has only " << grid.getNumOutputs() << " outputs!" << endl;
                cerr << " HINT: the outputs are indexed starting at zero!" << endl;
                return false;
            }
            if ((ref_output == -1) && (grid.getNumOutputs() > 1)){
                cerr << "ERROR: must specify a refinement output with -ref_output option!" << endl;
                return false;
            }else if (ref_output == -1) ref_output = 0;
            grid.setSurplusRefinement(tolerance, ref_output, llimits);
        }
    }
    writeGrid();
    return true;
}
bool TasgridWrapper::cancelRefine(){
    grid.clearRefinement();
    return true;
}
bool TasgridWrapper::mergeRefine(){
    grid.mergeRefinement();
    return true;
}
bool TasgridWrapper::getPoly(){
    if ((grid.isGlobal()) || (grid.isSequence())){
        int num_d = grid.getNumDimensions();
        bool integrate = ((depth_type == type_iptotal) || (depth_type == type_ipcurved) || (depth_type == type_iptensor) || (depth_type == type_iphyperbolic));
        std::vector<int> poly = grid.getGlobalPolynomialSpace(integrate);
        std::vector<double> double_poly(poly.size());
        std::transform(poly.begin(), poly.end(), double_poly.begin(), [](int x)->double{ return static_cast<double>(x); });
        writeMatrix(outfilename, (int) double_poly.size() / num_d, num_d, double_poly.data());
        printMatrix((int) double_poly.size() / num_d, num_d, double_poly.data());
    }else{
        cerr << "ERROR: cannot call -getpoly for a grid that is neither Global nor Sequence" << endl;
        return false;
    }
    return true;
}
bool TasgridWrapper::getSummary(){
    grid.printStats();
    return true;
}

bool TasgridWrapper::getEvalHierarchyDense(){
    auto x = readMatrix(xfilename);
    if (x.getStride() != (size_t) grid.getNumDimensions()){
        cerr << "ERROR: grid is set for " << grid.getNumDimensions() << " dimensions, but " << xfilename << " specifies " << x.getStride() << endl;
        return false;
    }
    if (x.empty()){
        cerr << "ERROR: no points specified in " << xfilename << endl;
        return false;
    }
    int num_p = grid.getNumPoints();
    int num_x = x.getNumStrips();
    auto result = grid.evaluateHierarchicalFunctions(x.release());

    writeMatrix(outfilename, num_x, ((grid.isFourier()) ? 2 * num_p : num_p), result.data());
    printMatrix(num_x, (int) num_p, result.data(), grid.isFourier());

    return true;
}
template<bool iomode>
void writeSparseMatrix(int cols, std::vector<int> const &pntr, std::vector<int> const &indx, std::vector<double> const &vals, std::ostream &os){
    int nnz = static_cast<int>(indx.size());
    int rows = static_cast<int>(pntr.size() - 1);
    IO::writeNumbers<iomode, IO::pad_line>(os, rows, cols, nnz);
    IO::writeVector<iomode, IO::pad_line>(pntr, os);
    IO::writeVector<iomode, IO::pad_line>(indx, os);
    IO::writeVector<iomode, IO::pad_line>(vals, os);
}
bool TasgridWrapper::getEvalHierarchySparse(){
    auto x = readMatrix(xfilename);
    if (x.getStride() != (size_t) grid.getNumDimensions()){
        cerr << "ERROR: grid is set for " << grid.getNumDimensions() << " dimensions, but " << xfilename << " specifies " << x.getStride() << endl;
        return false;
    }
    if (x.empty()){
        cerr << "ERROR: no points specified in " << xfilename << endl;
        return false;
    }
    std::vector<int> pntr , indx;
    std::vector<double> vals;
    grid.evaluateSparseHierarchicalFunctions(x.release(), pntr, indx, vals);
    int num_p = grid.getNumPoints();
    if (!outfilename.empty()){
        if (useASCII){
            std::ofstream ofs(outfilename);
            ofs << std::scientific; ofs.precision(17);
            writeSparseMatrix<mode_ascii>(num_p, pntr, indx, vals, ofs);
            ofs.close();
        }else{
            std::ofstream ofs(outfilename, std::ios::out | std::ios::binary);
            char charTSG[3] = {'T', 'S', 'G'};
            ofs.write(charTSG, 3 * sizeof(char));
            writeSparseMatrix<mode_binary>(num_p, pntr, indx, vals, ofs);
            ofs.close();
        }
    }
    if (printCout){
        cout << std::scientific; cout.precision(17);
        writeSparseMatrix<mode_ascii>(num_p, pntr, indx, vals, cout);
    }
    return true;
}
bool TasgridWrapper::setHierarchy(){
    auto vals = readMatrix(valsfilename);
    if (vals.getNumStrips() != grid.getNumPoints()){
        cerr << "ERROR: grid is awaiting " << grid.getNumPoints() << " hierarchical surpluses, but " << valsfilename
             << " specifies " << vals.getNumStrips() << endl;
        return false;
    }
    if (!(grid.isFourier()) && (vals.getStride() != (size_t) grid.getNumOutputs())){
        cerr << "ERROR: grid is set for " << grid.getNumOutputs() << " outputs, but " << valsfilename << " specifies " << vals.getStride() << endl;
        return false;
    }else if (grid.isFourier() && (vals.getStride() != (size_t) (2 * grid.getNumOutputs()))){
        cerr << "ERROR: fourier grid is set for " << grid.getNumOutputs() << " outputs, but " << valsfilename << " specifies " << (vals.getStride() / 2) << endl;
        return false;
    }
    grid.setHierarchicalCoefficients(vals.release());
    return true;
}

bool TasgridWrapper::getPointsIndexes(){
    const int *p = grid.getPointsIndexes();
    int num_p = grid.getNumPoints();
    int num_d = grid.getNumDimensions();
    Data2D<double> pv(num_d, num_p);
    std::transform(p, p + Utils::size_mult(num_p, num_d), pv.getStrip(0), [](int i)->double{ return double(i); });

    writeMatrix(outfilename, num_p, num_d, pv.getStrip(0));
    printMatrix(num_p, num_d, pv.getStrip(0));

    return true;
}
bool TasgridWrapper::getNeededIndexes(){
    const int *p = grid.getNeededIndexes();
    int num_p = grid.getNumNeeded();
    int num_d = grid.getNumDimensions();
    Data2D<double> pv(num_d, num_p);
    std::transform(p, p + Utils::size_mult(num_p, num_d), pv.getStrip(0), [](int i)->double{ return double(i); });

    writeMatrix(outfilename, num_p, num_d, pv.getStrip(0));
    printMatrix(num_p, num_d, pv.getStrip(0));

    return true;
}

std::vector<int> TasgridWrapper::readAnisotropicFile(int num_weights) const{
    if (anisofilename.empty()) return std::vector<int>();
    auto mat = readMatrix(anisofilename);
    if (mat.getNumStrips() != 1)
        throw std::runtime_error("ERROR: anisotropy file must contain only one row");
    if (mat.getStride() != ((size_t) num_weights)){
        cerr << "ERROR: anisotropy file has wrong number of entries, " << num_weights << " expected " << mat.getStride() << " found." << endl;
        throw std::runtime_error("ERROR: anisotropy file has wrong number of entries");
    }

    std::vector<int> weights(mat.getTotalEntries());
    std::transform(mat.begin(), mat.end(), weights.begin(), [](double x)->int{ return static_cast<int>(x); });
    return weights;
}
std::pair<std::vector<double>, std::vector<double>> TasgridWrapper::readTransform() const{
    auto mat = readMatrix(transformfilename);
    if (mat.getStride() != 2){
        cerr << "ERROR: file " << transformfilename << " must have exactly two columns." << endl;
        throw std::runtime_error("ERROR: incorrect format for the transform file");
    }
    if (mat.getNumStrips() != num_dimensions){
        cerr << "ERROR: file " << transformfilename << " has " << mat.getNumStrips() << " rows, instead of the number of dimensions " << num_dimensions << endl;
        throw std::runtime_error("ERROR: incorrect format for the transform file");
    }
    std::vector<double> transa((size_t) num_dimensions);
    std::vector<double> transb((size_t) num_dimensions);
    for(int i=0; i<num_dimensions; i++){
        transa[i] = mat.getStrip(i)[0];
        transb[i] = mat.getStrip(i)[1];
    }
    return std::make_pair(transa, transb);
}
std::vector<int> TasgridWrapper::readLevelLimits(int num_weights) const{
    if (levellimitfilename.empty()) return std::vector<int>();
    auto mat = readMatrix(levellimitfilename);
    if (mat.getNumStrips() != 1)
        throw std::runtime_error("ERROR: level limits file must contain only one row");

    if (mat.getStride() != ((size_t)  num_weights)){
        cerr << "ERROR: level limits file has wrong number of entries, " << num_weights << " expected " << mat.getStride() << " found." << endl;
        throw std::runtime_error("ERROR: level limits file has incorrect format.");
    }

    std::vector<int> llimits(mat.getTotalEntries());
    std::transform(mat.begin(), mat.end(), llimits.begin(), [](double x)->int{ return static_cast<int>(x); });
    return llimits;
}

template<typename iomode>
Data2D<double> readMatrixFromOpen(std::istream &is){
    int rows = IO::readNumber<iomode, int>(is);
    int cols = IO::readNumber<iomode, int>(is);
    return Data2D<double>(cols, rows, IO::readVector<iomode, double>(is, Utils::size_mult(cols, rows)));
}

Data2D<double> TasgridWrapper::readMatrix(std::string const &filename){
    Data2D<double> matrix;
    if (filename.empty()) return matrix;
    std::ifstream ifs;
    ifs.open(filename, std::ios::in | std::ios::binary);
    if (!(ifs.good())){
        cerr << "ERROR: could not open file " << filename << endl;
        ifs.close();
        return matrix;
    }
    char tsg[3] = {'A', 'A', 'A'};
    ifs.read(tsg, 3*sizeof(char));
    if ((tsg[0] == 'T') && (tsg[1] == 'S') && (tsg[2] == 'G')){
        matrix = readMatrixFromOpen<IO::mode_binary_type>(ifs);
    }else{ // not a binary file
        ifs.close();
        ifs.open(filename);
        matrix = readMatrixFromOpen<IO::mode_ascii_type>(ifs);
    }
    if (matrix.empty())
        cerr << "WARNING: empty file " << filename << endl;
    ifs.close();
    return matrix;
}
void TasgridWrapper::writeMatrix(std::string const &filename, int rows, int cols, const double mat[]) const{
    if (filename.empty()) return;
    size_t cols_t = (size_t) cols;
    std::ofstream ofs;
    if (useASCII){
        Utils::Wrapper2D<const double> matrix(cols, mat);
        ofs.open(filename);
        ofs << rows << " " << cols << endl;
        ofs.precision(17);
        ofs << std::scientific;
        for(int i=0; i<rows; i++){
            double const * r = matrix.getStrip(i);
            ofs << setw(25) << r[0];
            for(size_t j=1; j<cols_t; j++){
                ofs << " " << setw(25) << r[j];
            }
            ofs << "\n";
        }
    }else{
        ofs.open(filename, std::ios::out | std::ios::binary);
        char tsg[3] = {'T', 'S', 'G'};
        ofs.write(tsg, 3*sizeof(char));
        ofs.write((char*) &rows, sizeof(int));
        ofs.write((char*) &cols, sizeof(int));
        ofs.write((char*) mat, Utils::size_mult(rows, cols) * sizeof(double));
    }
    ofs.close();
}
void TasgridWrapper::printMatrix(int rows, int cols, const double mat[], bool isComplex) const{
    if (!printCout) return;
    cout << rows << " " << cols << endl;
    cout.precision(17);
    cout << std::scientific;
    size_t cols_t = (size_t) cols;
    Utils::Wrapper2D<const double> matrix(cols, mat);
    for(int i=0; i<rows; i++){
        double const * r = matrix.getStrip(i);
        if (isComplex){
            cout << setw(50) << std::complex<double>(r[0], r[1]);
            for(size_t j=1; j<cols_t; j++)
                cout << setw(50) << std::complex<double>(r[2*j], r[2*j + 1]);
        }else{
            cout << setw(25) << r[0];
            for(size_t j=1; j<cols_t; j++)
                cout << " " << setw(25) << r[j];
        }
        cout << '\n';
    }
    cout << endl;
}

bool TasgridWrapper::executeCommand(){
    if (!checkSane()) return false;
    if (command == command_makeglobal){
        createGlobalGird();
    }else if (command == command_makesequence){
        createSequenceOrFourierGird();
    }else if (command == command_makelocalp){
        createLocalPolynomialGird();
    }else if (command == command_makewavelet){
        createWaveletGird();
    }else if (command == command_makefourier){
        createSequenceOrFourierGird();
    }else if (command == command_makequadrature){
        createQuadrature();
        outputQuadrature();
    }else if (command == command_makeexoquad){
        createExoticQuadrature();
        outputExoticQuadrature();
    }
    if (isCreateCommand(command)){
        if (command != command_makequadrature && command != command_makeexoquad){
            outputPoints(false);
            if (!gridfilename.empty()) writeGrid();
        }
    }else{
        if (!readGrid()){
            cerr << "ERROR: could not read the grid in file: " << gridfilename << ", please specify -gridfile" << endl;
            return false;
        }
    }

    if (command == command_update){
        if (!updateGrid()){
            cerr << "ERROR: could not update the grid" << endl;
            return false;
        }
    }else if (command == command_setconformal){
        if (setConformalTransformation()){
            writeGrid();
        }else{
            cerr << "ERROR: could not set the conformal grid" << endl;
            return false;
        }
    }else if (command == command_getcoefficients){
        outputHierarchicalCoefficients();
    }else if (command == command_getquadrature){
        outputQuadrature();
    }else if (command == command_getpoints){
        outputPoints(false);
    }else if (command == command_getneeded){
        outputPoints(true);
    }else if (command == command_gethsupport){
        outputHierachicalSupport();
    }else if (command == command_loadvalues){
        if (loadValues()){
            writeGrid();
        }else{
            cerr << "ERROR: values could not be loaded!" << endl;
        }
    }else if (command == command_getinterweights){
        if (!getInterWeights()){
            cerr << "ERROR: could not generate interpolation weights" << endl;
            return false;
        }
    }else if (command == command_evaluate){
        if (!getEvaluate()){
            cerr << "ERROR: could not evaluate the grid" << endl;
            return false;
        }
    }else if (command == command_integrate){
        if (!getIntegrate()){
            cerr << "ERROR: could not integrate the grid" << endl;
            return false;
        }
    }else if (command == command_getanisocoeff){
        getAnisoCoeff();
    }else if ((command == command_refine)||(command == command_refine_aniso)||(command == command_refine_surp)){
        refineGrid();
        outputPoints(true);
    }else if (command == command_refine_clear){
        if (cancelRefine()){
            writeGrid();
        }else{
            cerr << "ERROR: could not clear the refinement" << endl;
        }
    }else if (command == command_refine_merge){
        if (mergeRefine()){
            writeGrid();
        }else{
            cerr << "ERROR: could not merge the refinement" << endl;
        }
    }else if (command == command_getpoly){
        if (!getPoly()){
            cerr << "ERROR: could not get polynomial basis" << endl;
            return false;
        }
    }else if (command == command_summary){
        getSummary();
    }else if (command == command_evalhierarchical_dense){
        if (!getEvalHierarchyDense()){
            cerr << "ERROR: could not evaluate the (dense) hierarchical basis functions" << endl;
            return false;
        }
    }else if (command == command_evalhierarchical_sparse){
        if (!getEvalHierarchySparse()){
            cerr << "ERROR: could not evaluate the (sparse) hierarchical basis functions" << endl;
            return false;
        }
    }else if (command == command_setcoefficients){
        if (setHierarchy()){
            writeGrid();
        }else{
            cerr << "ERROR: could not set the hierarchical coefficients" << endl;
            return false;
        }
    }else if (command == command_getpointsindex){
        if (!getPointsIndexes()){
            cerr << "ERROR: could not get the indexes" << endl;
        }
    }else if (command == command_getneededindex){
        if (!getNeededIndexes()){
            cerr << "ERROR: could not get the indexes" << endl;
        }
    }

    return true;
}

#endif
