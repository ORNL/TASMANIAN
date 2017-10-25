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

#include "tasgridWrapper.hpp"

TasgridWrapper::TasgridWrapper() : grid(0), command(command_none), num_dimensions(0), num_outputs(-1), depth(-1), order(1),
    depth_type(type_none), rule(rule_none),
    conformal(conformal_none), alpha(0.0), beta(0.0), set_alpha(false), set_beta(false), tolerance(0.0), set_tolerance(false),
    ref_output(-1), min_growth(-1), tref(refine_fds), set_tref(false),
    gridfilename(0), outfilename(0), valsfilename(0), xfilename(0), anisofilename(0), transformfilename(0), conformalfilename(0), customfilename(0),
    printCout(false), useASCII(false)
{}
TasgridWrapper::~TasgridWrapper(){
    if (grid != 0) delete grid;
}

TypeConformalMap TasgridWrapper::getConfromalType(const char* name){
    if (strcmp(name, "asin") == 0){
        return conformal_asin;
    }else{
        return conformal_none;
    }
}

bool TasgridWrapper::isCreateCommand(TypeCommand com){
    return ( (com == command_makeglobal) || (com == command_makesequence) || (com == command_makelocalp) || (com == command_makewavelet) || (com == command_makequadrature) );
}

void TasgridWrapper::setCommand(TypeCommand com){ command = com; }
TypeCommand TasgridWrapper::getCommand() const{ return command; }
void TasgridWrapper::setNumDimensions(int dim){ num_dimensions = dim; }
void TasgridWrapper::setNumOutputs(int out){ num_outputs = out; }
void TasgridWrapper::setNumDepth(int d){ depth = d; }
void TasgridWrapper::setOrder(int o){ order = o; }
void TasgridWrapper::setDepthType(TypeDepth dt){ depth_type = dt; }
void TasgridWrapper::setConformalType(TypeConformalMap con){ conformal = con; }
void TasgridWrapper::setRule(TypeOneDRule r){ rule = r;}
TypeOneDRule TasgridWrapper::getRule() const{ return rule; }
void TasgridWrapper::setAlpha(double a){ alpha = a;  set_alpha = true; }
void TasgridWrapper::setBeta(double b){ beta = b; set_beta = true; }
void TasgridWrapper::setTolerance(double tol){ tolerance = tol;  set_tolerance = true; }
void TasgridWrapper::setRefOutput(int out){ ref_output = out; }
void TasgridWrapper::setMinGrowth(int mg){ min_growth = mg; }
void TasgridWrapper::setTypeRefinement(TypeRefinement rt){ tref = rt; set_tref = true; }
void TasgridWrapper::setGridFilename(const char *filename){ gridfilename = filename; }
void TasgridWrapper::setOutFilename(const char *filename){ outfilename = filename; }
void TasgridWrapper::setValsFilename(const char *filename){ valsfilename = filename; }
void TasgridWrapper::setXFilename(const char *filename){ xfilename = filename; }
void TasgridWrapper::setAnisoFilename(const char *filename){ anisofilename = filename; }
void TasgridWrapper::setTransformFilename(const char *filename){ transformfilename = filename; }
void TasgridWrapper::setConformalFilename(const char *filename){ conformalfilename = filename; }
void TasgridWrapper::setCustomFilename(const char *filename){ customfilename = filename; }
void TasgridWrapper::setPrintPoints(bool pp){ printCout = pp; }
void TasgridWrapper::setUseASCII(bool ascii){ useASCII = ascii; }

bool TasgridWrapper::checkSane() const{
    bool pass = true;
    if (command == command_none){
        cerr << "ERROR: no command specified" << endl;  return false;
    }else if (command == command_makeglobal){
        if (num_dimensions < 1){ cerr << "ERROR: must specify number of dimensions" << endl; pass = false; }
        if (num_outputs < 0){ cerr << "ERROR: must specify number of outputs (could use zero)" << endl; pass = false; }
        if (depth < 0){ cerr << "ERROR: must specify depth (e.g., level or polynomial degree)" << endl; pass = false; }
        if (depth_type == type_none){ cerr << "ERROR: must specify depth_type (e.g., select levels or polynomial basis)" << endl; pass = false; }
        if (rule == rule_none){ cerr << "ERROR: must specify rule to use (e.g., clenshaw-curtis)" << endl; pass = false; }
        if (!(OneDimensionalMeta::isGlobal(rule))){ cerr << "ERROR: cannot use global grids with rule: " << OneDimensionalMeta::getIORuleString(rule) << endl; pass = false; }
        if ((rule == rule_gaussgegenbauer) || (rule == rule_gausslaguerre) || (rule == rule_gausshermite) || (rule == rule_gaussgegenbauerodd) || (rule == rule_gausshermiteodd) ){
            if (!set_alpha){ cerr << "ERROR: one dimensional rule " << OneDimensionalMeta::getIORuleString(rule) << " requires alpha parameter" << endl; pass = false; }
        }else if (rule == rule_gaussjacobi){
            if (!set_alpha){ cerr << "ERROR: one dimensional rule " << OneDimensionalMeta::getIORuleString(rule) << " requires alpha parameter" << endl; pass = false; }
            if (!set_beta ){ cerr << "ERROR: one dimensional rule " << OneDimensionalMeta::getIORuleString(rule) << " requires beta parameter" << endl; pass = false; }
        }else{
            if (set_alpha){ cerr << "WARNING: alpha parameter set, but one dimensional rule " << OneDimensionalMeta::getIORuleString(rule) << " doesn't depend on alpha" << endl; }
            if (set_beta ){ cerr << "WARNING: beta parameter set, but one dimensional rule " << OneDimensionalMeta::getIORuleString(rule) << " doesn't depend on beta" << endl; }
        }
        if ((rule == rule_customtabulated) && (customfilename == 0)){
            cerr << "ERROR: custom-tabulated rule specified, but no -customflile given" << endl; pass = false;
        }
        if ((rule != rule_customtabulated) && (customfilename != 0)){
            cerr << "WARNING: -customflile given, but is only valid for the custom-tabulated rule and not rule " << OneDimensionalMeta::getIORuleString(rule) << endl;
        }
        if ((gridfilename == 0) && (outfilename == 0) && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -gridfile, -outfile or -print" << endl; pass = false;
        }
        if ((conformal != conformal_none) ^ (conformalfilename != 0)){
            cerr << "WARNING: conformal transform requires both -conformaltype and -conformalfile, ignoring conformal mapping" << endl;
        }
        return pass;
    }else if (command == command_makesequence){
        if (num_dimensions < 1){ cerr << "ERROR: must specify number of dimensions" << endl; pass = false; }
        if (num_outputs < 0){ cerr << "ERROR: must specify number of outputs (could use zero)" << endl; pass = false; }
        if (depth < 0){ cerr << "ERROR: must specify depth (e.g., level or polynomial degree)" << endl; pass = false; }
        if (depth_type == type_none){ cerr << "ERROR: must specify depth_type (e.g., select levels or polynomial basis)" << endl; pass = false; }
        if (rule == rule_none){ cerr << "ERROR: must specify rule to use (e.g., rleja)" << endl; pass = false; }
        if (!(OneDimensionalMeta::isSequence(rule))){ cerr << "ERROR: rule is set to " << OneDimensionalMeta::getIORuleString(rule) << " which is not a sequence rule (e.g., leja, rleja, min/max-lebesgue)" << endl; pass = false; }
        if ((gridfilename == 0) && (outfilename == 0) && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -gridfile, -outfile or -print" << endl; pass = false;
        }
        if ((conformal != conformal_none) ^ (conformalfilename != 0)){
            cerr << "WARNING: conformal transform requires both -conformaltype and -conformalfile, ignoring conformal mapping" << endl;
        }
        return pass;
    }else if (command == command_makelocalp){
        if (num_dimensions < 1){ cerr << "ERROR: must specify number of dimensions" << endl; pass = false; }
        if (num_outputs < 0){ cerr << "ERROR: must specify number of outputs (could use zero)" << endl; pass = false; }
        if (depth < 0){ cerr << "ERROR: must specify depth (e.g., level or polynomial degree)" << endl; pass = false; }
        if (order < -1){ cerr << "ERROR: the maximum order cannot be less than -1"; pass = false; }
        if (rule == rule_none){ cerr << "ERROR: must specify rule to use (e.g., localp)" << endl; pass = false; }
        if (!(OneDimensionalMeta::isLocalPolynomial(rule))){ cerr << "ERROR: cannot use a local polynomial grid with rule: " << OneDimensionalMeta::getIORuleString(rule) << endl; pass = false; }
        if ((gridfilename == 0) && (outfilename == 0) && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -gridfile, -outfile or -print" << endl; pass = false;
        }
        if ((conformal != conformal_none) ^ (conformalfilename != 0)){
            cerr << "WARNING: conformal transform requires both -conformaltype and -conformalfile, ignoring conformal mapping" << endl;
        }
        return pass;
    }else if (command == command_makewavelet){
        if (num_dimensions < 1){ cerr << "ERROR: must specify number of dimensions" << endl; pass = false; }
        if (num_outputs < 0){ cerr << "ERROR: must specify number of outputs (could use zero)" << endl; pass = false; }
        if (depth < 0){ cerr << "ERROR: must specify depth (e.g., level or polynomial degree)" << endl; pass = false; }
        if ((order != 1) && (order != 3)){ cerr << "ERROR: the order must be either 1 or 3"; pass = false; }
        if ((gridfilename == 0) && (outfilename == 0) && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -gridfile, -outfile or -print" << endl; pass = false;
        }
        if ((conformal != conformal_none) ^ (conformalfilename != 0)){
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
                if (!set_alpha){ cerr << "ERROR: one dimensional rule " << OneDimensionalMeta::getIORuleString(rule) << " requires alpha parameter" << endl; pass = false; }
            }
            if (rule == rule_gaussjacobi){
                if (!set_alpha){ cerr << "ERROR: one dimensional rule " << OneDimensionalMeta::getIORuleString(rule) << " requires alpha parameter" << endl; pass = false; }
                if (!set_beta ){ cerr << "ERROR: one dimensional rule " << OneDimensionalMeta::getIORuleString(rule) << " requires beta parameter" << endl; pass = false; }
            }
        }else if (OneDimensionalMeta::isLocalPolynomial(rule)){
            if (order < -1){ cerr << "ERROR: the maximum order cannot be less than -1"; pass = false; }
        }else if (OneDimensionalMeta::isWavelet(rule)){
            if ((order != 1) && (order != 3)){ cerr << "ERROR: the order must be either 1 or 3"; pass = false; }
        }else{
            if (rule == rule_none){
                cerr << "ERROR: must specify rule to use (e.g., clenshaw-curtis or localp)" << endl; pass = false;
            }else{
                cerr << "ERROR: cannot make a quadrature with rule " << OneDimensionalMeta::getIORuleString(rule) << endl; pass = false;
            }
        }
        if ((outfilename == 0) && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -outfile or -print" << endl; pass = false;
        }
        if (gridfilename != 0){  cerr << "WARNING: quadrature does not output a -gridfile, if you need a gridfile use -makeglobal/-makelocalpoly commands followed by -getquadrature" << endl; }
        if ((conformal != conformal_none) ^ (conformalfilename != 0)){
            if ((conformal != conformal_none)) cout << "conf not none" << endl;
            if ((conformalfilename != 0)) cout << "conf file not 0" << endl;
            cerr << "WARNING: conformal transform requires both -conformaltype and -conformalfile, ignoring conformal mapping" << endl;
        }
        return pass;
    }else if (command == command_update){
        if (gridfilename == 0){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        if (depth < 0){ cerr << "ERROR: must specify depth (e.g., level or polynomial degree)" << endl; pass = false; }
        if (depth_type == type_none){ cerr << "ERROR: must specify depth_type (e.g., select levels or polynomial basis)" << endl; pass = false; }
    }else if (command == command_setconformal){
        if (gridfilename == 0){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        if (conformal == conformal_none){  cerr << "ERROR: must specify valid -conformaltype" << endl; pass = false;  }
        if (conformalfilename == 0){ cerr << "ERROR: must specify valid -conformalfile" << endl; pass = false; }
    }else if ((command == command_getquadrature) || (command == command_getpoints) || (command == command_getneeded)){
        if (gridfilename == 0){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        if ((outfilename == 0) && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -outfile or -print" << endl; pass = false;
        }
        return pass;
    }else if (command == command_loadvalues){
        if (gridfilename == 0){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        if (valsfilename == 0){ cerr << "ERROR: must specify valid -valsfile" << endl; pass = false; }
        return pass;
    }else if (command == command_sethierarchical){
        if (gridfilename == 0){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        if (valsfilename == 0){ cerr << "ERROR: must specify valid -valsfile" << endl; pass = false; }
        return pass;
    }else if ((command == command_getinterweights) || (command == command_evaluate) || (command == command_evalhierarchical)){
        if (gridfilename == 0){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        if (xfilename == 0){ cerr << "ERROR: must specify valid -pointsfile" << endl; pass = false; }
        if ((outfilename == 0) && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -outfile or -print" << endl; pass = false;
        }
    }else if (command == command_integrate){
        if (gridfilename == 0){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        if ((outfilename == 0) && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -outfile or -print" << endl; pass = false;
        }
    }else if (command == command_getanisocoeff){
        if (depth_type == type_none){
            cerr << "ERROR: must specify type of coefficients with valid -type!" << endl;
            return false;
        }
        if (gridfilename == 0){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        if ((outfilename == 0) && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -outfile or -print" << endl; pass = false;
        }
    }else if (command == command_refine){
        if (gridfilename == 0){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
    }else if (command == command_refine_aniso){
        if (gridfilename == 0){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
    }else if (command == command_refine_surp){
        if (gridfilename == 0){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
    }else if (command == command_refine_clear){
        if (gridfilename == 0){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
    }else if (command == command_getrefcoeff){
        if (gridfilename == 0){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
    }else if (command == command_getpoly){
        if (gridfilename == 0){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        if (depth_type == type_none){ cerr << "ERROR: must specify depth_type (e.g., select levels or polynomial basis)" << endl; pass = false; }
        if ((depth_type != type_iptotal) && (depth_type != type_ipcurved) && (depth_type != type_iptensor) && (depth_type != type_iphyperbolic) &&
             (depth_type != type_qptotal) && (depth_type != type_qpcurved) && (depth_type != type_qptensor) && (depth_type != type_qphyperbolic)){
            cerr << "ERROR: the type here must start with either i or q indicating whether we seek the polynomils for integration or interpolation." << endl; pass = false;
        }
        if ((outfilename == 0) && (printCout == false)){
            cerr << "ERROR: no means of output are specified, you should specify -outfile or -print" << endl; pass = false;
        }
    }else if (command == command_summary){
        if (gridfilename == 0){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
    }else if (command == command_getsurpluses){
        if (gridfilename == 0){ cerr << "ERROR: must specify valid -gridfile" << endl; pass = false; }
        // ask for an output file
    }

    return pass;
}

void TasgridWrapper::createGlobalGird(){
    int *weights = 0;
    if ((depth_type == type_curved) || (depth_type == type_ipcurved) || (depth_type == type_qpcurved)){
        weights = readAnisotropicFile(2*num_dimensions);
    }else{
        weights = readAnisotropicFile(num_dimensions);
    }
    if (grid == 0) grid = new TasmanianSparseGrid();
    grid->makeGlobalGrid(num_dimensions, num_outputs, depth, depth_type, rule, weights, alpha, beta, customfilename);
    if (weights == 0) delete[] weights;
    if (transformfilename != 0){
        double *transforms = readTransform();
        grid->setDomainTransform(transforms, &(transforms[num_dimensions]));
        delete[] transforms;
    }
    if ((conformal != conformal_none) && (conformalfilename != 0)){
        if (!setConformalTransformation()){
            cerr << "ERROR: could not set conformal transform" << endl;
        }
    }
}
void TasgridWrapper::createSequenceGird(){
    int *weights = 0;
    if ((depth_type == type_curved) || (depth_type == type_ipcurved) || (depth_type == type_qpcurved)){
        weights = readAnisotropicFile(2*num_dimensions);
    }else{
        weights = readAnisotropicFile(num_dimensions);
    }
    if (grid == 0) grid = new TasmanianSparseGrid();
    grid->makeSequenceGrid(num_dimensions, num_outputs, depth, depth_type, rule, weights);
    if (weights == 0) delete[] weights;
    if (transformfilename != 0){
        double *transforms = readTransform();
        grid->setDomainTransform(transforms, &(transforms[num_dimensions]));
        delete[] transforms;
    }
    if ((conformal != conformal_none) && (conformalfilename != 0)){
        if (!setConformalTransformation()){
            cerr << "ERROR: could not set conformal transform" << endl;
        }
    }
}
void TasgridWrapper::createLocalPolynomialGird(){
    if (grid == 0) grid = new TasmanianSparseGrid();
    grid->makeLocalPolynomialGrid(num_dimensions, num_outputs, depth, order, rule);
    if (transformfilename != 0){
        double *transforms = readTransform();
        grid->setDomainTransform(transforms, &(transforms[num_dimensions]));
        delete[] transforms;
    }
    if ((conformal != conformal_none) && (conformalfilename != 0)){
        if (!setConformalTransformation()){
            cerr << "ERROR: could not set conformal transform" << endl;
        }
    }
}
void TasgridWrapper::createWaveletGird(){
    if (grid == 0) grid = new TasmanianSparseGrid();
    grid->makeWaveletGrid(num_dimensions, num_outputs, depth, order);
    if (transformfilename != 0){
        double *transforms = readTransform();
        grid->setDomainTransform(transforms, &(transforms[num_dimensions]));
        delete[] transforms;
    }
    if ((conformal != conformal_none) && (conformalfilename != 0)){
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
    }else{
        cerr << "ERROR: createQuadrature" << endl;
    }
}
bool TasgridWrapper::updateGrid(){
    if (!(grid->isGlobal() || grid->isSequence())){
        cerr << "ERROR: -makeupdate can be called only for Global and Sequence grids" << endl;
        return false;
    }
    int nd = grid->getNumDimensions();
    int *weights = 0;
    if ((depth_type == type_curved) || (depth_type == type_ipcurved) || (depth_type == type_qpcurved)){
        weights = readAnisotropicFile(2*nd);
    }else{
        weights = readAnisotropicFile(nd);
    }
    if (grid->isGlobal()){
        grid->updateGlobalGrid(depth, depth_type, weights);
    }else{
        grid->updateSequenceGrid(depth, depth_type, weights);
    }
    return true;
}
void TasgridWrapper::writeGrid() const{ grid->write(gridfilename, !useASCII); }
bool TasgridWrapper::readGrid(){
    if (grid == 0) grid = new TasmanianSparseGrid();
    char s[3];
    std::ifstream ifs;
    ifs.open(gridfilename, std::ios::in | std::ios::binary);
    ifs.read(s, 3 * sizeof(char));
    ifs.close();
    if ((s[0] == 'T') && (s[1] == 'S') && (s[2] == 'G')){
        return grid->read(gridfilename, true);
    }else{
        return grid->read(gridfilename, false);
    }
}
void TasgridWrapper::outputPoints(bool useNeeded) const{
    int num_p, num_d = grid->getNumDimensions();
    double *points;
    if ((outfilename == 0) && (!printCout)) return;
    if (useNeeded){
        num_p = grid->getNumNeeded();
        points = grid->getNeededPoints();
    }else{
        num_p = grid->getNumPoints();
        points = grid->getPoints();
    }
    if (outfilename != 0) writeMatrix(outfilename, num_p, num_d, points, useASCII);
    if (printCout) printMatrix(num_p, num_d, points);
    delete[] points;
}
void TasgridWrapper::outputQuadrature() const{
    double *points, *weights, *combined;
    if ((outfilename == 0) && (!printCout)) return;
    int num_p = grid->getNumPoints();
    int num_d = grid->getNumDimensions();
    int offset = num_d + 1;
    points  = grid->getPoints();
    weights = grid->getQuadratureWeights();
    combined = new double[num_p * offset];
    for(int i=0; i<num_p; i++){
        combined[i * offset] = weights[i];
        for(int j=0; j<num_d; j++) combined[i * offset + j + 1] = points[i*num_d + j];
    }
    if (outfilename != 0){
        writeMatrix(outfilename, num_p, offset, combined, useASCII);
    }
    if (printCout){
        printMatrix(num_p, offset, combined);
    }
    delete[] combined;
    delete[] weights;
    delete[] points;
}
bool TasgridWrapper::setConformalTransformation(){
    if (conformal == conformal_asin){
        int rows, cols;
        double* mat;
        readMatrix(conformalfilename, rows, cols, mat);
        if (rows != 1){
            cerr << "ERROR: the conformal file for asin should contain only one row" << endl;
            return false;
        }
        if (cols != grid->getNumDimensions()){
            cerr << "ERROR: the conformal file for asin should contain " << grid->getNumDimensions() << " columns, instead it has " << cols << endl;
            return false;
        }
        int *coef = new int[cols];
        for(int j=0; j<cols; j++) coef[j] = (int) mat[j];
        grid->setConformalTransformASIN(coef);
        delete[] mat;
        delete[] coef;
        return true;
    }else{
        cerr << "ERROR: conformal type not implemented" << endl;
        return false;
    }
}
bool TasgridWrapper::loadValues(){
    int rows, cols;
    double *vals = 0;
    readMatrix(valsfilename, rows, cols, vals);
    if (rows != grid->getNumNeeded()){
        cerr << "ERROR: grid is awaiting " << grid->getNumNeeded() << " new values, but " << valsfilename << " specifies " << rows << endl;
        return false;
    }
    if (cols != grid->getNumOutputs()){
        cerr << "ERROR: grid is set for " << grid->getNumOutputs() << " outputs, but " << valsfilename << " specifies " << cols << endl;
        return false;
    }
    grid->loadNeededPoints(vals);
    delete[] vals;
    return true;
}
bool TasgridWrapper::getInterWeights(){
    int rows, cols;
    double *x = 0, *res;
    readMatrix(xfilename, rows, cols, x);
    if (cols != grid->getNumDimensions()){
        cerr << "ERROR: grid is set for " << grid->getNumDimensions() << " dimensions, but " << xfilename << " specifies " << cols << endl;
        return false;
    }
    if (rows < 1){
        cerr << "ERROR: no points specified in " << xfilename << endl;
        return false;
    }
    int num_p = grid->getNumPoints();
    res = new double[num_p * rows];
    #pragma omp parallel for
    for(int i=0; i<rows; i++){
        double *r = grid->getInterpolationWeights(&(x[i*cols]));
        std::copy(r, r + num_p, &(res[i * num_p]));
        delete[] r;
    }
    if (outfilename != 0){
        writeMatrix(outfilename, rows, num_p, res, useASCII);
    }
    if (printCout){
        printMatrix(rows, num_p, res);
    }
    delete[] x;
    delete[] res;
    return true;
}
bool TasgridWrapper::getEvaluate(){
    if (grid->getNumLoaded() == 0){
        cerr << "ERROR: no values loaded in the grid, cannot evaluate!" << endl;
        return false;
    }
    if (grid->getNumOutputs() == 0){
        cerr << "ERROR: no outputs set for the grid, nothing to evaluate!" << endl;
        return false;
    }
    int rows, cols;
    double *x = 0, *res;
    readMatrix(xfilename, rows, cols, x);
    if (cols != grid->getNumDimensions()){
        cerr << "ERROR: grid is set for " << grid->getNumDimensions() << " dimensions, but " << xfilename << " specifies " << cols << endl;
        return false;
    }
    if (rows < 1){
        cerr << "ERROR: no points specified in " << xfilename << endl;
        return false;
    }
    int num_out = grid->getNumOutputs();
    res = new double[num_out * rows];

    grid->evaluateBatch(x, rows, res);

    if (outfilename != 0){
        writeMatrix(outfilename, rows, num_out, res, useASCII);
    }
    if (printCout){
        printMatrix(rows, num_out, res);
    }
    delete[] x;
    delete[] res;
    return true;
}
bool TasgridWrapper::getIntegrate(){
    if (grid->getNumLoaded() == 0){
        cerr << "ERROR: no values loaded in the grid, cannot evaluate!" << endl;
        return false;
    }
    if (grid->getNumOutputs() == 0){
        cerr << "ERROR: no outputs set for the grid, nothing to evaluate!" << endl;
        return false;
    }
    int num_out = grid->getNumOutputs();
    double *q = new double[num_out];
    grid->integrate(q);
    if (outfilename != 0){
        writeMatrix(outfilename, 1, num_out, q, useASCII);
    }
    if (printCout){
        printMatrix(1, num_out, q);
    }
    delete[] q;
    return true;
}
bool TasgridWrapper::getAnisoCoeff(){
    if (grid->getNumOutputs() == 0){
        cerr << "ERROR: cannot estimate coefficients with no outputs!" << endl;
        return false;
    }
    if (grid->getNumLoaded() == 0){
        cerr << "ERROR: cannot estimate coefficients for a grid with no loaded values!" << endl;
        return false;
    }
    int *ab;
    if (grid->isSequence()){
        ab = grid->estimateAnisotropicCoefficients(depth_type, ref_output);
    }else{
        if (ref_output >= grid->getNumOutputs()){
            cerr << "ERROR: -ref_output " << ref_output << " is specified, however, the grid has only " << grid->getNumOutputs() << " outputs!" << endl;
            cerr << " HINT: the outputs are indexed starting at zero!" << endl;
            return false;
        }
        if ((ref_output == -1) && (grid->getNumOutputs() > 1)){
            cerr << "ERROR: must specify a refinement output with -ref_output option!" << endl;
            return false;
        }else if (ref_output == -1) ref_output = 0;
        ab = grid->estimateAnisotropicCoefficients(depth_type, ref_output);
    }
    double *coeff;
    int num_coeff = grid->getNumDimensions();
    if ((depth_type == type_curved) || (depth_type == type_ipcurved) || (depth_type == type_qpcurved)){
        num_coeff *= 2;
    }
    coeff = new double[num_coeff];
    for(int i=0; i<num_coeff; i++){
        coeff[i] = (double) ab[i];
    }
    if (outfilename != 0){
        writeMatrix(outfilename, num_coeff, 1, coeff, useASCII);
    }
    if (printCout){
        printMatrix(num_coeff, 1, coeff);
    }
    delete[] coeff;
    delete[] ab;

    return true;
}

bool TasgridWrapper::refineGrid(){
    // put the sanity check code here, since many of the parameters of the refinement depend on the type of grid being used
    if (grid->getNumOutputs() == 0){
        cerr << "ERROR: cannot refine a grid with no outputs!" << endl;
        return false;
    }
    if (grid->getNumLoaded() == 0){
        cerr << "ERROR: cannot refine a grid with no loaded values!" << endl;
        return false;
    }
    TypeCommand effective_command = command;
    if (command == command_refine){
        if (grid->isGlobal() || grid->isSequence()){
            effective_command = command_refine_aniso;
        }else{
            effective_command = command_refine_surp;
        }
    }
    if (effective_command == command_refine_aniso){
        if ((!grid->isGlobal()) && (!grid->isSequence())){
            cerr << "ERROR: anisotropic refinement can be used only for global and sequence grids!" << endl;
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
        if (grid->isSequence()){
            grid->setAnisotropicRefinement(depth_type, min_growth, ref_output);
        }else{
            if (ref_output >= grid->getNumOutputs()){
                cerr << "ERROR: -ref_output " << ref_output << " is specified, however, the grid has only " << grid->getNumOutputs() << " outputs!" << endl;
                cerr << " HINT: the outputs are indexed starting at zero!" << endl;
                return false;
            }
            if ((ref_output == -1) && (grid->getNumOutputs() > 1)){
                cerr << "ERROR: must specify a refinement output with -ref_output option!" << endl;
                return false;
            }else if (ref_output == -1) ref_output = 0;
            grid->setAnisotropicRefinement(depth_type, min_growth, ref_output);
        }
    }else{ // using surplus refinement
        if (!set_tolerance){
            cerr << "ERROR: must specify -tolerance for surplus refinement!" << endl;
            return false;
        }
        if ((grid->isLocalPolynomial()) || (grid->isWavelet())){
            if (!set_tref){
                cerr << "ERROR: must specify -reftype option!" << endl;
                return false;
            }
            grid->setSurplusRefinement(tolerance, tref, ref_output);
        }else if (grid->isSequence()){
            grid->setSurplusRefinement(tolerance, ref_output);
        }else{
            if (ref_output >= grid->getNumOutputs()){
                cerr << "ERROR: -ref_output " << ref_output << " is specified, however, the grid has only " << grid->getNumOutputs() << " outputs!" << endl;
                cerr << " HINT: the outputs are indexed starting at zero!" << endl;
                return false;
            }
            if ((ref_output == -1) && (grid->getNumOutputs() > 1)){
                cerr << "ERROR: must specify a refinement output with -ref_output option!" << endl;
                return false;
            }else if (ref_output == -1) ref_output = 0;
            grid->setSurplusRefinement(tolerance, ref_output);
        }
    }
    writeGrid();
    return true;
}
bool TasgridWrapper::cancelRefine(){
    grid->clearRefinement();
    return true;
}
bool TasgridWrapper::getPoly(){
    if ((grid->isGlobal()) || (grid->isSequence())){
        int num_d = grid->getNumDimensions();
        bool integrate = ((depth_type == type_iptotal) || (depth_type == type_ipcurved) || (depth_type == type_iptensor) || (depth_type == type_iphyperbolic));
        int n, *poly = 0;
        grid->getGlobalPolynomialSpace(integrate, n, poly);
        double *double_poly = new double[n * num_d];
        for(int i=0; i<num_d * n; i++) double_poly[i] = (double) poly[i];
        if (outfilename != 0){
            writeMatrix(outfilename, n, num_d, double_poly, useASCII);
        }
        if (printCout){
            printMatrix(n, num_d, double_poly);
        }
        delete[] poly;
        delete[] double_poly;
    }else{
        cerr << "ERROR: cannot call -getpoly for a grid that is neither Glboal nor Sequence" << endl;
        return false;
    }
    return true;
}
bool TasgridWrapper::getSummary(){
    grid->printStats();
    return true;
}
bool TasgridWrapper::getSurpluses(){
    const double *surp = grid->getHierarchicalCoefficients();
    int num_p = grid->getNumLoaded();
    int num_o = grid->getNumOutputs();
    if (outfilename != 0){
        writeMatrix(outfilename, num_p, num_o, surp, useASCII);
    }
    if (printCout){
        printMatrix(num_p, num_o, surp);
    }
    return true;
}

bool TasgridWrapper::getEvalHierarchy(){
    int rows, cols;
    double *x = 0, *res;
    readMatrix(xfilename, rows, cols, x);
    if (cols != grid->getNumDimensions()){
        cerr << "ERROR: grid is set for " << grid->getNumDimensions() << " dimensions, but " << xfilename << " specifies " << cols << endl;
        return false;
    }
    if (rows < 1){
        cerr << "ERROR: no points specified in " << xfilename << endl;
        return false;
    }
    int num_p = grid->getNumPoints();
    res = new double[num_p * rows];
    #pragma omp parallel for
    for(int i=0; i<rows; i++){
        double *r = grid->evalHierarchicalFunctions(&(x[i*cols]));
        std::copy(r, r + num_p, &(res[i * num_p]));
        delete[] r;
    }
    if (outfilename != 0){
        writeMatrix(outfilename, rows, num_p, res, useASCII);
    }
    if (printCout){
        printMatrix(rows, num_p, res);
    }
    delete[] x;
    delete[] res;
    return true;
}
bool TasgridWrapper::setHierarchy(){
    int rows, cols;
    double *vals = 0;
    readMatrix(valsfilename, rows, cols, vals);
    if (rows != grid->getNumPoints()){
        cerr << "ERROR: grid is awaiting " << grid->getNumPoints() << " hierarchical surpluses, but " << valsfilename << " specifies " << rows << endl;
        return false;
    }
    if (cols != grid->getNumOutputs()){
        cerr << "ERROR: grid is set for " << grid->getNumOutputs() << " outputs, but " << valsfilename << " specifies " << cols << endl;
        return false;
    }
    cout << "Setting up hierach" << endl;
    grid->setHierarchicalCoefficients(vals);
    delete[] vals;
    return true;
}

bool TasgridWrapper::getPointsIndexes(){
    const int *p = grid->getPointsIndexes();
    int num_p = grid->getNumPoints();
    //num_p = (grid->getNumNeeded() > 0) ? grid->getNumNeeded() : num_p;
    int num_d = grid->getNumDimensions();
    double *pv = new double[num_p * num_d];
    for(int i=0; i<num_p*num_d; i++) pv[i] = (double) (p[i]);
    if (outfilename != 0){
        writeMatrix(outfilename, num_p, num_d, pv, useASCII);
    }
    if (printCout){
        printMatrix(num_p, num_d, pv);
    }
    delete[] pv;
    return true;
}
bool TasgridWrapper::getNeededIndexes(){
    const int *p = grid->getNeededIndexes();
    int num_p = grid->getNumNeeded();
    //num_p = (grid->getNumNeeded() > 0) ? grid->getNumNeeded() : num_p;
    int num_d = grid->getNumDimensions();
    double *pv = new double[num_p * num_d];
    for(int i=0; i<num_p*num_d; i++) pv[i] = (double) (p[i]);
    if (outfilename != 0){
        writeMatrix(outfilename, num_p, num_d, pv, useASCII);
    }
    if (printCout){
        printMatrix(num_p, num_d, pv);
    }
    delete[] pv;
    return true;
}

int* TasgridWrapper::readAnisotropicFile(int num_weights) const{
    if (anisofilename == 0) return 0;
    int rows, cols;
    double* mat;
    readMatrix(anisofilename, rows, cols, mat);
    if (rows != 1){
        cerr << "ERROR: anisotropy file must contain only one row" << endl;
        delete[] mat;
        exit(1);
    }
    if (cols != num_weights){
        cerr << "ERROR: anisotropy file has wrong number of entries, " << num_weights << " expected " << cols << " found." << endl;
        delete[] mat;
        exit(1);
    }
    int *weights = new int[num_weights];
    for(int i=0; i<num_weights; i++) weights[i] = (int)(mat[i] + 0.3);
    return weights;
}
double* TasgridWrapper::readTransform() const{
    int rows, cols;
    double* mat;
    readMatrix(transformfilename, rows, cols, mat);
    if (cols != 2){
        cerr << "ERROR: file " << transformfilename << " must have exactly two columns." << endl;
        delete mat;
        exit(1);
    }
    if (rows != num_dimensions){
        cerr << "ERROR: file " << transformfilename << " has " << rows << " rows, instead of the number of dimensions " << num_dimensions << endl;
        delete mat;
        exit(1);
    }
    double *transforms = new double[2*num_dimensions];
    for(int i=0; i<num_dimensions; i++){
        transforms[i] = mat[2*i];
        transforms[num_dimensions + i] = mat[2*i+1];
    }
    delete[] mat;
    return transforms;
}

void TasgridWrapper::readMatrix(const char *filename, int &rows, int &cols, double* &mat){
    if (mat == 0) delete[] mat;
    std::ifstream ifs;
    ifs.open(filename, std::ios::in | std::ios::binary);
    if (!(ifs.good())){
        cerr << "ERROR: could not open file " << filename << endl;
        mat = 0;
        ifs.close();
        return;
    }
    char tsg[3] = {'A', 'A', 'A'};
    ifs.read(tsg, 3*sizeof(char));
    if ((tsg[0] == 'T') && (tsg[1] == 'S') && (tsg[2] == 'G')){
        ifs.read((char*) &rows, sizeof(int));
        ifs.read((char*) &cols, sizeof(int));
        if ((rows == 0) && (cols == 0)){
            cerr << "WARNING: empty file " << filename << endl;
            mat = 0;
            ifs.close();
            return;
        }
        mat = new double[rows*cols];
        ifs.read((char*) mat, rows * cols * sizeof(double));
    }else{ // not a binary file
        ifs.close();
        ifs.open(filename);
        ifs >> rows >> cols;
        if ((rows == 0) && (cols == 0)){
            cerr << "WARNING: empty file " << filename << endl;
            mat = 0;
            ifs.close();
            return;
        }
        mat = new double[rows*cols];
        for(int i=0; i<rows; i++){
            for(int j=0; j<cols; j++){
                ifs >> mat[i*cols + j];
            }
        }
    }
    ifs.close();
}
void TasgridWrapper::writeMatrix(const char *filename, int rows, int cols, const double mat[], bool ascii){
    std::ofstream ofs;
    if (ascii){
        ofs.open(filename);
        ofs << rows << " " << cols << endl;
        ofs.precision(17);
        ofs << std::scientific;
        for(int i=0; i<rows; i++){
            for(int j=0; j<cols; j++){
                ofs << setw(25) << mat[i*cols + j] << " ";
            }
            ofs << endl;
        }
    }else{
        ofs.open(filename, std::ios::out | std::ios::binary);
        char tsg[3] = {'T', 'S', 'G'};
        ofs.write(tsg, 3*sizeof(char));
        ofs.write((char*) &rows, sizeof(int));
        ofs.write((char*) &cols, sizeof(int));
        ofs.write((char*) mat, rows * cols * sizeof(double));
    }
    ofs.close();
}
void TasgridWrapper::printMatrix(int rows, int cols, const double mat[]){
    cout << rows << " " << cols << endl;
    cout.precision(17);
    cout << std::scientific;
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
            cout << setw(25) << mat[i*cols + j] << " ";
        }
        cout << endl;
    }
}

bool TasgridWrapper::executeCommand(){
    if (!checkSane()) return false;
    if (command == command_makeglobal){
        createGlobalGird();
    }else if (command == command_makesequence){
        createSequenceGird();
    }else if (command == command_makelocalp){
        createLocalPolynomialGird();
    }else if (command == command_makewavelet){
        createWaveletGird();
    }else if (command == command_makequadrature){
        createQuadrature();
        outputQuadrature();
    }
    if (isCreateCommand(command)){
        if (command != command_makequadrature){
            outputPoints(false);
            if (gridfilename != 0) writeGrid();
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
    }else if (command == command_getquadrature){
        outputQuadrature();
    }else if (command == command_getpoints){
        outputPoints(false);
    }else if (command == command_getneeded){
        outputPoints(true);
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
        cancelRefine();
    }else if (command == command_getpoly){
        if (!getPoly()){
            cerr << "ERROR: could not get polynomial basis" << endl;
            return false;
        }
    }else if (command == command_summary){
        getSummary();
    }else if (command == command_evalhierarchical){
        if (!getEvalHierarchy()){
            cerr << "ERROR: could not evaluate the hierarchical basis functions" << endl;
            return false;
        }
    }else if (command == command_sethierarchical){
        if (setHierarchy()){
            writeGrid();
        }else{
            cerr << "ERROR: could not set the hierarchical surpluses" << endl;
            return false;
        }
    }else if (command == command_getsurpluses){
        if (!getSurpluses()){
            cerr << "ERROR: could not get the surpluses" << endl;
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
