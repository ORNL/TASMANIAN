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


// helper class
struct internal_sparse_matrix{
    int num_cols;
    std::vector<int> pntr, indx;
    std::vector<double> vals;
    internal_sparse_matrix(int cols) : num_cols(cols){}
    template<bool iomode>
    void writeSparseMatrix(std::ostream &os) const{
        int nnz = static_cast<int>(indx.size());
        int rows = static_cast<int>(pntr.size() - 1);
        IO::writeNumbers<iomode, IO::pad_line>(os, rows, num_cols, nnz);
        IO::writeVector<iomode, IO::pad_line>(pntr, os);
        IO::writeVector<iomode, IO::pad_line>(indx, os);
        IO::writeVector<iomode, IO::pad_line>(vals, os);
    }
    void write(std::string const &filename, bool use_ascii) const{
        if (not filename.empty()){
            if (use_ascii){
                std::ofstream ofs(filename);
                ofs << std::scientific; ofs.precision(17);
                writeSparseMatrix<mode_ascii>(ofs);
            }else{
                std::ofstream ofs(filename, std::ios::out | std::ios::binary);
                char charTSG[3] = {'T', 'S', 'G'};
                ofs.write(charTSG, 3 * sizeof(char));
                writeSparseMatrix<mode_binary>(ofs);
            }
        }
    }
    void write(bool print_to_cout) const{
        if (print_to_cout){
            cout << std::scientific; cout.precision(17);
            writeSparseMatrix<mode_ascii>(cout);
        }
    }
};

TasgridWrapper::TasgridWrapper() : command(command_none), num_dimensions(0), num_outputs(-1), depth(-1), order(1),
    depth_type(type_none), rule(rule_none),
    conformal(conformal_none), alpha(0.0), beta(0.0), set_alpha(false), set_beta(false), tolerance(0.0), set_tolerance(false),
    ref_output(-1), min_growth(-1), tref(refine_fds), set_tref(false),
    printCout(false), useASCII(false), set_gpuid(-2), shift(0.0), set_shift(false)
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
            {"-getdiffweights", command_getdiffweights}, {"-gd",   command_getdiffweights},
            {"-getpoints",  command_getpoints},  {"-gp", command_getpoints},
            {"-getneeded",  command_getneeded},  {"-gn", command_getneeded},
            {"-loadvalues", command_loadvalues}, {"-l",  command_loadvalues},
            {"-evaluate",   command_evaluate},   {"-e",  command_evaluate},
            {"-integrate",  command_integrate},  {"-i",  command_integrate},
            {"-differentiate",  command_differentiate},  {"-d",  command_differentiate},
            {"-evalhierarchyd", command_evalhierarchical_dense},  {"-ehd", command_evalhierarchical_dense},
            {"-evalhierarchys", command_evalhierarchical_sparse}, {"-ehs", command_evalhierarchical_sparse},
            {"-gethsupport", command_gethsupport}, {"-ghsup", command_gethsupport},
            {"-getanisotropy", command_getanisocoeff}, {"-ga", command_getanisocoeff},
            {"-refinesurp",    command_refine_surp},   {"-rs", command_refine_surp},
            {"-refineaniso",   command_refine_aniso},  {"-ra", command_refine_aniso},
            {"-refine",        command_refine},        {"-r",  command_refine},
            {"-cancelrefine",  command_refine_clear},  {"-cr",   command_refine_clear},
            {"-mergerefine",   command_refine_merge},  {"-mr",   command_refine_merge},
            {"-using-construct", command_using_construct},
            {"-getconstructpnts", command_get_candidate_construction}, {"-gcp", command_get_candidate_construction},
            {"-loadconstructed", command_load_construction}, {"-lcp", command_load_construction},
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

struct command_tester{
    TypeCommand command;
    bool inside(std::vector<TypeCommand> const &command_list1, std::vector<TypeCommand> const &command_list2 = {}){
        return (std::any_of(command_list1.begin(), command_list1.end(), [&](TypeCommand c)->bool{ return (c == command); })
            or std::any_of(command_list2.begin(), command_list2.end(), [&](TypeCommand c)->bool{ return (c == command); }));
    }
};
struct test_result_wrapper{
    bool pass = true;
    operator bool() const{ return pass; }
    void fail_if(bool condition, const char *text){
        if (condition){
            cerr << "ERROR: " << text << "\n";
            pass = false;
        }
    }
    void worry_if(bool condition, const char *text){
        if (condition) cerr << "WARNING: " << text << "\n";
    }
};

bool TasgridWrapper::checkSane() const{
    if (command == command_none){
        cerr << "ERROR: no command specified\n";
        return false;
    }

    command_tester com{command};
    std::vector<TypeCommand> makecoms = {command_makeglobal, command_makesequence, command_makelocalp, command_makewavelet, command_makefourier};

    test_result_wrapper test;
    // adopt the signature for the checks, problem and commands that have the problem
    // e.g., num_dimensions < 1 is problem when using make-grid command or make quadrature
    test.fail_if(num_dimensions < 1 and com.inside(makecoms, {command_makequadrature}),
                 "must specify number of dimensions (e.g., number of model inputs)");
    test.fail_if(num_outputs < 1 and com.inside(makecoms),
                 "must specify number of outputs (could be zero)");
    test.fail_if(depth < 0 and com.inside(makecoms, {command_makequadrature, command_makeexoquad, command_update}),
                 "must specify depth (e.g., level or polynomial degree)");
    test.fail_if(order < -1 and (command == command_makelocalp
                                 or (command == command_makequadrature and OneDimensionalMeta::isLocalPolynomial(rule))),
                 "the maximum order for local polynomial cannot be less than -1");
    test.fail_if(order != 1 and order != 3 and (command == command_makewavelet
                                                or (command == command_makequadrature and OneDimensionalMeta::isWavelet(rule))),
                 "the wavelet order must be either 1 or 3");

    test.fail_if(depth_type == type_none and
                 (com.inside({command_makeglobal, command_makesequence, command_makefourier,
                              command_update, command_getanisocoeff, command_getpoly})
                  or (command == command_makequadrature and (OneDimensionalMeta::isGlobal(rule) or OneDimensionalMeta::isFourier(rule)))),
                 "must specify depth_type (e.g., select levels or polynomial basis)");
    test.fail_if(rule == rule_none and com.inside({command_makeglobal, command_makesequence, command_makelocalp, command_makequadrature}),
                 "must specify rule to use (e.g., clenshaw-curtis or localp)");

    if (command == command_makeglobal or command == command_makequadrature){
        bool needs_alpha = (rule == rule_gaussgegenbauer or rule == rule_gausslaguerre or rule == rule_gausshermite or
                            rule == rule_gaussgegenbauerodd or rule == rule_gausshermiteodd or rule == rule_gaussjacobi);
        bool needs_beta = (rule == rule_gaussjacobi);

        test.fail_if(not set_alpha and needs_alpha,
                     (std::string("one dimensional rule ") + IO::getRuleString(rule) + " requires alpha parameter").c_str());
        test.fail_if(not set_beta and needs_beta,
                     (std::string("one dimensional rule ") + IO::getRuleString(rule) + " requires alpha parameter").c_str());

        test.worry_if(set_alpha and not needs_alpha,
                      (std::string("alpha parameter set, but one dimensional rule ") + IO::getRuleString(rule) + " doesn't depend on alpha").c_str());
        test.worry_if(set_beta and not needs_beta,
                      (std::string("beta parameter set, but one dimensional rule ") + IO::getRuleString(rule) + " doesn't depend on beta").c_str());

        test.fail_if(customfilename.empty() and rule == rule_customtabulated,
                     "ustom-tabulated rule specified, but no -customflile given");
        test.worry_if(not customfilename.empty() and rule != rule_customtabulated,
                      "custom-tabulated rule specified, but no -customflile given");
    }

    // not output at all
    test.fail_if(gridfilename.empty() and outfilename.empty() and not printCout and com.inside(makecoms),
                 "no means of output are specified, you should specify -gridfile, -outfile or -print");
    // cannot output to gridfile and outfile and print are not set
    test.fail_if(outfilename.empty() and not printCout
                 and com.inside({command_makequadrature, command_makeexoquad, command_getinterweights, command_evaluate,
                                 command_evalhierarchical_dense, command_evalhierarchical_sparse, command_differentiate,
                                 command_get_candidate_construction, command_getquadrature, command_getpoints,
                                 command_getneeded, command_getcoefficients, command_gethsupport, command_integrate,
                                 command_getanisocoeff, command_getpoly}),
                 "no means of output are specified, you should specify -outfile or -print");

    test.fail_if(conformalfilename.empty() and conformal != conformal_none,
                 "conformal transform requires both -conformaltype and -conformalfile");

    test.fail_if(conformal == conformal_none and command == command_setconformal,
                 "must specify valid -conformaltype");
    test.fail_if(conformalfilename.empty() and command == command_setconformal,
                 "must specify valid -conformalfile");

    test.fail_if(gridfilename.empty() and not com.inside(makecoms, {command_makequadrature, command_makeexoquad}),
                 "must specify valid -gridfile");
    test.fail_if(xfilename.empty()
                 and com.inside({command_getinterweights, command_evaluate, command_evalhierarchical_dense,
                                 command_evalhierarchical_sparse, command_differentiate, command_load_construction}),
                 "must specify valid -pointsfile");

    test.fail_if(valsfilename.empty() and com.inside({command_loadvalues, command_setcoefficients, command_load_construction}),
                 "must specify valid -valsfile");

    // handle special cases per command
    switch(command){
        case command_makeglobal:
            test.fail_if(not OneDimensionalMeta::isGlobal(rule),
                         (std::string("cannot use global grids with rule: ") + IO::getRuleString(rule)).c_str());
            break;
        case command_makesequence:
            test.fail_if(not OneDimensionalMeta::isSequence(rule),
                         (std::string("rule is set to ") + IO::getRuleString(rule) + " which is not a sequence rule (e.g., leja, rleja, min/max-lebesgue)").c_str());
            break;
        case command_makelocalp:
            test.fail_if(not OneDimensionalMeta::isLocalPolynomial(rule),
                         (std::string("cannot use local polynomial grids with rule: ") + IO::getRuleString(rule)).c_str());
            break;
        case command_makeexoquad:
            test.fail_if(not set_shift, "must specify shift parameter");
            test.fail_if(weightfilename.empty(), "must specify shift parameter");
            test.fail_if(description.empty(), "must specify description string");
            break;
        case command_getpoly:
            test.fail_if(depth_type == type_level or depth_type == type_curved or depth_type == type_hyperbolic,
                         "the type here must start with either i or q indicating whether we seek the polynomils for integration or interpolation");
            break;
        default:
            break;
    }

    test.worry_if(num_outputs != -1 and command == command_makequadrature,
                  "ignoring the -outputs specified for the -makequadrature command");
    test.worry_if(not gridfilename.empty() and command == command_makequadrature,
                  "quadrature does not output a -gridfile, if you need a gridfile use -makeglobal/-makelocalpoly commands followed by -getquadrature");
    test.worry_if(not conformalfilename.empty() and conformal == conformal_none,
                  "conformal transform requires both -conformaltype and -conformalfile, ignoring conformal mapping");

    return test;
}
bool TasgridWrapper::checkSanePostRead() const{
    command_tester com{command};
    test_result_wrapper test;

    test.fail_if(grid.getNumLoaded() == 0 and
                 com.inside({command_evaluate, command_differentiate, command_integrate, command_getanisocoeff,
                             command_getcoefficients, command_refine, command_refine_aniso, command_refine_surp}),
                 "the grid has no loaded data");
    test.fail_if(grid.getNumOutputs() == 0 and
                 com.inside({command_evaluate, command_differentiate, command_integrate, command_getanisocoeff,
                             command_getcoefficients, command_refine, command_refine_aniso, command_refine_surp,
                             command_get_candidate_construction}),
                 "the grid has no outputs");

    test.fail_if(ref_output >= grid.getNumOutputs() and
                 com.inside({command_getanisocoeff, }),
                 "-refout outside of the range, note the outputs are indexed from zero");
    test.fail_if(ref_output == -1 and grid.getNumOutputs() > 1 and grid.isGlobal() and
                 com.inside({command_getanisocoeff}),
                 "-refout cannot use -1 when working with Global grids");

    test.fail_if((grid.isLocalPolynomial() or grid.isWavelet()) and command == command_refine_aniso,
                 "anisotropic operations are available only for Global, Sequence and Fourier grids");
    test.fail_if((grid.isFourier() and command == command_refine_surp),
                 "surplus refinement cannot be applied to Fourier grids");

    bool is_refine = (command == command_refine or command == command_get_candidate_construction);
    if (command == command_refine_aniso or (is_refine and (grid.isGlobal() or grid.isSequence() or grid.isFourier()))){
        test.fail_if(grid.isLocalPolynomial() or grid.isWavelet(),
                     "anisotropic refinement can be applied only to Global, Sequence and Fourier grids");
        test.fail_if(depth_type == type_none, "anisotropic refinement requires depth type with -tt");
        test.fail_if(ref_output == -1 and grid.getNumOutputs() > 1 and grid.isGlobal() and
                     com.inside({command_getanisocoeff}),
                     "-refout cannot use -1 when working with Global grids");

        test.worry_if(set_tolerance, "anisotropic refinement (and grids Global, Sequence, Fourier) ignores the -tolerance option");
        test.worry_if(set_tref, "anisotropic refinement (and grids Global, Sequence, Fourier) ignores the -reftype option");
    }

    if (command == command_refine_surp or (is_refine and (grid.isLocalPolynomial() or grid.isWavelet()))){
        test.fail_if(not set_tolerance, "must specify -tolerance for surplus refinement");
        test.fail_if(not set_tref, "must specify -reftype option");

        test.worry_if(not valsfilename.empty() and (grid.isGlobal() or grid.isSequence() or grid.isFourier()),
                      "the scale factors are not used with Global, Sequence and Fourier grids");
    }

    test.fail_if((grid.isLocalPolynomial() or grid.isWavelet()) and command == command_getpoly,
                 "cannot call -getpoly for a grid that is neither Global nor Sequence");

    return test;
}
void TasgridWrapper::iassert(bool condition, const char* text) const{
    if (not condition){
        cerr << "ERROR: " << text << "\n";
        pass_flag = false;
    }
}

bool TasgridWrapper::executeCommand(){
    if (not checkSane()) return false;
    pass_flag = true;
    if (command == command_makeexoquad){
        createExoticQuadrature();
        return pass_flag;
    }
    if (set_gpuid > -2){
        try{
            grid.enableAcceleration(accel_gpu_cuda, set_gpuid);
        }catch(std::runtime_error &e){
            cerr << "WARNING: setting the GPU failed with the following message: \n";
            cerr << e.what() << "\n";
        }
    }

    command_tester com{command};
    std::vector<TypeCommand> makecoms = {command_makeglobal, command_makesequence, command_makelocalp, command_makewavelet,
                                         command_makefourier, command_makequadrature};
    std::vector<TypeCommand> constcoms = {
        command_getquadrature, command_getinterweights, command_getdiffweights, command_getpoints,
        command_getneeded, command_evaluate, command_integrate, command_differentiate, command_getanisocoeff,
        command_getpoly, command_summary, command_getcoefficients, command_evalhierarchical_sparse,
        command_evalhierarchical_dense, command_gethsupport, command_getpointsindex, command_getneededindex
    };
    // read grid or make a new grid
    if (not com.inside(makecoms, {command_makeexoquad})){
        if (not readGridfile()) return false;
    }

    if (not checkSanePostRead()) return false;

    if (command == command_makequadrature) num_outputs = 0;

    if (com.inside(makecoms)){
        auto llimits = readLimits();
        auto aniso = readAnisotropic();

        if (command == command_makeglobal or (command == command_makequadrature and OneDimensionalMeta::isGlobal(rule))){
            grid.makeGlobalGrid(num_dimensions, num_outputs, depth, depth_type, rule, aniso, alpha, beta, customfilename.c_str(), llimits);
        }else if (command == command_makesequence){
            grid.makeSequenceGrid(num_dimensions, num_outputs, depth, depth_type, rule, aniso, llimits);
        }else if (command == command_makefourier or (command == command_makequadrature and rule == rule_fourier)){
            grid.makeFourierGrid(num_dimensions, num_outputs, depth, depth_type, aniso, llimits);
        }else if (command == command_makelocalp or (command == command_makequadrature and OneDimensionalMeta::isGlobal(rule))){
            grid.makeLocalPolynomialGrid(num_dimensions, num_outputs, depth, order, rule, llimits);
        }else{ // wavelets
            grid.makeWaveletGrid(num_dimensions, num_outputs, depth, order, llimits);
        }
        setTransform();
    }

    if (com.inside(makecoms) or command == command_setconformal)
        setConformal();

    switch(command){
        case command_update:
            grid.updateGrid(depth, depth_type, readAnisotropic());
            break;
        case command_getdiffweights:
        case command_getinterweights:
        case command_evalhierarchical_dense:
        case command_evalhierarchical_sparse:
        case command_evaluate:
        case command_differentiate:
            processEvalLike();
            break;
        case command_integrate:
        case command_gethsupport:
        case command_getanisocoeff:
            processOutputLike();
            break;
        case command_getcoefficients:
            outputHierarchicalCoefficients();
            break;
        case command_loadvalues:
        case command_load_construction:
            loadComputedValues();
            break;
        case command_setcoefficients:
            setHierarchy();
            break;
        case command_refine_clear:
            grid.clearRefinement();
            if (grid.isUsingConstruction())
                grid.finishConstruction();
            break;
        case command_refine_merge:
            grid.mergeRefinement();
            break;
        case command_using_construct:
            cout << "dynamic construction: " << ((grid.isUsingConstruction()) ? "enabled" : "disabled") << "\n";
            break;
        case command_summary:
            grid.printStats();
            break;
        case command_getneededindex:
        case command_getpointsindex:
            outputIndexes((command == command_getneededindex) ? output_points_mode::needed : output_points_mode::regular);
            break;
        case command_getpoly:
            getPoly();
            break;
        case command_refine:
        case command_refine_aniso:
        case command_refine_surp:
            refineGrid();
            break;
        case command_get_candidate_construction:
            getConstructedPoints();
            break;
        default:
            break;
    }

    if (com.inside(makecoms) or command == command_getpoints)
        outputPoints(output_points_mode::regular);
    if (com.inside({command_getneeded, command_refine, command_refine_aniso, command_refine_surp}))
        outputPoints(output_points_mode::needed);
    if (command == command_makequadrature or command == command_getquadrature)
        outputQuadrature();

    if (not com.inside(constcoms))
        writeGrid();

    return pass_flag;
}
void TasgridWrapper::createExoticQuadrature(){
    TasGrid::TasmanianSparseGrid weight_surrogate;
    weight_surrogate.read(weightfilename.c_str());
    iassert(weight_surrogate.getNumDimensions() == 1, "the weight function surrogate must be one-dimensional");
    iassert(weight_surrogate.getNumLoaded() > 0, "the weight function surrogate must have loaded values for interpolation");
    ct = TasGrid::getExoticQuadrature(depth, shift, weight_surrogate, description.c_str(), is_symmetric_weight_function);
    if (not outfilename.empty()){
        std::ofstream ofs(outfilename, std::ios::out | std::ios::trunc);
        ct.write<mode_ascii>(ofs);
    }
    if (printCout)
        ct.write<mode_ascii>(cout);
}
bool TasgridWrapper::readGridfile(){
    try{
        grid.read(gridfilename);
    }catch(std::runtime_error &e){
        cerr << e.what() << "\n";
        return false;
    }
    num_dimensions = grid.getNumDimensions();
    num_outputs = grid.getNumOutputs();
    rule = grid.getRule();
    return true;
}
std::vector<int> TasgridWrapper::readLimits() const{
    if (levellimitfilename.empty()) return std::vector<int>();
    auto mat = readMatrix(levellimitfilename);
    iassert(mat.getNumStrips() == 1, "level limits file must contain only one row");
    iassert(static_cast<int>(mat.getStride()) == num_dimensions,
            (std::string("level limits file has wrong number of entries, expected: ") + std::to_string(num_dimensions) +
             " but found " + std::to_string(mat.getStride())).c_str());

    std::vector<int> llimits(num_dimensions);
    std::transform(mat.begin(), mat.end(), llimits.begin(), [](double x)->int{ return static_cast<int>(x); });
    return llimits;
}
std::vector<int> TasgridWrapper::readAnisotropic() const{
    if (anisofilename.empty()) return std::vector<int>();
    auto mat = readMatrix(anisofilename);
    iassert(mat.getNumStrips() == 1, "anisotropy file must contain only one row");
    size_t expected_size = static_cast<size_t>(
        (OneDimensionalMeta::getControurType(depth_type) == type_curved) ? 2*num_dimensions : num_dimensions);
    iassert(mat.getStride() == expected_size,
            (std::string("level limits file has wrong number of entries, expected: ") + std::to_string(expected_size) +
             " but found " + std::to_string(mat.getStride())).c_str());
    std::vector<int> weights(expected_size);
    std::transform(mat.begin(), mat.end(), weights.begin(), [](double x)->int{ return static_cast<int>(x); });
    return weights;
}
void TasgridWrapper::setTransform(){
    if (transformfilename.empty()) return;
    auto mat = readMatrix(transformfilename);
    iassert(mat.getStride() == 2, "the matrix in the transform file must have exactly two columns");
    iassert(mat.getNumStrips() == num_dimensions,
            (std::string("the domain transform expects ") + std::to_string(num_dimensions) +
             " rows but found " + std::to_string(mat.getNumStrips()) + " in the file: " + transformfilename).c_str());
    std::vector<double> transa((size_t) num_dimensions);
    std::vector<double> transb((size_t) num_dimensions);
    for(int i=0; i<num_dimensions; i++){
        transa[i] = mat.getStrip(i)[0];
        transb[i] = mat.getStrip(i)[1];
    }
    grid.setDomainTransform(transa, transb);
}
void TasgridWrapper::setConformal(){
    if (conformal == conformal_asin){
        auto mat = readMatrix(conformalfilename);
        iassert(mat.getNumStrips() == 1, "the conformal file for asin should contain only one row");
        iassert(static_cast<int>(mat.getStride()) == num_dimensions,
                (std::string("conformal file for asin wrong number of entries, expected: ") + std::to_string(num_dimensions) +
                 " but found " + std::to_string(mat.getStride())).c_str());
        std::vector<int> coeff(mat.getTotalEntries());
        std::transform(mat.begin(), mat.end(), coeff.begin(), [](double x)->int{ return static_cast<int>(x); });
        grid.setConformalTransformASIN(coeff);
    }
}
void TasgridWrapper::outputPoints(output_points_mode mode) const{
    if (outfilename.empty() and not printCout) return;
    int num_points = (mode == output_points_mode::needed) ? grid.getNumNeeded() : grid.getNumPoints();
    auto points = (mode == output_points_mode::needed) ? grid.getNeededPoints() : grid.getPoints();
    writeMatrix(outfilename, num_points, num_dimensions, points.data());
    printMatrix(num_points, num_dimensions, points.data());
}
void TasgridWrapper::outputIndexes(output_points_mode mode) const{
    const int *p = (mode == output_points_mode::needed) ? grid.getNeededIndexes() : grid.getPointsIndexes();
    int num_points = (mode == output_points_mode::needed) ? grid.getNumNeeded() : grid.getNumPoints();
    Data2D<double> pv(num_dimensions, num_points);
    std::transform(p, p + Utils::size_mult(num_points, num_dimensions), pv.getStrip(0), [](int i)->double{ return double(i); });

    writeMatrix(outfilename, num_points, num_dimensions, pv.getStrip(0));
    printMatrix(num_points, num_dimensions, pv.getStrip(0));
}
void TasgridWrapper::outputQuadrature() const{
    if (outfilename.empty() and not printCout) return;
    int num_points = grid.getNumPoints();
    auto pnts  = grid.getPoints();
    Utils::Wrapper2D<double> points(num_dimensions, pnts.data());
    auto weights = grid.getQuadratureWeights();
    Data2D<double> combined(num_dimensions + 1, num_points);
    for(int i=0; i<num_points; i++){
        combined.getStrip(i)[0] = weights[i];
        std::copy_n(points.getStrip(i), num_dimensions, &(combined.getStrip(i)[1]));
    }
    writeMatrix(outfilename, num_points, num_dimensions+1, combined.getStrip(0));
    printMatrix(num_points, num_dimensions+1, combined.getStrip(0));
}
void TasgridWrapper::writeGrid() const{
    if (gridfilename.empty()) return;
    grid.write(gridfilename.c_str(), (useASCII) ? mode_ascii : mode_binary);
}
Data2D<double> TasgridWrapper::verifiedRead(std::string const& filename, int expected_stride) const{
    if (filename.empty() or expected_stride == 0) return Data2D<double>();
    Data2D<double> x = readMatrix(filename);
    iassert(x.getStride() == static_cast<size_t>(expected_stride),
            (std::string("the matrix in file ") + filename + " has " + std::to_string(x.getStride()) +
             " rows, but it should have " + std::to_string(expected_stride)).c_str());
    return x;
}

void TasgridWrapper::processEvalLike() const{
    int num_points = grid.getNumPoints();
    if (not pass_flag) return;

    auto x = verifiedRead(xfilename, num_dimensions);

    Data2D<double> result;
    switch(command){
        case command_evaluate:
            result = Data2D<double>(num_outputs, x.getNumStrips());
            grid.evaluateBatch(x.data(), x.getNumStrips(), result.data());
            break;
        case command_differentiate:
            result = Data2D<double>(num_outputs * num_dimensions, x.getNumStrips());
            #pragma omp parallel for
            for (int i=0; i<x.getNumStrips(); i++)
                grid.differentiate(x.getStrip(i), result.getStrip(i));
            break;
        case command_getinterweights:
            result = Data2D<double>(num_points, x.getNumStrips());
            #pragma omp parallel for
            for(int i=0; i<x.getNumStrips(); i++)
                grid.getInterpolationWeights(x.getStrip(i), result.getStrip(i));
            break;
        case command_getdiffweights:
            result = Data2D<double>(num_points * num_dimensions, x.getNumStrips());
            #pragma omp parallel for
            for(int i=0; i<x.getNumStrips(); i++)
                grid.getDifferentiationWeights(x.getStrip(i), result.getStrip(i));
            break;
        case command_evalhierarchical_dense:
            result = Data2D<double>(num_points * (grid.isFourier() ? 2 : 1), x.getNumStrips());
            grid.evaluateHierarchicalFunctions(x.data(), x.getNumStrips(), result.data());
            break;
        case command_evalhierarchical_sparse:
        { // scope is needed to declare sparse variables
            internal_sparse_matrix matrix(num_points);
            grid.evaluateSparseHierarchicalFunctions(x.release(), matrix.pntr, matrix.indx, matrix.vals);
            matrix.write(outfilename, useASCII);
            matrix.write(printCout);
            return; // after this, there is output section but invalid for this command
        }
            break;
        default:
            throw std::runtime_error("ERROR: internal problem, processEvalLike() called with wrong command");
            break;
    }

    writeMatrix(outfilename, result);
    printMatrix(result);
}
void TasgridWrapper::processOutputLike() const{
    // returns data and updates the above through capture
    Data2D<double> result;
    switch(command){
        case command_integrate:
            if (not pass_flag) return;
            result = Data2D<double>(1, grid.getNumOutputs(), grid.integrate());
            break;
        case command_gethsupport:
            result = Data2D<double>(grid.getNumDimensions(), grid.getNumPoints(), grid.getHierarchicalSupport());
            break;
        case command_getanisocoeff:
            if (not pass_flag) return;
            {
                std::vector<int> ab = grid.estimateAnisotropicCoefficients(depth_type,
                                                                           (grid.isGlobal() and ref_output == -1) ? 0 : ref_output);
                std::vector<double> temp(ab.size());
                std::transform(ab.begin(), ab.end(), temp.begin(), [](int i)->double{ return double(i); });
                result = Data2D<double>(1, grid.getNumDimensions() * ((OneDimensionalMeta::getControurType(depth_type) == type_curved) ? 2 : 1), std::move(temp));
            }
            break;
        default:
            throw std::runtime_error("invalid command for processOutputLike() method");
    }
    writeMatrix(outfilename, result);
    printMatrix(result);
}
void TasgridWrapper::outputHierarchicalCoefficients() const{
    const double *coeff = grid.getHierarchicalCoefficients();
    int num_points = grid.getNumPoints();
    if (grid.isFourier()){
        // use interwoven format for complex coefficients
        Data2D<double> coeff_fourier(2 * num_outputs, num_points);
        Utils::Wrapper2D<const double> real(num_outputs, coeff);
        Utils::Wrapper2D<const double> imag(num_outputs, real.getStrip(num_points));
        for(int p=0; p<num_points; p++){
            double *c = coeff_fourier.getStrip(p);
            double const *r = real.getStrip(p);
            double const *i = imag.getStrip(p);

            for(size_t j=0; j<size_t(num_outputs); j++){
                c[2*j    ] = r[j];
                c[2*j + 1] = i[j];
            }
        }
        writeMatrix(outfilename, num_points, 2 * num_outputs, coeff_fourier.getStrip(0));
        printMatrix(num_points, num_outputs, coeff_fourier.getStrip(0), true);
    }else{
        writeMatrix(outfilename, num_points, num_outputs, coeff);
        printMatrix(num_points, num_outputs, coeff);
    }
}
void TasgridWrapper::loadComputedValues(){
    auto v = verifiedRead(valsfilename, num_outputs);
    if (command == command_loadvalues){
        grid.loadNeededValues(v.release());
    }else{
        if (not grid.isUsingConstruction())
            grid.beginConstruction();
        auto x = verifiedRead(xfilename, num_dimensions);
        grid.loadConstructedPoints(x.release(), v.release());
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

void TasgridWrapper::refineGrid(){
    auto llimits = readLimits();
    TypeCommand effective_command = command;
    if (command == command_refine){
        if (grid.isGlobal() || grid.isSequence()){
            effective_command = command_refine_aniso;
        }else{
            effective_command = command_refine_surp;
        }
    }
    if (effective_command == command_refine_aniso){
        if (min_growth < 1) min_growth = 1;
        if (grid.isGlobal() and ref_output == -1) ref_output = 0;
        grid.setAnisotropicRefinement(depth_type, min_growth, ref_output, llimits);
    }else{ // using surplus refinement
        Data2D<double> scale;
        if (not valsfilename.empty() and (grid.isLocalPolynomial() or grid.isWavelet())){
            scale = readMatrix(valsfilename);
            iassert(scale.getNumStrips() == grid.getNumPoints(), "the number of weights must match the number of points");
            if (ref_output == -1)
                iassert(scale.getStride() == 1, "the number of weights must match the number of outputs");
            if (ref_output > -1)
                iassert(scale.getStride() == (size_t) grid.getNumOutputs(), "there must be one weight per output");
        }
        if (not pass_flag) return;
        if (grid.isGlobal() and ref_output == -1) ref_output = 0;
        if (grid.isGlobal() and grid.isSequence()){
            grid.setSurplusRefinement(tolerance, ref_output, llimits);
        }else{
            grid.setSurplusRefinement(tolerance, tref, ref_output, llimits, scale.release());
        }
    }
}
void TasgridWrapper::getConstructedPoints(){
    if (!grid.isUsingConstruction())
        grid.beginConstruction();

    auto llimits = readLimits();

    std::vector<double> points;
    if (grid.isLocalPolynomial() || grid.isWavelet()){
        Data2D<double> scale;
        if (!valsfilename.empty()){
            scale = readMatrix(valsfilename);
            iassert(scale.getNumStrips() == grid.getNumPoints(), "the number of weights must match the number of points");
            if (ref_output == -1)
                iassert(scale.getStride() == (size_t) grid.getNumOutputs(), "the number of weights must match the number of outputs");
            else
                iassert(scale.getStride() == 1, "there must be one weight per output");
            if (not pass_flag) return;
        }
        points = grid.getCandidateConstructionPoints(tolerance, tref, ref_output, llimits, scale.release());
    }else{
        if (not anisofilename.empty()){
            auto weights = readAnisotropic();
            points = grid.getCandidateConstructionPoints(depth_type, weights, llimits);
        }else{
            points = grid.getCandidateConstructionPoints(depth_type, ref_output, llimits);
        }
    }
    writeMatrix(outfilename, points.size() / num_dimensions, num_dimensions, points.data());
    printMatrix(points.size() / num_dimensions, num_dimensions, points.data());
}

void TasgridWrapper::getPoly(){
    bool integrate = ((depth_type == type_iptotal) || (depth_type == type_ipcurved) || (depth_type == type_iptensor) || (depth_type == type_iphyperbolic));
    std::vector<int> poly = grid.getGlobalPolynomialSpace(integrate);
    std::vector<double> double_poly(poly.size());
    std::transform(poly.begin(), poly.end(), double_poly.begin(), [](int x)->double{ return static_cast<double>(x); });
    writeMatrix(outfilename, (int) double_poly.size() / num_dimensions, num_dimensions, double_poly.data());
    printMatrix((int) double_poly.size() / num_dimensions, num_dimensions, double_poly.data());
}
void TasgridWrapper::setHierarchy(){
    auto vals = readMatrix(valsfilename);
    iassert(vals.getNumStrips() == grid.getNumPoints(),
            (std::string("grid is awaiting ") + std::to_string(grid.getNumPoints()) + " hierarchical surpluses, but "
             + valsfilename + " specifies " + std::to_string(vals.getNumStrips())).c_str());
    if (grid.isFourier())
        iassert((vals.getStride() == (size_t) (2 * grid.getNumOutputs())),
                (std::string("fourier grid is set for ") + std::to_string(grid.getNumOutputs()) + " outputs, but "
                 + valsfilename + " specifies " + std::to_string(vals.getStride())).c_str());
    else
        iassert((vals.getStride() == (size_t) grid.getNumOutputs()),
                (std::string("grid is set for ") + std::to_string(grid.getNumOutputs()) + " outputs, but "
                 + valsfilename + " specifies " + std::to_string(vals.getStride())).c_str());
    if (not pass_flag) return;
    grid.setHierarchicalCoefficients(vals.release());
}

template<typename iomode>
Data2D<double> readMatrixFromOpen(std::istream &is){
    int rows = IO::readNumber<iomode, int>(is);
    int cols = IO::readNumber<iomode, int>(is);
    return Data2D<double>(cols, rows, IO::readVector<iomode, double>(is, Utils::size_mult(cols, rows)));
}

Data2D<double> TasgridWrapper::readMatrix(std::string const &filename) const{
    Data2D<double> matrix;
    if (filename.empty()) return matrix;
    std::ifstream ifs;
    ifs.open(filename, std::ios::in | std::ios::binary);
    iassert(ifs.good(), (std::string("could not open file ") + filename).c_str());
    if (not pass_flag) return matrix;
    char tsg[3] = {'A', 'A', 'A'};
    ifs.read(tsg, 3*sizeof(char));
    if ((tsg[0] == 'T') && (tsg[1] == 'S') && (tsg[2] == 'G')){
        matrix = readMatrixFromOpen<IO::mode_binary_type>(ifs);
    }else{ // not a binary file
        ifs.close();
        ifs.open(filename);
        matrix = readMatrixFromOpen<IO::mode_ascii_type>(ifs);
    }
    if (matrix.empty()) cerr << "WARNING: empty file " << filename << "\n";
    return matrix;
}
void TasgridWrapper::writeMatrix(std::string const &filename, int rows, int cols, const double mat[]) const{
    if (filename.empty()) return;
    size_t cols_t = (size_t) cols;
    std::ofstream ofs;
    if (useASCII){
        Utils::Wrapper2D<const double> matrix(cols, mat);
        ofs.open(filename);
        ofs << rows << " " << cols << "\n";
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
    cout << rows << " " << cols << "\n";
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
    cout << "\n";
}

#endif
