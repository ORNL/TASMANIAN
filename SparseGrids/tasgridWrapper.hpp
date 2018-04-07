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

#ifndef __TASGRID_WRAPPER_HPP
#define __TASGRID_WRAPPER_HPP

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <string.h>
#include <math.h>

#include "TasmanianSparseGrid.hpp"

using std::cout;
using std::endl;
using std::setw;

using namespace TasGrid;

enum TypeCommand{
    command_none,
    command_makeglobal,
    command_makesequence,
    command_makelocalp,
    command_makewavelet,
    command_makequadrature,
    command_update,

    command_setconformal,

    command_getquadrature,
    command_getinterweights,
    command_getpoints,
    command_getneeded,
    command_getrefcoeff,

    command_loadvalues,

    command_evaluate,
    command_integrate,

    command_getanisocoeff,
    command_refine_surp,
    command_refine_aniso,
    command_refine,
    command_refine_clear,
    command_refine_merge,

    command_getpoly,

    command_summary,

    command_getcoefficients, // ML section
    command_setcoefficients, // ML section
    command_evalhierarchical_sparse, // ML section
    command_evalhierarchical_dense, // ML section

    command_getpointsindex,
    command_getneededindex
};

enum TypeConformalMap{
    conformal_none,
    conformal_asin
};

class TasgridWrapper{
public:
    TasgridWrapper();
    ~TasgridWrapper();

    static TypeConformalMap getConfromalType(const char* name);

    void setCommand(TypeCommand com);
    TypeCommand getCommand() const;

    void setNumDimensions(int dim);
    void setNumOutputs(int out);
    void setNumDepth(int d);
    void setOrder(int o);
    void setDepthType(TypeDepth dt);
    void setConformalType(TypeConformalMap con);
    void setRule(TypeOneDRule r);
    TypeOneDRule getRule() const;
    void setAlpha(double a);
    void setBeta(double b);

    void setTolerance(double tol);
    void setRefOutput(int out);
    void setMinGrowth(int mg);
    void setTypeRefinement(TypeRefinement rt);

    void setGridFilename(const char *filename);
    void setOutFilename(const char *filename);
    void setValsFilename(const char *filename);
    void setXFilename(const char *filename);
    void setAnisoFilename(const char *filename);
    void setTransformFilename(const char *filename);
    void setConformalFilename(const char *filename);
    void setCustomFilename(const char *filename);
    void setLevelLimitsFilename(const char *filename);

    void setPrintPoints(bool pp);
    void setUseASCII(bool ascii);
    void setGPID(int gpuid);

    bool executeCommand();

    static bool isCreateCommand(TypeCommand com);

protected:
    bool checkSane() const;

    void createGlobalGird();
    void createSequenceGird();
    void createLocalPolynomialGird();
    void createWaveletGird();
    void createQuadrature();
    bool updateGrid();
    void writeGrid() const;
    bool readGrid();

    void outputPoints(bool useNeeded) const;
    void outputQuadrature() const;
    void outputHierarchicalCoefficients() const;

    bool setConformalTransformation();

    bool loadValues();

    bool getInterWeights();
    bool getEvaluate();
    bool getIntegrate();
    bool getAnisoCoeff();

    bool refineGrid();
    bool cancelRefine();
    bool mergeRefine();

    bool setHierarchy();
    bool getEvalHierarchyDense();
    bool getEvalHierarchySparse();

    bool getPoly();

    bool getSummary();

    bool getSurpluses();
    bool getPointsIndexes();
    bool getNeededIndexes();

    int* readAnisotropicFile(int num_weights) const;
    double* readTransform() const;
    int* readLevelLimits(int num_weights) const;

    static void readMatrix(const char *filename, size_t &rows, size_t &cols, double* &mat);
    static void writeMatrix(const char *filename, int rows, int cols, const double mat[], bool ascii);
    static void printMatrix(int rows, int cols, const double mat[]);

private:
    TasmanianSparseGrid *grid;

    TypeCommand command;

    int num_dimensions, num_outputs, depth, order;
    TypeDepth depth_type;
    TypeOneDRule rule;
    TypeConformalMap conformal;

    double alpha, beta;
    bool set_alpha, set_beta;

    double tolerance;
    bool set_tolerance;

    int ref_output, min_growth;
    TypeRefinement tref;
    bool set_tref;

    const char *gridfilename;
    const char *outfilename;
    const char *valsfilename;
    const char *xfilename;
    const char *anisofilename;
    const char *transformfilename;
    const char *conformalfilename;
    const char *customfilename;
    const char *levellimitfilename;

    bool printCout;
    bool useASCII;
    int set_gpuid;
};

#endif
