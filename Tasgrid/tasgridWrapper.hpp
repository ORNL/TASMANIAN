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

#include "gridtestCLICommon.hpp"

enum TypeCommand{
    command_none,
    command_makeglobal,
    command_makesequence,
    command_makelocalp,
    command_makewavelet,
    command_makefourier,
    command_makequadrature,
    command_makeexoquad,
    command_update,

    command_setconformal,

    command_getquadrature,
    command_getinterweights,
    command_getdiffweights,

    command_getpoints,
    command_getneeded,

    command_loadvalues,

    command_evaluate,
    command_integrate,
    command_differentiate,

    command_getanisocoeff,
    command_refine_surp,
    command_refine_aniso,
    command_refine,
    command_refine_clear,
    command_refine_merge,

    command_using_construct,
    command_get_candidate_construction,
    command_load_construction,

    command_getpoly,

    command_summary,

    command_getcoefficients, // ML section
    command_setcoefficients, // ML section
    command_evalhierarchical_sparse, // ML section
    command_evalhierarchical_dense, // ML section
    command_gethsupport, // ML section

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

    static TypeCommand hasCommand(std::string const &s);
    static TypeConformalMap getConfromalType(const char* name);

    void setCommand(TypeCommand com){ command = com; }
    TypeCommand getCommand() const{ return command; }
    void setNumDimensions(int dim){ num_dimensions = dim; }
    void setNumOutputs(int out){ num_outputs = out; }
    void setNumDepth(int d){ depth = d; }
    void setOrder(int o){ order = o; }
    void setDepthType(TypeDepth dt){ depth_type = dt; }
    void setConformalType(TypeConformalMap con){ conformal = con; }
    void setRule(TypeOneDRule r){ rule = r;}
    TypeOneDRule getRule() const{ return rule; }
    void setAlpha(double a){ alpha = a;  set_alpha = true; }
    void setBeta(double b){ beta = b; set_beta = true; }
    void setTolerance(double tol){ tolerance = tol;  set_tolerance = true; }
    void setRefOutput(int out){ ref_output = out; }
    void setMinGrowth(int mg){ min_growth = mg; }
    void setTypeRefinement(TypeRefinement rt){ tref = rt; set_tref = true; }
    void setPrintPoints(bool pp){ printCout = pp; }
    void setUseASCII(bool ascii){ useASCII = ascii; }
    void setGPID(int gpuid){ set_gpuid = gpuid; }

    void setGridFilename(std::string const &filename){ gridfilename = filename; }
    void setOutFilename(std::string const &filename){ outfilename = filename; }
    void setValsFilename(std::string const &filename){ valsfilename = filename; }
    void setXFilename(std::string const &filename){ xfilename = filename; }
    void setAnisoFilename(std::string const &filename){ anisofilename = filename; }
    void setTransformFilename(std::string const &filename){ transformfilename = filename; }
    void setConformalFilename(std::string const &filename){ conformalfilename = filename; }
    void setCustomFilename(std::string const &filename){ customfilename = filename; }
    void setLevelLimitsFilename(std::string const &filename){ levellimitfilename = filename; }

    void setShift(double s){ shift = s; set_shift = true; }
    void setWeightFilename(std::string const &filename){ weightfilename = filename; }
    void setDescription(std::string const &desc){ description = desc; }
    void setIsSymmetric(bool b) { is_symmetric_weight_function = b; }

    bool executeCommand();

    static bool isCreateCommand(TypeCommand com);

protected:
    enum class output_points_mode{
        regular, needed
    };

    bool checkSane() const;
    bool checkSanePostRead() const;

    void iassert(bool, const char*) const;
    bool readGridfile();
    std::vector<int> readLimits() const;
    std::vector<int> readAnisotropic() const;

    void createExoticQuadrature();

    void setTransform();
    void setConformal();

    void outputPoints(output_points_mode mode) const;
    void outputIndexes(output_points_mode mode) const;
    void outputQuadrature() const;

    void writeGrid() const;

    void outputPoints(bool useNeeded) const;

    void outputHierarchicalCoefficients() const;

    Data2D<double> verifiedRead(std::string const& filename, int expected_stride) const;
    void processEvalLike() const;
    void processOutputLike() const;

    void loadComputedValues();
    void setHierarchy();
    void refineGrid();

    void getPoly();

    void getConstructedPoints();

    Data2D<double> readMatrix(std::string const &filename) const;
    void writeMatrix(std::string const &filename, int rows, int cols, const double mat[]) const;
    void printMatrix(int rows, int cols, const double mat[], bool isComplex = false) const;

    // Overloads.
    template<typename T>
    inline void writeMatrix(std::string const &filename, const Data2D<T> &mat) const { writeMatrix(filename, mat.getNumStrips(), mat.getStride(), mat.data()); }
    template<typename T>
    inline void printMatrix(const Data2D<T> &mat, bool isComplex = false) const { printMatrix(mat.getNumStrips(), mat.getStride(), mat.data(), isComplex); }

  private:
    TasmanianSparseGrid grid;
    CustomTabulated ct;

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

    std::string gridfilename;
    std::string outfilename;
    std::string valsfilename;
    std::string xfilename;
    std::string anisofilename;
    std::string transformfilename;
    std::string conformalfilename;
    std::string customfilename;
    std::string levellimitfilename;

    bool printCout;
    bool useASCII;
    int set_gpuid;

    double shift;
    bool set_shift;
    std::string weightfilename;
    std::string description;
    bool is_symmetric_weight_function;

    mutable bool pass_flag; // report internal errors
};

#endif
