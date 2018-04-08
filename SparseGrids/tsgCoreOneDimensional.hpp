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

#ifndef __TSG_CORE_ONE_DIMENSIONAL_HPP
#define __TSG_CORE_ONE_DIMENSIONAL_HPP

#include "tsgEnumerates.hpp"
#include "tsgLinearSolvers.hpp"

#include <string>
#include <math.h>

namespace TasGrid{

class CustomTabulated{
public:
    CustomTabulated(std::ostream *os = 0);
    CustomTabulated(const char* filename, std::ostream *os = 0);
    CustomTabulated(const CustomTabulated &custom, std::ostream *os = 0);
    ~CustomTabulated();

    // I/O subroutines
    void write(std::ofstream &ofs) const;
    bool read(std::ifstream &ifs);
    void writeBinary(std::ofstream &ofs) const;
    bool readBinary(std::ifstream &ifs);
    void copyRule(const CustomTabulated *custom);

    int getNumLevels() const;
    int getNumPoints(int level) const;
    int getIExact(int level) const;
    int getQExact(int level) const;

    void getWeightsNodes(int level, double* &w, double* &x) const;
    const char* getDescription() const;

protected:
    void reset();

private:
    int num_levels;
    int *num_nodes;
    int *precision;
    int *offsets;
    double *nodes;
    double *weights;
    std::string *description;
    std::ostream *logstream;
};

class OneDimensionalMeta{
public:
    OneDimensionalMeta();
    OneDimensionalMeta(const CustomTabulated *ccustom);
    ~OneDimensionalMeta();

    int getNumPoints(int level, TypeOneDRule rule) const;
    int getIExact(int level, TypeOneDRule rule) const;
    int getQExact(int level, TypeOneDRule rule) const;

    const CustomTabulated *getCustom() const;

    static bool isNonNested(TypeOneDRule rule);
    static bool isSequence(TypeOneDRule rule);
    static bool isGlobal(TypeOneDRule rule);
    static bool isSingleNodeGrowth(TypeOneDRule rule);
    static bool isLocalPolynomial(TypeOneDRule rule);
    static bool isWavelet(TypeOneDRule rule);

    static TypeOneDRule getIORuleString(const char *name);
    static const char* getIORuleString(TypeOneDRule rule);
    static const char* getHumanString(TypeOneDRule rule);
    static TypeOneDRule getIORuleInt(int index);
    static int getIORuleInt(TypeOneDRule rule);

    static TypeDepth getIOTypeString(const char *name);
    static TypeDepth getIOTypeInt(int type);

    static TypeRefinement getIOTypeRefinementString(const char *name);
    static TypeRefinement getIOTypeRefinementInt(int ref);

private:
    // add the custom class here, but alias to something created by GlobalGrid
    const CustomTabulated *custom;
};

namespace OneDimensionalNodes{
    // non-nested rules
    void getGaussLegendre(int m, double* &w, double* &x);
    void getChebyshev(int m, double* &w, double* &x);
    void getGaussChebyshev1(int m, double* &w, double* &x);
    void getGaussChebyshev2(int m, double* &w, double* &x);
    void getGaussJacobi(int m, double* &w, double* &x, double alpha, double beta);
    void getGaussHermite(int m, double* &w, double* &x, double alpha);
    void getGaussLaguerre(int m, double* &w, double* &x, double alpha);

    // nested rules
    double* getClenshawCurtisNodes(int level);
    double getClenshawCurtisWeight(int level, int point);

    double* getClenshawCurtisNodesZero(int level); // assuming zero boundary
    double getClenshawCurtisWeightZero(int level, int point); // assuming zero boundary

    double* getFejer2Nodes(int level);
    double getFejer2Weight(int level, int point);

    double* getRLeja(int n);
    double* getRLejaCentered(int n);
    double* getRLejaShifted(int n);
}

}

#endif
