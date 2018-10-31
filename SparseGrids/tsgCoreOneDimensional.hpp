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
    CustomTabulated();
    ~CustomTabulated();

    // I/O subroutines
    void write(std::ofstream &ofs) const;
    void read(const char* filename);
    void read(std::ifstream &ifs);
    void writeBinary(std::ofstream &ofs) const;
    void readBinary(std::ifstream &ifs);

    int getNumLevels() const;
    int getNumPoints(int level) const;
    int getIExact(int level) const;
    int getQExact(int level) const;

    void getWeightsNodes(int level, std::vector<double> &w, std::vector<double> &x) const;
    const char* getDescription() const;

protected:
    void reset();

private:
    int num_levels;
    std::vector<int> num_nodes;
    std::vector<int> precision;

    std::vector<std::vector<double>> nodes;
    std::vector<std::vector<double>> weights;
    std::string description;
};

namespace OneDimensionalMeta{
    int getNumPoints(int level, TypeOneDRule rule);
    int getIExact(int level, TypeOneDRule rule);
    int getQExact(int level, TypeOneDRule rule);

    bool isNonNested(TypeOneDRule rule);
    bool isSequence(TypeOneDRule rule);
    bool isGlobal(TypeOneDRule rule);
    bool isSingleNodeGrowth(TypeOneDRule rule);
    bool isLocalPolynomial(TypeOneDRule rule);
    bool isWavelet(TypeOneDRule rule);
    bool isFourier(TypeOneDRule rule);

    TypeOneDRule getIORuleString(const char *name);
    const char* getIORuleString(TypeOneDRule rule);
    const char* getHumanString(TypeOneDRule rule);
    TypeOneDRule getIORuleInt(int index);
    int getIORuleInt(TypeOneDRule rule);

    bool isTypeCurved(TypeDepth type);

    TypeDepth getIOTypeString(const char *name);
    TypeDepth getIOTypeInt(int type);

    TypeRefinement getIOTypeRefinementString(const char *name);
    TypeRefinement getIOTypeRefinementInt(int ref);
}

namespace OneDimensionalNodes{
    // non-nested rules
    void getGaussLegendre(int m, std::vector<double> &w, std::vector<double> &x);
    void getChebyshev(int m, std::vector<double> &w, std::vector<double> &x);
    void getGaussChebyshev1(int m, std::vector<double> &w, std::vector<double> &x);
    void getGaussChebyshev2(int m, std::vector<double> &w, std::vector<double> &x);
    void getGaussJacobi(int m, std::vector<double> &w, std::vector<double> &x, double alpha, double beta);
    void getGaussHermite(int m, std::vector<double> &w, std::vector<double> &x, double alpha);
    void getGaussLaguerre(int m, std::vector<double> &w, std::vector<double> &x, double alpha);

    // nested rules
    void getClenshawCurtisNodes(int level, std::vector<double> &nodes);
    double getClenshawCurtisWeight(int level, int point);

    void getClenshawCurtisNodesZero(int level, std::vector<double> &nodes); // assuming zero boundary
    double getClenshawCurtisWeightZero(int level, int point); // assuming zero boundary

    void getFejer2Nodes(int level, std::vector<double> &nodes);
    double getFejer2Weight(int level, int point);

    void getRLeja(int n, std::vector<double> &nodes);
    void getRLejaCentered(int n, std::vector<double> &nodes);
    void getRLejaShifted(int n, std::vector<double> &nodes);

    void getFourierNodes(int level, std::vector<double> &nodes);
}

}

#endif
