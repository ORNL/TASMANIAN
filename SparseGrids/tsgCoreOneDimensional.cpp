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

#ifndef __TSG_CORE_ONE_DIMENSIONAL_CPP
#define __TSG_CORE_ONE_DIMENSIONAL_CPP

#include "tsgCoreOneDimensional.hpp"
#include "tsgIOHelpers.hpp"

namespace TasGrid{

template<bool iomode> void CustomTabulated::write(std::ostream &ofs) const{
    if (iomode == mode_ascii){ // format is messy here
        ofs << "description: " << description.c_str() << std::endl;
        ofs << "levels: " << num_levels << std::endl;
        for(int i=0; i<num_levels; i++){
            ofs << num_nodes[i] << " " << precision[i] << std::endl;
        }
        ofs << std::scientific; ofs.precision(17);
        for(int l=0; l<num_levels; l++){
            auto x = nodes[l].begin();
            for(auto w : weights[l]) ofs << w << " " << *x++ << std::endl;
        }
    }else{
        int num_description = (int) description.size();
        ofs.write((char*) &num_description, sizeof(int));
        ofs.write(description.c_str(), num_description * sizeof(char));
        ofs.write((char*) &num_levels, sizeof(int));
        ofs.write((char*) num_nodes.data(), num_levels * sizeof(int));
        ofs.write((char*) precision.data(), num_levels * sizeof(int));
        for(int l=0; l<num_levels; l++){
            ofs.write((char*) weights[l].data(), weights[l].size() * sizeof(double));
            ofs.write((char*) nodes[l].data(), nodes[l].size() * sizeof(double));
        }
    }
}

template<bool iomode> void CustomTabulated::read(std::istream &is){
    if (iomode == mode_ascii){ // messy format chosen for better human readability, hard to template and generalize
        std::string T;
        char dummy;
        is >> T;
        if (!(T.compare("description:") == 0)){ throw std::invalid_argument("ERROR: wrong file format of custom tables on line 1"); }
        is.get(dummy);
        description = std::string();
        getline(is, description);

        is >> T;
        if (!(T.compare("levels:") == 0)){ throw std::invalid_argument("ERROR: wrong file format of custom tables on line 2"); }
        is >> num_levels;

        num_nodes.resize(num_levels);
        precision.resize(num_levels);
        for(int i=0; i<num_levels; i++)
            is >> num_nodes[i] >> precision[i];

        nodes.resize(num_levels);
        weights.resize(num_levels);
        for(int l=0; l<num_levels; l++){
            nodes[l].resize(num_nodes[l]);
            weights[l].resize(num_nodes[l]);
            auto x = nodes[l].begin();
            for(auto &w : weights[l]) is >> w >> *x++;
        }
    }else{
        int num_description = 0;
        is.read((char*) &num_description, sizeof(int));
        std::vector<char> desc((size_t) num_description+1);
        is.read(desc.data(), num_description);
        desc[num_description] = '\0';
        description = desc.data();

        is.read((char*) &num_levels, sizeof(int));
        num_nodes.resize(num_levels);
        precision.resize(num_levels);
        is.read((char*) num_nodes.data(), num_levels * sizeof(int));
        is.read((char*) precision.data(), num_levels * sizeof(int));

        nodes.resize(num_levels);
        weights.resize(num_levels);
        for(int l=0; l<num_levels; l++){
            nodes[l].resize(num_nodes[l]);
            weights[l].resize(num_nodes[l]);
            is.read((char*) weights[l].data(), num_nodes[l] * sizeof(double));
            is.read((char*) nodes[l].data(), num_nodes[l] * sizeof(double));
        }
    }
}

template void CustomTabulated::write<mode_ascii>(std::ostream &) const;
template void CustomTabulated::write<mode_binary>(std::ostream &) const;
template void CustomTabulated::read<mode_ascii>(std::istream &is);
template void CustomTabulated::read<mode_binary>(std::istream &is);

void CustomTabulated::read(const char* filename){
    std::ifstream ifs; ifs.open(filename);
    if (!ifs){
        std::string message = "Could not open the custom rule file: ";
        message += filename;
        throw std::invalid_argument(message);
    }
    read<mode_ascii>(ifs);
    ifs.close();
}

void CustomTabulated::getWeightsNodes(int level, std::vector<double> &w, std::vector<double> &x) const{
    w = weights[level];
    x = nodes[level];
}
void CustomTabulated::getWeightsNodes(int level, double w[], double x[]) const{
    std::copy(weights[level].begin(), weights[level].end(), w);
    std::copy(nodes[level].begin(), nodes[level].end(), x);
}

const char* CustomTabulated::getDescription() const{ return description.c_str(); }

int OneDimensionalMeta::getNumPoints(int level, TypeOneDRule rule){
    int lcc;
    switch (rule){
        case rule_chebyshev:
        case rule_gausslegendre:
        case rule_leja:
        case rule_rleja:
        case rule_rlejashifted:
        case rule_maxlebesgue:
        case rule_minlebesgue:
        case rule_mindelta:
        case rule_gausschebyshev1:
        case rule_gausschebyshev2:
        case rule_gaussgegenbauer:
        case rule_gaussjacobi:
        case rule_gausslaguerre:
        case rule_gausshermite:       return level+1;

        case rule_chebyshevodd:
        case rule_gausslegendreodd:
        case rule_lejaodd:
        case rule_rlejaodd:
        case rule_maxlebesgueodd:
        case rule_minlebesgueodd:
        case rule_mindeltaodd:
        case rule_gausschebyshev1odd:
        case rule_gausschebyshev2odd:
        case rule_gaussgegenbauerodd:
        case rule_gaussjacobiodd:
        case rule_gausslaguerreodd:
        case rule_gausshermiteodd:    return 2*level+1;

        case rule_rlejashiftedeven:   return 2*(level + 1);

        case rule_rlejashifteddouble: return Maths::pow2(level+1);

        case rule_clenshawcurtis0:
        case rule_gausspatterson:     return Maths::pow2(level+1) - 1;

        case rule_clenshawcurtis:     return (level == 0) ? 1 : (Maths::pow2(level) + 1);
        case rule_rlejadouble2:       if (level < 3){ return getNumPoints(level, rule_clenshawcurtis); }; lcc = 2 + (level-3)/2; return (getNumPoints(lcc,rule_clenshawcurtis) + ((getNumPoints(lcc+1,rule_clenshawcurtis)-getNumPoints(lcc,rule_clenshawcurtis))/2) * ((level-3)%2 +1));
        case rule_rlejadouble4:       if (level < 3){ return getNumPoints(level, rule_clenshawcurtis); }; lcc = 2 + (level-3)/4; return (getNumPoints(lcc,rule_clenshawcurtis) + ((getNumPoints(lcc+1,rule_clenshawcurtis)-getNumPoints(lcc,rule_clenshawcurtis))/4) * ((level-3)%4 +1));
        case rule_fejer2:             return Maths::pow2(level+1) - 1;

        case rule_fourier:            return Maths::pow3(level);

        default:
            return level; // should not be called, but compiler complains for the lack of return/default
    }
}
int OneDimensionalMeta::getIExact(int level, TypeOneDRule rule){
    switch (rule){
        case rule_chebyshev:
        case rule_gausslegendre:
        case rule_leja:
        case rule_rleja:
        case rule_rlejashifted:
        case rule_maxlebesgue:
        case rule_minlebesgue:
        case rule_mindelta:
        case rule_gausschebyshev1:
        case rule_gausschebyshev2:
        case rule_gaussgegenbauer:
        case rule_gaussjacobi:
        case rule_gausslaguerre:
        case rule_gausshermite:       return level;

        case rule_chebyshevodd:
        case rule_gausslegendreodd:
        case rule_lejaodd:
        case rule_rlejaodd:
        case rule_maxlebesgueodd:
        case rule_minlebesgueodd:
        case rule_mindeltaodd:
        case rule_gausschebyshev1odd:
        case rule_gausschebyshev2odd:
        case rule_gaussgegenbauerodd:
        case rule_gaussjacobiodd:
        case rule_gausslaguerreodd:
        case rule_gausshermiteodd:    return 2*level;

        case rule_rlejashiftedeven:   return 2*level + 1;

        case rule_rlejashifteddouble: return Maths::pow2(level+1) -1;

        case rule_gausspatterson:
        case rule_fejer2:             return Maths::pow2(level+1) -2;

        case rule_clenshawcurtis:     return (level > 0) ? Maths::pow2(level) : 0;
        case rule_clenshawcurtis0:    return Maths::pow2(level+1) +1;
        case rule_rlejadouble2:       return getNumPoints(level,rule_rlejadouble2)-1;
        case rule_rlejadouble4:       return getNumPoints(level,rule_rlejadouble4)-1;
        case rule_fourier:            return (Maths::pow3(level)-1)/2;
        default:
            return level; // should not be called, but compiler complains for the lack of return/default
    }
}
int OneDimensionalMeta::getQExact(int level, TypeOneDRule rule){
    switch (rule){
        case rule_chebyshevodd:
        case rule_gausslegendre:
        case rule_gausschebyshev1:
        case rule_gausschebyshev2:
        case rule_gaussgegenbauer:
        case rule_gaussjacobi:
        case rule_gausslaguerre:
        case rule_gausshermite:       return 2*level+1;

        case rule_gausschebyshev1odd:
        case rule_gausschebyshev2odd:
        case rule_gaussgegenbauerodd:
        case rule_gausshermiteodd:
        case rule_gaussjacobiodd:
        case rule_gausslaguerreodd:
        case rule_gausslegendreodd:   return 4*level+1;

        case rule_rleja:
        case rule_rlejashifted:       return level;
        case rule_leja:
        case rule_maxlebesgue:
        case rule_minlebesgue:
        case rule_mindelta:           return ((level == 0) ? 1 : level) + ((level == 2) ? 1 : 0);

        case rule_lejaodd:
        case rule_maxlebesgueodd:
        case rule_minlebesgueodd:
        case rule_mindeltaodd:        return ((level == 0) ? 1 : 2*level) + ((level == 1) ? 1 : 0);
        case rule_rlejaodd:           return 2*level;

        case rule_rlejashiftedeven:   return 2*level;

        case rule_rlejashifteddouble: return Maths::pow2(level+1) -1;

        case rule_gausspatterson:     return (level == 0) ? 1 : (3 * Maths::pow2(level) - 1);
        case rule_clenshawcurtis:     return (level == 0) ? 1 : (Maths::pow2(level) + 1);
        case rule_clenshawcurtis0:    return (level == 0) ? 1 : (Maths::pow2(level+1) + 1);
        case rule_chebyshev:          return level+1;
        case rule_rlejadouble2:       return getNumPoints(level,rule_rlejadouble2);
        case rule_rlejadouble4:       return getNumPoints(level,rule_rlejadouble4)-1;
        case rule_fejer2:             return Maths::pow2(level+1) - 1;
        case rule_fourier:            return (Maths::pow3(level)-1)/2;
        default:
            return level; // should not be called, but compiler complains for the lack of return/default
    }
}

bool OneDimensionalMeta::isNonNested(TypeOneDRule rule){
    return  ((rule == rule_chebyshev) || (rule == rule_chebyshevodd) || (rule == rule_gausslegendre) || (rule == rule_gausslegendreodd) || (rule == rule_gausschebyshev1)
           || (rule == rule_gausschebyshev1odd) || (rule == rule_gausschebyshev2) || (rule == rule_gausschebyshev2odd) || (rule == rule_gaussgegenbauer)
           || (rule == rule_gaussgegenbauerodd) || (rule == rule_gaussjacobi) || (rule == rule_gaussjacobiodd) || (rule == rule_gausslaguerre) || (rule == rule_gausslaguerreodd)
           || (rule == rule_gausshermite) || (rule == rule_gausshermiteodd) || (rule == rule_customtabulated));
}

bool OneDimensionalMeta::isSequence(TypeOneDRule rule){
    return ((rule == rule_leja) || (rule == rule_rleja) || (rule == rule_rlejashifted) || (rule == rule_maxlebesgue) || (rule == rule_minlebesgue) || (rule == rule_mindelta));
}
bool OneDimensionalMeta::isGlobal(TypeOneDRule rule){
    return !((rule == rule_semilocalp) || (rule == rule_localp0) || (rule == rule_localp) || (rule == rule_localpb)
          || (rule == rule_wavelet) || (rule == rule_fourier) || (rule == rule_none));
}
bool OneDimensionalMeta::isSingleNodeGrowth(TypeOneDRule rule){
    return ((rule == rule_leja) || (rule == rule_rleja) || (rule == rule_rlejashifted) || (rule == rule_maxlebesgue) || (rule == rule_minlebesgue) || (rule == rule_mindelta) ||
         (rule == rule_gausslegendre) || (rule == rule_gausschebyshev1) || (rule == rule_gausschebyshev2) || (rule == rule_gaussgegenbauer) ||
         (rule == rule_gaussjacobi) || (rule == rule_gausslaguerre) || (rule == rule_gausshermite));
}
bool OneDimensionalMeta::isLocalPolynomial(TypeOneDRule rule){
    return ((rule == rule_localp) || (rule == rule_localp0) || (rule == rule_semilocalp) || (rule == rule_localpb));
}
bool OneDimensionalMeta::isWavelet(TypeOneDRule rule){ // used by the tasgrid wrapper
    return (rule == rule_wavelet);
}
bool OneDimensionalMeta::isFourier(TypeOneDRule rule){ // used by the tasgrid wrapper
    return (rule == rule_fourier);
}

const char* OneDimensionalMeta::getHumanString(TypeOneDRule rule){
    switch (rule){
        case rule_clenshawcurtis:     return "Clenshaw-Curtis";
        case rule_clenshawcurtis0:    return "Clenshaw-Curtis zero boundary conditions";
        case rule_chebyshev:          return "Chebyshev";
        case rule_chebyshevodd:       return "Chebyshev, odd rules only";
        case rule_gausslegendre:      return "Gauss-Legendre";
        case rule_gausslegendreodd:   return "Gauss-Legendre, odd rules only";
        case rule_gausspatterson:     return "Gauss-Patterson";
        case rule_leja:               return "Leja";
        case rule_lejaodd:            return "Leja, odd rules only";
        case rule_rleja:              return "R-Leja";
        case rule_rlejaodd:           return "R-Leja odd rules only";
        case rule_rlejadouble2:       return "R-Leja doubling every 2 levels";
        case rule_rlejadouble4:       return "R-Leja doubling every 4 levels";
        case rule_rlejashifted:       return "R-Leja shifted by 1/2";
        case rule_rlejashiftedeven:   return "R-Leja shifted by 1/2, even rules only";
        case rule_rlejashifteddouble: return "R-Leja shifted by 1/2, doubling rule";
        case rule_maxlebesgue:        return "Maximum of the Lebesgue function";
        case rule_maxlebesgueodd:     return "Maximum of the Lebesgue function, odd rules only";
        case rule_minlebesgue:        return "Minimum (greedy) of the Lebesgue constant";
        case rule_minlebesgueodd:     return "Minimum (greedy) of the Lebesgue constant, odd rules only";
        case rule_mindelta:           return "Minimum surplus operator norm";
        case rule_mindeltaodd:        return "Minimum surplus operator norm, odd rules only";
        case rule_gausschebyshev1:    return "Gauss-Chebyshev of type 1";
        case rule_gausschebyshev1odd: return "Gauss-Chebyshev of type 1, odd rules only";
        case rule_gausschebyshev2:    return "Gauss-Chebyshev of type 2";
        case rule_gausschebyshev2odd: return "Gauss-Chebyshev of type 2, odd rules only";
        case rule_fejer2:             return "Fejer type 2";
        case rule_gaussgegenbauer:    return "Gauss-Gegenbauer";
        case rule_gaussgegenbauerodd: return "Gauss-Gegenbauer, odd rules only";
        case rule_gaussjacobi:        return "Gauss-Jacobi";
        case rule_gaussjacobiodd:     return "Gauss-Jacobi, odd rules only";
        case rule_gausslaguerre:      return "Gauss-Laguerre";
        case rule_gausslaguerreodd:   return "Gauss-Laguerre, odd rules only";
        case rule_gausshermite:       return "Gauss-Hermite";
        case rule_gausshermiteodd:    return "Gauss-Hermite, odd rules only";
        case rule_customtabulated:    return "Custom rule";
        case rule_localp:             return "Local polynomials";
        case rule_localp0:            return "Local polynomials zero boundary conditions";
        case rule_localpb:            return "Local polynomials focused nodes on the boundary";
        case rule_semilocalp:         return "Semi-Local polynomials";
        case rule_wavelet:            return "Wavelets";
        case rule_fourier:            return "Fourier / trigonometric";
        default:
            return "unknown";
    }
}

// Gauss-Legendre
void OneDimensionalNodes::getGaussLegendre(int m, std::vector<double> &w, std::vector<double> &x){
    double mu0 = 2.0;
    std::vector<double> diag(m, 0.0), off_diag(m-1);
    for(int i=0; i<m-1; i++){
        off_diag[i] = std::sqrt((double) ((i+1)*(i+1)) / ((double) (4*(i+1)*(i+1) - 1)));
    }
    TasmanianTridiagonalSolver::decompose(diag, off_diag, mu0, x, w);
}

// Chebyshev
void OneDimensionalNodes::getChebyshev(int m, std::vector<double> &w, std::vector<double> &x){
    // get Clanshaw-Curtis quadrature points
    w.resize(m);
    x.resize(m);
    int i, j;
    double b;
    if (m == 1){
        w[0] = 2.0; x[0] = 0.0;
        return;
    };

    for(i=0; i<m; i++){
        x[i] = std::cos(((double) (m-i-1)) * Maths::pi / ((double) (m-1)));
    };

    // may also have to set the mid-point to 0.0
    x[0] = -1.0; x[m-1] = 1.0;

    for(i=0; i<m; i++){
        w[i] = 1.0;
        double theta = ((double) i) * Maths::pi / ((double) (m-1));
        for(j=1; j<=(m-1)/2; j++){
            if (2*j == (m-1)){
                b = 1.0;
            }else{
                b = 2.0;
            };
            w[i] = w[i] - b * std::cos(2.0 * j * theta) / ((double) (4*j*j - 1));
        };
    };

    w[0] = w[0] / ((double) (m-1));
    for(i=1; i<m-1; i++){
        w[i] = 2.0 * w[i] / ((double) (m-1));
    };
    w[m-1] = w[m-1] / ((double) (m-1));
}

// get Gauss-Chebyshev type 1 quadrature points
void OneDimensionalNodes::getGaussChebyshev1(int m, std::vector<double> &w, std::vector<double> &x){
    w.resize(m);
    x.resize(m);

    for(int i=0; i<m; i++){
        x[m-i-1] = std::cos(Maths::pi*(2*i+1) / (2*((double)m)));
        w[i] = Maths::pi / m;
    }
}
// get Gauss-Chebyshev-type2 quadrature points
void OneDimensionalNodes::getGaussChebyshev2(int m, std::vector<double> &w, std::vector<double> &x){
    w.resize(m);
    x.resize(m);

    for(int i=0; i<m; i++){
        double theta = Maths::pi*((double)(i+1))/((double)(m+1));
        x[m-i-1] = std::cos(theta);
        w[i] = (Maths::pi / ((double)(m+1))) * std::sin(theta)* std::sin(theta);
    }
}
// get Gauss-Jacobi quadrature points
void OneDimensionalNodes::getGaussJacobi(int m, std::vector<double> &w, std::vector<double> &x, double alpha, double beta){
    double ab = alpha + beta;
    double mu0 = pow(2.0, 1.0 + ab) * tgamma(alpha + 1.0) * tgamma(beta + 1.0) / tgamma(2.0 + ab);
    std::vector<double> diag(m), off_diag(m-1);
    diag[0] = (beta - alpha) / (2.0 + ab);
    if (m > 1) {
        off_diag[0] =  std::sqrt(4.0 * (1.0 + alpha) * (1.0 + beta) / ((3.0 + ab) * (2.0 + ab) * (2.0 + ab)));
        for(int i=1; i<m-1; i++){
            double di = (double) (i+1);
            diag[i] =  (beta * beta - alpha * alpha) / ((2.0 * di + ab -2.0) * (2.0 * di + ab));
            off_diag[i] = std::sqrt(4.0 * di * (di + alpha) * (di + beta) * (di + ab) /
                                    (((2.0 * di + ab) * (2.0 * di + ab) - 1.0) * (2.0 * di + ab) * (2.0 * di + ab)));
        }
        double dm = (double) m;
        diag[m-1] = (beta * beta - alpha * alpha) / ((2.0 * dm + ab - 2.0) * (2.0 * dm + ab));;
    }
    TasmanianTridiagonalSolver::decompose(diag, off_diag, mu0, x, w);
}
// get Gauss-Hermite quadrature points
void OneDimensionalNodes::getGaussHermite(int m, std::vector<double> &w, std::vector<double> &x, double alpha){
    double mu0 = tgamma(0.5 * (alpha + 1.0));
    std::vector<double> diag(m, 0.0), off_diag(m-1);
    for(int i=0; i<m-1; i++){
        double di = (double) (i+1);
        off_diag[i] = std::sqrt(0.5 * (di + alpha * ((double) ((i + 1) % 2))));
    }
    TasmanianTridiagonalSolver::decompose(diag, off_diag, mu0, x, w);
}
// get Gauss-Laguerre quadrature points
void OneDimensionalNodes::getGaussLaguerre(int m, std::vector<double> &w, std::vector<double> &x, double alpha){
    double mu0 = tgamma(alpha + 1.0);
    std::vector<double> diag(m), off_diag(m-1);
    for(int i=0; i<m-1; i++){
        double di = (double) (i+1);
        diag[i] = 2.0 * di - 1.0 + alpha;
        off_diag[i] = std::sqrt(di * (di + alpha));
    }
    double dm = (double) m;
    diag[m-1] = 2.0 * dm - 1.0 + alpha;
    TasmanianTridiagonalSolver::decompose(diag, off_diag, mu0, x, w);
}

// Clenshaw-Curtis
std::vector<double> OneDimensionalNodes::getClenshawCurtisNodes(int level){
    int n = OneDimensionalMeta::getNumPoints(level, rule_clenshawcurtis);
    std::vector<double> nodes(n, 0.0);
    if (level > 0){
        nodes[1] = -1.0; nodes[2] = 1.0;
        int count = 3;
        for(int l=2; l<=level; l++){
            n = OneDimensionalMeta::getNumPoints(l, rule_clenshawcurtis);
            for(int i=1; i<n; i+=2){
                  nodes[count++] = std::cos(Maths::pi * ((double) (n-i-1)) / ((double) (n - 1)));
            }
        }
    }
    return nodes;
}
double OneDimensionalNodes::getClenshawCurtisWeight(int level, int point){
    int ieffective, n = OneDimensionalMeta::getNumPoints(level, rule_clenshawcurtis);
    if (level == 0){ return 2.0; }
    if (point == 0){
        ieffective = (n-1) / 2;
    }else if (point == 1){
        ieffective = 0;
    }else if (point == 2){
        ieffective = n-1;
    }else{
        int z = point - 1, l = 1;
        while (z >>= 1){ l++; };
        ieffective = (1 + 2 *(point - OneDimensionalMeta::getNumPoints(l-1, rule_clenshawcurtis))) * (n-1) / (1 << l);
    }

    double weight = 1.0;
    double theta = ((double) ieffective) * Maths::pi / ((double) (n-1));
    for(int j=1; j<(n-1)/2; j++){
        weight -= 2.0 * std::cos(2.0 * j * theta) / ((double) (4*j*j - 1));
    }
    weight -= std::cos(2.0 * (n-1) * theta / 2.0) / ((double) (n*n - 2*n));
    weight /= (double) (n-1);

    if ((point != 1) && (point != 2)){ weight *= 2.0; }

    return weight;
}
// Clenshaw-Curtis-Zero
std::vector<double> OneDimensionalNodes::getClenshawCurtisNodesZero(int level){
    int n = OneDimensionalMeta::getNumPoints(level+1, rule_clenshawcurtis);
    std::vector<double> nodes(n-2, 0.0);
    if (level > 0){
        int count = 1;
        for(int l=2; l<=level+1; l++){
            n = OneDimensionalMeta::getNumPoints(l, rule_clenshawcurtis);
            for(int i=1; i<n; i+=2){
                nodes[count++] = std::cos(Maths::pi * ((double) (n-i-1)) / ((double) (n - 1)));
            }
        }
    }
    return nodes;
}
double OneDimensionalNodes::getClenshawCurtisWeightZero(int level, int point){
    // this should be equivalent to  return getClenshawCurtisWeight(level + 1, ((point == 0) ? 0 : point +2));
    int ieffective, n = OneDimensionalMeta::getNumPoints(level+1, rule_clenshawcurtis);
    if (level == 0){ return 4.0/3.0; }
    if (point == 0){
        ieffective = (n-1) / 2;
    }else{
        int z = point +1, l = 1;
        while (z >>= 1){ l++; };
        ieffective = (1 + 2 *(point + 2 - OneDimensionalMeta::getNumPoints(l-1, rule_clenshawcurtis))) * (n-1) / (1 << l);
    }

    double weight = 1.0;
    double theta = ((double) ieffective) * Maths::pi / ((double) (n-1));
    for(int j=1; j<(n-1)/2; j++){
        weight -= 2.0 * std::cos(2.0 * j * theta) / ((double) (4*j*j - 1));
    }
    weight -= std::cos(2.0 * (n-1) * theta / 2.0) / ((double) (n*n - 2*n));
    weight /= (double) (n-1);
    weight *= 2.0;

    return weight;
}
// Fejer-2
std::vector<double> OneDimensionalNodes::getFejer2Nodes(int level){
    int n = OneDimensionalMeta::getNumPoints(level, rule_fejer2);
    std::vector<double> nodes(n, 0.0);
    if (level > 0){
        int count = 1;
        for(int l=2; l<=level+1; l++){
            n = OneDimensionalMeta::getNumPoints(l, rule_clenshawcurtis);
            for(int i=1; i<n; i+=2){
                nodes[count++] = std::cos(Maths::pi * ((double) (n-i-1)) / ((double) (n - 1)));
            }
        }
    }
    return nodes;
}
double OneDimensionalNodes::getFejer2Weight(int level, int point){
    if (level == 0){ return 2.0; }
    int ieffective, n = OneDimensionalMeta::getNumPoints(level, rule_fejer2);
    if (point == 0){
        ieffective = (n-1) / 2;
    }else{
        int z = point + 1, l = 0;
        while (z >>= 1){ l++; };
        z = (1 << (level-l))-1;
        ieffective = z + (point - OneDimensionalMeta::getNumPoints(l-1, rule_fejer2)) * ((n-1) / (1 << (l))+1);
    }

    double weight = 1.0;
    double theta = ((double) (n-ieffective)) * Maths::pi / ((double) (n+1));
    for(int j=1; j<=(n-1)/2; j++){
        weight -= 2.0 * std::cos(2.0 * j * theta) / ((double) (4*j*j - 1));
    }
    weight -= std::cos(((double) (n+1)) * theta) / ((double) (n));
    weight *= 2.0 / ((double) (n+1));

    return weight;
}

std::vector<double> OneDimensionalNodes::getRLeja(int n){
    std::vector<double> nodes(n, 0.0);
    if (n > 1){ nodes[1] = Maths::pi; }
    if (n > 2){ nodes[2] = 0.5 * Maths::pi; }
    for(int i=3; i<n; i++){
        if (i % 2 == 0){
            nodes[i] = nodes[i-1] + Maths::pi;
        }else{
            nodes[i] = 0.5 * nodes[(i+1)/2];
        }
    }
    for(int i=0; i<n; i++){  nodes[i] = std::cos(nodes[i]);  }
    if (n > 2){ nodes[2] = 0.0; } // not sure which version is better, starting at 0 or starting at 1
    return nodes;
}
std::vector<double> OneDimensionalNodes::getRLejaCentered(int n){
    std::vector<double> nodes = getRLeja(n);
    nodes[0] = 0.0;
    if (n > 1){ nodes[1] = 1.0; }
    if (n > 2){ nodes[2] = -1.0; }
    return nodes;
}
std::vector<double> OneDimensionalNodes::getRLejaShifted(int n){
    std::vector<double> nodes(n, -0.5);
    if (n > 1){ nodes[1] = 0.5; }
    for(int i=2; i<n; i++){
        if (i % 2 == 0){
            nodes[i] = std::sqrt(0.5 * (nodes[i/2] + 1));
        }else{
            nodes[i] = -nodes[i-1];
        }
    }
    return nodes;
}

std::vector<double> OneDimensionalNodes::getFourierNodes(int level) {
    int n = OneDimensionalMeta::getNumPoints(level, rule_fourier);
    std::vector<double> nodes(n, 0.0);
    if(level > 0){ nodes[1] = 1.0/3.0; nodes[2] = 2.0/3.0; }
    size_t c = 3;
    for(int l=2; l<=level; l++) {
        n = OneDimensionalMeta::getNumPoints(l, rule_fourier);
        for(int i=1; i<n; i+=3) {
            nodes[c]   = (double) i / (double) n;
            nodes[c+1] = (double) (i+1) / (double) n;
            c += 2;
        }
    }
    return nodes;
}

}

#endif
