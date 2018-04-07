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

namespace TasGrid{

CustomTabulated::CustomTabulated(std::ostream *os) : num_levels(0), num_nodes(0), precision(0), offsets(0), nodes(0), weights(0), description(0), logstream(os){}
CustomTabulated::CustomTabulated(const char* filename, std::ostream *os) : num_levels(0), num_nodes(0), precision(0), offsets(0), nodes(0),
                                                                                  weights(0), description(0), logstream(os)
{
    std::ifstream ifs; ifs.open(filename);
    if (!read(ifs)){
        if (logstream != 0) (*logstream) << "ERROR: could not read the custom rule! Check the file format." << endl;
    }
    ifs.close();
}
CustomTabulated::CustomTabulated(const CustomTabulated &custom, std::ostream *os) : num_levels(0), num_nodes(0), precision(0), offsets(0), nodes(0),
                                                                                    weights(0), description(0), logstream(os){
    copyRule(&custom);
}
CustomTabulated::~CustomTabulated(){
    reset();
}
void CustomTabulated::reset(){
    if (num_nodes != 0){ delete[] num_nodes; num_nodes = 0; }
    if (precision != 0){ delete[] precision; precision = 0; }
    if (offsets != 0){ delete[] offsets; offsets = 0; }
    if (nodes != 0){ delete[] nodes; nodes = 0; }
    if (weights != 0){ delete[] weights; weights = 0; }
    if (description != 0){ delete description; description = 0; }
    num_levels = 0;
}

// I/O subroutines
void CustomTabulated::write(std::ofstream &ofs) const{
    ofs << "description: " << description->c_str() << endl;
    ofs << "levels: " << num_levels << endl;
    for(int i=0; i<num_levels; i++){
        ofs << num_nodes[i] << " " << precision[i] << endl;
    }
    ofs << std::scientific; ofs.precision(17);
    for(int l=0; l<num_levels; l++){
        for(int i=offsets[l]; i<offsets[l] + num_nodes[l]; i++){
            ofs << weights[i] << " " << nodes[i] << endl;
        }
    }
}
void CustomTabulated::writeBinary(std::ofstream &ofs) const{
    int num_description = description->size();
    ofs.write((char*) &num_description, sizeof(int));
    ofs.write(description->c_str(), num_description * sizeof(char));
    ofs.write((char*) &num_levels, sizeof(int));
    ofs.write((char*) num_nodes, num_levels * sizeof(int));
    ofs.write((char*) precision, num_levels * sizeof(int));
    int total_points = offsets[num_levels-1] + num_nodes[num_levels-1];
    ofs.write((char*) weights, total_points * sizeof(double));
    ofs.write((char*) nodes, total_points * sizeof(double));
}
bool CustomTabulated::read(std::ifstream &ifs){
    reset();

    std::string T;
    char dummy;
    ifs >> T;
    if (!(T.compare("description:") == 0)){ if (logstream != 0) (*logstream) << "ERROR: wrong file format of custom tables on line 1" << endl; ifs.close(); return false; }
    ifs.get(dummy);
    description = new std::string;
    getline(ifs, *description);

    ifs >> T;
    if (!(T.compare("levels:") == 0)){ if (logstream != 0) (*logstream) << "ERROR: wrong file format of custom tables on line 2" << endl; ifs.close(); return false; }
    ifs >> num_levels;

    num_nodes = new int[num_levels];
    precision = new int[num_levels];
    for(int i=0; i<num_levels; i++){
        ifs >> num_nodes[i] >> precision[i];
    }

    int total_points = 0;
    offsets = new int[num_levels];
    for(int i=0; i<num_levels; i++){
        offsets[i] = total_points;
        total_points += num_nodes[i];
    }

    nodes   = new double[total_points];
    weights = new double[total_points];
    for(int i=0; i<total_points; i++){
        ifs >> weights[i] >> nodes[i];
    }

    return true;
}
bool CustomTabulated::readBinary(std::ifstream &ifs){
    reset();

    int num_description = 0;
    ifs.read((char*) &num_description, sizeof(int));
    char *desc = new char[num_description+1];
    ifs.read(desc, num_description);
    desc[num_description] = '\0';
    description = new std::string;
    *description = desc;
    delete[] desc;

    ifs.read((char*) &num_levels, sizeof(int));
    num_nodes = new int[num_levels];
    precision = new int[num_levels];
    ifs.read((char*) num_nodes, num_levels * sizeof(int));
    ifs.read((char*) precision, num_levels * sizeof(int));

    int total_points = 0;
    offsets = new int[num_levels];
    for(int i=0; i<num_levels; i++){
        offsets[i] = total_points;
        total_points += num_nodes[i];
    }

    nodes   = new double[total_points];
    weights = new double[total_points];
    ifs.read((char*) weights, total_points * sizeof(double));
    ifs.read((char*) nodes, total_points * sizeof(double));

    return true;
}

void CustomTabulated::copyRule(const CustomTabulated *custom){
    description = new std::string(*(custom->description));

    num_levels = custom->num_levels;

    num_nodes = new int[num_levels]; std::copy(custom->num_nodes, custom->num_nodes+custom->num_levels, num_nodes);
    precision = new int[num_levels]; std::copy(custom->precision, custom->precision+custom->num_levels, precision);

    int total_points = 0;
    offsets = new int[num_levels];
    for(int i=0; i<num_levels; i++){
        offsets[i] = total_points;
        total_points += num_nodes[i];
    }

    nodes   = new double[total_points]; std::copy(custom->nodes, custom->nodes + total_points, nodes);
    weights = new double[total_points]; std::copy(custom->weights, custom->weights + total_points, weights);
}

int CustomTabulated::getNumLevels() const{ return num_levels; }
int CustomTabulated::getNumPoints(int level) const{
    if (level >= num_levels){
        if (logstream != 0) (*logstream) << "ERROR: requested custom rule with level " << level << " the tabulated rules end at " << num_levels-1 << endl;
        return num_nodes[num_levels-1];
    }
    return num_nodes[level];
}
int CustomTabulated::getIExact(int level) const{
    if (level >= num_levels){
        if (logstream != 0) (*logstream) << "ERROR: requested custom rule with level " << level << " the tabulated rules end at " << num_levels-1 << endl;
        return 1 << 30;
    }
    return num_nodes[level] -1;
}
int CustomTabulated::getQExact(int level) const{
    if (level >= num_levels){
        if (logstream != 0) (*logstream) << "ERROR: requested custom rule with level " << level << " the tabulated rules end at " << num_levels-1 << endl;
        return 1 << 30;
    }
    return precision[level];
}

void CustomTabulated::getWeightsNodes(int level, double* &w, double* &x) const{
    if (w != 0) delete[] w;
    if (x != 0) delete[] x;
    w = new double[num_nodes[level]];
    x = new double[num_nodes[level]];
    int off = offsets[level];
    for(int i=0; i<num_nodes[level]; i++){
        x[i] = nodes[off + i];
        w[i] = weights[off + i];
    }
}
const char* CustomTabulated::getDescription() const{  return description->c_str();  }

OneDimensionalMeta::OneDimensionalMeta(): custom(0) {}
OneDimensionalMeta::OneDimensionalMeta(const CustomTabulated *ccustom) : custom(ccustom) {}
OneDimensionalMeta::~OneDimensionalMeta(){}
const CustomTabulated* OneDimensionalMeta::getCustom() const{  return custom;  }

int OneDimensionalMeta::getNumPoints(int level, TypeOneDRule rule) const{
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

        case rule_rlejashifteddouble: return (1 << (level+1));

        case rule_clenshawcurtis0:
        case rule_gausspatterson:     return ((1 << (level+1)) - 1);

        case rule_clenshawcurtis:     return (level == 0) ? 1 : ((1 << level) + 1);
        case rule_rlejadouble2:       if (level < 3){ return getNumPoints(level, rule_clenshawcurtis); }; lcc = 2 + (level-3)/2; return (getNumPoints(lcc,rule_clenshawcurtis) + ((getNumPoints(lcc+1,rule_clenshawcurtis)-getNumPoints(lcc,rule_clenshawcurtis))/2) * ((level-3)%2 +1));
        case rule_rlejadouble4:       if (level < 3){ return getNumPoints(level, rule_clenshawcurtis); }; lcc = 2 + (level-3)/4; return (getNumPoints(lcc,rule_clenshawcurtis) + ((getNumPoints(lcc+1,rule_clenshawcurtis)-getNumPoints(lcc,rule_clenshawcurtis))/4) * ((level-3)%4 +1));
        case rule_fejer2:             return ((1 << (level+1)) - 1);

        case rule_customtabulated:    return custom->getNumPoints(level);

        default:
            return level;
    }
}
int OneDimensionalMeta::getIExact(int level, TypeOneDRule rule) const{
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

        case rule_rlejashifteddouble: return ((1 << (level+1)) -1);

        case rule_gausspatterson:
        case rule_fejer2:             return ((1 << (level+1)) - 2);

        case rule_clenshawcurtis:     return (level == 0) ? 0 : ((1 << level));
        case rule_clenshawcurtis0:    return ((1 << (level+1)) +1);
        case rule_rlejadouble2:       return getNumPoints(level,rule_rlejadouble2)-1;
        case rule_rlejadouble4:       return getNumPoints(level,rule_rlejadouble4)-1;
        case rule_customtabulated:    return custom->getIExact(level);
        default:
            return level;
    }
}
int OneDimensionalMeta::getQExact(int level, TypeOneDRule rule) const{
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

        case rule_rlejashifteddouble: return ((1 << (level+1)) -1);

        case rule_gausspatterson:     return (level == 0) ? 1 : (3*(1 << level) - 1);
        case rule_clenshawcurtis:     return (level == 0) ? 1 : ((1 << level) + 1);
        case rule_clenshawcurtis0:    return (level == 0) ? 1 : ((1 << (level+1)) + 1);
        case rule_chebyshev:          return level+1;
        case rule_rlejadouble2:       return getNumPoints(level,rule_rlejadouble2);
        case rule_rlejadouble4:       return getNumPoints(level,rule_rlejadouble4)-1;
        case rule_fejer2:             return ((1 << (level+1)) - 1);
        case rule_customtabulated:    return custom->getQExact(level);
        default:
            return level;
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
    return !((rule == rule_semilocalp) || (rule == rule_localp0) || (rule == rule_localp) || (rule == rule_wavelet) || (rule == rule_none));
}
bool OneDimensionalMeta::isSingleNodeGrowth(TypeOneDRule rule){
    return ((rule == rule_leja) || (rule == rule_rleja) || (rule == rule_rlejashifted) || (rule == rule_maxlebesgue) || (rule == rule_minlebesgue) || (rule == rule_mindelta) ||
         (rule == rule_gausslegendre) || (rule == rule_gausschebyshev1) || (rule == rule_gausschebyshev2) || (rule == rule_gaussgegenbauer) ||
         (rule == rule_gaussjacobi) || (rule == rule_gausslaguerre) || (rule == rule_gausshermite));
}
bool OneDimensionalMeta::isLocalPolynomial(TypeOneDRule rule){
    return ((rule == rule_localp) || (rule == rule_localp0) || (rule == rule_semilocalp));
}
bool OneDimensionalMeta::isWavelet(TypeOneDRule rule){
    return (rule == rule_wavelet);
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
        case rule_semilocalp:         return "Semi-Local polynomials";
        case rule_wavelet:            return "Wavelets";
        default:
            return "unknown";
    }
}
const char* OneDimensionalMeta::getIORuleString(TypeOneDRule rule){
    switch (rule){
        case rule_clenshawcurtis:     return "clenshaw-curtis";
        case rule_clenshawcurtis0:    return "clenshaw-curtis-zero";
        case rule_chebyshev:          return "chebyshev";
        case rule_chebyshevodd:       return "chebyshev-odd";
        case rule_gausslegendre:      return "gauss-legendre";
        case rule_gausslegendreodd:   return "gauss-legendre-odd";
        case rule_gausspatterson:     return "gauss-patterson";
        case rule_leja:               return "leja";
        case rule_lejaodd:            return "leja-odd";
        case rule_rleja:              return "rleja";
        case rule_rlejaodd:           return "rleja-odd";
        case rule_rlejadouble2:       return "rleja-double2";
        case rule_rlejadouble4:       return "rleja-double4";
        case rule_rlejashifted:       return "rleja-shifted";
        case rule_rlejashiftedeven:   return "rleja-shifted-even";
        case rule_rlejashifteddouble: return "rleja-shifted-double";
        case rule_maxlebesgue:        return "max-lebesgue";
        case rule_maxlebesgueodd:     return "max-lebesgue-odd";
        case rule_minlebesgue:        return "min-lebesgue";
        case rule_minlebesgueodd:     return "min-lebesgue-odd";
        case rule_mindelta:           return "min-delta";
        case rule_mindeltaodd:        return "min-delta-odd";
        case rule_gausschebyshev1:    return "gauss-chebyshev1";
        case rule_gausschebyshev1odd: return "gauss-chebyshev1-odd";
        case rule_gausschebyshev2:    return "gauss-chebyshev2";
        case rule_gausschebyshev2odd: return "gauss-chebyshev2-odd";
        case rule_fejer2:             return "fejer2";
        case rule_gaussgegenbauer:    return "gauss-gegenbauer";
        case rule_gaussgegenbauerodd: return "gauss-gegenbauer-odd";
        case rule_gaussjacobi:        return "gauss-jacobi";
        case rule_gaussjacobiodd:     return "gauss-jacobi-odd";
        case rule_gausslaguerre:      return "gauss-laguerre";
        case rule_gausslaguerreodd:   return "gauss-laguerre-odd";
        case rule_gausshermite:       return "gauss-hermite";
        case rule_gausshermiteodd:    return "gauss-hermite-odd";
        case rule_customtabulated:    return "custom-tabulated";
        case rule_localp:             return "localp";
        case rule_localp0:            return "localp-zero";
        case rule_semilocalp:         return "semi-localp";
        case rule_wavelet:            return "wavelet";
        default:
            return "unknown";
    }
}
TypeOneDRule OneDimensionalMeta::getIORuleString(const char *name){
    if (strcmp(name, "clenshaw-curtis") == 0){
        return rule_clenshawcurtis;
    }else if (strcmp(name, "clenshaw-curtis-zero") == 0){
        return rule_clenshawcurtis0;
    }else if (strcmp(name, "chebyshev") == 0){
        return rule_chebyshev;
    }else if (strcmp(name, "chebyshev-odd") == 0){
        return rule_chebyshevodd;
    }else if (strcmp(name, "gauss-legendre") == 0){
        return rule_gausslegendre;
    }else if (strcmp(name, "gauss-legendre-odd") == 0){
        return rule_gausslegendreodd;
    }else if (strcmp(name, "gauss-patterson") == 0){
        return rule_gausspatterson;
    }else if (strcmp(name, "leja") == 0){
        return rule_leja;
    }else if (strcmp(name, "leja-odd") == 0){
        return rule_lejaodd;
    }else if (strcmp(name, "rleja") == 0){
        return rule_rleja;
    }else if (strcmp(name, "rleja-double2") == 0){
        return rule_rlejadouble2;
    }else if (strcmp(name, "rleja-double4") == 0){
        return rule_rlejadouble4;
    }else if (strcmp(name, "rleja-odd") == 0){
        return rule_rlejaodd;
    }else if (strcmp(name, "rleja-shifted") == 0){
        return rule_rlejashifted;
    }else if (strcmp(name, "rleja-shifted-even") == 0){
        return rule_rlejashiftedeven;
    }else if (strcmp(name, "rleja-shifted-double") == 0){
        return rule_rlejashifteddouble;
    }else if (strcmp(name, "max-lebesgue") == 0){
        return rule_maxlebesgue;
    }else if (strcmp(name, "max-lebesgue-odd") == 0){
        return rule_maxlebesgueodd;
    }else if (strcmp(name, "min-lebesgue") == 0){
        return rule_minlebesgue;
    }else if (strcmp(name, "min-lebesgue-odd") == 0){
        return rule_minlebesgueodd;
    }else if (strcmp(name, "min-delta") == 0){
        return rule_mindelta;
    }else if (strcmp(name, "min-delta-odd") == 0){
        return rule_mindeltaodd;
    }else if (strcmp(name, "gauss-chebyshev1") == 0){
        return rule_gausschebyshev1;
    }else if (strcmp(name, "gauss-chebyshev1-odd") == 0){
        return rule_gausschebyshev1odd;
    }else if (strcmp(name, "gauss-chebyshev2") == 0){
        return rule_gausschebyshev2;
    }else if (strcmp(name, "gauss-chebyshev2-odd") == 0){
        return rule_gausschebyshev2odd;
    }else if (strcmp(name, "fejer2") == 0){
        return rule_fejer2;
    }else if (strcmp(name, "gauss-gegenbauer") == 0){
        return rule_gaussgegenbauer;
    }else if (strcmp(name, "gauss-gegenbauer-odd") == 0){
        return rule_gaussgegenbauerodd;
    }else if (strcmp(name, "gauss-jacobi") == 0){
        return rule_gaussjacobi;
    }else if (strcmp(name, "gauss-jacobi-odd") == 0){
        return rule_gaussjacobiodd;
    }else if (strcmp(name, "gauss-laguerre") == 0){
        return rule_gausslaguerre;
    }else if (strcmp(name, "gauss-laguerre-odd") == 0){
        return rule_gausslaguerreodd;
    }else if (strcmp(name, "gauss-hermite") == 0){
        return rule_gausshermite;
    }else if (strcmp(name, "gauss-hermite-odd") == 0){
        return rule_gausshermiteodd;
    }else if (strcmp(name, "custom-tabulated") == 0){
        return rule_customtabulated;
    }else  if (strcmp(name, "localp") == 0){
        return rule_localp;
    }else  if (strcmp(name, "localp-zero") == 0){
        return rule_localp0;
    }else  if (strcmp(name, "semi-localp") == 0){
        return rule_semilocalp;
    }else if (strcmp(name, "wavelet") == 0){
        return rule_wavelet;
    }else{
        return rule_none;
    }
}
TypeOneDRule OneDimensionalMeta::getIORuleInt(int index){
    switch (index){
        case  1: return rule_clenshawcurtis;
        case  2: return rule_clenshawcurtis0;
        case  3: return rule_chebyshev;
        case  4: return rule_chebyshevodd;
        case  5: return rule_gausslegendre;
        case  6: return rule_gausslegendreodd;
        case  7: return rule_gausspatterson;
        case  8: return rule_leja;
        case  9: return rule_lejaodd;
        case 10: return rule_rleja;
        case 11: return rule_rlejaodd;
        case 12: return rule_rlejadouble2;
        case 13: return rule_rlejadouble4;
        case 14: return rule_rlejashifted;
        case 15: return rule_rlejashiftedeven;
        case 16: return rule_rlejashifteddouble;
        case 17: return rule_maxlebesgue;
        case 18: return rule_maxlebesgueodd;
        case 19: return rule_minlebesgue;
        case 20: return rule_minlebesgueodd;
        case 21: return rule_mindelta;
        case 22: return rule_mindeltaodd;
        case 23: return rule_gausschebyshev1;
        case 24: return rule_gausschebyshev1odd;
        case 25: return rule_gausschebyshev2;
        case 26: return rule_gausschebyshev2odd;
        case 27: return rule_fejer2;
        case 28: return rule_gaussgegenbauer;
        case 29: return rule_gaussgegenbauerodd;
        case 30: return rule_gaussjacobi;
        case 31: return rule_gaussjacobiodd;
        case 32: return rule_gausslaguerre;
        case 33: return rule_gausslaguerreodd;
        case 34: return rule_gausshermite;
        case 35: return rule_gausshermiteodd;
        case 36: return rule_customtabulated;
        case 37: return rule_localp;
        case 38: return rule_localp0;
        case 39: return rule_semilocalp;
        case 40: return rule_wavelet;
        default:
            return rule_none;
    }
}
int OneDimensionalMeta::getIORuleInt(TypeOneDRule rule){
    switch (rule){
        case rule_clenshawcurtis:     return  1;
        case rule_clenshawcurtis0:    return  2;
        case rule_chebyshev:          return  3;
        case rule_chebyshevodd:       return  4;
        case rule_gausslegendre:      return  5;
        case rule_gausslegendreodd:   return  6;
        case rule_gausspatterson:     return  7;
        case rule_leja:               return  8;
        case rule_lejaodd:            return  9;
        case rule_rleja:              return 10;
        case rule_rlejaodd:           return 11;
        case rule_rlejadouble2:       return 12;
        case rule_rlejadouble4:       return 13;
        case rule_rlejashifted:       return 14;
        case rule_rlejashiftedeven:   return 15;
        case rule_rlejashifteddouble: return 16;
        case rule_maxlebesgue:        return 17;
        case rule_maxlebesgueodd:     return 18;
        case rule_minlebesgue:        return 19;
        case rule_minlebesgueodd:     return 20;
        case rule_mindelta:           return 21;
        case rule_mindeltaodd:        return 22;
        case rule_gausschebyshev1:    return 23;
        case rule_gausschebyshev1odd: return 24;
        case rule_gausschebyshev2:    return 25;
        case rule_gausschebyshev2odd: return 26;
        case rule_fejer2:             return 27;
        case rule_gaussgegenbauer:    return 28;
        case rule_gaussgegenbauerodd: return 29;
        case rule_gaussjacobi:        return 30;
        case rule_gaussjacobiodd:     return 31;
        case rule_gausslaguerre:      return 32;
        case rule_gausslaguerreodd:   return 33;
        case rule_gausshermite:       return 34;
        case rule_gausshermiteodd:    return 35;
        case rule_customtabulated:    return 36;
        case rule_localp:             return 37;
        case rule_localp0:            return 38;
        case rule_semilocalp:         return 39;
        case rule_wavelet:            return 40;
        default:
            return 0;
    }
}

TypeDepth OneDimensionalMeta::getIOTypeString(const char *name){
    if (strcmp(name, "level") == 0){
        return type_level;
    }else if (strcmp(name, "curved") == 0){
        return type_curved;
    }else if (strcmp(name, "iptotal") == 0){
        return type_iptotal;
    }else if (strcmp(name, "ipcurved") == 0){
        return type_ipcurved;
    }else if (strcmp(name, "qptotal") == 0){
        return type_qptotal;
    }else if (strcmp(name, "qpcurved") == 0){
        return type_qpcurved;
    }else if (strcmp(name, "hyperbolic") == 0){
        return type_hyperbolic;
    }else if (strcmp(name, "iphyperbolic") == 0){
        return type_iphyperbolic;
    }else if (strcmp(name, "qphyperbolic") == 0){
        return type_qphyperbolic;
    }else if (strcmp(name, "tensor") == 0){
        return type_tensor;
    }else if (strcmp(name, "iptensor") == 0){
        return type_iptensor;
    }else if (strcmp(name, "qptensor") == 0){
        return type_qptensor;
    }else{
        //cerr << "ERROR: " << argv[k] << " is not a valid type!!!  For help see: ./tasgrid -help" << endl << endl;
        return type_none;
    }
}
TypeDepth OneDimensionalMeta::getIOTypeInt(int type){
    switch (type){
        case  1: return type_level;
        case  2: return type_curved;
        case  3: return type_iptotal;
        case  4: return type_ipcurved;
        case  5: return type_qptotal;
        case  6: return type_qpcurved;
        case  7: return type_hyperbolic;
        case  8: return type_iphyperbolic;
        case  9: return type_qphyperbolic;
        case 10: return type_tensor;
        case 11: return type_iptensor;
        case 12: return type_qptensor;
        default:
            return type_none;
    }
}

TypeRefinement OneDimensionalMeta::getIOTypeRefinementString(const char *name){
    if (strcmp(name, "classic") == 0){
        return refine_classic;
    }else if (strcmp(name, "parents") == 0){
        return refine_parents_first;
    }else if (strcmp(name, "direction") == 0){
        return refine_direction_selective;
    }else if (strcmp(name, "fds") == 0){
        return refine_fds;
    }else{
        return refine_none;
    }
}
TypeRefinement OneDimensionalMeta::getIOTypeRefinementInt(int ref){
    switch (ref){
        case  1: return refine_classic;
        case  2: return refine_parents_first;
        case  3: return refine_direction_selective;
        case  4: return refine_fds;
        default:
            return refine_none;
    }
}

// Gauss-Legendre
void OneDimensionalNodes::getGaussLegendre(int m, double* &w, double* &x){
    if (w != 0){ delete[] w; }; w = new double[m];
    if (x != 0){ delete[] x; }; x = new double[m];

    double *s = new double[m];
    for(int i=0; i<m; i++){ x[i] = w[i] = s[i] = 0.0; }

    for(int i=0; i<m; i++){
        s[i] = sqrt((double) ((i+1)*(i+1)) / ((double) (4*(i+1)*(i+1) - 1)));
    }
    w[0] = sqrt(2.0);

    TasmanianTridiagonalSolver::decompose(m, x, s, w);

    delete[] s;
}

// Chebyshev
void OneDimensionalNodes::getChebyshev(int m, double* &w, double* &x){
    // get Clanshaw-Curtis quadrature points
    if (w != 0){ delete[] w; }
    if (x != 0){ delete[] x; }
    w = new double[m];
    x = new double[m];
    int i, j;
    double b;
    if (m == 1){
        w[0] = 2.0; x[0] = 0.0;
        return;
    };

    for(i=0; i<m; i++){
        x[i] = cos(((double) (m-i-1)) * M_PI / ((double) (m-1)));
    };

    // may also have to set the mid-point to 0.0
    x[0] = -1.0; x[m-1] = 1.0;

    for(i=0; i<m; i++){
        w[i] = 1.0;
        double theta = ((double) i) * M_PI / ((double) (m-1));
        for(j=1; j<=(m-1)/2; j++){
            if (2*j == (m-1)){
                b = 1.0;
            }else{
                b = 2.0;
            };
            w[i] = w[i] - b * cos(2.0 * j * theta) / ((double) (4*j*j - 1));
        };
    };

    w[0] = w[0] / ((double) (m-1));
    for(i=1; i<m-1; i++){
        w[i] = 2.0 * w[i] / ((double) (m-1));
    };
    w[m-1] = w[m-1] / ((double) (m-1));
}

// get Gauss-Chebyshev type 1 quadrature points
void OneDimensionalNodes::getGaussChebyshev1(int m, double* &w, double* &x){
    if (w != 0){ delete[] w; }
    if (x != 0){ delete[] x; }
    w = new double[m];
    x = new double[m];

    for(int i=0; i<m; i++){
        x[m-i-1] = cos(M_PI*(2*i+1) / (2*((double)m)));
        w[i] = M_PI / m;
    }
}
// get Gauss-Chebyshev-type2 quadrature points
void OneDimensionalNodes::getGaussChebyshev2(int m, double* &w, double* &x){
    if (w != 0){ delete[] w; }
    if (x != 0){ delete[] x; }
    w = new double[m];
    x = new double[m];

    for(int i=0; i<m; i++){
        double theta = M_PI*((double)(i+1))/((double)(m+1));
        x[m-i-1] = cos(theta);
        w[i] = (M_PI / ((double)(m+1))) * sin(theta)* sin(theta);
    }
}
// get Gauss-Jacobi quadrature points
void OneDimensionalNodes::getGaussJacobi(int m, double* &w, double* &x, double alpha, double beta){
    if (w != 0){ delete[] w; }
    if (x != 0){ delete[] x; }
    w = new double[m];
    x = new double[m];

    double *s = new double[m];

    for(int i=0; i<m; i++){ x[i] = w[i] = s[i] = 0.0; }

    double ab = alpha + beta;

    w[0] = sqrt(pow(2.0, 1.0 + ab) * tgamma(alpha + 1.0) * tgamma(beta + 1.0) / tgamma(2.0 + ab));

    x[0] = (beta - alpha) / (2.0 + ab);
    s[0] = sqrt(4.0 * (1.0 + alpha) * (1.0 + beta) / ((3.0 + ab) * (2.0 + ab) * (2.0 + ab)));
    for(int i=1; i<m; i++){
        double di = (double) (i+1);
        x[i] = (beta*beta - alpha*alpha) / ((2.0*di + ab -2.0)*(2.0*di + ab));
        s[i] = sqrt(4.0 * di * (di + alpha) * (di + beta) * (di + ab)/ (((2.0*di + ab)*(2.0*di + ab) - 1.0) * (2.0*di + ab) * (2.0*di + ab)));
    }
    s[m-1] = 0.0;

    TasmanianTridiagonalSolver::decompose(m, x, s, w);

    delete[] s;
}
// get Gauss-Hermite quadrature points
void OneDimensionalNodes::getGaussHermite(int m, double* &w, double* &x, double alpha){
    if (w != 0){ delete[] w; }
    if (x != 0){ delete[] x; }
    w = new double[m];
    x = new double[m];

    double *s = new double[m];

    for(int i=0; i<m; i++){ x[i] = w[i] = s[i] = 0.0; }

    w[0] = sqrt(tgamma(0.5 * (alpha + 1.0)));

    for(int i=0; i<m; i++){
        double di = (double) (i+1);
        s[i] = sqrt(0.5 * (di + alpha * ((double) ((i+1)%2))));
    }
    s[m-1] = 0.0;

    TasmanianTridiagonalSolver::decompose(m, x, s, w);

    delete[] s;
}
// get Gauss-Laguerre quadrature points
void OneDimensionalNodes::getGaussLaguerre(int m, double* &w, double* &x, double alpha){
    if (w != 0){ delete[] w; }
    if (x != 0){ delete[] x; }
    w = new double[m];
    x = new double[m];

    double *s = new double[m];

    for(int i=0; i<m; i++){ x[i] = w[i] = s[i] = 0.0; }

    w[0] = sqrt(tgamma(alpha + 1.0));

    for(int i=0; i<m; i++){
        double di = (double) (i+1);
        x[i] = 2.0 * di - 1.0 + alpha;
        s[i] = sqrt(di * (di + alpha));
    }
    s[m-1] = 0.0;

    TasmanianTridiagonalSolver::decompose(m, x, s, w);

    delete[] s;
}

// Clenshaw-Curtis
double* OneDimensionalNodes::getClenshawCurtisNodes(int level){
    OneDimensionalMeta meta;
    int n = meta.getNumPoints(level, rule_clenshawcurtis);
    double* x = new double[n];
    x[0] = 0.0;
    if (level > 0){
        x[1] = -1.0; x[2] = 1.0;
        int count = 3;
        for(int l=2; l<=level; l++){
            n = meta.getNumPoints(l, rule_clenshawcurtis);
            for(int i=1; i<n; i+=2){
                  x[count++] = cos(M_PI * ((double) (n-i-1)) / ((double) (n - 1)));
            }
        }
    }
    return x;
}
double OneDimensionalNodes::getClenshawCurtisWeight(int level, int point){
    OneDimensionalMeta meta;
    int ieffective, n = meta.getNumPoints(level, rule_clenshawcurtis);
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
        ieffective = (1 + 2 *(point - meta.getNumPoints(l-1, rule_clenshawcurtis))) * (n-1) / (1 << l);
    }

    double weight = 1.0;
    double theta = ((double) ieffective) * M_PI / ((double) (n-1));
    for(int j=1; j<(n-1)/2; j++){
        weight -= 2.0 * cos(2.0 * j * theta) / ((double) (4*j*j - 1));
    }
    weight -= cos(2.0 * (n-1) * theta / 2.0) / ((double) (n*n - 2*n));
    weight /= (double) (n-1);

    if ((point != 1) && (point != 2)){ weight *= 2.0; }

    return weight;
}
// Clenshaw-Curtis-Zero
double* OneDimensionalNodes::getClenshawCurtisNodesZero(int level){
    OneDimensionalMeta meta;
    int n = meta.getNumPoints(level+1, rule_clenshawcurtis);
    double* x = new double[n-2];
    x[0] = 0.0;
    if (level > 0){
        int count = 1;
        for(int l=2; l<=level+1; l++){
            n = meta.getNumPoints(l, rule_clenshawcurtis);
            for(int i=1; i<n; i+=2){
                  x[count++] = cos(M_PI * ((double) (n-i-1)) / ((double) (n - 1)));
            }
        }
    }
    return x;
}
double OneDimensionalNodes::getClenshawCurtisWeightZero(int level, int point){
    // this should be equivalent to  return getClenshawCurtisWeight(level + 1, ((point == 0) ? 0 : point +2));
    OneDimensionalMeta meta;
    int ieffective, n = meta.getNumPoints(level+1, rule_clenshawcurtis);
    if (level == 0){ return 4.0/3.0; }
    if (point == 0){
        ieffective = (n-1) / 2;
    }else{
        int z = point +1, l = 1;
        while (z >>= 1){ l++; };
        ieffective = (1 + 2 *(point + 2 - meta.getNumPoints(l-1, rule_clenshawcurtis))) * (n-1) / (1 << l);
    }

    double weight = 1.0;
    double theta = ((double) ieffective) * M_PI / ((double) (n-1));
    for(int j=1; j<(n-1)/2; j++){
        weight -= 2.0 * cos(2.0 * j * theta) / ((double) (4*j*j - 1));
    }
    weight -= cos(2.0 * (n-1) * theta / 2.0) / ((double) (n*n - 2*n));
    weight /= (double) (n-1);
    weight *= 2.0;

    return weight;
}
// Fejer-2
double* OneDimensionalNodes::getFejer2Nodes(int level){
    OneDimensionalMeta meta;
    int n = meta.getNumPoints(level, rule_fejer2);
    double *x = new double[n];
    x[0] = 0.0;
    if (level > 0){
        int count = 1;
        for(int l=2; l<=level+1; l++){
            n = meta.getNumPoints(l, rule_clenshawcurtis);
            for(int i=1; i<n; i+=2){
                  x[count++] = cos(M_PI * ((double) (n-i-1)) / ((double) (n - 1)));
            }
        }
    }
    return x;
}
double OneDimensionalNodes::getFejer2Weight(int level, int point){
    if (level == 0){ return 2.0; }
    OneDimensionalMeta meta;
    int ieffective, n = meta.getNumPoints(level, rule_fejer2);
    if (point == 0){
        ieffective = (n-1) / 2;
    }else{
        int z = point + 1, l = 0;
        while (z >>= 1){ l++; };
        z = (1 << (level-l))-1;
        ieffective = z + (point - meta.getNumPoints(l-1, rule_fejer2)) * ((n-1) / (1 << (l))+1);
    }


    double weight = 1.0;
    double theta = ((double) (n-ieffective)) * M_PI / ((double) (n+1));
    for(int j=1; j<=(n-1)/2; j++){
        weight -= 2.0 * cos(2.0 * j * theta) / ((double) (4*j*j - 1));
    }
    weight -= cos(((double) (n+1)) * theta) / ((double) (n));
    weight *= 2.0 / ((double) (n+1));

    return weight;
}

double* OneDimensionalNodes::getRLeja(int n){
    double* nodes = new double[n];
    nodes[0] = 0.0;
    if (n > 1){ nodes[1] = M_PI; }
    if (n > 2){ nodes[2] = 0.5 * M_PI; }
    for(int i=3; i<n; i++){
        if (i % 2 == 0){
            nodes[i] = nodes[i-1] + M_PI;
        }else{
            nodes[i] = 0.5 * nodes[(i+1)/2];
        }
    }
    for(int i=0; i<n; i++){  nodes[i] = cos(nodes[i]);  }
    if (n > 2){ nodes[2] = 0.0; } // not sure which version is better, starting at 0 or starting at 1
    return nodes;
}
double* OneDimensionalNodes::getRLejaCentered(int n){
    double* nodes = getRLeja(n);
    nodes[0] = 0.0;
    if (n > 1){ nodes[1] = 1.0; }
    if (n > 2){ nodes[2] = -1.0; }
    return nodes;
}
double* OneDimensionalNodes::getRLejaShifted(int n){
    double* nodes = new double[n];
    nodes[0] = -0.5;
    if (n > 1){ nodes[1] = 0.5; }
    for(int i=2; i<n; i++){
        if (i % 2 == 0){
            nodes[i] = sqrt(0.5 * (nodes[i/2] + 1));
        }else{
            nodes[i] = -nodes[i-1];
        }
    }
    return nodes;
}

}

#endif
