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

#ifndef __TASMANIAN_IOHELPERS_HPP
#define __TASMANIAN_IOHELPERS_HPP

#include "tsgEnumerates.hpp"

/*!
 * \internal
 * \file tsgIOHelpers.hpp
 * \brief Templates to simply file I/O.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianIO
 *
 * Several templates that simplify the I/O of Tasmanian.
 * Commonly used operations are lumped into templates with simple binary/ascii switch.
 * \endinternal
 */

/*!
 * \internal
 * \ingroup TasmanianSG
 * \addtogroup TasmanianIO Templates for common I/O methods, simple binary/ascii switch provided for all templates
 * \endinternal
 */

namespace TasGrid{

/*!
 * \ingroup TasmanianSG
 * \brief Constant allowing for more expressive selection of ascii and binary mode in IO methods.
 */
constexpr bool mode_ascii = false;

/*!
 * \ingroup TasmanianSG
 * \brief Constant allowing for more expressive selection of ascii and binary mode in IO methods.
 */
constexpr bool mode_binary = true;


/*!
 * \internal
 * \ingroup TasmanianIO
 * \brief Collection of I/O handling templates.
 */
namespace IO{

/*!
 * \ingroup TasmanianIO
 * \brief Type indicating ascii I/O mode.
 */
struct mode_ascii_type{};

/*!
 * \ingroup TasmanianIO
 * \brief Type indicating binary I/O mode.
 */
struct mode_binary_type{};

/*!
 * \internal
 * \ingroup TasmanianIO
 * \brief Indicate the type of padding to use, none, space, new-line, etc.
 *
 * Ascii formats are human readable (for small data), but to achieve that
 * new lines and spaces must be added in the right place using the IOPad enumerate.
 * \endinternal
 */
enum IOPad{
    //! \brief Do not add padding.
    pad_none,
    //! \brief Pad with space at the end.
    pad_rspace,
    //! \brief Pad with space at the beginning.
    pad_lspace,
    //! \brief Pad with new line.
    pad_line,
    //! \brief Pad with space if the flag is true, newline if false.
    pad_auto
};

/*!
 * \ingroup TasmanianIO
 * \brief Creates a map with \b std::string rule names (used by C/Python/CLI) mapped to \b TypeOneDRule enums.
 */
inline std::map<std::string, TypeOneDRule> getStringRuleMap(){
    return std::initializer_list<std::pair<std::string const, TypeOneDRule>>{
        {"none",                 rule_none},
        {"clenshaw-curtis",      rule_clenshawcurtis},
        {"clenshaw-curtis-zero", rule_clenshawcurtis0},
        {"chebyshev",            rule_chebyshev},
        {"chebyshev-odd",        rule_chebyshevodd},
        {"gauss-legendre",       rule_gausslegendre},
        {"gauss-legendre-odd",   rule_gausslegendreodd},
        {"gauss-patterson",      rule_gausspatterson},
        {"leja",                 rule_leja},
        {"leja-odd",             rule_lejaodd},
        {"rleja",                rule_rleja},
        {"rleja-double2",        rule_rlejadouble2},
        {"rleja-double4",        rule_rlejadouble4},
        {"rleja-odd",            rule_rlejaodd},
        {"rleja-shifted",        rule_rlejashifted},
        {"rleja-shifted-even",   rule_rlejashiftedeven},
        {"rleja-shifted-double", rule_rlejashifteddouble},
        {"max-lebesgue",         rule_maxlebesgue},
        {"max-lebesgue-odd",     rule_maxlebesgueodd},
        {"min-lebesgue",         rule_minlebesgue},
        {"min-lebesgue-odd",     rule_minlebesgueodd},
        {"min-delta",            rule_mindelta},
        {"min-delta-odd",        rule_mindeltaodd},
        {"gauss-chebyshev1",     rule_gausschebyshev1},
        {"gauss-chebyshev1-odd", rule_gausschebyshev1odd},
        {"gauss-chebyshev2",     rule_gausschebyshev2},
        {"gauss-chebyshev2-odd", rule_gausschebyshev2odd},
        {"fejer2",               rule_fejer2},
        {"gauss-gegenbauer",     rule_gaussgegenbauer},
        {"gauss-gegenbauer-odd", rule_gaussgegenbauerodd},
        {"gauss-jacobi",         rule_gaussjacobi},
        {"gauss-jacobi-odd",     rule_gaussjacobiodd},
        {"gauss-laguerre",       rule_gausslaguerre},
        {"gauss-laguerre-odd",   rule_gausslaguerreodd},
        {"gauss-hermite",        rule_gausshermite},
        {"gauss-hermite-odd",    rule_gausshermiteodd},
        {"custom-tabulated",     rule_customtabulated},
        {"localp",               rule_localp},
        {"localp-zero",          rule_localp0},
        {"localp-boundary",      rule_localpb},
        {"semi-localp",          rule_semilocalp},
        {"wavelet",              rule_wavelet},
        {"fourier",              rule_fourier}};
}

/*!
 * \ingroup TasmanianIO
 * \brief Map the string rule name to the enumerate, used in ASCII I/O, command line and Python.
 */
inline TypeOneDRule getRuleString(std::string const &name){
    try{
        return getStringRuleMap().at(name);
    }catch(std::out_of_range &){
        return rule_none;
    }
}

/*!
 * \ingroup TasmanianIO
 * \brief Map the enumerate to a string, used in ASCII I/O, command line and Python.
 */
inline std::string getRuleString(TypeOneDRule rule){
    auto smap = getStringRuleMap();
    return std::find_if(smap.begin(), smap.end(),
                        [&](std::pair<std::string, TypeOneDRule> r)->bool{ return (r.second == rule); })->first;
}

/*!
 * \ingroup TasmanianIO
 * \brief Creates a map with \b int (used by Fortran and binary I/O) mapped to \b TypeOneDRule enums.
 */
inline std::vector<TypeOneDRule> getIntRuleMap(){
    return {rule_none, rule_clenshawcurtis, rule_clenshawcurtis0,
        rule_chebyshev, rule_chebyshevodd, rule_gausslegendre, rule_gausslegendreodd,
        rule_gausspatterson, rule_leja, rule_lejaodd,
        rule_rleja, rule_rlejadouble2, rule_rlejadouble4, rule_rlejaodd,
        rule_rlejashifted, rule_rlejashiftedeven, rule_rlejashifteddouble,
        rule_maxlebesgue, rule_maxlebesgueodd, rule_minlebesgue, rule_minlebesgueodd,
        rule_mindelta, rule_mindeltaodd, rule_gausschebyshev1, rule_gausschebyshev1odd,
        rule_gausschebyshev2, rule_gausschebyshev2odd, rule_fejer2,
        rule_gaussgegenbauer, rule_gaussgegenbauerodd, rule_gaussjacobi, rule_gaussjacobiodd,
        rule_gausslaguerre, rule_gausslaguerreodd, rule_gausshermite, rule_gausshermiteodd,
        rule_customtabulated,
        rule_localp, rule_localp0, rule_semilocalp,
        rule_wavelet, rule_fourier, rule_localpb};
}

/*!
 * \ingroup TasmanianIO
 * \brief Map the int rule index to the enumerate, used in Fortran and binary IO.
 */
inline TypeOneDRule getRuleInt(int r){
    auto rmap = getIntRuleMap();
    return ((size_t) r < rmap.size()) ? rmap[(size_t) r] : rule_none;
}

/*!
 * \ingroup TasmanianIO
 * \brief Map the enumerate to an int, used in Fortran and binary IO.
 */
inline int getRuleInt(TypeOneDRule rule){
    auto rmap = getIntRuleMap();
    return (int) std::distance(rmap.begin(), std::find_if(rmap.begin(), rmap.end(),
                               [&](TypeOneDRule r)->bool{ return (r == rule); }));
}

/*!
 * \ingroup TasmanianIO
 * \brief Creates a map with \b std::string rule names (used by C/Python/CLI) mapped to \b TypeDepth enums.
 */
inline std::map<std::string, TypeDepth> getStringToDepthMap(){
    return std::initializer_list<std::pair<std::string const, TypeDepth>>{
        {"level",        type_level},
        {"curved",       type_curved},
        {"iptotal",      type_iptotal},
        {"ipcurved",     type_ipcurved},
        {"qptotal",      type_qptotal},
        {"qpcurved",     type_qpcurved},
        {"hyperbolic",   type_hyperbolic},
        {"iphyperbolic", type_iphyperbolic},
        {"qphyperbolic", type_qphyperbolic},
        {"tensor",       type_tensor},
        {"iptensor",     type_iptensor},
        {"qptensor",     type_qptensor}};
}
/*!
 * \ingroup TasmanianIO
 * \brief Map the string to the enumerate multi-index selection strategy, used in command line and Python.
 */
inline TypeDepth getDepthTypeString(std::string const &name){
    try{
        return getStringToDepthMap().at(name);
    }catch(std::out_of_range &){
        return type_none;
    }
}

/*!
 * \ingroup TasmanianIO
 * \brief Map the integer to the enumerate multi-index selection strategy, used in Fortran.
 */
inline TypeDepth getDepthTypeInt(int t){
    std::vector<TypeDepth> imap = {type_none, type_level, type_curved, type_iptotal,
        type_ipcurved, type_qptotal, type_qpcurved, type_hyperbolic, type_iphyperbolic,
        type_qphyperbolic, type_tensor, type_iptensor, type_qptensor};
    return ((size_t) t < imap.size()) ? imap[(size_t) t] : type_none;
}

/*!
 * \ingroup TasmanianIO
 * \brief Creates a map with \b std::string rule names (used by C/Python/CLI) mapped to \b TypeRefinement enums.
 */
inline std::map<std::string, TypeRefinement> getStringToRefinementMap(){
    return std::initializer_list<std::pair<std::string const, TypeRefinement>>{
        {"classic",   refine_classic},
        {"parents",   refine_parents_first},
        {"direction", refine_direction_selective},
        {"fds",       refine_fds},
        {"stable",    refine_stable}};
}
/*!
 * \ingroup TasmanianIO
 * \brief Map the string to the enumerate hierarchical refinement strategy, used in command line and Python.
 */
inline TypeRefinement getTypeRefinementString(std::string const &name){
    try{
        return getStringToRefinementMap().at(name);
    }catch(std::out_of_range &){
        return refine_none;
    }
}
/*!
 * \ingroup TasmanianIO
 * \brief Map the integer to the enumerate hierarchical refinement strategy, used by Fortran.
 */
inline TypeRefinement getTypeRefinementInt(int refinement){
    std::vector<TypeRefinement> imap = {refine_none, refine_classic,
                                        refine_parents_first, refine_direction_selective,
                                        refine_fds, refine_stable};
    return ((size_t) refinement < imap.size()) ? imap[(size_t) refinement] : refine_none;
}

/*!
 * \ingroup TasmanianIO
 * \brief Write the flag to file, ascii uses 0 and 1, binary uses characters y and n (counter intuitive, I know).
 */
template<bool iomode, IOPad pad>
void writeFlag(bool flag, std::ostream &os){
    if (iomode == mode_ascii){
        os << ((flag) ? "1" : "0");
        if ((pad == pad_rspace) || ((pad == pad_auto) && flag)) os << " ";
        if ((pad == pad_line) || ((pad == pad_auto) && !flag)) os << std::endl;
    }else{
        char cflag = ((flag) ? 'y' : 'n');
        os.write(&cflag, sizeof(char));
    }
}

/*!
 * \ingroup TasmanianIO
 * \brief Read a flag, ascii uses 0 and 1, binary uses characters y and n (counter intuitive, I know).
 */
template<typename iomode>
bool readFlag(std::istream &os){
    if (std::is_same<iomode, mode_ascii_type>::value){
        int flag;
        os >> flag;
        return (flag != 0);
    }else{
        char cflag;
        os.read(&cflag, sizeof(char));
        return (cflag == 'y');
    }
}

/*!
 * \ingroup TasmanianIO
 * \brief Write the vector to the stream, the vector cannot be empty.
 */
template<bool iomode, IOPad pad, typename VecType>
void writeVector(const std::vector<VecType> &x, std::ostream &os){
    if (iomode == mode_ascii){
        if (pad == pad_lspace)
            for(auto i : x) os << " " << i;
        if (pad == pad_rspace)
            for(auto i : x) os << i << " ";
        if ((pad == pad_none) || (pad == pad_line)){
            os << x[0];
            for(size_t i = 1; i < x.size(); i++) os << " " << x[i];
            if (pad == pad_line) os << std::endl;
        }
    }else{
        os.write((char*) x.data(), x.size() * sizeof(VecType));
    }
}

/*!
 * \ingroup TasmanianIO
 * \brief Read the vector from the stream, the size must already be set.
 */
template<typename iomode, typename VecType>
void readVector(std::istream &is, std::vector<VecType> &x){
    if (std::is_same<iomode, mode_ascii_type>::value){
        for(auto &i : x) is >> i;
    }else{
        is.read((char*) x.data(), x.size() * sizeof(VecType));
    }
}

/*!
 * \ingroup TasmanianIO
 * \brief Read the vector with the specified size.
 */
template<typename iomode, typename VecType, typename SizeType>
std::vector<VecType> readVector(std::istream &is, SizeType num_entries){
    std::vector<VecType> x((size_t) num_entries);
    readVector<iomode, VecType>(is, x);
    return x;
}

/*!
 * \ingroup TasmanianIO
 * \brief Write a bunch of numbers with the same type.
 */
template<bool iomode, IOPad pad, typename... Vals>
void writeNumbers(std::ostream &os, Vals... vals){
    std::vector<typename std::tuple_element<0, std::tuple<Vals...>>::type> values = {vals...};
    writeVector<iomode, pad>(values, os);
}

/*!
 * \ingroup TasmanianIO
 * \brief Read a single number, used to read ints (and potentially cast to size_t) or read a double.
 */
template<typename iomode, typename Val>
Val readNumber(std::istream &is){
    Val v;
    if (std::is_same<iomode, mode_ascii_type>::value){
        is >> v;
    }else{
        is.read((char*) &v, sizeof(Val));
    }
    return v;
}

/*!
 * \ingroup TasmanianIO
 * \brief Write a rule.
 */
template<bool iomode>
void writeRule(TypeOneDRule rule, std::ostream &os){
    if (iomode == mode_ascii){
        os << getRuleString(rule) << std::endl;
    }else{
        int r = getRuleInt(rule);
        os.write((char*) &r, sizeof(int));
    }
}

/*!
 * \ingroup TasmanianIO
 * \brief Read a rule.
 */
template<typename iomode>
TypeOneDRule readRule(std::istream &is){
    if (std::is_same<iomode, mode_ascii_type>::value){
        std::string T;
        is >> T;
        return getRuleString(T);
    }else{
        return getRuleInt(readNumber<mode_binary_type, int>(is));
    }
}

}

}

#endif
