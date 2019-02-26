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

#include <tuple>

#include "tsgCoreOneDimensional.hpp"

//! \internal
//! \file tsgIOHelpers.hpp
//! \brief Templates to simply file I/O.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianIO
//!
//! Several templates that simplify the I/O of Tasmanian.
//! Commonly used operations are lumped into templates with simple binary/ascii switch.

/*!
 * \internal
 * \ingroup TasmanianSG
 * \addtogroup TasmanianIO Templates for common I/O methods, simple binary/ascii switch provided for all templates
 * \endinternal
 */

namespace TasGrid{

namespace IO{

//! \internal
//! \brief Indicate the type of padding to use, none, space, new-line, etc.
//! \ingroup TasmanianIO
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

//! \internal
//! \brief Write the flag to file, ascii uses 0 and 1, binary uses characters y and n (counter intuitive, I know).
//! \ingroup TasmanianIO
template<bool useAscii, IOPad pad>
void writeFlag(bool flag, std::ostream &os){
    if (useAscii){
        os << ((flag) ? "1" : "0");
        if ((pad == pad_rspace) || ((pad == pad_auto) && flag)) os << " ";
        if ((pad == pad_line) || ((pad == pad_auto) && !flag)) os << std::endl;
    }else{
        char cflag = ((flag) ? 'y' : 'n');
        os.write(&cflag, sizeof(char));
    }
}

//! \internal
//! \brief Read a flag, ascii uses 0 and 1, binary uses characters y and n (counter intuitive, I know).
//! \ingroup TasmanianIO
template<bool useAscii>
bool readFlag(std::istream &os){
    if (useAscii){
        int flag;
        os >> flag;
        return (flag == 1);
    }else{
        char cflag;
        os.read(&cflag, sizeof(char));
        return (cflag == 'y');
    }
}

//! \internal
//! \brief Write the vector to the stream.
//! \ingroup TasmanianIO
template<bool useAscii, IOPad pad, typename VecType>
void writeVector(const std::vector<VecType> &x, std::ostream &os){
    if (useAscii){
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

//! \internal
//! \brief Read the vector from the stream.
//! \ingroup TasmanianIO
template<bool useAscii, typename VecType>
void readVector(std::istream &os, std::vector<VecType> &x){
    if (useAscii){
        for(auto &i : x) os >> i;
    }else{
        os.read((char*) x.data(), x.size() * sizeof(VecType));
    }
}

//! \internal
//! \brief Write a bunch of numbers with the same type.
//! \ingroup TasmanianIO
template<bool useAscii, IOPad pad, typename... Vals>
void writeNumbers(std::ostream &os, Vals... vals){
    std::vector<typename std::tuple_element<0, std::tuple<Vals...>>::type> values = {vals...};
    writeVector<useAscii, pad>(values, os);
}

//! \internal
//! \brief Read a single number, used to read ints (and potentially cast to size_t) or read a double.
//! \ingroup TasmanianIO
template<bool useAscii, typename Val>
Val readNumber(std::istream &os){
    Val v;
    if (useAscii){
        os >> v;
    }else{
        os.read((char*) &v, sizeof(Val));
    }
    return v;
}

//! \internal
//! \brief Write a rule.
//! \ingroup TasmanianIO
template<bool useAscii>
void writeRule(TypeOneDRule rule, std::ostream &os){
    if (useAscii){
        os << OneDimensionalMeta::getIORuleString(rule) << std::endl;
    }else{
        int r = OneDimensionalMeta::getIORuleInt(rule);
        os.write((char*) &r, sizeof(int));
    }
}

//! \internal
//! \brief Read a rule.
//! \ingroup TasmanianIO
template<bool useAscii>
TypeOneDRule readRule(std::istream &is){
    if (useAscii){
        std::string T;
        is >> T;
        return OneDimensionalMeta::getIORuleString(T.c_str());
    }else{
        return OneDimensionalMeta::getIORuleInt(readNumber<false, int>(is));
    }
}

}

}

#endif
