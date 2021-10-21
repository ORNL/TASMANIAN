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

#ifndef __TASGRID_COMMON_HPP
#define __TASGRID_COMMON_HPP

/*!
 * \internal
 * \file tasgridCLICommon.hpp
 * \brief Common executable includes and templates.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianCLI
 *
 * Defines common includes for the various executables and templates for managing command line arguments.
 * \endinternal
 */

#include <random>
#include <deque>
#include <cctype>

#include "TasmanianSparseGrid.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using namespace TasGrid;

/*!
 * \internal
 * \ingroup TasmanianCLI
 * \brief Creates a std::deque of strings from the CLI arguments, skips the first argument.
 *
 * Converts the CLI arguments to strings, the first argument is the name of the executable and it is omitted.
 * \endinternal
 */
inline std::deque<std::string> stringArgs(int argc, const char** argv){
    std::deque<std::string> args;
    for(int i = 1; i < argc; i++)
        args.push_back(std::string(argv[i]));
    return args;
}

/*!
 * \internal
 * \ingroup TasmanianCLI
 * \brief Returns \b true if the string contains a sub-string with the word "help" (case insensitive).
 *
 * \endinternal
 */
inline bool hasHelp(std::string const &arg){
    std::string lower(arg.size(), ' ');
    std::transform(arg.begin(), arg.end(), lower.begin(),
        [](char c)->char{
            return static_cast<char>(std::tolower(static_cast<int>(c)));
        });

    auto pos = lower.find("help");
    if ((pos < lower.size()) && (pos  + 4 <= lower.size()))
        return (lower.substr(pos, pos + 4).compare("help") == 0);

    return false;
}

/*!
 * \internal
 * \ingroup TasmanianCLI
 * \brief Returns \b true if the string contains a request for version information.
 *
 * Accepted strings are "-v", "version", "verbose", and "info" with "-" or "--".
 * \endinternal
 */
inline bool hasInfo(std::string const &s){
    std::map<std::string, bool> accpetable = {
        {"-v",   true},
        {"version", true},
        {"-version", true},
        {"--version",   true},
        {"verbose",   true},
        {"-verbose",   true},
        {"--verbose",   true},
        {"info",  true},
        {"-info",  true},
        {"--info",  true},
    };

    try{
        return accpetable.at(s);
    }catch(std::out_of_range &){
        return false;
    }
}

/*!
 * \internal
 * \ingroup TasmanianCLI
 * \brief Returns \b true if the string contains "random", "-random", "rand", or "-rand"
 *
 * \endinternal
 */
inline bool hasRandom(std::string const &s){
    return ((s == "random") || (s == "-random") || (s == "rand") || (s == "-rand"));
}

/*!
 * \internal
 * \ingroup TasmanianCLI
 * \brief Returns \b true if the string contains "gpu", "-gpu", "gpuid", or "-gpuid"
 *
 * \endinternal
 */
inline bool hasGpuID(std::string const &s){
    return ((s == "gpu") or (s == "-gpu") or (s == "gpuid") or (s == "-gpuid"));
}

/*!
 * \ingroup TasmanianCLI
 * \brief Returns the specified GPU id (could be -1) or throws
 */
inline int getGpuID(std::deque<std::string> const &args){
    int gpuid = (args.empty()) ? -2 : std::stoi(args.front());
    if ((gpuid < -1) || (gpuid >= TasmanianSparseGrid::getNumGPUs())){
        cerr << "ERROR: -gpuid requires a valid gpuid!" << endl;
        cerr << "      see ./tasgrid -v for a list of detected GPUs." << endl;
        throw std::invalid_argument("Invalid GPU ID!");
    }
    return gpuid;
}


#endif
