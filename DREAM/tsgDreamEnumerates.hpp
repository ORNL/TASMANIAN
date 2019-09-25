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

#ifndef __TASMANIAN_DREAM_ENUMERATES_HPP
#define __TASMANIAN_DREAM_ENUMERATES_HPP

#include <random>

#include "TasmanianSparseGrid.hpp"

/*!
 * \internal
 * \file tsgDreamEnumerates.hpp
 * \brief The enumerated types used in the DREAM module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianDREAM
 *
 * Defines the enumerated types used throughout the DREAM module.
 * The file is included in every other DREAM header.
 * \endinternal
 */

/*!
 * \ingroup TasmanianDREAM
 * \addtogroup DREAMEnumerates Enumerated types
 *
 * The Enumerate types used in the external and internal API of the DREAM module.
 */

namespace TasDREAM{

/*!
 * \internal
 * \ingroup TasmanianDREAM
 * \addtogroup DREAMUtils Miscellaneous utility templates, borrowed from TasGrid.
 *
 * Several templates borrowed from the Sparse Grids module and included here for easier usage.
 * \endinternal
 */

/*!
 * \internal
 * \ingroup DREAMUtils
 *
 * Similar to TasGrid::Utils this holds aliases to size_mult(), TasGrid::Utils::Wrapper2D, and copyArray().
 * \endinternal
 */
namespace Utils{
    //! \brief See TasGrid::Utils::size_mult().
    using TasGrid::Utils::size_mult;
    //! \brief See TasGrid::Utils::Wrapper2D.
    using TasGrid::Utils::Wrapper2D;
    //! \brief See TasGrid::Utils::copyArray().
    using TasGrid::Utils::copyArray;
}

//! \brief Describes whether sampling should be done with the regular or logarithm form of the probability density.
//! \ingroup DREAMEnumerates

//! Probability distributions with very localized support will have values close to zero over most of the sampling domain.
//! Sampling using the logarithm of the probability density can be more stable, when operations with very small numbers cause problems.
//! In most cases, it is up to the user to provide the appropriate values, but the sampling algorithm needs to know which form is being used.
enum TypeSamplingForm{
    //! \brief Use the standard form for the probability density.
    regform,

    //! \brief Use the logarithm form for the probability density.
    logform
};

//! \brief Indicates a specific probability distribution for the associated function.
//! \ingroup DREAMEnumerates

//! Used to instantiate the \b getDensity() variadric template. See the template documentation for the specific formulas.
enum TypeDistribution{
    //! \brief Uniform distribution.
    dist_uniform,

    //! \brief Gaussian or Normal distribution defined by mean and variance.
    dist_gaussian,

    //! \brief Exponential distribution (i.e., special case of the Gamma distribution).
    dist_exponential,

    //! \brief Beta distribution, corresponds to Gauss-Jacobi sparse grid rule \b TasGrid::rule_gaussjacobi.
    dist_beta,

    //! \brief Gamma distribution, corresponds to Gauss-Laguerre sparse grid rule \b TasGrid::rule_gausslaguerre.
    dist_gamma,

    //! \brief Indicates a no-distribution (no correction).
    dist_none,

    //! \brief Indicates I/O error or unspecified distribution.
    dist_null
};

/*!
 * \internal
 * \ingroup DREAMUtils
 * \brief Similar to TasGrid::IO this holds conversion methods between enums and int/string types.
 *
 * \endinternal
 */
namespace IO{
    /*!
     * \internal
     * \ingroup TasmanianIO
     * \brief Converts an integer to TypeSamplingForm, synced with the Python interface.
     *
     * \endinternal
     */
    inline TypeSamplingForm intToForm(int form){ return (form == 0) ? regform : logform; }

    /*!
     * \internal
     * \ingroup TasmanianIO
     * \brief Returns the map from a string to a distribution type.
     *
     * \endinternal
     */
    inline std::map<std::string, TypeDistribution> getStringRuleMap(){
        return std::initializer_list<std::pair<std::string const, TypeDistribution>>{
            {"null",        dist_null},
            {"none",        dist_none},
            {"uniform",     dist_uniform},
            {"gaussian",    dist_gaussian},
            {"exponential", dist_exponential},
            {"beta",        dist_beta},
            {"gamma",       dist_gamma}};
    }

    /*!
     * \internal
    * \ingroup TasmanianIO
    * \brief Map the string distribution name to the enumerate, used by python.
    *
    * \endinternal
    */
    inline TypeDistribution getDistributionString(std::string const &name){
        try{
            return getStringRuleMap().at(name);
        }catch(std::out_of_range &){
            return dist_null;
        }
    }
}

namespace DreamMaths{
/*!
 * \internal
 * \ingroup DREAMEnumerates
 * \brief Dream copy of \b TasGrid::Maths::pi.
 *
 * \endinternal
 */
constexpr double pi = TasGrid::Maths::pi;
}

}

#endif
