/*
 * Copyright (c) 2022, Miroslav Stoyanov & Weiwei Kong
 *
 * This file is part of
 * Toolkit for Adaptive Stochastic Modeling And Non-Intrusive ApproximatioN: TASMANIAN
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following
 * conditions are met:
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
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * UT-BATTELLE, LLC AND THE UNITED STATES GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL WARRANTIES, BOTH EXPRESSED AND
 * IMPLIED. THERE ARE NO EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF
 * THE SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL
 * ACCOMPLISH THE INTENDED RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE. THE USER ASSUMES
 * RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS, CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING
 * FROM OR ARISING OUT OF, IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.
 */

#ifndef __TASMANIAN_OPTIM_ENUMERATES_HPP
#define __TASMANIAN_OPTIM_ENUMERATES_HPP

#include "TasmanianDREAM.hpp"

/*!
 * \internal
 * \file tsgOptimizationUtils.hpp
 * \brief Utility functions and aliases in the optimization module.
 * \author Weiwei Kong & Miroslav Stoyanov
 * \ingroup TasmanianOptimization
 *
 * Defines functions and type aliases that are used in the Tasmanian Optimization module. The file is included in every other
 * TasOptimization header.
 * \endinternal
 */

/*!
 * \ingroup TasmanianOptimization
 * \addtogroup OptimizationUtil Miscellaneous utility functions and aliases
 *
 * Several type aliases and utility functions based on similar ones in the DREAM module.
 */

namespace TasOptimization {

/*! \internal
 * \ingroup OptimizationUtil
 *
 * Checks if a variable size \b var_name associated with \b var_name inside \b method_name matches an expected size \b exp_size.
 * If it does not match, a runtime error is thrown.
 * \endinternal
 */
inline void checkVarSize(const std::string method_name, const std::string var_name, const int var_size, const int exp_size) {
    if (var_size != exp_size) {
        throw std::runtime_error("Size of " + var_name + " (" + std::to_string(var_size) + ") in the function " + method_name +
                                 "() is not equal to its expected value of (" + std::to_string(exp_size) + ")");
    }
}

// Functions used in optimization.
using GenericBatchedFunction = std::function<void(const std::vector<double> &x_batch, std::vector<double> &y_batch)>;

/*! \ingroup OptimizationUtil
 * \brief Generic non-batched objective function signature.
 *
 * Accepts a single input \b x and returns the evaluation of the function on the point \b x.
 */
using ObjectiveFunctionSingle = std::function<double(const std::vector<double> &x)>;

/*! \ingroup OptimizationUtil
 * \brief Generic batched objective function signature.
 *
 * Batched version of TasOptimization::ObjectiveFunctionSingle.
 * Accepts multiple points \b x_batch and writes their corresponding values into
 * \b fval_batch. It is expected that the size of \b x_batch is a multiple of the size of the output.
 */
using ObjectiveFunction = GenericBatchedFunction;

/*! \ingroup OptimizationUtil
 * \brief Creates a TasOptimization::ObjectiveFunction object from a TasOptimization::ObjectiveFunctionSingle object.
 *
 * Given a TasOptimization::ObjectiveFunctionSingle \b f_single and the size of its input \b num_dimensions,
 * returns a TasOptimization::ObjectiveFunction that evaluates
 * a batch of points \f$ x_1,\ldots,x_k \f$ to \f$ {\rm f\_single}(x_1),\ldots, {\rm f\_single}(x_k) \f$.
 */
inline ObjectiveFunction makeObjectiveFunction(const int num_dimensions, const ObjectiveFunctionSingle f_single) {
    return [=](const std::vector<double> &x_values, std::vector<double> &fval_values)->void {
        int num_points = x_values.size() / num_dimensions;
        std::vector<double> x(num_dimensions);
        for (int i=0; i<num_points; i++) {
            std::copy_n(x_values.begin() + i * num_dimensions, num_dimensions, x.begin());
            fval_values[i] = f_single(x);
        }
    };
}

/*! \ingroup OptimizationUtil
 * \brief Generic non-batched gradient function signature.
 *
 * Accepts a single input \b x and returns the gradient at the point \b x.
 */
using GradientFunctionSingle = std::function<std::vector<double>(const std::vector<double> &x_batch)>;

/*! \ingroup OptimizationUtil
 * \brief Generic non-batched gradient function signature.
 *
 * Batched version of TasOptimization::GradientFunctionSingle.
 * Accepts multiple points \b x_batch and writes their corresponding gradients into
 * \b grad_batch. It is expected that the size of the output is a multiple of the size of \b x_batch.
 */
using GradientFunction = GenericBatchedFunction;

/*! \ingroup OptimizationUtil
 * \brief Creates a TasOptimization::GradientFunction object from a TasOptimization::GradientFunctionSingle object.
 *
 * Given a TasOptimization::GradientFunctionSingle \b grad_single, the size of the domain of the function \b num_dimensions,
 * and the size of the codomain of the function \b num_outputs, returns a TasOptimization::GradientFunction that evaluates
 * a batch of points \f$ x_1,\ldots,x_k \f$ to \f$ {\rm grad\_single}(x_1),\ldots, {\rm grad\_single}(x_k) \f$.
 */
inline GradientFunction makeGradientFunction(const int num_dimensions, const int num_outputs, const GradientFunctionSingle grad_single) {
    return [=](const std::vector<double> &x_values, std::vector<double> &grad_values)->void {
        int num_points = x_values.size() / num_dimensions;
        std::vector<double> x(num_dimensions), jacobian(num_dimensions * num_outputs);
        for (int i=0; i<num_points; i++) {
            std::copy_n(x_values.begin() + i * num_dimensions, num_dimensions, x.begin());
            jacobian = grad_single(x);
            std::copy_n(jacobian.begin(), num_dimensions * num_outputs, grad_values.begin() + i * num_dimensions * num_outputs);
        }
    };
}

/*! \ingroup OptimizationUtil
 * \brief Generic non-batched projection function signature.
 *
 * Accepts a single input \b x and returns the projection of \b x onto a user-specified domain.
 */
using ProjectionFunctionSingle = GradientFunctionSingle; // Same function prototype.

/*! \ingroup OptimizationUtil
 * \brief Generic non-batched projection function signature.
 *
 * Batched version of TasOptimization::ProjectionFunctionSingle.
 * Accepts multiple points \b x_batch and writes their corresponding projections into
 * \b grad_batch. It is expected that the size of the output is equal the size of \b x_batch.
 */
using ProjectionFunction = GenericBatchedFunction;

/*! \ingroup OptimizationUtil
 * \brief Creates a TasOptimization::ProjectionFunction object from a TasOptimization::ProjectionFunctionSingle object.
 *
 * Given a TasOptimization::GradientFunctionSingle \b grad_single, the size of the domain of the function \b num_dimensions,
 * returns a TasOptimization::GradientFunction that evaluates a batch of points \f$ x_1,\ldots,x_k \f$ to \f$
 * {\rm proj\_single}(x_1),\ldots, {\rm proj\_single}(x_k) \f$.
 */
inline ProjectionFunction makeProjectionFunction(const int num_dimensions, const ProjectionFunctionSingle proj_single) {
    // Since ProjectionFunctionSingle and ProjectionFunction are the same GradientFunctionSingle and GradientFunction, we can
    // just call the existing batch making function for gradients.
    return makeGradientFunction(num_dimensions, 1, proj_single);
}


} // End namespace

#endif
