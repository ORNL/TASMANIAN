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

#ifndef __TASMANIAN_SPARSE_GRID_SEQUENCE_OPTIMIZER_CPP
#define __TASMANIAN_SPARSE_GRID_SEQUENCE_OPTIMIZER_CPP

#include "tsgSequenceOptimizer.hpp"

namespace TasGrid{

namespace Optimizer{
/*!
 * \ingroup TasmanianSequenceOpt
 * \brief Computes the coefficients needed for fast evaluation of the Lagrange polynomials.
 */
std::vector<double> makeCoefficients(std::vector<double> const &nodes){
    size_t num_nodes = nodes.size();
    std::vector<double> coeffs(num_nodes);
    for(size_t i=0; i<num_nodes; i++){
        double c = 1.0;
        for(size_t j=0; j<i; j++){
            c *= (nodes[i] - nodes[j]);
        }
        for(size_t j=i+1; j<num_nodes; j++){
            c *= (nodes[i] - nodes[j]);
        }
        coeffs[i] = c;
    }
    return coeffs;
}
/*!
 * \ingroup TasmanianSequenceOpt
 * \brief Computes the values of the Lagrange polynomials at \b x, used in most functionals.
 */
std::vector<double> evalLagrange(std::vector<double> const &nodes, std::vector<double> const &coeffs, double x){
    int num_nodes = (int) nodes.size();
    std::vector<double> lag(nodes.size());
    lag[0] = 1.0;
    for(int i=0; i<num_nodes-1; i++){
        lag[i+1] = (x - nodes[i]) * lag[i];
    }
    double w = 1.0;
    lag[num_nodes-1] /= coeffs[num_nodes-1];
    for(int i= num_nodes-2; i>=0; i--){
        w *= (x - nodes[i+1]);
        lag[i] *= w / coeffs[i];
    }
    return lag;
}
/*!
 * \ingroup TasmanianSequenceOpt
 * \brief Computes the derivative of the \b inode basis functions at \b x.
 */
double differentiateBasis(std::vector<double> const &nodes, std::vector<double> const &coeffs, size_t inode, double x){
    size_t num_nodes = nodes.size();
    double s = 1.0;
    double p = 1.0;
    double n = (inode != 0) ? (x - nodes[0]) : (x - nodes[1]);

    for(size_t j=1; j<inode; j++){
        p *= n;
        n = (x - nodes[j]);
        s *= n;
        s += p;
    }
    for(size_t j = ((inode == 0) ? 2 : inode+1); j<num_nodes; j++){
        p *= n;
        n = (x - nodes[j]);
        s *= n;
        s += p;
    }
    return s / coeffs[inode];
}

/*!
 * \ingroup TasmanianSequenceOpt
 * \brief Data needed for the functional associated with the sequence rule, specialized for each case.
 *
 * In most cases, we are working with Lagrange polynomials and associated coefficients.
 */
template<TypeOneDRule> struct CurrentNodes{
    //! \brief Default constructor, retain a copy of the nodes and compute the coefficients.
    CurrentNodes(std::vector<double> const &cnodes)
        : nodes(cnodes), coeff(makeCoefficients(cnodes)){}
    //! \brief Constructor that combines the \b cnodes with the \b new_node.
    CurrentNodes(std::vector<double> const &cnodes, double new_node)
            : nodes(cnodes){
        nodes.push_back(new_node);
        coeff = makeCoefficients(nodes);
    }
    //! \brief Current set of nodes.
    std::vector<double> nodes;
    //! \brief Coefficients cache.
    std::vector<double> coeff;
};
/*!
 * \ingroup TasmanianSequenceOpt
 * \brief Specialization for \b rule_leja, no need for coefficients.
 */
template<> struct CurrentNodes<rule_leja>{
    //! \brief Retain a copy of the nodes.
    CurrentNodes(std::vector<double> const &cnodes) : nodes(cnodes){}
    //! \brief Current set of nodes.
    std::vector<double> nodes;
};
/*!
 * \ingroup TasmanianSequenceOpt
 * \brief Specialization for \b rule_mindeltaodd, requires two levels.
 */
template<> struct CurrentNodes<rule_mindeltaodd>{
    //! \brief Construct two levels using the \b cnodes and \b cnodes with added \b new_node.
    CurrentNodes(std::vector<double> const &cnodes, double new_node)
            : nodes(cnodes), nodes_less1(cnodes), coeff_less1(makeCoefficients(cnodes)){
        nodes.push_back(new_node);
        coeff = makeCoefficients(nodes);
    }
    //! \brief Nodes for the current level.
    std::vector<double> nodes;
    //! \brief Nodes for the previous level.
    std::vector<double> nodes_less1;
    //! \brief Coefficients cache current level.
    std::vector<double> coeff;
    //! \brief Coefficients cache for previous level.
    std::vector<double> coeff_less1;
};

/*!
 * \ingroup TasmanianSequenceOpt
 * \brief Indicates whether a \b rule has associated derivative, most do.
 */
template<TypeOneDRule rule> struct HasDerivative{
    //! \brief Indicates whether the functional is differentiable.
    static constexpr bool value = true;
};
/*!
 * \ingroup TasmanianSequenceOpt
 * \brief Specialization for \b rule_minlebesgue which uses a min-max problem and cannot be differentiated.
 */
template<> struct HasDerivative<rule_minlebesgue>{
    //! \brief Indicates non-differentiable.
    static constexpr bool value = false;
};
/*!
 * \ingroup TasmanianSequenceOpt
 * \brief Specialization for \b rule_mindelta which uses a min-max problem and cannot be differentiated.
 */
template<> struct HasDerivative<rule_mindelta>{
    //! \brief Indicates non-differentiable.
    static constexpr bool value = false;
};

/*!
 * \ingroup TasmanianSequenceOpt
 * \brief Computes the value of the functional for given \b x, specialized for each sequence.
 */
template<TypeOneDRule rule> double getValue(CurrentNodes<rule> const&, double){ return 0.0; }

/*!
 * \ingroup TasmanianSequenceOpt
 * \brief Computes the derivative of the functional for given \b x, specialized for each sequence with \b HasDerivative<rule>::value \b = \b true.
 */
template<TypeOneDRule rule> double getDerivative(CurrentNodes<rule> const&, double){ return 0.0; }

/*!
 * \ingroup TasmanianSequenceOpt
 * \brief Uses the secant method to find local maximum of the  functional of the current nodes, uses \b left and \b right as starting points.
 *
 * Uses the secant method to find the zero of the derivative of the functional associated with the \b rule.
 * The method converges fast but may converge to the wrong answer if the initial guess is not close to the zero.
 * When called from \b computeLocalMaximum(), the result of a coarse pattern search is used as initial guess.
 */
template<TypeOneDRule rule>
OptimizerResult performSecantSearch(CurrentNodes<rule> const& current, double left, double right){
    auto func = [&](double x)->OptimizerResult{ return {x, getDerivative<rule>(current, x)}; };

    OptimizerResult past    = func(left);
    OptimizerResult present = func(right);

    if (std::abs(past.value) < std::abs(present.value))
        std::swap(past, present);

    int iterations = 0;
    while((std::abs(present.value) > 3.0 * Maths::num_tol) && (iterations < 1000)){ // 1000 guards against stagnation
        OptimizerResult future = func(present.node - present.value * (present.node - past.node) / (present.value - past.value));
        past    = present;
        present = future;
        iterations++;
    }

    return present;
}

/*!
 * \ingroup TasmanianSequenceOpt
 * \brief Finds the maximum of the functional of the current nodes in the interval between \b left_node and \b right_node.
 *
 * Uses pattern search with simple left, middle, and right points.
 * If the functional of the \b rule is differentiable, then the pattern search is used as an initial guess
 * to secant optimization method.
 */
template<TypeOneDRule rule>
OptimizerResult computeLocalMaximum(CurrentNodes<rule> const& current, double left_node, double right_node){
    auto func = [&](double x)->OptimizerResult{ return {x, getValue<rule>(current, x)}; };

    double pattern = 0.5 * (right_node - left_node); // pattern width
    OptimizerResult left   = func(left_node);
    OptimizerResult middle = func(left_node + pattern);
    OptimizerResult right  = func(right_node);

    // if differentiable, use coarse tolerance since we will post-process.
    double tolerance = (HasDerivative<rule>::value) ? Maths::num_tol * 1.E-3 : Maths::num_tol;

    while(pattern > tolerance){
        if (middle.value >= std::max(left.value, right.value)){ // middle is largest, shrink
            pattern /= 2.0;
            left  = func(middle.node - pattern);
            right = func(middle.node + pattern);
        }else if (left.value >= std::max(middle.value, right.value)){ // left is the largest
            if (left.node - pattern < left_node){ // if going out of bounds
                pattern /= 2.0;
                right  = middle;
                middle = func(left.node + pattern);
            }else{ // shift left
                right  = middle;
                middle = left;
                left  = func(middle.node - pattern);
            }
        }else{ // right must be the largest
            if (right.node + pattern > right_node){ // if going out of bounds
                pattern /= 2.0;
                left   = middle;
                middle = func(right.node - pattern);
            }else{ // shift right
                left   = middle;
                middle = right;
                right  = func(middle.node + pattern);
            }
        }
    }

    if (HasDerivative<rule>::value){
        middle = performSecantSearch(current, left.node, right.node); // this returns the derivative
        middle = func(middle.node);
    }

    return middle;
}

/*!
 * \ingroup TasmanianSequenceOpt
 * \brief Finds the maximum of the functional over the interval (-1, 1).
 *
 * Given the \b current set of nodes, construct the functional for the given \b rule
 * and perform a local optimization on every interval, i.e., compute the local maximum
 * between every two adjacent nodes in \b current.
 * The work in done in parallel and the global maximum is reported as the largest among
 * the local maximums.
 */
template<TypeOneDRule rule> OptimizerResult computeMaximum(CurrentNodes<rule> const& current){
    std::vector<double> sorted = current.nodes;
    std::sort(sorted.begin(), sorted.end());
    int num_intervals = (int) sorted.size() - 1;

    auto func = [&](double x)->OptimizerResult{ return {x, getValue<rule>(current, x)}; };

    OptimizerResult max_result = func(-1.0);

    OptimizerResult right_result = func(1.0);

    if (right_result.value > max_result.value)
        max_result = right_result;

    #pragma omp parallel
    {
        OptimizerResult thread_max = max_result, thread_result;

        #pragma omp for schedule(dynamic)
        for(int i=0; i<num_intervals; i++){
            thread_result = computeLocalMaximum(current, sorted[i], sorted[i+1]);
            if ((!std::isnan(thread_result.value)) && (thread_result.value > thread_max.value))
                thread_max = thread_result;
        }

        #pragma omp critical
        {
            if (thread_max.value > max_result.value){
                max_result = thread_max;
            }
        }
    }

    return max_result;
}

/*!
 * \ingroup TasmanianSequenceOpt
 * \brief The \b rule_leja functional.
 */
template<> double getValue<rule_leja>(CurrentNodes<rule_leja> const& current, double x){
    double p = 1.0;
    for(auto n : current.nodes) p *= (x - n);
    return std::abs(p);
}
/*!
 * \ingroup TasmanianSequenceOpt
 * \brief The \b rule_maxlebesgue functional.
 */
template<> double getValue<rule_maxlebesgue>(CurrentNodes<rule_maxlebesgue> const& current, double x){
    std::vector<double> lag = evalLagrange(current.nodes, current.coeff, x);
    double sum = 0.0;
    for(auto l : lag) sum += std::abs(l);
    return sum;
}
/*!
 * \ingroup TasmanianSequenceOpt
 * \brief The \b rule_mindeltaodd functional (indicates companion to \b rule_mindelta for the min-max problem).
 */
template<> double getValue<rule_mindeltaodd>(CurrentNodes<rule_mindeltaodd> const& current, double x){
    auto lag       = evalLagrange(current.nodes, current.coeff, x);
    auto lag_less1 = evalLagrange(current.nodes_less1, current.coeff_less1, x);

    return std::abs(lag.back()) + std::inner_product(lag_less1.begin(), lag_less1.end(), lag.begin(), 0.0,
                                     std::plus<double>(), [](double a, double b)->double{ return std::abs(a - b); });
}

/*!
 * \ingroup TasmanianSequenceOpt
 * \brief The \b rule_minlebesgue functional, uses the \b rule_maxlebesgue functions in min-max problem.
 */
template<> double getValue<rule_minlebesgue>(CurrentNodes<rule_minlebesgue> const& current, double x){
    for(auto n : current.nodes) if (std::abs(x - n) < 10 * Maths::num_tol) return -1.E+100;

    CurrentNodes<rule_maxlebesgue> companion(current.nodes, x);
    return - computeMaximum(companion).value;
}

/*!
 * \ingroup TasmanianSequenceOpt
 * \brief The \b rule_mindelta functional, uses the \b rule_mindeltaodd functions in min-max problem.
 */
template<> double getValue<rule_mindelta>(CurrentNodes<rule_mindelta> const& current, double x){
    for(auto n : current.nodes) if (std::abs(x - n) < 10 * Maths::num_tol) return -1.E+100;

    CurrentNodes<rule_mindeltaodd> companion(current.nodes, x);
    return - computeMaximum(companion).value;
}

/*!
 * \ingroup TasmanianSequenceOpt
 * \brief The derivative of the \b rule_leja functional.
 */
template<> double getDerivative<rule_leja>(CurrentNodes<rule_leja> const& current, double x){
    double s = 1.0, p = 1.0, n = (x - current.nodes[0]);
    for(size_t j=1; j<current.nodes.size(); j++){
        p *= n;
        n = (x - current.nodes[j]);
        s *= n;
        s += p;
    }
    return s;
}

/*!
 * \ingroup TasmanianSequenceOpt
 * \brief The derivative of the \b rule_maxlebesgue functional.
 */
template<> double getDerivative<rule_maxlebesgue>(CurrentNodes<rule_maxlebesgue> const& current, double x){
    std::vector<double> lag = evalLagrange(current.nodes, current.coeff, x);
    double sum = 0.0;
    for(size_t i=0; i<lag.size(); i++)
        sum += Maths::sign(lag[i]) * differentiateBasis(current.nodes, current.coeff, i, x);
    return sum;
}

/*!
 * \ingroup TasmanianSequenceOpt
 * \brief The derivative of the \b rule_mindeltaodd functional.
 */
template<> double getDerivative<rule_mindeltaodd>(CurrentNodes<rule_mindeltaodd> const& current, double x){
    auto lag       = evalLagrange(current.nodes,       current.coeff,       x);
    auto lag_less1 = evalLagrange(current.nodes_less1, current.coeff_less1, x);

    double sum = 0.0;
    for(size_t i=0; i<lag_less1.size(); i++){
        sum += Maths::sign(lag[i] - lag_less1[i])
                * (differentiateBasis(current.nodes,        current.coeff,       i, x)
                   - differentiateBasis(current.nodes_less1, current.coeff_less1, i, x));
    }
    return sum + Maths::sign(lag.back()) * differentiateBasis(current.nodes, current.coeff, lag.size() - 1, x);
}

std::vector<double> getPrecomputedMinLebesgueNodes(){
    return        { 0.00000000000000000e+00,
                    1.00000000000000000e+00,
                   -1.00000000000000000e+00,
                    5.77350269189625731e-01,
                   -6.82983516382591915e-01,
                    3.08973641709742008e-01,
                   -8.88438098101029028e-01,
                    9.01204348257380272e-01,
                   -3.79117423735989223e-01,
                    7.68752974570125480e-01,
                   -5.49318186832010835e-01,
                   -9.64959809152743819e-01,
                    9.67682365886019302e-01,
                   -1.28290369854639041e-01,
                    4.45625686400429988e-01,
                   -8.04184165478832425e-01,
                    1.67074539986499709e-01,
                    8.44556006683382599e-01,
                   -4.67815101214832718e-01,
                   -9.88223645152207397e-01,
                    6.75111201618315282e-01,
                   -2.31781319689011223e-01,
                    9.90353577025220644e-01,
                   -8.50249643305543756e-01,
                    3.78518581322005110e-01,
                   -7.34179511891844605e-01,
                    9.37348529794296281e-01,
                    7.43693755708555865e-02,
                   -9.43140751144033618e-01,
                    7.22894448368867515e-01,
                   -3.13966521558236233e-01,
                    5.21332302672803283e-01,
                   -6.22914216661596742e-01,
                    9.81207150951301732e-01,
                   -9.96112540596109541e-01,
                    2.27115897487506796e-01,
                    8.71343953281246475e-01,
                   -5.83881096023192936e-01,
                   -7.07087869617637754e-02,
                   -9.21222090369298474e-01,
                    6.35275009690072778e-01,
                    9.97425667442937813e-01,
                   -7.74028611711955805e-01,
                    1.19609717340386237e-01,
                    8.06079215959057738e-01,
                   -9.77828412798477653e-01,
                   -2.81002869364477770e-01,
                    4.81831205470734381e-01,
                   -4.33236587065744694e-01,
                    9.51233414543255273e-01};
}

std::vector<double> getPrecomputedMinDeltaNodes(){
    return        { 0.00000000000000000e+00,
                    1.00000000000000000e+00,
                   -1.00000000000000000e+00,
                    5.85786437626666157e-01,
                   -7.11391303688287735e-01,
                   -3.65515299787630255e-01,
                    8.73364116971659943e-01,
                   -9.09934766576546594e-01,
                    3.48149066530255846e-01,
                    7.61022503813320816e-01,
                   -5.63728220989981099e-01,
                    1.61287729090120513e-01,
                   -9.70195195938610255e-01,
                    9.71760208817352256e-01,
                   -2.12010912937840634e-01,
                   -8.19500788126849566e-01,
                    9.28317867167342214e-01,
                    4.63681477197408431e-01,
                   -4.75500758437281013e-01,
                    6.83991900744136849e-01,
                   -9.89565716162813747e-01,
                   -9.47327399278441035e-02,
                    9.91020313545615483e-01,
                   -7.70894179555761561e-01,
                    2.57418053836716398e-01,
                   -9.38159793780884654e-01,
                    8.23375672488699584e-01,
                   -2.90168112727829441e-01,
                    5.28028188455664571e-01,
                   -6.43769278853295157e-01,
                    9.51854409585821237e-01,
                   -8.71963189117939352e-01,
                    8.85485426816430832e-02,
                    7.23563624423755547e-01,
                   -9.96440136176172664e-01,
                   -1.50233510276863047e-01,
                    9.96883045450636107e-01,
                   -5.22343767211049914e-01,
                    4.02868136231700480e-01,
                   -8.46787710296258989e-01,
                    8.96281851550158049e-01,
                   -4.11624139689389490e-01,
                    6.30422054421239331e-01,
                   -9.57396180946457842e-01,
                    2.09670340519686721e-01,
                    8.47208170728777743e-01,
                   -6.80096440009749670e-01,
                   -4.07635282203260146e-02,
                    9.81787150881257009e-01,
                   -9.81715548483999001e-01};
}

inline std::vector<double> getPrecomputed(TypeOneDRule rule){
    if (rule == rule_leja){
        return {0.0, 1.0, -1.0, std::sqrt(1.0/3.0)};
    }else if (rule == rule_maxlebesgue){
        return {0.0, 1.0, -1.0, 0.5};
    }else if (rule == rule_minlebesgue){
        return getPrecomputedMinLebesgueNodes();
    }else{ // rule_mindelta
        return getPrecomputedMinDeltaNodes();
    }
}

template<TypeOneDRule rule> double getNextNode(std::vector<double> const &nodes){
    return computeMaximum(CurrentNodes<rule>(nodes)).node;
}

template double getNextNode<rule_leja>(std::vector<double> const &nodes);
template double getNextNode<rule_maxlebesgue>(std::vector<double> const &nodes);
template double getNextNode<rule_minlebesgue>(std::vector<double> const &nodes);
template double getNextNode<rule_mindelta>(std::vector<double> const &nodes);

template<TypeOneDRule rule>
std::vector<double> getGreedyNodes(int n){
    // load the first few precomputed nodes
    auto precomputed = getPrecomputed(rule);
    size_t usefirst = std::min(precomputed.size(), (size_t) n);
    std::vector<double> nodes(precomputed.begin(), precomputed.begin() + usefirst);
    if (n > (int) precomputed.size()){
        nodes.reserve((size_t) n);
        for(int i = (int) precomputed.size(); i<n; i++)
            nodes.push_back(getNextNode<rule>(nodes));
    }

    return nodes;
}

template std::vector<double> getGreedyNodes<rule_leja>(int n);
template std::vector<double> getGreedyNodes<rule_maxlebesgue>(int n);
template std::vector<double> getGreedyNodes<rule_minlebesgue>(int n);
template std::vector<double> getGreedyNodes<rule_mindelta>(int n);

}

}

#endif
