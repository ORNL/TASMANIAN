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

#ifndef __TSG_ONE_DIMENSIONAL_WRAPPER_HPP
#define __TSG_ONE_DIMENSIONAL_WRAPPER_HPP

#include "tsgHardCodedTabulatedRules.hpp"

/*!
 * \internal
 * \file tsgOneDimensionalWrapper.hpp
 * \brief Cache for one dimensional rules.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianCoreOneDimensional
 *
 * A class to cache one dimensional rules, nodes, weight, etc.
 * \endinternal
 */

namespace TasGrid{

/*!
 * \ingroup TasmanianCoreOneDimensional
 * \brief A class to cache one dimensional rules, nodes, weight, meta-data, etc.
 */
class OneDimensionalWrapper{
public:
    //! \brief Default constructor, create an emptry cache.
    OneDimensionalWrapper() : isNonNested(false), num_levels(0), rule(rule_none){}
    //! \brief Default destructor.
    ~OneDimensionalWrapper() = default;

    /*!
     * \brief Load the cache.
     *
     * Load \b max_level of cache for a rule given by \b crule. If the rule is not
     * custom tabulated, then \b custom is not used (could be empty).
     * Similarly, \b alpha and \b beta are used only by rules that use the corresponding
     * transformation parameters.
     */
    OneDimensionalWrapper(const CustomTabulated &custom, int max_level, TypeOneDRule crule, double alpha, double beta) :
        isNonNested(OneDimensionalMeta::isNonNested(crule)),
        num_levels(max_level + 1),
        rule(crule),
        num_points(num_levels),
        pntr(num_levels + 1, 0), // the first entry must be set to zero
        weights(num_levels),
        nodes(num_levels),
        coeff(num_levels)
    {
        if (crule == rule_customtabulated){
            if (max_level + 1 > custom.getNumLevels()){
                std::string message = "ERROR: custom-tabulated rule needed with levels ";
                message += std::to_string(num_levels);
                message += ", but only ";
                message += std::to_string(custom.getNumLevels());
                message += " are provided.";
                throw std::runtime_error(message);
            }
        }
        // find the points per level and the cumulative pointers
        if (rule != rule_customtabulated){
            for(int l=0; l<num_levels; l++){
                num_points[l] = OneDimensionalMeta::getNumPoints(l, rule);
                pntr[l+1] = pntr[l] + num_points[l];
            }
        }else{
            for(int l=0; l<num_levels; l++){
                num_points[l] = custom.getNumPoints(l);
                pntr[l+1] = pntr[l] + num_points[l];
            }
        }
        int num_total = pntr[num_levels];

        if (isNonNested){
            indx.reserve(num_total);
            unique.reserve(num_total);

            for(int l=0; l<num_levels; l++){
                int n = num_points[l];
                if ((rule == rule_chebyshev) || (rule == rule_chebyshevodd)){
                    OneDimensionalNodes::getChebyshev(n, weights[l], nodes[l]);
                }else if ((rule == rule_gausslegendre) || (rule == rule_gausslegendreodd)){
                    OneDimensionalNodes::getGaussLegendre(n, weights[l], nodes[l]);
                }else if ((rule == rule_gausschebyshev1) || (rule == rule_gausschebyshev1odd)){
                    OneDimensionalNodes::getGaussChebyshev1(n, weights[l], nodes[l]);
                }else if ((rule == rule_gausschebyshev2) || (rule == rule_gausschebyshev2odd)){
                    OneDimensionalNodes::getGaussChebyshev2(n, weights[l], nodes[l]);
                }else if ((rule == rule_gaussgegenbauer) || (rule == rule_gaussgegenbauerodd)){
                    OneDimensionalNodes::getGaussJacobi(n, weights[l], nodes[l], alpha, alpha);
                }else if ((rule == rule_gausshermite) || (rule == rule_gausshermiteodd)){
                    OneDimensionalNodes::getGaussHermite(n, weights[l], nodes[l], alpha);
                }else if ((rule == rule_gaussjacobi) || (rule == rule_gaussjacobiodd)){
                    OneDimensionalNodes::getGaussJacobi(n, weights[l], nodes[l], alpha, beta);
                }else if ((rule == rule_gausslaguerre) || (rule == rule_gausslaguerreodd)){
                    OneDimensionalNodes::getGaussLaguerre(n, weights[l], nodes[l], alpha);
                }else if (rule == rule_customtabulated){
                    custom.getWeightsNodes(l, weights[l], nodes[l]);
                }

                for(auto x : nodes[l]){
                    int point = -1;
                    for(int j=0; j<(int) unique.size(); j++){
                        if (std::abs(x - unique[j]) < Maths::num_tol){
                            point = j;
                            break;
                        }
                    }
                    if (point == - 1){ // new point found
                        unique.push_back(x);
                        point = (int) (unique.size() - 1);
                    }
                    indx.push_back(point);
                }
            }

            unique.shrink_to_fit();

            // make coefficients
            for(int l=0; l<num_levels; l++){
                int n = num_points[l];
                coeff[l].resize(n);
                double *x = nodes[l].data();
                for(int i=0; i<n; i++){
                    double c = 1.0;
                    for(int j=0; j<i; j++){
                        c /= (x[i] - x[j]);
                    }
                    for(int j=i+1; j<n; j++){
                        c /= (x[i] - x[j]);
                    }
                    coeff[l][i] = c;
                }
            }
        }else{
            if (rule == rule_clenshawcurtis){
                unique = OneDimensionalNodes::getClenshawCurtisNodes(max_level);
            }else if (rule == rule_clenshawcurtis0){
                unique = OneDimensionalNodes::getClenshawCurtisNodesZero(max_level);
            }else if (rule == rule_fejer2){
                unique = OneDimensionalNodes::getFejer2Nodes(max_level);
            }else if (rule == rule_gausspatterson){
                TableGaussPatterson gp;
                if (num_levels > gp.getNumLevels()){
                    std::string message = "ERROR: gauss-patterson rule needed with level ";
                    message += std::to_string(max_level);
                    message += ", but only ";
                    message += std::to_string(gp.getNumLevels());
                    message += " are hardcoded.";
                    throw std::runtime_error(message);
                }
                unique = gp.getNodes(max_level);

                // move here to avoid a warning
                for(int l=0; l<num_levels; l++){
                    weights[l].resize(num_points[l]);
                    for(int i=0; i<num_points[l]; i++){
                        weights[l][i] = gp.getWeight(l, i);
                    }
                }
            }else if (rule == rule_rleja){
                unique = OneDimensionalNodes::getRLeja(OneDimensionalMeta::getNumPoints(max_level,rule));
            }else if ((rule == rule_rlejaodd) || (rule == rule_rlejadouble2) || (rule == rule_rlejadouble4)){
                unique = OneDimensionalNodes::getRLejaCentered(OneDimensionalMeta::getNumPoints(max_level,rule));
            }else if ((rule == rule_rlejashifted) || (rule == rule_rlejashiftedeven) || (rule == rule_rlejashifteddouble)){
                unique = OneDimensionalNodes::getRLejaShifted(OneDimensionalMeta::getNumPoints(max_level,rule));
            }else if ((rule == rule_leja) || (rule == rule_lejaodd)){
                unique = Optimizer::getGreedyNodes<rule_leja>(OneDimensionalMeta::getNumPoints(max_level, rule));
            }else if ((rule == rule_maxlebesgue) || (rule == rule_maxlebesgueodd)){
                unique = Optimizer::getGreedyNodes<rule_maxlebesgue>(OneDimensionalMeta::getNumPoints(max_level, rule));
            }else if ((rule == rule_minlebesgue) || (rule == rule_minlebesgueodd)){
                unique = Optimizer::getGreedyNodes<rule_minlebesgue>(OneDimensionalMeta::getNumPoints(max_level, rule));
            }else if ((rule == rule_mindelta) || (rule == rule_mindeltaodd)){
                unique = Optimizer::getGreedyNodes<rule_mindelta>(OneDimensionalMeta::getNumPoints(max_level, rule));
            }else{ // if (rule==rule_fourier)
                unique = OneDimensionalNodes::getFourierNodes(max_level);
            }

            for(int l=0; l<num_levels; l++){
                int n = num_points[l];
                weights[l].resize(n);
                coeff[l].resize(n);
                for(int i=0; i<n; i++){
                    double c = 1.0;
                    for(int j=0; j<i; j++){
                        c /= (unique[i] - unique[j]);
                    }
                    for(int j=i+1; j<n; j++){
                        c /= (unique[i] - unique[j]);
                    }
                    if (rule == rule_clenshawcurtis0){
                        c /= (unique[i] - 1.0) * (unique[i] + 1.0);
                    }
                    coeff[l][i] = c;
                }
            }

            if (rule == rule_clenshawcurtis){
                for(int l=0; l<num_levels; l++){
                    for(int i=0; i<num_points[l]; i++){
                        weights[l][i] = OneDimensionalNodes::getClenshawCurtisWeight(l, i);
                    }
                }
            }else if (rule == rule_clenshawcurtis0){
                for(int l=0; l<num_levels; l++){
                    for(int i=0; i<num_points[l]; i++){
                        weights[l][i] = OneDimensionalNodes::getClenshawCurtisWeightZero(l, i);
                    }
                }
            }else if (rule == rule_fejer2){
                for(int l=0; l<num_levels; l++){
                    for(int i=0; i<num_points[l]; i++){
                        weights[l][i] = OneDimensionalNodes::getFejer2Weight(l, i);
                    }
                }
            }else if ((rule == rule_rleja) || (rule == rule_rlejaodd) || (rule == rule_rlejadouble2) || (rule == rule_rlejadouble4) || (rule == rule_leja) || (rule == rule_lejaodd) ||
                    (rule == rule_rlejashifted) || (rule == rule_rlejashiftedeven) || (rule == rule_rlejashifteddouble) ||
                (rule == rule_maxlebesgue) || (rule == rule_maxlebesgueodd) || (rule == rule_minlebesgue) || (rule == rule_minlebesgueodd) || (rule == rule_mindelta) || (rule == rule_mindeltaodd)){
                // if missing the quadrature weights
                int n = 1 + num_points[max_level] / 2; // number of Gauss-Legendre points needed to integrate the basis functions
                std::vector<double> lag_x, lag_w;
                OneDimensionalNodes::getGaussLegendre(n, lag_w, lag_x);
                for(int l=0; l<num_levels; l++){
                    std::fill(weights[l].begin(), weights[l].end(), 0.0);
                    int npl = num_points[l];
                    std::vector<double> v(npl);
                    for(int i=0; i<n; i++){
                        v[0] = 1.0;
                        for(int j=0; j<npl-1; j++){
                            v[j+1] = (lag_x[i] - unique[j]) * v[j];
                        }
                        v[npl-1] *= coeff[l][npl-1];
                        double s = 1.0;
                        for(int j=npl-2; j>=0; j--){
                            s *= (lag_x[i] - unique[j+1]);
                            v[j] *= s * coeff[l][j];
                        }
                        for(int j=0; j<npl; j++){
                            weights[l][j] += lag_w[i] * v[j];
                        }
                    }
                }
            }
        }
    }

    //! \brief Overload that skips the custom rule altogether, more convenient in Fourier grids.
    OneDimensionalWrapper(int max_level, TypeOneDRule crule, double alpha, double beta) :
        OneDimensionalWrapper(CustomTabulated(), max_level, crule, alpha, beta){}

    //! \brief Get the number of points for the \b level, can only querry loaded levels.
    int getNumPoints(int level) const{ return num_points[level]; }
    //! \brief Get the global index of point \b j on \b level (used only by non-nested rules).
    int getPointIndex(int level, int j) const{ return indx[pntr[level] + j]; }

    //! \brief Get number of loaded nodes.
    int getNumNodes() const{ return (int) unique.size(); }
    //! \brief Get the canonical coordinate of the node with global index \b j.
    double getNode(int j) const{ return unique[j]; }
    //! \brief Get the quadrature weight of the \b j-th node on the \b level (for non-nested rules, using index local to the level)
    double getWeight(int level, int j) const{ return weights[level][j]; }

    //! \brief Get an array of all nodes associated with the \b level (the array is ordered by local index).
    const double* getNodes(int level) const{ return (isNonNested) ? nodes[level].data() : unique.data(); }
    //! \brief Get an array of the Lagrange coefficients associated with the \b level, the order matches \b getNodes().
    const double* getCoefficients(int level) const{ return coeff[level].data(); }

    //! \brief Get a reference to the offsets of all point counts used by the \b CacheLagrange class.
    const std::vector<int>& getPointsCount() const{ return pntr; }

    //! \brief Get the loaded rule.
    TypeOneDRule getType() const{ return rule; }
    //! \brief Get the number of cached levels.
    int getNumLevels() const{ return num_levels; }

    //! \brief Return reference to all unique nodes.
    const std::vector<double>& getUnique() const{ return unique; }
    //! \brief Return all nodes concatenated per level.
    std::vector<double> getAllNodes() const{ return Utils::mergeVectors(nodes); }
    //! \brief Return all coefficients concatenated per level.
    std::vector<double> getAllCoeff() const{ return Utils::mergeVectors(coeff); }
    //! \brief Return reference to the number of nodes per level.
    const std::vector<int>& getNumNodesPerLevel() const{ return num_points; }
    //! \brief Return cumulative points per level.
    std::vector<int> getOffsetNodesPerLevel() const{
        std::vector<int> offsets(num_points.size());
        offsets[0] = 0;
        for(size_t i=1; i<num_points.size(); i++)
            offsets[i] = offsets[i-1] + num_points[i-1];
        return offsets;
    }

private:
    bool isNonNested;
    int num_levels;
    TypeOneDRule rule;

    std::vector<int> num_points; // give the number of points per level
    std::vector<int> pntr; // gives the cumulative offset of each level
    std::vector<int> indx; // gives a list of the points associated with each level

    std::vector<std::vector<double>> weights; // contains the weight associated with each level
    std::vector<std::vector<double>> nodes; // contains all nodes for each level
    std::vector<double> unique; // contains the x-coordinate of each sample point

    std::vector<std::vector<double>> coeff; // the coefficients of the Lagrange
};

}

#endif
