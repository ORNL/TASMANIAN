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

#ifndef __TSG_RULE_LOCAL_POLYNOMIAL_HPP
#define __TSG_RULE_LOCAL_POLYNOMIAL_HPP

#include "tsgCoreOneDimensional.hpp"

namespace TasGrid{

#ifndef __TASMANIAN_DOXYGEN_SKIP
class BaseRuleLocalPolynomial{
public:
    BaseRuleLocalPolynomial() : max_order(0){}
    virtual ~BaseRuleLocalPolynomial() = default;

    virtual int getMaxOrder() const = 0;
    virtual void setMaxOrder(int order) = 0;

    virtual TypeOneDRule getType() const = 0;

    virtual int getNumPoints(int level) const = 0; // for building the initial grid

    virtual double getNode(int point) const = 0;

    virtual int getLevel(int point) const = 0;
    virtual double getSupport(int point) const = 0;

    virtual int getMaxNumKids() const = 0;
    virtual int getMaxNumParents() const = 0;

    virtual int getParent(int point) const = 0;
    virtual int getStepParent(int point) const = 0;
    virtual int getKid(int point, int kid_number) const = 0;

    virtual double evalRaw(int point, double x) const = 0; // normalizes x (i.e., (x-node) / support), but it does not check the support
    virtual double evalSupport(int point, double x, bool &isSupported) const = 0; // // normalizes x (i.e., (x-node) / support) and checks if x is within the support

    // Analogous to the eval functions, but differentiates at x instead.
    virtual double diffRaw(int point, double x) const = 0;
    virtual double diffSupport(int point, double x, bool &isSupported) const = 0;

    virtual double getArea(int point, std::vector<double> const &w, std::vector<double> const &x) const = 0;
    // integrate the function associated with the point, constant to cubic are known analytically, higher order need a 1-D quadrature rule

    // constructs a Vandemond matrix for the first num_rows points
    virtual void van_matrix(int num_rows, std::vector<int> &pntr, std::vector<int> &indx, std::vector<double> &vals) = 0;

protected:
    int max_order;
};


template<TypeOneDRule rule, bool isZeroOrder>
class templRuleLocalPolynomial : public BaseRuleLocalPolynomial{
public:
    templRuleLocalPolynomial() = default;
    ~templRuleLocalPolynomial() = default;

    int getMaxOrder() const override{ return max_order; }
    void setMaxOrder(int order) override{ max_order = order; }

    TypeOneDRule getType() const override{ return rule; }

    int getNumPoints(int level) const override{
        if (isZeroOrder){
            int n = 1;
            while (level-- > 0) n *= 3;
            return n;
        }else{
            if ((rule == rule_localp) || (rule == rule_semilocalp)){
                return (level == 0) ? 1 : ((1 << level) + 1);
            }else if (rule == rule_localp0){
                return (1 << (level+1)) -1;
            }else{ // rule == rule_localpb
                return ((1 << level) + 1);
            }
        }
    }

    double getNode(int point) const override{
        if (isZeroOrder){
            return -2.0 + (1.0 / ((double) Maths::int3log3(point))) * (3*point + 2 - point % 2);
        }else{
            if ((rule == rule_localp) || (rule == rule_semilocalp)){
                if (point == 0) return  0.0;
                if (point == 1) return -1.0;
                if (point == 2) return  1.0;
                return ((double)(2*point - 1)) / ((double) Maths::int2log2(point - 1)) - 3.0;
            }else if (rule == rule_localp0){
                return ((double)(2*point +3) ) / ((double) Maths::int2log2(point + 1) ) - 3.0;
            }else{ // rule == rule_localpb
                if (point == 0) return -1.0;
                if (point == 1) return  1.0;
                if (point == 2) return  0.0;
                return ((double)(2*point - 1)) / ((double) Maths::int2log2(point - 1)) - 3.0;
            }
        }
    }

    int getLevel(int point) const override{
        if (isZeroOrder){
            int level = 0;
            while(point >= 1){ point /= 3; level += 1; }
            return level;
        }else{
            if ((rule == rule_localp) || (rule == rule_semilocalp)){
                return (point == 0) ? 0 : (point == 1) ? 1 : (Maths::intlog2(point - 1) + 1);
            }else if (rule == rule_localp0){
                return Maths::intlog2(point + 1);
            }else{ // rule == rule_localpb
                return ((point == 0) || (point == 1)) ? 0 : (Maths::intlog2(point - 1) + 1);
            }
        }
    }
    double getSupport(int point) const override{
        if (isZeroOrder){
            return 1.0 / (double) Maths::int3log3(point);
        }else{
            if ((rule == rule_localp) || (rule == rule_semilocalp)){
                return (point == 0) ? 1.0 : 1.0 / ((double) Maths::int2log2(point - 1));
            }else if (rule == rule_localp0){
                return 1.0 / ((double) Maths::int2log2(point + 1));
            }else{ // rule == rule_localpb
                return ((point == 0) || (point == 1)) ? 2.0 : 1.0 / ((double) Maths::int2log2(point - 1));
            }
        }
    }

    int getMaxNumKids() const override{ return (isZeroOrder) ? 4 : 2; }
    int getMaxNumParents() const override{ return ((isZeroOrder || (rule == rule_semilocalp) || (rule == rule_localpb)) ? 2 : 1); }

    int getParent(int point) const override{
        if (isZeroOrder){
            return (point == 0) ? -1 : point / 3;
        }else{
            if ((rule == rule_localp) || (rule == rule_semilocalp)){
                int dad = (point + 1) / 2;
                if (point < 4) dad--;
                return dad;
            }else if (rule == rule_localp0){
                if (point == 0) return -1;
                return (point - 1) / 2;
            }else{ // rule == rule_localpb
                return (point < 2) ? -1 : ((point + 1) / 2);
            }
        }
    }

    int getStepParent(int point) const override{
        if (isZeroOrder){
            int i3l3 = Maths::int3log3(point);
            if (point == i3l3/3) return -1;
            if (point == i3l3-1) return -1;
            if ((point % 3 == 2) && (point % 2 == 0)) return point / 3 + 1;
            if ((point % 3 == 0) && (point % 2 == 1)) return point / 3 - 1;
            return -1;
        }else{
            if (rule == rule_semilocalp){
                if (point == 3) return 2;
                if (point == 4) return 1;
            }else if (rule == rule_localpb){
                if (point == 2) return 0;
            }
            return -1;
        }
    }

    int getKid(int point, int kid_number) const override{
        if (isZeroOrder){
            if (point == 0) return (kid_number == 0) ? 1 : (kid_number==1) ? 2 : -1;
            if (kid_number == 3){
                int i3l3 = Maths::int3log3(point);
                if (point == i3l3/3) return -1;
                if (point == i3l3-1) return -1;
                return (point % 2 == 0) ? 3*point + 3 : 3*point - 1;
            }
            return 3*point + kid_number;
        }else if ((rule == rule_localp) || (rule == rule_semilocalp)){
            if (kid_number == 0){
                if (point == 0) return 1;
                if (point == 1) return 3;
                if (point == 2) return 4;
                return 2*point-1;
            }else{
                if (point == 0) return 2;
                if ((point == 1) || (point == 2)) return -1;
                return 2*point;
            }
        }else if (rule == rule_localp0){
            return 2*point + ((kid_number == 0) ? 1 : 2);
        }else{ // localpb
            if ((point == 0) || (point == 1)) return (kid_number == 0) ? 2 : -1;
            return 2*point - ((kid_number == 0) ? 1 : 0);
        }
    }

    double evalRaw(int point, double x) const override{
        if (isZeroOrder){
            if (std::abs(x - getNode(point)) > getSupport(point)) return 0.0; // maybe can speed this up
            return 1.0;
        }else{
            if ((rule == rule_localp) || (rule == rule_semilocalp)){
                if (point == 0) return 1.0;
                if (rule == rule_semilocalp){
                    if (point == 1) return 0.5 * x * (x - 1.0);
                    if (point == 2) return 0.5 * x * (x + 1.0);
                }
            }
            double xn = scaleX(point, x);
            if (rule != rule_semilocalp) if (max_order == 1) return 1.0 - std::abs(xn);
            if (max_order == 2) return evalPWQuadratic(point, xn);
            if (max_order == 3) return evalPWCubic(point, xn);
            return evalPWPower(point, xn);
        }
    }

    double evalSupport(int point, double x, bool &isSupported) const override{
        if (isZeroOrder){
            double distance = std::abs(x - getNode(point)), support = getSupport(point); // can speed this up, see above
            isSupported = (distance <= (2.0) * support);
            return (distance > support) ? 0.0 : 1.0;
        }else{
            isSupported = true;
            if ((rule == rule_localp) || (rule == rule_semilocalp)){
                if (point == 0) return 1.0;
                if (rule == rule_semilocalp){
                    if (point == 1) return 0.5 * x * (x - 1.0);
                    if (point == 2) return 0.5 * x * (x + 1.0);
                }
            }
            double xn = scaleX(point, x);
            if (std::abs(xn) <= 1.0){
                if (rule != rule_semilocalp) if (max_order == 1) return 1.0 - std::abs(xn);
                if (max_order == 2) return evalPWQuadratic(point, xn);
                if (max_order == 3) return evalPWCubic(point, xn);
                return evalPWPower(point, xn);
            }else{
                isSupported = false;
                return 0.0;
            }
        }
    }

    double diffRaw(int point, double x) const override{
        if (isZeroOrder) {
            return 0.0;
        } else {
            if (rule == rule_localp or rule == rule_semilocalp) {
                if (point == 0) return 0.0;
                if (rule == rule_semilocalp) {
                    if (point == 1) return x - 0.5;
                    if (point == 2) return x + 0.5;
                }
            }
            double xn = scaleX(point, x);
            double an = scaleDiffX(point);
            if (rule != rule_semilocalp and max_order == 1) {
                if (x == 1.0) {
                    if (rule == rule_localp0 and point == 0)
                        return  -1.0;
                    else if (rule == rule_localpb and point == 1)
                        return 0.5;
                    else if (rule == rule_localpb and point == 2)
                        return -1.0;
                    else if (point == 2)
                        return an;
                }
                return (xn >= 0 ? -1.0 : 1.0) * an;
            }
            if (max_order == 2) return an * diffPWQuadratic(point, xn);
            if (max_order == 3) return an * diffPWCubic(point, xn);
            return an * diffPWPower(point, xn);
        }
    }

    double diffSupport(int point, double x, bool &isSupported) const override{
        if (isZeroOrder) {
            isSupported = false;
            return 0.0;
        } else {
            // We still need (isSupported == true) at point = 0 so that TasGrid::GridLocalPolynomial::walkTree<3>()
            // walks through the children of this point.
           isSupported = true;
            if (rule == rule_localp or rule == rule_semilocalp) {
                if (point == 0) {
                    return 0.0;
                }
                if (rule == rule_semilocalp) {
                    if (point == 1) return x - 0.5;
                    if (point == 2) return x + 0.5;
                }
            }
            // Note that the support will be of the form [a, b), except for x = +1.0, i.e., the rightmost node in the domain. This
            // is to avoid double-counting derivatives at points of discontinuity.
            double xn = scaleX(point, x);
            isSupported = (-1.0 <= xn and xn < 1.0) or (x == 1.0 and xn == 1.0);
            if (isSupported) {
                double an = scaleDiffX(point);
                if (rule != rule_semilocalp and max_order == 1) {
                    // Edge case due to the logic of avoiding double-counting.
                    if (x == 1.0 and rule == rule_localp and point == 2) return an;
                    return (xn >= 0 ? -1.0 : 1.0) * an;
                }
                if (max_order == 2) return an * diffPWQuadratic(point, xn);
                if (max_order == 3) return an * diffPWCubic(point, xn);
                return an * diffPWPower(point, xn);
            } else {
                return 0.0;
            }
        }
    }

    double getArea(int point, std::vector<double> const &w, std::vector<double> const &x) const override{
        if (isZeroOrder){
            return 2.0 * getSupport(point);
        }else{
            if (rule == rule_localp){
                if (point == 0) return 2.0;
                if ((point == 1) || (point == 2)) return 0.5;
                if (max_order == 1) return getSupport(point);
                if ((max_order == 2) || (max_order == 3) || (point <= 8))  return (4.0/3.0) * getSupport(point);
            }else if (rule == rule_semilocalp){
                if (point == 0) return 2.0;
                if ((point == 1) || (point == 2)) return 1.0/3.0;
                if ((max_order == 2) || (max_order == 3) || (point <= 4))  return (4.0/3.0) * getSupport(point);
            }else if (rule == rule_localpb){
                if ((point == 0) || (point == 1)) return 1.0;
                if (max_order == 1) return getSupport(point);
                if ((max_order == 2) || (max_order == 3) || (point <= 4))  return (4.0/3.0) * getSupport(point);
            }else if (rule == rule_localp0){
                if (max_order == 1) return getSupport(point);
                if ((max_order == 2) || (max_order == 3) || (point <= 2))  return (4.0/3.0) * getSupport(point);
            }
            double sum = 0.0;
            for(size_t i=0; i<w.size(); i++) sum += w[i] * evalPWPower(point, x[i]);
            return sum * getSupport(point);
        }
    }

    void van_matrix(int num_rows, std::vector<int> &pntr, std::vector<int> &indx, std::vector<double> &vals) override{
        int max_level = getLevel(num_rows);

        if (isZeroOrder) {
            switch(num_rows) {
                case 0:
                    indx = {0,};
                    vals = {1.0,};
                    pntr = {0, 1};
                    break;
                case 1:
                    indx = {0, 0, 1};
                    vals = {1.0, 1.0, 1.0};
                    pntr = {0, 1, 3};
                    break;
                case 2:
                    indx = {0, 0, 1, 0, 2};
                    vals = {1.0, 1.0, 1.0, 1.0, 1.0};
                    pntr = {0, 1, 3, 5};
                    break;
                default:
                    indx.clear();
                    indx.reserve(num_rows * max_level);
                    vals.clear();
                    vals.reserve(num_rows * max_level);
                    pntr = std::vector<int>(num_rows + 1, 0);
                    for(auto i : std::array<int, 5>{0, 0, 1, 0, 2})
                        indx.push_back(i);
                    for(size_t i=0; i<5; i++)
                        vals.push_back(1.0);
                    pntr[1] = 1;
                    pntr[2] = 3;
                    {
                        std::vector<int> ancestors;
                        ancestors.reserve(max_level);
                        for(int r=3; r<num_rows; r++) {
                            pntr[r] = static_cast<int>(indx.size());
                            ancestors.clear();
                            int kid = r;
                            int dad = kid / 3;
                            while(dad != 0) {
                                if (dad % 2 == 0) {
                                    if (kid != 3 * dad)
                                        ancestors.push_back(dad);
                                } else {
                                    if (kid != 3 * dad + 2)
                                        ancestors.push_back(dad);
                                }
                                kid = dad;
                                dad = kid / 3;
                            }
                            indx.push_back(0);
                            indx.insert(indx.end(), ancestors.rbegin(), ancestors.rend());
                            indx.push_back(r);
                            for(size_t i=0; i<ancestors.size() + 2; i++)
                                vals.push_back(1.0);
                        }
                        pntr.back() = static_cast<int>(indx.size());
                    }
                    break;
            };
        } else {
            if (rule == rule_localp) {
                switch(num_rows) {
                    case 0:
                        indx = {0,};
                        vals = {1.0,};
                        pntr = {0, 1};
                        break;
                    default:
                        indx.clear();
                        indx.reserve(num_rows * max_level);
                        vals.clear();
                        vals.reserve(num_rows * max_level);
                        pntr = std::vector<int>(num_rows + 1, 0);
                        indx.push_back(0);
                        vals.push_back(1.0);
                        {
                            std::vector<int> ancestors;
                            std::vector<double> ancestors_vals;
                            ancestors.reserve(max_level);
                            ancestors_vals.reserve(max_level);
                            for(int r=1; r<num_rows; r++) {
                                double x = getNode(r);
                                pntr[r] = static_cast<int>(indx.size());
                                ancestors.clear();
                                ancestors_vals.clear();
                                int kid = r;
                                int dad = (kid + 1) / 2;
                                if (kid < 4) dad--;
                                while(dad != 0) {
                                    ancestors.push_back(dad);
                                    ancestors_vals.push_back( evalRaw(dad, x) );
                                    kid = dad;
                                    dad = (kid + 1) / 2;
                                    if (kid < 4) dad--;
                                }
                                indx.push_back(0);
                                indx.insert(indx.end(), ancestors.rbegin(), ancestors.rend());
                                indx.push_back(r);
                                vals.push_back(1.0);
                                vals.insert(vals.end(), ancestors_vals.rbegin(), ancestors_vals.rend());
                                vals.push_back(1.0);
                            }
                            pntr.back() = static_cast<int>(indx.size());
                        }
                        break;
                };
            } else if (rule == rule_semilocalp) {
                switch(num_rows) {
                    case 0:
                        indx = {0,};
                        vals = {1.0,};
                        pntr = {0, 1};
                        break;
                    case 1:
                        indx = {0, 0, 1};
                        vals = {1.0, 1.0, 1.0};
                        pntr = {0, 1, 3};
                        break;
                    case 2:
                        indx = {0, 0, 1, 0, 2};
                        vals = {1.0, 1.0, 1.0, 1.0, 1.0};
                        pntr = {0, 1, 3, 5};
                        break;
                    default:
                        indx.clear();
                        indx.reserve(num_rows * max_level);
                        vals.clear();
                        vals.reserve(num_rows * max_level);
                        pntr = std::vector<int>(num_rows + 1, 0);
                        for(auto i : std::array<int, 5>{0, 0, 1, 0, 2})
                            indx.push_back(i);
                        for(int i=0; i<5; i++)
                            vals.push_back(1.0);
                        pntr[1] = 1;
                        pntr[2] = 3;
                        {
                            std::vector<int> ancestors;
                            std::vector<double> ancestors_vals;
                            ancestors.reserve(max_level);
                            ancestors_vals.reserve(max_level);
                            for(int r=3; r<num_rows; r++) {
                                pntr[r] = static_cast<int>(indx.size());
                                double x = getNode(r);
                                ancestors.clear();
                                ancestors_vals.clear();
                                int kid = r;
                                int dad = (kid + 1) / 2;
                                while(dad > 2) {
                                    ancestors.push_back(dad);
                                    ancestors_vals.push_back( evalRaw(dad, x) );
                                    kid = dad;
                                    dad = (kid + 1) / 2;
                                }
                                indx.push_back(0);
                                indx.push_back(1);
                                indx.push_back(2);
                                indx.insert(indx.end(), ancestors.rbegin(), ancestors.rend());
                                indx.push_back(r);
                                vals.push_back(1.0);
                                vals.push_back( evalRaw(1, x) );
                                vals.push_back( evalRaw(2, x) );
                                vals.insert(vals.end(), ancestors_vals.rbegin(), ancestors_vals.rend());
                                vals.push_back(1.0);
                            }
                            pntr.back() = static_cast<int>(indx.size());
                        }
                        break;
                };
            } else if (rule == rule_localp0) {
                switch(num_rows) {
                    case 0:
                        indx = {0,};
                        vals = {1.0,};
                        pntr = {0, 1};
                        break;
                    default:
                        indx.clear();
                        indx.reserve(num_rows * max_level);
                        vals.clear();
                        vals.reserve(num_rows * max_level);
                        pntr = std::vector<int>(num_rows + 1, 0);
                        indx.push_back(0);
                        vals.push_back(1.0);
                        {
                            std::vector<int> ancestors;
                            std::vector<double> ancestors_vals;
                            ancestors.reserve(max_level);
                            ancestors_vals.reserve(max_level);
                            for(int r=1; r<num_rows; r++) {
                                double x = getNode(r);
                                pntr[r] = static_cast<int>(indx.size());
                                ancestors.clear();
                                ancestors_vals.clear();
                                int kid = r;
                                int dad = (kid - 1) / 2;
                                while(dad != 0) {
                                    ancestors.push_back(dad);
                                    ancestors_vals.push_back( evalRaw(dad, x) );
                                    kid = dad;
                                    dad = (kid - 1) / 2;
                                }
                                indx.push_back(0);
                                indx.insert(indx.end(), ancestors.rbegin(), ancestors.rend());
                                indx.push_back(r);
                                vals.push_back( evalRaw(0, x) );
                                vals.insert(vals.end(), ancestors_vals.rbegin(), ancestors_vals.rend());
                                vals.push_back(1.0);
                            }
                            pntr.back() = static_cast<int>(indx.size());
                        }
                        break;
                };
            } else { // if (rule == rule_localpb) {
                switch(num_rows) {
                    case 0:
                        indx = {0,};
                        vals = {1.0,};
                        pntr = {0, 1};
                        break;
                    case 1:
                        indx = {0, 1};
                        vals = {1.0, 1.0};
                        pntr = {0, 1, 2};
                        break;
                    default:
                        indx.clear();
                        indx.reserve(num_rows * max_level);
                        vals.clear();
                        vals.reserve(num_rows * max_level);
                        pntr = std::vector<int>(num_rows + 1, 0);
                        indx.push_back(0);
                        indx.push_back(1);
                        vals.push_back(1.0);
                        vals.push_back(1.0);
                        pntr[1] = 1;
                        {
                            std::vector<int> ancestors;
                            std::vector<double> ancestors_vals;
                            ancestors.reserve(max_level);
                            ancestors_vals.reserve(max_level);
                            for(int r=2; r<num_rows; r++) {
                                pntr[r] = static_cast<int>(indx.size());
                                double x = getNode(r);
                                ancestors.clear();
                                ancestors_vals.clear();
                                int kid = r;
                                int dad = (kid + 1) / 2;
                                while(dad > 1) {
                                    ancestors.push_back(dad);
                                    ancestors_vals.push_back( evalRaw(dad, x) );
                                    kid = dad;
                                    dad = (kid + 1) / 2;
                                }
                                indx.push_back(0);
                                indx.push_back(1);
                                indx.insert(indx.end(), ancestors.rbegin(), ancestors.rend());
                                indx.push_back(r);
                                vals.push_back( evalRaw(0, x) );
                                vals.push_back( evalRaw(1, x) );
                                vals.insert(vals.end(), ancestors_vals.rbegin(), ancestors_vals.rend());
                                vals.push_back(1.0);
                            }
                            pntr.back() = static_cast<int>(indx.size());
                        }
                        break;
                };
            }
        }
    }

protected:
    double scaleX(int point, double x) const{
        if (rule == rule_localp0){
            if (point == 0) return x;
            return ((double) Maths::int2log2(point + 1) * (x + 3.0) - 3.0 - (double) (2*point));
        }else{
            if (rule == rule_localp){
                if (point == 0) return x;
                if (point == 1) return (x + 1.0);
                if (point == 2) return (x - 1.0);
            }else if (rule == rule_localpb){
                if (point == 0) return (x + 1.0) / 2.0;
                if (point == 1) return (x - 1.0) / 2.0;
                if (point == 2) return x;
            }
            return ((double) Maths::int2log2(point - 1) * (x + 3.0) + 1.0 - (double) (2*point));
        }
    }

    double scaleDiffX(int point) const{
        if (rule == rule_localp0) {
            if (point == 0) return 1.0;
            return (double) Maths::int2log2(point + 1);
        } else {
            if (rule == rule_localp) {
                if (point == 0 or point == 1 or point == 2) return 1.0;
            } else if (rule == rule_localpb) {
                if (point == 0 or point == 1) return 0.5;
                if (point == 2) return 1.0;
            }
            return (double) Maths::int2log2(point - 1);
        }
    }

    double evalPWQuadratic(int point, double x) const{
        if (rule == rule_localp){
            if (point == 1) return 1.0 - x;
            if (point == 2) return 1.0 + x;
        }else if (rule == rule_localpb){
            if (point == 0) return 1.0 - x;
            if (point == 1) return 1.0 + x;
        }
        return (1.0 - x) * (1.0 + x);
    }
    double evalPWCubic(int point, double x) const{
        if (rule == rule_localp){
            if (point == 0) return 1.0;
            if (point == 1) return 1.0 - x;
            if (point == 2) return 1.0 + x;
            if (point <= 4) return (1.0 - x) * (1.0 + x);
        }else if (rule == rule_localpb){
            if (point == 0) return 1.0 - x;
            if (point == 1) return 1.0 + x;
            if (point == 2) return (1.0 - x) * (1.0 + x);
        }else if (rule == rule_localp0){
            if (point == 0) return (1.0 - x) * (1.0 + x);
        }
        return (point % 2 == 0) ? (1.0 - x) * (1.0 + x) * (3.0 + x) / 3.0 : (1.0 - x) * (1.0 + x) * (3.0 - x) / 3.0;
    }
    double evalPWPower(int point, double x) const{
        if (rule == rule_localp)     if (point <= 8) return evalPWCubic(point, x); // if order is cubic or less, use the hard-coded functions
        if (rule == rule_semilocalp) if (point <= 4) return evalPWCubic(point, x);
        if (rule == rule_localpb)    if (point <= 4) return evalPWCubic(point, x);
        if (rule == rule_localp0)    if (point <= 2) return evalPWCubic(point, x);
        int level = getLevel(point);
        int most_turns = 1;
        double value = (1.0 - x)*(1.0 + x), phantom_distance = 1.0;
        int max_ancestors; // maximum number of ancestors to consider, first set the possible max then constrain by max_order
        if (rule == rule_localp)     max_ancestors = level-2; // number of ancestors to consider (two counted in the initialization of value)
        if (rule == rule_semilocalp) max_ancestors = level-1; // semi-localp and localpb add one more ancestor due to the global support at levels 1 and 0 (respectively)
        if (rule == rule_localpb)    max_ancestors = level-1;
        if (rule == rule_localp0)    max_ancestors = level; // localp0 adds two ancestors, the assumed zero nodes at the boundary
        if (max_order > 0) max_ancestors = std::min(max_ancestors, max_order - 2); // use the minimum of the available ancestors or the order restriction

        for(int j=0; j < max_ancestors; j++){
            // Lagrange polynomial needs to be constructed using normalized x (basis support (-1, 1)).
            // The support of the basis is equal to 2, thus we use units of "half-support"
            // The first two nodes are used in the initialization of value, those are the nearest ancestors (in spacial distance, not hierarchy level).
            // The other nodes are "phantoms" that lay strictly outside of the support [-1, 1] (i.e., not on the edge)
            // The walking distance more than doubles and is an odd number (due to the half-support)
            // The walk through the ancestors can take a left or right turn at each step,
            //   most_turns is the total number of turns possible for the current ancestor, turns (or most_turns - 1 - turns) is the actual number of turns
            // Every time we turn, we backtrack and we lose 2 units, the phantom node is at maximum distance minus 2 time the number of turns (to the left or right)
            most_turns *= 2;
            phantom_distance = 2.0 * phantom_distance + 1.0;
            int turns = (rule == rule_localp0) ? ((point+1) % most_turns) : ((point-1) % most_turns);
            double node = (turns < most_turns / 2) ? (phantom_distance - 2.0 * ((double) turns)) : (-phantom_distance + 2.0 * ((double) (most_turns - 1 - turns)));
            value *= ( x - node ) / ( - node);
        }
        return value;
    }

    // Same logic as the corresponding eval functions.
    double diffPWQuadratic(int point, double x) const{
        if (rule == rule_localp) {
            if (point == 1) return -1.0;
            if (point == 2) return 1.0;
        } else if (rule == rule_localpb) {
            if (point == 0) return -1.0;
            if (point == 1) return 1.0;
        }
        return -2.0 * x;
    }
    double diffPWCubic(int point, double x) const{
        if (rule == rule_localp) {
            if (point == 0) return 0.0;
            if (point == 1) return -1.0;
            if (point == 2) return 1.0;
            if (point <= 4) return -2.0 * x;
        } else if (rule == rule_localpb) {
            if (point == 0) return -1.0;
            if (point == 1) return 1.0;
            if (point == 2) return -2.0 * x;
        } else if (rule == rule_localp0) {
            if (point == 0) return -2.0 * x;
        }
        return (point % 2 == 0) ? 1.0 / 3.0 - x * (x + 2.0) : -1.0 / 3.0 + x * (x - 2.0);
    }
    double diffPWPower(int point, double x) const{
        if (rule == rule_localp     and point <= 8) return diffPWCubic(point, x);
        if (rule == rule_semilocalp and point <= 4) return diffPWCubic(point, x);
        if (rule == rule_localpb    and point <= 4) return diffPWCubic(point, x);
        if (rule == rule_localp0    and point <= 2) return diffPWCubic(point, x);
        int level = getLevel(point);
        int max_ancestors;
        if (rule == rule_localp)     max_ancestors = level-2;
        if (rule == rule_semilocalp) max_ancestors = level-1;
        if (rule == rule_localpb)    max_ancestors = level-1;
        if (rule == rule_localp0)    max_ancestors = level;
        if (max_order > 0) max_ancestors = std::min(max_ancestors, max_order - 2);

        // This lambda captures most_turns and phantom_distance by reference, uses those as internal state variables, and on each
        // call it returns the next normalized ancestor node.
        int most_turns = 1;
        double phantom_distance = 1.0;
        auto update_and_get_next_node = [&]() {
            most_turns *= 2;
            phantom_distance = 2.0 * phantom_distance + 1.0;
            int turns = (rule == rule_localp0) ? ((point+1) % most_turns) : ((point-1) % most_turns);
            return (turns < most_turns / 2) ?
                    (phantom_distance - 2.0 * ((double) turns)) :
                    (-phantom_distance + 2.0 * ((double) (most_turns - 1 - turns)));
        };

        // This lambda is the inverse transform of update_and_get_next_node() above, and returns the previous descendant node.
        auto rollback_and_get_prev_node = [&]() {
            most_turns /= 2;
            phantom_distance = 0.5 * (phantom_distance - 1.0);
            int turns = (rule == rule_localp0) ? ((point+1) % most_turns) : ((point-1) % most_turns);
            return (turns < most_turns / 2) ?
                    (phantom_distance - 2.0 * ((double) turns)) :
                    (-phantom_distance + 2.0 * ((double) (most_turns - 1 - turns)));
        };

        // Does not include the additional factor (1-x) * (1+x) and the Lagrange coefficient.
        std::vector<double> left_prods(max_ancestors);
        left_prods[0] = 1.0;
        double node = update_and_get_next_node();
        double coeff = 1.0 / (-node);
        for(int j=1; j<max_ancestors; j++) {
            left_prods[j] = left_prods[j-1] * (x - node);
            node = update_and_get_next_node();
            coeff *= 1.0 / (-node);
        }
        double right_prod = 1.0;
        double derivative = left_prods[max_ancestors-1];
        for (int j=max_ancestors-2; j>=0; j--) {
            right_prod *= x - node;
            derivative += right_prod * left_prods[j];
            node = rollback_and_get_prev_node();
        }

        // Adjust for the additional factor (1-x) * (1+x) and the Lagrange coefficient.
        derivative = derivative * (1.0 - x) * (1.0 + x) + right_prod * (x - node) * (-2.0) * x;
        derivative *= coeff;

        return derivative;
    }
};

inline std::unique_ptr<BaseRuleLocalPolynomial> makeRuleLocalPolynomial(TypeOneDRule rule, int order){
    if (order == 0){
        return Utils::make_unique<templRuleLocalPolynomial<rule_localp, true>>();
    }else{
        std::unique_ptr<BaseRuleLocalPolynomial> rule_pntr = [=]()->std::unique_ptr<BaseRuleLocalPolynomial>{
            if (rule == rule_localp){
                return Utils::make_unique<templRuleLocalPolynomial<rule_localp, false>>();
            }else if (rule == rule_semilocalp){
                return Utils::make_unique<templRuleLocalPolynomial<rule_semilocalp, false>>();
            }else if (rule == rule_localp0){
                return Utils::make_unique<templRuleLocalPolynomial<rule_localp0, false>>();
            }else{ // must be (rule == rule_localpb)
                return Utils::make_unique<templRuleLocalPolynomial<rule_localpb, false>>();
            }
        }();
        rule_pntr->setMaxOrder(order);
        return rule_pntr;
    }
}

#endif // __TASMANIAN_DOXYGEN_SKIP


}

#endif
