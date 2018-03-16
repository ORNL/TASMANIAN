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

#include "tsgEnumerates.hpp"

namespace TasGrid{

class BaseRuleLocalPolynomial{
public:
    BaseRuleLocalPolynomial();
    ~BaseRuleLocalPolynomial();

    virtual int getMaxOrder() const;
    virtual void setMaxOrder(int order);

    virtual TypeOneDRule getType() const;

    virtual int getNumPoints(int level) const = 0; // for building the initial grid

    virtual double getNode(int point) const;

    virtual int getLevel(int point) const;
    virtual double getSupport(int point) const;

    virtual int getMaxNumKids() const;
    virtual int getMaxNumParents() const;

    virtual int getParent(int point) const;
    virtual int getStepParent(int point) const;
    virtual int getKid(int point, int kid_number) const;

    virtual double evalRaw(int point, double x) const; // normalizes x (i.e., (x-node) / support), but it does not check the support
    virtual double evalSupport(int point, double x, bool &isSupported) const; // // normalizes x (i.e., (x-node) / support) and checks if x is within the support

    virtual double getArea(int point, int n, const double w[], const double x[]) const;
    // integrate the function associated with the point, constant to cubic are known analytically, higher order need a 1-D quadrature rule

    static int intlog2(int i);
    static int int2log2(int i);
    static int int3log3(int i);

protected:
    int max_order;
};


template<TypeOneDRule rule, bool isZeroOrder>
class templRuleLocalPolynomial : public BaseRuleLocalPolynomial{
public:
    templRuleLocalPolynomial(){}
    ~templRuleLocalPolynomial(){}

    int getMaxOrder() const{ return max_order; }
    void setMaxOrder(int order){ max_order = order; }

    TypeOneDRule getType() const{ return rule; }

    int getNumPoints(int level) const{
        if (isZeroOrder){
            int n = 1;
            while (level-- > 0) n *= 3;
            return n;
        }else{
            if ((rule == rule_localp) || (rule == rule_semilocalp)){
                return (level == 0) ? 1 : ((1 << level) + 1);
            }else if (rule == rule_localp0){
                return (1 << (level+1)) -1;
            }
        }
    }

    double getNode(int point) const{
        if (isZeroOrder){
            return -2.0 + (1.0 / ((double) int3log3(point))) * (3*point + 2 - point % 2);
        }else{
            if (point == 0) return  0.0;
            if ((rule == rule_localp) || (rule == rule_semilocalp)){
                if (point == 1) return -1.0;
                if (point == 2) return  1.0;
                return ((double)(2*point - 1)) / ((double) int2log2(point - 1)) - 3.0;
            }else if (rule == rule_localp0){
                return ((double)(2*point +3) ) / ((double) int2log2(point + 1) ) - 3.0;
            }
        }
    }

    int getLevel(int point) const{
        if (isZeroOrder){
            int level = 0;
            while(point >= 1){ point /= 3; level += 1; }
            return level;
        }else{
            if ((rule == rule_localp) || (rule == rule_semilocalp)){
                return (point == 0) ? 0 : (point == 1) ? 1 : intlog2(point - 1) + 1;
            }else if (rule == rule_localp0){
                return intlog2(point + 1);
            }
        }
    }
    double getSupport(int point) const{
        if (isZeroOrder){
            return 1.0 / (double) int3log3(point);
        }else{
            if ((rule == rule_localp) || (rule == rule_semilocalp)){
                return (point == 0) ? 1.0 : 1.0 / ((double) int2log2(point - 1));
            }else if (rule == rule_localp0){
                return 1.0 / ((double) int2log2(point + 1));
            }
        }
    }

    int getMaxNumKids() const{ return (isZeroOrder) ? 4 : 2; }
    int getMaxNumParents() const{ return ((isZeroOrder || (rule == rule_semilocalp)) ? 2 : 1); }

    int getParent(int point) const{
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
            }
        }
    }

    int getStepParent(int point) const{
        if (isZeroOrder){
            int i3l3 = int3log3(point);
            if (point == i3l3/3) return -1;
            if (point == i3l3-1) return -1;
            if ((point % 3 == 2) && (point % 2 == 0)) return point / 3 + 1;
            if ((point % 3 == 0) && (point % 2 == 1)) return point / 3 - 1;
            return -1;
        }else{
            if (rule == rule_semilocalp){
                if (point == 3) return 2;
                if (point == 4) return 1;
            }
            return -1;
        }
    }

    int getKid(int point, int kid_number) const{
        if (isZeroOrder){
            if (point == 0) return (kid_number == 0) ? 1 : (kid_number==1) ? 2 : -1;
            if (kid_number == 3){
                int i3l3 = int3log3(point);
                if (point == i3l3/3) return -1;
                if (point == i3l3-1) return -1;
                return (point % 2 == 0) ? 3*point + 3 : 3*point - 1;
            }
            return 3*point + kid_number;
        }else{
            if ((rule == rule_localp) || (rule == rule_semilocalp)){
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
            }
        }
    }

    double evalRaw(int point, double x) const{
        if (isZeroOrder){
            if (fabs(x - getNode(point)) > getSupport(point)) return 0.0; // maybe can speed this up
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
            if (rule != rule_semilocalp) if (max_order == 1) return 1.0 - fabs(xn);
            if (max_order == 2) return evalPWQuadratic(point, xn);
            if (max_order == 3) return evalPWCubic(point, xn);
            return evalPWPower(point, xn);
        }
    }

    double evalSupport(int point, double x, bool &isSupported) const{
        if (isZeroOrder){
            double distance = fabs(x - getNode(point)), support = getSupport(point); // can speed this up, see above
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
            if (fabs(xn) <= 1.0){
                if ((rule == rule_localp) || (rule == rule_localp0)) if (max_order == 1) return 1.0 - fabs(xn);
                if (max_order == 2) return evalPWQuadratic(point, xn);
                if (max_order == 3) return evalPWCubic(point, xn);
                return evalPWPower(point, xn);
            }else{
                isSupported = false;
                return 0.0;
            }
        }
    }
    double getArea(int point, int n, const double w[], const double x[]) const{
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
            }else if (rule == rule_localp0){
                if (max_order == 1) return getSupport(point);
                if ((max_order == 2) || (max_order == 3) || (point <= 2))  return (4.0/3.0) * getSupport(point);
            }
            double sum = 0.0;
            for(int i=0; i<n; i++) sum += w[i] * evalPWPower(point, x[i]);
            return sum * getSupport(point);
        }
    }

protected:
    double scaleX(int point, double x) const{
        if (rule == rule_localp0){
            if (point == 0) return x;
            return ((double) int2log2(point + 1) * (x + 3.0) - 3.0 - (double) (2*point));
        }
        if ((rule == rule_localp) || (rule == rule_semilocalp)){
            if (rule == rule_localp){
                if (point == 0) return x;
                if (point == 1) return (x + 1.0);
                if (point == 2) return (x - 1.0);
            }
            return ((double) int2log2(point - 1) * (x + 3.0) + 1.0 - (double) (2*point));
        }
    }

    double evalPWQuadratic(int point, double x) const{
        if (rule == rule_localp){
            if (point == 1) return 1.0 - x;
            if (point == 2) return 1.0 + x;
        }
        return (1.0 - x) * (1.0 + x);
    }
    double evalPWCubic(int point, double x) const{
        if (rule == rule_localp){
            if (point == 0) return 1.0;
            if (point == 1) return 1.0 - x;
            if (point == 2) return 1.0 + x;
            if (point <= 4) return (1.0 - x) * (1.0 + x);
        }
        if (rule == rule_localp0) if (point == 0) return (1.0 - x) * (1.0 + x);
        return (point % 2 == 0) ? (1.0 - x) * (1.0 + x) * (3.0 + x) / 3.0 : (1.0 - x) * (1.0 + x) * (3.0 - x) / 3.0;
    }
    double evalPWPower(int point, double x) const{
        if (rule == rule_localp)     if (point <= 8) return evalPWCubic(point, x);
        if (rule == rule_semilocalp) if (point <= 4) return evalPWCubic(point, x);
        if (rule == rule_localp0)    if (point <= 2) return evalPWCubic(point, x);
        int level = getLevel(point);
        int mod = 1;
        double value = (1.0 - x)*(1.0 + x), offset = 1.0;
        int imax;
        if (rule == rule_localp)     imax = (max_order < 0) ? level-2 : ((max_order   < level) ? max_order -2 : level-2);
        if (rule == rule_semilocalp) imax = (max_order < 0) ? level-1 : ((max_order-1 < level) ? max_order -2 : level-1);
        if (rule == rule_localp0)    imax = (max_order < 0) ? level   : ((max_order-2 < level) ? max_order -2 : level);
        for(int j=0; j < imax; j++){
            mod *= 2;
            offset = 2.0 * offset + 1.0;
            int z;
            if ((rule == rule_localp) || (rule == rule_semilocalp)){
                z = (point-1) % mod;
            }else{
                z = (point+1) % mod;
            }
            if (z < mod / 2){
                value *= (x - offset + 2.0 * ((double) z)) / (- offset + 2.0 * ((double) z));
            }else{
                value *= (x + offset - 2.0 * ((double) (mod - 1 - z))) / (offset - 2.0 * ((double) (mod - 1 - z)));
            }
        }
        return value;
    }
};


}

#endif
