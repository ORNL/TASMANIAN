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

#ifndef __TSG_RULE_LOCAL_POLYNOMIAL_CPP
#define __TSG_RULE_LOCAL_POLYNOMIAL_CPP

#include "tsgRuleLocalPolynomial.hpp"

namespace TasGrid{

BaseRuleLocalPolynomial::BaseRuleLocalPolynomial() : max_order(0){}
BaseRuleLocalPolynomial::~BaseRuleLocalPolynomial(){}

int BaseRuleLocalPolynomial::getMaxOrder() const{ return max_order; }
void BaseRuleLocalPolynomial::setMaxOrder(int order){ max_order = order; }
TypeOneDRule BaseRuleLocalPolynomial::getType() const{ return rule_none; }

double BaseRuleLocalPolynomial::getNode(int) const{ return 0.0; }
int BaseRuleLocalPolynomial::getLevel(int) const{ return 0; }
double BaseRuleLocalPolynomial::getSupport(int) const{ return 0; }
//bool BaseRuleLocalPolynomial::isSemiLocal() const{  return false;  }
int BaseRuleLocalPolynomial::getMaxNumKids() const{ return 0; }
int BaseRuleLocalPolynomial::getMaxNumParents() const{ return 0; }

int BaseRuleLocalPolynomial::getParent(int) const{ return 0; }
int BaseRuleLocalPolynomial::getStepParent(int) const{ return 0; }
//int BaseRuleLocalPolynomial::getKidLeft(int point) const{  return 0;  }
//int BaseRuleLocalPolynomial::getKidRight(int point) const{  return 0;  }
int BaseRuleLocalPolynomial::getKid(int, int) const{ return 0; }

double BaseRuleLocalPolynomial::evalRaw(int, double) const{ return 0.0; }
double BaseRuleLocalPolynomial::evalSupport(int, double, bool&) const{ return 0.0; }
double BaseRuleLocalPolynomial::getArea(int, int, const double*, const double*) const{ return 0.0; }

int BaseRuleLocalPolynomial::intlog2(int i){ // this is effectively: floor(log_2(i))
    int result = 0;
    while (i >>= 1){ result++; }
    return result;
}
int BaseRuleLocalPolynomial::int2log2(int i){ // this is effectively: 2^(floor(log_2(i))), when i == 0, this returns 1
    int result = 1;
    while (i >>= 1){ result <<= 1; }
    return result;
}
int BaseRuleLocalPolynomial::int3log3(int i){ // this is effectively: 3^(ceil(log_3(i+1)))
    int result = 1;
    while(i >= 1){ i /= 3; result *= 3; }
    return result;
}

RuleLocalPolynomial::RuleLocalPolynomial(){}
RuleLocalPolynomial::~RuleLocalPolynomial(){}

int RuleLocalPolynomial::getMaxOrder() const{ return max_order; }
void RuleLocalPolynomial::setMaxOrder(int order){ max_order = order; }
TypeOneDRule RuleLocalPolynomial::getType() const{  return rule_localp; }
int RuleLocalPolynomial::getNumPoints(int level) const{
    return (level == 0) ? 1 : ((1 << level) + 1);
}
double RuleLocalPolynomial::getNode(int point) const{
    if (point == 0) return  0.0;
    if (point == 1) return -1.0;
    if (point == 2) return  1.0;
    return ((double)(2*point - 1) ) / ((double) int2log2(point - 1) ) - 3.0;
}
int RuleLocalPolynomial::getLevel(int point) const{
    return (point == 0) ? 0 : (point == 1) ? 1 : intlog2(point - 1) + 1;
}

int RuleLocalPolynomial::getMaxNumKids() const{ return 2; }
int RuleLocalPolynomial::getMaxNumParents() const{ return 1; }

int RuleLocalPolynomial::getParent(int point) const{
    int dad = (point + 1) / 2;
    if (point < 4) dad--;
    return dad;
}
int RuleLocalPolynomial::getStepParent(int) const{ return -1; }

int RuleLocalPolynomial::getKid(int point, int kid_number) const{
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
}
double RuleLocalPolynomial::evalRaw(int point, double x) const{
    if ((max_order == 0) || (point == 0)) return 1.0;
    double xn = scaleX(point, x);
    if (max_order == 1) return 1.0 - fabs(xn);
    if (max_order == 2) return evalPWQuadratic(point, xn);
    if (max_order == 3) return evalPWCubic(point, xn);
    return evalPWPower(point, xn);
}

double RuleLocalPolynomial::evalSupport(int point, double x, bool &isSupported) const{
    isSupported = true;
    if (point == 0) return 1.0;
    double xn = scaleX(point, x);
    if (fabs(xn) <= 1.0){
        if (max_order == 0) return 1.0;
        if (max_order == 1) return 1.0 - fabs(xn);
        if (max_order == 2) return evalPWQuadratic(point, xn);
        if (max_order == 3) return evalPWCubic(point, xn);
        return evalPWPower(point, xn);
    }else{
        isSupported = false;
        return 0.0;
    }
}
double RuleLocalPolynomial::getArea(int point, int n, const double w[], const double x[]) const{
    if (point == 0) return 2.0;
    if (max_order == 0){
        if ((point == 1) || (point == 2)) return 1.0;
        return 2.0 * getSupport(point);
    }
    if ((point == 1) || (point == 2)) return 0.5;
    if (max_order == 1) return getSupport(point);
    if ((max_order == 2) || (max_order == 3) || (point <= 8))  return (4.0/3.0) * getSupport(point);
    double sum = 0.0;
    for(int i=0; i<n; i++) sum += w[i] * evalPWPower(point, x[i]);
    return sum * getSupport(point);
}
double RuleLocalPolynomial::getSupport(int point) const{
    return (point == 0) ? 1.0 : 1.0 / ((double) int2log2(point - 1));
}
double RuleLocalPolynomial::scaleX(int point, double x) const{
    if (point == 0) return x;
    if (point == 1) return (x + 1.0);
    if (point == 2) return (x - 1.0);
    return ((double) int2log2(point - 1) * (x + 3.0) + 1.0 - (double) (2*point));
}
double RuleLocalPolynomial::evalPWQuadratic(int point, double x) const{
    if (point == 1) return 1.0 - x;
    if (point == 2) return 1.0 + x;
    return (1.0 - x) * (1.0 + x);
}
double RuleLocalPolynomial::evalPWCubic(int point, double x) const{
    if (point == 0) return 1.0;
    if (point == 1) return 1.0 - x;
    if (point == 2) return 1.0 + x;
    if (point <= 4) return (1.0 - x) * (1.0 + x);
    return (point % 2 == 0) ? (1.0 - x) * (1.0 + x) * (3.0 + x) / 3.0 : (1.0 - x) * (1.0 + x) * (3.0 - x) / 3.0;
}
double RuleLocalPolynomial::evalPWPower(int point, double x) const{
    if (point <= 8) return evalPWCubic(point, x); //cout << " x = " << x << endl;
    int level = getLevel(point);
    int mod = 1;
    double value = (1.0 - x)*(1.0 + x), offset = 1.0;
    int imax = (max_order < 0) ? level-2 : ((max_order < level) ? max_order -2 : level-2);
    for(int j=0; j < imax; j++){
        mod *= 2;
        offset = 2.0 * offset + 1.0;
        int z = (point-1) % mod;
        if (z < mod / 2){
            value *= (x - offset + 2.0 * ((double) z)) / (- offset + 2.0 * ((double) z));
        }else{
            value *= (x + offset - 2.0 * ((double) (mod - 1 - z))) / (offset - 2.0 * ((double) (mod - 1 - z)));
        }
    }
    return value;
}

RuleSemiLocalPolynomial::RuleSemiLocalPolynomial(){}
RuleSemiLocalPolynomial::~RuleSemiLocalPolynomial(){}

int RuleSemiLocalPolynomial::getMaxOrder() const{ return max_order; }
void RuleSemiLocalPolynomial::setMaxOrder(int order){ max_order = order; }
TypeOneDRule RuleSemiLocalPolynomial::getType() const{  return rule_semilocalp; }
int RuleSemiLocalPolynomial::getNumPoints(int level) const{
    return (level == 0) ? 1 : ((1 << level) + 1);
}
double RuleSemiLocalPolynomial::getNode(int point) const{
    if (point == 0){ return  0.0; }else
    if (point == 1){ return -1.0; }else
    if (point == 2){ return  1.0; }
    return ((double)(2*point - 1) ) / ((double) int2log2(point - 1) ) - 3.0;
}
int RuleSemiLocalPolynomial::getLevel(int point) const{
    return (point == 0) ? 0 : (point == 1) ? 1 : intlog2(point - 1) + 1;
}
//bool RuleSemiLocalPolynomial::isSemiLocal() const{  return true;  }
int RuleSemiLocalPolynomial::getMaxNumKids() const{  return 2;  }
int RuleSemiLocalPolynomial::getMaxNumParents() const{  return 2;  }

int RuleSemiLocalPolynomial::getParent(int point) const{
    int dad = (point + 1) / 2;
    if (point < 4) dad--;
    return dad;
}
int RuleSemiLocalPolynomial::getStepParent(int point) const{
    if (point == 3) return 2;
    if (point == 4) return 1;
    return -1;
}

int RuleSemiLocalPolynomial::getKid(int point, int kid_number) const{
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
}
double RuleSemiLocalPolynomial::evalRaw(int point, double x) const{
    if (point == 0) return 1.0;
    if (point == 1) return 0.5 * x * (x - 1.0);
    if (point == 2) return 0.5 * x * (x + 1.0);
    double xn = scaleX(point, x);
    if (max_order == 2) return evalPWQuadratic(point, xn);
    if (max_order == 3) return evalPWCubic(point, xn);
    return evalPWPower(point, xn);
}

double RuleSemiLocalPolynomial::evalSupport(int point, double x, bool &isSupported) const{
    isSupported = true;
    if (point == 0) return 1.0;
    if (point == 1) return 0.5 * x * (x - 1.0);
    if (point == 2) return 0.5 * x * (x + 1.0);
    double xn = scaleX(point, x);
    if (fabs(xn) <= 1.0){
        if (max_order == 2) return evalPWQuadratic(point, xn);
        if (max_order == 3) return evalPWCubic(point, xn);
        return evalPWPower(point, xn);
    }else{
        isSupported = false;
        return 0.0;
    }
}
double RuleSemiLocalPolynomial::getArea(int point, int n, const double w[], const double x[]) const{
    if (point == 0) return 2.0;
    if (max_order == 0){
        if ((point == 1) || (point == 2)) return 1.0;
        return 2.0 * getSupport(point);
    }
    if ((point == 1) || (point == 2)) return 1.0/3.0;
    if ((max_order == 2) || (max_order == 3) || (point <= 4))  return (4.0/3.0) * getSupport(point);
    double sum = 0.0;
    for(int i=0; i<n; i++) sum += w[i] * evalPWPower(point, x[i]);
    return sum * getSupport(point);
}
double RuleSemiLocalPolynomial::getSupport(int point) const{
    return (point == 0) ? 1.0 : 1.0 / ((double) int2log2(point - 1));
}
double RuleSemiLocalPolynomial::scaleX(int point, double x) const{
    return ((double) int2log2(point - 1) * (x + 3.0) + 1.0 - (double) (2*point));
}
double RuleSemiLocalPolynomial::evalPWQuadratic(int, double x) const{
    return (1.0 - x) * (1.0 + x);
}
double RuleSemiLocalPolynomial::evalPWCubic(int point, double x) const{
    return (point % 2 == 0) ? (1.0 - x) * (1.0 + x) * (3.0 + x) / 3.0 : (1.0 - x) * (1.0 + x) * (3.0 - x) / 3.0;
}
double RuleSemiLocalPolynomial::evalPWPower(int point, double x) const{
    if (point <= 4) return evalPWCubic(point, x); //cout << " x = " << x << endl;
    int level = getLevel(point);
    int mod = 1;
    double value = (1.0 - x)*(1.0 + x), offset = 1.0;
    int imax = (max_order < 0) ? level-1 : ((max_order-1 < level) ? max_order -2 : level-1);
    for(int j=0; j < imax; j++){
        mod *= 2;
        offset = 2.0 * offset + 1.0;
        int z = (point-1) % mod;
        if (z < mod / 2){
            value *= (x - offset + 2.0 * ((double) z)) / (- offset + 2.0 * ((double) z));
        }else{
            value *= (x + offset - 2.0 * ((double) (mod - 1 - z))) / (offset - 2.0 * ((double) (mod - 1 - z)));
        }
    }
    return value;
}

RuleLocalPolynomialZero::RuleLocalPolynomialZero(){}
RuleLocalPolynomialZero::~RuleLocalPolynomialZero(){}

int RuleLocalPolynomialZero::getMaxOrder() const{ return max_order; }
void RuleLocalPolynomialZero::setMaxOrder(int order){ max_order = order; }
TypeOneDRule RuleLocalPolynomialZero::getType() const{  return rule_localp0; }
int RuleLocalPolynomialZero::getNumPoints(int level) const{
    return (1 << (level+1)) -1;
}
double RuleLocalPolynomialZero::getNode(int point) const{
    if (point == 0) return  0.0;
    return ((double)(2*point +3) ) / ((double) int2log2(point + 1) ) - 3.0;
}
int RuleLocalPolynomialZero::getLevel(int point) const{
    return intlog2(point + 1);
}
//bool RuleLocalPolynomialZero::isSemiLocal() const{  return false;  }
int RuleLocalPolynomialZero::getMaxNumKids() const{  return 2;  }
int RuleLocalPolynomialZero::getMaxNumParents() const{  return 1;  }

int RuleLocalPolynomialZero::getParent(int point) const{
    if (point == 0) return -1;
    return (point - 1) / 2;
}
int RuleLocalPolynomialZero::getStepParent(int) const{  return -1;  }

int RuleLocalPolynomialZero::getKid(int point, int kid_number) const{
    return 2*point + ((kid_number == 0) ? 1 : 2);
}
double RuleLocalPolynomialZero::evalRaw(int point, double x) const{
    if (max_order == 0) return 1.0;
    double xn = scaleX(point, x);
    if (max_order == 1) return 1.0 - fabs(xn);
    if (max_order == 2) return (1.0 - xn) * (1.0 + xn);
    if (max_order == 3) return evalPWCubic(point, xn);
    return evalPWPower(point, xn);
}

double RuleLocalPolynomialZero::evalSupport(int point, double x, bool &isSupported) const{
    isSupported = true;
    double xn = scaleX(point, x);
    if (fabs(xn) <= 1.0){
        if (max_order == 0) return 1.0;
        if (max_order == 1) return 1.0 - fabs(xn);
        if (max_order == 2) return (1.0 - xn) * (1.0 + xn);
        if (max_order == 3) return evalPWCubic(point, xn);
        return evalPWPower(point, xn);
    }else{
        isSupported = false;
        return 0.0;
    }
}
double RuleLocalPolynomialZero::getArea(int point, int n, const double w[], const double x[]) const{
    if (max_order == 0){
        if (point == 0) return 2.0;
        return 2.0 * getSupport(point);
    }
    if (max_order == 1) return getSupport(point);
    if ((max_order == 2) || (max_order == 3) || (point <= 2))  return (4.0/3.0) * getSupport(point);
    double sum = 0.0;
    for(int i=0; i<n; i++) sum += w[i] * evalPWPower(point, x[i]);
    return sum * getSupport(point);
}
double RuleLocalPolynomialZero::getSupport(int point) const{
    return 1.0 / ((double) int2log2(point + 1));
}
double RuleLocalPolynomialZero::scaleX(int point, double x) const{
    if (point == 0) return x;
    return ((double) int2log2(point + 1) * (x + 3.0) - 3.0 - (double) (2*point));
}
double RuleLocalPolynomialZero::evalPWCubic(int point, double x) const{
    if (point == 0) return (1.0 - x) * (1.0 + x);
    return (point % 2 == 0) ? (1.0 - x) * (1.0 + x) * (3.0 + x) / 3.0 : (1.0 - x) * (1.0 + x) * (3.0 - x) / 3.0;
}
double RuleLocalPolynomialZero::evalPWPower(int point, double x) const{
    if (point <= 2) return evalPWCubic(point, x);
    int level = getLevel(point);
    int mod = 1;
    double value = (1.0 - x)*(1.0 + x), offset = 1.0;
    int imax = (max_order < 0) ? level : ((max_order-2 < level) ? max_order -2 : level);
    for(int j=0; j < imax; j++){
        mod *= 2;
        offset = 2.0 * offset + 1.0;
        int z = (point+1) % mod;
        if (z < mod / 2){
            value *= (x - offset + 2.0 * ((double) z)) / (- offset + 2.0 * ((double) z));
        }else{
            value *= (x + offset - 2.0 * ((double) (mod - 1 - z))) / (offset - 2.0 * ((double) (mod - 1 - z)));
        }
    }
    return value;
}

RuleLocalPolynomialConstant::RuleLocalPolynomialConstant(){}
RuleLocalPolynomialConstant::~RuleLocalPolynomialConstant(){}

int RuleLocalPolynomialConstant::getMaxOrder() const{ return max_order; }
void RuleLocalPolynomialConstant::setMaxOrder(int order){ max_order = order; }

TypeOneDRule RuleLocalPolynomialConstant::getType() const{ return rule_localp; }

int RuleLocalPolynomialConstant::getNumPoints(int level) const{
    int n = 1; while (level-- > 0) n *= 3; return n;
}
double RuleLocalPolynomialConstant::getNode(int point) const{
    /*if (point == 0) return 0.0;
    int i3l3 = int3log3(point);
    double supp = 1.0 / ((double) i3l3);
    if (point == i3l3/3) return -1.0 + supp;
    if (point == i3l3-1) return 1.0 - supp;
    i3l3 = (point - i3l3/3 + 1) / 2; // pair index
    double x = -1.0 + supp * ((double) (6*i3l3-1));
    if (point % 2 == 1) x += 2.0 * supp;
    return x;*/
    //double supp = 1.0 / ((double) int3log3(point));
    //return (point % 2 == 0) ? -2.0 + supp*(3*point+2) : -2.0 + supp*(3*point+1);
    return -2.0 + (1.0 / ((double) int3log3(point))) * (3*point + 2 - point % 2);
}
int RuleLocalPolynomialConstant::getLevel(int point) const{
    int level = 0;
    while(point >= 1){ point /= 3; level += 1; }
    return level;
}
double RuleLocalPolynomialConstant::getSupport(int point) const{
    return 1.0 / (double) int3log3(point);
}
int RuleLocalPolynomialConstant::getMaxNumKids() const{ return 4; }
int RuleLocalPolynomialConstant::getMaxNumParents() const{ return 2; }

int RuleLocalPolynomialConstant::getParent(int point) const{
    return (point == 0) ? -1 : point / 3;
}
int RuleLocalPolynomialConstant::getStepParent(int point) const{
    int i3l3 = int3log3(point);
    if (point == i3l3/3) return -1;
    if (point == i3l3-1) return -1;
    if ((point % 3 == 2) && (point % 2 == 0)) return point / 3 + 1;
    if ((point % 3 == 0) && (point % 2 == 1)) return point / 3 - 1;
    return -1;
}
int RuleLocalPolynomialConstant::getKid(int point, int kid_number) const{
    if (point == 0) return (kid_number == 0) ? 1 : (kid_number==1) ? 2 : -1;
    if (kid_number == 3){
        int i3l3 = int3log3(point);
        if (point == i3l3/3) return -1;
        if (point == i3l3-1) return -1;
        return (point % 2 == 0) ? 3*point + 3 : 3*point - 1;
    }
    return 3*point + kid_number;
}
double RuleLocalPolynomialConstant::evalRaw(int point, double x) const{
    if (fabs(x - getNode(point)) > getSupport(point)) return 0.0; // can speed this up
    return 1.0;
}
double RuleLocalPolynomialConstant::evalSupport(int point, double x, bool &isSupported) const{
    double distance = fabs(x - getNode(point)), support = getSupport(point); // can speed this up, see above
    isSupported = (distance <= (2.0) * support);
    return (distance > support) ? 0.0 : 1.0;
}
double RuleLocalPolynomialConstant::getArea(int point, int, const double*, const double*) const{
    return 2.0 * getSupport(point);
}

}

#endif
