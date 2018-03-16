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
int BaseRuleLocalPolynomial::getMaxNumKids() const{ return 0; }
int BaseRuleLocalPolynomial::getMaxNumParents() const{ return 0; }

int BaseRuleLocalPolynomial::getParent(int) const{ return 0; }
int BaseRuleLocalPolynomial::getStepParent(int) const{ return 0; }
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

}

#endif
