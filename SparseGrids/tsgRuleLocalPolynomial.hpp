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

    //virtual bool isSemiLocal() const; // second level for quadratic functions, two ways
    virtual int getMaxNumKids() const;
    virtual int getMaxNumParents() const;

    virtual int getParent(int point) const;
    virtual int getStepParent(int point) const;
    //virtual int getKidLeft(int point) const;
    //virtual int getKidRight(int point) const;
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

class RuleLocalPolynomial : public BaseRuleLocalPolynomial{
public:
    RuleLocalPolynomial();
    ~RuleLocalPolynomial();

    int getMaxOrder() const;
    void setMaxOrder(int order);

    TypeOneDRule getType() const;

    int getNumPoints(int level) const;

    double getNode(int point) const;

    int getLevel(int point) const;
    double getSupport(int point) const;

    //bool isSemiLocal() const;
    int getMaxNumKids() const;
    int getMaxNumParents() const;

    int getParent(int point) const;
    int getStepParent(int point) const;
    int getKidLeft(int point) const;
    int getKidRight(int point) const;
    int getKid(int point, int kid_number) const;

    double evalRaw(int point, double x) const;
    double evalSupport(int point, double x, bool &isSupported) const;
    double getArea(int point, int n, const double w[], const double x[]) const;

protected:
    double scaleX(int point, double x) const;
    // those take a normalized x
    double evalPWQuadratic(int point, double x) const;
    double evalPWCubic(int point, double x) const;
    double evalPWPower(int point, double x) const;
};

class RuleSemiLocalPolynomial : public BaseRuleLocalPolynomial{
public:
    RuleSemiLocalPolynomial();
    ~RuleSemiLocalPolynomial();

    int getMaxOrder() const;
    void setMaxOrder(int order);

    TypeOneDRule getType() const;

    int getNumPoints(int level) const;

    double getNode(int point) const;

    int getLevel(int point) const;
    double getSupport(int point) const;

    //bool isSemiLocal() const;
    int getMaxNumKids() const;
    int getMaxNumParents() const;

    int getParent(int point) const;
    int getStepParent(int point) const;
    //int getKidLeft(int point) const;
    //int getKidRight(int point) const;
    int getKid(int point, int kid_number) const;

    double evalRaw(int point, double x) const;
    double evalSupport(int point, double x, bool &isSupported) const;
    double getArea(int point, int n, const double w[], const double x[]) const;

protected:
    double scaleX(int point, double x) const;
    // those take a normalized x
    double evalPWQuadratic(int point, double x) const;
    double evalPWCubic(int point, double x) const;
    double evalPWPower(int point, double x) const;
};

class RuleLocalPolynomialZero : public BaseRuleLocalPolynomial{
public:
    RuleLocalPolynomialZero();
    ~RuleLocalPolynomialZero();

    int getMaxOrder() const;
    void setMaxOrder(int order);

    TypeOneDRule getType() const;

    int getNumPoints(int level) const;

    double getNode(int point) const;

    int getLevel(int point) const;
    double getSupport(int point) const;

    //bool isSemiLocal() const;
    int getMaxNumKids() const;
    int getMaxNumParents() const;

    int getParent(int point) const;
    int getStepParent(int point) const;
    //int getKidLeft(int point) const;
    //int getKidRight(int point) const;
    int getKid(int point, int kid_number) const;

    double evalRaw(int point, double x) const;
    double evalSupport(int point, double x, bool &isSupported) const;
    double getArea(int point, int n, const double w[], const double x[]) const;

protected:
    double scaleX(int point, double x) const;
    // those take a normalized x
    double evalPWCubic(int point, double x) const;
    double evalPWPower(int point, double x) const;
};

class RuleLocalPolynomialConstant : public BaseRuleLocalPolynomial{
public:
    RuleLocalPolynomialConstant();
    ~RuleLocalPolynomialConstant();

    int getMaxOrder() const;
    void setMaxOrder(int order);

    TypeOneDRule getType() const;

    int getNumPoints(int level) const;
//
    double getNode(int point) const;
//
    int getLevel(int point) const;
    double getSupport(int point) const;
//
//    //bool isSemiLocal() const;
    int getMaxNumKids() const;
    int getMaxNumParents() const;
//
    int getParent(int point) const;
    int getStepParent(int point) const;
//    int getKidLeft(int point) const;
//    int getKidRight(int point) const;
    int getKid(int point, int kid_number) const;
    double evalRaw(int point, double x) const; // must check support here too
    double evalSupport(int point, double x, bool &isSupported) const; // bool should hold extended support, i.e., 4/3 * supp
    double getArea(int point, int n, const double w[], const double x[]) const;

protected:
    //double scaleX(int point, double x) const;
    // those take a normalized x
    //double evalPWCubic(int point, double x) const;
    //double evalPWPower(int point, double x) const;
};

}

#endif
