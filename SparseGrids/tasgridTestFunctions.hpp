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

#ifndef __TASMANIAN_TASGRID_FUNCTIONS_HPP
#define __TASMANIAN_TASGRID_FUNCTIONS_HPP

class BaseFunction{
public:
    BaseFunction();
    ~BaseFunction();

    virtual int getNumInputs() const = 0;
    virtual int getNumOutputs() const = 0;
    virtual const char* getDescription() const = 0;
    virtual void eval(const double x[], double y[]) const = 0;
    virtual void getIntegral(double y[]) const = 0;
};

class OneOneP0: public BaseFunction{
public:
    OneOneP0(); ~OneOneP0();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class OneOneP3: public BaseFunction{
public:
    OneOneP3(); ~OneOneP3();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class OneOneP4: public BaseFunction{
public:
    OneOneP4(); ~OneOneP4();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class OneOneExpMX: public BaseFunction{
public:
    OneOneExpMX(); ~OneOneExpMX();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class TwoOneP4: public BaseFunction{
public:
    TwoOneP4(); ~TwoOneP4();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class TwoOneP5: public BaseFunction{
public:
    TwoOneP5(); ~TwoOneP5();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class TwoOneExpNX2: public BaseFunction{
public:
    TwoOneExpNX2(); ~TwoOneExpNX2();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class ThreeOneExpNX2: public BaseFunction{
public:
    ThreeOneExpNX2(); ~ThreeOneExpNX2();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class TwoOneCos: public BaseFunction{
public:
    TwoOneCos(); ~TwoOneCos();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class TwoOneSinSin: public BaseFunction{
public:
    TwoOneSinSin(); ~TwoOneSinSin();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class TwoOneCosCos: public BaseFunction{
public:
    TwoOneCosCos(); ~TwoOneCosCos();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class TwoOneExpSinCos : public BaseFunction{
public:
    TwoOneExpSinCos(); ~TwoOneExpSinCos();
    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class TwoOneSinCosAxis : public BaseFunction{
public:
    TwoOneSinCosAxis(); ~TwoOneSinCosAxis();
    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class TwoTwoSinCos : public BaseFunction{
public:
    TwoTwoSinCos(); ~TwoTwoSinCos();
    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class TwoOneExpm40: public BaseFunction{
public:
    TwoOneExpm40(); ~TwoOneExpm40();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class FiveOneExpSum: public BaseFunction{
public:
    FiveOneExpSum(); ~FiveOneExpSum();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class SixOneExpSum: public BaseFunction{
public:
    SixOneExpSum(); ~SixOneExpSum();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class EightOneCosSum: public BaseFunction{
public:
    EightOneCosSum(); ~EightOneCosSum();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class ThreeOneUnitBall: public BaseFunction{
public:
    ThreeOneUnitBall(); ~ThreeOneUnitBall();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

// Function to test special integration rules, i.e. integrals against weight functions
class TwoOneConstGC1: public BaseFunction{
public:
    TwoOneConstGC1(); ~TwoOneConstGC1();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};
class TwoOneConstGC2: public BaseFunction{
public:
    TwoOneConstGC2(); ~TwoOneConstGC2();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};
class TwoOneConstGG: public BaseFunction{
public:
    TwoOneConstGG(); ~TwoOneConstGG();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};
class TwoOneConstGJ: public BaseFunction{
public:
    TwoOneConstGJ(); ~TwoOneConstGJ();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};
class TwoOneConstGGL: public BaseFunction{
public:
    TwoOneConstGGL(); ~TwoOneConstGGL();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};
class TwoOneConstGH: public BaseFunction{
public:
    TwoOneConstGH(); ~TwoOneConstGH();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class TwoOneENX2aniso: public BaseFunction{
public:
    TwoOneENX2aniso(); ~TwoOneENX2aniso();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};
class TwoTwoExpAsym: public BaseFunction{
public:
    TwoTwoExpAsym(); ~TwoTwoExpAsym();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};
class SixteenOneActive3: public BaseFunction{
public:
    SixteenOneActive3(); ~SixteenOneActive3();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class TwoOneDivisionAnisotropic: public BaseFunction{
public:
    TwoOneDivisionAnisotropic(); ~TwoOneDivisionAnisotropic();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};
class TwoOne1DCurved: public BaseFunction{
public:
    TwoOne1DCurved(); ~TwoOne1DCurved();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};
class TwoOneExpShiftedDomain: public BaseFunction{
public:
    TwoOneExpShiftedDomain(); ~TwoOneExpShiftedDomain();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};
class OneOneConformalOne: public BaseFunction{
public:
    OneOneConformalOne(); ~OneOneConformalOne();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};
class TwoOneConformalOne: public BaseFunction{
public:
    TwoOneConformalOne(); ~TwoOneConformalOne();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class Two3KExpSinCos: public BaseFunction{
public:
    Two3KExpSinCos(); ~Two3KExpSinCos();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

class TwoOneC1C2Periodic: public BaseFunction{
public:
    TwoOneC1C2Periodic(); ~TwoOneC1C2Periodic();

    int getNumInputs() const; int getNumOutputs() const; const char* getDescription() const; void eval(const double x[], double y[]) const; void getIntegral(double y[]) const;
};

#endif
