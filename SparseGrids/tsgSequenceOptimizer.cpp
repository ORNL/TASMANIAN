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

VectorFunctional::VectorFunctional(){}
VectorFunctional::~VectorFunctional(){}

OptimizerResult argMaxGlobal(const VectorFunctional &F){
    std::vector<double> sorted;
    F.getIntervals(sorted);
    int num_intervals = (int) sorted.size() - 1;

    OptimizerResult MaxResult;

    MaxResult.xmax = -1.0;
    MaxResult.fmax = F.getValue(MaxResult.xmax);

    double x = 1.0;
    double f = F.getValue(x);
    if (f > MaxResult.fmax){
        MaxResult.xmax = x;
        MaxResult.fmax = f;
    }

    #pragma omp parallel
    {
        OptimizerResult LocalMaxResult = MaxResult, LocalResult;

        #pragma omp for schedule(dynamic)
        for(int i=0; i<num_intervals; i++){
            LocalResult = argMaxLocalPattern(F, sorted[i], sorted[i+1]);
            if (LocalResult.fmax > LocalMaxResult.fmax){
                LocalMaxResult = LocalResult;
            }

            #pragma omp critical
            {
                if (LocalMaxResult.fmax > MaxResult.fmax){
                    MaxResult = LocalMaxResult;
                }
            }
        }
    }

    return MaxResult;
}

OptimizerResult argMaxLocalPattern(const VectorFunctional &F, double left, double right){
    double xl = left;
    double xr = right;
    double xm = 0.5 * (left + right);
    double fl = F.getValue(xl);
    double fm = F.getValue(xm);
    double fr = F.getValue(xr);
    double d = xr - xm;

    double tol = (F.hasDerivative()) ? TSG_NUM_TOL * 1.E-3 : TSG_NUM_TOL;

    while(d > tol){
        if ((fm >= fl) && (fm >= fr)){ // shrink the pattern
            d /= 2.0;
            xl = xm - d;
            xr = xm + d;
            fr = F.getValue(xr);
            fl = F.getValue(xl);
        }else if ((fl >= fm) && (fl >= fr)){ // if the left point is the maximum
            if (xl - d < left){ // if shifting left would get us outside the interval
                d /= 2.0;
                xm = xl + d;
                xr = xm + d;
                fm = F.getValue(xm);
                fr = F.getValue(xr);
            }else{ // shift left
                xr = xm;
                fr = fm;
                xm = xl;
                fm = fl;
                xl -= d;
                fl = F.getValue(xl);
            }
        }else{ // the right point is the largest
            if (xr + d > right){ // if shifting right would get us outside the interval
                d /= 2.0;
                xm = xr - d;
                xl = xm - d;
                fm = F.getValue(xm);
                fl = F.getValue(xl);
            }else{ // shift right
                xl = xm;
                fl = fm;
                xm = xr;
                fm = fr;
                xr += d;
                fr = F.getValue(xr);
            }
        }
    }
    // the post-processing here can give you a digit or two of extra accuracy
    // however, the secant method by itself is unstable when the number of poitns grows
    OptimizerResult R;
    if (F.hasDerivative()){
        R.xmax = argMaxLocalSecant(F, xl, xr);
        R.fmax = F.getValue(R.xmax);
    }else{
        R.fmax = fm;
        R.xmax = xm;
    }
    return R;
}

double argMaxLocalSecant(const VectorFunctional &F, double left, double right){
    double xm = left;
    double dm = F.getDiff(xm);
    double x = right;
    double d = F.getDiff(x);

    if (fabs(d) > fabs(dm)){
        double t = x; x = xm; xm = t;
        t = d; d = dm; dm = t;
    }
    int itr = 0;
    while((fabs(d) > 3*TSG_NUM_TOL) && (itr < TSG_MAX_SECANT_ITERATIONS)){
        double xp = x - d * (x - xm) / (d - dm);
        xm = x; dm = d; x = xp;

        d = F.getDiff(x);

        itr++;
    }
    return (fabs(d) < fabs(dm)) ? x : xm;
}

void getPrecomputedMinLebesgueNodes(std::vector<double> &precomputed){
    precomputed = { 0.00000000000000000e+00,
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

void getPrecomputedMinDeltaNodes(std::vector<double> &precomputed){
    precomputed = { 0.00000000000000000e+00,
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

}

}

#endif
