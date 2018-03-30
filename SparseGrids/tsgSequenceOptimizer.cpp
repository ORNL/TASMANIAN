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

GreedySequences::GreedySequences(){}
GreedySequences::~GreedySequences(){}

double* GreedySequences::getLejaNodes(int n) const{
    double *nodes = new double[n];
    nodes[0] = 0.0;
    if (n > 1) nodes[1] =  1.0;
    if (n > 2) nodes[2] = -1.0;
    if (n > 3) nodes[3] = sqrt(1.0 / 3.0);
    for(int i=4; i<n; i++){
        Residual g(i, nodes);
        OptimizerResult R = Optimizer::argMaxGlobal(&g);
        nodes[i] = R.xmax;
    }
    return nodes;
}
double* GreedySequences::getMaxLebesgueNodes(int n) const{
    double *nodes = new double[n];
    nodes[0] = 0.0;
    if (n > 1) nodes[1] =  1.0;
    if (n > 2) nodes[2] = -1.0;
    if (n > 3) nodes[3] = 0.5;
    for(int i=4; i<n; i++){
        MaxLebesgue g(i, nodes);
        OptimizerResult R = Optimizer::argMaxGlobal(&g);
        nodes[i] = R.xmax;
        //nodes[i] = argMaxGlobal(functional_lebesgue_function, i, nodes);
    }
    return nodes;
}
double* GreedySequences::getMinLebesgueNodes(int n) const{
    double *nodes = new double[n];
    int stored = getNumMinLebesgueStored();
    int usefirst = (n > stored) ? stored : n;
    for(int i=0; i<usefirst; i++){
        nodes[i] = getMinLebesgueStored(i);
    }
    if (n > stored){
        for(int i=stored; i<n; i++){
            MinLebesgue g(i, nodes);
            OptimizerResult R = Optimizer::argMaxGlobal(&g);
            nodes[i] = R.xmax;
        }
    }
    return nodes;
}
double* GreedySequences::getMinDeltaNodes(int n) const{
    double *nodes = new double[n];
    int stored = getNumMinDeltaStored();
    int usefirst = (n > stored) ? stored : n;
    for(int i=0; i<usefirst; i++){
        nodes[i] = getMinDeltaStored(i);
    }
    if (n > stored){
        cout << "Computing extra Min-Delta nodes" << endl;
        for(int i=stored; i<n; i++){
            MinDelta g(i, nodes);
            OptimizerResult R = Optimizer::argMaxGlobal(&g);
            nodes[i] = R.xmax;
        }
    }
    return nodes;
}
double GreedySequences::findLebesgueConstant(int n, const double nodes[]) const{
    OptimizerResult R;
    MaxLebesgue g(n, nodes);
    R = Optimizer::argMaxGlobal(&g);
    return R.fmax;
}
int GreedySequences::getNumMinLebesgueStored() const{ return 50; }
double GreedySequences::getMinLebesgueStored(int i) const{
    switch(i){
        case  0: return  0.00000000000000000e+00;
        case  1: return  1.00000000000000000e+00;
        case  2: return -1.00000000000000000e+00;
        case  3: return  5.77350269189625731e-01;
        case  4: return -6.82983516382591915e-01;
        case  5: return  3.08973641709742008e-01;
        case  6: return -8.88438098101029028e-01;
        case  7: return  9.01204348257380272e-01;
        case  8: return -3.79117423735989223e-01;
        case  9: return  7.68752974570125480e-01;
        case 10: return -5.49318186832010835e-01;
        case 11: return -9.64959809152743819e-01;
        case 12: return  9.67682365886019302e-01;
        case 13: return -1.28290369854639041e-01;
        case 14: return  4.45625686400429988e-01;
        case 15: return -8.04184165478832425e-01;
        case 16: return  1.67074539986499709e-01;
        case 17: return  8.44556006683382599e-01;
        case 18: return -4.67815101214832718e-01;
        case 19: return -9.88223645152207397e-01;
        case 20: return  6.75111201618315282e-01;
        case 21: return -2.31781319689011223e-01;
        case 22: return  9.90353577025220644e-01;
        case 23: return -8.50249643305543756e-01;
        case 24: return  3.78518581322005110e-01;
        case 25: return -7.34179511891844605e-01;
        case 26: return  9.37348529794296281e-01;
        case 27: return  7.43693755708555865e-02;
        case 28: return -9.43140751144033618e-01;
        case 29: return  7.22894448368867515e-01;
        case 30: return -3.13966521558236233e-01;
        case 31: return  5.21332302672803283e-01;
        case 32: return -6.22914216661596742e-01;
        case 33: return  9.81207150951301732e-01;
        case 34: return -9.96112540596109541e-01;
        case 35: return  2.27115897487506796e-01;
        case 36: return  8.71343953281246475e-01;
        case 37: return -5.83881096023192936e-01;
        case 38: return -7.07087869617637754e-02;
        case 39: return -9.21222090369298474e-01;
        case 40: return  6.35275009690072778e-01;
        case 41: return  9.97425667442937813e-01;
        case 42: return -7.74028611711955805e-01;
        case 43: return  1.19609717340386237e-01;
        case 44: return  8.06079215959057738e-01;
        case 45: return -9.77828412798477653e-01;
        case 46: return -2.81002869364477770e-01;
        case 47: return  4.81831205470734381e-01;
        case 48: return -4.33236587065744694e-01;
        case 49: return  9.51233414543255273e-01;
        default :
            return 0.0;
    }
}
int GreedySequences::getNumMinDeltaStored() const{ return 50; }
double GreedySequences::getMinDeltaStored(int i) const{
    switch(i){
        case  0: return  0.00000000000000000e+00;
        case  1: return  1.00000000000000000e+00;
        case  2: return -1.00000000000000000e+00;
        case  3: return  5.85786437626666157e-01;
        case  4: return -7.11391303688287735e-01;
        case  5: return -3.65515299787630255e-01;
        case  6: return  8.73364116971659943e-01;
        case  7: return -9.09934766576546594e-01;
        case  8: return  3.48149066530255846e-01;
        case  9: return  7.61022503813320816e-01;
        case 10: return -5.63728220989981099e-01;
        case 11: return  1.61287729090120513e-01;
        case 12: return -9.70195195938610255e-01;
        case 13: return  9.71760208817352256e-01;
        case 14: return -2.12010912937840634e-01;
        case 15: return -8.19500788126849566e-01;
        case 16: return  9.28317867167342214e-01;
        case 17: return  4.63681477197408431e-01;
        case 18: return -4.75500758437281013e-01;
        case 19: return  6.83991900744136849e-01;
        case 20: return -9.89565716162813747e-01;
        case 21: return -9.47327399278441035e-02;
        case 22: return  9.91020313545615483e-01;
        case 23: return -7.70894179555761561e-01;
        case 24: return  2.57418053836716398e-01;
        case 25: return -9.38159793780884654e-01;
        case 26: return  8.23375672488699584e-01;
        case 27: return -2.90168112727829441e-01;
        case 28: return  5.28028188455664571e-01;
        case 29: return -6.43769278853295157e-01;
        case 30: return  9.51854409585821237e-01;
        case 31: return -8.71963189117939352e-01;
        case 32: return  8.85485426816430832e-02;
        case 33: return  7.23563624423755547e-01;
        case 34: return -9.96440136176172664e-01;
        case 35: return -1.50233510276863047e-01;
        case 36: return  9.96883045450636107e-01;
        case 37: return -5.22343767211049914e-01;
        case 38: return  4.02868136231700480e-01;
        case 39: return -8.46787710296258989e-01;
        case 40: return  8.96281851550158049e-01;
        case 41: return -4.11624139689389490e-01;
        case 42: return  6.30422054421239331e-01;
        case 43: return -9.57396180946457842e-01;
        case 44: return  2.09670340519686721e-01;
        case 45: return  8.47208170728777743e-01;
        case 46: return -6.80096440009749670e-01;
        case 47: return -4.07635282203260146e-02;
        case 48: return  9.81787150881257009e-01;
        case 49: return -9.81715548483999001e-01;
        default : //cerr << "ERROR: GreedySequences::getMinDeltaStored called with wrong i, how did this happen?" << endl;
            return 0.0;
    }
}

Functional::Functional(){}
Functional::~Functional(){}
double Functional::getDiff(double) const{ return 0.0; }
double* Functional::makeCoeff(int num_nodes, const double nodes[]){
    double* coeff = new double[num_nodes];
    for(int i=0; i<num_nodes; i++){
        coeff[i] = 1.0;
        for(int j=0; j<i; j++){
            coeff[i] *= (nodes[i] - nodes[j]);
        }
        for(int j=i+1; j<num_nodes; j++){
            coeff[i] *= (nodes[i] - nodes[j]);
        }
    }
    return coeff;
}
double* Functional::evalLag(int num_nodes, const double nodes[], const double coeff[], double x){
    double *lag = new double[num_nodes];
    lag[0] = 1.0;
    for(int i=0; i<num_nodes-1; i++){
        lag[i+1] = (x - nodes[i]) * lag[i];
    }
    double w = 1.0;
    lag[num_nodes-1] /= coeff[num_nodes-1];
    for(int i=num_nodes-2; i>=0; i--){
        w *= (x - nodes[i+1]);
        lag[i] *= w / coeff[i];
    }
    return lag;
}
double Functional::basisDx(int num_nodes, const double nodes[], const double coeff[], int inode, double x){
    double s = 1.0;
    double p = 1.0;
    double n = (inode != 0) ? (x - nodes[0]) : (x - nodes[1]);

    for(int j=1; j<inode; j++){
        p *= n;
        n = (x - nodes[j]);
        s *= n;
        s += p;
    }
    for(int j= ((inode == 0) ? 2 : inode+1); j<num_nodes; j++){
        p *= n;
        n = (x - nodes[j]);
        s *= n;
        s += p;
    }
    return s / coeff[inode];
}
double* Functional::sortIntervals(int num_nodes, const double nodes[]){
    double *sorted = new double[num_nodes];
    std::copy(nodes, nodes + num_nodes, sorted);
    qsort(sorted, num_nodes, sizeof(double), &nodeCompar);
    return sorted;
}
int Functional::nodeCompar(const void * a, const void * b){ return (*(double*) a < *(double*) b) ? -1 : 1; }

Residual::Residual(int cnum_nodes, const double cnodes[]) : num_nodes(cnum_nodes){
    nodes = new double[num_nodes];
    std::copy(cnodes, cnodes + num_nodes, nodes);
}
Residual::~Residual(){  delete[] nodes;  }
double Residual::getValue(double x) const{
    double prod = 1.0;
    for(int j=0; j<num_nodes; j++){ prod *= (x - nodes[j]); }
    return fabs(prod);
}
bool Residual::hasDerivative() const{  return true;  }
double Residual::getDiff(double x) const{
    double s = 1.0, p = 1.0, n = (x - nodes[0]);
    for(int j=1; j<num_nodes; j++){
        p *= n;
        n = (x - nodes[j]);
        s *= n;
        s += p;
    }
    return s;
}
int Residual::getNumIntervals() const{  return num_nodes-1;  }
double* Residual::getIntervals() const{
    return sortIntervals(num_nodes, nodes);
}

MaxLebesgue::MaxLebesgue(int cnum_nodes, const double cnodes[], double new_node) : num_nodes(cnum_nodes+1){
    nodes = new double[num_nodes];
    std::copy(cnodes, cnodes + num_nodes-1, nodes);
    nodes[num_nodes-1] = new_node;

    coeff = makeCoeff(num_nodes, nodes);
}
MaxLebesgue::MaxLebesgue(int cnum_nodes, const double cnodes[]) : num_nodes(cnum_nodes){
    nodes = new double[num_nodes];
    std::copy(cnodes, cnodes + num_nodes, nodes);

    coeff = makeCoeff(num_nodes, nodes);
}
MaxLebesgue::~MaxLebesgue(){  delete[] nodes;  delete[] coeff;  }
double MaxLebesgue::getValue(double x) const{
    double *lag = evalLag(num_nodes, nodes, coeff, x);
    double sum = fabs(lag[0]); for(int i=1; i<num_nodes; i++) sum += fabs(lag[i]);
    delete[] lag;
    return sum;
}
bool MaxLebesgue::hasDerivative() const{  return true;  }
double MaxLebesgue::getDiff(double x) const{
    double *lag = evalLag(num_nodes, nodes, coeff, x);
    double sum = 0.0;
    for(int i=0; i<num_nodes; i++) sum += (lag[i] > 0.0) ? basisDx(num_nodes, nodes, coeff, i, x) : -basisDx(num_nodes, nodes, coeff, i, x);
    delete[] lag;
    return sum;
}
int MaxLebesgue::getNumIntervals() const{  return num_nodes-1;  }
double* MaxLebesgue::getIntervals() const{
    return sortIntervals(num_nodes, nodes);
}


MinLebesgue::MinLebesgue(int cnum_nodes, const double cnodes[]) : num_nodes(cnum_nodes) {
    nodes = new double[num_nodes];
    std::copy(cnodes, cnodes + num_nodes, nodes);
}
MinLebesgue::~MinLebesgue(){  delete[] nodes;  }

double MinLebesgue::getValue(double x) const{
    for(int i=0; i<num_nodes; i++){
        if (fabs(x - nodes[i]) < 10*TSG_NUM_TOL){ // points almost overlap, will cause an explosion in the Lebesgue constant
            return -1.E+100;
        }
    }
    MaxLebesgue M(num_nodes, nodes, x);
    OptimizerResult R = Optimizer::argMaxGlobal(&M);
    return -R.fmax;
}
bool MinLebesgue::hasDerivative() const{  return false;  }
int MinLebesgue::getNumIntervals() const{ return num_nodes - 1; }
double* MinLebesgue::getIntervals() const{
    return sortIntervals(num_nodes, nodes);
}

MaxDelta::MaxDelta(int cnum_nodes, const double cnodes[], double new_node) : num_nodes(cnum_nodes + 1){
    nodes = new double[num_nodes];
    std::copy(cnodes, cnodes + num_nodes-1, nodes);
    nodes[num_nodes-1] = new_node;

    coeff_plus  = makeCoeff(num_nodes,   nodes);
    coeff_minus = makeCoeff(num_nodes-1, nodes);
}
MaxDelta::MaxDelta(int cnum_nodes, const double cnodes[]) : num_nodes(cnum_nodes){
    nodes = new double[num_nodes];
    std::copy(cnodes, cnodes + num_nodes, nodes);
    coeff_plus  = makeCoeff(num_nodes,   nodes);
    coeff_minus = makeCoeff(num_nodes-1, nodes);
}
MaxDelta::~MaxDelta(){  delete[] nodes;  delete[] coeff_plus;  delete[] coeff_minus;  }
double MaxDelta::getValue(double x) const{
    double* lag_plus  = evalLag(num_nodes, nodes, coeff_plus, x);
    double* lag_minus = evalLag(num_nodes-1, nodes, coeff_minus, x);
    double sum = 0.0;
    for(int i=0; i<num_nodes-1; i++) sum += fabs(lag_plus[i] - lag_minus[i]);
    sum += fabs(lag_plus[num_nodes-1]);
    delete[] lag_plus;  delete[] lag_minus;
    return sum;
}
bool MaxDelta::hasDerivative() const{  return true; }
double MaxDelta::getDiff(double x) const{
    double* lag_plus  = evalLag(num_nodes, nodes, coeff_plus, x);
    double* lag_minus = evalLag(num_nodes-1, nodes, coeff_minus, x);
    double sum = 0.0;
    for(int i=0; i<num_nodes-1; i++) sum += ((lag_plus[i] - lag_minus[i]) > 0.0) ?
            basisDx(num_nodes, nodes, coeff_plus, i, x) - basisDx(num_nodes-1, nodes, coeff_minus, i, x) :
            basisDx(num_nodes-1, nodes, coeff_minus, i, x) - basisDx(num_nodes, nodes, coeff_plus, i, x);
    sum += (lag_plus[num_nodes-1] > 0.0) ? basisDx(num_nodes, nodes, coeff_plus, num_nodes-1, x) :
                         -basisDx(num_nodes, nodes, coeff_plus, num_nodes-1, x);
    delete[] lag_plus;  delete[] lag_minus;
    return sum;
}
int MaxDelta::getNumIntervals() const{ return num_nodes - 1; }
double* MaxDelta::getIntervals() const{  return sortIntervals(num_nodes, nodes);  }

MinDelta::MinDelta(int cnum_nodes, const double cnodes[]) : num_nodes(cnum_nodes) {
    nodes = new double[num_nodes];
    std::copy(cnodes, cnodes + num_nodes, nodes);
}
MinDelta::~MinDelta(){  delete[] nodes;  }
double MinDelta::getValue(double x) const{
    for(int i=0; i<num_nodes; i++){
        if (fabs(x - nodes[i]) < 10*TSG_NUM_TOL){ // points almost overlap, will cause an explosion in the Lebesgue constant
            return -1.E+100;
        }
    }
    MaxDelta M(num_nodes, nodes, x);
    OptimizerResult R = Optimizer::argMaxGlobal(&M);
    return -R.fmax;
}
bool MinDelta::hasDerivative() const{  return false;  }
int MinDelta::getNumIntervals() const{  return num_nodes-1;  }
double* MinDelta::getIntervals() const{  return sortIntervals(num_nodes, nodes);  }


Optimizer::Optimizer(){}
Optimizer::~Optimizer(){}

OptimizerResult Optimizer::argMaxGlobal(const Functional *F){
    double *sorted = F->getIntervals();
    int num_intervals = F->getNumIntervals();

    OptimizerResult MaxResult;

    MaxResult.xmax = -1.0;
    MaxResult.fmax = F->getValue(MaxResult.xmax);

    double x = 1.0;
    double f = F->getValue(x);
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

    delete[] sorted;

    return MaxResult;
}
OptimizerResult Optimizer::argMaxLocalPattern(const Functional *F, double left, double right){
    double xl = left;
    double xr = right;
    double xm = 0.5 * (left + right);
    double fl = F->getValue(xl);
    double fm = F->getValue(xm);
    double fr = F->getValue(xr);
    double d = xr - xm;

    double tol = (F->hasDerivative()) ? TSG_NUM_TOL * 1.E-3 : TSG_NUM_TOL;

    while(d > tol){
        if ((fm >= fl) && (fm >= fr)){ // shrink the pattern
            d /= 2.0;
            xl = xm - d;
            xr = xm + d;
            fr = F->getValue(xr);
            fl = F->getValue(xl);
        }else if ((fl >= fm) && (fl >= fr)){ // if the left point is the maximum
            if (xl - d < left){ // if shifting left would get us outside the interval
                d /= 2.0;
                xm = xl + d;
                xr = xm + d;
                fm = F->getValue(xm);
                fr = F->getValue(xr);
            }else{ // shift left
                xr = xm;
                fr = fm;
                xm = xl;
                fm = fl;
                xl -= d;
                fl = F->getValue(xl);
            }
        }else{ // the right point is the largest
            if (xr + d > right){ // if shifting right would get us outside the interval
                d /= 2.0;
                xm = xr - d;
                xl = xm - d;
                fm = F->getValue(xm);
                fl = F->getValue(xl);
            }else{ // shift right
                xl = xm;
                fl = fm;
                xm = xr;
                fm = fr;
                xr += d;
                fr = F->getValue(xr);
            }
        }
    }
    // the post-processing here can give you a digit or two of extra accuracy
    // however, the secant method by itself is unstable when the number of poitns grows
    //if ((type == functional_residual) || (type == functional_lebesgue_function)){
    OptimizerResult R;
    if (F->hasDerivative()){
        R.xmax = argMaxLocalSecant(F, xl, xr);
        R.fmax = F->getValue(R.xmax);
    }else{
        R.fmax = fm;
        R.xmax = xm;
    }
    return R;
}
double Optimizer::argMaxLocalSecant(const Functional *F, double left, double right){
    double xm = left;
    double dm = F->getDiff(xm); // Functional1DDx(type, num_nodes, nodes, xm);
    double x = right;
    double d = F->getDiff(x); // Functional1DDx(type, num_nodes, nodes, x);

    if (fabs(d) > fabs(dm)){
        double t = x; x = xm; xm = t;
        t = d; d = dm; dm = t;
    }
    int itr = 0;
    while((fabs(d) > 3*TSG_NUM_TOL) && (itr < TSG_MAX_SECANT_ITERATIONS)){
        if (d == dm) std::cout << itr << "  " << fabs(d) << "  " << left << "  " << right << "  " << x << endl;
        double xp = x - d * (x - xm) / (d - dm);
        xm = x; dm = d; x = xp;

        d = F->getDiff(x); // Functional1DDx(type, num_nodes, nodes, x);

        itr++;
    }
    return (fabs(d) < fabs(dm)) ? x : xm;
}

}

#endif
