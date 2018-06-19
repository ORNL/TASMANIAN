
#ifndef __TASMANIAN_SPARSE_GRID_FOURIER_CPP
#define __TASMANIAN_SPARSE_GRID_FOURIER_CPP

#include <iostream>
#include <ctime>
#include "tsgGridFourier.hpp"
#include "tsgHiddenExternals.hpp"

namespace TasGrid{

GridFourier::GridFourier() : num_dimensions(0), num_outputs(0), wrapper(0), tensors(0), active_tensors(0), active_weights(0),
    points(0), needed(0), exponents(0), fourier_coefs(0), exponent_refs(0), tensor_refs(0), values(0), accel(0)
{}

GridFourier::GridFourier(const GridFourier &fourier) : num_dimensions(0), num_outputs(0), wrapper(0), tensors(0), active_tensors(0),
    active_weights(0), points(0), needed(0), exponents(0), fourier_coefs(0), exponent_refs(0), tensor_refs(0), values(0), accel(0) {
    copyGrid(&fourier);
}

GridFourier::~GridFourier() { reset(); }

void GridFourier::reset() {
    clearAccelerationData();
    if(exponent_refs != 0) { for(int i=0; i<active_tensors->getNumIndexes(); i++){ delete[] exponent_refs[i]; exponent_refs[i] = 0; } delete[] exponent_refs; exponent_refs = 0; }
    if(wrapper != 0) { delete wrapper; wrapper = 0; }
    if(tensors != 0) { delete tensors; tensors = 0; }
    if(active_tensors != 0) { delete active_tensors; active_tensors = 0; }
    if(active_weights != 0) { delete[] active_weights; active_weights = 0; }
    if(points != 0) { delete points; points = 0; }
    if(needed != 0) { delete needed; needed = 0; }
    if(exponents != 0) { delete exponents; exponents = 0; }
    if(values != 0) { delete values; values = 0; }
    if(fourier_coefs != 0) { delete[] fourier_coefs; fourier_coefs = 0; }
    num_dimensions = 0;
    num_outputs = 0;
}

void GridFourier::makeGrid(int cnum_dimensions, int cnum_outputs, int depth, TypeDepth type, const int* anisotropic_weights, const int* level_limits) {
    IndexManipulator IM(cnum_dimensions);

    IndexSet *tset = IM.selectTensors(depth, type, anisotropic_weights, rule_fourier);
    if (level_limits != 0){
        IndexSet *limited = IM.removeIndexesByLimit(tset, level_limits);
        if (limited != 0){
            delete tset;
            tset = limited;
        }
    }

    setTensors(tset, cnum_outputs);
}

void GridFourier::copyGrid(const GridFourier *fourier) {
    IndexSet *tset = new IndexSet(fourier->tensors);
    
    setTensors(tset, fourier->num_outputs);
    
    if ((num_outputs > 0) && (fourier->points != 0)){ // if there are values inside the source object
        loadNeededPoints(fourier->values->getValues(0));
    }
}

void GridFourier::setTensors(IndexSet* &tset, int cnum_outputs) {
    reset();
    num_dimensions = tset->getNumDimensions();
    num_outputs = cnum_outputs;
    
    tensors = tset;
    tset = 0;
    
    IndexManipulator IM(num_dimensions);
    
    OneDimensionalMeta meta(0);
    int* max_levels = new int[num_dimensions];
    int max_level; IM.getMaxLevels(tensors, max_levels, max_level);
    delete[] max_levels;
    wrapper = new OneDimensionalWrapper(&meta, max_level, rule_fourier);
    
    
    int* tensors_w = IM.makeTensorWeights(tensors);
    active_tensors = IM.nonzeroSubset(tensors, tensors_w);

    int nz_weights = active_tensors->getNumIndexes();

    active_weights = new int[nz_weights];
    tensor_refs = new int*[nz_weights];
    int count = 0;
    for(int i=0; i<tensors->getNumIndexes(); i++){ if (tensors_w[i] != 0) active_weights[count++] = tensors_w[i]; }

    delete[] tensors_w;
    
    needed = IM.generateNestedPoints(tensors, wrapper); // nested grids exploit nesting
    
    UnsortedIndexSet* exponents_unsorted = new UnsortedIndexSet(num_dimensions, needed->getNumIndexes());
    int *exponent = new int[num_dimensions];
    
    for (int i=0; i<needed->getNumIndexes(); i++) {
        for(int j=0; j<needed->getNumDimensions(); j++) {
            exponent[j] = (needed->getIndex(i)[j] % 2 == 0 ? -needed->getIndex(i)[j]/2 : (needed->getIndex(i)[j]+1)/2);
        }
        exponents_unsorted->addIndex(exponent); 
    }
    delete[] exponent;
    
    exponents = new IndexSet(exponents_unsorted);

    exponent_refs = new int*[nz_weights];
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<nz_weights; i++){
        exponent_refs[i] = referenceExponents(active_tensors->getIndex(i), exponents);
        tensor_refs[i] = IM.referenceNestedPoints(active_tensors->getIndex(i), wrapper, needed);
    }
    
    if (num_outputs == 0) {
        points = needed;
        needed = 0;
    } else {
        values = new StorageSet(num_outputs, needed->getNumIndexes());
    }
    
}

int* GridFourier::referenceExponents(const int levels[], const IndexSet *list) {
    
    // This is like IndexManipulator::referenceNestedPoints, but it references an ordered list of exponents for basis functions
    int *num_points = new int[num_dimensions];
    int num_total = 1;
    for(int j=0; j<num_dimensions; j++){  num_points[j] = wrapper->getNumPoints(levels[j]); num_total *= num_points[j];  }

    int* refs = new int[num_total];
    int *p = new int[num_dimensions];

    for(int i=0; i<num_total; i++){
        int t = i;
        for(int j=num_dimensions-1; j>=0; j--){
            int tmp = t % num_points[j];
            if(tmp <= (num_points[j]-1)/2) {
                p[j] = tmp;
            } else {
                p[j] = -num_points[j] + tmp;
            }
            t /= num_points[j];
        }
        refs[i] = list->getSlot(p);
    }

    delete[] p;
    delete[] num_points;
    
    return refs;
}

int GridFourier::getNumDimensions() const{ return num_dimensions; }
int GridFourier::getNumOutputs() const{ return num_outputs; }
TypeOneDRule GridFourier::getRule() const{ return rule_fourier; }

int GridFourier::getNumLoaded() const{ return (((points == 0) || (num_outputs == 0)) ? 0 : points->getNumIndexes()); }
int GridFourier::getNumNeeded() const{ return ((needed == 0) ? 0 : needed->getNumIndexes()); }
int GridFourier::getNumPoints() const{ return ((points == 0) ? getNumNeeded() : points->getNumIndexes()); }


void GridFourier::loadNeededPoints(const double *vals, TypeAcceleration acc) {
    if (accel != 0) accel->resetGPULoadedData();

    if(points == 0) { //setting points for the first time
        values->setValues(vals);
        points = needed;
        needed = 0;
    } else { //resetting the points
        values->setValues(vals);
    }
    //if we add anisotropic or surplus refinement, I'll need to add a third case here
    
    calculateFourierCoefficients();
}



double* GridFourier::getLoadedPoints() const{
    if (points == 0) return 0;
    int num_points = points->getNumIndexes();
    double *x = new double[num_dimensions * num_points];
    getLoadedPoints(x);
    return x;
}

void GridFourier::getLoadedPoints(double *x) const{
    int num_points = points->getNumIndexes();
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = points->getIndex(i);
        for(int j=0; j<num_dimensions; j++){
            x[i*num_dimensions + j] = wrapper->getNode(p[j]);
        }
    }
}

double* GridFourier::getNeededPoints() const{
    if (needed == 0) return 0;
    int num_points = needed->getNumIndexes();
    if (num_points == 0) return 0;
    double *x = new double[num_dimensions * num_points];
    getNeededPoints(x);
    return x;
}

void GridFourier::getNeededPoints(double *x) const{
    int num_points = needed->getNumIndexes();
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = needed->getIndex(i);
        for(int j=0; j<num_dimensions; j++){
            x[i*num_dimensions + j] = wrapper->getNode(p[j]);
        }
    }
}

double* GridFourier::getPoints() const{
    return ((points == 0) ? getNeededPoints() : getLoadedPoints());
}

void GridFourier::getPoints(double *x) const{
    if (points == 0){ getNeededPoints(x); }else{ getLoadedPoints(x); };
}

void GridFourier::calculateFourierCoefficients() {
    int num_nodes = getNumPoints();
    
    if(fourier_coefs != 0) { delete[] fourier_coefs; fourier_coefs = 0; }
    fourier_coefs = new std::complex<double>[num_outputs * num_nodes];
    std::fill(fourier_coefs, fourier_coefs + num_outputs*num_nodes, std::complex<double>(0,0));
    
    
    for(int k=0; k<num_outputs; k++) {
        for(int n=0; n<active_tensors->getNumIndexes(); n++) {  
            const int* levels = active_tensors->getIndex(n);
            int num_tensor_points = 1;
            int* num_oned_points = new int[num_dimensions];
            for(int j=0; j<num_dimensions; j++) {
                num_oned_points[j] = wrapper->getNumPoints(levels[j]);
                num_tensor_points *= num_oned_points[j];
            }
         
            // First we need to reorder from Tasmanian's internal indexing to get the point-grid in order of actual values
            // EXAMPLE: level-3 1D grid is [0, 1, 2, ..., 9] internally, but these indices correspond to
            //    [0, 1/3, 2/3, 1/9, 2/9, 4/9, 5/9, 7/9, 8/9] externally. We say external indexing is arranged
            //    in order of increasing x values.
         
         
            double* in = new double[num_tensor_points];
            std::complex<double>* out = new std::complex<double>[num_tensor_points];
         
            for(int i=0; i<num_tensor_points; i++) {
          
                // This i is external indexing, flattened according to standard C (row-major) indexing
                int t=i;
                int* p=new int[num_dimensions];
                for(int j=num_dimensions-1; j>=0; j--) {
                    p[j] = t % num_oned_points[j];
                    t /= num_oned_points[j];
                }
                // p[] stores the external indexing as a tensor address

                // Now we move from external to internal indexing
                int *p_internal = new int[num_dimensions];
                for(int j=0; j<num_dimensions; j++) {
                    int tmp = p[j]; 
                    if (tmp == 0) { 
                        p_internal[j] = 0; 
                    } else {
                        int go_back = 0;
                        while(tmp % 3 == 0 && tmp != 0) {
                            go_back += 1;
                            tmp /= 3;
                        }
                        
                        int num_prev_level = wrapper->getNumPoints(levels[j]-go_back-1);
                        p_internal[j] = num_prev_level + (tmp/3)*2 + tmp % 3 - 1;   
                    }
                }
          
                int fval_index = points->getSlot(p_internal);
                const double *v = values->getValues(fval_index);
                in[i] = v[k];
            }
         
         
            // Execute FFT
            TasmanianFourierTransform::dft(num_dimensions, num_oned_points, in, out);
        
            for(int i=0; i<num_tensor_points; i++) {
                fourier_coefs[k*num_nodes + exponent_refs[n][i]] += ((double) active_weights[n]) * out[i];    
                // Combine with tensor weights
            }
        }    
    }
}

std::complex<double>* GridFourier::getBasisFunctions(const double x[]) const {
    std::complex<double> unit_imag(0,1);    // this is sqrt(-1)
    std::complex<double> *weights = new std::complex<double>[exponents->getNumIndexes()];
    
    for(int i=0; i<exponents->getNumIndexes(); i++) {
        weights[i] = exp(2 * M_PI * unit_imag * ((double) exponents->getIndex(i)[0]) * x[0]);
        for(int j=1; j<exponents->getNumDimensions(); j++) {
            weights[i] *= exp(2 * M_PI * unit_imag * ((double) exponents->getIndex(i)[j]) * x[j]);
        }
    }
    return weights;
}

void GridFourier::getBasisFunctions(const double x[], double weights[]) const {
    
    // weights has length 2*num_nodes here
    std::complex<double> *tmp = getBasisFunctions(x);
    for(int i=0; i<exponents->getNumIndexes(); i++) {
        weights[2*i] = tmp[i].real();
        weights[2*i+1] = tmp[i].imag();
    }
    delete[] tmp;
}

void GridFourier::getBasisFunctions(const double x[], std::complex<double> weights[]) const {
    
    // weights has length num_nodes here
    weights = getBasisFunctions(x);
}

double* GridFourier::getInterpolationWeights(const double x[]) const { double *tmp = new double; *tmp = 0.0; return tmp; }
void GridFourier::getInterpolationWeights(const double x[], double weights[]) const { 
    std::fill(weights, weights + exponents->getNumIndexes(), 0.0); 
}

double* GridFourier::getQuadratureWeights() const {
    int num_points = getNumPoints();
    double *w = new double[num_points]; 
    getQuadratureWeights(w);
    return w;
}

void GridFourier::getQuadratureWeights(double weights[]) const{
    
    /* 
     * When integrating the Fourier series on a tensored grid, all the 
     * nonzero modes vanish, and we're left with the normalized Fourier 
     * coeff for e^0 (sum of the data divided by number of points)
     */
    
    int num_points = getNumPoints();
    std::fill(weights, weights+num_points, 0.0);
    for(int n=0; n<active_tensors->getNumIndexes(); n++) {
        const int *levels = active_tensors->getIndex(n);
        int num_tensor_points = 1;
        for(int j=0; j<num_dimensions; j++) {
            num_tensor_points *= wrapper->getNumPoints(levels[j]);
        }
        for(int i=0; i<num_tensor_points; i++) {
            weights[tensor_refs[n][i]] += active_weights[n]/num_tensor_points;
        }
        delete[] levels;
    }
}

void GridFourier::evaluate(const double x[], double y[]) const{
    std::complex<double> *w = getBasisFunctions(x);
    TasBLAS::setzero(num_outputs, y);
    for(int k=0; k<num_outputs; k++){
        for(int i=0; i<points->getNumIndexes(); i++){
            y[k] += (w[i] * fourier_coefs[k*(points->getNumIndexes())+i]).real();
        }
    }
    delete[] w;
}

void GridFourier::evaluateBatch(const double x[], int num_x, double y[]) const{
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        evaluate(&(x[((size_t) i) * ((size_t) num_dimensions)]), &(y[((size_t) i) * ((size_t) num_outputs)]));
    }
}

void GridFourier::evaluateFastCPUblas(const double x[], double y[]) const{
    #ifdef Tasmanian_ENABLE_BLAS
    std::cout<<"BLAS enabled"<<std::endl;
    std::complex<double> *w = getBasisFunctions(x);
    std::complex<double> *y_tmp = new std::complex<double>[num_outputs];
    TasBLAS::zgemv(num_outputs, points->getNumIndexes(), fourier_coefs, w, y_tmp);
    
    #pragma omp parallel for
    for(int i=0; i<num_outputs; i++) {
        y[i] = y_tmp[i].real();
    }
    
    delete[] w;
    delete[] y_tmp;
    #else
    std::cout<<"BLAS disabled"<<std::endl;
    evaluate(x,y);
    #endif // Tasmanian_ENABLE_BLAS
}

void GridFourier::evaluateFastGPUcublas(const double x[], double y[], std::ostream *os) const{
    evaluateFastCPUblas(x,y);
}

void GridFourier::evaluateFastGPUcuda(const double x[], double y[], std::ostream *os) const{
    evaluateFastCPUblas(x,y);
}

void GridFourier::evaluateFastGPUmagma(const double x[], double y[], std::ostream *os) const{
    evaluateFastCPUblas(x,y);
}

void GridFourier::evaluateBatchCPUblas(const double x[], int num_x, double y[]) const{
    #ifdef Tasmanian_ENABLE_BLAS
    int num_points = points->getNumIndexes();
    std::complex<double> *y_tmp = new std::complex<double>[num_outputs * num_x];
    double *weights = new double[num_points * num_x];
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        getBasisFunctions(&(x[((size_t) i) * ((size_t) num_dimensions)]), &(weights[((size_t) i) * ((size_t) num_points)]));
    }

    TasBLAS::zgemm(num_outputs, num_x, num_points, std::complex<double>(1.0,0.0), fourier_coefs, weights, std::complex<double>(0.0,0.0), y_tmp);

    #pragma omp parallel for
    for(int i=0; i<num_outputs*num_x; i++) {
        y[i] = y_tmp[i].real();
    }
    
    delete[] y_tmp;
    delete[] weights;
    #else
    evaluateBatch(x, num_x, y);
    #endif // Tasmanian_ENABLE_BLAS
}


void GridFourier::evaluateBatchGPUcublas(const double x[], int num_x, double y[], std::ostream *os) const {
    evaluateBatchCPUblas(x, num_x, y);
}

void GridFourier::evaluateBatchGPUcuda(const double x[], int num_x, double y[], std::ostream *os) const {
    evaluateBatchCPUblas(x, num_x, y);
}

void GridFourier::evaluateBatchGPUmagma(const double x[], int num_x, double y[], std::ostream *os) const {
    evaluateBatchCPUblas(x, num_x, y);
}

void GridFourier::integrate(double q[], double *conformal_correction) const{
    double *w = getQuadratureWeights();
    if(conformal_correction != 0) { for(int i=0; i<points->getNumIndexes(); i++) w[i] *= conformal_correction[i]; }
    std::fill(q, q+num_outputs, 0.0);
    #pragma omp parallel for schedule(static)
    for(int k=0; k<num_outputs; k++) {
        for(int i=0; i<points->getNumIndexes(); i++) {
            const double *v = values->getValues(i);
            q[k] += w[i] * v[k];
        }
    }
    delete[] w;
}

void GridFourier::evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const{
    // y must be of size num_x * num_nodes * 2
    int num_points = getNumPoints();
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<num_x; i++) {
        getBasisFunctions(&(x[((size_t) i) * ((size_t) num_dimensions)]), &(y[((size_t) i) * ((size_t) num_points)]));
    }
}

void GridFourier::setHierarchicalCoefficients(const double c[], TypeAcceleration acc, std::ostream *os) {
    
    // takes c to be length 2*num_outputs*num_nodes
    // first two entries are real and imag parts of the Fourier coef for the first basis function and first output dimension
    // second two entries are real/imag parts of the Fourier coef for the second basis function and first output dim
    // and so on
    if(accel != 0) accel->resetGPULoadedData();
    for(int i=0; i<num_outputs*getNumPoints(); i++) {
        fourier_coefs[i] = std::complex<double>(c[2*i], c[2*i+1]);
    }
}

void GridFourier::clearAccelerationData() {
    if(accel != 0) {
        delete accel;
        accel = 0;
    }
}

void GridFourier::clearRefinement() { return; }     // to be expanded later
void GridFourier::mergeRefinement() { return; }     // to be expanded later

const int* GridFourier::getPointIndexes() const{
    return ((points == 0) ? needed->getIndex(0) : points->getIndex(0));
}

const IndexSet* GridFourier::getExponents() const{
    return exponents;
}

const double* GridFourier::getFourierCoefs() const{
    double* fc = new double[2*(exponents->getNumIndexes())];
    for(int i=0; i<exponents->getNumIndexes(); i++) {
        fc[2*i] = fourier_coefs[i].real();
        fc[2*i+1] = fourier_coefs[i].imag();
    }
    return fc;
}


} // end TasGrid


#endif
