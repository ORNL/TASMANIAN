
#ifndef __TASMANIAN_SPARSE_GRID_FOURIER_HPP
#define __TASMANIAN_SPARSE_GRID_FOURIER_HPP

#include <cstdlib>
#include <math.h>
#include <complex>

#include "tsgEnumerates.hpp"
#include "tsgIndexSets.hpp"
#include "tsgCoreOneDimensional.hpp"
#include "tsgIndexManipulator.hpp"
#include "tsgLinearSolvers.hpp"
#include "tsgOneDimensionalWrapper.hpp"
#include "tsgGridCore.hpp"

#include "tsgAcceleratedDataStructures.hpp"

namespace TasGrid{

class GridFourier : public BaseCanonicalGrid {
public:
    GridFourier();
    GridFourier(const GridFourier &fourier);
    ~GridFourier();

    void makeGrid(int cnum_dimensions, int cnum_outputs, int depth, TypeDepth type, const int* anisotropic_weights = 0, const int* level_limits = 0);
    void copyGrid(const GridFourier *fourier);

    void setTensors(IndexSet* &tset, int cnum_outputs);
    int* referenceExponents(const int levels[], const IndexSet *list);

    int getNumDimensions() const;
    int getNumOutputs() const;
    TypeOneDRule getRule() const;

    int getNumLoaded() const;
    int getNumNeeded() const;
    int getNumPoints() const; // returns the number of loaded points unless no points are loaded, then returns the number of needed points
    
    void loadNeededPoints(const double *vals, TypeAcceleration acc = accel_none);
    
    double* getLoadedPoints() const;
    void getLoadedPoints(double *x) const;
    double* getNeededPoints() const;
    void getNeededPoints(double *x) const;
    double* getPoints() const;
    void getPoints(double *x) const; // returns the loaded points unless no points are loaded, then returns the needed points
    
    double* getInterpolationWeights(const double x[]) const;
    void getInterpolationWeights(const double x[], double weights[]) const;
    
    double* getQuadratureWeights() const;
    void getQuadratureWeights(double weights[]) const;

    void evaluate(const double x[], double y[]) const;
    void evaluateBatch(const double x[], int num_x, double y[]) const;
    
    void evaluateFastCPUblas(const double x[], double y[]) const;
    void evaluateFastGPUcublas(const double x[], double y[], std::ostream *os) const;
    void evaluateFastGPUcuda(const double x[], double y[], std::ostream *os) const;
    void evaluateFastGPUmagma(const double x[], double y[], std::ostream *os) const;
    
    void evaluateBatchCPUblas(const double x[], int num_x, double y[]) const;
    void evaluateBatchGPUcublas(const double x[], int num_x, double y[], std::ostream *os) const;
    void evaluateBatchGPUcuda(const double x[], int num_x, double y[], std::ostream *os) const;
    void evaluateBatchGPUmagma(const double x[], int num_x, double y[], std::ostream *os) const;
    
    void integrate(double q[], double *conformal_correction) const;
    
    void evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const;
    void setHierarchicalCoefficients(const double c[], TypeAcceleration acc, std::ostream *os);

    void clearAccelerationData();
    void clearRefinement();
    void mergeRefinement();

    const int* getPointIndexes() const;
    const IndexSet* getExponents() const;
    const double* getFourierCoefs() const;
    bool blasEnabled() const;

protected:
    void reset();
    void calculateFourierCoefficients();

    std::complex<double>* getBasisFunctions(const double x[]) const;    
    void getBasisFunctions(const double x[], std::complex<double> weights[]) const;
    void getBasisFunctions(const double x[], double weights[]) const;   // for evaluateHierarchicalFunctions (parallelized)
    
    int convertIndexes(const int i, const int levels[]) const;

private:
    int num_dimensions, num_outputs;

    OneDimensionalWrapper *wrapper;

    IndexSet *tensors;
    IndexSet *active_tensors;
    int *active_weights;
    IndexSet *points;
    IndexSet *needed;
    IndexSet *exponents;
    
    std::complex<double> *fourier_coefs;
    int **exponent_refs;
    int **tensor_refs;

    StorageSet *values;
    
    mutable BaseAccelerationData *accel;

};

}

#endif
