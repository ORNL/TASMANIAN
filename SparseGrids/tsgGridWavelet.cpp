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

#ifndef __TASMANIAN_SPARSE_GRID_WAVELET_CPP
#define __TASMANIAN_SPARSE_GRID_WAVELET_CPP

#include "tsgGridWavelet.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace TasGrid{

GridWavelet::GridWavelet() : num_dimensions(0), num_outputs(0), order(1), coefficients(0), points(0), needed(0), values(0), inter_matrix(0)
{}
GridWavelet::GridWavelet(const GridWavelet &wav) : num_dimensions(0), num_outputs(0), order(1), coefficients(0), points(0), needed(0), values(0), inter_matrix(0)
{  copyGrid(&wav);  }
GridWavelet::~GridWavelet(){ reset(); }

void GridWavelet::reset(){
    if (coefficients != 0){ delete[] coefficients; coefficients = 0; }
    if (points != 0){ delete points; points = 0; }
    if (needed != 0){ delete needed; needed = 0; }
    if (values != 0){ delete values; values = 0; }
    if (inter_matrix != 0){ delete inter_matrix; inter_matrix = 0; }
}

void GridWavelet::write(std::ofstream &ofs) const{
    ofs << std::scientific; ofs.precision(17);
    ofs << num_dimensions << " " << num_outputs << " " << order << endl;
    if (num_dimensions > 0){
        if (points == 0){
            ofs << "0" << endl;
        }else{
            ofs << "1 ";
            points->write(ofs);
        }
        if (coefficients == 0){
            ofs << "0" << endl;
        }else{
            ofs << "1 ";
            for(size_t i=0; i<((size_t) points->getNumIndexes()) * ((size_t) num_outputs); i++){ ofs << " " << coefficients[i]; } ofs << endl;
        }
        if (needed == 0){
            ofs << "0" << endl;
        }else{
            ofs << "1 ";
            needed->write(ofs);
        }
        if (num_outputs > 0) values->write(ofs);
    }
}
void GridWavelet::writeBinary(std::ofstream &ofs) const{
    int dims[3];
    dims[0] = num_dimensions;
    dims[1] = num_outputs;
    dims[2] = order;
    ofs.write((char*) dims, 3 * sizeof(int));
    if (num_dimensions > 0){
        char flag;
        if (points == 0){
            flag = 'n'; ofs.write(&flag, sizeof(char));
        }else{
            flag = 'y'; ofs.write(&flag, sizeof(char));
            points->writeBinary(ofs);
        }
        if (needed == 0){
            flag = 'n'; ofs.write(&flag, sizeof(char));
        }else{
            flag = 'y'; ofs.write(&flag, sizeof(char));
            needed->writeBinary(ofs);
        }
        if (coefficients == 0){
            flag = 'n'; ofs.write(&flag, sizeof(char));
        }else{
            flag = 'y'; ofs.write(&flag, sizeof(char));
            ofs.write((char*) coefficients, ((size_t) num_outputs) * ((size_t) points->getNumIndexes()) * sizeof(double));
        }
        if (num_outputs > 0) values->writeBinary(ofs);
    }
}
void GridWavelet::read(std::ifstream &ifs){
    reset();
    ifs >> num_dimensions >> num_outputs >> order;
    if (num_dimensions > 0){
        int flag;
        rule1D.updateOrder(order);

        ifs >> flag;
        if (flag == 1){
            points = new IndexSet(num_dimensions);
            points->read(ifs);
        }
        ifs >> flag;
        if (flag == 1){
            size_t size_coeff = ((size_t) points->getNumIndexes()) * ((size_t) num_outputs);
            coefficients = new double[size_coeff];
            for(size_t i=0; i<size_coeff; i++) ifs >> coefficients[i];
        }
        ifs >> flag;
        if (flag == 1){
            needed = new IndexSet(num_dimensions);
            needed->read(ifs);
        }

        if (num_outputs > 0){  values = new StorageSet(0, 0); values->read(ifs);  }
    }
    buildInterpolationMatrix();
}
void GridWavelet::readBinary(std::ifstream &ifs){
    int dims[3];
    ifs.read((char*) dims, 3 * sizeof(int));
    num_dimensions = dims[0];
    num_outputs = dims[1];
    order = dims[2];
    if (num_dimensions > 0){
        char flag;
        rule1D.updateOrder(order);

        ifs.read((char*) &flag, sizeof(char)); if (flag == 'y'){ points = new IndexSet(num_dimensions); points->readBinary(ifs); }
        ifs.read((char*) &flag, sizeof(char)); if (flag == 'y'){ needed = new IndexSet(num_dimensions); needed->readBinary(ifs); }

        ifs.read((char*) &flag, sizeof(char));
        if (flag == 'y'){
            size_t size_coeff = ((size_t) num_outputs) * ((size_t) points->getNumIndexes());
            coefficients = new double[size_coeff];
            ifs.read((char*) coefficients, size_coeff * sizeof(double));
        }

        if (num_outputs > 0){ values = new StorageSet(0, 0); values->readBinary(ifs); }
    }
    buildInterpolationMatrix();
}

void GridWavelet::makeGrid(int cnum_dimensions, int cnum_outputs, int depth, int corder, const int *level_limits){
    reset();
    num_dimensions = cnum_dimensions;
    num_outputs = cnum_outputs;
    order = corder;

    rule1D.updateOrder(order);

    IndexManipulator IM(num_dimensions);
    UnsortedIndexSet* deltas = IM.getToalDegreeDeltas(depth);
    if (level_limits != 0){
        UnsortedIndexSet *limited = IM.removeIndexesByLimit(deltas, level_limits);
        if (limited != 0){
            delete deltas;
            deltas = limited;
        }
    }

    WaveletLevels wavelevels(order); // describes the wavelet hierarchy, namely points per delta

    needed = IM.generatePointsFromDeltas(deltas, &wavelevels);
    delete deltas;

    if (num_outputs == 0){
        points = needed;
        needed = 0;
    }else{
        values = new StorageSet(num_outputs, needed->getNumIndexes());
    }

    buildInterpolationMatrix();
}
void GridWavelet::copyGrid(const GridWavelet *wav){
    reset();
    num_dimensions = wav->num_dimensions;
    num_outputs = wav->num_outputs;
    order = wav->order;

    rule1D.updateOrder(order);

    if (wav->points != 0) points = new IndexSet(wav->points);
    if (wav->needed != 0) needed = new IndexSet(wav->needed);

    if (wav->values != 0) values = new StorageSet(wav->values);

    buildInterpolationMatrix();

    if ((points != 0) && (num_outputs > 0)){ // points are loaded
        recomputeCoefficients();
    }
}

void GridWavelet::setNodes(IndexSet* &nodes, int cnum_outputs, int corder){
    reset();
    num_dimensions = nodes->getNumDimensions();
    num_outputs = cnum_outputs;
    order = corder;

    rule1D.updateOrder(order);

    needed = nodes; nodes = 0;

    if (num_outputs == 0){
        points = needed;
        needed = 0;
    }else{
        values = new StorageSet(num_outputs, needed->getNumIndexes());
    }

    buildInterpolationMatrix();
}

int GridWavelet::getNumDimensions() const{ return num_dimensions;  }
int GridWavelet::getNumOutputs() const{ return num_outputs;  }
TypeOneDRule GridWavelet::getRule() const{ return rule_wavelet;  }
int GridWavelet::getOrder() const{ return order;  }

int GridWavelet::getNumLoaded() const{ return (((points == 0) || (num_outputs == 0)) ? 0 : points->getNumIndexes()); }
int GridWavelet::getNumNeeded() const{ return ((needed == 0) ? 0 : needed->getNumIndexes()); }
int GridWavelet::getNumPoints() const{ return ((points == 0) ? getNumNeeded() : points->getNumIndexes()); }

double* GridWavelet::getLoadedPoints() const{
    if (points == 0) return 0;
    int num_points = points->getNumIndexes();
    if (num_points == 0) return 0;
    double *x = new double[num_dimensions * num_points];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = points->getIndex(i);
        for(int j=0; j<num_dimensions; j++){
            x[i*num_dimensions + j] = rule1D.getNode(p[j]);
        }
    }
    return x;
}
void GridWavelet::getLoadedPoints(double *x) const{
    int num_points = points->getNumIndexes();
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = points->getIndex(i);
        for(int j=0; j<num_dimensions; j++){
            x[i*num_dimensions + j] = rule1D.getNode(p[j]);
        }
    }
}
double* GridWavelet::getNeededPoints() const{
    if (needed == 0) return 0;
    int num_points = needed->getNumIndexes();
    if (num_points == 0) return 0;
    double *x = new double[num_dimensions * num_points];
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = needed->getIndex(i);
        for(int j=0; j<num_dimensions; j++){
            x[i*num_dimensions + j] = rule1D.getNode(p[j]);
        }
    }
    return x;
}
void GridWavelet::getNeededPoints(double *x) const{
    int num_points = needed->getNumIndexes();
    #pragma omp parallel for schedule(static)
    for(int i=0; i<num_points; i++){
        const int *p = needed->getIndex(i);
        for(int j=0; j<num_dimensions; j++){
            x[i*num_dimensions + j] = rule1D.getNode(p[j]);
        }
    }
}
double* GridWavelet::getPoints() const{
    return ((points == 0) ? getNeededPoints() : getLoadedPoints());
}
void GridWavelet::getPoints(double *x) const{
    if (points == 0){ getNeededPoints(x); }else{ getLoadedPoints(x); }
}

double* GridWavelet::getQuadratureWeights() const{
    IndexSet *work = (points == 0) ? needed : points;
    int num_points = work->getNumIndexes();
	double *weights = new double[num_points];
	getQuadratureWeights(weights);
	return weights;
}
void GridWavelet::getQuadratureWeights(double *weights) const{
    IndexSet *work = (points == 0) ? needed : points;
    int num_points = work->getNumIndexes();
	#pragma omp parallel for
	for(int i=0; i<num_points; i++){
		weights[i] = evalIntegral(work->getIndex(i));
	}
	solveTransposed(weights);
}
double* GridWavelet::getInterpolationWeights(const double x[]) const{
    IndexSet *work = (points == 0) ? needed : points;
    int num_points = work->getNumIndexes();
	double *weights = new double[num_points];
	getInterpolationWeights(x, weights);
	return weights;
}
void GridWavelet::getInterpolationWeights(const double x[], double *weights) const{
    IndexSet *work = (points == 0) ? needed : points;
    int num_points = work->getNumIndexes();
	#pragma omp parallel for
	for(int i=0; i<num_points; i++){
        weights[i] = evalBasis(work->getIndex(i), x);
	}
	solveTransposed(weights);
}
void GridWavelet::loadNeededPoints(const double *vals, TypeAcceleration){
    if (points == 0){
        values->setValues(vals);
        points = needed;
        needed = 0;
    }else if (needed == 0){
        values->setValues(vals);
    }else{
        values->addValues(points, needed, vals);
        points->addIndexSet(needed);
        delete needed; needed = 0;
        buildInterpolationMatrix();
    }
    recomputeCoefficients();
}
void GridWavelet::mergeRefinement(){
    if (needed == 0) return; // nothing to do
    int num_all_points = getNumLoaded() + getNumNeeded();
    size_t size_vals = ((size_t) num_all_points) * ((size_t) num_outputs);
    double *vals = new double[size_vals];
    std::fill(vals, vals + size_vals, 0.0);
    values->setValuesPointer(vals, num_all_points);
    if (points == 0){
        points = needed;
        needed = 0;
    }else{
        points->addIndexSet(needed);
        delete needed; needed = 0;
        buildInterpolationMatrix();
    }
    if (coefficients != 0) delete[] coefficients;
    coefficients = new double[size_vals];
    std::fill(coefficients, coefficients + size_vals, 0.0);
}
void GridWavelet::evaluate(const double x[], double y[]) const{
    int num_points = points->getNumIndexes();
	double *basis_values = new double[num_points];
	#pragma omp parallel for
	for(int i=0; i<num_points; i++){
        basis_values[i] = evalBasis(points->getIndex(i), x);
	}
	for(int j=0; j<num_outputs; j++){
        double sum = 0.0;
        #pragma omp parallel for reduction(+ : sum)
        for(int i=0; i<num_points; i++){
            sum += basis_values[i] * coefficients[((size_t) i) * ((size_t) num_outputs) + j];
        }
        y[j] = sum;
	}
	delete[] basis_values;
}

void GridWavelet::evaluateFastCPUblas(const double x[], double y[]) const{ evaluate(x, y); }
void GridWavelet::evaluateFastGPUcublas(const double x[], double y[], std::ostream*) const{ evaluate(x, y); }
void GridWavelet::evaluateFastGPUcuda(const double x[], double y[], std::ostream*) const{ evaluate(x, y); }
void GridWavelet::evaluateFastGPUmagma(const double x[], double y[], std::ostream*) const{ evaluate(x, y); }

void GridWavelet::evaluateBatch(const double x[], int num_x, double y[]) const{
    #pragma omp parallel for
    for(int i=0; i<num_x; i++){
        evaluate(&(x[i*num_dimensions]), &(y[i*num_outputs]));
    }
}
void GridWavelet::evaluateBatchCPUblas(const double x[], int num_x, double y[]) const{
    evaluateBatch(x, num_x, y);
}
void GridWavelet::evaluateBatchGPUcublas(const double x[], int num_x, double y[], std::ostream*) const{
    evaluateBatch(x, num_x, y);
}
void GridWavelet::evaluateBatchGPUcuda(const double x[], int num_x, double y[], std::ostream*) const{
    evaluateBatch(x, num_x, y);
}
void GridWavelet::evaluateBatchGPUmagma(const double x[], int num_x, double y[], std::ostream*) const{
    evaluateBatch(x, num_x, y);
}

void GridWavelet::integrate(double q[], double *conformal_correction) const{
    int num_points = points->getNumIndexes();

    if (conformal_correction == 0){
        double *basis_integrals = new double[num_points];
        #pragma omp parallel for
        for(int i=0; i<num_points; i++){
            basis_integrals[i] = evalIntegral(points->getIndex(i));
        }
        for(int j=0; j<num_outputs; j++){
            double sum = 0.0;
            #pragma omp parallel for reduction(+ : sum)
            for(int i=0; i<num_points; i++){
                sum += basis_integrals[i] * coefficients[((size_t) i) * ((size_t) num_outputs) + j];
            }
            q[j] = sum;
        }
        delete[] basis_integrals;
    }else{
        std::fill(q, q + num_outputs, 0.0);
        double *w = getQuadratureWeights();
        for(int i=0; i<num_points; i++){
            w[i] *= conformal_correction[i];
            const double *vals = values->getValues(i);
            for(int k=0; k<num_outputs; k++){
                q[k] += w[i] * vals[k];
            }
        }
        delete[] w;
    }
}

double GridWavelet::evalBasis(const int p[], const double x[]) const{
	/*
	 * Evaluates the wavelet basis given at point p at the coordinates given by x.
	 */
	double v = 1.0;
	for(int i = 0; i < num_dimensions; i++){
		v *= rule1D.eval(p[i], x[i]);
		if (v == 0.0){ break; }; // MIRO: reduce the expensive wavelet evaluations
	}
	return v;
}
double GridWavelet::evalIntegral(const int p[]) const{
	/*
	 * For a given node p, evaluate the integral of the associated wavelet.
	 */
	double v = 1.0;
	for(int i = 0; i < num_dimensions; i++){
		v *= rule1D.getWeight(p[i]);
		if (v == 0.0){ break; }; // MIRO: reduce the expensive wavelet evaluations
	}
	return v;
}

void GridWavelet::buildInterpolationMatrix(){
    // change the logic here, build the matrix only if num_outputs are zero or if there is need (load surpluses etc)
	/*
	 * Using the IndexSet points, constructs the interpolation matrix needed for methods
	 * such as recomputeCoefficients, solveTransposed, etc.
	 */
	IndexSet *work = (points == 0) ? needed : points;
	if(inter_matrix != 0) { delete inter_matrix; }

	int num_points = work->getNumIndexes();

	TasSparse::TsgSparseCOO coo_mat(num_points, num_points);

#ifdef _OPENMP
	int num_threads = omp_get_max_threads();
	TasSparse::TsgSparseCOO *coos = new TasSparse::TsgSparseCOO[num_threads];
	#pragma omp parallel
#endif
	{
		double *xs = new double[num_dimensions];

    // This is O(n^2) but Wavelets don't have a nice rule of support to use monkeys and graphs!
#ifdef _OPENMP
		int thread_id = omp_get_thread_num();
		#pragma omp for
#endif
		for(int i = 0; i < num_points; i++){ /* Loop over points */
			const int *point = work->getIndex(i);

			for(int k = 0; k < num_dimensions; k++){ /* Fill in x values of point */
				xs[k] = rule1D.getNode(point[k]);
			}

			for(int j = 0; j < num_points; j++){ /* Loop over wavelets */
				const int *wavelet = work->getIndex(j);
				double v = 1.;

				for(int k = 0; k < num_dimensions; k++){ /* Loop over dimensions */
					v *= rule1D.eval(wavelet[k], xs[k]);
					if (v == 0.0){ break; }; // MIRO: evaluating the wavelets is expensive, stop if any one of them is zero
				} /* End for dimensions */

				if(v != 0.0){
				/*
				 * Testing for equal to zero is safe since v*0 is the only way
				 * for this to happen.
				 */
#ifdef _OPENMP
					coos[thread_id].addPoint(i, j, v);
#else
					coo_mat.addPoint(i,j,v);
#endif
				}


			} /* End for wavelets */

		} /* End for points */

		delete[] xs;

#ifdef _OPENMP
		#pragma omp master
		{ /* Master thread combines the sub-thread work */

			for(int i = 0; i < num_threads; i++){
				coo_mat.combine(coos[i]);
			}
			delete[] coos;
		} /* End master section */
#endif
	} /* End parallel section */
	inter_matrix = new TasSparse::SparseMatrix(coo_mat);
}

void GridWavelet::recomputeCoefficients(){
	/*
	 * Recalculates the coefficients to interpolate the values in points.
	 * Make sure buildInterpolationMatrix has been called since the list was updated.
	 */

	int num_points = points->getNumIndexes();
	if (coefficients != 0){ delete[] coefficients; }
	coefficients = new double[((size_t) num_points) * ((size_t) num_outputs)];

	if ((inter_matrix == 0) || (inter_matrix->getNumRows() != num_points)){
        buildInterpolationMatrix();
	}

	double *workspace = new double[2*num_points];
	double *b = workspace;
	double *x = workspace + num_points;

	for(int output = 0; output < num_outputs; output++){
		// Copy relevant portion
		std::fill(workspace, workspace + 2*num_points, 0.0);

		// Populate RHS
		for(int i = 0; i < num_points; i++){
			b[i] = values->getValues(i)[output];
		}

		// Solve system
		inter_matrix->solve(b, x);

		// Populate surplus
		for(int i = 0; i < num_points; i++){
			coefficients[((size_t) num_outputs) * ((size_t) i) + output] = x[i];
		}
	}

	delete[] workspace;
}

void GridWavelet::solveTransposed(double w[]) const{
	/*
	 * Solves the system A^T * w = y. Used to calculate interpolation and integration
	 * weights. RHS values should be passed in through w. At exit, w will contain the
	 * required weights.
	 */
	int num_points = inter_matrix->getNumRows();

	double *y = new double[num_points];

	std::copy(w, w + num_points, y);

	inter_matrix->solve(y, w, true);

	delete[] y;
}

double* GridWavelet::getNormalization() const{
    double* norm = new double[num_outputs];  std::fill(norm, norm + num_outputs, 0.0);
    for(int i=0; i<points->getNumIndexes(); i++){
        const double *v = values->getValues(i);
        for(int j=0; j<num_outputs; j++){
            if (norm[j] < fabs(v[j])) norm[j] = fabs(v[j]);
        }
    }
    return norm;
}
int* GridWavelet::buildUpdateMap(double tolerance, TypeRefinement criteria, int output) const{
    int num_points = points->getNumIndexes();
    int *pmap = new int[num_points * num_dimensions];  std::fill(pmap, pmap + num_points * num_dimensions, 0);

    double *norm = getNormalization();

    if ((criteria == refine_classic) || (criteria == refine_parents_first)){
        // classic refinement
        #pragma omp parallel for
        for(int i=0; i<num_points; i++){
            bool small = true;
            if (output == -1){
                for(size_t k=0; k<((size_t) num_outputs); k++){
                    if (small && ((fabs(coefficients[((size_t) i) * ((size_t) num_outputs) + k]) / norm[k]) > tolerance)) small = false;
                }
            }else{
                small = !((fabs(coefficients[((size_t) i) * ((size_t) num_outputs) + output]) / norm[output]) > tolerance);
            }
            if (!small){
                std::fill(&(pmap[i*num_dimensions]), &(pmap[i*num_dimensions]) + num_dimensions, 1);
            }
        }
    }else{
        SplitDirections split(points);

        for(int s=0; s<split.getNumJobs(); s++){
            int d = split.getJobDirection(s);
            int nump = split.getJobNumPoints(s);
            const int *pnts = split.getJobPoints(s);

            int *global_to_pnts = new int[num_points];

            int active_outputs = (output == -1) ? num_outputs : 1;

            double *vals = new double[nump * active_outputs];
            int *indexes = new int[nump * num_dimensions];

            for(int i=0; i<nump; i++){
                const double* v = values->getValues(pnts[i]);
                if (output == -1){
                    std::copy(v, v + num_outputs, &(vals[((size_t) i) * ((size_t) num_outputs)]));
                }else{
                    vals[i] = v[output];
                }
                const int *p = points->getIndex(pnts[i]);
                std::copy(p, p + num_dimensions, &(indexes[((size_t) i) * ((size_t) num_dimensions)]));
                global_to_pnts[pnts[i]] = i;
            }
            IndexSet *pointset = new IndexSet(num_dimensions, nump, indexes);

            GridWavelet direction_grid;
            direction_grid.setNodes(pointset, active_outputs, order);
            direction_grid.loadNeededPoints(vals, accel_none);
            const double *coeff = direction_grid.coefficients;

            for(size_t i=0; i<((size_t) nump); i++){
                bool small = true;
                if (output == -1){
                    for(size_t k=0; k<((size_t) num_outputs); k++){
                        if (small && ((fabs(coefficients[pnts[i]*num_outputs + k]) / norm[k]) > tolerance) && ((fabs(coeff[i*num_outputs + k]) / norm[k]) > tolerance)) small = false;
                    }
                }else{
                    if (((fabs(coefficients[pnts[i]*num_outputs + output]) / norm[output]) > tolerance) && ((fabs(vals[i]) / norm[output]) > tolerance)) small = false;
                }
                pmap[pnts[i]*num_dimensions + d] = (small) ? 0 : 1;;
            }

            delete[] vals;
            delete[] global_to_pnts;
        }
    }

    delete[] norm;

    return pmap;
}

bool GridWavelet::addParent(const int point[], int direction, GranulatedIndexSet *destination, IndexSet *exclude) const{
    int *dad = new int[num_dimensions];  std::copy(point, point + num_dimensions, dad);
    bool added = false;
    dad[direction] = rule1D.getParent(point[direction]);
    if (dad[direction] == -2){
        for(int c=0; c<rule1D.getNumPoints(0); c++){
            dad[direction] = c;
            if (exclude->getSlot(dad) == -1){
                destination->addIndex(dad);
                added = true;
            }
        }
    }else if (dad[direction] >= 0){
        if (exclude->getSlot(dad) == -1){
            destination->addIndex(dad);
            added = true;
        }
    }

    delete[] dad;
    return added;
}
void GridWavelet::addChild(const int point[], int direction, GranulatedIndexSet *destination, IndexSet *exclude) const{
    int *kid = new int[num_dimensions];  std::copy(point, point + num_dimensions, kid);
    int L, R; rule1D.getChildren(point[direction], L, R);
    kid[direction] = L;
    if ((kid[direction] != -1) && (exclude->getSlot(kid) == -1)){
        destination->addIndex(kid);
    }
    kid[direction] = R;
    if ((kid[direction] != -1) && (exclude->getSlot(kid) == -1)){
        destination->addIndex(kid);
    }
    delete[] kid;
}
void GridWavelet::addChildLimited(const int point[], int direction, GranulatedIndexSet *destination, IndexSet *exclude, const int *level_limits) const{
    int *kid = new int[num_dimensions];  std::copy(point, point + num_dimensions, kid);
    int L, R; rule1D.getChildren(point[direction], L, R);
    kid[direction] = L;
    if ((kid[direction] != -1) && (rule1D.getLevel(L) <= level_limits[direction]) && (exclude->getSlot(kid) == -1)){
        destination->addIndex(kid);
    }
    kid[direction] = R;
    if ( (kid[direction] != -1) && (rule1D.getLevel(R) <= level_limits[direction]) && (exclude->getSlot(kid) == -1)){
        destination->addIndex(kid);
    }
    delete[] kid;
}

void GridWavelet::clearRefinement(){
    if (needed != 0){ delete needed;  needed = 0; }
}
const double* GridWavelet::getSurpluses() const{
    return coefficients;
}

void GridWavelet::evaluateHierarchicalFunctions(const double x[], int num_x, double y[]) const{
    IndexSet *work = (points == 0) ? needed : points;
    int num_points = work->getNumIndexes();
    for(int i=0; i<num_x; i++){
        double *this_y = &(y[i*num_points]);
        const double *this_x = &(x[i*num_dimensions]);
        for(int j=0; j<num_points; j++){
            const int* p = work->getIndex(j);
            double v = 1.0;
            for(int k=0; k<num_dimensions; k++){
                v *= rule1D.eval(p[k], this_x[k]);
                if (v == 0.0){ break; }; // MIRO: evaluating the wavelets is expensive, stop if any one of them is zero
            }
            this_y[j] = v;
        }
    }
}

void GridWavelet::setHierarchicalCoefficients(const double c[], TypeAcceleration acc, std::ostream *os){
    double *vals = 0;
    bool aliased = false;
    if (points != 0){
        clearRefinement();
        vals = values->aliasValues();
        aliased = true;
    }else{
        points = needed;
        needed = 0;
        vals = new double[((size_t) points->getNumIndexes()) * ((size_t) num_outputs)];
    }
    int num_ponits = points->getNumIndexes();
    size_t size_coeff = ((size_t) num_ponits) * ((size_t) num_outputs);
    if (coefficients != 0) delete[] coefficients;
    coefficients = new double[size_coeff];
    std::copy(c, c + size_coeff, coefficients);
    double *x = getPoints();
    if (acc == accel_cpu_blas){
        evaluateBatchCPUblas(x, points->getNumIndexes(), vals);
    }else if (acc == accel_gpu_cublas){
        evaluateBatchGPUcublas(x, points->getNumIndexes(), vals, os);
    }else if (acc == accel_gpu_cuda){
        evaluateBatchGPUcuda(x, points->getNumIndexes(), vals, os);
    }else{
        evaluateBatch(x, points->getNumIndexes(), vals);
    }
    delete[] x;
    if (! aliased) values->setValuesPointer(vals, num_ponits);
}

const int* GridWavelet::getPointIndexes() const{
    return ((points == 0) ? needed->getIndex(0) : points->getIndex(0));
}
void GridWavelet::setSurplusRefinement(double tolerance, TypeRefinement criteria, int output, const int *level_limits){
    clearRefinement();

    int *pmap = buildUpdateMap(tolerance, criteria, output);

    bool useParents = (criteria == refine_fds) || (criteria == refine_parents_first);

    GranulatedIndexSet *refined = new GranulatedIndexSet(num_dimensions);

    int num_points = points->getNumIndexes();

    if (level_limits == 0){
        for(int i=0; i<num_points; i++){
            for(int j=0; j<num_dimensions; j++){
                if (pmap[i*num_dimensions+j] == 1){ // if this dimension needs to be refined
                    if (!(useParents && addParent(points->getIndex(i), j, refined, points))){
                        addChild(points->getIndex(i), j, refined, points);
                    }
                }
            }
        }
    }else{
        for(int i=0; i<num_points; i++){
            for(int j=0; j<num_dimensions; j++){
                if (pmap[i*num_dimensions+j] == 1){ // if this dimension needs to be refined
                    if (!(useParents && addParent(points->getIndex(i), j, refined, points))){
                        addChildLimited(points->getIndex(i), j, refined, points, level_limits);
                    }
                }
            }
        }
    }

    if (refined->getNumIndexes() > 0){
        needed = new IndexSet(refined);
    }
    delete refined;

    delete[] pmap;
}

void GridWavelet::clearAccelerationData(){}

}

#endif
