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

#ifndef __TASMANIAN_ADDONS_CMANAGER_HPP
#define __TASMANIAN_ADDONS_CMANAGER_HPP

/*!
 * \internal
 * \file tsgCandidateManager.hpp
 * \brief Manager for manipulations of candidate construction points.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianAddonsCommon
 *
 * A class that manipulate candidate construction points.
 * \endinternal
 */

#include "tsgAddonsCommon.hpp"

/*!
 * \ingroup TasmanianAddons
 * \addtogroup TasmanianAddonsConstruct Automated Surrogate Construction Procedure
 *
 */

namespace TasGrid{

/*!
 * \internal
 * \ingroup TasmanianAddonsConstruct
 * \brief Manages candidate points.
 *
 * Class for book keeping of candidate construction points,
 * e.g., started jobs, finished jobs, available jobs.
 * \endinternal
 */
class CandidateManager{
protected:
    enum TypeStatus{ free, running, done };
public:
    //! \brief Constructor, accepts number of dimensions as a constant parameter.
    template<typename IntType>
    CandidateManager(IntType dimensions) : num_dimensions((size_t) dimensions),
        num_candidates(0), num_running(0), num_done(0){}
    //! \brief Default destructor.
    ~CandidateManager(){}

    /*!
     * \brief Assign a new set of candidate points.
     *
     * The new candidates are provided in the order of importance,
     * but a lexicographical sorted list will also be created for searching purposes.
     * The new set of points is also assumed to contain the currently running jobs.
     */
    void operator=(std::vector<double> &&new_candidates){
        candidates = std::move(new_candidates);
        num_candidates = candidates.size() / num_dimensions;
        num_done = 0;

        sort_candidates();
        status.resize(num_candidates);
        std::fill(status.begin(), status.end(), free);
        for(auto const &p : running_jobs){
            size_t i = find(p);
            if (i < num_candidates) status[sorted[i]] = running;
        }
    }

    //! \brief Mark a point as "complete".
    void complete(std::vector<double> const &p){
        num_done++;
        num_running--;
        auto i = find(p);
        status[sorted[i]] = done;

        for(auto ip = running_jobs.before_begin(), it = running_jobs.begin();
            it != running_jobs.end();
            ip++){
            if (match(p, *it)){
                running_jobs.erase_after(ip);
                it = running_jobs.end(); // exit the for-loop
            }else{
                it++;
            }
        }
    }

    //! \brief Returns the next best point to compute, returns empty vector if no points are available.
    std::vector<double> next(){
        size_t i = 0;
        while((i < num_candidates) && (status[i] != free)) i++;
        if (i == num_candidates) return std::vector<double>();
        num_running++;
        status[i] = running;
        running_jobs.push_front(std::vector<double>(&candidates[i*num_dimensions], &candidates[i*num_dimensions] + num_dimensions));
        return running_jobs.front();
    }

    //! \brief Returns the number of running jobs.
    size_t getNumRunning() const{ return num_running; }

    //! \brief Returns the number of complete jobs.
    size_t getNumDone() const{ return num_done; }

    //! \brief Returns the number of all candidate jobs.
    size_t getNumCandidates() const{ return num_candidates; }

protected:
    //! \brief Returns \b true if entries in \b a and \b b match to \b Maths::num_tol, assumes sizes match already.
    bool match(std::vector<double> const &a, std::vector<double> const &b) const{
        for(auto ia = a.begin(), ib = b.begin(); ia != a.end(); ia++, ib++)
            if (std::abs(*ia - *ib) > Maths::num_tol) return false;
        return true;
    }

    //! \brief Returns \b true if the lexicographical order of \b a is before \b b.
    bool compare(double const a[], double const b[]) const{
        for(size_t i=0; i<num_dimensions; i++){
            if (a[i] < b[i] - Maths::num_tol) return true;
            if (a[i] > b[i] + Maths::num_tol) return false;
        }
        return false;
    }

    //! \brief Creates a sorted list of all candidates.
    void sort_candidates(){
        sorted.resize(num_candidates);
        std::iota(sorted.begin(), sorted.end(), 0);
        std::sort(sorted.begin(), sorted.end(), [&](size_t a, size_t b)->bool{ return compare(&candidates[a*num_dimensions], &candidates[b*num_dimensions]); });
    }

    //! \brief Find the index of the \b point within the \b canidates vector, returns \b num_candidates if not found.
    size_t find(std::vector<double> const &point){
        size_t sstart = 0, send = num_candidates - 1, current = (sstart + send) / 2;
        while (sstart <= send){
            const double *cp = &candidates[sorted[current]*num_dimensions];
            if (compare(cp, point.data())){
                sstart = current + 1;
            }else if (compare(point.data(), cp)){
                send = current - 1;
            }else{
                return current;
            }
            current = (sstart + send) / 2;
        }
        return num_candidates;
    }

private:
    size_t const num_dimensions;
    size_t num_candidates, num_running, num_done;
    std::vector<double> candidates;
    std::vector<size_t> sorted;
    std::vector<TypeStatus> status;
    std::forward_list<std::vector<double>> running_jobs;
};

/*!
 * \internal
 * \ingroup TasmanianAddonsConstruct
 * \brief Stores complete set of points before adding to the sparse grid.
 *
 * Constructing a grid one point at a time has high computational cost of at least O(N^2),
 * adding points in batches can amortize the cost to O(N log(N)).
 * This class stores the points temporarily before they can be loaded.
 * \endinternal
 */
class CompleteStorage{
public:
    //! \brief Constructor, accepts number of dimensions as a constant parameter.
    template<typename IntType>
    CompleteStorage(IntType dimensions) : num_dimensions((size_t) dimensions){}
    //! \brief Default destructor.
    ~CompleteStorage(){}

    //! \brief Add a point to the stored list.
    void add(std::vector<double> const &x, std::vector<double> const &y){
        points.insert(points.end(), x.begin(), x.end());
        values.insert(values.end(), y.begin(), y.end());
    }

    //! \brief Move the stored points into the grid.
    void load(TasmanianSparseGrid &grid){
        if (points.empty()) return;
        grid.loadConstructedPoints(points, values);
        points.clear();
        values.clear();
    }

    //! \brief Returns the number of stored points.
    size_t getNumStored() const{ return points.size() / num_dimensions; }

private:
    size_t const num_dimensions;
    std::vector<double> points, values;
};

}

#endif
