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
 * Templates that perform an automated surrogate construction to a model defined
 * by a lambda expression.
 * This allows for a one-click (or one function call) grid construction where
 * the samples can be computed either sequentially or in parallel.
 * The procedure is carried out until either a target accuracy is reached
 * or a computational budget is exhausted.
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
    //! \brief Jobs are divided into already checked-out (running), checked-in (done), and available (free).
    enum TypeStatus{ free, running, done };
public:
    //! \brief Constructor, accepts number of dimensions as a constant parameter.
    template<typename IntTypeDim, typename IntTypeBatch>
    CandidateManager(IntTypeDim dimensions, IntTypeBatch batch_size) : num_dimensions(static_cast<size_t>(dimensions)),
        num_batch(batch_size), num_candidates(0), num_running(0), num_done(0){}
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
        if (num_candidates == 0) return;

        sort_candidates();
        status.resize(num_candidates);
        std::fill(status.begin(), status.end(), free);
        for(auto const &p : running_jobs){
            size_t i = find(p.data());
            if (i < num_candidates) status[sorted[i]] = running;
        }
    }

    //! \brief Mark a point as "complete".
    void complete(std::vector<double> const &p){
        size_t num_complete = p.size() / num_dimensions;
        num_done += num_complete;
        num_running -= num_complete;

        for(auto ip = p.begin(); ip != p.end(); std::advance(ip, num_dimensions)){
            auto i = find(&*ip);
            // there is a scenario where a point is checked out
            // then while computing another point is done which changes the surpluses
            // and then the checked out point is no longer a candidate
            if (i < num_candidates) status[sorted[i]] = done;
        }

        // takes an iterator and returns an iterator to the next entry
        auto inext = [](std::forward_list<std::vector<double>>::iterator ib)->
            std::forward_list<std::vector<double>>::iterator{
                return ++ib;
            };

        for(size_t i=0; i<num_complete; i++){
            auto ib = running_jobs.before_begin();
            while(not match(&p[i * num_dimensions], inext(ib)->data())) ib++;
            running_jobs.erase_after(ib);
        }
    }

    //! \brief Returns the next best point to compute, returns empty vector if no points are available.
    std::vector<double> next(size_t remaining_budget){
        size_t this_batch = std::min(remaining_budget, num_batch);
        size_t i = 0;
        while((i < num_candidates) && (status[i] != free)) i++;
        if (i == num_candidates) return std::vector<double>();
        num_running++;
        status[i] = running;
        std::vector<double> result(&candidates[i*num_dimensions], &candidates[i*num_dimensions] + num_dimensions);
        running_jobs.push_front(result);

        size_t num_next = 1;
        while((i < num_candidates) && (num_next < this_batch)){
            while((i < num_candidates) && (status[i] != free)) i++;
            if (i < num_candidates){
                running_jobs.push_front(std::vector<double>(&candidates[i*num_dimensions], &candidates[i*num_dimensions] + num_dimensions));
                result.insert(result.end(), &candidates[i*num_dimensions], &candidates[i*num_dimensions] + num_dimensions);
                num_next++;
                num_running++;
                status[i++] = running;
            }
        }
        return result;
    }

    //! \brief Returns the number of running jobs.
    size_t getNumRunning() const{ return num_running; }

    //! \brief Returns the number of complete jobs.
    size_t getNumDone() const{ return num_done; }

    //! \brief Returns the number of all candidate jobs.
    size_t getNumCandidates() const{ return num_candidates; }

protected:
    //! \brief Returns \b true if entries in \b a and \b b match to \b Maths::num_tol, assumes sizes match already.
    bool match(double const a[], double const b[]) const{
        for(size_t i=0; i<num_dimensions; i++)
            if (std::abs(a[i] - b[i]) > Maths::num_tol) return false;
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
    size_t find(double const point[]) const{
        // if the point precedes the entire list, then send will craw back to -1
        // if send reaches -1 in unsigned arithmetic the method will segfault
        int sstart = 0, send = (int) num_candidates - 1, current = (sstart + send) / 2;
        while (sstart <= send){
            const double *cp = &candidates[sorted[current]*num_dimensions];
            if (compare(cp, point)){
                sstart = current + 1;
            }else if (compare(point, cp)){
                send = current - 1;
            }else{
                return (size_t) current;
            }
            current = (sstart + send) / 2;
        }
        return num_candidates;
    }

private:
    size_t const num_dimensions, num_batch;
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

    //! \brief Write the stored samples to a stream
    void write(std::ostream &os) const{
        IO::writeNumbers<mode_binary, IO::pad_auto>(os, points.size(), values.size());
        IO::writeVector<mode_binary, IO::pad_auto>(points, os);
        IO::writeVector<mode_binary, IO::pad_auto>(values, os);
    }

    //! \brief Read the stored samples from the stream
    void read(std::istream &is){
        points.resize(IO::readNumber<IO::mode_binary_type, size_t>(is));
        values.resize(IO::readNumber<IO::mode_binary_type, size_t>(is));
        IO::readVector<IO::mode_binary_type>(is, points);
        IO::readVector<IO::mode_binary_type>(is, values);
    }

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
