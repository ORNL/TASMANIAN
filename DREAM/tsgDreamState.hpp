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

#ifndef __TASMANIAN_DREAM_STATE_HPP
#define __TASMANIAN_DREAM_STATE_HPP

#include "tsgDreamEnumerates.hpp"

//! \internal
//! \file tsgDreamState.hpp
//! \brief The container class holding the DREAM history.
//! \author Miroslav Stoyanov
//! \ingroup TasmanianDREAM
//!
//! Defines the TasmanianDREAM class which stores the history of the MCMC samples.

/*!
 * \ingroup TasmanianDREAM
 * \addtogroup DREAMState TasmanianDREAM State
 *
 * Keeps the current state of the DREAM chains and stores the history
 * of a MCMC sampling.
 */

namespace TasDREAM{

//! \brief Contains the current state and the history of the DREAM chains.
//! \ingroup DREAMState

//! \par DREAM State
//! The DREAM state consists of \b num_chains vectors of size \b num_dimensions,
//! where the dimension is either specified or inferred from the given sparse grid.
//! Each vector is also associated with a
//! value of the corresponding probability density function (pdf).
//! The history of the state consists of multiple snapshots of the vectors and
//! the pdf values.
//!
//! \par Container class
//! The TasmanianDREAM class is a container for the current state and history,
//! where the vectors are stored using std::vector data structures.
//! Methods are provided for file I/O, direct access to the data and
//! several handy queries for inference or optimization problems.
class TasmanianDREAM{
public:
    //! \brief Constructor for a null DREAM state, no chains and no dimensions (used for MPI purposes).
    TasmanianDREAM();
    //! \brief Constructor for a DREAM state with the number of chains and dimensions.
    TasmanianDREAM(int cnum_chains, int cnum_dimensions);
    //! \brief Constructor for a DREAM state with the number of chains and the dimension of the sparse grid.
    TasmanianDREAM(int cnum_chains, const TasGrid::TasmanianSparseGrid &grid);
    //! \brief Default destructor, release all used memory.
    ~TasmanianDREAM();

    //! \brief Return human readable string with the version.
    static const char* getVersion(){ return TASMANIAN_VERSION_STRING; }
    //! \brief Return human readable sting with the license, refer to the LICENSE file for details.
    static const char* getLicense(){ return TASMANIAN_LICENSE; }
    //! \brief Return the major version of the library.
    static int getVersionMajor(){ return TASMANIAN_VERSION_MAJOR; }
    //! \brief Return the minor version of the library.
    static int getVersionMinor(){ return TASMANIAN_VERSION_MINOR; }

    //! \brief Return the number of dimensions.
    int getNumDimensions() const{ return (int) num_dimensions; }
    //! \brief Return the number of chains.
    int getNumChains() const{ return (int) num_chains; }
    //! \brief Return the number of saved vectors in the history.
    size_t getNumHistory() const{ return pdf_history.size(); }

    //! \brief Return \b true if the state has already been initialized with \b setState().
    bool isStateReady() const{ return init_state; }
    //! \brief Return \b true if the pdf values have already been initialized with \b setPDF().
    bool isPDFReady() const{ return init_values; }

    //! \brief Set the current DREAM state to the \b new_state.

    //! The \b new_state must have size equal to num_chains times num_dimensions, the first
    //! chain should be stored in the first num_dimensions entries, the second chain in the
    //! second batch of num_dimensions entries, and so on.
    void setState(const std::vector<double> &new_state);

    //! \brief Set the current DREAM state by calling \b update_state on each vector.

    //! The \b update_state function will be given each state vector as an array and
    //! the function should overwrite the entries with the desired current state.
    //! For example, \b update_state can overwrite the array with random numbers sampled
    //! according to a desired probability distribution.
    void setState(std::function<void(double *)> update_state);

    //! \brief Set the current set of pdf values to the \b new_values.

    //! The \b new_values must have size equal to the number of chains,
    //! the current set of pdf values are overwritten with the new ones.
    //! Note: that there is no restriction on the sign of the values,
    //! since log-form sampling can result in negative values.
    void setPDFvalues(const std::vector<double> &new_values);

    //! \brief Set the current set of pdf values using the given \b probability_distribution.

    //! Calls the \b probability_distribution function wit the current state (the state must be already initialized).
    //! The \b probability_distribution is the same function as in the \b DREAM sampling call,
    //! in order to avoid discrepancy, it may be best to leave the PDF values uninitialized and
    //! let the \b DREAM sampler to call \b setPDFvalues() automatically.
    void setPDFvalues(std::function<void(const std::vector<double> &state, std::vector<double> &values)> probability_distribution);

    //! \brief Erase the currently stored pdf values, useful when switching between log and non-log form sampling.
    void clearPDFvalues();

    //! \brief Used by the \b DREAM sampler, \f$ x = s_i + w ( s_k - s_j) \f$, where \b s indicates the current i, j, and k-th state vectors.

    //! Used by the \b DREAM sampler and probably should not be called by the user.
    //! This returns the vector of the \b i-th chain corrected by the difference
    //! between the \b k-th and \b j-th chain and weighted by \b w.
    //! The result is written to the array \b x, which must have size at least num_dimensions.
    void getIJKdelta(size_t i, size_t j, size_t k, double w, std::vector<double> &x) const;

    //! \brief Return the state of the \b i-th chain, store an array \b x of size at least num_dimensions.

    //! Used by the \b DREAM sampler and probably should not be called by the user.
    void getChainState(size_t i, double *x) const{ std::copy_n(state.begin() + i * num_dimensions, num_dimensions, x); }

    //! \brief Return a const reference to the internal state vector.
    const std::vector<double>& getChainState() const{ return state; }

    //! \brief Return the value of the probability_distribution of the \b i-th chain.

    //! Used by the \b DREAM sampler and probably should not be called by the user.
    double getPDFvalue(size_t i) const{ return pdf_values[i]; }

    //! \brief Allocate (expand) internal storage for the history snapshots, avoids reallocating data when saving a snapshot.

    //! Used by the \b DREAM sampler and probably should not be called by the user.
    void expandHistory(int num_snapshots);

    //! \brief Appends the current state to the history.

    //! Used by the \b DREAM sampler and probably should not be called by the user; \b num_accepted is the number of new chains in the state.
    void saveStateHistory(size_t num_accepted);

    //! \brief Return a const reference to the internal state vector.
    const std::vector<double>& getHistory() const{ return history; }

    //! \brief Return a const reference to the internal state vector.
    const std::vector<double>& getHistoryPDF() const{ return pdf_history; }

    //! \brief Compute the means and variance of the saved history.
    void getHistoryMeanVariance(std::vector<double> &mean, std::vector<double> &var) const;

    //! \brief Return the sample with highest probability, searchers within the history.
    void getApproximateMode(std::vector<double> &mode) const;

    //! \brief Overload that returns the vector.
    std::vector<double> getApproximateMode() const{
        std::vector<double> mode;
        getApproximateMode(mode);
        return mode;
    }

    //! \brief Clear the stored history (does not touch the state).
    void clearHistory();

    //! \brief Returns the acceptance rate of the current history.
    double getAcceptanceRate() const{ return ((pdf_history.empty()) ? 0 : ((double) accepted) / ((double) pdf_history.size())); }

    // file I/O
private:
    size_t num_chains, num_dimensions;
    bool init_state, init_values;

    size_t accepted;

    std::vector<double> state, history;
    std::vector<double> pdf_values, pdf_history;
};

}

#endif
