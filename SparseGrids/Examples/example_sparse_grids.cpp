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

/*!
 * \internal
 * \file example_sparse_grids.cpp
 * \brief Examples for the Tasmanian Sparse Grid module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianSGExamples
 *
 * The main file that runs all sparse grid examples.
 */

/*!
 * \defgroup TasmanianSGExamples Examples for the Tasmanian Sparse Grid module
 *
 * Several examples are included that demonstrate the proper usage of
 * the TasmanianSparseGrid class.
 * The examples include quadrature and interpolation tasks,
 * as well as constructing surrogates from random data.
 */

#ifndef __TASMANIAN_DOXYGEN_SKIP
void sparse_grids_example_01();
void sparse_grids_example_02();
void sparse_grids_example_03();
void sparse_grids_example_04();
void sparse_grids_example_05();
void sparse_grids_example_06();
void sparse_grids_example_07();
void sparse_grids_example_08();
void sparse_grids_example_09();
void sparse_grids_example_10();

int main(int argc, const char**){

    // some of the examples take long time, for post-install testing purposes
    // optionally limit examples to the first 3
    bool limit_examples = (argc > 1);

    sparse_grids_example_01();
    sparse_grids_example_02();
    sparse_grids_example_03();
    sparse_grids_example_04();

    if (limit_examples) return 0; // skip the more expensive examples

    sparse_grids_example_05();
    sparse_grids_example_06();
    sparse_grids_example_07();
    sparse_grids_example_08();
    sparse_grids_example_09();
    sparse_grids_example_10();

    return 0;
}

#endif // __TASMANIAN_DOXYGEN_SKIP
