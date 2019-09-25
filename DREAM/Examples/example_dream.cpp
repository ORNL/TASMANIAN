#include "Tasmanian.hpp"

/*!
 * \internal
 * \file example_dream.cpp
 * \brief Examples for the Tasmanian DREAM module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianDREAMExamples
 *
 * The main header required to gain access to the DREAM capabilities of Tasmanian.
 * The header will include all files needed by the DREAM module including
 * the TasmanianSparseGrid.hpp header.
 * \endinternal
 */

//! \defgroup TasmanianDREAMExamples Examples for the Tasmanian DREAM module
//!
//! Several examples are included that demonstrate the proper usage of the DREAM interface.
//! The examples include random sampling from an arbitrary probability distribution,
//! Bayesian inference using a custom model or sparse grids interpolant,
//! as well as simple optimization problem (i.e., search for the mode of a probability distribution).

using namespace std;

#ifndef __TASMANIAN_DOXYGEN_SKIP
void dream_example_01();
void dream_example_02();
void dream_example_03();
void dream_example_04();
void dream_example_05();

int main(int argc, const char**){
/*
 * The purpose of this file is to demonstrate the proper way to call
 * functions from the TASMANIAN DREAM Module.
 */
    dream_example_01();
    dream_example_02();

    if (argc > 1) return 0; // fast testing used to check if the library linked correctly

    dream_example_03();
    dream_example_04();
    dream_example_05();

    cout << "\n" << "---------------------------------------------------------------------------------------------------\n\n";

    return 0;
}

#endif
