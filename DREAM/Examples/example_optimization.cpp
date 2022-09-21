#include "Tasmanian.hpp"

/*!
 * \internal
 * \file example_optimization.cpp
 * \brief Examples for the Tasmanian Optimization module.
 * \author Miroslav Stoyanov
 * \ingroup TasmanianOPTExamples
 *
 * The main file for all optimization examples.
 * \endinternal
 */

/*!
 * \defgroup TasmanianOPTExamples Examples for the Tasmanian Optimization module
 *
 * Several examples are included that demonstrate the proper usage of the Optimization interface.
 */

using namespace std;

#ifndef __TASMANIAN_DOXYGEN_SKIP
void optimizaiton_example_01();

int main(int argc, const char**){
/*
 * The purpose of this file is to demonstrate the proper way to call
 * functions from the Tasmanian Optimization Module.
 */
    optimizaiton_example_01();

    if (argc > 1) return 0; // fast testing used to check if the library linked correctly

    cout << "\n" << "---------------------------------------------------------------------------------------------------\n\n";

    return 0;
}

#endif
