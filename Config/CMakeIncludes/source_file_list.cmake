########################################################################
# list all source files
########################################################################

set(Tasmanian_source_tasgrid
                SparseGrids/tasgrid_main.cpp
                SparseGrids/TasmanianSparseGrid.hpp
                SparseGrids/tasgridExternalTests.hpp
                SparseGrids/tasgridExternalTests.cpp
                SparseGrids/tasgridTestFunctions.hpp
                SparseGrids/tasgridTestFunctions.cpp
                SparseGrids/tasgridWrapper.hpp
                SparseGrids/tasgridWrapper.cpp)



set(Tasmanian_source_libsparsegrid_cuda
               SparseGrids/tsgCudaKernels.cu
               SparseGrids/tsgCudaBasisEvaluations.hpp
               SparseGrids/tsgCudaLinearAlgebra.hpp)

set(Tasmanian_source_libdream
               DREAM/TasmanianDREAM.hpp
               DREAM/TasmanianDREAM.cpp
               DREAM/tdrEnumerates.hpp
               DREAM/tdrCorePDF.hpp
               DREAM/tdrCorePDF.cpp)

set(Tasmanian_source_tasdream
               DREAM/tasdream_main.cpp
               DREAM/tasdreamExternalTests.hpp
               DREAM/tasdreamExternalTests.cpp
               DREAM/tasdreamTestPDFs.hpp
               DREAM/tasdreamTestPDFs.cpp
               DREAM/tasdreamBenchmark.hpp
               DREAM/tasdreamBenchmark.cpp)

set(Tasmanian_source_matlabsg
               tsgCancelRefine.m
               tsgCleanTempFiles.m
               tsgCopyGrid.m
               tsgDeleteGridByName.m
               tsgDeleteGrid.m
               tsgEstimateAnisotropicCoefficients.m
               tsgEvaluateHierarchy.m
               tsgEvaluate.m
               tsgExample.m
               tsgGetHCoefficients.m
               tsgGetInterpolationWeights.m
               tsgGetNeededIndexes.m
               tsgGetNeededPoints.m
               tsgGetPointsIndexes.m
               tsgGetPoints.m
               tsgGetPolynomialSpaces.m
               tsgGetQuadrature.m
               tsgIntegrate.m
               tsgListGridsByName.m
               tsgLoadHCoefficients.m
               tsgLoadValues.m
               tsgMakeFilenames.m
               tsgMakeGlobal.m
               tsgMakeLocalPolynomial.m
               tsgMakeQuadrature.m
               tsgMakeSequence.m
               tsgMakeWavelet.m
               tsgMergeRefine.m
               tsgPlotPoints2D.m
               tsgReadMatrix.m
               tsgRefineAnisotropic.m
               tsgRefineSurplus.m
               tsgReloadGrid.m
               tsgSummary.m
               tsgCoreTests.m
               tsgWriteCustomRuleFile.m
               tsgWriteMatrix.m)

set(Tasmanian_source_libsparsegrid_fortran
               InterfaceFortran/TasmanianSG.f90
               InterfaceFortran/tsgC2FortranBridge.f90
               InterfaceFortran/tsgC2Fortran.cpp)

set(Tasmanian_source_fortester
               Testing/fortester.f90)
