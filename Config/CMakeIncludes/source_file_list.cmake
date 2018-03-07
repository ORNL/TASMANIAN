########################################################################
# list all source files
########################################################################

set(tasgrid_src SparseGrids/tasgrid_main.cpp
                SparseGrids/TasmanianSparseGrid.hpp
                SparseGrids/tasgridExternalTests.hpp
                SparseGrids/tasgridExternalTests.cpp
                SparseGrids/tasgridTestFunctions.hpp
                SparseGrids/tasgridTestFunctions.cpp
                SparseGrids/tasgridWrapper.hpp
                SparseGrids/tasgridWrapper.cpp)

set(libtsg_src SparseGrids/TasmanianSparseGrid.hpp
               SparseGrids/TasmanianSparseGrid.cpp
               SparseGrids/tsgAcceleratedDataStructures.hpp
               SparseGrids/tsgAcceleratedDataStructures.cpp
               SparseGrids/tsgCacheLagrange.hpp
               SparseGrids/tsgCacheLagrange.cpp
               SparseGrids/tsgCoreOneDimensional.hpp
               SparseGrids/tsgCoreOneDimensional.cpp
               SparseGrids/tsgEnumerates.hpp
               SparseGrids/tsgGridCore.hpp
               SparseGrids/tsgGridCore.cpp
               SparseGrids/tsgGridGlobal.hpp
               SparseGrids/tsgGridGlobal.cpp
               SparseGrids/tsgGridWavelet.hpp
               SparseGrids/tsgGridWavelet.cpp
               SparseGrids/tsgHardCodedTabulatedRules.hpp
               SparseGrids/tsgHardCodedTabulatedRules.cpp
               SparseGrids/tsgHiddenExternals.hpp
               SparseGrids/tsgGridLocalPolynomial.hpp
               SparseGrids/tsgGridLocalPolynomial.cpp
               SparseGrids/tsgGridSequence.hpp
               SparseGrids/tsgGridSequence.cpp
               SparseGrids/tsgIndexManipulator.hpp
               SparseGrids/tsgIndexManipulator.cpp
               SparseGrids/tsgIndexSets.hpp
               SparseGrids/tsgIndexSets.cpp
               SparseGrids/tsgLinearSolvers.hpp
               SparseGrids/tsgLinearSolvers.cpp
               SparseGrids/tsgOneDimensionalWrapper.hpp
               SparseGrids/tsgOneDimensionalWrapper.cpp
               SparseGrids/tsgRuleLocalPolynomial.hpp
               SparseGrids/tsgRuleLocalPolynomial.cpp
               SparseGrids/tsgRuleWavelet.hpp
               SparseGrids/tsgRuleWavelet.cpp
               SparseGrids/tsgSequenceOptimizer.hpp
               SparseGrids/tsgSequenceOptimizer.cpp)

set(libtdr_src DREAM/TasmanianDREAM.hpp
               DREAM/TasmanianDREAM.cpp
               DREAM/tdrEnumerates.hpp
               DREAM/tdrCorePDF.hpp
               DREAM/tdrCorePDF.cpp
               SparseGrids/TasmanianSparseGrid.hpp)

set(tasdre_src DREAM/tasdream_main.cpp
               DREAM/tasdreamExternalTests.hpp
               DREAM/tasdreamExternalTests.cpp
               DREAM/tasdreamTestPDFs.hpp
               DREAM/tasdreamTestPDFs.cpp
               DREAM/tasdreamBenchmark.hpp
               DREAM/tasdreamBenchmark.cpp)

set(tasmat_src tsgCancelRefine.m
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

set(libfor_src InterfaceFortran/TasmanianSG.f90
               InterfaceFortran/tsgC2FortranBridge.f90
               InterfaceFortran/tsgC2Fortran.cpp)

set(tasfor_src Testing/fortester.f90)
