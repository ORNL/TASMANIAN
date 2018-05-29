cmake_minimum_required(VERSION 3.5)

########################################################################
# MATLAB interface, examples and tests
########################################################################

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

# set MATLAB for the build folder (stage 1)
set(Tasmanian_tasgrid_path "${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/tasgrid")
set(Tasmanian_string_matlab_work ${CMAKE_CURRENT_BINARY_DIR}/tempMATLAB/)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/tempMATLAB/)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tsgGetPaths.in.m" "${CMAKE_CURRENT_BINARY_DIR}/matlab/tsgGetPaths.m")

foreach(Tasmanian_matlab_file ${Tasmanian_source_matlabsg})
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${Tasmanian_matlab_file}" "${CMAKE_CURRENT_BINARY_DIR}/matlab/${Tasmanian_matlab_file}" COPYONLY)
endforeach()

set(Tasmanian_tasgrid_path "${CMAKE_INSTALL_PREFIX}/bin/tasgrid")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tsgGetPaths.in.m" "${CMAKE_CURRENT_BINARY_DIR}/configured_for_install/tsgGetPaths.m")

# installed files
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
        DESTINATION "matlab"
        FILES_MATCHING PATTERN "*.m"
        PATTERN "tsgGetPaths.*" EXCLUDE)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured_for_install/tsgGetPaths.m"
        DESTINATION "matlab"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
