########################################################################
# MATLAB interface, examples and tests
########################################################################

# copy the source files to the build folder to be used pre-install (stage 1)
file(COPY  tsgCancelRefine.m
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
           tsgMakeFourier.m
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
           tsgWriteMatrix.m
           DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/matlab/")

# set MATLAB (tasgrid and temp folder) for the build folder (stage 1)
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tempMATLAB/") # put temp files if running matlab from the build location
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tsgGetPaths.build.m" "${CMAKE_CURRENT_BINARY_DIR}/matlab/tsgGetPaths.m" @ONLY)

# set MATLAB (tasgrid and temp folder) for the install folder (stage 2)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tsgGetPaths.install.m" "${CMAKE_CURRENT_BINARY_DIR}/configured/tsgGetPaths.m" @ONLY)

# create the MATLAB work-folder on install
install(CODE "file(MAKE_DIRECTORY \"${Tasmanian_MATLAB_WORK_FOLDER}\")")

# installed files
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
        DESTINATION "share/Tasmanian/matlab/"
        FILES_MATCHING PATTERN "*.m"
        PATTERN "tsgGetPaths.*" EXCLUDE)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/configured/tsgGetPaths.m"
        DESTINATION "share/Tasmanian/matlab/"
        PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ)
