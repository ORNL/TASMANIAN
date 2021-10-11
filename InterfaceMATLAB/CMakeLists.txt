########################################################################
# MATLAB interface, examples and tests
########################################################################

set(Tasmanian_matlab_source_files tsgCancelRefine.m
                                  tsgCleanTempFiles.m
                                  tsgCopyGrid.m
                                  tsgDeleteGridByName.m
                                  tsgDeleteGrid.m
                                  tsgEstimateAnisotropicCoefficients.m
                                  tsgEvaluateHierarchy.m
                                  tsgEvaluate.m
                                  tsgExample.m
                                  tsgGetHCoefficients.m
                                  tsgGetHSupport.m
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
                                  tsgMakeExoticQuadrature.m
                                  tsgMakeFilenames.m
                                  tsgMakeFourier.m
                                  tsgMakeGlobal.m
                                  tsgMakeLocalPolynomial.m
                                  tsgMakeQuadrature.m
                                  tsgMakeSequence.m
                                  tsgMakeWavelet.m
                                  tsgMergeRefine.m
                                  tsgPlotPoints2D.m
                                  tsgReadCustomRuleFile.m
                                  tsgReadMatrix.m
                                  tsgRefineAnisotropic.m
                                  tsgRefineSurplus.m
                                  tsgReloadGrid.m
                                  tsgSummary.m
                                  tsgCoreTests.m
                                  tsgWriteCustomRuleFile.m
                                  tsgWriteMatrix.m)

# configure_file() command causes update to the entire CMake environment, if the file changes
# using copy will copy the files on change but will not reconfigure the environment
foreach(Tasmanian_matlab_source_file ${Tasmanian_matlab_source_files})
    add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/matlab/${Tasmanian_matlab_source_file}"
                       COMMAND "${CMAKE_COMMAND}"
                       ARGS -E copy "${CMAKE_CURRENT_SOURCE_DIR}/${Tasmanian_matlab_source_file}" "${CMAKE_CURRENT_BINARY_DIR}/matlab/${Tasmanian_matlab_source_file}"
                       DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${Tasmanian_matlab_source_file}"
                       COMMENT "Copying ${CMAKE_CURRENT_SOURCE_DIR}/${Tasmanian_matlab_source_file}")
    list(APPEND Tasmanian_matlab_pre_install_files "${CMAKE_CURRENT_BINARY_DIR}/matlab/${Tasmanian_matlab_source_file}")
endforeach()
add_custom_target(Tasmanian_matlab_interface ALL DEPENDS "${Tasmanian_matlab_pre_install_files}")

# set MATLAB (tasgrid and temp folder) for the build folder (stage 1)
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/tempMATLAB/") # place for temp files if running matlab from the build location
set(Tasmanian_matlab_tasgrid "${CMAKE_CURRENT_BINARY_DIR}/../SparseGrids/tasgrid")
string(REPLACE " " "\\ " Tasmanian_matlab_tasgrid "${Tasmanian_matlab_tasgrid}")
set(Tasmanian_matlab_workdir "${CMAKE_CURRENT_BINARY_DIR}/tempMATLAB/")
string(REPLACE " " "\\ " Tasmanian_matlab_workdir "${Tasmanian_matlab_workdir}")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/tsgGetPaths.build.m" "${CMAKE_CURRENT_BINARY_DIR}/matlab/tsgGetPaths.m" @ONLY)

# set MATLAB (tasgrid and temp folder) for the install folder (stage 2)
set(Tasmanian_matlab_tasgrid "${Tasmanian_final_install_path}/bin/tasgrid")
string(REPLACE " " "\\ " Tasmanian_matlab_tasgrid "${Tasmanian_matlab_tasgrid}")
message("${Tasmanian_matlab_tasgrid}")
set(Tasmanian_matlab_workdir "${Tasmanian_MATLAB_WORK_FOLDER}/")
string(REPLACE " " "\\ " Tasmanian_matlab_workdir "${Tasmanian_matlab_workdir}")
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
