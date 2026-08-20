cmake_minimum_required(VERSION 3.20)
project(SAMsrcV5 VERSION 5.1.0 LANGUAGES C)

include(CTest)
include(GNUInstallDirs)

if(WIN32)
  add_compile_options(-include "${PROJECT_SOURCE_DIR}/include/sam_platform.h")
endif()

option(SAM_ENABLE_OPENMP "Build parallel SAM programs when OpenMP is available" ON)
option(SAM_BUILD_TESTING "Build native SAM test helpers" ${BUILD_TESTING})
set(SAM_DEPS_ROOT "$ENV{SAM_DEPS_ROOT}" CACHE PATH "Prefix containing FFTW and GSL")

if(SAM_DEPS_ROOT)
  list(PREPEND CMAKE_PREFIX_PATH "${SAM_DEPS_ROOT}")
endif()

find_path(GSL_INCLUDE_DIR gsl/gsl_math.h HINTS "${SAM_DEPS_ROOT}/include")
find_library(GSL_LIBRARY NAMES gsl libgsl HINTS "${SAM_DEPS_ROOT}/lib" "${SAM_DEPS_ROOT}/lib64")
find_library(GSL_CBLAS_LIBRARY NAMES gslcblas libgslcblas HINTS "${SAM_DEPS_ROOT}/lib" "${SAM_DEPS_ROOT}/lib64")
if(NOT GSL_INCLUDE_DIR OR NOT GSL_LIBRARY OR NOT GSL_CBLAS_LIBRARY)
  message(FATAL_ERROR "GSL was not found; set SAM_DEPS_ROOT or CMAKE_PREFIX_PATH")
endif()
add_library(GSL::gsl UNKNOWN IMPORTED)
set_target_properties(GSL::gsl PROPERTIES
  IMPORTED_LOCATION "${GSL_LIBRARY}"
  INTERFACE_INCLUDE_DIRECTORIES "${GSL_INCLUDE_DIR}"
)
add_library(GSL::gslcblas UNKNOWN IMPORTED)
set_target_properties(GSL::gslcblas PROPERTIES
  IMPORTED_LOCATION "${GSL_CBLAS_LIBRARY}"
  INTERFACE_INCLUDE_DIRECTORIES "${GSL_INCLUDE_DIR}"
)
find_path(FFTW3_INCLUDE_DIR fftw3.h HINTS "${SAM_DEPS_ROOT}/include")
find_library(FFTW3_LIBRARY NAMES fftw3 libfftw3 HINTS "${SAM_DEPS_ROOT}/lib" "${SAM_DEPS_ROOT}/lib64")
if(NOT FFTW3_INCLUDE_DIR OR NOT FFTW3_LIBRARY)
  message(FATAL_ERROR "FFTW3 was not found; set SAM_DEPS_ROOT or CMAKE_PREFIX_PATH")
endif()
add_library(FFTW3::fftw3 UNKNOWN IMPORTED)
set_target_properties(FFTW3::fftw3 PROPERTIES
  IMPORTED_LOCATION "${FFTW3_LIBRARY}"
  INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}"
)

set(SAM_LIBRARY_SOURCES
  Libs/LFULIB/Coil2B.c
  Libs/LFULIB/Coil2B_plus.c
  Libs/LFULIB/DIPLeadField.c
  Libs/LFULIB/DIPLeadField2D.c
  Libs/LFULIB/DIPSolve.c
  Libs/LFULIB/DIPSolve2D.c
  Libs/LFULIB/DIPSolveFwd.c
  Libs/LFULIB/ECDLeadField.c
  Libs/LFULIB/ECDSolve.c
  Libs/LFULIB/ECDSolveFwd.c
  Libs/LFULIB/FwdOPM.c
  Libs/LFULIB/FwdSensIntPnt.c
  Libs/LFULIB/FwdSensor.c
  Libs/LFULIB/FwdVec.c
  Libs/LFULIB/GetBasis.c
  Libs/LFULIB/GetCoeffs.c
  Libs/LFULIB/GetHull.c
  Libs/LFULIB/GetPatchVertices.c
  Libs/LFULIB/LUsolve.c
  Libs/LFULIB/LtoB.c
  Libs/LFULIB/LtoB_old.c
  Libs/LFULIB/MtoB.c
  Libs/LFULIB/MtoB_plus.c
  Libs/LFULIB/NLeadField.c
  Libs/LFULIB/NSolveFwd.c
  Libs/LFULIB/NSolveMoments.c
  Libs/LFULIB/NSolveWts.c
  Libs/LFULIB/ROISolveWts.c
  Libs/LFULIB/RandomMesh.c
  Libs/LFULIB/SolveMoment.c
  Libs/LFULIB/geig22.c
  Libs/LFULIB/rnynm.c
  Libs/LFULIB/xbd.c
  Libs/SAMLIB/AdaptiveMeshScatter.c
  Libs/SAMLIB/Atlas2Mesh.c
  Libs/SAMLIB/Atlas2RandomMesh.c
  Libs/SAMLIB/CtoS.c
  Libs/SAMLIB/Dflip.c
  Libs/SAMLIB/DipoleNM.c
  Libs/SAMLIB/EtoR.c
  Libs/SAMLIB/FIELD2Vector.c
  Libs/SAMLIB/FitDip.c
  Libs/SAMLIB/FitOPM.c
  Libs/SAMLIB/FitOPM_plus.c
  Libs/SAMLIB/GetAtlas.c
  Libs/SAMLIB/GetCov.c
  Libs/SAMLIB/GetDsData.c
  Libs/SAMLIB/GetDsInfo.c
  Libs/SAMLIB/GetFilePath.c
  Libs/SAMLIB/GetGIFTIWts.c
  Libs/SAMLIB/GetGiiAtlas.c
  Libs/SAMLIB/GetHDM.c
  Libs/SAMLIB/GetMarkers.c
  Libs/SAMLIB/GetNIFTI.c
  Libs/SAMLIB/GetNIFTIMask.c
  Libs/SAMLIB/GetNIFTIWts.c
  Libs/SAMLIB/GetNoise.c
  Libs/SAMLIB/GetROIVertices.c
  Libs/SAMLIB/GetSurfer.c
  Libs/SAMLIB/GetTargets.c
  Libs/SAMLIB/GetVoxels.c
  Libs/SAMLIB/GetVoxelsFS.c
  Libs/SAMLIB/MagFitDipole.c
  Libs/SAMLIB/MagFitSensor.c
  Libs/SAMLIB/MeshScatter.c
  Libs/SAMLIB/MeshScatterMS.c
  Libs/SAMLIB/MoveField.c
  Libs/SAMLIB/MoveFrame.c
  Libs/SAMLIB/MoveMEG.c
  Libs/SAMLIB/MoveSensor.c
  Libs/SAMLIB/OPMFIELD2Vector.c
  Libs/SAMLIB/OrthoPlane.c
  Libs/SAMLIB/P5toP6.c
  Libs/SAMLIB/ParseFileName.c
  Libs/SAMLIB/ParsePath.c
  Libs/SAMLIB/PutCov.c
  Libs/SAMLIB/PutGIFTISurf.c
  Libs/SAMLIB/PutGIFTIWts.c
  Libs/SAMLIB/PutNIFTIImage.c
  Libs/SAMLIB/PutNIFTITargets.c
  Libs/SAMLIB/PutNIFTIWts.c
  Libs/SAMLIB/PutWts.c
  Libs/SAMLIB/RandomGauss.c
  Libs/SAMLIB/SAMsolve.c
  Libs/SAMLIB/SensorNM.c
  Libs/SAMLIB/SetHullOrigin.c
  Libs/SAMLIB/SetSphereOrigins.c
  Libs/SAMLIB/StoC.c
  Libs/SAMLIB/Vector2FIELD.c
  Libs/SAMLIB/Vector2OPMFIELD.c
  Libs/SAMLIB/Vector_FIELD_convert.c
  Libs/SAMLIB/gifti_util.c
  Libs/SAMLIB/gsl_diag.c
  Libs/SAMLIB/gsl_vector_norm.c
  Libs/SAMLIB/read_segfile.c
  Libs/SIGLIB/Butterworth.c
  Libs/SIGLIB/FFTfilter.c
  Libs/SIGLIB/Hanning.c
  Libs/SIGLIB/Levy.c
  Libs/SIGLIB/bdiir.c
  Libs/SIGLIB/debase.c
  Libs/SIGLIB/demean.c
  Libs/SIGLIB/detrend.c
  Libs/SIGLIB/hilbert.c
  Libs/SIGLIB/kendall.c
  Libs/SIGLIB/mkfilter.c
  Libs/SIGLIB/mkiir.c
  Libs/SIGLIB/mknotch.c
  Libs/SIGLIB/mode.c
  Libs/SIGLIB/pinv.c
  Libs/SIGLIB/power.c
  Libs/SIGLIB/response.c
  Libs/SIGLIB/spectrum.c
  Libs/SIGLIB/whiten.c
  Libs/samutil/cleanup.c
  Libs/samutil/file.c
  Libs/samutil/mem.c
  Libs/samutil/platform.c
  Libs/samutil/sam_param.c
  Libs/samutil/sam_parse.c
  Libs/samutil/shuffle.c
)

add_library(samcore STATIC ${SAM_LIBRARY_SOURCES})
target_compile_features(samcore PUBLIC c_std_17)
target_compile_definitions(samcore PUBLIC HAVE_INLINE=1 GSL_RANGE_CHECK_OFF=1)
target_include_directories(samcore PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(samcore PUBLIC GSL::gsl GSL::gslcblas FFTW3::fftw3)
if(NOT WIN32)
  target_link_libraries(samcore PUBLIC m)
endif()

set(SAM_OPENMP_TARGETS sam_wts sam_3d sam_4d sam_power)
if(SAM_ENABLE_OPENMP)
  find_package(OpenMP QUIET COMPONENTS C)
endif()

function(sam_add_program target source output_name)
  add_executable(${target} ${source})
  target_compile_features(${target} PRIVATE c_std_17)
  target_include_directories(${target} PRIVATE "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(${target} PRIVATE samcore)
  set_target_properties(${target} PROPERTIES OUTPUT_NAME "${output_name}")
  if(OpenMP_C_FOUND AND target IN_LIST SAM_OPENMP_TARGETS)
    target_link_libraries(${target} PRIVATE OpenMP::OpenMP_C)
  endif()
  if(MINGW)
    target_link_options(${target} PRIVATE -static-libgcc)
  endif()
  install(TARGETS ${target} RUNTIME DESTINATION samsrcv5/_bin COMPONENT python)
endfunction()

sam_add_program(sam_cov Mains/SAMcov/SAMcov.c sam_cov)
sam_add_program(sam_wts Mains/SAMwts/SAMwts.c sam_wts)
sam_add_program(sam_3d Mains/SAM3d/SAM3d.c sam_3d)
sam_add_program(sam_4d Mains/SAM4d/SAM4d.c sam_4d)
sam_add_program(sam_ers Mains/SAMers/SAMers.c sam_ers)
sam_add_program(sam_epi Mains/SAMepi/SAMepi.c sam_epi)
sam_add_program(sam_power Mains/SAMpower/SAMpower.c sam_power)
sam_add_program(sam_simulate Mains/SAMsimulate/SAMsimulate.c sam_simulate)
sam_add_program(opm_sim Mains/OPMsim/OPMsim.c OPMsim)

add_executable(stats_1d Mains/3dNormalize/1dstats.c)
set_target_properties(stats_1d PROPERTIES OUTPUT_NAME 1dstats)
target_compile_features(stats_1d PRIVATE c_std_17)
if(NOT WIN32)
  target_link_libraries(stats_1d PRIVATE m)
endif()
install(TARGETS stats_1d RUNTIME DESTINATION samsrcv5/_bin COMPONENT python)

add_executable(meshnorm
  Mains/orthohull/meshnorm.c
  Mains/orthohull/err.c
  Mains/orthohull/mem.c
)
target_compile_features(meshnorm PRIVATE c_std_17)
target_include_directories(meshnorm PRIVATE Mains/orthohull)
if(NOT WIN32)
  target_link_libraries(meshnorm PRIVATE m)
endif()
install(TARGETS meshnorm RUNTIME DESTINATION samsrcv5/_bin COMPONENT python)

set(SAM_LEGACY_PYTHON
  Mains/FSnormals/FSnormals.py
  Mains/FSnormals/mkGiiAtlas.py
  Mains/orthohull/fid.py
  Mains/orthohull/fiddist.py
  Mains/orthohull/orthohull.py
  Mains/orthohull/plothull.py
  Mains/orthohull/runscript.py
  Mains/orthohull/samutil.py
  Mains/orthohull/thd_atr.py
)
install(FILES ${SAM_LEGACY_PYTHON} DESTINATION samsrcv5/_legacy COMPONENT python)
install(FILES Mains/orthohull/ply2fid
  DESTINATION samsrcv5/_legacy RENAME ply2fid.py COMPONENT python)
install(FILES Mains/ROIbuilder/ROIbuilder
  DESTINATION samsrcv5/_legacy RENAME ROIbuilder.py COMPONENT python)
install(FILES
  Mains/orthohull/master+orig.BRIK.gz
  Mains/orthohull/master+orig.HEAD
  Mains/orthohull/null.tag
  Mains/ROIbuilder/ROIbuilder.ui
  Mains/ROIbuilder/SessionAtlases.niml
  DESTINATION samsrcv5/data COMPONENT python
)
install(FILES LICENSE LICENSE.package THIRD_PARTY_NOTICES.md
  DESTINATION samsrcv5/licenses COMPONENT python)
install(FILES
  "${SAM_DEPS_ROOT}/share/samsrcv5-licenses/GSL-GPL-3.0.txt"
  "${SAM_DEPS_ROOT}/share/samsrcv5-licenses/FFTW-GPL-2.0-or-later.txt"
  DESTINATION samsrcv5/licenses COMPONENT python OPTIONAL
)

if(SAM_BUILD_TESTING)
  add_executable(sam_test_core test/unit/test_core.c)
  target_include_directories(sam_test_core PRIVATE "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(sam_test_core PRIVATE samcore)
  add_test(NAME sam_test_core COMMAND sam_test_core)

  add_executable(sam_inspect_ctf test/integration/inspect_ctf.c)
  target_include_directories(sam_inspect_ctf PRIVATE "${PROJECT_SOURCE_DIR}/include")
  target_link_libraries(sam_inspect_ctf PRIVATE samcore)
  install(TARGETS sam_inspect_ctf RUNTIME DESTINATION samsrcv5/_bin COMPONENT python)
endif()
