# =========================================================================
# HOPR CMake configuration
# =========================================================================
# NOTE: use minimum CMake version required by tools/libs (Paraview, HDF5, FFTW)
CMAKE_MINIMUM_REQUIRED(VERSION 3.5.2)
# Policies (see https://cmake.org/cmake/help/v3.13/manual/cmake-policies.7.html)
# VERSION_GREATER_EQUAL is incompatible with CMAKE versions < 3.7
IF(NOT(${CMAKE_VERSION} VERSION_LESS "3.6.0"))
  SET(GITSHALLOW "GIT_SHALLOW ON")
ENDIF()
IF(NOT(${CMAKE_VERSION} VERSION_LESS "3.12.0"))
  CMAKE_POLICY(SET CMP0074 NEW)
ENDIF()
# Folder for custom CMake source files
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
INCLUDE(ExternalProject)
INCLUDE(CMakeDependentOption)

# =========================================================================
# Check whether we are on Github or at IAG
# =========================================================================
EXECUTE_PROCESS(COMMAND git ls-remote --get-url OUTPUT_VARIABLE GIT_ORIGIN)
STRING(FIND "${GIT_ORIGIN}" ".iag.uni-stuttgart.de" HOMERUN)
STRING(FIND "${GIT_ORIGIN}" "git@" USESSH)
IF (HOMERUN GREATER -1)
  IF(USESSH EQUAL 0)
    SET(HOPR_DLPATH "git@gitlab.iag.uni-stuttgart.de:")
  ELSE()
    SET(HOPR_DLPATH "https://gitlabext.iag.uni-stuttgart.de/")
  ENDIF()
  SET(HDF5DOWNLOAD ${HOPR_DLPATH}libs/hdf5.git )
  SET(CGNSDOWNLOAD ${HOPR_DLPATH}libs/CGNS.git )
ELSE()
  SET(CGNSDOWNLOAD https://github.com/CGNS/CGNS.git )
  SET(HDF5DOWNLOAD https://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git )
ENDIF()
SET (HDF5TAG "hdf5-1_10_4" CACHE STRING "HDF5 version tag")
SET (CGNSTAG "v3.3.1"      CACHE STRING "CGNS version tag")

# =========================================================================
# Detect machine environements
# =========================================================================
SITE_NAME(CMAKE_HOSTNAME)
STRING(FIND "${CMAKE_HOSTNAME}" "eslogin" ONHORNET)
# Set options that need to be enabled BEFORE ANY PROJECT OR LANGUAGE COMMAND
IF (ONHORNET GREATER -1)
  SET(CMAKE_C_COMPILER       cc)
  SET(CMAKE_CXX_COMPILER     CC)
  SET(CMAKE_Fortran_COMPILER ftn)
  SET(HDF5_USE_STATIC_LIBRARIES TRUE)
  MARK_AS_ADVANCED(C_PATH CXX_PATH Fortran_PATH)
  MESSAGE(STATUS "Compiling on Hazelhen")
  ADD_DEFINITIONS(-DLUSTRE)
ENDIF()

# =========================================================================
# Project definitions and output paths
# =========================================================================
PROJECT(Hopr)
ENABLE_LANGUAGE(C Fortran)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# Set custom install dir (needs to be done after project has been defined!)
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Prefix prepended to install directories" FORCE)
ENDIF()

# =========================================================================
# Options
# =========================================================================

MACRO(SUBDIRLIST result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      LIST(APPEND dirlist ${child})
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

OPTION(HOPR_MPI "Compile SINGLE or MPI version" ON)
OPTION(HOPR_BUILD_YPLUSESTIMATOR "Compile and build Y+-Estimator" OFF)
OPTION(HOPR_BUILD_HDF5 "Compile and build hdf5 library" ON)
OPTION(HOPR_BUILD_CGNS "Compile and build CGNS library" ON)
SET(HOPR_CGNS_INT "32" CACHE STRING "integer type in CGNS lib")

# =========================================================================
# USERBLOCK + Preproc_flags
# =========================================================================

# A function to get all user defined variables with a specified prefix
function (getListOfVarsStartingWith _prefix _varResult)
    GET_CMAKE_PROPERTY(_vars CACHE_VARIABLES)
    STRING(REGEX MATCHALL "(^|;)${_prefix}[A-Za-z0-9_]*" _matchedVars "${_vars}")
    SET(${_varResult} ${_matchedVars} PARENT_SCOPE)
endfunction()

# A function to get all user defined variables with a specified prefix
SET(configuration ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/configuration.cmake)
FILE(WRITE ${configuration} "")
GETLISTOFVARSSTARTINGWITH("HOPR_" hoprVars)
FOREACH (_var IN LISTS hoprVars)
  GET_PROPERTY(currentHelpString CACHE "${_var}" PROPERTY HELPSTRING)
  SET(boolian (${${_var}} STREQUAL "ON" OR ${${_var}} STREQUAL "OFF")) 
  IF(${boolian})
    FILE(APPEND ${configuration} "OPTION(${_var} \"${currentHelpString}\" ${${_var}})\n")
  ELSE()
    IF(${_var})
      STRING(REPLACE "\\" "\\\\" ${_var} ${${_var}})
    ENDIF()
    FILE(APPEND ${configuration} "SET(${_var} \"${${_var}}\" CACHE STRING \"${currentHelpString}\")\n")
  ENDIF()
ENDFOREACH()

# =========================================================================
# Build type
# =========================================================================
# make sure that the default is a RELEASE
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug Release Profile."
      FORCE)
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release Profile)
ENDIF (NOT CMAKE_BUILD_TYPE)

IF (CMAKE_BUILD_TYPE MATCHES "Debug")
  ADD_DEFINITIONS("-DDEBUG")
ENDIF()

# =========================================================================
# Location of binary and filenames
# =========================================================================

# append relative filename-macro for __FILENAME__ in Stamp of abort function (see hopr.h)
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D__FILENAME__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")

# add basedir for location of corresponding userblock-file
ADD_DEFINITIONS("-DBASEDIR='\"${CMAKE_CURRENT_BINARY_DIR}/\"'")

# =========================================================================
# COMPILER FLAGS
# =========================================================================

# FFLAGS depend on the compiler
GET_FILENAME_COMPONENT (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -fbackslash -ffree-line-length-0 -DGNU")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}     -O3 -march=native  -finline-functions")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -pg -O3 -march=native  -finline-functions")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -ggdb3 -fbounds-check -finit-real=nan -fbacktrace  -Wall")
  # add flags only for compiling not linking!
  SET (HOPR_COMPILE_FLAGS "-xf95-cpp-input -fPIC")
ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -r8 -i4 -traceback -warn all -shared-intel -DINTEL")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}    -O2 -xHost -qopt-report0 -qopt-report-phase=vec")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -p -O2 -xHost -qopt-report0 -qopt-report-phase=vec")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -fpe0 -traceback -check all,noarg_temp_created,noformat,nooutput_conversion,pointer,uninit -init=snan -init=arrays")
  SET (CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS}    -O2 -xHost -ip")
  SET (CMAKE_C_FLAGS_PROFILE   "${CMAKE_C_FLAGS} -p -O2 -xHost -ip")
  SET (CMAKE_C_FLAGS_DEBUG     "${CMAKE_C_FLAGS} -g")
  # add flags only for compiling not linking!
  SET (HOPR_COMPILE_FLAGS "-fpp -allow nofpp_comments -assume bscc")

ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -f free -s real64 -em -DCRAY")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O2 -hfp3 -p . -rm")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -O2 -hfp3 -h profile_generate -p . -rm")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -O0 -eD -rm")
  # add flags only for compiling not linking!
  SET (HOPR_COMPILE_FLAGS "-F")

ELSE()
  MESSAGE(SEND_ERROR "Unknown compiler")
ENDIF()

# =========================================================================
# Add the libraries
# =========================================================================

IF(HOPR_MPI)
  SET(HOPR_EXTERNAL_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${CMAKE_Fortran_COMPILER_ID}-MPI)
ELSE()
  SET(HOPR_EXTERNAL_LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${CMAKE_Fortran_COMPILER_ID})
ENDIF()
MARK_AS_ADVANCED(HOPR_EXTERNAL_LIB_DIR)

# =========================================================================
# HDF5 library
# =========================================================================
SET(LIB_TYPE STATIC) # OR SHARED
STRING(TOLOWER ${LIB_TYPE} SEARCH_TYPE)

# We support two methods for finding HDF5:
# the version built using configure scripts and using CMake
# Support for CMake-built HDF5 is limited to version >1.10.2
# which require at CMake >3.10
# Recommendations: use configure version if possible

IF(NOT HOPR_BUILD_HDF5)
  SET(HOPR_HDF5_CMAKE TRUE) # try to find the CMake version
  FIND_PACKAGE(HDF5 NAMES hdf5 COMPONENTS C Fortran ${SEARCH_TYPE} QUIET PATH_SUFFIXES share/cmake )
  IF (HDF5_FOUND)
    IF(${CMAKE_VERSION} VERSION_LESS "3.10.0")
      MESSAGE (WARNING "HDF5 built with CMake was found! This feature is only supported for CMake 3.10+ and HDF5 1.10.2+ !")
    ENDIF()
    MESSAGE (STATUS "HDF5 C libs:${HDF5_FOUND} static:${HDF5_static_C_FOUND} and shared:${HDF5_shared_C_FOUND}")
    MESSAGE (STATUS "HDF5 Fortran libs: static:${HDF5_static_Fortran_FOUND} and shared:${HDF5_shared_Fortran_FOUND}")
    MESSAGE (STATUS "h5diff found:${HDF5_DIFF_EXECUTABLE}")
  ELSE()
    SET(HOPR_HDF5_CMAKE FALSE) # try to find the configure version
    FIND_PACKAGE(HDF5 COMPONENTS C Fortran)
  ENDIF()
ELSE()
  # define here the path to download from

  SET(HOPR_HDF5_DLDIR ${HOPR_EXTERNAL_LIB_DIR}/HDF5)
  SET(HOPR_HDF5_DIR  ${HOPR_HDF5_DLDIR}/build)
  IF (NOT EXISTS "${HOPR_HDF5_DIR}/lib/libhdf5.a")
    IF(HOPR_MPI)
      SET(HOPR_HDF5PARALLEL --enable-parallel)
      SET(HOPR_HDF5FC ${MPI_Fortran_COMPILER})
      SET(HOPR_HDF5CC ${MPI_C_COMPILER})
    ELSE()
      UNSET(HOPR_HDF5PARALLEL)
      SET(HOPR_HDF5FC ${CMAKE_Fortran_COMPILER})
      SET(HOPR_HDF5CC ${CMAKE_C_COMPILER} )
    ENDIF()
    #NOTE: old 1.8.X versions of HDF5 require additional flag --enable-fortran2003
    ExternalProject_Add(HDF5
      GIT_REPOSITORY ${HDF5DOWNLOAD}
      GIT_TAG ${HDF5TAG}
      ${${GITSHALLOW}}
      PREFIX ${HOPR_HDF5_DIR}
      UPDATE_COMMAND ""
      CONFIGURE_COMMAND FC=${HOPR_HDF5FC} CC=${HOPR_HDF5CC} ${HOPR_HDF5_DIR}/src/HDF5/configure --prefix=${HOPR_HDF5_DIR} --with-pic --enable-fortran ${HOPR_HDF5PARALLEL} --libdir=${HOPR_HDF5_DIR}/lib
      BUILD_COMMAND ${MAKE}
    )
    LIST(APPEND SELFBUILTEXTERNALS HDF5)
    # CMake HDF5 is fast but not that reliable
    #    ExternalProject_Add(HDF5
    #      GIT_REPOSITORY ${HDF5DOWNLOAD}
    #  GIT_TAG ${HDF5TAG}
    #      PREFIX ${HOPR_HDF5_DIR}
    #      UPDATE_COMMAND ""
    #  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${HOPR_HDF5_DIR} -DCMAKE_C_COMPILER=${HOPR_HDF5CC} -DCMAKE_Fortran_COMPILER=${HOPR_HDF5FC} -DBUILD-TESTING=OFF -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_BUILD_TOOLS=OFF -DHDF5_BUILD_FORTRAN=ON -DHDF5_ENABLE_PARALLEL=ON
    #      BUILD_COMMAND ${MAKE}
    #    )
  ENDIF()

  SET(HOPR_HDF5_CMAKE FALSE)
  SET(HDF5_hdf5_LIBRARY_RELEASE         ${HOPR_HDF5_DIR}/lib/libhdf5.a)
  SET(HDF5_hdf5_fortran_LIBRARY_RELEASE ${HOPR_HDF5_DIR}/lib/libhdf5_fortran.a)
  SET(HDF5_C_INCLUDE_DIR                ${HOPR_HDF5_DIR}/include)
  SET(HDF5_Fortran_INCLUDE_DIR          ${HOPR_HDF5_DIR}/include)
  UNSET(HDF5_LIBRARIES)
  UNSET(HDF5_INCLUDE_DIR_FORTRAN)
  UNSET(HDF5_INCLUDE_DIR)
  LIST(APPEND HDF5_INCLUDE_DIR_FORTRAN ${HDF5_Fortran_INCLUDE_DIR} ${HDF5_C_INCLUDE_DIR})
  LIST(APPEND HDF5_INCLUDE_DIR  ${HDF5_C_INCLUDE_DIR})
  FIND_LIBRARY(HDF5_z_LIBRARY_RELEASE z)
  LIST(APPEND HDF5_LIBRARIES ${HDF5_hdf5_fortran_LIBRARY_RELEASE} ${HDF5_hdf5_LIBRARY_RELEASE} ${HDF5_z_LIBRARY_RELEASE} -ldl)
ENDIF()

IF(HOPR_HDF5_CMAKE)
  INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR} ${HDF5_INCLUDE_DIR_FORTRAN})
  LIST(APPEND linkedlibs ${HDF5_C_${LIB_TYPE}_LIBRARY} ${HDF5_FORTRAN_${LIB_TYPE}_LIBRARY} )
ELSE()
  INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR_FORTRAN} ${HDF5_INCLUDE_DIR})
  MESSAGE(STATUS "Found HDF5 Libs: " ${HDF5_LIBRARIES})
  LIST(APPEND linkedlibs ${HDF5_LIBRARIES} )
ENDIF()

# =========================================================================
# LAPACK/MKL library
# =========================================================================
# use Lapack/Blas for GNU
MESSAGE(STATUS "Using BLAS/Lapack library")
IF (HOPR_MKL)
  SET(BLA_VENDOR "Intel") #search only for Intel BLAS (=MKL)
ENDIF()
FIND_PACKAGE(LAPACK QUIET)
IF (LAPACK_FOUND)
  LIST(APPEND linkedlibs ${LAPACK_LIBRARIES})
ENDIF()

# =========================================================================
# MPI
# =========================================================================
IF(HOPR_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  IF (NOT MPI_Fortran_NO_INTERROGATE)
    FOREACH(DIR ${MPI_INCLUDE_PATH})
      INCLUDE_DIRECTORIES(${DIR})
    ENDFOREACH()
    FOREACH(DIR ${MPI_Fortran_INCLUDE_PATH})
      INCLUDE_DIRECTORIES(${DIR})
    ENDFOREACH()
    LIST(APPEND linkedlibs ${MPI_Fortran_LIBRARIES})
  ENDIF()
  MARK_AS_ADVANCED(MPI_LIBRARY MPI_EXTRA_LIBRARY) # don't show these ones
  MESSAGE(STATUS "Building HOPR with MPI")
  MESSAGE(STATUS "MPI Compiler: " ${MPI_Fortran_COMPILER})
  ADD_DEFINITIONS(-DUSE_MPI=1)
ELSE()
  ADD_DEFINITIONS(-DUSE_MPI=0)
ENDIF()

# =========================================================================
# CGNS library
# =========================================================================
ADD_DEFINITIONS(-DPP_CGNS_INT=${HOPR_CGNS_INT})

IF (HOPR_BUILD_CGNS)
  SET(HOPR_CGNS_DLDIR ${HOPR_EXTERNAL_LIB_DIR}/CGNS)
  SET(HOPR_CGNS_DIR   ${HOPR_CGNS_DLDIR}/build)
  IF (NOT EXISTS "${HOPR_CGNS_DIR}/lib/libcgns.a")
    STRING(COMPARE EQUAL ${HOPR_CGNS_INT} "64" HOPR_CGNS_64BIT)

    ExternalProject_Add(cgns
    GIT_REPOSITORY ${CGNSDOWNLOAD}
    GIT_TAG ${CGNSTAG}
    #PATCH_COMMAND sed -i s/^c/!/g ${HOPR_CGNS_DIR}/src/cgns/src/cgnslib_f.h.in
    PREFIX ${HOPR_CGNS_DIR}
    CMAKE_ARGS  -DCMAKE_INSTALL_PREFIX=${HOPR_CGNS_DIR}/build /
    -DCGNS_ENABLE_FORTRAN=ON   -DCGNS_ENABLE_64BIT=${HOPR_CGNS_64BIT} -DCGNS_BUILD_SHARED=OFF -DCGNS_USE_SHARED=OFF /
    -DCMAKE_BUILD_TYPE=Release -DCGNS_BUILD_CGNSTOOLS=OFF -DCGNS_ENABLE_HDF5=OFF -DCGNS_ENABLE_PARALLEL=OFF -DCGNS_ENABLE_TESTS=OFF -DCMAKE_SKIP_RPATH=ON
    )
  LIST(APPEND SELFBUILTEXTERNALS cgns)
  ENDIF()
  LIST(APPEND linkedlibs ${HOPR_CGNS_DIR}/build/lib/libcgns.a)
  INCLUDE_DIRECTORIES(${HOPR_CGNS_DIR}/build/include)
ELSE()
  FIND_PACKAGE(CGNS)
  IF (CGNS_FOUND)
    MESSAGE(STATUS "CGNS include dir: " ${CGNS_INCLUDE_DIR})
    LIST(APPEND linkedlibs ${CGNS_LIBRARIES})
    INCLUDE_DIRECTORIES (${CGNS_INCLUDE_DIR})
 ELSE()
    MESSAGE(ERROR "CGNS not found")
  ENDIF()
ENDIF()

# =========================================================================
# HOPR sources & targets
# =========================================================================

function(add_exec target)
  ADD_EXECUTABLE(${target} ${ARGN})
  TARGET_LINK_LIBRARIES(${target} hoprlib stacksizelib ${linkedlibs})
  SET_TARGET_PROPERTIES(${target} PROPERTIES COMPILE_FLAGS ${HOPR_COMPILE_FLAGS}) 
  SET_TARGET_PROPERTIES(${target} PROPERTIES LINKER_LANGUAGE Fortran)
  INSTALL(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${target} DESTINATION bin)
  STRING(TOUPPER ${target} TARGETUP)
  ADD_CUSTOM_COMMAND(TARGET ${target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "SUCCESS: ${TARGETUP} BUILD COMPLETE!")
endfunction()

# include hopr.h
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src)


# collect all F90 files needed for hoprlib
FILE(GLOB_RECURSE hoprF90 ./src/basis/*.f90
                          ./src/mesh/*.f90
                          ./src/output/*.f90
                          ./src/readin/*.f90
                          ./src/search/*.f90
                          ./src/sfc/*.f90
                          globals.f90
                          headers.f90
                          io_hdf5.f90
                          sortingtools.f90)

FILE(GLOB_RECURSE hoprC   ./src/sfc/*.c
                          ./src/globals/stacksize.c)

# define libs
ADD_LIBRARY(hoprlibF90 OBJECT ${hoprF90})
SET_TARGET_PROPERTIES(hoprlibF90 PROPERTIES COMPILE_FLAGS ${HOPR_COMPILE_FLAGS}) 
ADD_LIBRARY(hoprlibC OBJECT ${hoprC})
SET_TARGET_PROPERTIES(hoprlibC PROPERTIES LINKER_LANGUAGE C)
SET_TARGET_PROPERTIES(hoprlibC PROPERTIES COMPILE_FLAGS "-fPIC")
ADD_LIBRARY(hoprlib STATIC $<TARGET_OBJECTS:hoprlibF90> $<TARGET_OBJECTS:hoprlibC>)
IF(SELFBUILTEXTERNALS)
  ADD_DEPENDENCIES(hoprlibF90 ${SELFBUILTEXTERNALS})
ENDIF()

ADD_LIBRARY(stacksizelib STATIC ./src/stacksize.c)
SET_TARGET_PROPERTIES(stacksizelib PROPERTIES LINKER_LANGUAGE C)
SET_TARGET_PROPERTIES(stacksizelib PROPERTIES COMPILE_FLAGS "-fPIC")

add_exec(hopr ./src/hopr.f90)
ADD_DEPENDENCIES(hopr hoprlib stacksizelib)

# =========================================================================
# Y+-Estimator Tool
# =========================================================================
IF (HOPR_BUILD_YPLUSESTIMATOR)
  add_exec(yplusestimator ./tools/yplusestimator/yplusestimator.f90)
ENDIF()

# =========================================================================
# Ctags
# =========================================================================
FIND_PROGRAM(CTAGS_PATH ctags-exuberant)
IF(CTAGS_PATH)
  SET(EXCLUDES "")
  ADD_CUSTOM_COMMAND(OUTPUT ctags.txt COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && ${CTAGS_PATH} -R --fortran-kinds=+i+L -f ctags.txt ${EXCLUDES} src)
  ADD_CUSTOM_TARGET(tags DEPENDS ctags.txt)
  MESSAGE(STATUS "Found ctags: ${CTAGS_PATH}")
  ADD_DEPENDENCIES(hopr tags)
ENDIF(CTAGS_PATH)
