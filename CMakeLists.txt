CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
IF ("${CMAKE_VERSION}" VERSION_GREATER 3.0.0)
  cmake_policy (SET CMP0046 OLD)
ENDIF()
# =========================================================================
# Detect machine environements
# =========================================================================
SITE_NAME(CMAKE_HOSTNAME)
STRING(FIND "${CMAKE_HOSTNAME}" "eslogin" ONHORNET)
# Set options that need to be enabled BEFORE ANY PROJECT OR LANGUAGE COMMAND
IF (ONHORNET GREATER -1)
  SET(CMAKE_C_COMPILER       cc)
  SET(CMAKE_Fortran_COMPILER ftn)
  MESSAGE(STATUS "Compiling on Hazelhen")
  ADD_DEFINITIONS(-DLUSTRE)
ENDIF()

# =========================================================================
# Machine environements and specific library paths
# =========================================================================
SET (CGNSDOWNLOAD https://github.com/CGNS/CGNS.git )
SET (HOPR_CGNSTAG "v3.2.1" CACHE STRING "CGNS version tag")
SET (HDF5DOWNLOAD https://gitlabext.iag.uni-stuttgart.de/libs/HDF5.git )
SET (HOPR_HDF5TAG "1.8.17" CACHE STRING "HDF5 version tag")

IF (ONHORNET GREATER -1)
  ADD_DEFINITIONS(-DLUSTRE)
  SET (CGNSDOWNLOAD git@129.69.43.151:libs/CGNS.git )
ENDIF()

# =========================================================================
# Project definitions and output paths
# =========================================================================
PROJECT(Hopr)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
ENABLE_LANGUAGE(C Fortran)
INCLUDE(ExternalProject)

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  SET(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Prefix prepended to install directories" FORCE)
ENDIF()

# =========================================================================
# Options
# =========================================================================
OPTION(HOPR_MPI "Compile SINGLE or MPI version" ON)
OPTION(HOPR_HDF5F90 "enable old hdf5 interface" OFF)
SET(HOPR_CGNS_INT "32" CACHE STRING "integer type in CGNS lib")
ADD_DEFINITIONS(-DPP_CGNS_INT=${HOPR_CGNS_INT})

# =========================================================================
# USERBLOCK + Preproc_flags
# =========================================================================

# A function to get all user defined variables with a specified prefix
function (getListOfVarsStartingWith _prefix _varResult)
    get_cmake_property(_vars CACHE_VARIABLES)
    string (REGEX MATCHALL "(^|;)${_prefix}[A-Za-z0-9_]*" _matchedVars "${_vars}")
    set (${_varResult} ${_matchedVars} PARENT_SCOPE)
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

# =========================================================================
# HOPR sources & targets
# =========================================================================

# include hopr.h 
FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/src/hopr.h DESTINATION ${CMAKE_Fortran_MODULE_DIRECTORY})
INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})

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

FILE(GLOB_RECURSE hoprC   ./src/sfc/*.c)

# define hopr library 
ADD_LIBRARY(hoprlibF90 STATIC ${hoprF90})
SET_TARGET_PROPERTIES(hoprlibF90 PROPERTIES LINKER_LANGUAGE Fortran)

ADD_LIBRARY(hoprlibC STATIC ${hoprC})
SET_TARGET_PROPERTIES(hoprlibC PROPERTIES LINKER_LANGUAGE C)

# define hopr executable
ADD_EXECUTABLE(hopr ./src/hopr.f90)
TARGET_LINK_LIBRARIES(hopr hoprlibF90 hoprlibC) 
ADD_DEPENDENCIES(hopr HDF5 cgns)

# =========================================================================
# Additional tools for hopr
# =========================================================================
# Y+-Estimator
OPTION(HOPR_BUILDYPLUSESTIMATOR "Compile and build Y+-Estimator" OFF)
IF (HOPR_BUILDYPLUSESTIMATOR)
  ADD_EXECUTABLE(yplusestimator ./tools/yplusestimator/yplusestimator.f90)
  TARGET_LINK_LIBRARIES(yplusestimator hoprlibF90)
  SET_TARGET_PROPERTIES(yplusestimator PROPERTIES LINKER_LANGUAGE Fortran)
  INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/bin/yplusestimator DESTINATION bin)
  ADD_DEPENDENCIES(yplusestimator hoprlibF90)
ENDIF()

# =========================================================================
# Location of binary and filenames
# =========================================================================
# append relative filename-macro for __FILENAME__ in Stamp of abort function (see hopr.h)
ADD_DEFINITIONS("-D__FILENAME__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'") 
# add basedir for location of corresponding userblock-file
ADD_DEFINITIONS("-DBASEDIR='\"${CMAKE_CURRENT_BINARY_DIR}/\"'")

# =========================================================================
# COMPILER STUFF
# =========================================================================

# FFLAGS depend on the compiler
GET_FILENAME_COMPONENT (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # gfortran
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -fbackslash -ffree-line-length-0 -DGNU")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -march=native  -finline-functions")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -pg -O3 -march=native  -finline-functions")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -ggdb3 -fbounds-check -finit-real=nan -fbacktrace  -Wall")
  # add flags only for compiling not linking!
  SET_TARGET_PROPERTIES(hopr hoprlibF90 PROPERTIES COMPILE_FLAGS "-xf95-cpp-input")
ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
   # ifort
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -r8 -i4 -traceback -warn all -DINTEL")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O2 -xHost -qopt-report0 -qopt-report-phase=vec")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -p -O2 xHost -qopt-report0 -qopt-report-phase=vec")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -fpe0 -traceback -check all,noarg_temp_created,noformat,nooutput_conversion,pointer,uninit")
  SET (CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS} -O2 -xHost -ip")
  SET (CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS} -p -O2 -xHost -ip")
  SET (CMAKE_C_FLAGS_DEBUG     "${CMAKE_C_FLAGS} -g")
  # add flags only for compiling not linking!
  SET_TARGET_PROPERTIES(hopr hoprlibF90 PROPERTIES COMPILE_FLAGS "-fpp -assume bscc")
ELSE()
   MESSAGE(ERROR "Unknown compiler")
ENDIF()

# =========================================================================
# HDF5 library
# =========================================================================
#HDF5
SET(HOPR_HDF5DIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${CMAKE_Fortran_COMPILER_ID}/hdf5)
SET(HDF5_USE_STATIC_LIBRARIES TRUE)
FIND_PACKAGE (HDF5 NAMES hdf5 COMPONENTS C Fortran QUIET)
IF (NOT HDF5_FOUND)
  FIND_PACKAGE (HDF5 COMPONENTS C Fortran QUIET)
ENDIF()
IF (NOT HDF5_FOUND)
  OPTION(HOPR_BUILDHDF5 "Compile and build hdf5 library" ON)
ELSE()
  OPTION(HOPR_BUILDHDF5 "Compile and build hdf5 library" OFF)
ENDIF()

IF(HOPR_BUILDHDF5)

  SET(HDF5_DIR  ${HOPR_HDF5DIR}/build/lib)
  IF (NOT EXISTS "${HOPR_HDF5DIR}/build/lib/libhdf5.a")
    IF(HOPR_MPI)
      SET(HOPR_HDF5PARALLEL --enable-parallel)
      SET(HOPR_HDF5FC ${MPI_Fortran_COMPILER})
      SET(HOPR_HDF5CC ${MPI_C_COMPILER})
    ELSE()
      SET(HOPR_HDF5FC ${CMAKE_Fortran_COMPILER})
      SET(HOPR_HDF5CC ${CMAKE_C_COMPILER} )
    ENDIF()
    ExternalProject_Add(hdf5
      GIT_REPOSITORY ${HDF5DOWNLOAD}
      PREFIX ${HOPR_HDF5DIR}
      CONFIGURE_COMMAND FC=${HOPR_HDF5FC} CC=${HOPR_HDF5CC} ${HOPR_HDF5DIR}/src/hdf5/configure --prefix=${HOPR_HDF5DIR}/build --enable-fortran --enable-fortran2003 ${HOPR_HDF5PARALLEL} --disable-shared
      BUILD_COMMAND ${MAKE}
    )
  ENDIF()
  SET(HDF5_INCLUDE_DIR ${HOPR_HDF5DIR}/build/include)
  TARGET_LINK_LIBRARIES(hopr ${HDF5_DIR}/libhdf5_fortran.a ${HDF5_DIR}/libhdf5.a -ldl -lz)
  MESSAGE(STATUS "HDF5 include dir: " ${HDF5_INCLUDE_DIR})
  INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})

ELSE()
  SET(HOPR_HDF5_DIR $ENV{HDF5_DIR} CACHE PATH "Directory of HDF5")
  IF(EXISTS "${HOPR_HDF5_DIR}/include/")
    INCLUDE_DIRECTORIES (${HOPR_HDF5_DIR}/include)
    MESSAGE(STATUS "HDF5 include dir: " ${HOPR_HDF5_DIR}/include)
    SET(HOPR_HDF5_FOUND TRUE)
  ENDIF()
  IF(EXISTS "${HOPR_HDF5_DIR}/include/static")
    INCLUDE_DIRECTORIES (${HOPR_HDF5_DIR}/include/static)
    MESSAGE(STATUS "HDF5 include dir: " ${HOPR_HDF5_DIR}/include/static)
    SET(HOPR_HDF5_FOUND TRUE)
  ENDIF()
  IF(NOT HOPR_HDF5_FOUND)
    MESSAGE(ERROR "No hdf5 installation found!")
  ENDIF()
  FIND_LIBRARY(HOPR_HDF51 NAMES libhdf5_fortran.a hdf5_fortran PATHS  ${HOPR_HDF5_DIR}/lib)
  TARGET_LINK_LIBRARIES(hopr   ${HOPR_HDF51})
  IF(EXISTS "${HOPR_HDF5_DIR}/lib/libhdf5_hl_f90cstub.a")
    FIND_LIBRARY(HOPR_HDF52 NAMES libhdf5_f90cstub.a hdf5_f90cstub PATHS ${HOPR_HDF5_DIR}/lib)
    TARGET_LINK_LIBRARIES(hopr   ${HOPR_HDF52})
  ENDIF()
  FIND_LIBRARY(HOPR_HDF53 NAMES libhdf5.a hdf5 PATHS ${HOPR_HDF5_DIR}/lib)
  TARGET_LINK_LIBRARIES(hopr   ${HOPR_HDF53} -ldl -lz)
  MESSAGE(STATUS "Building hopr with parallel HDF5 (HyperSlab)") 
  IF (ONHORNET GREATER -1)
     # customized location of the tools on hazelhen
    ADD_DEFINITIONS(-DH5TOOLSDIR=\"${HOPR_HDF5_DIR}/../../../../hdf5/1.8.14/bin/\")
  ELSE()
    ADD_DEFINITIONS(-DH5TOOLSDIR=\"${HOPR_HDF5_DIR}/bin/\")
  ENDIF()

  IF (NOT HDF5_FOUND)
    MESSAGE(ERROR "HDF5 not found and HOPR_BUILD_HDF5=OFF!")
  ENDIF()

ENDIF()
IF (HOPR_HDF5F90)
  ADD_DEFINITIONS(-DHDF5_F90)
ENDIF()

# =========================================================================
# LAPACK/MKL library
# =========================================================================
IF (NOT HOPR_MKL)
  # use Lapack/Blas for GNU
  MESSAGE(STATUS "Using BLAS/Lapack library")
  FIND_PACKAGE(LAPACK)
  IF (LAPACK_FOUND)
    TARGET_LINK_LIBRARIES(hopr ${LAPACK_LIBRARIES}) 
  ENDIF()
ELSE()
  MESSAGE(STATUS "Using MKL library")
  # include custom FindMKL CMake file
  LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
  FIND_PACKAGE(MKL REQUIRED)
  INCLUDE_DIRECTORIES(${MKL_INCLUDE_DIR})
  LINK_DIRECTORIES(${MKL_LIBRARY})
  TARGET_LINK_LIBRARIES(hopr ${MKL_LIBRARY})
ENDIF()

# =========================================================================
# MPI
# =========================================================================
IF(HOPR_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  LIST(GET MPI_INCLUDE_PATH 0 MPI_INCLUDE_PATH)
  IF(EXISTS "${MPI_INCLUDE_PATH}")
    INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH} ${MPI_INCLUDE_PATH}/../lib)
    MESSAGE(STATUS "MPI include dir: " ${MPI_INCLUDE_PATH})
  ENDIF()
  ADD_DEFINITIONS(-DMPI)
  TARGET_LINK_LIBRARIES(hopr       ${MPI_Fortran_LIBRARIES} ${MPI_C_LIBRARIES})
  MESSAGE(STATUS "Building HOPR with MPI")
  MESSAGE(STATUS "MPI Compiler: "  ${MPI_Fortran_COMPILER})
ENDIF()

# =========================================================================
# CGNS library
# =========================================================================
SET(HOPR_CGNSDIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${HOPR_COMPILERNAME}/cgns)
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
FIND_PACKAGE(CGNS QUIET)
IF (NOT CGNS_FOUND)
  OPTION(HOPR_BUILDCGNS "Compile and build cgns library" ON)
ELSE()
  OPTION(HOPR_BUILDCGNS "Compile and build cgns library" OFF)
ENDIF()

STRING(COMPARE EQUAL ${HOPR_CGNS_INT} "64" HOPR_CGNS_64BIT)

IF (HOPR_BUILDCGNS)
  IF (NOT EXISTS "${HOPR_CGNSDIR}/build/lib/libcgns.a")
    ExternalProject_Add(cgns
    GIT_REPOSITORY ${CGNSDOWNLOAD}
    GIT_TAG ${HOPR_CGNSTAG}
    PATCH_COMMAND sed -i s/^c/!/g ${HOPR_CGNSDIR}/src/cgns/src/cgnslib_f.h.in
    PREFIX ${HOPR_CGNSDIR}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${HOPR_CGNSDIR}/build /
    -DCGNS_ENABLE_FORTRAN=ON   -DCGNS_ENABLE_64BIT=${HOPR_CGNS_64BIT} -DCGNS_BUILD_SHARED=OFF -DCGNS_USE_SHARED=OFF /
    -DCMAKE_BUILD_TYPE=Release -DCGNS_BUILD_CGNSTOOLS=OFF -DCGNS_ENABLE_TESTS=OFF -DCMAKE_SKIP_RPATH=ON
    )
  ENDIF()
  INCLUDE_DIRECTORIES(${HOPR_CGNSDIR}/build/include)
  TARGET_LINK_LIBRARIES(hopr ${HOPR_CGNSDIR}/build/lib/libcgns.a)
ELSE()
  INCLUDE_DIRECTORIES (${CGNS_INCLUDE_DIR})
  MESSAGE(STATUS "CGNS include dir: " ${CGNS_INCLUDE_DIR})
  TARGET_LINK_LIBRARIES(hopr ${CGNS_LIBRARIES})
ENDIF()

# =========================================================================
# Userblock
# =========================================================================
ADD_CUSTOM_COMMAND(TARGET hoprlibF90 PRE_BUILD COMMAND touch ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET hoprlibF90 PRE_BUILD COMMAND rm ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET hopr POST_BUILD COMMAND echo \"{[\( CMAKE \)]}\" > ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET hopr POST_BUILD COMMAND cat ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/configuration.cmake >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET hopr POST_BUILD COMMAND
  echo \"{[\( GIT REVISIONS \)]}\" >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET hopr POST_BUILD COMMAND
  cd ${CMAKE_CURRENT_SOURCE_DIR} && git rev-parse --abbrev-ref HEAD >>  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt &&
  git log --oneline -1  >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET hopr POST_BUILD COMMAND
  cd ${CMAKE_CURRENT_SOURCE_DIR} && echo "master: " >>  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt &&
  git log origin/master --oneline -1  >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET hopr POST_BUILD COMMAND
  echo \"{[\( GIT DIFF \)]}\" >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET hopr POST_BUILD COMMAND
  cd ${CMAKE_CURRENT_SOURCE_DIR} && git diff -p origin/master >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET hopr POST_BUILD COMMAND echo \"{[\( GIT URL \)]}\" >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)
ADD_CUSTOM_COMMAND(TARGET hopr POST_BUILD COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && git config --get remote.origin.url >> ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt)


# =========================================================================
# Install
# =========================================================================
INSTALL(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/hopr DESTINATION bin)
INSTALL(FILES  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/configuration.cmake ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/userblock.txt DESTINATION bin)

# =========================================================================
# Ctags
# =========================================================================
ADD_CUSTOM_COMMAND(OUTPUT ctags.txt COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && ctags -R -f ctags.txt src)
ADD_CUSTOM_TARGET(tags DEPENDS ctags.txt)
FIND_PROGRAM(CTAGS_PATH ctags)
IF(CTAGS_PATH)
  MESSAGE(STATUS "Found ctags: ${CTAGS_PATH}")
  ADD_DEPENDENCIES(hopr tags)
ENDIF(CTAGS_PATH)

# =========================================================================
# final message
# =========================================================================
ADD_CUSTOM_COMMAND(TARGET hopr POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "SUCCESS: HOPR BUILD COMPLETE!")
