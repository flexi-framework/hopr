CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
IF ("${CMAKE_VERSION}" VERSION_GREATER 3.0.0)
  cmake_policy (SET CMP0046 OLD)
ENDIF()
PROJECT(Hopr)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# =========================================================================
# Options
# =========================================================================
include(ExternalProject)

OPTION(HOPR_MPI "Compile SINGLE or MPI version" ON)
OPTION(HOPR_HDF5F90 "enable old hdf5 interface" OFF)
SET(HOPR_CGNS_INT "32" CACHE STRING "integer type in CGNS lib")
ADD_DEFINITIONS(-DPP_CGNS_INT=${HOPR_CGNS_INT})
ENABLE_LANGUAGE(C CXX Fortran)

# =========================================================================
# Machine environements and specific library paths
# =========================================================================
SET (CGNSDOWNLOAD https://github.com/CGNS/CGNS.git )
SET (HDF5DOWNLOAD http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-${HOPR_HDF5VERSION}/src/hdf5-${HOPR_HDF5VERSION}.tar.gz )

SITE_NAME(CMAKE_HOSTNAME)
STRING(SUBSTRING "${CMAKE_HOSTNAME}" 0 7 CMAKEHOST)

STRING(COMPARE EQUAL "${CMAKEHOST}" "eslogin" ONHORNET)
IF (ONHORNET OR ONNEHALEM)
  ADD_DEFINITIONS(-DLUSTRE)
  SET (CGNSDOWNLOAD git@129.69.43.151:libs/CGNS.git )
ENDIF()

#STRING(COMPARE EQUAL "${CMAKEHOST}" "juqueen" ONJUQUEEN)
#IF (ONJUQUEEN)
  #ADD_DEFINITIONS(-DGPFS)
#ENDIF()

# =========================================================================
# HOPR sources & targets
# =========================================================================

# include hopr.h 
FILE(COPY ${CMAKE_CURRENT_SOURCE_DIR}/src/defines.f90 DESTINATION ${CMAKE_Fortran_MODULE_DIRECTORY})
INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})

# collect all F90 files needed for hoprlib
FILE(GLOB_RECURSE hoprF90 ./src/basis/*.f90 ./src/mesh/*.f90 ./src/output/*.f90 ./src/readin/*.f90 ./src/search/*.f90 ./src/sfc/*.f90 globals.f90 headers.f90 io_hdf5.f90 sortingtools.f90)
FILE(GLOB_RECURSE hoprC ./src/sfc/*.c)

# define hopr library 
ADD_LIBRARY(hoprlibF90 STATIC ${hoprF90})
SET_TARGET_PROPERTIES(hoprlibF90 PROPERTIES LINKER_LANGUAGE Fortran)
TARGET_LINK_LIBRARIES(hoprlibF90 ${LINK_LIBS})
ADD_LIBRARY(hoprlibC STATIC ${hoprC})
SET_TARGET_PROPERTIES(hoprlibC PROPERTIES LINKER_LANGUAGE C)
TARGET_LINK_LIBRARIES(hoprlibC ${LINK_LIBS})

# define hopr executable
ADD_EXECUTABLE(hopr ./src/hopr.f90)
TARGET_LINK_LIBRARIES(hopr hoprlibF90 hoprlibC) 
ADD_DEPENDENCIES(hopr HDF5 cgns)
INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/bin/hopr DESTINATION bin)
INSTALL(FILES  ${CMAKE_CURRENT_BINARY_DIR}/bin/configuration.cmake DESTINATION bin)

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
# FLAGS
# =========================================================================

# append relative filename-macro for __FILENAME__ in Stamp of abort function (see globals.f90)
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D__FILENAME__='\"$(subst ${CMAKE_SOURCE_DIR}/,,$(abspath $<))\"'")

# =========================================================================
# COMPILER STUFF
# =========================================================================

# FFLAGS depend on the compiler
GET_FILENAME_COMPONENT (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # set compiler
  SET (HOPR_COMPILERNAME "GNU")
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -fbackslash -DGNU")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -march=native  -finline-functions")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -ggdb3 -fbounds-check -finit-real=nan -fbacktrace  -Wall")
  # add flags only for compiling not linking!
  SET_TARGET_PROPERTIES(hoprlibF90 PROPERTIES COMPILE_FLAGS "-xf95-cpp-input")
  SET_TARGET_PROPERTIES(hopr       PROPERTIES COMPILE_FLAGS "-xf95-cpp-input")
ELSEIF (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # set compiler
  SET (HOPR_COMPILERNAME "INTEL")
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -r8 -i4 -traceback -warn all -DINTEL -mkl=sequential")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O2 $(INSTRUCTION) -vec-report0 -no-prec-div")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -fpe0 -traceback -check all,noarg_temp_created,noformat,nooutput_conversion,pointer,uninit")
  SET (CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS} -O2 -xHost -ip")
  SET (CMAKE_C_FLAGS_DEBUG     "${CMAKE_C_FLAGS} -g")
  SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O2 -xHost -ip")
  SET (CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS} -g")
  # add flags only for compiling not linking!
  SET_TARGET_PROPERTIES(hoprlibF90 PROPERTIES COMPILE_FLAGS "-fpp -assume bscc")
  SET_TARGET_PROPERTIES(hopr       PROPERTIES COMPILE_FLAGS "-fpp -assume bscc")
ELSE()
   MESSAGE(ERROR "Unknown compiler")
ENDIF()

# =========================================================================
# MPI
# =========================================================================
IF(HOPR_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  LIST(GET MPI_INCLUDE_PATH 0 MPI_INCLUDE_PATH)
  #INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
  MESSAGE(STATUS "MPI include dir: " ${MPI_INCLUDE_PATH})
  ADD_DEFINITIONS(-DMPI)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${MPI_INCLUDE_PATH}")
  SET (LINK_LIBS ${LINK_LIBS} ${MPI_LIBRARIES})
  TARGET_LINK_LIBRARIES(hoprlibF90 ${MPI_Fortran_LIBRARIES})
  TARGET_LINK_LIBRARIES(hoprlibC   ${MPI_C_LIBRARIES})
  TARGET_LINK_LIBRARIES(hopr       ${MPI_Fortran_LIBRARIES})
  MESSAGE(STATUS "Building HOPR with MPI")
ENDIF()

# =========================================================================
# HDF5 library
# =========================================================================
#HDF5
SET(HOPR_HDF5DIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${HOPR_COMPILERNAME}/hdf5)
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

  SET(HOPR_HDF5VERSION 1.8.16 CACHE STRING "Sets the version of the hdf5 library")
  SET(HDF5_DIR  ${HOPR_HDF5DIR}/build/lib)
  IF ("${CMAKE_VERSION}" VERSION_GREATER 3.2.0) # Compile HDF5 using cmake
    IF (NOT EXISTS "${HOPR_HDF5DIR}/build/lib/libhdf5.a")
      ExternalProject_Add(hdf5
        DOWNLOAD_COMMAND ${HDF5DOWNLOAD}
        PREFIX ${HOPR_HDF5DIR}
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${HOPR_HDF5DIR}/build /
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS_RELEASE} /
        -DHDF5_BUILD_FORTRAN=ON -DHDF5_ENABLE_PARALLEL=${HOPR_MPI} -DHDF5_BUILD_CPP_LIB=OFF -DHDF5_BUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF /
        -DBUILD_SHARED_LIBS=OFF -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_BUILD_TOOLS=OFF -DHDF5_BUILD_HL_LIB=OFF /
        -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${HOPR_HDF5DIR}/build
      )
    ENDIF()
    SET(HDF5_INCLUDE_DIR ${HOPR_HDF5DIR}/build/include/static)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${HDF5_INCLUDE_DIR}")
    TARGET_LINK_LIBRARIES(hopr ${HDF5_DIR}/libhdf5_fortran.a ${HDF5_DIR}/libhdf5_f90cstub.a ${HDF5_DIR}/libhdf5.a -ldl -lz)
  ELSE()                                        # Compile HDF5 using automake
    IF (NOT EXISTS "${HOPR_HDF5DIR}/lib/libhdf5.a")
      ExternalProject_Add(hdf5
        DOWNLOAD_COMMAND ${HDF5DOWNLOAD}
        PREFIX ${HOPR_HDF5DIR}
        CONFIGURE_COMMAND FC=${CMAKE_FORTRAN_COMPULER} CC=${CMAKE_C_COMPILER} ${HOPR_HDF5DIR}/src/hdf5/configure --prefix=${HOPR_HDF5DIR}/build --enable-fortran --disable-shared
        BUILD_COMMAND ${MAKE}
      )
    ENDIF()
    SET(HDF5_INCLUDE_DIR ${HOPR_HDF5DIR}/build/include)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${HDF5_INCLUDE_DIR}")
    TARGET_LINK_LIBRARIES(hopr ${HDF5_DIR}/libhdf5_fortran.a ${HDF5_DIR}/libhdf5.a -ldl -lz)
  ENDIF()
  MESSAGE(STATUS "HDF5 include dir: " ${HDF5_INCLUDE_DIR})
  INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})

ELSE()

  IF(HDF5_DIR STREQUAL "HDF5_DIR-NOTFOUND")     # HDF5 built using standard configure script
    FIND_PACKAGE (HDF5 COMPONENTS C Fortran )
    LIST(GET HDF5_INCLUDE_DIR 0 HDF5_INCLUDE_DIR)
    MESSAGE(STATUS "HDF5 include dir: " ${HDF5_INCLUDE_DIR})
    INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${HDF5_INCLUDE_DIR}")
    TARGET_LINK_LIBRARIES(hopr ${HDF5_LIBRARIES})
    MESSAGE(STATUS "HDF5 libraries: " ${HDF5_LIBRARIES})
  ELSE()                                        # HDF5 built using cmake
    FIND_PACKAGE (HDF5 NAMES hdf5 COMPONENTS C Fortran static )
    LIST(GET HDF5_INCLUDE_DIR 0 HDF5_INCLUDE_DIR)
    MESSAGE(STATUS "HDF5 include dir: " ${HDF5_INCLUDE_DIR_FORTRAN})
    INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR_FORTRAN})
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${HDF5_INCLUDE_DIR_FORTRAN}")
    TARGET_LINK_LIBRARIES(hopr ${HDF5_FORTRAN_STATIC_LIBRARY} ${HDF5_LIBRARIES})
    MESSAGE(STATUS "HDF5 libraries: " ${HDF5_LIBRARIES} ${HDF5_FORTRAN_STATIC_LIBRARY})
  ENDIF()
  IF (NOT HDF5_FOUND)
    MESSAGE(ERROR "HDF5 not found and HOPR_BUILD_HDF5=OFF!")
  ENDIF()

ENDIF()
IF (HOPR_HDF5F90)
  ADD_DEFINITIONS(-DHDF5_F90)
ENDIF()

# =========================================================================
# CGNS library
# =========================================================================
SET(HOPR_CGNSDIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${HOPR_COMPILERNAME}/cgns)
FIND_PACKAGE(CGNS QUIET)
IF (NOT CGNS_FOUND)
  OPTION(HOPR_BUILDCGNS "Compile and build cgns library" ON)
ELSE()
  OPTION(HOPR_BUILDCGNS "Compile and build cgns library" OFF)
ENDIF()

STRING(COMPARE EQUAL ${HOPR_CGNS_INT} "64" HOPR_CGNS_64BIT)

IF (HOPR_BUILDCGNS)
  SET(HOPR_CGNSTAG "v3.2.1" CACHE STRING "CGNS version tag")
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
# LAPACK library
# =========================================================================
FIND_PACKAGE(LAPACK REQUIRED)
TARGET_LINK_LIBRARIES(hopr ${LAPACK_LIBRARIES}) 
TARGET_LINK_LIBRARIES(hoprlibF90 ${LAPACK_LIBRARIES})

# =========================================================================
# store configuratioPREFIX ${HOPR_HDF5DIR}n
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
   SET(boolian ${${_var}} STREQUAL "ON" OR ${${_var}} STREQUAL "OFF")
   IF(${boolian})
      FILE(APPEND ${configuration} "OPTION(${_var} \"${currentHelpString}\" ${${_var}})\n")
   ELSE()
      IF(${_var})
         STRING(REPLACE "\\" "\\\\" ${_var} ${${_var}})
      ENDIF()
      FILE(APPEND ${configuration} "SET(${_var} \"${${_var}}\" CACHE STRING \"${currentHelpString}\")\n")
   ENDIF()
ENDFOREACH()
