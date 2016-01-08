CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11)
PROJECT(Hopr)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
ENABLE_LANGUAGE(C Fortran)

# =========================================================================
# Options
# =========================================================================
OPTION(HOPR_MPI "Compile SINGLE or MPI version" ON)
OPTION(HOPR_HDF5F90 "enable old hdf5 interface" OFF)
SET(HOPR_CGNS_INT "32" CACHE STRING "integer type in CGNS lib")
ADD_DEFINITIONS(-DPP_CGNS_INT=${HOPR_CGNS_INT})
# =========================================================================
# Machine environements
# =========================================================================
STRING(SUBSTRING ${CMAKE_ROOT} 4 6 CMAKEHOST)
STRING(COMPARE EQUAL ${CMAKEHOST} "hornet" ONHORNET)
STRING(COMPARE EQUAL ${CMAKEHOST} "nehalem" ONNEHALEM)
IF (ONHORNET OR ONNEHALEM)
  ADD_DEFINITIONS(-DLUSTRE)
ENDIF()
STRING(COMPARE EQUAL ${CMAKEHOST} "juqueen" ONJUQUEEN)
IF (ONJUQUEEN)
  ADD_DEFINITIONS(-DGPFS)
ENDIF()

# =========================================================================
# HOPR sources & targets
# =========================================================================

# include hopr.h 
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/src)

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
ADD_DEPENDENCIES(hopr hdf5 cgns)


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

IF (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  # set compiler
  SET (HOPR_COMPILERNAME "GNU")
  IF (HOPR_MPI)
     SET (CMAKE_Fortran_COMPILER mpif90)
     SET (CMAKE_C_COMPILER mpicc)
     SET (CMAKE_CXX_COMPILER mpicc)
  ENDIF()
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -fdefault-double-8 -fbackslash -DGNU")
  SET (CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native  -finline-functions")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -ggdb3 -fbounds-check -finit-real=nan -fbacktrace  -Wall")
  # add flags only for compiling not linking!
  SET_TARGET_PROPERTIES(hoprlibF90 PROPERTIES COMPILE_FLAGS "-xf95-cpp-input")
  SET_TARGET_PROPERTIES(hopr PROPERTIES COMPILE_FLAGS "-xf95-cpp-input")
ELSEIF (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort
  # set compiler
  SET (HOPR_COMPILERNAME "INTEL")
  IF (HOPR_MPI)
     SET (CMAKE_Fortran_COMPILER mpif90)
     SET (CMAKE_C_COMPILER mpicc)
     SET (CMAKE_CXX_COMPILER mpicc)
  ENDIF()
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -r8 -i4 -traceback -warn all -DINTEL")
  SET (CMAKE_Fortran_FLAGS_RELEASE "-O2 $(INSTRUCTION) -vec-report0 -no-prec-div")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -fpe0 -traceback -check all,noarg_temp_created,noformat,nooutput_conversion,pointer,uninit")
  SET (CMAKE_C_FLAGS_RELEASE "-O2 -xHost -ip")
  SET (CMAKE_C_FLAGS_DEBUG "-g")
  SET (CMAKE_CXX_FLAGS_RELEASE "-O2 -xHost -ip")
  SET (CMAKE_CXX_FLAGS_DEBUG "-g")
  # add flags only for compiling not linking!
  SET_TARGET_PROPERTIES(hoprlibF90 PROPERTIES COMPILE_FLAGS "-fpp -assume bscc")
  SET_TARGET_PROPERTIES(hopr PROPERTIES COMPILE_FLAGS "-fpp -assume bscc")
ELSE()
   MESSAGE(ERROR "Unknown compiler")
ENDIF()

# =========================================================================
# MPI
# =========================================================================
IF(HOPR_MPI)
  FIND_PACKAGE(MPI REQUIRED)
  LIST(GET MPI_INCLUDE_PATH 0 MPI_INCLUDE_PATH)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
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
#OPTION(HOPR_BUILDHDF5 "Compile and build hdf5 library" OFF)
SET(HOPR_HDF5DIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${HOPR_COMPILERNAME}/hdf5)
FIND_PACKAGE (HDF5 NAMES hdf5 QUIET)
IF (NOT(HDF5_FOUND))
  OPTION(HOPR_BUILDHDF5 "Compile and build hdf5 library" ON)
ELSE()
  OPTION(HOPR_BUILDHDF5 "Compile and build hdf5 library" OFF)
ENDIF()
IF (HOPR_BUILDHDF5)
  SET(HOPR_HDF5VERSION 1.8.16 CACHE STRING "Sets the version of the hdf5 library")
  IF (EXISTS "${HOPR_HDF5DIR}/build/lib/libhdf5.a")
  ELSE()
    include(ExternalProject)
    ExternalProject_Add(hdf5
      URL http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-${HOPR_HDF5VERSION}/src/hdf5-${HOPR_HDF5VERSION}.tar.bz2
      PREFIX ${HOPR_HDF5DIR}
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${HOPR_HDF5DIR}/build /
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_C_FLAGS=${CMAKE_C_FLAGS_RELEASE} /
      -DHDF5_BUILD_FORTRAN=ON -DHDF5_ENABLE_PARALLEL=${HOPR_MPI} -DHDF5_BUILD_CPP_LIB=OFF -DHDF5_BUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF /
      -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${HOPR_HDF5DIR}/build
    )
  ENDIF()
  SET(HDF5_DIR  ${HOPR_HDF5DIR}/build/lib)
  SET(HDF5_INCLUDE_DIR ${HOPR_HDF5DIR}/build/include/static)
  MESSAGE(STATUS "HDF5 include dir: " ${HDF5_INCLUDE_DIR})
  INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${HDF5_INCLUDE_DIR}")
  TARGET_LINK_LIBRARIES(hopr ${HDF5_DIR}/libhdf5_fortran.a ${HDF5_DIR}/libhdf5_f90cstub.a ${HDF5_DIR}/libhdf5.a -ldl -lz)
ELSE()
  FIND_PACKAGE (HDF5 NAMES hdf5)
  LIST(GET HDF5_INCLUDE_DIR 0 HDF5_INCLUDE_DIR)
  MESSAGE(STATUS "HDF5 include dir: " ${HDF5_INCLUDE_DIR})
  INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${HDF5_INCLUDE_DIR}")
  TARGET_LINK_LIBRARIES(hopr ${HDF5_LIBRARIES})
ENDIF()
IF (HOPR_HDF5F90)
  ADD_DEFINITIONS(-DHDF5_F90)
ENDIF()

# =========================================================================
# CGNS library
# =========================================================================
SET(HOPR_CGNSDIR ${CMAKE_CURRENT_SOURCE_DIR}/share/${HOPR_COMPILERNAME}/cgns)
FIND_PACKAGE(CGNS QUIET)
IF (NOT(CGNS_FOUND))
  OPTION(HOPR_BUILDCGNS "Compile and build cgns library" ON)
ELSE()
  OPTION(HOPR_BUILDCGNS "Compile and build cgns library" OFF)
ENDIF()

IF (HOPR_BUILDCGNS)
  SET(HOPR_CGNSVERSION "3.2.1" CACHE STRING "CGNS version")
  IF (NOT(EXISTS "${HOPR_CGNSDIR}/build/lib/libcgns.a"))
    include(ExternalProject)
    ExternalProject_Add(cgns
    URL https://github.com/CGNS/CGNS/archive/v${HOPR_CGNSVERSION}.tar.gz
    PATCH_COMMAND sed -i s/^c/!/g ${HOPR_CGNSDIR}/src/cgns/src/cgnslib_f.h.in
    PREFIX ${HOPR_CGNSDIR}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${HOPR_CGNSDIR}/build /
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} /
    -DCMAKE_BUILD_TYPE=Release -DCGNS_ENABLE_FORTRAN=ON -DCMAKE_INSTALL_PREFIX=${HOPR_CGNSDIR}/build
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
# RECIPES library
# =========================================================================
#FIND_LIBRARY(HOPR_RECIPESLIB recipes ${CMAKE_CURRENT_SOURCE_DIR}/share/${HOPR_COMPILERNAME}/RECIPES/)
#TARGET_LINK_LIBRARIES(hopr ${HOPR_RECIPESLIB})

# =========================================================================
# LAPACK library
# =========================================================================
IF(NOT LAPACK_lapack_LIBRARY)
   FIND_PACKAGE(LAPACK)
ENDIF()
TARGET_LINK_LIBRARIES(hopr ${LAPACK_lapack_LIBRARY}) 
TARGET_LINK_LIBRARIES(hoprlibF90 ${LAPACK_lapack_LIBRARY}) 

