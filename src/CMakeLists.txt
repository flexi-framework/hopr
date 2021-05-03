# =========================================================================
# HOPR Options
# =========================================================================

# Build directory list
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

# =========================================================================
# HOPR sources & targets
# =========================================================================
FUNCTION(add_exec target)
  ADD_EXECUTABLE(${target} ${ARGN})
  TARGET_LINK_LIBRARIES(${target} hoprlib stacksizelib ${linkedlibs})
  SET_TARGET_PROPERTIES(${target} PROPERTIES COMPILE_FLAGS ${HOPR_COMPILE_FLAGS})
  SET_TARGET_PROPERTIES(${target} PROPERTIES LINKER_LANGUAGE Fortran)
  INSTALL(PROGRAMS ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${target} DESTINATION bin)
  STRING(TOUPPER ${target} TARGETUP)
  ADD_CUSTOM_COMMAND(TARGET ${target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --cyan "SUCCESS: ${TARGETUP} BUILD COMPLETE!")
ENDFUNCTION()

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
ADD_LIBRARY(          hoprlibF90 OBJECT ${hoprF90})
SET_TARGET_PROPERTIES(hoprlibF90 PROPERTIES COMPILE_FLAGS ${HOPR_COMPILE_FLAGS})
ADD_LIBRARY(          hoprlibC OBJECT ${hoprC})
SET_TARGET_PROPERTIES(hoprlibC PROPERTIES LINKER_LANGUAGE C)
SET_TARGET_PROPERTIES(hoprlibC PROPERTIES COMPILE_FLAGS "-fPIC")
ADD_LIBRARY(          hoprlib STATIC $<TARGET_OBJECTS:hoprlibF90> $<TARGET_OBJECTS:hoprlibC>)
IF(SELFBUILTEXTERNALS)
  ADD_DEPENDENCIES(hoprlibF90 ${SELFBUILTEXTERNALS})
ENDIF()

ADD_LIBRARY(stacksizelib STATIC ./src/stacksize.c)
SET_TARGET_PROPERTIES(stacksizelib PROPERTIES LINKER_LANGUAGE C)
SET_TARGET_PROPERTIES(stacksizelib PROPERTIES COMPILE_FLAGS "-fPIC")

ADD_EXEC(        hopr ./src/hopr.f90)
ADD_DEPENDENCIES(hopr hoprlib stacksizelib)

# =========================================================================
# Y+-Estimator Tool
# =========================================================================
OPTION(HOPR_YPLUSESTIMATOR "Compile and build Y+-Estimator" OFF)
IF (HOPR_YPLUSESTIMATOR)
  ADD_EXEC(yplusestimator ./tools/yplusestimator/yplusestimator.f90)
ENDIF()

# =========================================================================
# Ctags
# =========================================================================
MARK_AS_ADVANCED(FORCE CTAGS_PATH)
FIND_PROGRAM(CTAGS_PATH ctags-exuberant)
IF(CTAGS_PATH)
  SET(EXCLUDES "")
  ADD_CUSTOM_COMMAND(OUTPUT ctags.txt COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} && ${CTAGS_PATH} -R --fortran-kinds=+i+L -f ctags.txt ${EXCLUDES} src)
  ADD_CUSTOM_TARGET(tags DEPENDS ctags.txt)
  MESSAGE(STATUS "Found ctags: ${CTAGS_PATH}")
  ADD_DEPENDENCIES(hopr tags)
ENDIF(CTAGS_PATH)
