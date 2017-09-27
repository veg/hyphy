cmake_minimum_required(VERSION 3.0.0)
project(HyPhy)
cmake_policy(VERSION 3.0.0)
cmake_policy(SET CMP0026 OLD)

set(CMAKE_CONFIGURATION_TYPES Release)

#-------------------------------------------------------------------------------
# OPTIONS
#-------------------------------------------------------------------------------

option(NOAVX OFF)
option(NOSSE3 OFF)

#-------------------------------------------------------------------------------
# SSE MACROS
#-------------------------------------------------------------------------------

macro(PCL_CHECK_FOR_SSE3)
  include(CheckCXXSourceRuns)
  set(CMAKE_REQUIRED_FLAGS)

  if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
      set(CMAKE_REQUIRED_FLAGS "-msse3")
  endif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)

  check_cxx_source_runs("
      #include <pmmintrin.h>
      int main ()
      {
          __m128d a, b;
          double vals[2] = {0};
          a = _mm_loadu_pd (vals);
          b = _mm_hadd_pd (a,a);
          _mm_storeu_pd (vals, b);
          return (0);
      }"
      HAVE_SSE3_EXTENSIONS)
endmacro(PCL_CHECK_FOR_SSE3)

#-------------------------------------------------------------------------------
# SSE MACROS
#-------------------------------------------------------------------------------

macro(PCL_CHECK_FOR_AVX)
  include(CheckCXXSourceRuns)
  set(CMAKE_REQUIRED_FLAGS)

  if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
      set(CMAKE_REQUIRED_FLAGS "-march=corei7-avx -mtune=corei7-avx")
  endif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)

  check_cxx_source_runs("
      #include <immintrin.h>
      int main ()
      {
          __m256d a, b;
          double vals[4] = {0.};
          a = _mm256_loadu_pd (vals);
          b = _mm256_hadd_pd (a,a);
          _mm256_storeu_pd (vals, b);
          return (0);
      }"
      HAVE_AVX_EXTENSIONS)
endmacro(PCL_CHECK_FOR_AVX)


#-------------------------------------------------------------------------------
# set the CMAKE_MODULE_PATH so we can import OpenCL
#-------------------------------------------------------------------------------
set(CMAKE_MODULE_PATH cmake)
set(HYPHY_VERSION 2.1)

#-------------------------------------------------------------------------------
# setup the files we'll be using
#-------------------------------------------------------------------------------
file(GLOB SRC_CORE src/core/*.cpp)
file(GLOB SRC_GTESTS tests/gtests/*.cpp)
file(GLOB SRC_NEW src/new/*.cpp)

set(SRC_LINK src/lib/link/THyPhy.cpp)
set(SRC_SQLITE3 contrib/SQLite-3.8.2/sqlite3.c)
set(SRC_UNIXMAIN src/mains/unix.cpp)
set(SRC_UTILS src/utils/hyphyunixutils.cpp)

set(SRC_COMMON ${SRC_CORE} ${SRC_NEW} ${SRC_SQLITE3} ${SRC_UTILS})


#-------------------------------------------------------------------------------
# make it fast, scotty
#-------------------------------------------------------------------------------

set(DEFAULT_WARNING_FLAGS "-w")
set(DEFAULT_DEBUG_WARNING_FLAGS "-Wall -Wno-int-to-pointer-cast -Wno-conversion-null -Wno-sign-compare -Wno-maybe-uninitialized")



if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
        OUTPUT_VARIABLE GCC_VERSION
    )

    if (${GCC_VERSION} VERSION_LESS 4.9)
        message(FATAL_ERROR "Require at least gcc-4.9")
    endif()

    if(${GCC_VERSION} VERSION_LESS 6.0)
      set(DEFAULT_COMPILE_FLAGS "-fsigned-char -O3 -std=gnu++14")
    else(${GCC_VERSION} VERSION_LESS 6.0)
      set(DEFAULT_COMPILE_FLAGS "-fsigned-char -O3 -std=c++14")
    endif(${GCC_VERSION} VERSION_LESS 6.0)

    if(NOAVX)
        set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -mno-avx")
        if(NOSSE3)
            set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -mno-sse3")
        else(NOSSE)
            PCL_CHECK_FOR_SSE3()
            if(${HAVE_SSE3_EXTENSIONS})
                add_definitions (-D_SLKP_USE_SSE_INTRINSICS)
                set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -msse3 ")
            endif(${HAVE_SSE3_EXTENSIONS})
        endif(NOSSE3)
    else(NOAVX)
       PCL_CHECK_FOR_AVX()
       if(${HAVE_AVX_EXTENSIONS})
            set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -march=corei7-avx -mtune=corei7-avx")
            add_definitions (-D_SLKP_USE_AVX_INTRINSICS)
        else(${HAVE_AVX_EXTENSIONS})
            PCL_CHECK_FOR_SSE3()
            if(${HAVE_SSE3_EXTENSIONS})
                add_definitions (-D_SLKP_USE_SSE_INTRINSICS)
                set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -msse3 ")
            endif(${HAVE_SSE3_EXTENSIONS})
        endif(${HAVE_AVX_EXTENSIONS})
    endif(NOAVX)


    set(DEFAULT_LINK_FLAGS "${DEFAULT_COMPILE_FLAGS}")
endif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(DEFAULT_COMPILE_FLAGS "-fsigned-char -O3")
    PCL_CHECK_FOR_AVX()
    if(${HAVE_AVX_EXTENSIONS})
        set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -march=corei7-avx -mtune=corei7-avx")
        add_definitions (-D_SLKP_USE_AVX_INTRINSICS)
    else(${HAVE_AVX_EXTENSIONS})
        PCL_CHECK_FOR_SSE3()
        if(${HAVE_SSE3_EXTENSIONS})
            add_definitions (-D_SLKP_USE_SSE_INTRINSICS)
            set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -msse3 ")
        endif(${HAVE_SSE3_EXTENSIONS})
    endif (${HAVE_AVX_EXTENSIONS})

    set_property(
     SOURCE ${SRC_CORE}
     APPEND_STRING PROPERTY  COMPILE_FLAGS " -std=c++14 "
    )
endif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")

set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS}")
MESSAGE ("Set compiler flags to ${DEFAULT_COMPILE_FLAGS}")

if(NOT DEFINED DEFAULT_COMPILE_FLAGS)
    set(DEFAULT_COMPILE_FLAGS "")
endif(NOT DEFINED DEFAULT_COMPILE_FLAGS)

if(NOT DEFINED DEFAULT_LINK_FLAGS)
    set(DEFAULT_LINK_FLAGS "")
endif(NOT DEFINED DEFAULT_LINK_FLAGS)

if(NOT DEFINED DEFAULT_WARNING_FLAGS)
    set(DEFAULT_WARNING_FLAGS "")
endif(NOT DEFINED DEFAULT_WARNING_FLAGS)


#-------------------------------------------------------------------------------
# OpenMP support
#-------------------------------------------------------------------------------
find_package(OpenMP)

if(NOT ${OPENMP_FOUND})
    set(OpenMP_CXX_FLAGS "")
endif(NOT ${OPENMP_FOUND})

#-------------------------------------------------------------------------------
# default installation prefix
#-------------------------------------------------------------------------------
set(INSTALL_PREFIX /usr/local CACHE PATH "Installation prefix")
set(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX} CACHE INTERNAL "Installation prefix" FORCE)
set(DEFAULT_LIBRARIES dl pthread)

#-------------------------------------------------------------------------------
# LIBCurl support
#-------------------------------------------------------------------------------
find_package(CURL)
if(${CURL_FOUND} AND NOT APPLE)
    set(DEFAULT_LIBRARIES dl pthread crypto curl ssl)
    add_definitions (-D__HYPHYCURL__)
endif(${CURL_FOUND} AND NOT APPLE)

#-------------------------------------------------------------------------------
# gtest dependency
#-------------------------------------------------------------------------------
add_subdirectory(contrib/gtest-1.7.0)
set_property(
    TARGET gtest gtest_main
    PROPERTY EXCLUDE_FROM_ALL true
)


#-------------------------------------------------------------------------------
# uninstall target
#-------------------------------------------------------------------------------
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY
)
add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
)




#-------------------------------------------------------------------------------
# setup some source properties so that we don't spam stdout
#-------------------------------------------------------------------------------
set_property(
    SOURCE ${SRC_CORE}
    APPEND_STRING PROPERTY  COMPILE_FLAGS "${DEFAULT_WARNING_FLAGS}"
)

set_property(
    SOURCE ${SRC_SQLITE3}
    APPEND_STRING PROPERTY COMPILE_FLAGS " -Wno-int-to-pointer-cast -Wno-pointer-to-int-cast -std=c99"
)

set_property(
    SOURCE src/core/calcnode2.cpp
    APPEND_STRING PROPERTY COMPILE_FLAGS " -fstrict-aliasing -funroll-loops"
)

set_property(
  SOURCE ${SRC_CORE}
  APPEND_STRING PROPERTY COMPILE_FLAGS " -Wno-int-to-pointer-cast"
)

#-------------------------------------------------------------------------------
# shared compilation definitions and header includes
#-------------------------------------------------------------------------------
add_definitions(-D_SLKP_LFENGINE_REWRITE_ -D__AFYP_REWRITE_BGM__)

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    add_definitions(-D__HYPHY_64__)
endif(CMAKE_SIZEOF_VOID_P EQUAL 8)

include_directories(
    src/core/include
    src/contrib
    src/lib/Link src/new/include
    contrib/SQLite-3.8.2 # SQLite
    contrib/gtest-1.7.0/include
    tests/gtests
)


#-------------------------------------------------------------------------------
# shared hyphy hbl stdlib target
#-------------------------------------------------------------------------------
install(DIRECTORY res/ DESTINATION lib/hyphy)

#-------------------------------------------------------------------------------
# HYPHYMP target
#-------------------------------------------------------------------------------
add_executable(
    HYPHYMP
    EXCLUDE_FROM_ALL
    ${SRC_COMMON} ${SRC_UNIXMAIN}
)

target_link_libraries(HYPHYMP ${DEFAULT_LIBRARIES})

add_custom_target(MP DEPENDS HYPHYMP)

install(
    TARGETS HYPHYMP
    RUNTIME DESTINATION bin
    OPTIONAL
)

#-------------------------------------------------------------------------------
# hyphy OpenCL target
#-------------------------------------------------------------------------------
find_package(OpenCL)
if(${OPENCL_FOUND})

    add_executable(
        HYPHYOCL
        EXCLUDE_FROM_ALL
        ${SRC_COMMON} ${SRC_UNIXMAIN} ${CMAKE_CURRENT_BINARY_DIR}/src/ocl/opencl_kernels.h
    )

    include_directories(${OPENCL_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/src/ocl)
    target_link_libraries(HYPHYOCL ${DEFAULT_LIBRARIES} ${OPENCL_LIBRARIES})

    install(
        TARGETS HYPHYOCL
        RUNTIME DESTINATION bin
        OPTIONAL
    )

    add_custom_target(OCL DEPENDS HYPHYOCL)
    set_property(
        TARGET HYPHYOCL
        APPEND PROPERTY COMPILE_DEFINITIONS __MP__ __MP2__ __UNIX__ _HY_GPU_EXAMPLE_CALCULATOR MDSOCL
    )
    set_target_properties(
        HYPHYOCL
        PROPERTIES
        COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${DEFAULT_WARNING_FLAGS}"
        LINK_FLAGS "${DEFAULT_LINK_FLAGS}"
    )
    add_executable(
        cl2hdr
        EXCLUDE_FROM_ALL
        contrib/cl2hdr.c
    )
    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/src/ocl/opencl_kernels.h
        COMMAND cl2hdr ${CMAKE_CURRENT_SOURCE_DIR}/src/ocl/opencl_kernels.cl ${CMAKE_CURRENT_BINARY_DIR}/src/ocl/opencl_kernels.h
        DEPENDS cl2hdr ${CMAKE_CURRENT_SOURCE_DIR}/src/ocl/opencl_kernels.cl
    )
endif(${OPENCL_FOUND})


#-------------------------------------------------------------------------------
# hyphy mpi target
#-------------------------------------------------------------------------------
find_package(MPI)

if(${MPI_FOUND})

    add_executable(
        HYPHYMPI
        EXCLUDE_FROM_ALL
        ${SRC_COMMON} ${SRC_UNIXMAIN}
    )

    include_directories(${MPI_INCLUDE_PATH})

    target_link_libraries(HYPHYMPI ${DEFAULT_LIBRARIES} ${MPI_LIBRARIES})

    install(
        TARGETS HYPHYMPI
        RUNTIME DESTINATION bin
        OPTIONAL
    )

    add_custom_target(MPI DEPENDS HYPHYMPI)

    # set definitions here 'cause "s_t_p" below doesn't like lists
    set_property(
        TARGET HYPHYMPI
        APPEND PROPERTY COMPILE_DEFINITIONS __HYPHYMPI__ __UNIX__ _HYPHY_LIBDIRECTORY_="${CMAKE_INSTALL_PREFIX}/lib/hyphy"
    )

    set_target_properties(
        HYPHYMPI
        PROPERTIES
        COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${DEFAULT_WARNING_FLAGS} ${MPI_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS}"
        LINK_FLAGS "${DEFAULT_LINK_FLAGS} ${MPI_LINK_FLAGS} ${OpenMP_CXX_FLAGS}"
    )
endif(${MPI_FOUND})


#-------------------------------------------------------------------------------
# hyphy debug target
#-------------------------------------------------------------------------------
add_executable(
    HYPHYDEBUG
    EXCLUDE_FROM_ALL
    ${SRC_COMMON} ${SRC_UNIXMAIN}
)

target_link_libraries(HYPHYDEBUG ${DEFAULT_LIBRARIES})
add_custom_target(DEBUG DEPENDS HYPHYDEBUG)



set_target_properties(
    HYPHYDEBUG
    PROPERTIES
    COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${DEFAULT_DEBUG_WARNING_FLAGS} ${OpenMP_CXX_FLAGS} -g -fprofile-arcs -ftest-coverage -pg -O0 "
    LINK_FLAGS "${DEFAULT_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS} -g -fprofile-arcs -ftest-coverage -pg -O0 "
)



#-------------------------------------------------------------------------------
# hyphy gtest target
#-------------------------------------------------------------------------------
add_executable(
    HYPHYGTEST
    EXCLUDE_FROM_ALL
    ${SRC_COMMON} ${SRC_GTESTS} ${SRC_UNIXMAIN}
)
add_dependencies(HYPHYGTEST gtest)
target_link_libraries(HYPHYGTEST ${DEFAULT_LIBRARIES} gtest)
add_custom_target(GTEST DEPENDS HYPHYGTEST)

set_property(
    TARGET HYPHYGTEST
    APPEND PROPERTY COMPILE_DEFINITIONS __UNITTEST__
)
set_target_properties(
    HYPHYGTEST
    PROPERTIES
    COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${DEFAULT_DEBUG_WARNING_FLAGS} ${OpenMP_CXX_FLAGS} -msse3 -g -fprofile-arcs -ftest-coverage -pg"
    LINK_FLAGS "${DEFAULT_COMPILE_FLAGS} ${DEFAULT_DEBUG_WARNING_FLAGS} ${OpenMP_CXX_FLAGS} -msse3 -g -fprofile-arcs -ftest-coverage -pg"
)


#-------------------------------------------------------------------------------
# target-specific and occasionally-shared properties
#-------------------------------------------------------------------------------

if(UNIX)
    set_property(
        TARGET HYPHYMP HYPHYGTEST HYPHYDEBUG
        APPEND PROPERTY COMPILE_DEFINITIONS __UNIX__
    )
endif(UNIX)

set_property(
    TARGET HYPHYMP
    APPEND PROPERTY COMPILE_DEFINITIONS __MP__ __MP2__
)

set_property(
    TARGET HYPHYMP HYPHYGTEST HYPHYDEBUG
    APPEND PROPERTY COMPILE_DEFINITIONS _HYPHY_LIBDIRECTORY_="${CMAKE_INSTALL_PREFIX}/lib/hyphy"
)

set_property(
    TARGET HYPHYDEBUG HYPHYGTEST
    APPEND PROPERTY COMPILE_DEFINITIONS __HYPHYDEBUG__
)

set_target_properties(
    HYPHYMP
    PROPERTIES
    COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${DEFAULT_WARNING_FLAGS} ${OpenMP_CXX_FLAGS}"
    LINK_FLAGS "${DEFAULT_LINK_FLAGS} ${OpenMP_CXX_FLAGS}"
)

