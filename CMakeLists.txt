cmake_minimum_required(VERSION 2.8.1)
project(HyPhy)
cmake_policy(SET CMP0042 NEW)

set(CMAKE_CONFIGURATION_TYPES Release)

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
# figure out some system-stuff for compile-time workarounds
#-------------------------------------------------------------------------------


set(MACOSX_LION false)
set(MACOSX false)

if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    set(MACOSX true)
    if(${CMAKE_SYSTEM_VERSION} VERSION_GREATER 11.0.0 OR ${CMAKE_SYSTEM_VERSION} VERSION_EQUAL 11.0.0)
        set(MACOSX_LION true)
    endif(${CMAKE_SYSTEM_VERSION} VERSION_GREATER 11.0.0 OR ${CMAKE_SYSTEM_VERSION} VERSION_EQUAL 11.0.0)
    
    if(${CMAKE_SYSTEM_VERSION} VERSION_GREATER 12.0.0 OR ${CMAKE_SYSTEM_VERSION} VERSION_EQUAL 12.0.0)
        set(MACOSX_LION false)
    endif(${CMAKE_SYSTEM_VERSION} VERSION_GREATER 12.0.0 OR ${CMAKE_SYSTEM_VERSION} VERSION_EQUAL 12.0.0)
    
endif(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")

set (SSE_RELATED_FLAGS "")

#-------------------------------------------------------------------------------
# make it fast, scotty
#-------------------------------------------------------------------------------
if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
        OUTPUT_VARIABLE GCC_VERSION
    )

    set(GCC45 false)
    set(GCC46 false)

    if(${GCC_VERSION} VERSION_GREATER 4.5 OR ${GCC_VERSION} VERSION_EQUAL 4.5)
        set(GCC45 true)
    endif(${GCC_VERSION} VERSION_GREATER 4.5 OR ${GCC_VERSION} VERSION_EQUAL 4.5)

    if(${GCC_VERSION} VERSION_GREATER 4.6 OR ${GCC_VERSION} VERSION_EQUAL 4.6)
        set(GCC46 true)
    endif(${GCC_VERSION} VERSION_GREATER 4.6 OR ${GCC_VERSION} VERSION_EQUAL 4.6)

    #if(${MACOSX})
    #    set(DEFAULT_WARNING_FLAGS "-Wno-int-to-pointer-cast -Wno-conversion-null -Wno-dangling-else -Wno-logical-op-parentheses")
    #endif(${MACOSX})

    if(${GCC45})
        set(DEFAULT_COMPILE_FLAGS "-fsigned-char -O3 ")
    else(${GCC45})
        set(DEFAULT_COMPILE_FLAGS "-fsigned-char -O3  ")
    endif(${GCC45})

    if(${GCC46})
        set(DEFAULT_WARNING_FLAGS "-Wno-int-to-pointer-cast -Wno-conversion-null")
    endif(${GCC46})
    
    PCL_CHECK_FOR_AVX()
    if(${HAVE_AVX_EXTENSIONS})
        set(SSE_RELATED_FLAGS " -march=corei7-avx -mtune=corei7-avx ")
        add_definitions (-D_SLKP_USE_AVX_INTRINSICS)
    else(${HAVE_AVX_EXTENSIONS})
        PCL_CHECK_FOR_SSE3()
        if(${HAVE_SSE3_EXTENSIONS})
            add_definitions (-D_SLKP_USE_SSE_INTRINSICS)
            set(SSE_RELATED_FLAGS " -msse3 ")
        endif(${HAVE_SSE3_EXTENSIONS})
    endif (${HAVE_AVX_EXTENSIONS})
    
    set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${SSE_RELATED_FLAGS} ")    
    set(DEFAULT_LINK_FLAGS "${DEFAULT_COMPILE_FLAGS}")

endif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)

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

if(${MACOSX_LION} OR NOT ${OPENMP_FOUND})
    set(OpenMP_CXX_FLAGS "")
endif(${MACOSX_LION} OR NOT ${OPENMP_FOUND})


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
if(${CURL_FOUND})
    set(DEFAULT_LIBRARIES dl pthread crypto curl ssl)
    add_definitions (-D__HYPHYCURL__)
endif(${CURL_FOUND})

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
# setup the files we'll be using
#-------------------------------------------------------------------------------
file(GLOB SRC_CORE src/core/*.cpp)
file(GLOB SRC_GUI src/gui/*.cpp src/gui/Components/*.cpp)
file(GLOB SRC_GTESTS tests/gtests/core/*.cpp)
file(GLOB SRC_NEW src/new/*.cpp)
file(GLOB SRC_GRAMMAR src/grammar/*.cpp)

set(SRC_LINK src/lib/Link/THyPhy.cpp)
set(SRC_PREFS src/gui/preferences.cpp)
set(SRC_SQLITE3 contrib/SQLite-3.8.2/sqlite3.c)
set(SRC_UNIXMAIN src/mains/unix.cpp)
set(SRC_GTKMAIN src/mains/gtk.cpp)
set(SRC_MACMAIN src/mains/mac.cpp)
set(SRC_UTILS src/utils/hyphyunixutils.cpp)

set(SRC_COMMON ${SRC_CORE} ${SRC_NEW} ${SRC_SQLITE3} ${SRC_UTILS} ${SRC_GRAMMAR})
set (SRC_WORKING src/core/baseobj.cpp src/core/site.cpp src/core/strings.cpp src/core/error.cpp src/core/string_buffer.cpp src/core/helperfunctions.cpp src/core/regex.cpp src/core/globals.cpp src/core/trie.cpp)

#-------------------------------------------------------------------------------
# setup some source properties so that we don't spam stdout
#-------------------------------------------------------------------------------
set_property(
    SOURCE ${SRC_CORE}
    PROPERTY COMPILE_FLAGS "${DEFAULT_WARNING_FLAGS}"
)

set_property(
    SOURCE ${SRC_SQLITE3}
    PROPERTY COMPILE_FLAGS "-Wno-int-to-pointer-cast -Wno-pointer-to-int-cast"
)

set_property(
    SOURCE src/core/calcnode2.cpp
    PROPERTY COMPILE_FLAGS "-fstrict-aliasing -funroll-loops"
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
    src/core/templates
    src/lib/Link 
    src/new/include
    contrib/SQLite-3.7.11 # SQLite
    contrib/tinyxml
    contrib/pugixml/src/
    src/ocl/
    contrib/gtest-1.7.0/include
    contrib/gtest-1.7.0/gtest/include
    src/grammar
    tests/gtests
)


#-------------------------------------------------------------------------------
# shared hyphy hbl stdlib target
#-------------------------------------------------------------------------------
install(DIRECTORY res/ DESTINATION lib/hyphy)


#-------------------------------------------------------------------------------
# hyphy library target
#-------------------------------------------------------------------------------
add_library(
  hyphy_mp SHARED
  EXCLUDE_FROM_ALL
  ${SRC_COMMON} ${SRC_LINK} ${SRC_PREFS}
)
target_link_libraries(hyphy_mp ${DEFAULT_LIBRARIES})

#install(
#  TARGETS hyphy_mp
#  LIBRARY DESTINATION lib
#  OPTIONAL
#)

add_custom_target(LIB DEPENDS hyphy_mp)

set_property(
    TARGET hyphy_mp
    APPEND PROPERTY COMPILE_DEFINITIONS __HEADLESS__
)


#-------------------------------------------------------------------------------
# hyphy mp2 target
#-------------------------------------------------------------------------------
add_executable(
    HYPHYMP
    EXCLUDE_FROM_ALL
    ${SRC_COMMON} ${SRC_UNIXMAIN}
)

target_link_libraries(HYPHYMP ${DEFAULT_LIBRARIES})

#install(
#    TARGETS HYPHYMP
#    RUNTIME DESTINATION bin
#    OPTIONAL
#)

add_custom_target(MP2 DEPENDS HYPHYMP)


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

    #install(
    #    TARGETS HYPHYOCL
    #    RUNTIME DESTINATION bin
    #    OPTIONAL
    #)

    add_custom_target(OCL DEPENDS HYPHYOCL)

    set_property(
        TARGET HYPHYOCL
        APPEND PROPERTY COMPILE_DEFINITIONS __MP__ __MP2__ __UNIX__ _HY_GPU_EXAMPLE_CALCULATOR MDSOCL
    )

    set_target_properties(
        HYPHYOCL
        PROPERTIES
        COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${SSE_RELATED_FLAGS}"
        LINK_FLAGS "${DEFAULT_LINK_FLAGS} ${SSE_RELATED_FLAGS}"
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

    #install(
    #    TARGETS HYPHYMPI
    #    RUNTIME DESTINATION bin
    #    OPTIONAL
    #)

    add_custom_target(MPI DEPENDS HYPHYMPI)

    # set definitions here 'cause "s_t_p" below doesn't like lists
    set_property(
        TARGET HYPHYMPI
        APPEND PROPERTY COMPILE_DEFINITIONS __HYPHYMPI__ __UNIX__ _HYPHY_LIBDIRECTORY_="${CMAKE_INSTALL_PREFIX}/lib/hyphy"
    )

    set_target_properties(
        HYPHYMPI
        PROPERTIES
        COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS} ${SSE_RELATED_FLAGS}"
        LINK_FLAGS "${DEFAULT_LINK_FLAGS} ${MPI_LINK_FLAGS} ${SSE_RELATED_FLAGS}"
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
    COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${SSE_RELATED_FLAGS} ${OpenMP_CXX_FLAGS} -g -O0"
    LINK_FLAGS "${DEFAULT_COMPILE_FLAGS} ${SSE_RELATED_FLAGS} ${OpenMP_CXX_FLAGS} -g -O0"
)

set_property(
    TARGET HYPHYDEBUG
    APPEND PROPERTY COMPILE_DEFINITIONS __NEW_GRAMMAR__
)

#install(
#  TARGETS HYPHYDEBUG
#        RUNTIME DESTINATION bin
#        OPTIONAL
#    )



#-------------------------------------------------------------------------------
# hyphy gtest target
#-------------------------------------------------------------------------------
add_executable(
    HYPHYGTEST
    EXCLUDE_FROM_ALL
    ${SRC_WORKING} ${SRC_GTESTS} tests/gtests/main.cpp
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
    COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${SSE_RELATED_FLAGS} ${OpenMP_CXX_FLAGS} -g -fprofile-arcs -ftest-coverage -O0 "
    LINK_FLAGS "${DEFAULT_COMPILE_FLAGS} ${SSE_RELATED_FLAGS} ${OpenMP_CXX_FLAGS} -g -fprofile-arcs -ftest-coverage -O0 "
)


#-------------------------------------------------------------------------------
# target-specific and occasionally-shared properties
#-------------------------------------------------------------------------------
if(UNIX)
    set_property(
        TARGET HYPHYMP hyphy_mp HYPHYGTEST HYPHYDEBUG
        APPEND PROPERTY COMPILE_DEFINITIONS __UNIX__
    )
endif(UNIX)

set_property(
    TARGET hyphy_mp HYPHYMP
    APPEND PROPERTY COMPILE_DEFINITIONS __MP__ __MP2__
)

set_property(
    TARGET hyphy_mp HYPHYMP HYPHYGTEST HYPHYDEBUG
    APPEND PROPERTY COMPILE_DEFINITIONS _HYPHY_LIBDIRECTORY_="${CMAKE_INSTALL_PREFIX}/lib/hyphy"
)

set_property(
    TARGET HYPHYDEBUG HYPHYGTEST
    APPEND PROPERTY COMPILE_DEFINITIONS __HYPHYDEBUG__
)

set_target_properties(
    hyphy_mp HYPHYMP
    PROPERTIES
    COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS}"
    LINK_FLAGS "${DEFAULT_LINK_FLAGS} ${OpenMP_CXX_FLAGS}"
)

