cmake_minimum_required(VERSION 3.12)
project(HyPhy)
#cmake_policy(VERSION 3.0.0)

enable_testing()
#enable_language(Fortran)

set(CMAKE_CONFIGURATION_TYPES Release)

#-------------------------------------------------------------------------------
# OPTIONS
#-------------------------------------------------------------------------------

option(NOAVX OFF)
option(NOSSE4 OFF)
option(NONEON OFF)
option(NOZLIB OFF)
option(NOBLAS OFF)

if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
   set(NOAVX ON)
   #set(NOSSE4 ON)
   set(NONEON ON)
   set(NOZLIB ON)
   set(LIBRARY_PATH "/hyphy")
   #set(CMAKE_EXECUTABLE_SUFFIX ".html")
   add_definitions (-D_USE_EMSCRIPTEN_)
 endif()

#-------------------------------------------------------------------------------
# SSE MACROS
#-------------------------------------------------------------------------------

macro(PCL_CHECK_FOR_SSE4)
  include(CheckCXXSourceRuns)
  set(CMAKE_REQUIRED_FLAGS)

  #if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
      set(CMAKE_REQUIRED_FLAGS "-msse4.1")
  #endif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)

  check_cxx_source_runs("
      #include <smmintrin.h>
      #include <stdio.h>

      int main ()
      {
          __m128d a, b;
          double vals[2] = {0};
          a = _mm_loadu_pd (vals);
          b = _mm_hadd_pd (a,a);
          _mm_storeu_pd (vals, _mm_blend_pd (a,b,2));
          printf (\"%g\", vals[0]);
          return (0);
      }"
      HAVE_SSE4_EXTENSIONS)
endmacro(PCL_CHECK_FOR_SSE4)

#-------------------------------------------------------------------------------
# NEON MACROS
#-------------------------------------------------------------------------------

macro(PCL_CHECK_FOR_NEON)
	  include(CheckCXXSourceRuns)
	  set(CMAKE_REQUIRED_FLAGS)

	  check_cxx_source_runs("
		#include <stdio.h>
		#include <arm_neon.h>
		int main() {
			double array[4] = {4.,2.,3.,4};
			float64x2_t a1 = vld1q_f64 (array),
						a2 = vld1q_f64 (array+2);
			a1 = vmulq_f64  (a1,a2);
			vst1q_lane_f64 (array, a1, 0);
			printf (\"%g\", array[0]);
			return 0;
	   }"
       HAVE_NEON_EXTENSIONS)
endmacro(PCL_CHECK_FOR_NEON)


#-------------------------------------------------------------------------------
# SSE MACROS
#-------------------------------------------------------------------------------

macro(PCL_CHECK_FOR_AVX)
  include(CheckCXXSourceRuns)
  set(CMAKE_REQUIRED_FLAGS)

  #if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
      set(CMAKE_REQUIRED_FLAGS "-mavx2")
  #endif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)

  check_cxx_source_runs("
      #include <immintrin.h>
      #include <stdio.h>

      int main ()
      {
          __m256d a, b;
          double vals[4] = {0.};
          a = _mm256_loadu2_m128d (vals, vals+2);
          b = _mm256_permute4x64_pd (_mm256_hadd_pd (a,a), 2);
          _mm256_storeu_pd (vals, b);
          printf (\"%g\", vals[0]);
          return (0);
      }"
      HAVE_AVX_EXTENSIONS)
endmacro(PCL_CHECK_FOR_AVX)

#-------------------------------------------------------------------------------
# FMA3 MACROS
#-------------------------------------------------------------------------------

macro(PCL_CHECK_FOR_FMA3)
    include(CheckCXXSourceRuns)
    set(CMAKE_REQUIRED_FLAGS)

    #if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
        set(CMAKE_REQUIRED_FLAGS "-mavx2 -mfma ")
    #endif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)

    check_cxx_source_runs("
    #include <immintrin.h>
    #include <stdio.h>

    int main ()
    {
        __m256d a, b, c;
        double vals[4] = {2.};
        a = _mm256_loadu_pd (vals);
        b = _mm256_loadu_pd (vals);
        c = _mm256_loadu_pd (vals);
        a = _mm256_fmadd_pd      (a, b, c);
        _mm256_storeu_pd (vals, a);
        printf (\"%g\", vals[0]);
        return 0;
    }"
    HAVE_FMA3)
endmacro(PCL_CHECK_FOR_FMA3)


#-------------------------------------------------------------------------------
# set the CMAKE_MODULE_PATH so we can import OpenCL
#-------------------------------------------------------------------------------
set(CMAKE_MODULE_PATH cmake)
set(HYPHY_VERSION 2.1)

#-------------------------------------------------------------------------------
# setup the files we'll be using
#-------------------------------------------------------------------------------
file(GLOB SRC_CORE src/core/*.cpp)
file(GLOB SRC_NEW src/new/*.cpp)
file(GLOB SRC_CONTRIB src/contrib/*.cpp)

set(SRC_LINK src/lib/link/THyPhy.cpp)
set(SRC_SQLITE3 contrib/SQLite-3.8.2/sqlite3.c)
set(SRC_UNIXMAIN src/mains/unix.cpp)
set(SRC_UTILS src/utils/hyphyunixutils.cpp)


set(SRC_COMMON ${SRC_CORE} ${SRC_NEW} ${SRC_SQLITE3} ${SRC_UTILS} ${SRC_CONTRIB})
set_source_files_properties(${SRC_CORE} ${SRC_NEW} {SRC_UTILS} PROPERTIES COMPILE_FLAGS "-Weffc++ -Wextra -Wall")


#-------------------------------------------------------------------------------
# make it fast, scotty
#-------------------------------------------------------------------------------

set(DEFAULT_WARNING_FLAGS " -w -Weffc++ -Wextra -Wall ")
set(DEFAULT_DEBUG_WARNING_FLAGS "-Wall -Wno-int-to-pointer-cast -Wno-conversion-null -Wno-sign-compare")
set(ADDITIONAL_FLAGS "")

if (DEBUGFLAGS)
	set(ADDITIONAL_FLAGS "-g")
endif()


if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
     execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
        OUTPUT_VARIABLE GCC_VERSION
    )

    if (${GCC_VERSION} VERSION_LESS 4.9)
        message(FATAL_ERROR "Require at least gcc-4.9")
    endif()

    if(${GCC_VERSION} VERSION_LESS 6.0)
      set(DEFAULT_COMPILE_FLAGS "-fsigned-char -O3 -D_FORTIFY_SOURCE=2   -Wall -std=gnu++14")
    else(${GCC_VERSION} VERSION_LESS 6.0)
      set(DEFAULT_COMPILE_FLAGS "-fsigned-char -O3 -D_FORTIFY_SOURCE=2   -Wall -std=c++14")
    endif(${GCC_VERSION} VERSION_LESS 6.0)

    #set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS}  -fopt-info")

    if(NOAVX)
        set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -mno-avx")
        if(NOSSE4)
	        if(NOT CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
            	set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -mno-sse3")
            endif()
        else(NOSSE4)
            PCL_CHECK_FOR_SSE4()
            if(${HAVE_SSE4_EXTENSIONS})
                add_definitions (-D_SLKP_USE_SSE_INTRINSICS)
                set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -msse4.1 ")
                if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
				 	set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -msimd128 ")
			    endif()
            endif(${HAVE_SSE4_EXTENSIONS})
        endif(NOSSE4)
    else(NOAVX)
       PCL_CHECK_FOR_AVX()
       if(${HAVE_AVX_EXTENSIONS})
            set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -march=native -mtune=native -mavx2")
            add_definitions (-D_SLKP_USE_AVX_INTRINSICS)
            PCL_CHECK_FOR_FMA3()
            if (${HAVE_FMA3})
                set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS}  -mfma")
                add_definitions (-D_SLKP_USE_FMA3_INTRINSICS)
            endif (${HAVE_FMA3})
        else(${HAVE_AVX_EXTENSIONS})
            PCL_CHECK_FOR_SSE4()
            if(${HAVE_SSE4_EXTENSIONS})
                add_definitions (-D_SLKP_USE_SSE_INTRINSICS)
                set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -msse4.1 ")
                

            endif(${HAVE_SSE4_EXTENSIONS})
        endif(${HAVE_AVX_EXTENSIONS})
    endif(NOAVX)
    	
	if(NOT NONEON)
		PCL_CHECK_FOR_NEON()
		if(${HAVE_NEON_EXTENSIONS})
			add_definitions (-D_SLKP_USE_ARM_NEON)

		endif(${HAVE_NEON_EXTENSIONS})
	endif(NOT NONEON)

    set(DEFAULT_LINK_FLAGS "${DEFAULT_COMPILE_FLAGS}")
endif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    #set(DEFAULT_COMPILE_FLAGS "-fsigned-char -O3  -fsave-optimization-record=yaml ")
    set(DEFAULT_COMPILE_FLAGS "-fsigned-char -O3 -D_FORTIFY_SOURCE=2  -std=c++14 -Wall")
    if(NOAVX)
        if(NOSSE4)
	        if(NOT CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
            	set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -mno-sse3")
            endif()
         else(NOSSE4)
            PCL_CHECK_FOR_SSE4()
 			if("${HAVE_SSE4_EXTENSIONS}" OR ${CMAKE_SYSTEM_NAME} STREQUAL "Emscripten")
				add_definitions (-D_SLKP_USE_SSE_INTRINSICS)
				set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -msse4.1 ")
				if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
					set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -msimd128 ")
				endif()
			endif()  
        endif(NOSSE4)
    else(NOAVX)
        PCL_CHECK_FOR_AVX()
        if(${HAVE_AVX_EXTENSIONS})
            set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -march=native -mtune=native -mavx2")
            add_definitions (-D_SLKP_USE_AVX_INTRINSICS)
            PCL_CHECK_FOR_FMA3()
            if (${HAVE_FMA3})
                set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -mfma")
                add_definitions (-D_SLKP_USE_FMA3_INTRINSICS)
            endif (${HAVE_FMA3})
        else(${HAVE_AVX_EXTENSIONS})
            PCL_CHECK_FOR_SSE4()
            if(${HAVE_SSE4_EXTENSIONS})
                add_definitions (-D_SLKP_USE_SSE_INTRINSICS)
                set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -msse4.1 ")
            endif(${HAVE_SSE4_EXTENSIONS})
        endif (${HAVE_AVX_EXTENSIONS})
    endif(NOAVX)
    
	if(NOT NONEON)
		PCL_CHECK_FOR_NEON()
		if(${HAVE_NEON_EXTENSIONS})
			add_definitions (-D_SLKP_USE_ARM_NEON)
		endif(${HAVE_NEON_EXTENSIONS})
	endif(NOT NONEON)

    set_property(
     SOURCE ${SRC_CORE} ${SRC_NEW} ${SRC_UTILS} ${SRC_UNIXMAIN}
     APPEND_STRING PROPERTY  COMPILE_FLAGS " -std=c++14 "
    )
endif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")

#set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS}  -fopt-info -fopt-info-vec-missed")
set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${ADDITIONAL_FLAGS}")

if(NOT DEFINED DEFAULT_COMPILE_FLAGS)
    set(DEFAULT_COMPILE_FLAGS "")
endif(NOT DEFINED DEFAULT_COMPILE_FLAGS)

if(NOT DEFINED DEFAULT_LINK_FLAGS)
    set(DEFAULT_LINK_FLAGS "")
endif(NOT DEFINED DEFAULT_LINK_FLAGS)

set(DEFAULT_LINK_FLAGS "${DEFAULT_LINK_FLAGS}  ${ADDITIONAL_FLAGS}")

if(NOT DEFINED DEFAULT_WARNING_FLAGS)
    set(DEFAULT_WARNING_FLAGS "")
endif(NOT DEFINED DEFAULT_WARNING_FLAGS)

if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
   set(DEFAULT_COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -fwasm-exceptions ")
endif()

MESSAGE ("Set compiler flags to ${DEFAULT_COMPILE_FLAGS}")
MESSAGE ("Set linker flags to ${DEFAULT_LINK_FLAGS}")

set(BLA_VENDOR Apple)
find_package(BLAS)



#-------------------------------------------------------------------------------
# OpenMP support
#-------------------------------------------------------------------------------
find_package(OpenMP)

# emulate missing OpenMP for testing

#set (OPENMP_FOUND FALSE)
#set (OpenMP_CXX_FLAGS "")

#-------------------------------------------------------------------------------
# default installation prefix
#-------------------------------------------------------------------------------
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation prefix")
endif()
if(NOT DEFINED INSTALL_PREFIX)
    set(INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE INTERNAL "Installation prefix" FORCE)
endif()

if(CMAKE_SYSTEM_NAME MATCHES "FreeBSD")
	set(DEFAULT_LIBRARIES )
else(CMAKE_SYSTEM_NAME MATCHES "FreeBSD")
    set(DEFAULT_LIBRARIES dl)
endif(CMAKE_SYSTEM_NAME MATCHES "FreeBSD")

#-------------------------------------------------------------------------------
# LIBCurl support
#-------------------------------------------------------------------------------
find_package(CURL)
if(${CURL_FOUND})
    message(${CURL_LIBRARIES})
    list(APPEND DEFAULT_LIBRARIES ${CURL_LIBRARIES})
    add_definitions (-D__HYPHYCURL__)
endif(${CURL_FOUND})

if(NOT NOZLIB)
	find_package(ZLIB 1.2.9)
	if(${ZLIB_FOUND})
                list(APPEND DEFAULT_LIBRARIES ${ZLIB_LIBRARIES})
                include_directories(${ZLIB_INCLUDE_DIRS})
		add_definitions (-D__ZLIB__)
	endif(${ZLIB_FOUND})
endif(NOT NOZLIB)



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

find_program(NODE npm)
if(NOT NODE)
    message("Node not installed; API documentation will not be generated")
else(NOT NODE)
    add_custom_target(
        docs
        rm -rf ./docs && npm install -g documentation@4 && documentation build --polyglot -f md ./res/TemplateBatchFiles/libv3/**/*.bf -o ./docs.md
    )
endif(NOT NODE)


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


include_directories(
    src/core/include
    src/contrib
    src/lib/Link src/new/include
    contrib/SQLite-3.8.2 # SQLite
)

if (${BLAS_FOUND} AND NOT NOBLAS) 
	add_definitions (-D_SLKP_USE_APPLE_BLAS)
	set(DEFAULT_LIBRARIES "${DEFAULT_LIBRARIES}")
endif()

#-------------------------------------------------------------------------------
# shared hyphy hbl stdlib target
#-------------------------------------------------------------------------------
install(DIRECTORY res/ DESTINATION share/hyphy)

#-------------------------------------------------------------------------------
# HYPHYMP target
#-------------------------------------------------------------------------------
add_executable(
    hyphy
    ${SRC_COMMON} ${SRC_UNIXMAIN}
)

if(${OPENMP_FOUND})
    target_link_libraries(hyphy PRIVATE ${DEFAULT_LIBRARIES} OpenMP::OpenMP_CXX)
else(${OPENMP_FOUND})
    target_link_libraries(hyphy PRIVATE ${DEFAULT_LIBRARIES})
endif(${OPENMP_FOUND})

if(${BLAS_FOUND} AND NOT NOBLAS)
    target_link_libraries(hyphy PRIVATE ${BLAS_LIBRARIES})
endif ()

add_custom_target(MP DEPENDS hyphy)

install(
    TARGETS hyphy
    RUNTIME DESTINATION bin
    OPTIONAL
)

if (NOT DEFINED LIBRARY_PATH)
    set(LIBRARY_PATH "${CMAKE_INSTALL_PREFIX}/share/hyphy")
endif()

MESSAGE (${LIBRARY_PATH})

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

    if(${OPENMP_FOUND})
        target_link_libraries(HYPHYMPI ${DEFAULT_LIBRARIES} ${MPI_LIBRARIES} OpenMP::OpenMP_CXX)
    else(${OPENMP_FOUND})
        target_link_libraries(HYPHYMPI ${DEFAULT_LIBRARIES} ${MPI_LIBRARIES})
    endif(${OPENMP_FOUND})

	if(${BLAS_FOUND} AND NOT NOBLAS)
    	target_link_libraries(HYPHYMPI ${BLAS_LIBRARIES})
	endif ()

    install(
        TARGETS HYPHYMPI
        RUNTIME DESTINATION bin
        OPTIONAL
    )

    add_custom_target(MPI DEPENDS HYPHYMPI)

    # set definitions here 'cause "s_t_p" below doesn't like lists
    set_property(
        TARGET HYPHYMPI
        APPEND PROPERTY COMPILE_DEFINITIONS __HYPHYMPI__ __UNIX__ _HYPHY_LIBDIRECTORY_="${LIBRARY_PATH}"
    )

    message ("${DEFAULT_LINK_FLAGS} ${MPI_LINK_FLAGS} ${OpenMP_CXX_FLAGS}")

    set_target_properties(
        HYPHYMPI
        PROPERTIES
        if(${OPENMP_FOUND})
            COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${DEFAULT_WARNING_FLAGS} ${MPI_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS}"
            LINK_FLAGS "${DEFAULT_LINK_FLAGS} ${MPI_LINK_FLAGS} ${OpenMP_CXX_FLAGS}"
        else(${OPENMP_FOUND})
            COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${DEFAULT_WARNING_FLAGS} ${MPI_COMPILE_FLAGS}"
            LINK_FLAGS "${DEFAULT_LINK_FLAGS} ${MPI_LINK_FLAGS}"
        endif(${OPENMP_FOUND})
    )
endif(${MPI_FOUND})

#-------------------------------------------------------------------------------
# hyphy profile target
#-------------------------------------------------------------------------------
add_executable(
    HYPHYPROFILE
    EXCLUDE_FROM_ALL
    ${SRC_COMMON} ${SRC_UNIXMAIN}
)

if(${OPENMP_FOUND})
    target_link_libraries(HYPHYPROFILE ${DEFAULT_LIBRARIES} OpenMP::OpenMP_CXX)
else(${OPENMP_FOUND})
    target_link_libraries(HYPHYPROFILE ${DEFAULT_LIBRARIES})
endif(${OPENMP_FOUND})

add_custom_target(PROFILE DEPENDS HYPHYPROFILE)


set_target_properties(
    HYPHYPROFILE
    PROPERTIES
    if(${OPENMP_FOUND})
        COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS} -g -fprofile-generate -fanalyzer "
        LINK_FLAGS "${DEFAULT_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS} -g -fprofile-generate -fanalyzer "
    else(${OPENMP_FOUND})
        COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS}-g -fprofile-generate -fanalyzer "
        LINK_FLAGS "${DEFAULT_COMPILE_FLAGS} -g  -fprofile-generate -fanalyzer "
    endif(${OPENMP_FOUND})
)

#-------------------------------------------------------------------------------
 # hyphy profile target
 #-------------------------------------------------------------------------------
 add_executable(
     HYPHYGPERF
     EXCLUDE_FROM_ALL
     ${SRC_COMMON} ${SRC_UNIXMAIN}
 )

  if(${OPENMP_FOUND})
    target_link_libraries(HYPHYGPERF PRIVATE ${DEFAULT_LIBRARIES} OpenMP::OpenMP_CXX)
  else(${OPENMP_FOUND})
    target_link_libraries(HYPHYGPERF PRIVATE ${DEFAULT_LIBRARIES})
  endif(${OPENMP_FOUND})

  add_custom_target(GPERF DEPENDS HYPHYGPERF)

  set_target_properties(
     HYPHYGPERF
     PROPERTIES
         COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -g"
         LINK_FLAGS "${DEFAULT_LINK_FLAGS}  -g -lprofiler -ltcmalloc"
     endif(${OPENMP_FOUND})
 )


#-------------------------------------------------------------------------------
# hyphy debug target
#-------------------------------------------------------------------------------
add_executable(
    HYPHYDEBUG
    EXCLUDE_FROM_ALL
    ${SRC_COMMON} ${SRC_UNIXMAIN}
)

if(${OPENMP_FOUND})
    target_link_libraries(HYPHYDEBUG ${DEFAULT_LIBRARIES} OpenMP::OpenMP_CXX)
else(${OPENMP_FOUND})
    target_link_libraries(HYPHYDEBUG ${DEFAULT_LIBRARIES})
endif(${OPENMP_FOUND})

target_link_libraries(HYPHYDEBUG ${DEFAULT_LIBRARIES})

add_custom_target(DEBUG DEPENDS HYPHYDEBUG)

if(${BLAS_FOUND})
    target_link_libraries(HYPHYDEBUG ${BLAS_LIBRARIES})
endif ()

set_target_properties(
    HYPHYDEBUG
    PROPERTIES
    if(${OPENMP_FOUND})
        #COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${DEFAULT_DEBUG_WARNING_FLAGS} ${OpenMP_CXX_FLAGS} ${MPI_COMPILE_FLAGS} -g -fsanitize=address -O0 "
        COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${DEFAULT_DEBUG_WARNING_FLAGS} ${OpenMP_CXX_FLAGS} ${MPI_COMPILE_FLAGS} -g -fsanitize=address -O0 "
        LINK_FLAGS "${DEFAULT_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS} ${MPI_LINK_FLAGS}  -g -fsanitize=address -O0 "
    else(${OPENMP_FOUND})
        #COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${DEFAULT_DEBUG_WARNING_FLAGS} -g -fsanitize=address -O0 "
        COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${DEFAULT_DEBUG_WARNING_FLAGS} -g -fsanitize=address -O0 "
        LINK_FLAGS "${DEFAULT_COMPILE_FLAGS} -g -fsanitize=address -O0 "
    endif(${OPENMP_FOUND})
)



#-------------------------------------------------------------------------------
# target-specific and occasionally-shared properties
#-------------------------------------------------------------------------------

if(UNIX)
    set_property(
        TARGET hyphy HYPHYDEBUG HYPHYPROFILE HYPHYGPERF
        APPEND PROPERTY COMPILE_DEFINITIONS __UNIX__
    )
    if(${SWIG_FOUND})
    set_property(
        TARGET LIBRARY
        APPEND PROPERTY COMPILE_DEFINITIONS __UNIX__ __HEADLESS__
    )
    endif(${SWIG_FOUND})
endif(UNIX)

set_property(
    TARGET hyphy
    APPEND PROPERTY COMPILE_DEFINITIONS __MP__ __MP2__
)


set_property(
    TARGET hyphy HYPHYDEBUG HYPHYPROFILE
    APPEND PROPERTY COMPILE_DEFINITIONS _HYPHY_LIBDIRECTORY_="${LIBRARY_PATH}"
)


set_property(
    TARGET HYPHYDEBUG
    APPEND PROPERTY COMPILE_DEFINITIONS __HYPHYDEBUG__
)

set_target_properties(
    hyphy
    PROPERTIES
    if(${OPENMP_FOUND})
        COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${DEFAULT_WARNING_FLAGS} ${OpenMP_CXX_FLAGS}"
        LINK_FLAGS "${DEFAULT_LINK_FLAGS} ${OpenMP_CXX_FLAGS}"
    else(${OPENMP_FOUND})
        COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${DEFAULT_WARNING_FLAGS}"
        LINK_FLAGS "${DEFAULT_LINK_FLAGS}"
    endif(${OPENMP_FOUND})
)

macro(makeLink src dest target)
  add_custom_command(TARGET ${target} POST_BUILD
    COMMAND ln -sf ${dest} ${src} DEPENDS ${dest} COMMENT "mklink ${src} -> ${dest}")
endmacro()
makeLink(HYPHYMP hyphy hyphy)


#-------------------------------------------------------------------------------
# TESTS
#-------------------------------------------------------------------------------


add_test (NAME UNIT-TESTS COMMAND bash run_unit_tests.sh)
add_test (CODON HYPHYMP tests/hbltests/SimpleOptimizations/SmallCodon.bf)
add_test (PROTEIN HYPHYMP CPU=1 tests/hbltests/SimpleOptimizations/IntermediateProtein.bf)
add_test (MTCODON HYPHYMP tests/hbltests/libv3/mtDNA-code.wbf)
add_test (ALGAE HYPHYMP tests/hbltests/libv3/algae-mtDNA.wbf)
add_test (CILIATES HYPHYMP tests/hbltests/libv3/ciliate-code.wbf)
add_test (NAME SLAC COMMAND HYPHYMP tests/hbltests/libv3/SLAC.wbf)
add_test (NAME SLAC-PARTITIONED COMMAND HYPHYMP tests/hbltests/libv3/SLAC-partitioned.wbf)
add_test (NAME FEL COMMAND HYPHYMP tests/hbltests/libv3/FEL.wbf)
add_test (MEME HYPHYMP tests/hbltests/libv3/MEME.wbf)
#add_test (MEME_MPI mpirun -np 4 HYPHYMPI tests/hbltests/libv3/MEME.wbf)
add_test (MEME-PARTITIONED HYPHYMP tests/hbltests/libv3/MEME-partitioned.wbf)
add_test (BUSTED HYPHYMP tests/hbltests/libv3/BUSTED.wbf)
add_test (BUSTED-SRV HYPHYMP tests/hbltests/libv3/BUSTED-SRV.wbf)
add_test (RELAX HYPHYMP tests/hbltests/libv3/RELAX.wbf)
add_test (FUBAR HYPHYMP tests/hbltests/libv3/FUBAR.wbf)
add_test (BGM HYPHYMP tests/hbltests/libv3/BGM.wbf)
add_test (CONTRAST-FEL HYPHYMP tests/hbltests/libv3/CFEL.wbf)
add_test (GARD HYPHYMP tests/hbltests/libv3/GARD.wbf)
add_test (FADE HYPHYMP tests/hbltests/libv3/FADE.wbf)
add_test (ABSREL HYPHYMP tests/hbltests/libv3/ABSREL.wbf)
#add_test (FMM HYPHYMP tests/hbltests/libv3/FMM.wbf)
#add_test (NAME ABSREL-MH COMMAND HYPHYMP tests/hbltests/libv3/ABSREL-MH.wbf)

