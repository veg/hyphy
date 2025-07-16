cmake_minimum_required(VERSION 3.15)
project(HyPhy LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

enable_testing()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#------------------------------------------------------------------------------
# OPTIONS
#------------------------------------------------------------------------------

option(NOAVX "Disable AVX instruction set" OFF)
option(NOSSE4 "Disable SSE4 instruction set" OFF)
option(NONEON "Disable NEON instruction set" OFF)
option(NOZLIB "Disable ZLIB compression library" OFF)
option(NOBLAS "Disable BLAS linear algebra library" OFF)
option(RUN_CLANG_TIDY "Run clang-tidy on all files" OFF)

if(RUN_CLANG_TIDY)
  find_program(CLANG_TIDY_EXE NAMES clang-tidy)
  if(CLANG_TIDY_EXE)
    set(CMAKE_CXX_CLANG_TIDY ${CLANG_TIDY_EXE})
  else()
    message(SEND_ERROR "clang-tidy not found")
  endif()
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
   set(NOAVX ON)
   set(NONEON ON)
   set(NOZLIB ON)
   set(LIBRARY_PATH "/hyphy")
   add_definitions (-D_USE_EMSCRIPTEN_)
 endif()

#------------------------------------------------------------------------------
# Source Files
#------------------------------------------------------------------------------
# For larger projects, it is recommended to list all source files explicitly
# instead of using file(GLOB). This ensures that adding/removing files
# is explicitly captured by the build system.
file(GLOB SRC_CORE src/core/*.cpp)
file(GLOB SRC_NEW src/new/*.cpp)
file(GLOB SRC_CONTRIB src/contrib/*.cpp)

set(SRC_UNIXMAIN src/mains/unix.cpp)
set(SRC_UTILS src/utils/hyphyunixutils.cpp)

set(SRC_COMMON ${SRC_CORE} ${SRC_NEW} ${SRC_UTILS} ${SRC_CONTRIB})

#------------------------------------------------------------------------------
# Compiler specific flags
#------------------------------------------------------------------------------
include(CheckCXXCompilerFlag)
include(CheckCXXSourceCompiles)

set(COMMON_COMPILE_FLAGS "-fsigned-char;-O3;-D_FORTIFY_SOURCE=2;-Wall;-g;")

if(NOT NONEON)
    set(OLD_CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_FLAGS "")
    check_cxx_source_compiles("
#include <arm_neon.h>
int main() {
  float32x4_t vec = vdupq_n_f32(0.0f);
  (void)vec;
  return 0;
}
" HYPHY_HAVE_NEON)
    if(HYPHY_HAVE_NEON)
        add_compile_definitions(_SLKP_USE_ARM_NEON)
        if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            list(APPEND COMMON_COMPILE_FLAGS "-mcpu=native;-mtune=native")
        endif()
    endif()
    set(CMAKE_REQUIRED_FLAGS ${OLD_CMAKE_REQUIRED_FLAGS})
endif()

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(NOT NOAVX)
        set(OLD_CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
        set(CMAKE_REQUIRED_FLAGS "-mavx2")
        check_cxx_source_compiles("
#include <immintrin.h>
int main() {
  __m256d vec = _mm256_setzero_pd();
  (void)vec;
  return 0;
}
" HYPHY_HAVE_AVX2)
        if(HYPHY_HAVE_AVX2)
            list(APPEND COMMON_COMPILE_FLAGS "-mavx2")
            add_compile_definitions(_SLKP_USE_AVX_INTRINSICS)
            set(CMAKE_REQUIRED_FLAGS "-mfma")
            check_cxx_source_compiles("
#include <immintrin.h>
int main() {
  __m256d vec = _mm256_setzero_pd();
  (void)vec;
  return 0;
}
" HYPHY_HAVE_FMA)
            if(HYPHY_HAVE_FMA)
                list(APPEND COMMON_COMPILE_FLAGS "-mfma")
                add_compile_definitions(_SLKP_USE_FMA3_INTRINSICS)
            endif()
        endif()
        set(CMAKE_REQUIRED_FLAGS ${OLD_CMAKE_REQUIRED_FLAGS})
    elseif(NOT NOSSE4)
        set(OLD_CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
        set(CMAKE_REQUIRED_FLAGS "-msse4.1")
        check_cxx_source_compiles("
#include <smmintrin.h>
int main() {
  __m128i vec = _mm_setzero_si128();
  (void)vec;
  return 0;
}
" HYPHY_HAVE_SSE41)
        if(HYPHY_HAVE_SSE41)
            list(APPEND COMMON_COMPILE_FLAGS "-msse4.1")
            add_compile_definitions(_SLKP_USE_SSE_INTRINSICS)
        endif()
        set(CMAKE_REQUIRED_FLAGS ${OLD_CMAKE_REQUIRED_FLAGS})
    else()
        if(NOT CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
            list(APPEND COMMON_COMPILE_FLAGS "-mno-sse3")
        endif()
    endif()
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
   list(APPEND COMMON_COMPILE_FLAGS "-fwasm-exceptions")
   if(NOT NOSSE4)
       list(APPEND COMMON_COMPILE_FLAGS "-msimd128")
   endif()
endif()

#------------------------------------------------------------------------------
# Dependencies
#------------------------------------------------------------------------------
find_package(OpenMP)
find_package(BLAS)
find_package(CURL)
find_package(ZLIB)
find_package(MPI)

#------------------------------------------------------------------------------
# Installation
#------------------------------------------------------------------------------
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX "/usr/local" CACHE PATH "Installation prefix")
endif()

if (NOT DEFINED LIBRARY_PATH)
    set(LIBRARY_PATH "${CMAKE_INSTALL_PREFIX}/share/hyphy")
endif()

install(DIRECTORY res/ DESTINATION share/hyphy)

#------------------------------------------------------------------------------
# Uninstall Target
#------------------------------------------------------------------------------
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY
)
add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
)

#------------------------------------------------------------------------------
# Main Executable
#------------------------------------------------------------------------------
add_executable(hyphy ${SRC_COMMON} ${SRC_UNIXMAIN})
target_compile_options(hyphy PRIVATE ${COMMON_COMPILE_FLAGS})
target_compile_definitions(hyphy PRIVATE
    _SLKP_LFENGINE_REWRITE_
    __AFYP_REWRITE_BGM__
    __UNIX__
    __MP__
    __MP2__
    _HYPHY_LIBDIRECTORY_="${LIBRARY_PATH}"
)
target_include_directories(hyphy PRIVATE
    src/core/include
    src/contrib
    src/lib/Link
    src/new/include
)

if(CMAKE_SYSTEM_NAME MATCHES "FreeBSD")
    set(DEFAULT_LIBRARIES "")
else()
    set(DEFAULT_LIBRARIES dl)
endif()

if(CURL_FOUND)
    target_link_libraries(hyphy PRIVATE ${CURL_LIBRARIES})
    target_compile_definitions(hyphy PRIVATE __HYPHYCURL__)
endif()

if(ZLIB_FOUND AND NOT NOZLIB)
    target_link_libraries(hyphy PRIVATE ${ZLIB_LIBRARIES})
    target_include_directories(hyphy PRIVATE ${ZLIB_INCLUDE_DIRS})
    target_compile_definitions(hyphy PRIVATE __ZLIB__)
endif()

if(BLAS_FOUND AND NOT NOBLAS)
    target_link_libraries(hyphy PRIVATE ${BLAS_LIBRARIES})
    target_compile_definitions(hyphy PRIVATE _SLKP_USE_APPLE_BLAS)
endif()

if(OPENMP_FOUND)
    target_link_libraries(hyphy PRIVATE OpenMP::OpenMP_CXX)
endif()

target_link_libraries(hyphy PRIVATE ${DEFAULT_LIBRARIES})

install(TARGETS hyphy RUNTIME DESTINATION bin)

#------------------------------------------------------------------------------
# Other Executables
#------------------------------------------------------------------------------
if(MPI_FOUND)
    add_executable(HYPHYMPI EXCLUDE_FROM_ALL ${SRC_COMMON} ${SRC_UNIXMAIN})
    target_compile_options(HYPHYMPI PRIVATE ${COMMON_COMPILE_FLAGS})
    target_compile_definitions(HYPHYMPI PRIVATE
        _SLKP_LFENGINE_REWRITE_
        __AFYP_REWRITE_BGM__
        __UNIX__
        __HYPHYMPI__
        _HYPHY_LIBDIRECTORY_="${LIBRARY_PATH}"
    )
    target_include_directories(HYPHYMPI PRIVATE
        src/core/include
        src/contrib
        src/lib/Link
        src/new/include
        ${MPI_INCLUDE_PATH}
    )
    target_link_libraries(HYPHYMPI PRIVATE
        ${DEFAULT_LIBRARIES}
        ${MPI_LIBRARIES}
    )
    if(OPENMP_FOUND)
        target_link_libraries(HYPHYMPI PRIVATE OpenMP::OpenMP_CXX)
    endif()
    if(BLAS_FOUND AND NOT NOBLAS)
        target_link_libraries(HYPHYMPI PRIVATE ${BLAS_LIBRARIES})
    endif()
    install(TARGETS HYPHYMPI RUNTIME DESTINATION bin OPTIONAL)
endif()

add_executable(HYPHYDEBUG EXCLUDE_FROM_ALL ${SRC_COMMON} ${SRC_UNIXMAIN})
target_compile_options(HYPHYDEBUG PRIVATE
    ${COMMON_COMPILE_FLAGS}
    -g
    -fsanitize=address
    -O0
)
target_compile_definitions(HYPHYDEBUG PRIVATE
    _SLKP_LFENGINE_REWRITE_
    __AFYP_REWRITE_BGM__
    __UNIX__
    __HYPHYDEBUG__
    _HYPHY_LIBDIRECTORY_="${LIBRARY_PATH}"
)
target_include_directories(HYPHYDEBUG PRIVATE
    src/core/include
    src/contrib
    src/lib/Link
    src/new/include
)
target_link_libraries(HYPHYDEBUG PRIVATE ${DEFAULT_LIBRARIES})
if(OPENMP_FOUND)
    target_link_libraries(HYPHYDEBUG PRIVATE OpenMP::OpenMP_CXX)
endif()
if(BLAS_FOUND AND NOT NOBLAS)
    target_link_libraries(HYPHYDEBUG PRIVATE ${BLAS_LIBRARIES})
endif()

#------------------------------------------------------------------------------
# Tests
#------------------------------------------------------------------------------
add_test(NAME UNIT-TESTS COMMAND bash run_unit_tests.sh)
add_test(CODON hyphy tests/hbltests/SimpleOptimizations/SmallCodon.bf)
add_test(PROTEIN hyphy tests/hbltests/SimpleOptimizations/IntermediateProtein.bf)
add_test(MTCODON hyphy tests/hbltests/libv3/mtDNA-code.wbf)
add_test(ALGAE hyphy tests/hbltests/libv3/algae-mtDNA.wbf)
add_test(CILIATES hyphy tests/hbltests/libv3/ciliate-code.wbf)
add_test(NAME SLAC COMMAND hyphy tests/hbltests/libv3/SLAC.wbf)
add_test(NAME SLAC-PARTITIONED COMMAND hyphy tests/hbltests/libv3/SLAC-partitioned.wbf)
add_test(NAME FEL COMMAND hyphy tests/hbltests/libv3/FEL.wbf)
add_test(MEME hyphy tests/hbltests/libv3/MEME.wbf)
add_test(MEME-PARTITIONED hyphy tests/hbltests/libv3/MEME-partitioned.wbf)
add_test(BUSTED hyphy tests/hbltests/libv3/BUSTED.wbf)
add_test(BUSTED-SRV hyphy tests/hbltests/libv3/BUSTED-SRV.wbf)
add_test(RELAX hyphy tests/hbltests/libv3/RELAX.wbf)
add_test(FUBAR hyphy tests/hbltests/libv3/FUBAR.wbf)
add_test(BGM hyphy tests/hbltests/libv3/BGM.wbf)
add_test(CONTRAST-FEL hyphy tests/hbltests/libv3/CFEL.wbf)
add_test(FADE hyphy tests/hbltests/libv3/FADE.wbf)
add_test(NAME GARD COMMAND hyphy tests/hbltests/libv3/GARD.wbf "ENV=TOLERATE_NUMERICAL_ERRORS=1;")
add_test(ABSREL hyphy tests/hbltests/libv3/ABSREL.wbf)