cmake_minimum_required(VERSION 3.15)
project(HyPhy LANGUAGES CXX)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9")
        message(FATAL_ERROR "GCC version 9 or later is required.")
    endif()
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

enable_testing()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#------------------------------------------------------------------------------
# OPTIONS
#------------------------------------------------------------------------------

option(NOAVX "Disable AVX instruction set" OFF)
option(NOAVX2 "Disable AVX2 instruction set" OFF)
option(NOSSE4 "Disable SSE4 instruction set" OFF)
option(NONEON "Disable NEON instruction set" OFF)
option(NOZLIB "Disable ZLIB compression library" OFF)
option(NOBLAS "Disable BLAS linear algebra library" OFF)
option(RUN_CLANG_TIDY "Run clang-tidy on all files" OFF)
option(CODE_COVERAGE "Enable code coverage" OFF)
option(USE_ASAN "Enable Address Sanitizer" OFF)
option(USE_GCC_ANALYZER "Enable GCC static analyzer" OFF)
option(TARGET_NATIVE "Target the native CPU" ON)
option(INCLUDE_DEBUG_INFO "Include debug info in all builds" OFF)

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
include(CheckCXXSourceRuns)

set(COMMON_COMPILE_FLAGS "-fsigned-char;-O3;-D_FORTIFY_SOURCE=2;-Wall;")
set(HYPHY_SIMD_FLAGS "")
if(INCLUDE_DEBUG_INFO)
    list(APPEND COMMON_COMPILE_FLAGS "-g")
endif()

if(USE_GCC_ANALYZER)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        list(APPEND COMMON_COMPILE_FLAGS "-fanalyzer")
    endif()
endif()

if(TARGET_NATIVE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
        if(COMPILER_SUPPORTS_MARCH_NATIVE)
            list(APPEND COMMON_COMPILE_FLAGS "-march=native")
        endif()
        check_cxx_compiler_flag("-mtune=native" COMPILER_SUPPORTS_MTUNE_NATIVE)
        if(COMPILER_SUPPORTS_MTUNE_NATIVE)
            list(APPEND COMMON_COMPILE_FLAGS "-mtune=native")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        check_cxx_compiler_flag("-mtune=native" COMPILER_SUPPORTS_MTUNE_NATIVE)
        if(COMPILER_SUPPORTS_MTUNE_NATIVE)
            list(APPEND COMMON_COMPILE_FLAGS "-mtune=native")
        endif()
    endif()
endif()

if(NOT NONEON)
    set(OLD_CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_FLAGS "")
    check_cxx_source_runs("#include <arm_neon.h>\nint main() {\n  float32x4_t vec = vdupq_n_f32(0.0f);\n  (void)vec;\n  return 0;\n}\n" HYPHY_HAVE_NEON)
    if(HYPHY_HAVE_NEON)
        add_compile_definitions(_SLKP_USE_ARM_NEON)
    endif()
    set(CMAKE_REQUIRED_FLAGS ${OLD_CMAKE_REQUIRED_FLAGS})
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    list(APPEND COMMON_COMPILE_FLAGS "-Wno-overloaded-virtual")
    list(APPEND COMMON_COMPILE_FLAGS "-Wno-maybe-uninitialized")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(NOT HYPHY_HAVE_NEON AND NOT CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
    endif()
    if(NOT NOAVX)
        set(OLD_CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
        set(CMAKE_REQUIRED_FLAGS "-mavx")
        check_cxx_source_runs("#include <immintrin.h>\nint main() {\n  __m256d vec = _mm256_setzero_pd();\n  (void)vec;\n  return 0;\n}\n" HYPHY_HAVE_AVX)
        if(HYPHY_HAVE_AVX)
            list(APPEND HYPHY_SIMD_FLAGS "-mavx")
            add_compile_definitions(_SLKP_USE_AVX_INTRINSICS)
            if (NOT NOAVX2)
                set(CMAKE_REQUIRED_FLAGS "-mfma -mavx2")
                check_cxx_source_runs("#include <immintrin.h>\nint main() {\n  __m256d a = _mm256_setzero_pd(); __m256d b = _mm256_setzero_pd(); __m256d c = _mm256_setzero_pd(); a = _mm256_fmadd_pd(a, b, c); return 0;\n}\n" HYPHY_HAVE_FMA)
                if(HYPHY_HAVE_FMA)
                    list(APPEND HYPHY_SIMD_FLAGS "-mfma")
                    list(APPEND HYPHY_SIMD_FLAGS "-mavx2")
                    add_compile_definitions(_SLKP_USE_FMA3_INTRINSICS)
                endif()
            endif()
        endif()
        set(CMAKE_REQUIRED_FLAGS ${OLD_CMAKE_REQUIRED_FLAGS})
    elseif(NOT NOSSE4)
        set(OLD_CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
        set(CMAKE_REQUIRED_FLAGS "-msse4.1")
        check_cxx_source_runs("#include <smmintrin.h>\nint main() {\n  __m128i vec = _mm_setzero_si128();\n  (void)vec;\n  return 0;\n}\n" HYPHY_HAVE_SSE41)
        if(HYPHY_HAVE_SSE41)
            list(APPEND HYPHY_SIMD_FLAGS "-msse4.1")
            add_compile_definitions(_SLKP_USE_SSE_INTRINSICS)
        endif()
        set(CMAKE_REQUIRED_FLAGS ${OLD_CMAKE_REQUIRED_FLAGS})
    else()
        if(NOT CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
            list(APPEND COMMON_COMPILE_FLAGS "-mno-sse3")
        endif()
    endif()

endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    list(APPEND COMMON_COMPILE_FLAGS
        -Wextra
        -Wpedantic
        -Wshadow
        -Wformat=2
        -Wunused
        -Wmost
    )
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
   list(APPEND COMMON_COMPILE_FLAGS "-fwasm-exceptions")
   if(NOT NOSSE4)
       list(APPEND HYPHY_SIMD_FLAGS "-msimd128")
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
target_compile_options(hyphy PRIVATE ${COMMON_COMPILE_FLAGS} ${HYPHY_SIMD_FLAGS})
target_compile_definitions(hyphy PRIVATE
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
    add_executable(HYPHYMPI ${SRC_COMMON} ${SRC_UNIXMAIN})
    target_compile_options(HYPHYMPI PRIVATE ${COMMON_COMPILE_FLAGS} ${HYPHY_SIMD_FLAGS})
    target_compile_definitions(HYPHYMPI PRIVATE
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
    if(ZLIB_FOUND AND NOT NOZLIB)
        target_link_libraries(HYPHYMPI PRIVATE ${ZLIB_LIBRARIES})
        target_include_directories(HYPHYMPI PRIVATE ${ZLIB_INCLUDE_DIRS})
        target_compile_definitions(HYPHYMPI PRIVATE __ZLIB__)
    endif()
    install(TARGETS HYPHYMPI RUNTIME DESTINATION bin)
endif()

add_executable(HYPHYDEBUG EXCLUDE_FROM_ALL ${SRC_COMMON} ${SRC_UNIXMAIN})
target_compile_options(HYPHYDEBUG PRIVATE ${HYPHY_SIMD_FLAGS})
if(USE_ASAN)
    target_compile_options(HYPHYDEBUG PRIVATE -fsanitize=address -g)
    target_link_libraries(HYPHYDEBUG PRIVATE -fsanitize=address -g)
endif()
target_compile_definitions(HYPHYDEBUG PRIVATE
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
if(ZLIB_FOUND AND NOT NOZLIB)
    target_link_libraries(HYPHYDEBUG PRIVATE ${ZLIB_LIBRARIES})
    target_include_directories(HYPHYDEBUG PRIVATE ${ZLIB_INCLUDE_DIRS})
    target_compile_definitions(HYPHYDEBUG PRIVATE __ZLIB__)
endif()
if(CODE_COVERAGE)
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(HYPHYDEBUG PRIVATE -fprofile-instr-generate -fcoverage-mapping)
        target_link_libraries(HYPHYDEBUG PRIVATE -fprofile-instr-generate -fcoverage-mapping)
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(HYPHYDEBUG PRIVATE  -fprofile-arcs -ftest-coverage)
        target_link_libraries(HYPHYDEBUG PRIVATE  -fprofile-arcs -ftest-coverage)
    endif()
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