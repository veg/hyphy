cmake_minimum_required(VERSION 3.15)
project(HyPhy)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
# set the CMAKE_MODULE_PATH so we can import OpenCL
#-------------------------------------------------------------------------------
set(CMAKE_MODULE_PATH cmake)
set(HYPHY_VERSION 2.1)

#-------------------------------------------------------------------------------
# setup the files we'll be using
#-------------------------------------------------------------------------------
set(SRC_FILES
    src/core/associative_list.cpp
    src/core/avllist.cpp
    src/core/avllistxl.cpp
    src/core/batchlan.cpp
    src/core/batchlanguage.cpp
    src/core/bgm.cpp
    src/core/calcnode.cpp
    src/core/calcnode2.cpp
    src/core/category.cpp
    src/core/chardata.cpp
    src/core/chart.cpp
    src/core/codondata.cpp
    src/core/dataloader.cpp
    src/core/dataset.cpp
    src/core/dataset_filter.cpp
    src/core/defines.cpp
    src/core/formula.cpp
    src/core/function_templates.cpp
    src/core/global_things.cpp
    src/core/hyphy_strings.cpp
    src/core/hyphymatrix.cpp
    src/core/hyphytree.cpp
    src/core/likefunc.cpp
    src/core/likelihoodfunction.cpp
    src/core/list.cpp
    src/core/list_iterator.cpp
    src/core/matrix.cpp
    src/core/matrix_op_cmds.cpp
    src/core/mersenne_twister.cpp
    src/core/multihit.cpp
    src/core/operation.cpp
    src/core/optimizer.cpp
    src/core/polynoml.cpp
    src/core/scfg.cpp
    src/core/scfg_ml.cpp
    src/core/sequencedata.cpp
    src/core/site.cpp
    src/core/stack.cpp
    src/core/string_file_wrapper.cpp
    src/core/trie.cpp
    src/core/variable.cpp
    src/core/variable_container.cpp
    src/new/branch_attribute.cpp
    src/new/model.cpp
    src/new/parameter.cpp
    src/contrib/sqlite3.cpp
    src/lib/link/THyPhy.cpp
    src/mains/unix.cpp
    src/utils/hyphyunixutils.cpp
)
set_source_files_properties(${SRC_CORE} ${SRC_NEW} {SRC_UTILS} PROPERTIES COMPILE_FLAGS "-Weffc++ -Wextra -Wall")


#-------------------------------------------------------------------------------
# make it fast, scotty
#-------------------------------------------------------------------------------

set(DEFAULT_WARNING_FLAGS " -w -Weffc++ -Wextra -Wall ")
set(DEFAULT_DEBUG_WARNING_FLAGS "-Wall -Wno-int-to-pointer-cast -Wno-conversion-null -Wno-sign-compare")

if (DEBUGFLAGS)
	set(ADDITIONAL_FLAGS "-pg -g")
endif()

#include_directories("/Library/Developer/CommandLineTools/SDKs/MacOSX15.2.sdk/usr/include/c++/v1")




if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

target_compile_options(hyphy PRIVATE
  $<$<CONFIG:Release>:-O3>
  $<$<CONFIG:Debug>:-g>
  $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
)

if(NOT NOAVX)
  target_compile_features(hyphy PRIVATE avx2 fma)
endif()

if(NOT NOSSE4)
  target_compile_features(hyphy PRIVATE sse4.1)
endif()

if(NOT NONEON)
  target_compile_features(hyphy PRIVATE neon)
endif()

if(NOT NOZLIB)
  find_package(ZLIB)
  if(ZLIB_FOUND)
    target_link_libraries(hyphy PRIVATE ZLIB::ZLIB)
    target_compile_definitions(hyphy PRIVATE -D__ZLIB__)
  endif()
endif()

find_package(CURL)
if(CURL_FOUND)
  target_link_libraries(hyphy PRIVATE CURL::CURL)
  target_compile_definitions(hyphy PRIVATE -D__HYPHYCURL__)
endif()

find_package(OpenMP)
if(OpenMP_FOUND)
  target_link_libraries(hyphy PRIVATE OpenMP::OpenMP_CXX)
endif()




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

# SQLite removed

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
    ${SRC_FILES}
)




#-------------------------------------------------------------------------------
# TESTS
#-------------------------------------------------------------------------------

set(HYPHY_ENV_ARG "ENV=TOLERATE_NUMERICAL_ERRORS=1;")


add_test (NAME UNIT-TESTS COMMAND bash run_unit_tests.sh)
add_test (CODON HYPHYMP tests/hbltests/SimpleOptimizations/SmallCodon.bf)
add_test (PROTEIN HYPHYMP tests/hbltests/SimpleOptimizations/IntermediateProtein.bf)
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
add_test (FADE HYPHYMP tests/hbltests/libv3/FADE.wbf)
add_test (NAME GARD 
          COMMAND HYPHYMP 
          tests/hbltests/libv3/GARD.wbf 
          ${HYPHY_ENV_ARG}
          )
add_test (ABSREL HYPHYMP tests/hbltests/libv3/ABSREL.wbf)
#add_test (FMM HYPHYMP tests/hbltests/libv3/FMM.wbf)
#add_test (NAME ABSREL-MH COMMAND HYPHYMP tests/hbltests/libv3/ABSREL-MH.wbf)

