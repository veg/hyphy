cmake_minimum_required(VERSION 2.8)
project(HyPhy)

set(DEFAULT_COMPILE_FLAGS "-w -c -fsigned-char -fast -O3")

set(CMAKE_INSTALL_PREFIX /usr/local CACHE PATH "Installation prefix")

add_subdirectory(contrib/gtest-1.6.0)
set_property(
    TARGET gtest gtest_main
    PROPERTY EXCLUDE_FROM_ALL true
)

# file(GLOB_RECURSE HYPHY_STDLIB res/)
# add_custom_command(OUTPUT ${HYPHY_STDLIB}

install(DIRECTORY res/ DESTINATION lib/hyphy)

# uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY
)
add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake
)

file(GLOB SRC_CORE src/core/*.cpp)
file(GLOB SRC_NEW src/new/*.cpp)
file(GLOB SRC_GTESTS tests/gtests/*.cpp)

set(SRC_LINK src/lib/Link/THyPhy.cpp)
set(SRC_PREFS src/gui/preferences.cpp)
set(SRC_SQLITE3 contrib/SQLite-3.6.17/sqlite3.c)
set(SRC_UNIXMAIN src/mains/unix.cpp)
set(SRC_UTILS src/utils/hyphyunixutils.cpp)

set(SRC_COMMON ${SRC_CORE} ${SRC_NEW} ${SRC_SQLITE3} ${SRC_UTILS})

set_property(
    SOURCE ${SRC_SQLITE3}
    APPEND PROPERTY COMPILE_FLAGS "-Wno-int-to-pointer-cast -Wno-pointer-to-int-cast"
)

add_library(hyphy_mp SHARED
    EXCLUDE_FROM_ALL
    ${SRC_COMMON} ${SRC_LINK} ${SRC_PREFS}
)
target_link_libraries(hyphy_mp curl ssl crypto)
install(TARGETS hyphy_mp
    LIBRARY DESTINATION lib
    OPTIONAL
)
add_custom_target(LIB DEPENDS hyphy_mp)

add_executable(HYPHYMP
    EXCLUDE_FROM_ALL
    ${SRC_COMMON} ${SRC_UNIXMAIN}
)
target_link_libraries(HYPHYMP curl ssl crypto)
install(TARGETS HYPHYMP
    RUNTIME DESTINATION bin
    OPTIONAL
)
add_custom_target(MP2 DEPENDS HYPHYMP)

add_executable(HYPHYMPI
    EXCLUDE_FROM_ALL
    ${SRC_COMMON} ${SRC_UNIXMAIN}
)
target_link_libraries(HYPHYMPI curl ssl crypto)
install(TARGETS HYPHYMPI
    RUNTIME DESTINATION bin
    OPTIONAL
)
add_custom_target(MPI DEPENDS HYPHYMPI)

add_executable(HYPHYGTEST
    EXCLUDE_FROM_ALL
    ${SRC_COMMON} ${SRC_GTESTS} ${SRC_UNIXMAIN}
)
add_dependencies(HYPHYGTEST gtest)
target_link_libraries(HYPHYGTEST curl ssl crypto gtest)
add_custom_target(GTEST DEPENDS HYPHYGTEST)

include_directories(src/core/include src/gui/include src/lib/Link src/new/include contrib contrib/gtest-1.6.0/include tests/gtests)
add_definitions(-D_SLKP_LFENGINE_REWRITE_ -DINTPTR_TYPE=long)

if(UNIX)
    add_definitions(-D__UNIX__)
endif(UNIX)

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    add_definitions(-D__HYPHY_64__)
endif(CMAKE_SIZEOF_VOID_P EQUAL 8)

set_property(
    TARGET hyphy_mp
    APPEND PROPERTY COMPILE_DEFINITIONS __HEADLESS__
)

set_property(
    TARGET hyphy_mp HYPHYMP
    APPEND PROPERTY COMPILE_DEFINITIONS __MP__ __MP2__
)

set_property(
    TARGET HYPHYMPI
    PROPERTY COMPILE_DEFINITIONS __HYPHYMPI__
)

set_property(
    TARGET HYPHYGTEST
    APPEND PROPERTY COMPILE_DEFINITIONS __UNITTEST__
)

set_property(
    TARGET HYPHYMP HYPHYMPI HYPHYGTEST
    APPEND PROPERTY COMPILE_DEFINITIONS _HYPHY_LIBDIRECTORY_="${CMAKE_INSTALL_PREFIX}/lib/hyphy"
)

set_property(
    TARGET hyphy_mp HYPHYMP HYPHYMPI
    PROPERTY COMPILE_FLAGS ${DEFAULT_COMPILE_FLAGS}
)

set_property(
    TARGET HYPHYGTEST
    PROPERTY COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -g -fprofile-arcs -ftest-coverage"
)
