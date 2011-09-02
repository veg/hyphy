cmake_minimum_required(VERSION 2.8.1)
project(HyPhy)

# make it fast, scotty
set(DEFAULT_COMPILE_FLAGS "-c -fsigned-char -fast -O3 -flto")

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
file(GLOB SRC_GUI src/gui/*.cpp src/gui/Components/*.cpp)
file(GLOB SRC_GTESTS tests/gtests/*.cpp)
file(GLOB SRC_MACGUI src/gui/mac/*.cpp src/gui/mac/Components/HYPlatformList.cpp)
file(GLOB SRC_NEW src/new/*.cpp)

set(SRC_LINK src/lib/Link/THyPhy.cpp)
set(SRC_PREFS src/gui/preferences.cpp)
set(SRC_SQLITE3 contrib/SQLite-3.6.17/sqlite3.c)
set(SRC_UNIXMAIN src/mains/unix.cpp)
set(SRC_MACMAIN src/mains/mac.cpp)
set(SRC_UTILS src/utils/hyphyunixutils.cpp)

set(SRC_COMMON ${SRC_CORE} ${SRC_NEW} ${SRC_SQLITE3} ${SRC_UTILS})

set_property(
    SOURCE ${SRC_SQLITE3}
    PROPERTY COMPILE_FLAGS "-Wno-int-to-pointer-cast -Wno-pointer-to-int-cast"
)

set_property(
    SOURCE src/core/calcnode2.cpp
    PROPERTY COMPILE_FLAGS "-fstrict-aliasing -funroll-loops"
)

add_library(
    hyphy_mp SHARED
    EXCLUDE_FROM_ALL
    ${SRC_COMMON} ${SRC_LINK} ${SRC_PREFS}
)
target_link_libraries(hyphy_mp curl ssl crypto)
install(
    TARGETS hyphy_mp
    LIBRARY DESTINATION lib
    OPTIONAL
)
add_custom_target(LIB DEPENDS hyphy_mp)

add_executable(
    HYPHYMP
    EXCLUDE_FROM_ALL
    ${SRC_COMMON} ${SRC_UNIXMAIN}
)
target_link_libraries(HYPHYMP curl ssl crypto)
install(
    TARGETS HYPHYMP
    RUNTIME DESTINATION bin
    OPTIONAL
)
add_custom_target(MP2 DEPENDS HYPHYMP)

add_executable(
    HYPHYMPI
    EXCLUDE_FROM_ALL
    ${SRC_COMMON} ${SRC_UNIXMAIN}
)
target_link_libraries(HYPHYMPI curl ssl crypto)
install(
    TARGETS HYPHYMPI
    RUNTIME DESTINATION bin
    OPTIONAL
)
add_custom_target(MPI DEPENDS HYPHYMPI)

add_executable(
    HYPHYDEBUG
    EXCLUDE_FROM_ALL
    ${SRC_COMMON} ${SRC_UNIXMAIN}
)
target_link_libraries(HYPHYDEBUG curl ssl crypto)
add_custom_target(DEBUG DEPENDS HYPHYDEBUG)

add_executable(
    HYPHYGTEST
    EXCLUDE_FROM_ALL
    ${SRC_COMMON} ${SRC_GTESTS} ${SRC_UNIXMAIN}
)
add_dependencies(HYPHYGTEST gtest)
target_link_libraries(HYPHYGTEST curl ssl crypto gtest)
add_custom_target(GTEST DEPENDS HYPHYGTEST)

include_directories(src/core/include src/lib/Link src/new/include contrib src/gui/include src/gui/include/Components src/gui/include/WindowClasses contrib/gtest-1.6.0/include tests/gtests)
add_definitions(-D_SLKP_LFENGINE_REWRITE_ -DINTPTR_TYPE=long)

if(APPLE AND FALSE)
    set(CMAKE_OSX_ARCHITECTURES i386)
    add_executable(HYPHYMACOSX
        MACOSX_BUNDLE
        EXCLUDE_FROM_ALL
        ${SRC_CORE} ${SRC_GUI} ${SRC_MACGUI} ${SRC_MACMAIN}
    )
    include_directories(src/gui/mac/include src/gui/mac/include/Components /Developer/Headers/FlatCarbon)
    target_link_libraries(HYPHYMACOSX curl ssl crypto)
    add_custom_target(MAC DEPENDS HYPHYMACOSX)
    set_property(
        TARGET HYPHYMACOSX
        APPEND PROPERTY COMPILE_DEFINITIONS __HYPHYXCODE__ __MAC__ __MP__ __MP2__
    )
    set_property(
        SOURCE src/gui/preferences.cpp src/gui/mac/hydialogs.cpp
        PROPERTY COMPILE_FLAGS "-fpascal-strings"
    )
    set_target_properties(
        HYPHYMACOSX
        PROPERTIES
        MACOSX_BUNDLE_INFO_PLIST src/gui/mac/Info.plist
        MACOSX_BUNDLE_ICON_FILE src/gui/mac/HYPHY.icns
    )
endif(APPLE AND FALSE)

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
    TARGET HYPHYMP HYPHYMPI HYPHYGTEST HYPHYDEBUG
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

set_property(
    TARGET HYPHYDEBUG
    PROPERTY COMPILE_FLAGS "-g -O0"
)
