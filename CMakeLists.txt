cmake_minimum_required(VERSION 2.8.1)
project(HyPhy)


set(CMAKE_CONFIGURATION_TYPES Release)


#-------------------------------------------------------------------------------
# set the CMAKE_MODULE_PATH so we can import OpenCL 
#-------------------------------------------------------------------------------
set(CMAKE_MODULE_PATH cmake)

find_package(OpenCL)

if(${OPENCL_FOUND})
    message("OpenCL FOUND!!!")
endif(${OPENCL_FOUND})


#-------------------------------------------------------------------------------
# figure out some system-stuff for compile-time workarounds
#-------------------------------------------------------------------------------
set(MACOSX_LION false)
if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    if(${CMAKE_SYSTEM_VERSION} VERSION_GREATER 11.0.0 OR ${CMAKE_SYSTEM_VERSION} VERSION_EQUAL 11.0.0)
        set(MACOSX_LION true)
    endif(${CMAKE_SYSTEM_VERSION} VERSION_GREATER 11.0.0 OR ${CMAKE_SYSTEM_VERSION} VERSION_EQUAL 11.0.0)
endif(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")


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

    if(${MACOSX_LION})
        set(GCC45 true)
    endif(${MACOSX_LION})

    if(${GCC45})
        set(DEFAULT_COMPILE_FLAGS "-fsigned-char -O3 -flto")
    else(${GCC45})
        set(DEFAULT_COMPILE_FLAGS "-fsigned-char -O3")
    endif(${GCC45})

    if(${GCC46})
        set(DEFAULT_WARNING_FLAGS "-Wno-int-to-pointer-cast -Wno-conversion-null")
    endif(${GCC46})

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


#-------------------------------------------------------------------------------
# gtest dependency
#-------------------------------------------------------------------------------
add_subdirectory(contrib/gtest-1.6.0)
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
file(GLOB SRC_GTESTS tests/gtests/*.cpp)
file(GLOB SRC_NEW src/new/*.cpp)

set(SRC_LINK src/lib/Link/THyPhy.cpp)
set(SRC_PREFS src/gui/preferences.cpp)
set(SRC_SQLITE3 contrib/SQLite-3.6.17/sqlite3.c)
set(SRC_UNIXMAIN src/mains/unix.cpp)
set(SRC_MACMAIN src/mains/mac.cpp)
set(SRC_UTILS src/utils/hyphyunixutils.cpp)

set(SRC_COMMON ${SRC_CORE} ${SRC_NEW} ${SRC_SQLITE3} ${SRC_UTILS})


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
add_definitions(-D_SLKP_LFENGINE_REWRITE_ -DINTPTR_TYPE=long)

if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    add_definitions(-D__HYPHY_64__)
endif(CMAKE_SIZEOF_VOID_P EQUAL 8)

include_directories(
    src/core/include
    src/lib/Link src/new/include
    contrib # SQLite
    src/gui/include
    src/gui/include/Components
    src/gui/include/WindowClasses
    contrib/gtest-1.6.0/include
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
target_link_libraries(hyphy_mp crypto curl pthread ssl)
install(
    TARGETS hyphy_mp
    LIBRARY DESTINATION lib
    OPTIONAL
)
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
target_link_libraries(HYPHYMP crypto curl pthread ssl)
install(
    TARGETS HYPHYMP
    RUNTIME DESTINATION bin
    OPTIONAL
)
add_custom_target(MP2 DEPENDS HYPHYMP)


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
    target_link_libraries(HYPHYMPI crypto curl pthread ssl ${MPI_LIBRARIES})
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
        COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}"
        LINK_FLAGS "${DEFAULT_LINK_FLAGS} ${MPI_LINK_FLAGS}"
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
target_link_libraries(HYPHYDEBUG crypto curl pthread ssl)
add_custom_target(DEBUG DEPENDS HYPHYDEBUG)

set_target_properties(
    HYPHYDEBUG
    PROPERTIES
    COMPILE_FLAGS "${OpenMP_CXX_FLAGS} -g -O0"
    LINK_FLAGS "${OpenMP_CXX_FLAGS} -g -O0"
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
target_link_libraries(HYPHYGTEST crypto curl gtest pthread ssl)
add_custom_target(GTEST DEPENDS HYPHYGTEST)

set_property(
    TARGET HYPHYGTEST
    APPEND PROPERTY COMPILE_DEFINITIONS __UNITTEST__
)
set_target_properties(
    HYPHYGTEST
    PROPERTIES
    COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS} -g -fprofile-arcs -ftest-coverage"
    LINK_FLAGS "${DEFAULT_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS} -g -fprofile-arcs -ftest-coverage"
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

set_target_properties(
    hyphy_mp HYPHYMP
    PROPERTIES
    COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} ${OpenMP_CXX_FLAGS}"
    LINK_FLAGS "${DEFAULT_LINK_FLAGS} ${OpenMP_CXX_FLAGS}"
)


#-------------------------------------------------------------------------------
# hyphy MacOSX gui target
#-------------------------------------------------------------------------------
if(APPLE)
    #---------------------------------------------------------------------------
    # set this or blow up
    #---------------------------------------------------------------------------
    set(CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.6.sdk/")


    #---------------------------------------------------------------------------
    # MacOSX gui files
    #---------------------------------------------------------------------------
    set(RES_MACGUI src/gui/mac/HyPhy.rsrc src/gui/mac/HYPHYMP.icns)
    file(
        GLOB SRC_MACGUI
        src/gui/mac/*.cpp
        src/gui/mac/Components/*.cpp
        src/gui/mac/WindowClasses/*.cpp
    )

    #---------------------------------------------------------------------------
    # MacOSX gui executable and directives
    #---------------------------------------------------------------------------
    add_executable(HyPhy
        MACOSX_BUNDLE
        EXCLUDE_FROM_ALL
        ${SRC_CORE} ${SRC_NEW} ${SRC_SQLITE3} ${SRC_GUI} ${SRC_MACGUI} ${SRC_MACMAIN}
        ${RES_MACGUI}
    )
    include_directories(src/gui/mac/include src/gui/mac/include/Components /Developer/Headers/FlatCarbon)
    find_library(CARBON NAMES Carbon)
    find_library(QUICKTIME NAMES Quicktime)
    target_link_libraries(HyPhy crypto curl pthread ssl ${CARBON} ${QUICKTIME})
    add_custom_target(MAC DEPENDS HyPhy)

    #---------------------------------------------------------------------------
    # install the hyphy hbl stdlib into the .app
    #---------------------------------------------------------------------------
    get_target_property(HYPHY_APPEXE HyPhy LOCATION)
    get_filename_component(HYPHY_APPDIR ${HYPHY_APPEXE} PATH)
    add_custom_command(
        TARGET HyPhy
        PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${HYPHY_APPDIR}/../Resources/HBL
    )
    add_custom_command(
        TARGET HyPhy
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/res ${HYPHY_APPDIR}/../Resources/HBL
    )

    #---------------------------------------------------------------------------
    # MacOSX gui-specific properties
    #---------------------------------------------------------------------------
    set_property(
        TARGET HyPhy
        APPEND PROPERTY COMPILE_DEFINITIONS __HYPHYXCODE__ __MAC__ __MP__ __MP2__ __BUNDLE__
    )
    set_target_properties(
        HyPhy
        PROPERTIES
        COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -Wno-char-subscripts -Wno-deprecated -Wno-deprecated-declarations -Wno-sign-compare"
        LINK_FLAGS "${DEFAULT_COMPILE_FLAGS} -Wno-char-subscripts -Wno-deprecated -Wno-deprecated-declarations -Wno-sign-compare"
    )
    set_property(
        SOURCE
        src/gui/mac/hydialogs.cpp
        src/gui/mac/HYPlatformComponent.cpp
        src/gui/mac/HYPlatformUtils.cpp
        src/gui/mac/HYPlatformWindow.cpp
        src/gui/mac/Components/HYPlatformCheckBox.cpp
        src/gui/mac/Components/HYPlatformPullDown.cpp
        src/gui/mac/WindowClasses/HYPlatformChartWindow.cpp
        src/gui/mac/WindowClasses/HYPlatformConsoleWindow.cpp
        src/gui/mac/WindowClasses/HYPlatformDataPanel.cpp
        src/gui/mac/WindowClasses/HYPlatformModelWindow.cpp
        src/gui/mac/WindowClasses/HYPlatformParameterTable.cpp
        src/gui/mac/WindowClasses/HYPlatformPWindow.cpp
        src/gui/mac/WindowClasses/HYPlatformTreePanel.cpp
        src/gui/preferences.cpp
        src/mains/mac.cpp
        PROPERTY COMPILE_FLAGS "-fpascal-strings"
    )
    set_source_files_properties(
        ${RES_MACGUI}
        PROPERTIES MACOSX_PACKAGE_LOCATION Resources
    )
    set_target_properties(
        HyPhy
        PROPERTIES MACOSX_BUNDLE_INFO_PLIST src/gui/mac/Info.plist
    )
    set_target_properties(
        HyPhy MAC
        PROPERTIES
        OSX_ARCHITECTURES i386
    )
endif(APPLE)
