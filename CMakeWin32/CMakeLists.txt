#-------------------------------------------------------------------------------
# hyphy Win32 target for compiling with OS X hosted mingw-gcc
# Download the mingw toolchain from and install into /usr/local
# http://crossgcc.rts-software.org/
#-------------------------------------------------------------------------------


cmake_minimum_required(VERSION 2.8)
project(HyPhyWin32 NONE)

set(CMAKE_CROSSCOMPILING true)
set(CMAKE_CONFIGURATION_TYPES Release)


# this one is important
SET(CMAKE_SYSTEM_NAME Generic)
#this one not so much
SET(CMAKE_SYSTEM_VERSION 1)

# search for programs in the build host directories
SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
# for libraries and headers in the target directories
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)


# where is the target environment 
SET(CMAKE_FIND_ROOT_PATH  /usr/local/i386-mingw32-4.3.0/)

# specify the cross compiler
SET(CMAKE_C_COMPILER   ${CMAKE_FIND_ROOT_PATH}/bin/i386-mingw32-gcc)
ENABLE_LANGUAGE(C)
SET(CMAKE_CXX_COMPILER ${CMAKE_FIND_ROOT_PATH}/bin/i386-mingw32-g++)
ENABLE_LANGUAGE(CXX)

set(DEFAULT_COMPILE_FLAGS "-fsigned-char -O3 ")
add_definitions(-D_SLKP_LFENGINE_REWRITE_ -D__AFYP_REWRITE_BGM__)


set(CMAKE_RC_COMPILER_INIT ${CMAKE_FIND_ROOT_PATH}/bin/i386-mingw32-windres)
ENABLE_LANGUAGE(RC)
SET(CMAKE_RC_COMPILE_OBJECT
"<CMAKE_RC_COMPILER> <FLAGS> <DEFINES> -O coff -o <OBJECT> -i <SOURCE> ")

#-------------------------------------------------------------------------------
# setup the files we'll be using
#-------------------------------------------------------------------------------
file(GLOB SRC_CORE ../src/core/*.cpp)
file(GLOB SRC_GUI ../src/gui/*.cpp ../src/gui/Components/*.cpp)
file(GLOB SRC_GTESTS ../tests/gtests/*.cpp)
file(GLOB SRC_NEW ../src/new/*.cpp)
file(GLOB SRC_WIN32_GUI ../src/gui/win/*.cpp ../src/gui/win/Components/*.cpp ../src/gui/win/WindowClasses/*.cpp)

set(SRC_SQLITE3 ../contrib/SQLite-3.8.2/sqlite3.c)
set(SRC_WINMAIN ../src/mains/win.cpp)
set(SRC_WINRES  ../src/gui/res/Windows/Win.rc)

set(SRC_COMMON ${SRC_CORE} ${SRC_NEW} ${SRC_SQLITE3})

include_directories(
	../src/core/include
	../src/lib/Link 
	../src/new/include
	../contrib/SQLite-3.8.2
	../src/gui/include
	../src/gui/include/Components
	../src/gui/win/include
	../src/gui/win/include/Components
	../src/gui/include/WindowClasses
)

#-------------------------------------------------------------------------------
# hyphy win32 GUI target
#-------------------------------------------------------------------------------
add_executable(
	HYPHY.EXE
	WIN32
	${SRC_WINRES}
	${SRC_WINMAIN}
	${SRC_WIN32_GUI}
	${SRC_COMMON} 
	${SRC_GUI}
)
	
target_link_libraries(HYPHY.EXE comctl32 winspool wininet msimg32)


add_custom_target(WIN32 DEPENDS HYPHY.EXE
				  COMMENT "Win32 GUI HyPhy"	
)

find_program (HAZ_NSIS "makensis" 
			HINTS /usr/bin/ /usr/local/bin ENV PATH)


if (${HAZ_NSIS} MATCHES ".+makensis$")
	add_custom_command(
		TARGET HYPHY.EXE
		POST_BUILD
		COMMAND makensis ${CMAKE_CURRENT_SOURCE_DIR}/../installers/Win/hyphy.nsi
		)	
endif (${HAZ_NSIS} MATCHES ".+makensis$")
 

if (${CMAKE_SYSTEM} MATCHES "^Darwin.*")
	set (CMAKE_CXX_LINK_FLAGS "")
endif (${CMAKE_SYSTEM} MATCHES "^Darwin.*")

set_target_properties(
	HYPHY.EXE
	PROPERTIES
	COMPILE_DEFINITIONS "__WINDOZE__;WINVER=0x0500;__HYPHY_NO_CURL__;"
	COMPILE_FLAGS "${DEFAULT_COMPILE_FLAGS} -mwindows "
	LINK_FLAGS "${DEFAULT_COMPILE_FLAGS} -mwindows -v -static "
	LINKER_LANGUAGE CXX
	WIN32_EXECUTABLE TRUE
)

