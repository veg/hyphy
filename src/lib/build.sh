#!/bin/sh


TARGET_NAME="BLANK"
LIBRARY_BINDINGS=""

if [ $# -ne 1 -a $# -ne 2 ]
then
	TARGET_NAME="HELP";
else
	if [ $1 != "SP" -a $1 != "MP" -a $1 != "MP2"  -a $1 != "MPI" -a $1 != "DEBUG" -a $1 != "LIBRARY" ] 
	then
		$TARGET_NAME = "HELP"
	else
		TARGET_NAME=$1
	fi
fi

if [ $TARGET_NAME = "HELP" ] 
then
	echo "Usage: build.sh package_name"
	echo "  LIBRARY [Python|R]: multi-threaded library version with optional wrappers for Python or R."
	exit 1
fi

if [ $TARGET_NAME = "LIBRARY" -a $# -eq 2 ]
then
	if [ $2 != "R" -a $2 != "Python" ] 
	then
		echo "Library binding options must be one of the following:"
		echo "  LIBRARY [Python|R]: multi-threaded library version with optional wrappers for Python or R."
		exit 1
	else
		echo "Library bindings for $2 will be linked into the library. See README for details"
		LIBRARY_BINDINGS=$2
	fi
fi


# MODIFY THESE BASED ON YOUR SYSTEM
# DEFAULT SETTINGS ARE FOR GCC
	
COMPILER="g++";
COMPILERC="gcc";
sysName=`uname`;
echo $sysName;
CURL_LINKER_LIBS=" -lssl -lcrypto -lcurl";

if [ $sysName == "Darwin" ]
then
	machName=`machine`;
	if [ $machName == "ppc7450" ] 
	then
		COMPILER_FLAGS=" -D __UNIX__ -w -c -fsigned-char -fast -mcpu=7450 -fpermissive -I`pwd`/../Core -I`pwd`/../NewerFunctionality -I`pwd`/../../SQLite/trunk -D SQLITE_PTR_SIZE=sizeof(long) "
	else
		COMPILER_FLAGS=" -D __UNIX__ -w -c -fsigned-char -fast -fpermissive -I`pwd`/../Core -I`pwd`/../NewerFunctionality -I`pwd`/../../SQLite/trunk -D SQLITE_PTR_SIZE=sizeof(long) "	
	fi
	COMPILER_LINK_FLAGS=" -w -fsigned-char ";
else
	if [ $sysName == "AIX" ]
	then
		COMPILER="xlC";
		COMPILERC="xlc";
		COMPILER_FLAGS=" -c -qchar=signed -O3 -D SQLITE_PTR_SIZE=sizeof(long) -D __UNIX__ -I`pwd`/../Core -I`pwd`/../NewerFunctionality -I`pwd`/../../SQLite/trunk ";
		COMPILER_LINK_FLAGS="  -qchar=signed ";
	else
		COMPILER_LINK_FLAGS=" -w -fsigned-char ";
		COMPILER_FLAGS=" -w -c -fsigned-char -O3 -fpermissive -I`pwd`/../Core -I`pwd`/../NewerFunctionality -I`pwd`/../../SQLite/trunk -D SQLITE_PTR_SIZE=sizeof(long) -D __UNIX__ ";
	fi
fi


# END MODIFY

echo "Checking for curl";
echo "#include <curl/curl.h>\nint main(void) {return 0;}" > curl_check.cpp

if `$COMPILER -o curl_check -w $CURL_LINKER_LIBS curl_check.cpp`
then
 echo "Curl seems to be present"
else
	echo "Curl seems to be absent (setting up compiler options skip CURL code)";
	CURL_LINKER_LIBS="";
	COMPILER_FLAGS=$COMPILER_FLAGS" -D__HYPHY_NO_CURL__";
	COMPILER_LINK_FLAGS=$COMPILER_LINK_FLAGS" -D__HYPHY_NO_CURL__";
fi

rm -rf curl_check*

	
makedir () {
	if [ -f $1 ] 
	then
		echo "Insufficient permissions to create an object directory";
		exit 1;
	fi
	
	if [ ! -d $1  ]
	then
		if [ `mkdir $1` ]
		then
			echo "Failed to create directory $1";
			exit 1;
		fi
	fi
}

compileAll () {
	cd $1
	for fileName in *$2
	do
	  obj_file=$3/$OBJ_DIR_NAME/${fileName}.o;
	  if [ $obj_file -nt $fileName ]
	  then
		echo File "$fileName" is up to date
	  else
		  echo Building "$fileName";
		  if [ $2 = "c" ]
		  then 
		  	ccd=$COMPILERC
		  else
		  	ccd=$COMPILER
		  fi
		  if `$ccd -o $obj_file $COMPILER_FLAGS -fvisibility=hidden $fileName `
		   then
			 echo Complete
		   else
				echo Error during compilation;
				exit 1;
		   fi
	  fi
	done
	cd $3
}

OBJ_DIR_NAME="obj_$TARGET_NAME"

if [ -f $OBJ_DIR_NAME ] 
then
	rm -rf $OBJ_DIR_NAME;
fi

makedir $OBJ_DIR_NAME


if [ $1 = "LIBRARY" ] 
then
	LINKER_FLAGS=$CURL_LINKER_LIBS" -ldl -lm -lpthread ";
	echo "+-----------------------------------------------------------+"
	echo "|Building a multi-threaded HYPHY library version            |"
	echo "+-----------------------------------------------------------+"
	if [ $sysName == "Darwin" ]
	then
		TARGET_NAME="libhyphy.so";
		COMPILER_FLAGS=$COMPILER_FLAGS" -fno-strict-aliasing -D __MP__ -D __MP2__ -D __HEADLESS__ -fPIC -I`pwd`/Link "
		COMPILER_LINK_FLAGS=$COMPILER_LINK_FLAGS" -bundle -flat_namespace -undefined suppress "	
	else
		COMPILER_FLAGS=$COMPILER_FLAGS" -D __MP__ -D __MP2__ -D __HEADLESS__ -fPIC -I`pwd`/Link "
		COMPILER_LINK_FLAGS=$COMPILER_LINK_FLAGS" -Wl,-shared "
	fi
	
fi

COMPILER_FLAGS=$COMPILER_FLAGS" -I `pwd`/../Source"

echo "COMPILER=$COMPILER, $COMPILERC";
echo "COMPILER_FLAGS=$COMPILER_FLAGS";

compileAll ../Source cpp ../Library
compileAll ../GUI preferences.cpp ../Library
compileAll ../Source/SQLite c ../../Library

if [ $LIBRARY_BINDINGS = "R" ]
then	
	echo $COMPILER_FLAGS
	compileAll Link THyPhy.cpp ../
	echo Linking HyPhy.so
	cppf="PKG_CPPFLAGS=\"-I`pwd`/Link/\""
	#echo  $cppf R CMD SHLIB -o LibraryModules/R/HyPhy.so $OBJ_DIR_NAME/*.o    Link/Source/THyPhy_R.cpp $LINKER_FLAGS
	export $cppf; R CMD SHLIB -o LibraryModules/R/HyPhy.so $OBJ_DIR_NAME/*.o  SWIGWrappers/THyPhy_R.cpp	$LINKER_FLAGS	
	echo R library written to `pwd`/LibraryModules/R/HyPhy.so
fi

echo Finished


