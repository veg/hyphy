#!/bin/sh


TARGET_NAME="BLANK"
LIBRARY_BINDINGS="NONE"

export PATH=$PATH:/opt/ibmcmp/xlc/ssc/0.9/bin/

if [ $# -ne 1 -a $# -ne 2 ]
then
	TARGET_NAME="HELP";
else
	if [ $1 != "SP" -a $1 != "GTEST" -a $1 != "MP" -a $1 != "MP2"  -a $1 != "MPI" -a $1 != "DEBUG" -a $1 != "LIBRARY" -a $1 != "DEV" -a $1 != "PS3" -a $1 != "DMALLOC" ] 
	then
		$TARGET_NAME = "HELP"
	else
		TARGET_NAME=$1
	fi
fi

if [ $TARGET_NAME = "HELP" ] 
then
	echo "Usage: build.sh package_name"
	echo "  where package_name is one of the following:"
	echo "  SP : for a single threaded build"
	echo "  MP : for a multi-threaded build with pthreads."
	echo "  MP2 : for a multi-threaded build with pthreads which support setconcurrency function."
	echo "  DEV : developmental OpenMP build with likelihood function speedups."
	echo "  PS3 : developmental PS3 OpenMP build with likelihood function speedups."
	echo "  MPI : for a single-threaded build with MPI message passing support."
	echo "  DEBUG: single threaded debug version."
	echo "  GTEST: a unit test version for coverage and regression testing."
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
# NEEDS gcc, g++ 4 and > for fopemp flags on DEV version
	
COMPILER="g++";
COMPILERC="gcc";

CURL_LINKER_LIBS=" -lssl -lcrypto -lcurl";

if [ $TARGET_NAME = "PS3" ]
then
	sysName="AIX"
else
	sysName=`uname`;
fi
echo $sysName;


if [ $sysName == "Darwin" ]
then
	machName=`machine`;
	if [ $machName == "ppc7450" ] 
	then
		COMPILER_FLAGS=" -D __UNIX__ -w -c -fsigned-char -fast -mcpu=7450 -fpermissive -I`pwd`/Source -I`pwd`/Source/SQLite -D INTPTR_TYPE=long "
	else
		COMPILER_FLAGS=" -D __UNIX__ -w -c -fsigned-char -fast -fpermissive -I`pwd`/Source -I`pwd`/Source/SQLite -D INTPTR_TYPE=long "	
	fi
	COMPILER_LINK_FLAGS=" -w -fsigned-char ";
else
	if [ $sysName == "AIX" ]
	then
		if [ $TARGET_NAME = "PS3" ]
		then
			COMPILER="cbexlc++"
			COMPILERC="cbexlc"
			COMPILER_FLAGS=" -c -qchar=signed -O3 -D INTPTR_TYPE=long -D __UNIX__ ";
			COMPILER_LINK_FLAGS="  -qchar=signed ";
		else
			COMPILER="xlC";
			COMPILERC="xlc";
			COMPILER_FLAGS=" -qsmp=omp -c -D _SLKP_LFENGINE_REWRITE_ -qchar=signed -O3 -D INTPTR_TYPE=long -D __UNIX__ -qreport  -qarch=auto -qtune=auto ";
			COMPILER_LINK_FLAGS=" -qsmp=omp -qchar=signed -lxlsmp  -bmaxdata:0x80000000 ";
		fi
	else
		COMPILER_LINK_FLAGS=" -w -fsigned-char ";
		COMPILER_FLAGS=" -w -c -fsigned-char -O3 -fpermissive -I`pwd`/Source -I`pwd`/Source/SQLite -D INTPTR_TYPE=long -D __UNIX__ ";
	fi
fi


# END MODIFY

COMPILER_FLAGS=$COMPILER_FLAGS" -D _SLKP_LFENGINE_REWRITE_ ";

echo "Checking for curl";

rm -rf curl_check*

(echo "#include <curl/curl.h>"; echo "int main(void) {return 0;}") | cat - > curl_check.cpp

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

OBJ_DIR_NAME="obj_$TARGET_NAME"

if [ -f $OBJ_DIR_NAME ] 
then
	rm -rf $OBJ_DIR_NAME;
fi

makedir $OBJ_DIR_NAME

if [ $1 = "SP" ] 
then
	TARGET_NAME="HYPHY";
	LINKER_FLAGS=$CURL_LINKER_LIBS" -lm -ldl -lpthread ";
	echo "+--------------------------------------+"
	echo "|Building a single threaded HYPHYKernel|"
	echo "+--------------------------------------+"
fi

if [ $1 = "MP" ] 
then
	TARGET_NAME="HYPHYMP";
	LINKER_FLAGS=$CURL_LINKER_LIBS" -lm -lpthread -ldl ";
	echo "+---------------------------------------+"
	echo "|Building a multi-threaded HYPHYKernelMP|"
	echo "+---------------------------------------+"
	COMPILER_FLAGS=$COMPILER_FLAGS" -D __MP__"
fi

if [ $1 = "MP2" ] 
then
	TARGET_NAME="HYPHYMP";
	LINKER_FLAGS=$CURL_LINKER_LIBS" -lm -lpthread -fopenmp -ldl ";
	echo "+-----------------------------------------------------------+"
	echo "|Building a multi-threaded HYPHYKernelMP with setconcurrency|"
	echo "+-----------------------------------------------------------+"
	COMPILER_FLAGS=$COMPILER_FLAGS" -D __MP__ -D __MP2__ -fopenmp "
fi

if [ $1 = "PS3" ] 
then
	TARGET_NAME="HYPHY_PS3";
	LINKER_FLAGS=$CURL_LINKER_LIBS" -lm ";
	echo "+-----------------------------------------------------------+"
	echo "|Building a OpenMP/PS3      developmental version of HyPhy  |"
	echo "+-----------------------------------------------------------+"
	COMPILER_FLAGS=$COMPILER_FLAGS" -D __MP__ -D __MP2__ -D _SLKP_LFENGINE_REWRITE_ -qtune=cell -qarch=cell -O5 "
fi

if [ $1 = "MPI" ] 
then
	if [ $sysName == "Darwin" ]
	then
		COMPILER="mpic++";
		COMPILERC="mpicc";
		COMPILER_FLAGS=$COMPILER_FLAGS" -D _SLKP_LFENGINE_REWRITE_ "
		LINKER_FLAGS=$CURL_LINKER_LIBS" -lm -ldl ";
	else
		LINKER_FLAGS=$CURL_LINKER_LIBS" -lpthread -lm -lmpich -ldl ";	
	fi 
	
	if [ $sysName == "AIX" ]
	then
		COMPILER="mpCC";
		COMPILERC="mpcc";
		LINKER_FLAGS=$CURL_LINKER_LIBS" -lm ";
	fi
	
	TARGET_NAME="HYPHYMPI";
	echo "+-----------------------------------------------------------+"
	echo "|Building a single-threaded HYPHYKernelMPI for MPI          |"
	echo "+-----------------------------------------------------------+"
	COMPILER_FLAGS=$COMPILER_FLAGS" -D __HYPHYMPI__ -D _SLKP_LFENGINE_REWRITE_ "
fi

if [ $1 = "DEBUG" ]
then
    TARGET_NAME="HYPHYDebug";
    LINKER_FLAGS=$CURL_LINKER_LIBS" -g -lm -fopenmp ";
    echo "+---------------------------------------+"
    echo "|Building a debug version HYPHYDebug    |"
    echo "+---------------------------------------+"
    COMPILER_FLAGS=" -w -c -g -fsigned-char  -fpermissive -D __UNIX__ -D _SLKP_LFENGINE_REWRITE_ -D INTPTR_TYPE=long -I`pwd`/Source "
fi

if [ $1 = "GTEST" ]
then
    TARGET_NAME="HYPHYTest";
    LINKER_FLAGS=$CURL_LINKER_LIBS" -g -fprofile-arcs -ftest-coverage -lm -fopenmp ../UnitTests/libgtest.a";
    echo "+---------------------------------------+"
    echo "|Building a debug version HYPHYDebug    |"
    echo "+---------------------------------------+"
    COMPILER_FLAGS=" -w -c -g -fprofile-arcs -ftest-coverage -fsigned-char  -fpermissive -D __UNIX__ -D__UNITTEST__ -D _SLKP_LFENGINE_REWRITE_ -D INTPTR_TYPE=long -I`pwd`/Source/  -I`pwd`/UnitTests "

	echo "COMPILER=$COMPILER, $COMPILERC";
	echo "COMPILER_FLAGS=$COMPILER_FLAGS";
	
	cd Source
	for fileName in *.cpp main-unix.cxx SQLite/*c ../UnitTests/*cpp
		do
			name=${fileName%%.[a-z0-9]*}
			name=${name##*/}
			extension=${fileName##*.}
			#echo $name $extension
			obj_file="$name.o";
			#echo Building "$name";
			
			if [ $obj_file -nt $fileName ]
			then
				echo File "$fileName" is up to date
	  		else
				if [ $extension = "c" ]
				then
					cmd="$COMPILERC $COMPILER_FLAGS $fileName"
				else
					cmd="$COMPILER $COMPILER_FLAGS $fileName"
				fi
				echo $cmd
				if `$cmd`
				then
					echo Complete
				else
					echo Error during compilation;
					exit 1;
				fi
			fi
		done
	  
	  
	echo Linking $TARGET_NAME
	echo $COMPILER $COMPILER_LINK_FLAGS -o ../$TARGET_NAME *.o $LINKER_FLAGS
    `$COMPILER $COMPILER_LINK_FLAGS -o ../$TARGET_NAME *.o $LINKER_FLAGS`
 
	exit 0 
fi

if [ $1 = "DMALLOC" ]
then
    TARGET_NAME="HYPHYDmalloc";
    LINKER_FLAGS=$CURL_LINKER_LIBS" -g -lm -ldmalloc ";
    echo "+---------------------------------------+"
    echo "|Building a debug version HYPHYDebug    |"
    echo "+---------------------------------------+"
    COMPILER_FLAGS=" -w -c -g -fsigned-char -D__HYPHYDMALLOC__ -fpermissive -D __UNIX__ -D _SLKP_LFENGINE_REWRITE_ -D INTPTR_TYPE=long ";
fi

if [ $1 = "LIBRARY" ] 
then
	LINKER_FLAGS=$CURL_LINKER_LIBS" -lm -lpthread ";
	echo "+-----------------------------------------------------------+"
	echo "|Building a multi-threaded HYPHY library version            |"
	echo "+-----------------------------------------------------------+"
	if [ $sysName == "Darwin" ]
	then
		TARGET_NAME="libhyphy.so";
		COMPILER_FLAGS=$COMPILER_FLAGS" -fno-strict-aliasing -D __MP__ -D __MP2__ -D __HEADLESS__ -fPIC -I`pwd`/Source/Link "
		COMPILER_LINK_FLAGS=$COMPILER_LINK_FLAGS" -bundle -flat_namespace -undefined suppress "	
	else
		COMPILER_FLAGS=$COMPILER_FLAGS" -D __MP__ -D __MP2__ -D __HEADLESS__ -fPIC -I`pwd`/Source/Link "
		COMPILER_LINK_FLAGS=$COMPILER_LINK_FLAGS" -Wl,-shared "
	fi
	
fi

#COMPILER_FLAGS=$COMPILER_FLAGS" -D __AFYP_REWRITE_BGM__ "

echo "COMPILER=$COMPILER, $COMPILERC";
echo "COMPILER_FLAGS=$COMPILER_FLAGS";


cd Source

if [ $1 = "LIBRARY" ] 
then
	for fileName in *.cpp
	do
	  obj_file=../$OBJ_DIR_NAME/${fileName}.o;
	  if [ $obj_file -nt $fileName ]
	  then
	  	echo File "$fileName" is up to date
	  else
		  echo Building "$fileName";
		  if `$COMPILER -o $obj_file $COMPILER_FLAGS $fileName `
		   then
	      	 echo Complete
		   else
		  		echo Error during compilation;
		  		exit 1;
		   fi
	  fi
	done
cd Link
	for fileName in *.cpp
	do
	  obj_file=../../$OBJ_DIR_NAME/${fileName}.o;
	  if [ $obj_file -nt $fileName ]
	  then
			echo File "$fileName" is up to date
	  else
		  echo Building "$fileName";
		  if `$COMPILER -o $obj_file $COMPILER_FLAGS $fileName `
		   then
			 echo Complete
		   else
				echo Error during compilation;
				exit 1;
		   fi
	  fi
	done
cd ..
# not a library build
else 
	for fileName in *.cpp main-unix.cxx
	do
	  obj_file=../$OBJ_DIR_NAME/${fileName}.o;
	  if [ $obj_file -nt $fileName ]
	  then
	  	echo File "$fileName" is up to date
	  else
		  echo Building "$fileName";
		  if `$COMPILER -o $obj_file $COMPILER_FLAGS $fileName `
		   then
	      	 echo Complete
		   else
		  		echo Error during compilation;
		  		exit 1;
		   fi
	  fi
	done
fi

cd SQLite

for fileName in *.c
do
  obj_file=../../$OBJ_DIR_NAME/${fileName}.o;
  if [ $obj_file -nt $fileName ]
  then 
  	echo SQLite File "$fileName" is up to date
  else
	  echo Building "SQLite file $fileName";
	  if `$COMPILERC -o $obj_file $COMPILER_FLAGS $fileName `
	   then
	                echo Complete
	   else
	                echo Error during compilation;
	                exit 1;
	   fi
  fi
done

cd ../..




if [ $1 != "LIBRARY" -o $LIBRARY_BINDINGS = "NONE" ]
then
	echo Linking $TARGET_NAME
	echo $COMPILER $COMPILER_LINK_FLAGS -o $TARGET_NAME $OBJ_DIR_NAME/*.o  $LINKER_FLAGS
    `$COMPILER $COMPILER_LINK_FLAGS -o $TARGET_NAME $OBJ_DIR_NAME/*.o  $LINKER_FLAGS`
else
	if [ $LIBRARY_BINDINGS = "R" ]
	then
		echo Linking HyPhy.so
		cppf="PKG_CPPFLAGS=\"-I`pwd`/Source/Link/\""
		echo  $cppf R CMD SHLIB -o LibraryModules/R/HyPhy.so $OBJ_DIR_NAME/*.o  LibraryModules/Source/THyPhy_R.cpp $LINKER_FLAGS
		export $cppf; R CMD SHLIB -o LibraryModules/R/HyPhy.so $OBJ_DIR_NAME/*.o  LibraryModules/Source/THyPhy_R.cpp	$LINKER_FLAGS	
		echo R library written to `pwd`/LibraryModules/R/HyPhy.so
	fi
fi
echo Finished


