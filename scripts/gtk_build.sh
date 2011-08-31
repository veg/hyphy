#!/bin/sh

# MODIFY THESE BASED ON YOUR SYSTEM
# DEFAULT SETTINGS ARE FOR GCC
	
CURL_LINKER_LIBS=" -lssl -lcrypto -lcurl";

COMPILER="g++";
COMPILERC="gcc";
sysName=`uname`;
echo $sysName;

if [ $sysName == "Darwin" ]
then
	machName=`machine`;
	if [ $machName == "ppc7450" ] 
	then
		COMPILER_FLAGS=" -w -c -fsigned-char -fast -mcpu=7450 -fpermissive -I`pwd`/GUI -I`pwd`/Source -I`pwd`/Source/SQLite -D INTPTR_TYPE=long `pkg-config gtk+-2.0 --cflags`  -D GDK_PIXBUF_ENABLE_BACKEND -D __HYPHY_GTK__ "
	else
		COMPILER_FLAGS=" -w -c -fsigned-char -fast -fpermissive -I`pwd`/GUI -I`pwd`/Source -I`pwd`/Source/SQLite -D INTPTR_TYPE=long `pkg-config gtk+-2.0 --cflags`  -D GDK_PIXBUF_ENABLE_BACKEND  -D __HYPHY_GTK__ "
	fi
	COMPILER_LINK_FLAGS=" -w -fsigned-char -D __HYPHY_GTK__ -D GDK_PIXBUF_ENABLE_BACKEND ";
else
	COMPILER_LINK_FLAGS=" -w -fsigned-char ";
	COMPILER_FLAGS=" -w -c -O3 -D INTPTR_TYPE=long -fsigned-char -fpermissive -I`pwd`/GUI -I`pwd`/Source -I`pwd`/Source/SQLite `pkg-config gtk+-2.0 --cflags`  -D GDK_PIXBUF_ENABLE_BACKEND -D __HYPHY_GTK__ -D _SLKP_LFENGINE_REWRITE_ ";
fi

echo "Checking for curl";
echo -e "#include <curl/curl.h>\nint main(void) {return 0;}" > curl_check.cpp

if `$COMPILER -o curl_check -w $CURL_LINKER_LIBS curl_check.cpp`
then
 echo "Curl seems to be present"
else
	echo "Curl seems to be absent (setting up compiler options skip CURL code)";
	CURL_LINKER_LIBS="";
	COMPILER_FLAGS=$COMPILER_FLAGS" -D__HYPHY_NO_CURL__";
	COMPILER_LINK_FLAGS=$COMPILER_LINK_FLAGS" -D__HYPHY_NO_CURL__";
fi

rm -rf curl_check.*

# END MODIFY


TARGET_NAME="BLANK"

if [ $# -ne 1 ]
then
	TARGET_NAME="HELP";
else
	if [ $1 != "SP" -a $1 != "MP" -a $1 != "MP2"  -a $1 != "MPI" -a $1 != "DEBUG" ] 
	then
		TARGET_NAME="HELP"
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
	echo "  MPI : for a single-threaded build with MPI message passing support."
	echo "  DEBUG: single threaded debug version."
	exit 1
fi
	

OBJ_DIR_NAME="objGTK_$TARGET_NAME"

if [ -f $OBJ_DIR_NAME ] 
then
	rm -rf $OBJ_DIR_NAME;
fi

if [ -f $OBJ_DIR_NAME ] 
then
	echo "Insufficient permissions to create an object directory";
	exit 1;
fi

if [ ! -d $OBJ_DIR_NAME  ]
then
	if [ `mkdir $OBJ_DIR_NAME` ]
	then
		echo "Failed to create object directory";
		exit 1;
	fi
fi

if [ $1 = "SP" ] 
then
	TARGET_NAME="HYPHY_GTK";
	LINKER_FLAGS=$CURL_LINKER_LIBS" -lm ";
	echo "+--------------------------------------+"
	echo "|Building a single threaded HYPHYKernel|"
	echo "+--------------------------------------+"
fi

if [ $1 = "MP" ] 
then
	TARGET_NAME="HYPHYMP_GTK";
	LINKER_FLAGS=$CURL_LINKER_LIBS" -lm -lpthread ";
	echo "+---------------------------------------+"
	echo "|Building a multi-threaded HYPHYKernelMP|"
	echo "+---------------------------------------+"
	COMPILER_FLAGS=$COMPILER_FLAGS" -D __MP__ "
fi

if [ $1 = "MP2" ] 
then
	TARGET_NAME="HYPHYMP_GTK";
	LINKER_FLAGS=$CURL_LINKER_LIBS" -lm -lpthread -fopenmp ";
	echo "+-----------------------------------------------------------+"
	echo "|Building a multi-threaded HYPHYKernelMP with setconcurrency|"
	echo "+-----------------------------------------------------------+"
	COMPILER_FLAGS=$COMPILER_FLAGS" -D __MP__ -D __MP2__ -fopenmp "
fi

if [ $1 = "MPI" ] 
then
	TARGET_NAME="HYPHYMPI_GTK";
	LINKER_FLAGS=$CURL_LINKER_LIBS" -lm -lmpich ";
	echo "+-----------------------------------------------------------+"
	echo "|Building a single-threaded HYPHYKernelMPI for MPI          |"
	echo "+-----------------------------------------------------------+"
	COMPILER_FLAGS=$COMPILER_FLAGS" -D __HYPHYMPI__ "
fi

if [ $1 = "DEBUG" ]
then
        TARGET_NAME="HYPHYDebug";
        LINKER_FLAGS=" -lssl -lcrypto -lcurl -lm";
        echo "+---------------------------------------+"
        echo "|Building a debug version HYPHYDebug    |"
        echo "+---------------------------------------+"
        (echo "#define __HYPHY_GTK__") | cat > Source/platform.h
        COMPILER_FLAGS=" -w -c -g -fpermissive ";
fi

echo "COMPILER=$COMPILER, $COMPILERC";
echo "COMPILER_FLAGS=$COMPILER_FLAGS";

cd Source

for fileName in *.cpp main*GTK.cxx
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

cd "../../GUI/"

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

cd ..

echo Linking $TARGET_NAME
GTK_LIBS=`pkg-config gtk+-2.0 --libs`;
echo $GTK_LIBS $LINKER_FLAGS
`$COMPILER $COMPILER_LINK_FLAGS -o $TARGET_NAME $OBJ_DIR_NAME/*.o  $GTK_LIBS $LINKER_FLAGS`
echo Finished



