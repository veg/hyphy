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

if [ $# -ne 1 ]
then
	installDirectory=HYPHY_MinGW
else
	installDirectory=$1
fi

cd ../..

makedir $installDirectory
makedir $installDirectory/Source
makedir $installDirectory/Win32GUI
makedir $installDirectory/Source/SQLite
cp trunk/Core/*.{h,cp,cpp} $installDirectory/Source/
rm -f $installDirectory/Source/preferences.cpp
cp trunk/Core/preferences.cpp $installDirectory/Win32GUI

cp trunk/NewerFunctionality/*.{h,cpp} $installDirectory/Source/
cp SQLite/trunk/*.{c,h} $installDirectory/Source/SQLite/
cp trunk/Scripts/*.sh $installDirectory/

cp 'trunk/GUIElements/Shared Source/'*.cpp $installDirectory/Win32GUI/
cp 'trunk/GUIElements/Shared Source/Components/'*.cpp $installDirectory/Win32GUI/
cp 'trunk/GUIElements/Shared Include/'*.h $installDirectory/Win32GUI/
cp 'trunk/GUIElements/Shared Include/'{Components,WindowClasses}/*.h $installDirectory/Win32GUI/
cp  'trunk/GUIElements/Platform Include/Windows/'*.h $installDirectory/Win32GUI/
cp  'trunk/GUIElements/Platform Include/Windows/Components/'*.h $installDirectory/Win32GUI/
cp  'trunk/GUIElements/Platform Source/Windows/'*.cpp $installDirectory/Win32GUI/
cp  'trunk/GUIElements/Platform Source/Windows/Components/'*.cpp $installDirectory/Win32GUI/
cp  'trunk/GUIElements/Platform Source/Windows/WindowClasses/'*.cpp $installDirectory/Win32GUI/

cp  trunk/Mains/hyphyunixutils.cpp $installDirectory/Source/hyphyunixutils.cpp
cp  trunk/Mains/main-win.cpp $installDirectory/Source/main-win.cxx
cp  trunk/Mains/main-unix.cpp $installDirectory/Source/main-unix.cxx
cp  -R trunk/{ChartAddIns,DatapanelAddIns,GeneticCodes,Help,SubstitutionClasses,SubstitutionModels,TemplateBatchFiles,TopologyInference,TreeAddIns,UserAddIns} $installDirectory
cp  -R 'trunk/GUIElements/Resources/Windows' $installDirectory/Win32GUI/
rm -f $installDirectory/Win32GUI/Windows/*.{dll,lib}
cp  'trunk/GUIElements/Resources/Windows/pthreadGC2.dll' $installDirectory

cd $installDirectory
sh build_mingw.sh
sh build_mingw_gui.sh
