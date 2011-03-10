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
	installDirectory=HYPHY
else
	installDirectory=$1
fi

cd ../..

makedir $installDirectory
makedir $installDirectory/Source
makedir $installDirectory/Library
makedir $installDirectory/Source/Link
makedir $installDirectory/GUI
makedir $installDirectory/Source/SQLite
cp trunk/Core/*.{h,cp,cpp} $installDirectory/Source/
rm -f $installDirectory/Source/preferences.cpp
cp trunk/Core/preferences.cpp $installDirectory/GUI
cp trunk/HeadlessLink/*.{h,cpp} $installDirectory/Source/Link

cp trunk/NewerFunctionality/*.{h,cpp} $installDirectory/Source/
cp -R trunk/Library $installDirectory/
cp -R trunk/data $installDirectory/
cp SQLite/trunk/*.{c,h} $installDirectory/Source/SQLite/
cp trunk/Scripts/*.sh $installDirectory/
cp 'trunk/GUIElements/Shared Source/'*.cpp $installDirectory/GUI/
cp 'trunk/GUIElements/Shared Source/Components/'*.cpp $installDirectory/GUI/
cp 'trunk/GUIElements/Shared Include/'*.h $installDirectory/GUI/
cp 'trunk/GUIElements/Shared Include/'{Components,WindowClasses}/*.h $installDirectory/GUI/
cp  'trunk/GUIElements/Platform Include/GTK/'*.h $installDirectory/GUI/
cp  'trunk/GUIElements/Platform Include/GTK/Components/'*.h $installDirectory/GUI/
cp  'trunk/GUIElements/Platform Source/GTK/'*.cpp $installDirectory/GUI/
cp  'trunk/GUIElements/Platform Source/GTK/Components/'*.cpp $installDirectory/GUI/
cp  'trunk/GUIElements/Platform Source/GTK/WindowClasses/'*.cpp $installDirectory/GUI/

cp 	trunk/Mains/main-GTK.cpp  $installDirectory/Source/main-GTK.cxx
cp  trunk/Mains/main-unix.cpp $installDirectory/Source/main-unix.cxx
cp  trunk/Mains/hyphyunixutils.cpp $installDirectory/Source/hyphyunixutils.cpp
cp  -R trunk/{ChartAddIns,DatapanelAddIns,GeneticCodes,Help,SubstitutionClasses,SubstitutionModels,TemplateBatchFiles,TopologyInference,TreeAddIns,UserAddins} $installDirectory
makedir $installDirectory/GTKResources
cp  trunk/GUIElements/Resources/GTKResources/*.* $installDirectory/GTKResources/
makedir $installDirectory/GTKResources/theme
cp  trunk/GUIElements/Resources/GTKResources/theme/*.* $installDirectory/GTKResources/theme

cd $installDirectory
bash build.sh MP2
#bash build.sh DEV
bash gtk_build.sh MP2
