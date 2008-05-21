cp GUIElements/Shared\ Source/* ../UNIX_Build/GUI/
cp GUIElements/Shared\ Source/*/* ../UNIX_Build/GUI/
cp GUIElements/Shared\ Include/*/* ../UNIX_Build/GUI/
cp GUIElements/Shared\ Include/* ../UNIX_Build/GUI/

cp GUIElements/Platform\ Include/GTK/* ../UNIX_Build/GUI/
cp GUIElements/Platform\ Include/GTK/*/* ../UNIX_Build/GUI/
cp GUIElements/Platform\ Source/GTK/*/* ../UNIX_Build/GUI/
cp GUIElements/Platform\ Source/GTK/* ../UNIX_Build/GUI/

cp Mains/main-GTK.cpp  ../UNIX_Build/Source/main-GTK.cxx
cp Mains/main-unix.cpp ../UNIX_Build/Source/main-unix.cxx
cp Mains/hyphyunixutils.cpp ../UNIX_Build/Source/hyphyunixutils.cpp

cp HeadlessLink/* ../UNIX_Build/Source/Link/
cp Core/preferences* ../UNIX_Build/Source/Link/
cp SWIG/* ../UNIX_Build/Source/Link

cp Core/*.{h,cp,cpp} ../UNIX_Build/Source/
rm -f ../UNIX_Build/Source/preferences.*

cp Core/preferences* ../UNIX_Build/GUI/
cp NewerFunctionality/* ../UNIX_Build/Source/

cp ../SQLLite/*.h ../UNIX_Build/Source/SQLite/
cp ../SQLLite/*.c ../UNIX_Build/Source/SQLite/

