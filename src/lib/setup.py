#!/usr/bin/python

from distutils.core      import setup,Extension
from distutils.sysconfig import get_python_inc
from os                  import listdir,getcwd,path
from sys                 import maxsize
#incdir = get_python_inc(plat_specific=1)
#print incdir

#build the list of Source files

global sourceFiles, currentWDir

dirFiles    = ["../core", "../new", "../../contrib/SQLite-3.6.17/", "Link"]
sourceFiles = []
currentWDir = getcwd()
includeDirs = dirFiles[-2:]
includeDirs.append("../gui/include")

for aDir in dirFiles:
    sourceFiles.extend([path.normpath(path.join(currentWDir,aDir,aPath)) for aPath in listdir (aDir) if (aPath.endswith ('cpp') or aPath.endswith ('c'))])
    includeDirs.extend([path.normpath(path.join(currentWDir,aDir,aPath)) for aPath in listdir(aDir) if aPath == 'include'])

#sourceFiles.append (path.normpath(path.join(currentWDir,"../Mains/hyphyunixutils.cpp")))
sourceFiles.append (path.normpath(path.join(currentWDir,"SWIGWrappers/THyPhy_python.cpp")))
sourceFiles.append (path.normpath(path.join(currentWDir,"../gui/preferences.cpp")))

definitions = [
    ('SQLITE_PTR_SIZE','sizeof(long)'),
    ('__HYPHY_64__',None),
    ('__UNIX__',None),
    ('__MP__',None),
    ('__MP2__',None),
    ('_SLKP_LFENGINE_REWRITE_',None),
    ('__HEADLESS__',None)
]

# python.org recommended test for 64bit-edness
is_64bits = maxsize > 2**32
if is_64bits:
    definitions.append(('__HYPHY_64__',None))

setup(name='HyPhy',
      version='0.1',
      description   = 'HyPhy package interface library',
      author        = 'Sergei L Kosakovsky Pond',
      author_email  = 'spond@ucsd.edu',
      url           = 'http://www.hyphy.org/',
      package_dir   = {'': 'LibraryModules/Python'},
      packages      = ['HyPhy'],
      py_modules    = ['HyPhy.py'],
      ext_modules   = [Extension('_HyPhy',
            sourceFiles,
            include_dirs = includeDirs,
            define_macros = definitions,
            libraries = ['pthread','ssl','crypto','curl'],
            extra_compile_args = ['-w', '-fsigned-char', '-O3', '-fPIC']
          )]
     )
