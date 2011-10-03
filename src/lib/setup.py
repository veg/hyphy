#!/usr/bin/python

from distutils.core      import setup, Extension
from distutils.sysconfig import get_python_inc
from os                  import listdir, getcwd, path
from glob                import glob
import sys
#incdir = get_python_inc(plat_specific=1)
#print incdir

#build the list of Source files

scriptPath = path.realpath(path.dirname(sys.argv[0]))
srcPath, libDir = path.split(scriptPath)
hyphyPath, srcDir = path.split(srcPath)

contribPath = path.join(hyphyPath, 'contrib')
sqlitePath = path.join(contribPath, 'SQLite-3.6.17')

linkPath = path.join(scriptPath, 'Link')
coreSrcPath = path.join(srcPath, 'core')
newSrcPath = path.join(srcPath, 'new')
guiSrcPath = path.join(srcPath, 'gui')
prefFile = [path.join(guiSrcPath, 'preferences.cpp')]
swigFile = [path.join(scriptPath, 'SWIGWrappers', 'THyPhy_python.cpp')]

coreSrcFiles = glob(path.join(coreSrcPath, '*.cpp'))
newSrcFiles = glob(path.join(newSrcPath, '*.cpp'))
sqliteFiles = glob(path.join(sqlitePath, '*.c'))
linkFiles = glob(path.join(linkPath, '*.cpp')) # + glob(path.join(linkPath, '*.cxx'))
utilFiles = glob(path.join(srcPath, 'utils', '*.cpp'))

sourceFiles = coreSrcFiles + newSrcFiles +  sqliteFiles + prefFile + linkFiles + swigFile + utilFiles

includePaths =  [path.join(p, 'include') for p in [coreSrcPath, newSrcPath, guiSrcPath]]
includePaths += [linkPath, contribPath]


setup(
    name = 'HyPhy',
    version = '0.1',
    description = 'HyPhy package interface library',
    author = 'Sergei L Kosakovsky Pond',
    author_email = 'spond@ucsd.edu',
    url = 'http://www.hyphy.org/',
    package_dir = {'': 'LibraryModules/Python'},
    packages = ['HyPhy'],
    # py_modules = ['HyPhy'],
    ext_modules = [Extension('_HyPhy',
            sourceFiles,
            include_dirs = includePaths,
            define_macros = [('SQLITE_PTR_SIZE','sizeof(long)'),
                             ('__UNIX__', None),
                             ('__MP__', None),
                             ('__MP2__', None),
                             ('_SLKP_LFENGINE_REWRITE_', None),
                             ('__HEADLESS__', None)],
            libraries = ['pthread', 'ssl', 'crypto', 'curl'],
            extra_compile_args = [
#                     '-Wno-int-to-pointer-cast',
#                     '-Wno-pointer-to-int-cast',
                    '-Wno-char-subscripts',
                    '-Wno-sign-compare',
                    '-fsigned-char',
                    '-O3',
                    '-fpermissive',
                    '-fPIC'
            ]
    )]
)
