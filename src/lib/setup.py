#!/usr/bin/python

from distutils.core      import setup, Extension
from os                  import path
from glob                import glob
import sys

from platform import architecture, mac_ver

#incdir = get_python_inc(plat_specific=1)
#print incdir


#build the list of Source files

scriptPath = path.realpath(path.dirname(sys.argv[0]))
srcPath, libDir = path.split(scriptPath)
hyphyPath, srcDir = path.split(srcPath)
# with open('batchfiles.list') as fh:
#     resFiles = [(f, path.join(*(['..'] * 5 + f.split('/')))) for f in fh.read().split('\n') if f != '']

contribPath = path.join(hyphyPath, 'contrib')
sqlitePath = path.join(contribPath, 'SQLite-3.8.2')

linkPath = path.join(scriptPath, 'Link')
coreSrcPath = path.join(srcPath, 'core')
newSrcPath = path.join(srcPath, 'new')
guiSrcPath = path.join(srcPath, 'gui')
prefFile = [path.join(guiSrcPath, 'preferences.cpp')]

if sys.version_info >= (3,0,0):
    swigFile = [path.join(scriptPath, 'SWIGWrappers', 'THyPhy_py3.cpp')]
else:
    swigFile = [path.join(scriptPath, 'SWIGWrappers', 'THyPhy_python.cpp')]

coreSrcFiles = glob(path.join(coreSrcPath, '*.cpp'))
newSrcFiles = glob(path.join(newSrcPath, '*.cpp'))
sqliteFiles = glob(path.join(sqlitePath, '*.c'))
linkFiles = glob(path.join(linkPath, '*.cpp')) # + glob(path.join(linkPath, '*.cxx'))
utilFiles = glob(path.join(srcPath, 'utils', '*.cpp'))

sourceFiles = coreSrcFiles + newSrcFiles +  sqliteFiles + prefFile + linkFiles + swigFile + utilFiles

includePaths =  [path.join(p, 'include') for p in [coreSrcPath, newSrcPath, guiSrcPath]]
includePaths += [linkPath, contribPath, sqlitePath]

# check for 64bit and define as such
define_macros = [('__HYPHY_64__', None)] if '64' in architecture()[0] else []

# openmp on Mac OS X Lion is broken
major, minor = mac_ver()[0].split('.')
openmp = ['-fopenmp'] if int(major) < 10 or (int(major) == 10 and int(minor) < 7) else []

setup(
    name = 'HyPhy',
    version = '0.1.1',
    description = 'HyPhy package interface library',
    author = 'Sergei L Kosakovsky Pond',
    author_email = 'spond@temple.edu',
    url = 'http://www.hyphy.org/',
    packages = ['HyPhy'],
    package_dir = {'HyPhy': 'LibraryModules/Python/HyPhy'},
#    data_files = resFiles,
    # py_modules = ['HyPhy'],
    ext_modules = [Extension('_HyPhy',
            sourceFiles,
            include_dirs = includePaths,
            define_macros = [('SQLITE_PTR_SIZE','sizeof(long)'),
                             ('__UNIX__', None),
                             ('__MP__', None),
                             ('__MP2__', None),
                             ('_SLKP_LFENGINE_REWRITE_', None),
                             ('__AFYP_REWRITE_BGM__', None),
                             ('__HEADLESS__', None),
                             ('_HYPHY_LIBDIRECTORY_', '"/usr/local/lib/hyphy"')] + define_macros,
            libraries = ['pthread', 'ssl', 'crypto', 'curl'],
            extra_compile_args = [
                    '-Wno-int-to-pointer-cast',
                    # '-Wno-pointer-to-int-cast',
                    '-Wno-char-subscripts',
                    '-Wno-sign-compare',
                    '-Wno-parentheses',
                    '-Wno-uninitialized',
#                    '-Wno-conversion-null',
                    '-Wno-unused-variable',
#                    '-Wno-unused-but-set-variable',
                    '-Wno-shorten-64-to-32',
                    '-fsigned-char',
                    '-O3',
                    '-fpermissive',
                    '-fPIC',
            ] + openmp,
            extra_link_args = [
            ] + openmp
    )]
)
