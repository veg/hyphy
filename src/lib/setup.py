#!/usr/bin/python

from setuptools          import Extension, setup, find_packages
from os                  import path
from glob                import glob
import sys

from setuptools.command.build_py import build_py



#incdir = get_python_inc(plat_specific=1)
#print incdir

#build the list of Source files

scriptPath = path.realpath(path.dirname(sys.argv[0]))
srcPath, libDir = path.split(scriptPath)
hyphyPath, srcDir = path.split(srcPath)


# with open('batchfiles.list') as fh:
#     resFiles = [(f, path.join(*(['..'] * 5 + f.split('/')))) for f in fh.read().split('\n') if f != '']

contribPath = path.join(hyphyPath, 'contrib')

linkPath        = path.join(scriptPath)
coreSrcPath     = path.join(srcPath, 'core')
newSrcPath      = path.join(srcPath, 'new')
contribSrcPath  = path.join(srcPath, 'contrib')

swigFile = [path.join(scriptPath,  'THyPhy.i')]

coreSrcFiles = glob(path.join(coreSrcPath, '*.cpp'))
newSrcFiles = glob(path.join(newSrcPath, '*.cpp'))
sqliteFiles = [] # glob(path.join(sqlitePath, 'sqlite3.c'))
linkFiles = glob(path.join(linkPath, '*.cxx')) # + glob(path.join(linkPath, '*.cxx'))
utilFiles = glob(path.join(srcPath, 'utils', '*.cpp'))
contribFiles = glob(path.join(contribSrcPath, '*.cpp'))

sourceFiles = coreSrcFiles + newSrcFiles +  contribFiles  + linkFiles + swigFile + utilFiles

includePaths =  [path.join(p, 'include') for p in [coreSrcPath, newSrcPath]]
includePaths += [linkPath, contribSrcPath]

# Build extensions before python modules,
# or the generated SWIG python files will be missing.
class BuildPy(build_py):
    def run(self):
        self.run_command('build_ext')
        super(build_py, self).run()
        

setup(
    name = 'HyPhy',
    version = '0.2.0',
    description = 'HyPhy package interface library for Python',
    author = 'Sergei L Kosakovsky Pond',
    author_email = 'spond@temple.edu',
    url = 'http://www.hyphy.org/',
    py_modules=["HyPhy"],
    cmdclass={
                    'build_py': BuildPy,
              },
              
    package_dir={'': '.'},
    python_requires='>=3.4',
    ext_modules = [Extension(name='_HyPhy',
            sources=sourceFiles,
            include_dirs = includePaths,
            define_macros = [('__HYPHY_NO_SQLITE__',None),
                             ('__UNIX__', None),
                             ('__MP__', None),
                             ('__MP2__', None),
                             ('_SLKP_LFENGINE_REWRITE_', None),
                             ('__AFYP_REWRITE_BGM__', None),
                             ('__HEADLESS__', None),
                             ('_HYPHY_LIBDIRECTORY_', '"/usr/local/lib/hyphy"')] ,
            extra_compile_args = [
                    '-Wno-int-to-pointer-cast',
                    # '-Wno-pointer-to-int-cast',
                    '-Wno-char-subscripts',
                    '-Wno-sign-compare',
                    '-Wno-parentheses',
                    '-Wno-uninitialized',
                    '-Wno-unused-variable',
                    '-Wno-shorten-64-to-32',
                    '-fsigned-char',
                    '-O3',
                    '-fpermissive',
                    '-fPIC',
                    '-std=c++14'
            ],
            swig_opts = ['-c++']
    )]

)
