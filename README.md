[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/hyphy/README.html)

# HyPhy - Hypothesis testing using Phylogenies

HyPhy is an open-source software package for the analysis of genetic sequences using techniques in phylogenetics, molecular evolution, and machine learning. It features a complete graphical user interface (GUI) and a rich scripting language for limitless customization of analyses. Additionally, HyPhy features support for parallel computing environments (via message passing interface (MPI)) and it can be compiled as a shared library and called from other programming environments such as Python and R. HyPhy is the computational backbone powering datamonkey.org. Additional information is available at hyphy.org.

## Quick Start

#### Install
`conda install -c bioconda hyphy`

#### Running with Docker
You can also run HyPhy without having to install it on your system using the provided Dockerfile. Following the below instructions starts an interactive Docker container where HyPhy is already available. 

If you don't have Docker installed, first look at [Docker Desktop](https://docs.docker.com/desktop/).

**Please note you must change this code snippet to point to the appropriate location for your input data. This will be made available inside the container at `/hyphy/data`.**

```
git clone https://github.com/veg/hyphy.git
cd hyphy
docker build -t hyphy:latest .
docker run --rm -v [path-to-your-input-data]:/hyphy/data -it hyphy:latest
```

#### Run with Command Line Arguments
`hyphy <method_name> --alignment <path_to_alignment_file> <additional_method_specific_arguments>`  
+ _`<method_name>` is the name of the analysis you wish to run (can be: absrel, bgm, busted, fade, fel, fubar, gard, meme, relax or slac)_
+ _`<path_to_alignment_file>` is the relative or absolute path to a fasta, nexus or phylib file containing an alignment and tree_
+ _A list of the available `<additional_method_specific_arguments>` can be seen by running `hyphy <method_name> --help`_

or  

#### Run in Interactive Mode
`hyphy -i`  

## Building from Source

#### Requirements
* cmake >= 3.12
* gcc >= 4.9
* libcurl
* libpthread
* openmp (can be installed on mac via `brew install libomp`)

#### Download
You can download a specific release [here](https://github.com/veg/hyphy/releases) or clone this repo with

`git clone https://github.com/veg/hyphy.git`

Change your directory to the downloaded/cloned directory

`cd hyphy`

#### Build
`cmake .`

`make install`

## Additional Options for Building from Source

#### Build Systems
If you prefer to use other build systems, such as Xcode, configure using the -G switch

`cmake -G Xcode .`

CMake supports a number of build system generators, feel free to peruse these and use them if you wish.

If you are on an OS X platform, you can specify which OS X SDK to use

`cmake -DCMAKE_OSX_SYSROOT=/Developer/SDKs/MacOSX10.9.sdk/ .`

If building on a heterogeneous cluster with some nodes that do not support auto-vectorization  

`cmake -DNOAVX=ON .`.

If you're on a UNIX-compatible system, and you're comfortable with GNU make, then run `make` with one of the following build targets:

+   MP or hyphy - build a HyPhy executable (This used to be "HYPHYMP" but is now just "hyphy") using pthreads to do multiprocessing
+   MPI or HYPHYMPI - build a HyPhy executable (HYPHYMPI) using openMPI to do multiprocessing
+   LIB - build a HyPhy library (libhyphy_mp) using pthreads to do multiprocessing
-   GTEST - build HyPhy's gtest testing executable (HYPHYGTEST)

#### Example (MPI build of hyphy using openMPI)
Ensure that you have openmpi installed and available on your  path. You can check if this is the case after running `cmake .` you should see something similar to this in your output

`-- Found MPI_C: /opt/scyld/openmpi/1.6.3/gnu/lib/libmpi.so;/usr/lib64/libibverbs.so;/usr/lib64/libdat.so;/usr/lib64/librt.so;/usr/lib64/libnsl.so;/usr/lib64/libutil.so;/usr/lib64/libm.so;/usr/lib64/libtorque.so;/usr/lib64/libm.so;/usr/lib64/libnuma.so;/usr/lib64/librt.so;/usr/lib64/libnsl.so;/usr/lib64/libutil.so;/usr/lib64/libm.so `

`-- Found MPI_CXX: /opt/scyld/openmpi/1.6.3/gnu/lib/libmpi_cxx.so;/opt/scyld/openmpi/1.6.3/gnu/lib/libmpi.so;/usr/lib64/libibverbs.so;/usr/lib64/libdat.so;/usr/lib64/librt.so;/usr/lib64/libnsl.so;/usr/lib64/libutil.so;/usr/lib64/libm.so;/usr/lib64/libtorque.so;/usr/lib64/libm.so;/usr/lib64/libnuma.so;/usr/lib64/librt.so;/usr/lib64/libnsl.so;/usr/lib64/libutil.so;/usr/lib64/libm.so `

Then run 

`make HYPHYMPI`

And then run make install to install the software

`make install`

+   hyphy will be installed at  `/location/of/choice/bin`
+   libhyphy_mp.(so/dylib/dll) will be installed at `/location/of/choice/lib`
+   HyPhy's standard library of batchfiles will go into `/location/of/choice/lib/hyphy`

#### WebAssembly Build

HyPhy provides WebAssembly builds that can run in modern browsers without installation.

##### Using Pre-built WebAssembly Files

The latest WebAssembly build can be downloaded from the GitHub Actions artifacts. Go to the [GitHub Actions page](https://github.com/veg/hyphy/actions/workflows/wasm-build.yml), select the most recent successful workflow run, and download the `hyphy-wasm` artifact. The artifact contains:

```
hyphy.js
hyphy.wasm
hyphy.data
hyphy_resources.tar.gz
```

Extract the resources tarball (`hyphy_resources.tar.gz`) in the same directory as the other files. All files should be served from the same directory on your web server. For an example of HyPhy running in the browser, see https://observablehq.com/@spond/hyphy-biowasm.

##### Building WebAssembly Files Locally

To build the WebAssembly files yourself, you'll need [Emscripten](https://emscripten.org/docs/compiling/Building-Projects.html).

1. Configure the build - note that we don't use preload flags during configuration to avoid CMake test issues:
```
emcmake cmake -DCMAKE_EXE_LINKER_FLAGS="-sTOTAL_STACK=2097152 -O2 -sASSERTIONS=1 -sMODULARIZE=1 -sALLOW_MEMORY_GROWTH -sFORCE_FILESYSTEM=1 -sEXIT_RUNTIME=0 -s EXPORTED_RUNTIME_METHODS=[\"callMain\",\"FS\",\"PROXYFS\",\"WORKERFS\",\"UTF8ToString\",\"getValue\",\"AsciiToString\"] -lworkerfs.js -lproxyfs.js -s INVOKE_RUN=0 -s ENVIRONMENT=\"web,worker\" ${EM_FLAGS//-s /-s} -fwasm-exceptions"
``` 

2. Build the target:
```
emmake make -j hyphy 
```

3. This will create the following files:
```
hyphy.js
hyphy.wasm
hyphy.data
```

4. You'll also need to copy the resource files to serve alongside the WebAssembly files:
```
mkdir -p hyphy_resources
cp -r res hyphy_resources/
cp -r tests/hbltests hyphy_resources/tests
```

5. All files should be served from the same directory on your web server.

#### Testing

Use `make test` after running `cmake .`.

#### Benchmarks for CMake Tests

Benchmarks, using Github Actions, can be found at http://hyphy.org/bench

#### Executable Location

By default, HyPhy installs into `/usr/local` but it can be installed on any location of your system by providing an installation prefix

`cmake -DCMAKE_INSTALL_PREFIX:PATH=/location/of/choice`

For example, this configuration will install hyphy at /opt/hyphy

`mkdir -p /opt/hyphy`

`cmake -DCMAKE_INSTALL_PREFIX:PATH=/opt/hyphy .`

#### Building Documentation

```
make docs
cd docs
python3 -m http.server
```

## CLI notes

As noted in the documentation [here](http://www.hyphy.org/tutorials/CLI-tutorial/) 
`hyphy` can be run as a command line tool.
Indeed for many analyses the `hyphy` CLI will return useful 
help messages, showing which parameter values can be set to 
specify your analysis. For example, running `hyphy gard --help`

```{bash}
hyphy gard --help 

Available analysis command line options
---------------------------------------
Use --option VALUE syntax to invoke
If a [reqired] option is not provided on the command line, the analysis will prompt for its value
[conditionally required] options may or not be required based on the values of other options

type
	The type of data to perform screening on
	default value: nucleotide

code
	Genetic code to use (for codon alignments)
	default value: Universal
	applies to: Choose Genetic Code

alignment [required]
	Sequence alignment to screen for recombination

model
	The substitution model to use
	default value: JTT

rv
	Site to site rate variation
	default value: None

max-breakpoints
	Maximum number of breakpoints to consider
	default value: 10000

rate-classes
	How many site rate classes to use
	default value: 4

output
	Write the resulting JSON to this file (default is to save to the same path as the alignment file + 'GARD.json')
	default value: gard.defaultJsonFilePath [computed at run time]

mode
	Run mode (Normal or Faster)
	default value: Normal

output-lf
	Write the best fitting HyPhy analysis snapshot to (default is to save to the same path as the alignment file + 'best-gard')
	default value: gard.defaultFitFilePath [computed at run time]
```

will show you the options that can be set for the `gard` analysis.
So for instance one could specify a `gard` run on the command line
with the following command

```{bash}
hyphy gard --alignment /path/to/file --rv GDD --mode Faster --rate-classes 3
```

While this is a useful feature, it is not always the case that 
older analyses will have the same level of support for command line.
For instance, the `acd` analysis does not have CLI support and 
so if one runs the help command

```{bash}
hyphy acd --help 

Available analysis command line options
---------------------------------------
Use --option VALUE syntax to invoke
If a [reqired] option is not provided on the command line, the analysis will prompt for its value
[conditionally required] options may or not be required based on the values of other options

No annotated keyword arguments are available for this analysis
```

one will see that there are no options available. In this case, 
you can use a different CLI specification. Indeed the CLI will accept
all of the options that are asked for in an interactive session, as 
positional arguments. In this case I could run the `acd` analysis with

```{bash}
hyphy acd Universal <alignment.fa> MG94CUSTOMCF3X4 Global 012345 <treefile> Estimate

```

where the options are specified in the exact order that they are asked for in the interactive session. This will work for all `hyphy` analyses
and provides a less readable but more flexible way to run `hyphy` analyses.
