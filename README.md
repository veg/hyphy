[![Build Status](https://travis-ci.org/veg/hyphy.svg?branch=master)](https://travis-ci.org/veg/hyphy)

# HyPhy - Hypothesis testing using Phylogenies

HyPhy is an open-source software package for the analysis of genetic sequences using techniques in phylogenetics, molecular evolution, and machine learning. It features a complete graphical user interface (GUI) and a rich scripting language for limitless customization of analyses. Additionally, HyPhy features support for parallel computing environments (via message passing interface (MPI)) and it can be compiled as a shared library and called from other programming environments such as Python and R. HyPhy is the computational backbone powering datamonkey.org. Additional information is available at hyphy.org.

## Quick Start

#### Install  
`conda install hyphy`

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
+   MPI - build a HyPhy executable (HYPHYMPI) using MPI to do multiprocessing
+   HYPHYMPI - build a HyPhy executable (HYPHYMPI) using openMPI 
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

