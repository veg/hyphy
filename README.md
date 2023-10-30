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
