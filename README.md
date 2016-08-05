[![HyPhy.org](http://hyphy.org/w/skins/common/images/HyphyLogo_135x135.png | width=50)](http://hyphy.org) HyPhy - Hypothesis testing using Phylogenies
============================================

[![Build Status](https://travis-ci.org/veg/hyphy.svg)](https://travis-ci.org/veg/hyphy)

Introduction
------------
HyPhy is an open-source software package for the analysis of genetic sequences using techniques in phylogenetics, molecular evolution, and machine learning. It features a complete graphical user interface (GUI) and a rich scripting language for limitless customization of analyses. Additionally, HyPhy features support for parallel computing environments (via message passing interface (MPI)) and it can be compiled as a shared library and called from other programming environments such as Python and R. 

Installation
------------

Hyphy depends on CMake for its build system.

To install, make sure you have CMake >3.0 installed. Hyphy is dependent on other development libraries like
libcurl and libpthread. Libcurl requires development libraries such as crypto++ and openssl ( or gnutls depending on your configuration)
On Ubuntu these are libcurl-dev, libcrypto++-dev and libssl-dev.

You can download a specific release from this repository or if you prefer the master branch simply 
clone the repo with

`git clone git@github.com:veg/hyphy.git`

Change your directory to the newly cloned directory

`cd hyphy`

Configure the project from the source directory using CMake.

`cmake .`

If you prefer to use other build systems, such as Xcode,
configure using the -G switch

`cmake -G Xcode .`

CMake supports a number of build system generators,
feel free to peruse these and use them if you wish.

By default, HyPhy installs into `/usr/local`
but it can be installed on any location of your system 
by providing an installation prefix

`cmake -DINSTALL_PREFIX=/location/of/choice`

For example, this configuration will install hyphy at /opt/hyphy

`mkdir -p /opt/hyphy`

`cmake -DINSTALL_PREFIX=/opt/hyphy .`

If you are on an OS X platform, you can specify which OS X SDK to use

`cmake -DCMAKE_OSX_SYSROOT=/Developer/SDKs/MacOSX10.9.sdk/ .`

If you're on a UNIX-compatible system,
and you're comfortable with GNU make,
then run `make` with one of the following build targets:

+   MAC - build a Mac Carbon application
+   HYPHYGTK - HYPHY with GTK
+   SP - build a HyPhy executable (HYPHYSP) without multiprocessing
+   MP2 - build a HyPhy executable (HYPHYMP) using pthreads to do multiprocessing
+   MPI - build a HyPhy executable (HYPHYMPI) using MPI to do multiprocessing
+   HYPHYMPI - build a HyPhy executable (HYPHYMPI) using openMPI 
+   LIB - build a HyPhy library (libhyphy_mp) using pthreads to do multiprocessing
-   GTEST - build HyPhy's gtest testing executable (HYPHYGTEST)

For example to create a MPI build of HYPHY using openMPI ensure that you 
have openmpi installed and available on your  path. You can check if this
is the case after running 
`cmake .` you should see something similar to this in your output

`-- Found MPI_C: /opt/scyld/openmpi/1.6.3/gnu/lib/libmpi.so;/usr/lib64/libibverbs.so;/usr/lib64/libdat.so;/usr/lib64/librt.so;/usr/lib64/libnsl.so;/usr/lib64/libutil.so;/usr/lib64/libm.so;/usr/lib64/libtorque.so;/usr/lib64/libm.so;/usr/lib64/libnuma.so;/usr/lib64/librt.so;/usr/lib64/libnsl.so;/usr/lib64/libutil.so;/usr/lib64/libm.so `

`-- Found MPI_CXX: /opt/scyld/openmpi/1.6.3/gnu/lib/libmpi_cxx.so;/opt/scyld/openmpi/1.6.3/gnu/lib/libmpi.so;/usr/lib64/libibverbs.so;/usr/lib64/libdat.so;/usr/lib64/librt.so;/usr/lib64/libnsl.so;/usr/lib64/libutil.so;/usr/lib64/libm.so;/usr/lib64/libtorque.so;/usr/lib64/libm.so;/usr/lib64/libnuma.so;/usr/lib64/librt.so;/usr/lib64/libnsl.so;/usr/lib64/libutil.so;/usr/lib64/libm.so `

Then run 

    `make HYPHYMPI`

And then run make install to install the software

`make install`

+   HYPHYMP(I) will be installed at  `/location/of/choice/bin`
+   libhyphy_mp.(so/dylib/dll) will be installed at `/location/of/choice/lib`
+   HyPhy's standard library of batchfiles will go into `/location/of/choice/lib/hyphy`

HYPHYGTEST isn't installed normally,
because it serves no utility outside of testing.

To test HyPhy, build with the  GTEST target and run ./HYPHYGTEST from the source directory.
`make GTEST`
`./HYPHYGTEST`
