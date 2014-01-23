HyPhy - Hypothesis testing using Phylogenies
============================================

Installation
------------

HyPhy now uses CMake for its build system.
To install, make sure you have CMake installed,
then configure the project from the source directory using
`cmake .`
By default, CMake produces Makefiles (I think),
so if you prefer other build systems, such as Xcode,
configure using the -G switch, e.g.
`cmake -G Xcode .`
CMake has a number of build system generators,
feel free to peruse these and use them if you wish.

One should also be aware, HyPhy requires development libraries
for libcurl, its requirements, and libpthread.
Additionally, libcurl requires development libraries for
crypto++ and openssl (or gnutls, if your packages are so configured).
On Ubuntu, these are libcurl-dev, libcrypto++-dev, and libssl-dev.

By default, HyPhy installs into /usr/local,
but can be made to install anywhere by passing
`-DINSTALL_PREFIX=/wherever/you/want`
to cmake during the configuration, e.g.
`cmake -DINSTALL_PREFIX=/opt/hyphy .`.

Occasionally, you may have to specify which OSX SDK you are using. e.g.
`cmake -DCMAKE_OSX_SYSROOT=/Developer/SDKs/MacOSX10.9.sdk/ .`.

If you're on a UNIX-compatible system,
and you're comfortable with GNU make,
then just `make` away with one of the following targets:

+   MAC - build a Mac Carbon application
+   HYPHYGTK - HYPHY with GTK
+   SP - build a HyPhy executable (HYPHYSP) without multiprocessing
+   MP2 - build a HyPhy executable (HYPHYMP) using pthreads to do multiprocessing
+   MPI - build a HyPhy executable (HYPHYMPI) using MPI to do multiprocessing
+   LIB - build a HyPhy library (libhyphy_mp) using pthreads to do multiprocessing
-   GTEST - build HyPhy's gtest testing executable (HYPHYGTEST)

A subsequent `make install` should put everything where they belong:

+   HYPHYMP(I) goes into /installation/prefix/bin
+   libhyphy_mp.(so/dylib/dll) goes into /installation/prefix/lib
+   HyPhy's standard library of batchfiles goes into /installation/prefix/lib/hyphy

HYPHYGTEST isn't installed normally,
as it serves no utility outside of testing.
To test HyPhy,
build the GTEST target and run ./HYPHYGTEST from the source directory.
