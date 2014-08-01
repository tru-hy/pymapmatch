# pymapmatch - Map matching with Python interface

pymapmatch is a (work in progress) Python module for mapping
noisy trajectories to a map. Currently supports
openstreetmap data and is mostly usable with "noisy shapes"
with no time information (usually encountered in GTFS datasets).
Could be quite easily support also eg. GPS trajectories

The implementation itself is done in C++ for performance with
a swig wrapper for interfacing with Python.

### Dependencies
Probably works on some modern Unices (developed on Linux).
Won't work on Windows in the current form (mostly due to the Makefile and
ctypes wrappers).

#### Compilation
* A recent G++ (uses some C++0x/C++11 features, at least builds on 4.8.2)
* Boost graph library
* Boost geometry library
* TCMalloc
* SWIG

## Install
### Quick
Make sure you have the dependencies (see above), go to the source root and run
	
	make

After this there's a python module `osmmapmatch.py`, which requires also the
just built `_osmmapmatch.so` to work.

## Routematching

There's also a C++ module with a tiny C wrapper for matching GPS positions
to a predefined route.

### Dependencies
Probably works on some modern Unices (developed on Linux).
Won't work on Windows in the current form (mostly due to the Makefile and
ctypes wrappers).

#### Runtime
* Python 2.x (where x is about 7)
* NumPy
* [libspatialindex](http://libspatialindex.github.io/)

#### Compilation
* A recent G++ (uses some C++0x/C++11 features, at least builds on 4.8.0)
* [Eigen](http://eigen.tuxfamily.org) (included in the source)
* [libspatialindex](http://libspatialindex.github.io/)
* [LZZ](http://www.lazycplusplus.com/) if you want to make modifications.
  The repository includes the compiled .h and .cpp files (for the moment at least),
  so the Makefile should skip the lzz part.

### Building

The route matching library can be compiled using command:
	
	TODO

## Usage

TODO
