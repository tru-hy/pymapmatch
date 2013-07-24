# pymapmatch - Map matching with Python interface

pymapmatch is a (work in progress) Python module for mapping
noisy trajectories (from eg. GPS) to predefined routes
(and hopefully in future to general graphs/maps).

The implementation itself is done in C++ for performance with
a thin C-wrapper for Python/NumPy interfacing.

## Install
### Quick
Make sure you have the dependencies (see below), go to the source root and run
	
	git clone <thecloneurl>
	cd libspatialindex
	make

No real installer at the moment, and may stay so.

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


## Usage

TODO
