routematch.so: routematch.cpp routematch.hpp
	g++  -O3 -std=c++0x -lspatialindex -I. -shared -o routematch.so routematch.cpp
