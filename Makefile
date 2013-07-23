routematch.so: routematch.cpp routematch.hpp
	g++ -Ofast -std=c++0x -lspatialindex -I. -shared -o routematch.so routematch.cpp
