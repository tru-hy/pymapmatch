CC=clang++
CFLAGS=-O3 -fpic -std=c++11
ROUTEMATCHLIBS=-lspatialindex
MAPMATCHLIBS=-lreadosm -lproj

mapmatch_benchmark: mapmatch_benchmark.cpp osmmapmatch.hpp
	$(CC) $(CFLAGS) -I. $(MAPMATCHLIBS) -o $@ $<

routematch.so: routematch.o
	$(CC) $(CFLAGS) -I. $(ROUTEMAPLIBS) -shared -o $@ $<

routematch.cpp: routematch.lzz
	lzz $<

routematch.o: routematch.cpp routematch.h
	$(CC) $(CFLAGS) -I. -c -o $@ $<

.PHONY: memtest
memtest: routematch_memtest
	valgrind --tool=memcheck --leak-check=full ./$<

routematch_memtest: routematch_memtest.cpp routematch.o
	$(CC) $(CFLAGS) -I. $(LIBS) -o $@ $^
