CC=g++
CFLAGS=-O3 -fPIC -std=c++11 -fopenmp -Wl,--no-as-needed
ROUTEMATCHLIBS=-lspatialindex
MAPMATCHLIBS=-lreadosm -lproj -ltcmalloc
PYTHONCONF=`python2-config --includes`

DEBUG ?= 0
ifeq ($(DEBUG), 1)
        CFLAGS += -DDEBUG -g -O0
endif

_osmmapmatch.so: osmmapmatch_wrap.cxx
	$(CC) $(CFLAGS) $(MAPMATCHLIBS) $(PYTHONCONF) -shared -o $@ $<

osmmapmatch_wrap.cxx osmmapmatch.py: osmmapmatch.i osmmapmatch.hpp
	swig -python -threads -c++ $<

mapmatch_benchmark: mapmatch_benchmark.cpp osmmapmatch.hpp
	$(CC) -g $(CFLAGS) -I. $(MAPMATCHLIBS) -o $@ $<

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
