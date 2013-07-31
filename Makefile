CC=clang++
CFLAGS=-O3 -fpic
LIBS=-lspatialindex

routematch.so: routematch.o
	$(CC) $(CFLAGS) -std=c++0x -I. $(LIBS) -shared -o $@ $<

routematch.cpp: routematch.lzz
	lzz $<

routematch.o: routematch.cpp routematch.h
	$(CC) $(CFLAGS) -std=c++0x -I. -c -o $@ $<

.PHONY: memtest
memtest: routematch_memtest
	valgrind --tool=memcheck --leak-check=full ./$<

routematch_memtest: routematch_memtest.cpp routematch.o
	$(CC) $(CFLAGS) -std=c++0x -I. $(LIBS) -o $@ $^
