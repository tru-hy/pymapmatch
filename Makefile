CFLAGS=-Ofast

routematch.so: routematch.o
	g++ $(CFLAGS) -std=c++0x -lspatialindex -I. -shared -o $@ $<

routematch.cpp: routematch.lzz
	lzz $<

routematch.o: routematch.cpp routematch.h
	g++ $(CFLAGS) -std=c++0x -lspatialindex -I. -c -o $@ $<

.PHONY: memtest
memtest: routematch_memtest
	valgrind --tool=memcheck --leak-check=full ./$<

routematch_memtest: routematch_memtest.cpp routematch.o
	g++ $(CFLAGS) -std=c++0x -lspatialindex -I. -o $@ $^
