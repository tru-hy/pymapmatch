routematch.so: routematch.o
	g++ -Ofast -std=c++0x -lspatialindex -I. -shared -o $@ $<

routematch.cpp: routematch.lzz
	lzz $<

routematch.o: routematch.cpp routematch.h
	g++ -Ofast -std=c++0x -lspatialindex -I. -c -o $@ $<

.PHONY: memtest
memtest: routematch_memtest
	valgrind --tool=memcheck ./$<

routematch_memtest: routematch_memtest.cpp routematch.o
	g++ -Og -std=c++0x -lspatialindex -I. -o $@ $^
