CFLAGS=-O3 -std=c++11 #-msse -march=athlon

all:	test_nn knn

cover_tree.o: cover_tree.cc cover_tree.h point.h stack.h
	g++ -c -Wall $(CFLAGS) cover_tree.cc

point.o: point.cc point.h stack.h
	g++ -c -Wall $(CFLAGS) point.cc

test_nn:	point.o cover_tree.o stack.h point.h test_nn.cc
	g++ -Wall $(CFLAGS) -o test_nn test_nn.cc cover_tree.o point.o 

knn:	point.o cover_tree.o stack.h point.h knn.cc
	g++ -Wall $(CFLAGS) -o knn knn.cc cover_tree.o point.o 

clean:
	rm -f *.o test_nn knn
