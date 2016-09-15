CXX=g++

hub: main.o
	${CXX} -I${eig_path} main.o -o hub -O2 -w -std=c++11 

main.o: main.cpp include/hf.h include/hubbard.h include/schmidt.h include/read.h \
	include/hred.h include/ab_string.h include/dfci.h include/scf.h include/diis.h \
	include/matdeque.h include/dmet.h include/bootstrap.h include/nropt.h include/brent.h
	${CXX} -I${eig_path} -c main.cpp -O2 -w -std=c++11 

clean:
	rm -f *.o hub

run: hub
	./hub
