CXX=g++
CXXFLAGS=-Ofast -w
CXXINC=${eig_path}

hub: main.o
	${CXX} -I$(CXXINC) $(CXXFLAGS) main.o -o hub

main.o: main.cpp include/hf.h include/hubbard.h include/schmidt.h include/read.h \
	include/hred.h include/ab_string.h include/dfci.h include/scf.h include/diis.h \
	include/matdeque.h include/dmet.h include/bootstrap.h include/nropt.h include/brent.h \
	include/frag.h include/troyfci.h include/fciwrap.h 
	${CXX} -I${CXXINC} $(CXXFLAGS) -c main.cpp 

clean:
	rm -f *.o hub

run: hub
	./hub
