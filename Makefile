CXX=g++

hub: main.o
	${CXX} -I${eig_path} main.o -o hub -O2 -w

main.o: main.cpp include/hf.h include/hubbard.h include/schmidt.h include/read.h \
	include/hred.h
	${CXX} -I${eig_path} -c main.cpp -O2 -w 

clean:
	rm -f *.o hub

run: hub
	./hub
