CXX=g++

hub: main.o
	${CXX} -I${eig_path} main.o -o hub -O2

main.o: main.cpp include/hf.h include/hubbard.h include/matdeque.h include/diis.h
	${CXX} -I${eig_path} -c main.cpp -O2

clean:
	rm -f *.o hub

run: hub
	./hub
