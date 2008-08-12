
all:
	gcc -lm -Wall  -std=c99 -Wno-unused -lpthread -lz -lm -lhdf5 -lhdf5_hl  --pedantic Cluster.c -o cluster
	gcc -lm -Wall  -std=c99 -Wno-unused  -g -I/usr/include/python2.5/Include -I/usr/include/python2.5 -fpic -shared -o cluster.so

clean:
	rm -f cluster cluster.so
