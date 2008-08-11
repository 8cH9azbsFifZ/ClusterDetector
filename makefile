all:
	gcc -lm -Wall  -std=c99 -Wno-unused -lpthread -lz -lm -lhdf5 -lhdf5_hl  --pedantic Cluster.c -o cluster

clean:
	rm -f cluster
