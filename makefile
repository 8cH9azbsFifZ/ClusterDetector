all:
	gcc -lm -Wall  -std=c99 -Wno-unused --pedantic Cluster.c -o cluster

clean:
	rm -f cluster
