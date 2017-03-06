#Makefile

all:
	icc -lrt -qopenmp main.c -o main

run:
	numactl --cpunodebind=0 --preferred=0 ./main

clean:
	rm hw0_output main
