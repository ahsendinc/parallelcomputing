#Makefile

all:
	icc  -c  timer.c
	icc -lrt -qopenmp timer.o main.c -o main

run:
	numactl --cpunodebind=0 --preferred=0 ./main

clean:
	rm hw0_output main
