

phys: phys.c
	mpicc -g -Wall -std=c99 phys.c -o phys -lncurses -lm 
r:
	mpiexec -np 2 ./phys
 
oldphys:
	gcc oldphys.c -o oldphys -lncurses -lm 

clean:
	rm phys
