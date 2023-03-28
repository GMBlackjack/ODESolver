CC     = gcc
CFLAGS =   -O2 -march=native -fopenmp -Wall -Wno-unused-variable -std=gnu99
#CFLAGS =  -O2 -march=native -g -fopenmp -Wall -Wno-unused-variable -std=gnu99
#CFLAGS =  -O2 -g -Wall -Wno-unused-variable -Wno-unknown-pragmas -std=gnu99
#CFLAGS =   -O2 -march=native -fopenmp -Wall -Wno-unused-variable -std=gnu99
INCLUDEDIRS = 
all: ODESolverComplicated4  main.o
main.o: main.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c main.c -o main.o

ODESolverComplicated4:   main.o
	$(CC)   main.o -o ODESolverComplicated4 -lm -lgomp

clean:
	rm -f *.o */*.o *~ */*~ ./#* *.txt *.dat *.avi *.png ODESolverComplicated4
