CC     = gcc
CFLAGS =   -O2 -march=native -fopenmp -Wall -Wno-unused-variable -std=gnu99
#CFLAGS =  -O2 -march=native -g -fopenmp -Wall -Wno-unused-variable -std=gnu99
#CFLAGS =  -O2 -g -Wall -Wno-unused-variable -Wno-unknown-pragmas -std=gnu99
#CFLAGS =   -O2 -march=native -fopenmp -Wall -Wno-unused-variable -std=gnu99
INCLUDEDIRS = 
all: ODESolver  main.o
main.o: main.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c main.c -o main.o

ODESolver:   main.o
	$(CC)   main.o -o ODESolver -lm -lgomp

clean:
	rm -f *.o */*.o *~ */*~ ./#* *.txt *.dat *.avi *.png ODESolver
