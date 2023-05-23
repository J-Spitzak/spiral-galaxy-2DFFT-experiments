CFLAGS = -O 
CCFLAGS = -O 
FFLAGS = -O

all: StuffThatWorks

clean:
	rm -f *.o *.a core 2dfft

.f:
	g77 $(FFLAGS) -o $* $*.f -lm
	rm -f *.o

.c:
	cc -o $* $*.c -lm 
	rm -f *.o

.cpp:  
	g++ $(CCFLAGS) -o $* $*.cpp FileClass.cpp FileParse.cpp -lm
	rm -f *.o

.d:
	cc -o $* $*.c -L/usr/local/etc/pgplot -lpgplot -lX11 -lm 

StuffThatWorks:
	gcc -c fourn.c -o fourn.o
	gcc -c 2dfft.c -o 2dfft.o
	gcc 2dfft.o fourn.o -lm -o 2dfft

