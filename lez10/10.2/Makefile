CC = g++
CFLAGS = -Wall -O3 --std=c++11
N=4

main.exe : main.o random.o popolazione.o
	mpicxx random.o popolazione.o main.o -o main.exe 

main.o : main.cpp popolazione.h random.h
	mpicxx -c main.cpp -o main.o 

popolazione.o: popolazione.cpp popolazione.h
	g++ -c popolazione.cpp $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

compila :
	$(CC) -c random.cpp -o random.o $(CFLAGS)
	g++ -c popolazione.cpp $(CFLAGS)
	mpicxx -c main.cpp -o main.o $(CFLAGS)
	mpicxx random.o main.o popolazione.o -o main.exe

esegui:
	mpiexec -np $(N) --oversubscribe main.exe $(N)

clean :
	rm *.o main.exe seed.out
