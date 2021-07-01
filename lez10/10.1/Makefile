CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o popolazione.o
	$(CC) random.o popolazione.o main.o -o main.exe
main.o : main.cpp popolazione.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
popolazione.o: popolazione.cpp popolazione.h
	g++ -c popolazione.cpp
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
compila:
	$(CC) -c random.cpp -o random.o $(CFLAGS)
	g++ -c popolazione.cpp $(CFLAGS)
	$(CC) -c main.cpp -o main.o $(CFLAGS)
	$(CC) random.o popolazione.o main.o -o main.exe
clean :
	rm *.o main.exe seed.out
