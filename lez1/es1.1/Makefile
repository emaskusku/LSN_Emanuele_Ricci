CC = g++
CFLAGS = -Wall -O3 --std=c++11

lsn01.1 : lsn01.1.o random.o
	$(CC) random.o lsn01.1.o -o lsn01.1 

lsn01.1.o : lsn01.1.cpp
	$(CC) -c lsn01.1.cpp -o lsn01.1.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o lsn01.1 seed.out
