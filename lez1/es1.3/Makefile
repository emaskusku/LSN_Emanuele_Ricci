CC = g++
CFLAGS = -Wall -O3 --std=c++11

lsn01.3 : lsn01.3.o random.o
	$(CC) random.o lsn01.3.o -o lsn01.3
lsn01.3.o : lsn01.3.cpp
	$(CC) -c lsn01.3.cpp -o lsn01.3.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

compila: lsn01.3

esegui:
	./lsn01.3
clean :
	rm *.o lsn01.3 seed.out