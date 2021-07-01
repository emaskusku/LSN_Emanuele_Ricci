CC = g++
CFLAGS = -Wall -O3 --std=c++11

es3 : es3.o random.o
	$(CC) random.o es3.o -o es3 `root-config --cflags` `root-config --libs`

es3.o : es3.cpp
	$(CC) -c es3.cpp -o es3.o $(CFLAGS) `root-config --cflags` `root-config --libs`

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

compila:	es3
	
esegui:
	./es3

clean :
	rm *.o es3 seed.out
