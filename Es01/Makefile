CC = g++
CFLAGS = -Wall -O3 --std=c++11

es1.exe : es1.o random.o libremia.o
	$(CC) libremia.o random.o es1.o -o es1.exe
es1.o : es1.cpp
	$(CC) -c es1.cpp -o es1.o $(CFLAGS)
libremia.o : ../libremia.cpp ../libremia.h
	$(CC) -c ../libremia.cpp -o libremia.o $(CFLAGS)
random.o : ../random/random.cpp ../random/random.h
	$(CC) -c ../random/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o es1.exe ../random/seed.out
