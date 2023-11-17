CC = g++
CFLAGS = -Wall -O3 --std=c++11

es9.exe : es9.o random.o libremia.o
	$(CC) libremia.o random.o es9.o -o es9.exe
es9.o : es9.cpp es9pointers.h
	$(CC) -c es9.cpp -o es9.o $(CFLAGS)
libremia.o : ../libremia.cpp ../libremia.h
	$(CC) -c ../libremia.cpp -o libremia.o $(CFLAGS)
random.o : ../random/random.cpp ../random/random.h
	$(CC) -c ../random/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o es9.exe ../random/seed.out
