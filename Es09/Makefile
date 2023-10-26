CC = g++
CFLAGS = -Wall -O3 --std=c++11

es5.exe : es9.o random.o
	$(CC) random.o es9.o -o es9.exe
es3.o : es9.cpp 
	$(CC) -c es9.cpp -o es9.o $(CFLAGS)
random.o : ../random/random.cpp ../random/random.h
	$(CC) -c ../random/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o es9.exe ../random/seed.out
