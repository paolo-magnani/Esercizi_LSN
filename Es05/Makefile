CC = g++
CFLAGS = -Wall -O3 --std=c++11

es5.exe : es5.o random.o
	$(CC) random.o es5.o -o es5.exe
es3.o : es5.cpp 
	$(CC) -c es5.cpp -o es5.o $(CFLAGS)
random.o : ../random/random.cpp ../random/random.h
	$(CC) -c ../random/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o es5.exe ../random/seed.out
