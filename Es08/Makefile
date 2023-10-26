CC = g++
CFLAGS = -Wall -O3 --std=c++11

es5.exe : es8.o random.o
	$(CC) random.o es8.o -o es8.exe
es3.o : es8.cpp 
	$(CC) -c es8.cpp -o es8.o $(CFLAGS)
random.o : ../random/random.cpp ../random/random.h
	$(CC) -c ../random/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o es8.exe ../random/seed.out
