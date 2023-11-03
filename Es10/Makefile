CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

es10.exe : es10.o random.o
	$(CC) random.o es10.o -o es10.exe
es10.o : es10.cpp 
	$(CC) -c es10.cpp -o es10.o $(CFLAGS)
random.o : ../random/random.cpp ../random/random.h
	$(CC) -c ../random/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o es10.exe ../random/seed.out