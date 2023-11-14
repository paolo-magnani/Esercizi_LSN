CC = g++
CFLAGS = -Wall -O3 --std=c++11

es8.exe : es8.o random.o libremia.o
	$(CC) libremia.o random.o es8.o -o es8.exe
es8.o : es8.cpp es8.h
	$(CC) -c es8.cpp -o es8.o $(CFLAGS)
libremia.o : ../libremia.cpp ../libremia.h
	$(CC) -c ../libremia.cpp -o libremia.o $(CFLAGS)
random.o : ../random/random.cpp ../random/random.h
	$(CC) -c ../random/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o es8.exe ../random/seed.out
