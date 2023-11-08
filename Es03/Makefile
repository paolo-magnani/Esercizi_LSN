CC = g++
CFLAGS = -Wall -O3 --std=c++11

es3.exe : es3.o random.o libremia.o
	$(CC) libremia.o random.o es3.o -o es3.exe
es3.o : es3.cpp 
	$(CC) -c es3.cpp -o es3.o $(CFLAGS)
libremia.o : ../libremia.cpp ../libremia.h
	$(CC) -c ../libremia.cpp -o libremia.o $(CFLAGS)
random.o : ../random/random.cpp ../random/random.h
	$(CC) -c ../random/random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o es3.exe ../random/seed.out
