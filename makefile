CC = gcc
FLAGS = -std=c++11 
LIBS = -lm -lstdc++ -L$(HOME)/local/lib
INC = -I./include -I$(HOME)/local/include

all: program

program: main.o
	$(CC) -o program main.o $(FLAGS) $(LIBS)

main.o: source/main.cpp
	$(CC) -o main.o -c source/main.cpp $(FLAGS) $(INC)

clean:
	rm -f program
	rm -f *.o

