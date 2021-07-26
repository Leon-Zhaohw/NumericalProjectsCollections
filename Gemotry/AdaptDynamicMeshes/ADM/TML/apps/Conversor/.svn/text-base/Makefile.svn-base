include ../TMLmake.defs

CC = g++
FLAGS = -g -Wall

INC = $(TML_INC)
LIB = $(TML_LIB)

BIN = convert components boundary isolated 3ds2off

.cc.o: $<
	$(CC) -c $< -o $@ $(FLAGS) $(INC)

all: $(BIN)

convert: convert.o
	$(CC) -o $@ convert.o $(FLAGS) $(INC) $(LIB)

components: components.o
	$(CC) -o $@ components.o $(FLAGS) $(INC) $(LIB)

boundary: boundary.o
	$(CC) -o $@ boundary.o $(FLAGS) $(INC) $(LIB)

isolated: isolated.o
	$(CC) -o $@ isolated.o $(FLAGS) $(INC) $(LIB)

3ds2off: 3ds2off.o
	$(CC) -o $@ 3ds2off.o $(FLAGS) $(INC) $(LIB)

clean:
	rm -rf *.o *~ $(BIN)
