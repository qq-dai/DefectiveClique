# all: sdclique
# .PHONY : all
#CONFIG = -O3 -g -lpthread

CC = g++ "-std=c++11" 
CONFIG = -O3 -g

objects = algorithms.o polyEnum.o cliqueEnum.o main.o

.PHONY : clean

sdclique: $(objects)
	$(CC) -o sdclique $(objects) $(CONFIG)

%.o:%.cpp
	$(CC) -c $^ $(CONFIG)

clean:
	rm -f *.o sdclique