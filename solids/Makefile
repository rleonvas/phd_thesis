CC = gcc -g
OBJECTS = main.o files.o functions.o vector.o set_geometry.o computing.o minimization.o
DEPS = files.h functions.h vector.h set_geometry.h computing.h minimization.h
CPLEX_INCLUDE = /home/apps/cplex/12.6.1/cplex/include/
CPLEX_LIBS = /home/apps/cplex/12.6.1/cplex/lib/x86-64_linux/static_pic/
CFLAGS = -fPIC

all: main clean

%.o: %.c $(DEPS)
	$(CC) -pg -O2 -c $< -I$(CPLEX_INCLUDE)

main: $(OBJECTS)
	$(CC) -pg -O2 -o $@ $^ $(CFLAGS) -I$(CPLEX_INCLUDE) -L$(CPLEX_LIBS) -lcplex -lm -lpthread


#%.o: %.c $(DEPS)
#	$(CC) -c $< 

#main: $(OBJECTS)
#	$(CC) -o $@ $^ -lm

clean:
	rm *.o
