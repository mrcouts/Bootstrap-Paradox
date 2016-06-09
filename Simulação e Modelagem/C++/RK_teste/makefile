# program's name
NAME = RK_teste

# objects
OBJS = $(NAME).o RK.o

# compiler and linker
CC = g++
LD = g++

#debug
DEBUG =

# flags
CFLAGS = -O3 -std=c++11 -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

# libs
LIBS = -lm -larmadillo

# remove command
RM = rm -f

all: $(NAME)

$(NAME) : $(OBJS)
	$(LD) $(LFLAGS) $(OBJS) -o $(NAME) $(LIBS)

$(NAME).o: $(NAME).cpp RK.h
	$(CC) $(CFLAGS) $(NAME).cpp

RK.o: RK.h RK.cpp
	$(CC) $(CFLAGS) RK.cpp

clean:
	$(RM) *.o *~ $(NAME).txt $(NAME) txt2py.pyc