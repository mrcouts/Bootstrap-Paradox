# program's name
NAME = hello

# objects
OBJS = $(NAME).o

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

$(NAME).o: $(NAME).cpp
	$(CC) $(CFLAGS) $(NAME).cpp

clean:
	$(RM) *.o *~ $(NAME)