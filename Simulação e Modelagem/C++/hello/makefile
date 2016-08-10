# program's name
NAME = hello
CLASS = GNR

# objects
OBJS = $(NAME).o $(CLASS).o

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

$(NAME).o: $(NAME).cpp $(CLASS).h
	$(CC) $(CFLAGS) $(NAME).cpp

$(CLASS).o: $(CLASS).h $(CLASS).cpp
	$(CC) $(CFLAGS) $(CLASS).cpp

clean:
	$(RM) *.o *~ $(NAME)