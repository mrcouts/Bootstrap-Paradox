# program's name
NAME = 6R_sim

#extra functions
FUNC = SomeUtilities

#classes
CLASS = Serial

# objects
OBJS = $(NAME).o $(FUNC).o $(CLASS).o 6R.o

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

$(NAME).o: $(NAME).cpp $(FUNC).h $(CLASS).h 6R.h
	$(CC) $(CFLAGS) $(NAME).cpp

$(FUNC).o: $(FUNC).h $(FUNC).cpp
	$(CC) $(CFLAGS) $(FUNC).cpp

$(CLASS).o: $(CLASS).h $(CLASS).cpp $(FUNC).h
	$(CC) $(CFLAGS) $(CLASS).cpp

6R.o: 6R.h 6R.cpp
	$(CC) $(CFLAGS) 6R.cpp

clean:
	$(RM) *.o *~ $(NAME).txt $(NAME) txt2py.pyc