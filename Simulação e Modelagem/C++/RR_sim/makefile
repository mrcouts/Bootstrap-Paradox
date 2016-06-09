# program's name
NAME = RR_sim

#extra functions
FUNC = SomeUtilities

#classes
CLASS = Serial

# objects
OBJS = $(NAME).o $(FUNC).o $(CLASS).o RR.o

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

$(NAME).o: $(NAME).cpp $(FUNC).h $(CLASS).h RR.h
	$(CC) $(CFLAGS) $(NAME).cpp

$(FUNC).o: $(FUNC).h $(FUNC).cpp
	$(CC) $(CFLAGS) $(FUNC).cpp

$(CLASS).o: $(CLASS).h $(CLASS).cpp $(FUNC).h
	$(CC) $(CFLAGS) $(CLASS).cpp

RR.o: RR.h RR.cpp
	$(CC) $(CFLAGS) RR.cpp

clean:
	$(RM) *.o *~ $(NAME).txt $(NAME) txt2py.pyc