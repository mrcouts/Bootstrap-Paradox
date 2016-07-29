# program's name
NAME = PRRP_sim

#extra functions
FUNC = SomeUtilities

#classes
CLASS = Serial
CLASS2 = Dy

# objects
OBJS = $(NAME).o $(FUNC).o $(CLASS).o $(CLASS2).o PRRP.o

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

$(NAME).o: $(NAME).cpp $(FUNC).h $(CLASS).h PRRP.h
	$(CC) $(CFLAGS) $(NAME).cpp

$(FUNC).o: $(FUNC).h $(FUNC).cpp
	$(CC) $(CFLAGS) $(FUNC).cpp

$(CLASS).o: $(CLASS).h $(CLASS).cpp $(CLASS2).h $(FUNC).h
	$(CC) $(CFLAGS) $(CLASS).cpp

$(CLASS2).o: $(CLASS2).h $(CLASS2).cpp
	$(CC) $(CFLAGS) $(CLASS2).cpp

6R.o: PRRP.h PRRP.cpp
	$(CC) $(CFLAGS) PRRP.cpp

clean:
	$(RM) *.o *~ $(NAME).txt $(NAME) txt2py.pyc