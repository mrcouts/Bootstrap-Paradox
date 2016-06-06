# program's name
NAME = hello

#extra functions
FUNC = SomeUtilities

#classes
CLASS = Serial

# objects
OBJS = $(NAME).o $(FUNC).o $(CLASS).o RR.o 6R.o RK.o

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

$(NAME) : $(OBJS)
	$(LD) $(LFLAGS) $(OBJS) -o $(NAME) $(LIBS)

$(NAME).o: $(NAME).cpp $(FUNC).h $(CLASS).h RR.h 6R.h RK.h
	$(CC) $(CFLAGS) $(NAME).cpp

$(FUNC).o: $(FUNC).h $(FUNC).cpp
	$(CC) $(CFLAGS) $(FUNC).cpp

$(CLASS).o: $(CLASS).h $(CLASS).cpp $(FUNC).h
	$(CC) $(CFLAGS) $(CLASS).cpp

RR.o: RR.h RR.cpp
	$(CC) $(CFLAGS) RR.cpp

6R.o: 6R.h 6R.cpp
	$(CC) $(CFLAGS) 6R.cpp

RK.o: RK.h RK.cpp
	$(CC) $(CFLAGS) RK.cpp

clean:
	$(RM) *.o *~ $(NAME)