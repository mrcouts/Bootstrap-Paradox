# program's name
NAME = hello

#extra functions
FUNC = SomeUtilities

# objects
OBJS = $(NAME).o $(FUNC).o

# compiler and linker
CC = g++
LD = g++

#debug
DEBUG =

# flags
CFLAGS = -O3 -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

# libs
LIBS = -lm -larmadillo

# remove command
RM = rm -f

$(NAME) : $(OBJS)
	$(LD) $(LFLAGS) $(OBJS) -o $(NAME) $(LIBS)

$(NAME).o: $(NAME).cpp $(FUNC).h
	$(CC) $(CFLAGS) $(NAME).cpp

$(FUNC).o: $(FUNC).h $(FUNC).cpp
	$(CC) $(CFLAGS) $(FUNC).cpp

clean:
	$(RM) *.o *~ $(NAME)