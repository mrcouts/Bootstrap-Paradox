# program's name
NAME0 = hello
NAME1 = RK_teste
NAME2 = RR_sim
NAME3 = 6R_sim

#extra functions
FUNC = SomeUtilities

#classes
CLASS = Serial

# objects
OBJS0 = $(NAME0).o
OBJS1 = $(NAME1).o RK.o
OBJS2 = $(NAME2).o $(FUNC).o $(CLASS).o RR.o
OBJS3 = $(NAME3).o $(FUNC).o $(CLASS).o 6R.o

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

all: $(NAME0) $(NAME1) $(NAME2) $(NAME3)

$(NAME0) : $(OBJS0)
	$(LD) $(LFLAGS) $(OBJS0) -o $(NAME0) $(LIBS)

$(NAME1) : $(OBJS1)
	$(LD) $(LFLAGS) $(OBJS1) -o $(NAME1) $(LIBS)

$(NAME2) : $(OBJS2)
	$(LD) $(LFLAGS) $(OBJS2) -o $(NAME2) $(LIBS)

$(NAME3) : $(OBJS3)
	$(LD) $(LFLAGS) $(OBJS3) -o $(NAME3) $(LIBS)

$(NAME0).o: $(NAME0).cpp
	$(CC) $(CFLAGS) $(NAME0).cpp

$(NAME1).o: $(NAME1).cpp RK.h
	$(CC) $(CFLAGS) $(NAME1).cpp

$(NAME2).o: $(NAME2).cpp $(FUNC).h $(CLASS).h RR.h
	$(CC) $(CFLAGS) $(NAME2).cpp

$(NAME3).o: $(NAME3).cpp $(FUNC).h $(CLASS).h 6R.h
	$(CC) $(CFLAGS) $(NAME3).cpp

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
	$(RM) *.o *~ $(NAME1) $(NAME2) $(NAME3)