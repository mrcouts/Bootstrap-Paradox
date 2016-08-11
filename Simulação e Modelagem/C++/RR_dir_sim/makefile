# program's name
NAME = RR_sim

#extra functions
FUNC = SomeUtilities

#classes
CLASS = Serial
CLASS2 = Dy
CLASS3 = FLControlLaw
CLASS4 = Acceleration
CLASS5 = Parallel
CLASS6 = Reference
CLASS7 = GNR

# objects
OBJS = $(NAME).o $(FUNC).o $(CLASS).o $(CLASS2).o $(CLASS3).o $(CLASS4).o $(CLASS5).o $(CLASS6).o $(CLASS7).o RR.o RK.o

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

$(NAME).o: $(NAME).cpp $(FUNC).h $(CLASS).h $(CLASS2).h $(CLASS3).h $(CLASS4).h $(CLASS5).h $(CLASS6).h $(CLASS7).h RR.h RK.h
	$(CC) $(CFLAGS) $(NAME).cpp

$(FUNC).o: $(FUNC).h $(FUNC).cpp
	$(CC) $(CFLAGS) $(FUNC).cpp

$(CLASS).o: $(CLASS).h $(CLASS).cpp $(CLASS2).h $(FUNC).h
	$(CC) $(CFLAGS) $(CLASS).cpp

$(CLASS2).o: $(CLASS2).h $(CLASS2).cpp
	$(CC) $(CFLAGS) $(CLASS2).cpp

$(CLASS3).o: $(CLASS3).h $(CLASS3).cpp $(CLASS2).h $(CLASS5).h $(CLASS6).h
	$(CC) $(CFLAGS) $(CLASS3).cpp

$(CLASS4).o: $(CLASS4).h $(CLASS4).cpp $(CLASS).h $(CLASS2).h $(CLASS3).h $(CLASS5).h $(CLASS6).h
	$(CC) $(CFLAGS) $(CLASS4).cpp

$(CLASS5).o: $(CLASS5).h $(CLASS5).cpp $(CLASS).h $(FUNC).h RR.h
	$(CC) $(CFLAGS) $(CLASS5).cpp

$(CLASS6).o: $(CLASS6).h $(CLASS6).cpp
	$(CC) $(CFLAGS) $(CLASS6).cpp

$(CLASS7).o: $(CLASS7).h $(CLASS7).cpp
	$(CC) $(CFLAGS) $(CLASS7).cpp

RR.o: RR.h RR.cpp
	$(CC) $(CFLAGS) RR.cpp

RK.o: RK.h RK.cpp $(CLASS4).h $(CLASS7).h
	$(CC) $(CFLAGS) RK.cpp

clean:
	$(RM) *.o *~ $(NAME).txt $(NAME) txt2py.pyc