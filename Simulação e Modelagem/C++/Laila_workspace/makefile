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
CLASS8 = GNR2

# objects
OBJS = $(NAME).o $(FUNC).o $(CLASS).o $(CLASS2).o $(CLASS3).o $(CLASS4).o $(CLASS5).o $(CLASS6).o $(CLASS7).o $(CLASS8).o 6R.o PRRP.o RK.o

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

$(NAME).o: $(NAME).cpp $(FUNC).h $(CLASS).h $(CLASS2).h $(CLASS3).h $(CLASS4).h $(CLASS5).h $(CLASS6).h $(CLASS7).h $(CLASS8).h 6R.h PRRP.h RK.h
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

$(CLASS5).o: $(CLASS5).h $(CLASS5).cpp $(CLASS).h $(FUNC).h 6R.h PRRP.h
	$(CC) $(CFLAGS) $(CLASS5).cpp

$(CLASS6).o: $(CLASS6).h $(CLASS6).cpp
	$(CC) $(CFLAGS) $(CLASS6).cpp

$(CLASS7).o: $(CLASS7).h $(CLASS7).cpp $(CLASS5).h
	$(CC) $(CFLAGS) $(CLASS7).cpp

6R.o: 6R.h 6R.cpp
	$(CC) $(CFLAGS) 6R.cpp

PRRP.o: PRRP.h PRRP.cpp
	$(CC) $(CFLAGS) PRRP.cpp

RK.o: RK.h RK.cpp $(CLASS4).h $(CLASS7).h 
	$(CC) $(CFLAGS) RK.cpp

$(CLASS8).o: $(CLASS8).h $(CLASS8).cpp $(CLASS7).h RK.h
	$(CC) $(CFLAGS) $(CLASS8).cpp

clean:
	$(RM) *.o *~ $(NAME).txt $(NAME) txt2py.pyc