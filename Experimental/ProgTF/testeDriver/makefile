#nome do programa
NAME = testeDriver

#Arquivo externo
EXT = Utilities

#objetos
OBJS = $(NAME).o $(EXT).o

#compilador e linker
CC = gcc
LD = gcc

#flags
CFLAGS = -O3 -Wall -c
LFLAGS = -Wall

#libs
LIBS = -lm -lwiringPi

#comando remove
RM = rm -f

all: $(NAME)

$(NAME): $(OBJS)
	$(LD) $(LFLAGS) $(OBJS) -o $(NAME) $(LIBS)
	
$(NAME).o: $(NAME).c $(EXT).h
		$(CC) $(CFLAGS) $(NAME).c
		
$(EXT).o: $(EXT).h $(EXT).c
		$(CC) $(CFLAGS) $(EXT).c 
		
clean:
	$(RM) *.o *~ $(NAME)
