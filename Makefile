CC = gcc
CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LDFLAGS = -lm
OBJ = symnmf.o matrix.o
 
symnmf: $(OBJ)
	$(CC) $(CFLAGS) -o symnmf $(OBJ) $(LDFLAGS)
 
symnmf.o: symnmf.c symnmf.h matrix.h
	$(CC) $(CFLAGS) -c symnmf.c
 
matrix.o: matrix.c matrix.h
	$(CC) $(CFLAGS) -c matrix.c
 
clean:
	rm -f *.o symnmf

.PHONY: clean