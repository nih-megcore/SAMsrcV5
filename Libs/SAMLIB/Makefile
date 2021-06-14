include ../../config/Makefile.local

LIB = libsam
SRC = $(shell ls *.c)
OBJ = $(SRC:.c=.o)

include ../Makefile.lib

Dflip.o: Dflip.c
	$(CC) $(CFLAGS) -Wno-strict-aliasing -c -o Dflip.o Dflip.c

clean:
	rm -f *~ *.o *.a *.so core $(LIBDIR)/$(LIB).*
