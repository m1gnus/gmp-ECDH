IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR) -Wall

ODIR=obj

_DEPS = ecdh.h mini-gmp.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = ecdh.o mini-gmp.o main.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: ecdh

ecdh: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ ecdh
