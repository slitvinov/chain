.POSIX:
.SUFFIXES:
.SUFFIXES: .c .o

CHAIN_FLAGS = -static -Lsrc -lchain -lm -Isrc
M = \
double\
main\
tangle\
tbond\
tdihedral\

all: $M
.c:; $(CC) -o $@ $(CFLAGS) $< $(LDFLAGS) $(CHAIN_FLAGS)
clean:; -rm $M
