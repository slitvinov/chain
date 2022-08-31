.POSIX:
.SUFFIXES:
.SUFFIXES: .c .o

CFLAGS = -g -O2
L = libchain.a
O = \
energy.o\
force.o\

$L: $O; ar r $@ $O && ranlib $@
.c.o:; $(CC) $(CFLAGS) $< -c

clean:; -rm $L $O

force.o: force.h
force.o: energy.h
energy.o: energy.h
