.POSIX:
.SUFFIXES:
.SUFFIXES: .c .o

all: main double tangle tbond tdihedral
0 = main.o energy.o force.o
1 = double.o energy.o
2 = tbond.o energy.o force.o
3 = tangle.o energy.o force.o
4 = tdihedral.o energy.o force.o
main: $0; $(CC) $(LDFLAGS) $0 -o $@ -lm
double: $1; $(CC) $(LDFLAGS) $1 -o $@ -lm
tbond: $2; $(CC) $(LDFLAGS) $2 -o $@ -lm
tangle: $3; $(CC) $(LDFLAGS) $3 -o $@ -lm
tdihedral: $4; $(CC) $(LDFLAGS) $4 -o $@ -lm
.c.o:
	$(CC) $(CFLAGS) -c $<
clean:
	rm -f tanble tbond tdihedral main double $1 $2 $3 $$
