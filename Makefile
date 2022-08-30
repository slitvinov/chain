.POSIX:
.SUFFIXES:
.SUFFIXES: .c .o

all: main double tbond
0 = main.o energy.o force.o
1 = double.o energy.o
2 = tbond.o energy.o force.o
main: $0; $(CC) $(LDFLAGS) $0 -o $@ -lm
double: $1; $(CC) $(LDFLAGS) $1 -o $@ -lm
tbond: $2; $(CC) $(LDFLAGS) $2 -o $@ -lm
.c.o:
	$(CC) $(CFLAGS) -c $<
clean:
	rm -f tbond main double $0 $1 $2