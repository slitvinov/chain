.POSIX:
.SUFFIXES:
.SUFFIXES: .c .o

all: main double
0 = main.o energy.o force.o
1 = double.o energy.o
main: $0; $(CC) $(LDFLAGS) $0 -o $@ -lm
double: $1; $(CC) $(LDFLAGS) $1 -o $@ -lm
.c.o:
	$(CC) $(CFLAGS) -c $<
clean:
	rm -f main double $0 $1