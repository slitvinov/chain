.POSIX:
.SUFFIXES:
.SUFFIXES: .c .o

0 = main.o energy.o
main: $0
	$(CC) $(LDFLAGS) $0 -o $@ -lm
.c.o:
	$(CC) $(CFLAGS) -c $<
clean:
	rm -f main $0