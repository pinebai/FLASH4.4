CC = gcc
CFLAGS = -ansi -pedantic -Wall -O3

sedov: sedov.c
	${CC} ${CFLAGS} sedov.c -lm -o $@

clean:
	- rm *~
	- rm sedov
	- rm *.dat
	- rm *.png
	- rm this.param
	- rm plotfile
