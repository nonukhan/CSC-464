
all: baseline

baseline: timekeeper.o baseline.c
	gcc -Wall timekeeper.o baseline.c -o baseline

timekeeper.o: timekeeper.c timekeeper.h
	gcc -c -Wall timekeeper.c

clean:
	rm -rf timekeeper.o baseline
