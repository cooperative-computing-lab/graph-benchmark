CFLAGS = -Wall -fopenmp -O1

Match: clean Queue.o Vector.o Graph.o HopKarp.o Match.o 
	gcc $(CFLAGS) Queue.o Vector.o Graph.o HopKarp.o Match.o -o Match

Match.o: Match.c HopKarp.h Graph.h Vector.h Queue.h
	gcc $(CFLAGS) -c Match.c

HopKarp.o: HopKarp.c Graph.h Vector.h Queue.h
	gcc $(CFLAGS) -c HopKarp.c

Graph.o: Graph.c Vector.h Queue.h
	gcc $(CFLAGS) -c Graph.c

Vector.o: Vector.c
	gcc $(CFLAGS) -c Vector.c

Queue.o: Queue.c
	gcc $(CFLAGS) -c Queue.c

clean:
	rm -f Match Match.o HopKarp.o Vector.o Graph.o Queue.o

run:
	./Match $(args)