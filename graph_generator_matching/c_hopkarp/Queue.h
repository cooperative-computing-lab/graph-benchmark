#ifndef QUEUE_H
#define QUEUE_H

typedef struct Queue
{
  int *data;
  int front;
  int count;
  int total;
} Queue;

void queue_init( Queue *q, int *data, int total );
void enqueue( Queue *q, int k );
int dequeue( Queue *q );
int empty( Queue *q );
int full( Queue *q );
void queue_print( Queue *q );

#endif