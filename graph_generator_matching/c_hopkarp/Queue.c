#include "Queue.h"
#include <stdio.h>

void queue_init( Queue *q, int *data, int total )
{
  q->data = data;
  q->count = 0;
  q->front = 0;
  q->total = total;
}

void enqueue( Queue *q, int k )
{
  q->data[ (q->front + q->count++) % q->total ] = k;
}

int dequeue( Queue *q )
{
  int popped = q->data[q->front];
  q->front = ( q->front + 1 ) % q->total;
  q->count--;
  return popped;
}

int empty( Queue *q )
{
  return ( q->count == 0 );
}

int full( Queue *q )
{
  return ( q->count == q->total );
}

void queue_print( Queue *q )
{
  int counter = 0;
  int index = q->front;
  while( counter++ < q->count )
  {
    printf("%d ", q->data[index]);
    index = ( index + 1 ) % q->total;
  }
}
