#include <stdio.h>
#include <stdlib.h>
#include "Vector.h"

void vec_init( Vector *v, int cap )
{
  if( cap < 1 )
    cap = 1;
  v->capacity = cap;
  v->total = 0;
  v->items = malloc(sizeof(void *) * v->capacity);
}

int vec_total( Vector *v )
{
  return v->total;
}

int vec_empty( Vector *v )
{
  return (v->total == 0);
}

static void vec_resize( Vector *v, int capacity )
{
  void **items = realloc(v->items, sizeof(void *) * capacity);
  if (items)
  {
    v->items = items;
    v->capacity = capacity;
  }
}

void vec_add( Vector *v, void *item )
{
  if (v->capacity == v->total)
    vec_resize(v, v->capacity * 2);
  v->items[v->total++] = item;
}

void vec_add_int( Vector *v, int item )
{
  if (v->capacity == v->total)
    vec_resize(v, v->capacity * 2);
  int *tmp = malloc(sizeof(int));
  *tmp = item;
  v->items[v->total++] = tmp;
}

void vec_set( Vector *v, int index, void *item )
{
  if (index >= 0 && index < v->total)
    v->items[index] = item;
}

void* vec_get( Vector *v, int index )
{
  if (index >= 0 && index < v->total)
    return v->items[index];
  return NULL;
}

void* vec_back( Vector *v )
{   
  if(v->total > 0)
    return v->items[v->total - 1];
  return NULL;
}

void* vec_front( Vector *v )
{
  if(v->total > 0)
    return v->items[0];
  return NULL;
}

void vec_delete( Vector *v, int index )
{
  if (index < 0 || index >= v->total)
    return;

  free(v->items[index]);
  v->items[index] = NULL;
  int total = v->total - 1;

  for (int i = index; i < total; i++)
  {
     v->items[i] = v->items[i + 1];
  }
  v->items[v->total - 1] = NULL;
  v->total--;

  if (v->total > 0 && v->total == v->capacity / 4)
    vec_resize(v, v->capacity / 2);
}

/*
void vec_qpop( Vector *v )
{
  if( v->total > 0 )
  {
    free(v->items[0]);
    v->items[0] = NULL;
    v->items = &v->items[1];
  }
}
*/

void vec_free( Vector *v )
{
  for( int i = 0; i < v->total; i++ )
    free(v->items[i]);
  free(v->items);
}

void vec_print( Vector *v )
{
  for( int i = 0; i < v->total; i++ )
  {
    printf("%d ", *(int *)v->items[i]);
  }
  printf("\n");
}