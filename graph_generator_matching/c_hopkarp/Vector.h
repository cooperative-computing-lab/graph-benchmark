#ifndef VEC_H
#define VEC_H

#define VEC_INIT_CAP 2
#define vec_qpop(x) vec_delete(x, 0)

typedef struct Vector
{
	void **items;
	int capacity;
	int total;
} Vector;

void vec_add( Vector *, void * );
void vec_set( Vector *, int, void * );
void *vec_get( Vector *, int );
void *vec_back( Vector * );
void *vec_front( Vector * );
void vec_init( Vector *, int );
int vec_total( Vector * );
int vec_empty( Vector * );
void vec_add_int( Vector *, int );
void vec_delete( Vector *, int );
//void vec_qpop( Vector *v );
void vec_clear( Vector * );
void vec_free( Vector * );
void vec_print( Vector * );

#endif