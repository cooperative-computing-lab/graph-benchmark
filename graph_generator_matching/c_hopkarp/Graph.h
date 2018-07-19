#ifndef GRAPH_H
#define GRAPH_H

#include "Vector.h"

typedef struct Graph
{
	int l;    /* number of vertices in left set */
  int r;    /* number of vertices in right set */
	int m;    /* number of edges */
	Vector lverts;   /* Vector of left vertices */
	Vector alist;   /* Pointers to vectors of vertices in right set*/
} Graph;

Graph* graph_create();
void graph_destroy( Graph *g );
void graph_add_edge( Graph *g, int u, int v );
int graph_has_edge( Graph *g, int source, int dest );
void graph_add_lvert( Graph *g, void *item );
void graph_add_rverts( Graph *g, Vector *v );
void graph_print( Graph *g );
void graph_read_adj_list( Graph *g, char *filename );
void str_app( char *str, char ch );

#endif