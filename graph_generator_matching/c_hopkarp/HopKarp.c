#include "HopKarp.h"
#include "Vector.h"
#include "Graph.h"
#include "Queue.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <omp.h>


int hop_karp( Graph *graph, int **Uverts )
{
  // pairU[u] stores pair of u in matching where u
  // is a vertex on left side of Bipartite Graph.
  // If u doesn't have any pair, then pairU[u] is NIL
  //void **pairU = calloc(graph->l + 1, sizeof(void*));
  int *pairU = calloc(graph->l+1, sizeof(int));

  // pairV[v] stores pair of v in matching. If v
  // doesn't have any pair, then pairU[v] is NIL
  //void **pairV = calloc(graph->r + 1, sizeof(void*));
  int *pairV = calloc(graph->r+1, sizeof(int));

  // dist[u] stores distance of left side vertices
  // dist[u] is one more than dist[u'] if u is next
  // to u'in augmenting path
  int *dist = malloc(sizeof(int) * (graph->l + 1));

  int result = 0;

  while( bfs(graph, dist, pairU, pairV) )
  {
    int i;
    for( i = 1; i <= graph->l; i++ )
    {
      if( pairU[i] == 0 && dfs(i, graph, dist, pairU, pairV) )
        result++;
    }
  }

  *Uverts = pairU;
  free(dist);

  return result;
}


int bfs( Graph *graph, int *dist, int *pairU, int *pairV )
{
	Queue Q;
  int *data = malloc( sizeof(int) * (graph->l + 1) );
  queue_init( &Q, data, graph->l + 1 );

  int i;
  for( i = 1; i <= graph->l; i++ )
  {
    if( pairU[i] == 0 )
    {
      dist[i] = 0;
      enqueue(&Q, i);
    }
    else
       dist[i] = INT_MAX;
  }

  dist[0] = INT_MAX;
 
  while( !empty(&Q) )
  {
    int u = dequeue(&Q);

    if( dist[u] < dist[0] )
    {
      Vector *tmp = vec_get(&graph->alist, u - 1);

      for( int i = 0; i < tmp->total; i++ )
      {
        int v = *(int*)vec_get(tmp, i);
        
        if( dist[pairV[v]] == INT_MAX )
        {
          dist[pairV[v]] = dist[u] + 1;
          enqueue(&Q, pairV[v]);
        }
      }
    }
  }

  free(data);
  return dist[0] != INT_MAX;
}


int dfs( int u, Graph *graph, int *dist, int *pairU, int *pairV )
{
  if( u != 0 )
  {
    Vector *tmp = vec_get(&graph->alist, u - 1);

    for( int i = 0; i < tmp->total; i++ )
    {
      int v = *(int*)vec_get(tmp, i);
      if( dist[pairV[v]] == dist[u] + 1 )
      {
        if( dfs( pairV[v], graph, dist, pairU, pairV ) == 1 )
        {
          pairV[v] = u;
          pairU[u] = v;
          return 1;
        }
      }
    }

    dist[u] = INT_MAX;
    return 0;
  }

  return 1;
}


void print_matching( Graph *graph, int *Uverts )
{
  for(int i = 0; i < graph->l; i++)
  {
    if( i % 5 == 0 )
      printf("\n");
    printf("%-5d:", *(int*)vec_get(&graph->lverts, i));
    printf("%5d\t", Uverts[i+1]);
  }
  printf("\n\n");
}