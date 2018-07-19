#include "HopKarp.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main( int argc, char **argv )
{
  int *Uverts;
  clock_t start, end;
  double cpu_time_used;

  Graph *graph = graph_create( argv[1] );
  graph_print( graph );

  start = clock();
  int match_count = hop_karp( graph, &Uverts );
  end = clock();
  cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
  printf("Total time to complete matching: %f\n", cpu_time_used);

  printf("Total number of matches: %d\n", match_count);
  //print_matching( graph, Uverts );

  free(Uverts);
	graph_destroy( graph );
}
