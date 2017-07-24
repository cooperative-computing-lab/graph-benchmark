#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "graph.h"

// Compute the clustering coefficient for each node in a graph
double *cluster_coeff(struct Graph *graph) {
	int neighbor;
	double *coeff;
	int degree, present;
	struct AdjList *adjlist;

	coeff = calloc(graph->V, sizeof *coeff);

	#pragma omp parallel for private(degree, present)
	for(int nodei = 0; nodei < graph->V; nodei++) {
		adjlist = graph->array + nodei;

		degree = 0;	// race condition? N.B.
		present = 0;

		// For each neighbor...
		for(struct AdjListNode *node1 = adjlist->head; node1; node1 = node1->next) {
			neighbor = node1->dest;

			// ...find our neighbors...
			for(struct AdjListNode *node2 = adjlist->head; node2; node2 = node2->next) {
				if(node2->dest == neighbor)
					continue;

				// ...in its neighbors.
				for(struct AdjListNode *node3 = graph->array[neighbor].head;
					node3;
					node3 = node3->next) {
					if(node2->dest == node3->dest) {
						//this should be done with a reduction N.B. (make a global copy and then write in a critical section)
						present++;
						break;
					}
				}
			}

			// Count the number of our neighbors (for the denominator)
			degree++;
		}

		// Compute the local clustering coefficient for this node
		coeff[nodei] = (double) present/(degree*(degree - 1));
	}

	return coeff;
}

