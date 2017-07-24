#include<stdio.h>
#include <omp.h>
#include<stdlib.h>
#include "graph.h"
#include <string.h>

double jaccard(struct AdjListNode *n1, struct AdjListNode *n2);
double * jaccard_all_pairs(struct Graph *graph )
{
	int count;
	//return NULL;
	#pragma omp parallel for
	for(int x = 0; x < graph->V; x++)
		for(int y = x+1; y < graph->V; y++)
		{
			struct AdjListNode *first =  graph->array[x].head;
			struct AdjListNode *second =  graph->array[y].head;
			double result = jaccard(first, second);
			if(result != 0.0)
			{
			 count++;
			}
//			printf("jaccard was: %f for %d, %d\n", jaccard(first, second),x,y);
		}
		printf("count is  %d ", count);
	return NULL;
}

double jaccard(struct AdjListNode *n1, struct AdjListNode *n2temp)
{
	int unionSize = 0;
	int intersectSize = 0;
	while(n1)//while we still have neighbors
	{
		struct AdjListNode *n2 = n2temp;
		int neighborBoth = 0;
		int count = 0;
		while(n2)	
		{
			count++;
			if(n2->dest == n1->dest)
			{
				neighborBoth=1;
			}
			n2 = n2->next;
		}
		if(neighborBoth)
		{
			intersectSize++;
			unionSize++;
		}
		else
		{
			unionSize++;
		}
		n1 = n1->next;
	}
	return (double)intersectSize/(double)unionSize;
}
		
 

